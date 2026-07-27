#' Causal inference with a latent IRT confounder via BART and WALNUTS
#'
#' Runs a Gibbs sampler for a causal model in which one or more latent person
#' traits `theta`, each measured by its own two-parameter logistic (2PL)
#' item-response model, confound a binary treatment `z` and a continuous outcome
#' `y`. Within each scan:
#'
#' * Each trait gets one per-person random-walk Metropolis move (the `theta_j`
#'   are conditionally independent given the trees and item parameters, so each
#'   person is accepted or rejected on its own). The traits are visited in a
#'   systematic scan --- trait `k` is proposed and installed for every person
#'   before trait `k + 1` is touched. The accepted values are installed into both
#'   BART models by a joint per-observation sweep
#'   ([dbarts::updatePredictorPerObservationJointly()]): a person is kept only
#'   if its `theta` leaves every leaf non-empty in every tree of both models.
#' * The two BART surfaces --- `y ~ f(z, theta)` (response) and
#'   `z ~ f(theta)` (probit assignment) --- each take one tree-update step.
#' * The IRT item parameters `(alpha, beta)` of each bank are drawn conditional
#'   on that bank's trait by the WALNUTS sampler ([irt_item_sampler()]):
#'   adapting during burn-in, then frozen for the sampling phase. Conditional on
#'   the traits the banks are independent, so this is one sampler per bank, not
#'   one larger one.
#'
#' The treatment-effect estimand `E[f(1, theta) - f(0, theta)]` is accumulated
#' from the response surface each sampling scan.
#'
#' @param responses Either a `n_persons` x `n_items` matrix (or data frame) of
#'   0/1 item responses --- the single-trait case --- or a list of such matrices,
#'   one item bank per latent trait. Every item loads on exactly one trait; banks
#'   may differ in item count but must cover the same persons in the same order.
#' @param y Numeric outcome, length `n_persons`.
#' @param z Binary (0/1) treatment, length `n_persons`.
#' @param n_burnin Number of burn-in scans. During burn-in the proposal SD is
#'   adapted and empty-leaf `theta` moves are collapsed rather than rejected.
#' @param n_sampling Number of kept scans.
#' @param theta_sd Initial standard deviation of the per-person `theta` proposal;
#'   adapted toward `theta_accept_target` during burn-in and frozen afterward.
#'   One shared value across traits.
#' @param theta_accept_target Target per-person acceptance rate for the `theta`
#'   proposal adaptation.
#' @param warmup_start Burn-in scan after which the WALNUTS item samplers begin
#'   adapting (before this they are held at their initial `(alpha, beta)` so the
#'   tuning is not polluted by the not-yet-settled `theta`). Must be less than
#'   `n_burnin`, otherwise the item samplers would never adapt; if it is not, it
#'   is reduced to `n_burnin %/% 2` with a warning.
#' @param beta_sd Prior standard deviation for the IRT item difficulties.
#' @param step_size Initial WALNUTS leapfrog step size.
#' @param n_theta_cutpoints Number of cut points for each trait predictor in
#'   each BART model (interior quantiles of a standard normal).
#' @param n_trees Number of trees in each BART surface.
#' @param n_threads Number of threads for the BART updates.
#' @param seed Integer seed. This calls [set.seed()] internally, so it
#'   overwrites the global RNG state (`.Random.seed`) as a side effect; it also
#'   seeds each WALNUTS sampler's own (independent) RNG, bank `k` getting
#'   `seed + k - 1`. `y`, `z`, and `responses` may be numeric or logical, but not
#'   factors; `y` must not contain `NA`.
#' @param keep_theta If `TRUE`, also return the per-scan `theta` draws. Off by
#'   default to keep the result small.
#' @param verbose If `TRUE`, print progress periodically.
#'
#' @return A list of draws and diagnostics:
#'   \describe{
#'     \item{`ate`}{Length-`n_sampling` vector of treatment-effect draws.}
#'     \item{`alpha`, `beta`}{`n_sampling` x `n_items` matrices of item-parameter
#'       draws; with more than one trait, a list of one such matrix per bank.}
#'     \item{`sigma`}{Length-`n_sampling` vector of response residual SD draws.}
#'     \item{`theta`}{`n_sampling` x `n_persons` matrix of `theta` draws (a list
#'       of one such matrix per trait when there is more than one), if
#'       `keep_theta = TRUE`.}
#'     \item{`theta_accept`}{Per-scan `theta` acceptance rate, averaged over
#'       traits (length `n_burnin + n_sampling`).}
#'     \item{`theta_sd`}{Final (frozen) proposal SD.}
#'     \item{`theta_sd_trace`}{Per-scan adapted proposal SD.}
#'     \item{`tuning`}{Final WALNUTS tuning (see [irt_tuning()]); with more than
#'       one trait, a list of one per bank.}
#'   }
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 200, n_items = 30, ate = -0.2, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 100, n_sampling = 200, seed = 1)
#' mean(fit$ate)                      # posterior mean treatment effect (~ -0.2)
#' quantile(fit$ate, c(0.025, 0.975)) # 95% interval
#'
#' # two latent traits, two item banks: pass a list of response matrices
#' sim2 <- simulate_irt_causal(
#'   n_persons = 200, n_items = c(20, 30), n_traits = 2,
#'   ate = -0.2, prognostic = c(1.2, 0.8), seed = 1
#' )
#' fit2 <- irt_causal_bart(sim2$responses, sim2$y, sim2$z,
#'                         n_burnin = 100, n_sampling = 200, seed = 1)
#' mean(fit2$ate)
#' length(fit2$alpha)                 # one item-parameter matrix per bank
#'
#' @seealso [simulate_irt_causal()], [irt_item_sampler()]
#' @importFrom stats dnorm pnorm plogis qnorm reformulate rnorm rexp
#' @export
irt_causal_bart <- function(
  responses,
  y,
  z,
  n_burnin = 500L,
  n_sampling = 1000L,
  theta_sd = 0.6,
  theta_accept_target = 0.44,
  warmup_start = 50L,
  beta_sd = 10,
  step_size = 0.1,
  n_theta_cutpoints = 100L,
  n_trees = 75L,
  n_threads = 1L,
  seed = 1L,
  keep_theta = FALSE,
  verbose = FALSE
) {
  banks <- as_response_banks(responses)
  n_traits <- length(banks)
  n_persons <- nrow(banks[[1L]])
  n_items <- vapply(banks, ncol, integer(1L))
  y <- as.double(y)
  z <- as.double(z)
  if (length(y) != n_persons) {
    stop("'y' must have length n_persons = ", n_persons)
  }
  if (length(z) != n_persons) {
    stop("'z' must have length n_persons = ", n_persons)
  }
  if (anyNA(y)) {
    stop("'y' must not contain NA")
  }
  if (anyNA(z) || any(z != 0 & z != 1)) {
    stop("'z' must be a 0/1 treatment indicator")
  }
  n_burnin <- as.integer(n_burnin)
  n_sampling <- as.integer(n_sampling)
  if (n_burnin < 0L || n_sampling <= 0L) {
    stop("'n_burnin' must be >= 0 and 'n_sampling' > 0")
  }
  n_theta_cutpoints <- as.integer(n_theta_cutpoints)
  if (is.na(n_theta_cutpoints) || n_theta_cutpoints < 1L) {
    stop("'n_theta_cutpoints' must be >= 1")
  }
  if (as.integer(n_trees) < 1L) {
    stop("'n_trees' must be >= 1")
  }

  # warmup_start delays WALNUTS adaptation until theta has settled; if it leaves
  # no adapting scans the item sampler would freeze at its initial tuning, so
  # clamp it below n_burnin (a no-op when n_burnin == 0: there is no warm-up).
  warmup_start <- as.integer(warmup_start)
  if (is.na(warmup_start) || warmup_start < 0L) {
    stop("'warmup_start' must be >= 0")
  }
  if (n_burnin > 0L && warmup_start >= n_burnin) {
    new_warmup_start <- n_burnin %/% 2L
    warning(sprintf(
      paste0(
        "'warmup_start' (%d) >= 'n_burnin' (%d): the WALNUTS item sampler ",
        "would never adapt; reducing 'warmup_start' to %d."
      ),
      warmup_start,
      n_burnin,
      new_warmup_start
    ))
    warmup_start <- new_warmup_start
  }

  # One BART predictor column per trait. A single trait keeps the bare name
  # "theta", so the K = 1 model, formulas, and draws are exactly as before.
  theta_names <- if (n_traits == 1L) {
    "theta"
  } else {
    paste0("theta", seq_len(n_traits))
  }

  set.seed(seed)
  theta <- matrix(rnorm(n_persons * n_traits), n_persons, n_traits)

  cutpoints <- qnorm(seq(0, 1, length.out = n_theta_cutpoints + 2L)[
    -c(1L, n_theta_cutpoints + 2L)
  ])
  control <- dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = as.integer(n_threads),
    n.trees = as.integer(n_trees)
  )

  # --- response surface (BART): y ~ f(z, theta) ---------------------------
  response_model <- dbarts::dbarts(
    reformulate(c("z", theta_names), response = "y"),
    data.frame(y = y, z = z, trait_frame(theta, theta_names)),
    control = control
  )
  for (nm in theta_names) {
    response_model$setCutPoints(cutpoints, nm)
  }
  response_model$sampleTreesFromPrior()
  response_samples <- response_model$run(5L, 1L)

  # --- assignment model (BART): z ~ f(theta) ------------------------------
  assignment_model <- dbarts::dbarts(
    reformulate(theta_names, response = "z"),
    data.frame(z = z, trait_frame(theta, theta_names)),
    control = control
  )
  for (nm in theta_names) {
    assignment_model$setCutPoints(cutpoints, nm)
  }
  assignment_model$sampleTreesFromPrior()
  assignment_samples <- assignment_model$run(5L, 1L)

  # --- IRT item parameters: sampled by WALNUTS (alpha, beta | theta) ------
  # Conditional on the traits the banks factorize exactly, so each gets its own
  # independent sampler over its own response matrix.
  alpha <- lapply(n_items, function(p) rep(1, p)) # discrimination (positive)
  beta <- lapply(n_items, function(p) rnorm(p)) # difficulty
  ab_sampler <- lapply(seq_len(n_traits), function(k) {
    irt_item_sampler(
      banks[[k]],
      alpha[[k]],
      beta[[k]],
      beta_sd = beta_sd,
      step_size = step_size,
      seed = as.integer(seed) + (k - 1L)
    )
  })

  alpha_mat <- lapply(seq_len(n_traits), function(k) {
    matrix(alpha[[k]], n_persons, n_items[k], byrow = TRUE)
  })
  beta_mat <- lapply(seq_len(n_traits), function(k) {
    matrix(beta[[k]], n_persons, n_items[k], byrow = TRUE)
  })
  sign_mat <- lapply(banks, function(bank) 2 * bank - 1)

  response_fitted <- response_samples$train
  assignment_fitted <- assignment_samples$train

  # storage (sampling phase only)
  ate_draws <- numeric(n_sampling)
  alpha_draws <- lapply(n_items, function(p) matrix(NA_real_, n_sampling, p))
  beta_draws <- lapply(n_items, function(p) matrix(NA_real_, n_sampling, p))
  sigma_draws <- numeric(n_sampling)
  theta_draws <- if (keep_theta) {
    lapply(seq_len(n_traits), function(k) {
      matrix(NA_real_, n_sampling, n_persons)
    })
  } else {
    NULL
  }
  theta_accept <- numeric(n_burnin + n_sampling)
  theta_sd_trace <- numeric(n_burnin + n_sampling)

  if (n_burnin == 0L) {
    for (s in ab_sampler) {
      irt_freeze(s)
    }
  }

  n_total <- n_burnin + n_sampling
  for (i_sample in seq_len(n_total)) {
    sampling_phase <- i_sample > n_burnin

    ## --- theta: per-person random-walk Metropolis, trait by trait ---------
    # Systematic scan over traits: trait k is proposed, accepted, and installed
    # for every person before trait k + 1 is touched, so each move is an
    # ordinary Metropolis-within-Gibbs step conditional on the other traits at
    # their current values. The component-wise shape is forced by the joint
    # install below, which takes a single shared column; see
    # docs/design/multi-trait.md, "The one real constraint".
    accept_rate <- numeric(n_traits) # post-install, for reporting
    proposal_rate <- numeric(n_traits) # pre-install, for the SD adaptation
    fitted_cur <- response_fitted
    assignment_cur <- assignment_fitted

    for (k in seq_len(n_traits)) {
      theta_old_k <- theta[, k]
      theta_prop_k <- theta_old_k + rnorm(n_persons, sd = theta_sd)
      theta_prop <- theta
      theta_prop[, k] <- theta_prop_k
      prop_frame <- trait_frame(theta_prop, theta_names)

      response_prop <- response_model$predict(
        data.frame(z = z, prop_frame)
      )
      assignment_prop <- assignment_model$predict(prop_frame)

      # per-person IRT log-likelihood for bank k (sum over that person's items)
      p_k <- n_items[k]
      irt_ll_cur <- rowSums(plogis(
        sign_mat[[k]] *
          alpha_mat[[k]] *
          (matrix(theta_old_k, n_persons, p_k) - beta_mat[[k]]),
        log.p = TRUE
      ))
      irt_ll_prop <- rowSums(plogis(
        sign_mat[[k]] *
          alpha_mat[[k]] *
          (matrix(theta_prop_k, n_persons, p_k) - beta_mat[[k]]),
        log.p = TRUE
      ))

      # per-person log acceptance ratio (vector of length n_persons)
      lr <- (irt_ll_prop +
        dnorm(y, response_prop, response_samples$sigma, log = TRUE) +
        pnorm((2 * z - 1) * assignment_prop, log.p = TRUE) +
        dnorm(theta_prop_k, log = TRUE)) -
        (irt_ll_cur +
          dnorm(y, fitted_cur, response_samples$sigma, log = TRUE) +
          pnorm((2 * z - 1) * assignment_cur, log.p = TRUE) +
          dnorm(theta_old_k, log = TRUE))

      accept <- !is.na(lr) & (-rexp(n_persons)) <= lr
      theta[accept, k] <- theta_prop_k[accept]
      proposal_rate[k] <- mean(accept)
      accept_rate[k] <- proposal_rate[k]

      ## --- install trait k into both BART models --------------------------
      if (sampling_phase) {
        # A theta that empties a leaf has zero tree prior, so the move must be
        # rejected. The joint sweep installs each person's theta in BOTH models
        # or neither -- a systematic-scan Metropolis-within-Gibbs sweep that
        # keeps the two models' trait columns identical and empty-leaf-free.
        installed <- dbarts::updatePredictorPerObservationJointly(
          list(response_model, assignment_model),
          theta[, k],
          theta_names[k]
        )
        theta[!installed, k] <- theta_old_k[!installed]
        accept_rate[k] <- mean(theta[, k] != theta_old_k)
      } else {
        # Burn-in: collapse empty leaves so theta and the trees stay in sync
        # while the chain/adaptation move. These draws are discarded.
        response_model$setPredictor(
          theta[, k],
          theta_names[k],
          forceUpdate = TRUE
        )
        assignment_model$setPredictor(
          theta[, k],
          theta_names[k],
          forceUpdate = TRUE
        )
      }

      # The next trait's ratio conditions on this one's accepted values, so the
      # "current" BART terms have to be re-read at the updated theta. Skipped
      # for the last trait (and so entirely when there is only one).
      if (k < n_traits) {
        cur_frame <- trait_frame(theta, theta_names)
        fitted_cur <- response_model$predict(data.frame(z = z, cur_frame))
        assignment_cur <- assignment_model$predict(cur_frame)
      }
    }

    p_accept <- mean(proposal_rate)
    theta_accept[i_sample] <- mean(accept_rate)

    ## --- adapt theta proposal SD toward target (burn-in only) -------------
    # Robbins-Monro on log(theta_sd) with decreasing gain; frozen once sampling
    # begins, so kept draws come from a fixed kernel.
    if (!sampling_phase) {
      gamma_t <- i_sample^(-0.6)
      theta_sd <- exp(
        log(theta_sd) + gamma_t * (p_accept - theta_accept_target)
      )
      theta_sd <- min(max(theta_sd, 1e-3), 10)
    }
    theta_sd_trace[i_sample] <- theta_sd

    ## --- BART trees | theta (one Gibbs step each) -------------------------
    response_samples <- response_model$run(0L, 1L)
    response_fitted <- response_samples$train
    assignment_samples <- assignment_model$run(0L, 1L)
    assignment_fitted <- assignment_samples$train

    ## --- alpha, beta | theta via WALNUTS, one sampler per bank ------------
    for (k in seq_len(n_traits)) {
      ab <- if (sampling_phase) {
        irt_draw(ab_sampler[[k]], theta[, k])
      } else if (i_sample > warmup_start) {
        irt_warmup(ab_sampler[[k]], theta[, k])
      } else {
        NULL # settle: alpha/beta held at init
      }
      if (!is.null(ab)) {
        p_k <- n_items[k]
        alpha[[k]] <- ab[1:p_k]
        beta[[k]] <- ab[(p_k + 1):(2 * p_k)]
        alpha_mat[[k]] <- matrix(alpha[[k]], n_persons, p_k, byrow = TRUE)
        beta_mat[[k]] <- matrix(beta[[k]], n_persons, p_k, byrow = TRUE)
      }
      if (i_sample == n_burnin) {
        irt_freeze(ab_sampler[[k]])
      }
    }

    ## --- store draws + treatment-effect estimand --------------------------
    if (sampling_phase) {
      i_draw <- i_sample - n_burnin
      cur_frame <- trait_frame(theta, theta_names)
      f1 <- response_model$predict(data.frame(
        z = rep(1, n_persons),
        cur_frame
      ))
      f0 <- response_model$predict(data.frame(
        z = rep(0, n_persons),
        cur_frame
      ))
      ate_draws[i_draw] <- mean(f1 - f0)
      for (k in seq_len(n_traits)) {
        alpha_draws[[k]][i_draw, ] <- alpha[[k]]
        beta_draws[[k]][i_draw, ] <- beta[[k]]
        if (keep_theta) theta_draws[[k]][i_draw, ] <- theta[, k]
      }
      sigma_draws[i_draw] <- response_samples$sigma[1L]
    }

    if (verbose && (i_sample %% 100L == 0L || i_sample == n_total)) {
      message(sprintf(
        "scan %d/%d (%s), theta accept %.2f",
        i_sample,
        n_total,
        if (sampling_phase) "sampling" else "burn-in",
        theta_accept[i_sample]
      ))
    }
  }

  # A single trait unwraps to the bare matrices/vectors it has always returned.
  unwrap <- function(x) if (n_traits == 1L) x[[1L]] else x

  list(
    ate = ate_draws,
    alpha = unwrap(alpha_draws),
    beta = unwrap(beta_draws),
    sigma = sigma_draws,
    theta = if (keep_theta) unwrap(theta_draws) else NULL,
    theta_accept = theta_accept,
    theta_sd = theta_sd,
    theta_sd_trace = theta_sd_trace,
    tuning = unwrap(lapply(ab_sampler, irt_tuning))
  )
}
