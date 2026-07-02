#' Causal inference with a latent IRT confounder via BART and WALNUTS
#'
#' Runs a Gibbs sampler for a causal model in which a latent person trait
#' `theta`, measured by a two-parameter logistic (2PL) item-response model,
#' confounds a binary treatment `z` and a continuous outcome `y`. Within each
#' scan:
#'
#' * `theta` gets one per-person random-walk Metropolis move (the `theta_j` are
#'   conditionally independent given the trees and item parameters, so each
#'   person is accepted or rejected on its own). The accepted values are then
#'   installed into both BART models by a single joint per-observation sweep
#'   ([dbarts::updatePredictorPerObservationJointly()]): a person is kept only
#'   if its `theta` leaves every leaf non-empty in every tree of both models.
#' * The two BART surfaces --- `y ~ f(z, theta)` (response) and
#'   `z ~ f(theta)` (probit assignment) --- each take one tree-update step.
#' * The IRT item parameters `(alpha, beta)` are drawn conditional on `theta`
#'   by the WALNUTS sampler ([irt_item_sampler()]): adapting during burn-in,
#'   then frozen for the sampling phase.
#'
#' The treatment-effect estimand `E[f(1, theta) - f(0, theta)]` is accumulated
#' from the response surface each sampling scan.
#'
#' @param responses A `n_persons` x `n_items` matrix (or data frame) of 0/1 item
#'   responses.
#' @param y Numeric outcome, length `n_persons`.
#' @param z Binary (0/1) treatment, length `n_persons`.
#' @param n_burnin Number of burn-in scans. During burn-in the proposal SD is
#'   adapted and empty-leaf `theta` moves are collapsed rather than rejected.
#' @param n_sampling Number of kept scans.
#' @param theta_sd Initial standard deviation of the per-person `theta` proposal;
#'   adapted toward `theta_accept_target` during burn-in and frozen afterward.
#' @param theta_accept_target Target per-person acceptance rate for the `theta`
#'   proposal adaptation.
#' @param warmup_start Burn-in scan after which the WALNUTS item sampler begins
#'   adapting (before this it is held at its initial `(alpha, beta)` so the
#'   tuning is not polluted by the not-yet-settled `theta`). Must be less than
#'   `n_burnin`, otherwise the item sampler would never adapt; if it is not, it
#'   is reduced to `n_burnin %/% 2` with a warning.
#' @param beta_sd Prior standard deviation for the IRT item difficulties.
#' @param step_size Initial WALNUTS leapfrog step size.
#' @param n_theta_cutpoints Number of cut points for the `theta` predictor in
#'   each BART model (interior quantiles of a standard normal).
#' @param n_trees Number of trees in each BART surface.
#' @param n_threads Number of threads for the BART updates.
#' @param seed Integer seed. This calls [set.seed()] internally, so it
#'   overwrites the global RNG state (`.Random.seed`) as a side effect; it also
#'   seeds the WALNUTS sampler's own (independent) RNG. `y`, `z`, and
#'   `responses` may be numeric or logical, but not factors; `y` must not
#'   contain `NA`.
#' @param keep_theta If `TRUE`, also return the per-scan `theta` draws (a
#'   `n_sampling` x `n_persons` matrix). Off by default to keep the result small.
#' @param verbose If `TRUE`, print progress periodically.
#'
#' @return A list of draws and diagnostics:
#'   \describe{
#'     \item{`ate`}{Length-`n_sampling` vector of treatment-effect draws.}
#'     \item{`alpha`, `beta`}{`n_sampling` x `n_items` matrices of item-parameter
#'       draws.}
#'     \item{`sigma`}{Length-`n_sampling` vector of response residual SD draws.}
#'     \item{`theta`}{`n_sampling` x `n_persons` matrix of `theta` draws, if
#'       `keep_theta = TRUE`.}
#'     \item{`theta_accept`}{Per-scan `theta` acceptance rate (length
#'       `n_burnin + n_sampling`).}
#'     \item{`theta_sd`}{Final (frozen) proposal SD.}
#'     \item{`theta_sd_trace`}{Per-scan adapted proposal SD.}
#'     \item{`tuning`}{Final WALNUTS tuning (see [irt_tuning()]).}
#'   }
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 200, n_items = 30, ate = -0.2, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 100, n_sampling = 200, seed = 1)
#' mean(fit$ate)                      # posterior mean treatment effect (~ -0.2)
#' quantile(fit$ate, c(0.025, 0.975)) # 95% interval
#'
#' @seealso [simulate_irt_causal()], [irt_item_sampler()]
#' @importFrom stats dnorm pnorm plogis qnorm rnorm rexp
#' @export
irt_causal_bart <- function(responses, y, z,
                            n_burnin = 500L, n_sampling = 1000L,
                            theta_sd = 0.6, theta_accept_target = 0.44,
                            warmup_start = 50L,
                            beta_sd = 10, step_size = 0.1,
                            n_theta_cutpoints = 100L, n_trees = 75L,
                            n_threads = 1L, seed = 1L,
                            keep_theta = FALSE, verbose = FALSE) {
  responses <- as_response_matrix(responses)
  n_persons <- nrow(responses)
  n_items <- ncol(responses)
  y <- as.double(y)
  z <- as.double(z)
  if (length(y) != n_persons)
    stop("'y' must have length nrow(responses) = ", n_persons)
  if (length(z) != n_persons)
    stop("'z' must have length nrow(responses) = ", n_persons)
  if (anyNA(y))
    stop("'y' must not contain NA")
  if (anyNA(z) || any(z != 0 & z != 1))
    stop("'z' must be a 0/1 treatment indicator")
  n_burnin <- as.integer(n_burnin)
  n_sampling <- as.integer(n_sampling)
  if (n_burnin < 0L || n_sampling <= 0L)
    stop("'n_burnin' must be >= 0 and 'n_sampling' > 0")
  n_theta_cutpoints <- as.integer(n_theta_cutpoints)
  if (is.na(n_theta_cutpoints) || n_theta_cutpoints < 1L)
    stop("'n_theta_cutpoints' must be >= 1")
  if (as.integer(n_trees) < 1L)
    stop("'n_trees' must be >= 1")

  # warmup_start delays WALNUTS adaptation until theta has settled; if it leaves
  # no adapting scans the item sampler would freeze at its initial tuning, so
  # clamp it below n_burnin (a no-op when n_burnin == 0: there is no warm-up).
  warmup_start <- as.integer(warmup_start)
  if (is.na(warmup_start) || warmup_start < 0L)
    stop("'warmup_start' must be >= 0")
  if (n_burnin > 0L && warmup_start >= n_burnin) {
    new_warmup_start <- n_burnin %/% 2L
    warning(sprintf(
      paste0("'warmup_start' (%d) >= 'n_burnin' (%d): the WALNUTS item sampler ",
             "would never adapt; reducing 'warmup_start' to %d."),
      warmup_start, n_burnin, new_warmup_start))
    warmup_start <- new_warmup_start
  }

  set.seed(seed)
  theta <- rnorm(n_persons)

  cutpoints <- qnorm(seq(0, 1, length.out = n_theta_cutpoints + 2L)[
    -c(1L, n_theta_cutpoints + 2L)])
  control <- dbarts::dbartsControl(n.chains = 1L, n.threads = as.integer(n_threads),
                                   n.trees = as.integer(n_trees))

  # --- response surface (BART): y ~ f(z, theta) ---------------------------
  response_model <- dbarts::dbarts(y ~ z + theta, data.frame(y, z, theta),
                                   control = control)
  response_model$setCutPoints(cutpoints, "theta")
  response_model$sampleTreesFromPrior()
  response_samples <- response_model$run(5L, 1L)

  # --- assignment model (BART): z ~ f(theta) ------------------------------
  assignment_model <- dbarts::dbarts(z ~ theta, data.frame(z, theta),
                                     control = control)
  assignment_model$setCutPoints(cutpoints, "theta")
  assignment_model$sampleTreesFromPrior()
  assignment_samples <- assignment_model$run(5L, 1L)

  # --- IRT item parameters: sampled by WALNUTS (alpha, beta | theta) ------
  alpha <- rep(1, n_items)   # discrimination (positive)
  beta <- rnorm(n_items)     # difficulty
  ab_sampler <- irt_item_sampler(responses, alpha, beta,
                                 beta_sd = beta_sd, step_size = step_size,
                                 seed = seed)

  theta_mat <- matrix(theta, n_persons, n_items)
  alpha_mat <- matrix(alpha, n_persons, n_items, byrow = TRUE)
  beta_mat <- matrix(beta, n_persons, n_items, byrow = TRUE)

  response_fitted <- response_samples$train
  assignment_fitted <- assignment_samples$train

  # storage (sampling phase only)
  ate_draws <- numeric(n_sampling)
  alpha_draws <- matrix(NA_real_, n_sampling, n_items)
  beta_draws <- matrix(NA_real_, n_sampling, n_items)
  sigma_draws <- numeric(n_sampling)
  theta_draws <- if (keep_theta) matrix(NA_real_, n_sampling, n_persons) else NULL
  theta_accept <- numeric(n_burnin + n_sampling)
  theta_sd_trace <- numeric(n_burnin + n_sampling)

  if (n_burnin == 0L) irt_freeze(ab_sampler)

  n_total <- n_burnin + n_sampling
  for (i_sample in seq_len(n_total)) {
    sampling_phase <- i_sample > n_burnin

    ## --- theta: per-person random-walk Metropolis -------------------------
    theta_old <- theta
    theta_prop <- theta + rnorm(n_persons, sd = theta_sd)
    theta_prop_mat <- matrix(theta_prop, n_persons, n_items)

    response_prop <- response_model$predict(data.frame(z, theta = theta_prop))
    assignment_prop <- assignment_model$predict(data.frame(theta = theta_prop))

    # per-person IRT log-likelihood (sum over that person's items)
    irt_ll_cur <- rowSums(plogis(
      (2 * responses - 1) * alpha_mat * (theta_mat - beta_mat), log.p = TRUE))
    irt_ll_prop <- rowSums(plogis(
      (2 * responses - 1) * alpha_mat * (theta_prop_mat - beta_mat), log.p = TRUE))

    # per-person log acceptance ratio (vector of length n_persons)
    lr <- (irt_ll_prop +
      dnorm(y, response_prop, response_samples$sigma, log = TRUE) +
      pnorm((2 * z - 1) * assignment_prop, log.p = TRUE) +
      dnorm(theta_prop, log = TRUE)) -
      (irt_ll_cur +
        dnorm(y, response_fitted, response_samples$sigma, log = TRUE) +
        pnorm((2 * z - 1) * assignment_fitted, log.p = TRUE) +
        dnorm(theta, log = TRUE))

    accept <- !is.na(lr) & (-rexp(n_persons)) <= lr
    theta[accept] <- theta_prop[accept]
    theta_mat <- matrix(theta, n_persons, n_items)
    p_accept <- mean(accept)
    theta_accept[i_sample] <- p_accept

    ## --- install theta into both BART models ------------------------------
    if (sampling_phase) {
      # A theta that empties a leaf has zero tree prior, so the move must be
      # rejected. The joint sweep installs each person's theta in BOTH models or
      # neither -- a systematic-scan Metropolis-within-Gibbs sweep that keeps the
      # two models' 'theta' columns identical and empty-leaf-free.
      installed <- dbarts::updatePredictorPerObservationJointly(
        list(response_model, assignment_model), theta, "theta")
      theta[!installed] <- theta_old[!installed]
      theta_mat <- matrix(theta, n_persons, n_items)
      theta_accept[i_sample] <- mean(theta != theta_old)
    } else {
      # Burn-in: collapse empty leaves so theta and the trees stay in sync while
      # the chain/adaptation move. These draws are discarded.
      response_model$setPredictor(theta, "theta", forceUpdate = TRUE)
      assignment_model$setPredictor(theta, "theta", forceUpdate = TRUE)
    }

    ## --- adapt theta proposal SD toward target (burn-in only) -------------
    # Robbins-Monro on log(theta_sd) with decreasing gain; frozen once sampling
    # begins, so kept draws come from a fixed kernel.
    if (!sampling_phase) {
      gamma_t <- i_sample^(-0.6)
      theta_sd <- exp(log(theta_sd) + gamma_t * (p_accept - theta_accept_target))
      theta_sd <- min(max(theta_sd, 1e-3), 10)
    }
    theta_sd_trace[i_sample] <- theta_sd

    ## --- BART trees | theta (one Gibbs step each) -------------------------
    response_samples <- response_model$run(0L, 1L)
    response_fitted <- response_samples$train
    assignment_samples <- assignment_model$run(0L, 1L)
    assignment_fitted <- assignment_samples$train

    ## --- alpha, beta | theta via WALNUTS ----------------------------------
    ab <- if (sampling_phase) {
      irt_draw(ab_sampler, theta)
    } else if (i_sample > warmup_start) {
      irt_warmup(ab_sampler, theta)
    } else {
      NULL   # settle: alpha/beta held at init
    }
    if (!is.null(ab)) {
      alpha <- ab[1:n_items]
      beta <- ab[(n_items + 1):(2 * n_items)]
      alpha_mat <- matrix(alpha, n_persons, n_items, byrow = TRUE)
      beta_mat <- matrix(beta, n_persons, n_items, byrow = TRUE)
    }
    if (i_sample == n_burnin) irt_freeze(ab_sampler)

    ## --- store draws + treatment-effect estimand --------------------------
    if (sampling_phase) {
      k <- i_sample - n_burnin
      f1 <- response_model$predict(data.frame(z = rep(1, n_persons), theta = theta))
      f0 <- response_model$predict(data.frame(z = rep(0, n_persons), theta = theta))
      ate_draws[k] <- mean(f1 - f0)
      alpha_draws[k, ] <- alpha
      beta_draws[k, ] <- beta
      sigma_draws[k] <- response_samples$sigma[1L]
      if (keep_theta) theta_draws[k, ] <- theta
    }

    if (verbose && (i_sample %% 100L == 0L || i_sample == n_total))
      message(sprintf("scan %d/%d (%s), theta accept %.2f",
                      i_sample, n_total,
                      if (sampling_phase) "sampling" else "burn-in",
                      theta_accept[i_sample]))
  }

  list(
    ate = ate_draws,
    alpha = alpha_draws,
    beta = beta_draws,
    sigma = sigma_draws,
    theta = theta_draws,
    theta_accept = theta_accept,
    theta_sd = theta_sd,
    theta_sd_trace = theta_sd_trace,
    tuning = irt_tuning(ab_sampler)
  )
}
