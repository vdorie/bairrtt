# Multiple chains: the runner behind irt_causal_bart(n_chains > 1), plus the
# accessor and the between-chain diagnostic over the result.

# Run one chain per seed. Chains cannot share a process: each owns mutable
# dbarts samplers and an XPtr-held WALNUTS state (see the IrtSampler note in
# inst/include/bairrtt_types.h), so parallelism here is fork-based, not
# threaded, and is simply unavailable on Windows. No RNG brokering is needed --
# each chain calls set.seed() on its own seed, so sequential and parallel runs
# give identical draws.
run_chains <- function(spec, seeds, n_cores) {
  n_chains <- length(seeds)
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("'n_cores' must be a single integer >= 1")
  }
  n_cores <- min(n_cores, n_chains)
  forkable <- n_cores > 1L && .Platform$OS.type != "windows"
  if (n_cores > 1L && !forkable) {
    message(
      "'n_cores' > 1 needs fork-based parallelism, which this platform does ",
      "not provide; running ",
      n_chains,
      " chains sequentially."
    )
  }

  one <- function(i) {
    do.call(
      irt_causal_bart_chain,
      c(spec, list(seed = seeds[i], chain_id = i))
    )
  }
  fits <- if (forkable) {
    parallel::mclapply(seq_len(n_chains), one, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_chains), one)
  }

  # mclapply() returns a worker error as a try-error element rather than
  # raising it, so a silent partial result is the default. Refuse that.
  failed <- which(vapply(fits, inherits, logical(1L), what = "try-error"))
  if (length(failed) > 0L) {
    stop(
      "chain(s) ",
      paste(failed, collapse = ", "),
      " failed: ",
      conditionMessage(attr(fits[[failed[1L]]], "condition"))
    )
  }
  fits
}

#' Per-chain draws from a multi-chain fit
#'
#' Pulls one quantity out of every chain of an `"irt_causal_chains"` object and
#' stacks the chains into a single array, which is the shape both [irt_rhat()]
#' and any pooling want. Note the plural: [irt_draw()] is the unrelated
#' single-step advance of a standalone item sampler.
#'
#' @param x An `"irt_causal_chains"` object, from [irt_causal_bart()] with
#'   `n_chains > 1`.
#' @param what Which quantity to extract.
#'
#' @return For a per-draw scalar (`ate`, `sigma`), an `n_sampling` x `n_chains`
#'   matrix. For a per-draw vector (`alpha`, `beta`, `theta`), an `n_sampling` x
#'   `n_chains` x `n_par` array. With more than one trait, a list of one such
#'   array per item bank. Pool with `c()` or `as.vector()`.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 100, n_items = 15, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 20, n_sampling = 40,
#'                        warmup_start = 10, n_chains = 2, seed = 1)
#' dim(irt_chain_draws(fit, "ate"))   # 40 draws x 2 chains
#' mean(irt_chain_draws(fit, "ate"))  # pooled posterior mean
#'
#' @seealso [irt_rhat()], [irt_causal_bart()]
#' @export
irt_chain_draws <- function(
  x,
  what = c("ate", "sigma", "alpha", "beta", "theta")
) {
  check_chains(x)
  what <- match.arg(what)
  per_chain <- lapply(x$chains, `[[`, what)
  if (is.null(per_chain[[1L]])) {
    stop(
      "'",
      what,
      "' was not kept by this fit",
      if (what == "theta") "; refit with keep_theta = TRUE" else ""
    )
  }
  bind_chain_draws(per_chain)
}

# Stack a per-chain list into a chain-major array, recursing into the per-bank
# lists that more than one trait produces.
bind_chain_draws <- function(per_chain) {
  first <- per_chain[[1L]]
  if (is.list(first)) {
    return(lapply(seq_along(first), function(k) {
      bind_chain_draws(lapply(per_chain, `[[`, k))
    }))
  }
  if (is.matrix(first)) {
    # each chain is draws x parameters; unlist runs draw-fastest, so build
    # draws x parameters x chains and then swap the last two axes
    stacked <- array(
      unlist(per_chain, use.names = FALSE),
      dim = c(nrow(first), ncol(first), length(per_chain))
    )
    return(aperm(stacked, c(1L, 3L, 2L)))
  }
  matrix(unlist(per_chain, use.names = FALSE), ncol = length(per_chain))
}

#' Between-chain convergence diagnostic (R-hat)
#'
#' Rank-normalized split R-hat for every quantity a multi-chain fit kept. This
#' is the diagnostic a single chain cannot provide, and the one that exposes the
#' sign and label multimodality that IRT models invite: a chain that settles on
#' `-theta` rather than `theta` disagrees with its siblings, and R-hat says so.
#'
#' Each chain is split in half before comparison, so a within-chain trend
#' registers as between-chain disagreement. The reported value is the larger of
#' the plain and folded rank-normalized statistics, following Vehtari et al.
#' (2021), which is what `posterior::rhat()` computes. Values above roughly 1.01
#' indicate the chains have not mixed; a constant quantity gives `NA`, where
#' R-hat is undefined.
#'
#' @param x An `"irt_causal_chains"` object, from [irt_causal_bart()] with
#'   `n_chains > 1`.
#'
#' @return A named list with one entry per quantity (`ate` and `sigma` scalar,
#'   `alpha`/`beta`/`theta` one value per parameter, or a list of those per item
#'   bank when there is more than one trait), plus `max`, the largest R-hat over
#'   all of them.
#'
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and
#'   Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#'   improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*
#'   16(2), 667--718.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 100, n_items = 15, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 20, n_sampling = 40,
#'                        warmup_start = 10, n_chains = 2, seed = 1)
#' irt_rhat(fit)$ate
#' irt_rhat(fit)$max
#'
#' @seealso [irt_chain_draws()], [irt_causal_bart()]
#' @export
irt_rhat <- function(x) {
  check_chains(x)
  out <- list()
  for (nm in c("ate", "sigma", "alpha", "beta", "theta")) {
    if (is.null(x$chains[[1L]][[nm]])) {
      next
    }
    out[[nm]] <- rhat_of(irt_chain_draws(x, nm))
  }
  all_values <- unlist(out, use.names = FALSE)
  out$max <- if (all(is.na(all_values))) {
    NA_real_
  } else {
    max(all_values, na.rm = TRUE)
  }
  out
}

# Walk whatever irt_chain_draws() returned down to draws x chains matrices.
rhat_of <- function(draws) {
  if (is.list(draws)) {
    return(lapply(draws, rhat_of))
  }
  if (length(dim(draws)) == 3L) {
    # one draws x chains slice per parameter
    return(vapply(
      asplit(draws, 3L),
      rhat_matrix,
      numeric(1L),
      USE.NAMES = FALSE
    ))
  }
  rhat_matrix(draws)
}

# Rank-normalized split R-hat, as the larger of the plain and folded values.
# `m` is draws x chains.
rhat_matrix <- function(m) {
  m <- as.matrix(m)
  if (nrow(m) < 4L || anyNA(m) || !all(is.finite(m))) {
    return(NA_real_)
  }
  if (all(m == m[1L])) {
    return(NA_real_)
  }
  folded <- abs(m - stats::median(m))
  max(rhat_rank(m), rhat_rank(folded))
}

rhat_rank <- function(m) {
  rhat_basic(z_scale(split_chains(m)))
}

# Split every chain in half, so a within-chain trend shows up as between-chain
# disagreement. An odd number of draws drops the middle one.
split_chains <- function(m) {
  n <- nrow(m)
  if (n < 2L) {
    return(m)
  }
  half <- n / 2
  cbind(
    m[seq_len(floor(half)), , drop = FALSE],
    m[ceiling(half + 1):n, , drop = FALSE]
  )
}

# Blom rank normalization: rank across all chains at once, then to normal
# scores. This is what makes the statistic robust to heavy tails and to chains
# that agree on location but not on scale.
z_scale <- function(m) {
  r <- rank(m, ties.method = "average")
  matrix(stats::qnorm((r - 3 / 8) / (length(r) + 1 / 4)), nrow = nrow(m))
}

# The classic Gelman-Rubin between/within statistic, on already-transformed
# draws.
rhat_basic <- function(m) {
  n <- nrow(m)
  within <- mean(apply(m, 2L, stats::var))
  if (!is.finite(within) || within <= 0) {
    return(NA_real_)
  }
  between <- n * stats::var(colMeans(m))
  sqrt((between / within + n - 1) / n)
}

#' @export
print.irt_causal_chains <- function(x, ...) {
  ate <- irt_chain_draws(x, "ate")
  ci <- stats::quantile(ate, c(0.025, 0.975), names = FALSE)
  max_rhat <- irt_rhat(x)$max
  cat(sprintf(
    "<irt_causal_chains: %d chains x %d draws, %d trait%s>\n",
    x$n_chains,
    nrow(ate),
    x$n_traits,
    if (x$n_traits == 1L) "" else "s"
  ))
  cat(sprintf(
    "  ate       %.4f  95%% [%.4f, %.4f]\n",
    mean(ate),
    ci[1L],
    ci[2L]
  ))
  cat(sprintf(
    "  max R-hat %.4f%s\n",
    max_rhat,
    if (isTRUE(max_rhat > 1.01)) "   not converged (> 1.01)" else ""
  ))
  cat(sprintf("  seeds     %s\n", paste(x$seeds, collapse = ", ")))
  invisible(x)
}

check_chains <- function(x) {
  if (!inherits(x, "irt_causal_chains")) {
    stop(
      "'x' must be an 'irt_causal_chains' (see irt_causal_bart() with ",
      "n_chains > 1)"
    )
  }
}
