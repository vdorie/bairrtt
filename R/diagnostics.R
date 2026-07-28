# Extraction and convergence diagnostics over an "irt_causal_fit". Follows
# dbarts's own idiom (its R/diagnostics.R and R/hooks.R): summary() delegates
# to posterior when it is installed and degrades to quantiles when it is not,
# and posterior's as_draws_* generics are registered dynamically so posterior
# stays in Suggests.

# --- draws in posterior's (iteration, chain, variable) layout ----------------

# Flatten the requested quantities into one named array. Per-bank lists become
# variables named alpha[bank1,item3], so heterogeneous banks -- which cannot
# share a rectangular axis -- still land in a single array.
irt_draws_array <- function(object, vars) {
  n_traits <- object$n_traits
  pieces <- list()
  names_out <- character(0L)
  for (nm in vars) {
    value <- object[[nm]]
    if (is.null(value)) {
      next
    }
    banks <- if (is.list(value)) value else list(value)
    for (k in seq_along(banks)) {
      arr <- banks[[k]]
      if (length(dim(arr)) == 2L) {
        # a per-draw scalar: one variable
        pieces[[length(pieces) + 1L]] <- array(arr, c(dim(arr), 1L))
        names_out <- c(names_out, nm)
      } else {
        pieces[[length(pieces) + 1L]] <- arr
        label <- if (is.list(value) && n_traits > 1L) {
          sprintf("%s[bank%d,%d]", nm, k, seq_len(dim(arr)[3L]))
        } else {
          sprintf("%s[%d]", nm, seq_len(dim(arr)[3L]))
        }
        names_out <- c(names_out, label)
      }
    }
  }
  if (length(pieces) == 0L) {
    return(NULL)
  }
  out <- array(
    unlist(pieces, use.names = FALSE),
    dim = c(dim(pieces[[1L]])[1L:2L], length(names_out))
  )
  dimnames(out) <- list(NULL, NULL, names_out)
  out
}

as_draws_array.irt_causal_fit <- function(x, vars = irt_default_vars(x), ...) {
  posterior::as_draws_array(irt_draws_array(x, vars))
}

as_draws_df.irt_causal_fit <- function(x, vars = irt_default_vars(x), ...) {
  posterior::as_draws_df(irt_draws_array(x, vars))
}

irt_default_vars <- function(x) c("ate", "sigma")

posterior_available <- function() {
  requireNamespace("posterior", quietly = TRUE)
}

# --- extract ----------------------------------------------------------------

# dbarts owns the generic and stan4bart registers on it too. bairrtt Imports
# dbarts rather than Depending on it, so re-export the generic: without this
# `extract(fit)` only works for users who have separately attached dbarts.
#' @importFrom dbarts extract
#' @export
dbarts::extract


#' Extract draws from a fitted IRT-confounded causal model
#'
#' Method for \pkg{dbarts}'s `extract()` generic, so a `bairrtt` fit is read the
#' same way a \pkg{dbarts} or \pkg{stan4bart} fit is.
#'
#' @param object An `"irt_causal_fit"` from [irt_causal_bart()].
#' @param type Which quantity to extract.
#' @param combine_chains If `TRUE` (the default) the chain axis is collapsed, so
#'   a per-draw scalar comes back as a length-`n_sampling * n_chains` vector and
#'   a per-draw vector as a `(n_sampling * n_chains)` x `n_par` matrix. If
#'   `FALSE` the chain axis is kept.
#' @param ... Ignored.
#'
#' @return A vector, matrix, or array as described under `combine_chains`; with
#'   more than one trait, `alpha`, `beta`, `gamma`, and `theta` come back as a
#'   list with one element per item bank.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 100, n_items = 15, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 20, n_sampling = 40,
#'                        warmup_start = 10, n_chains = 2, seed = 1)
#' length(extract(fit, "ate"))                       # 40 x 2 draws, pooled
#' dim(extract(fit, "ate", combine_chains = FALSE))  # 40 x 2
#'
#' @seealso [irt_causal_bart()], [summary.irt_causal_fit()]
#' @importFrom dbarts extract
#' @method extract irt_causal_fit
#' @export
extract.irt_causal_fit <- function(
  object,
  type = c(
    "ate",
    "sigma",
    "alpha",
    "beta",
    "gamma",
    "theta",
    "theta_accept",
    "theta_sd_trace"
  ),
  combine_chains = TRUE,
  ...
) {
  type <- match.arg(type)
  value <- object[[type]]
  if (is.null(value)) {
    stop(
      "'",
      type,
      "' is not present in this fit",
      if (type == "theta") "; refit with keep_theta = TRUE" else "",
      if (type == "gamma") "; it needs covariates, so pass 'x'" else ""
    )
  }
  if (!combine_chains) {
    return(value)
  }
  collapse <- function(arr) {
    if (length(dim(arr)) == 2L) {
      return(as.vector(arr))
    }
    # draws x chains x par -> (draws * chains) x par
    matrix(aperm(arr, c(1L, 2L, 3L)), prod(dim(arr)[1L:2L]), dim(arr)[3L])
  }
  if (is.list(value)) lapply(value, collapse) else collapse(value)
}

# --- summary ----------------------------------------------------------------

#' Convergence summary for a fitted IRT-confounded causal model
#'
#' Posterior summaries with convergence diagnostics. When \pkg{posterior} is
#' installed this is [posterior::summarise_draws()], so it reports
#' rank-normalized split R-hat alongside bulk and tail effective sample size;
#' otherwise it degrades to means, SDs, and quantiles.
#'
#' `vars` defaults to `"ate"` and `"sigma"` deliberately. R-hat is conventionally
#' read against a 1.01 threshold, but the maximum of many R-hats has a null
#' distribution that is not centred at one: with a few hundred effective draws
#' per parameter, a *perfectly converged* fit exceeds 1.01 on at least one of a
#' few hundred item and person parameters with probability near one. Judge
#' convergence on the estimand and the residual scale, and read the item and
#' person parameters --- available via `vars` --- as informational, preferring
#' their effective sample sizes to their R-hats.
#'
#' @param object An `"irt_causal_fit"` from [irt_causal_bart()].
#' @param vars Which quantities to summarize: any of `"ate"`, `"sigma"`,
#'   `"alpha"`, `"beta"`, `"gamma"`, and `"theta"`.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.irt_causal_fit"`, with the per-variable
#'   statistics in `$stats`.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 100, n_items = 15, seed = 1)
#' fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
#'                        n_burnin = 20, n_sampling = 40,
#'                        warmup_start = 10, n_chains = 2, seed = 1)
#' summary(fit)
#' summary(fit, vars = c("ate", "alpha"))
#'
#' @seealso [irt_causal_bart()], [extract.irt_causal_fit()]
#' @method summary irt_causal_fit
#' @export
summary.irt_causal_fit <- function(object, vars = c("ate", "sigma"), ...) {
  arr <- irt_draws_array(object, vars)
  if (is.null(arr)) {
    stop("none of 'vars' is present in this fit")
  }
  have_posterior <- posterior_available()
  stats <- if (have_posterior) {
    posterior::summarise_draws(posterior::as_draws_array(arr))
  } else {
    quantile_summary(arr)
  }
  structure(
    list(
      call = object$call,
      stats = stats,
      posterior = have_posterior,
      n_chains = object$n_chains,
      n_sampling = object$n_sampling,
      n_traits = object$n_traits
    ),
    class = "summary.irt_causal_fit"
  )
}

# plain per-variable mean/sd/quantiles, pooling draws and chains; the degrade
# path when posterior is not installed
quantile_summary <- function(arr) {
  pooled <- matrix(arr, prod(dim(arr)[1L:2L]), dim(arr)[3L])
  out <- data.frame(
    variable = dimnames(arr)[[3L]],
    mean = colMeans(pooled),
    sd = apply(pooled, 2L, stats::sd),
    q2.5 = apply(pooled, 2L, stats::quantile, probs = 0.025, names = FALSE),
    median = apply(pooled, 2L, stats::quantile, probs = 0.5, names = FALSE),
    q97.5 = apply(pooled, 2L, stats::quantile, probs = 0.975, names = FALSE),
    row.names = NULL
  )
  # R-hat needs more than one chain to mean anything
  if (dim(arr)[2L] > 1L) {
    out$rhat <- vapply(
      asplit(arr, 3L),
      rhat_matrix,
      numeric(1L),
      USE.NAMES = FALSE
    )
  }
  out
}

#' @export
print.summary.irt_causal_fit <- function(x, ...) {
  cat(sprintf(
    "irt_causal_bart fit: %d chain%s x %d draws, %d trait%s\n",
    x$n_chains,
    if (x$n_chains == 1L) "" else "s",
    x$n_sampling,
    x$n_traits,
    if (x$n_traits == 1L) "" else "s"
  ))
  if (!is.null(x$call)) {
    cat("call: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  }
  print(as.data.frame(x$stats), row.names = FALSE)
  if (!x$posterior) {
    cat("\ninstall 'posterior' for R-hat and effective sample size\n")
  } else if (x$n_chains == 1L) {
    cat(
      "\nwith one chain R-hat is within-chain (split-half) only; refit with ",
      "n_chains > 1 for a between-chain diagnostic\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @export
print.irt_causal_fit <- function(x, ...) {
  ate <- as.vector(x$ate)
  ci <- stats::quantile(ate, c(0.025, 0.975), names = FALSE)
  cat(sprintf(
    "<irt_causal_fit: %d chain%s x %d draws, %d trait%s>\n",
    x$n_chains,
    if (x$n_chains == 1L) "" else "s",
    x$n_sampling,
    x$n_traits,
    if (x$n_traits == 1L) "" else "s"
  ))
  cat(sprintf("  ate  %.4f  95%% [%.4f, %.4f]\n", mean(ate), ci[1L], ci[2L]))
  cat(sprintf("  seeds %s\n", paste(x$seeds, collapse = ", ")))
  if (x$n_chains > 1L) {
    cat("  summary() for R-hat and effective sample size\n")
  } else {
    cat("  n_chains > 1 for between-chain convergence diagnostics\n")
  }
  invisible(x)
}

# --- rank-normalized split R-hat (the no-posterior degrade path) -------------

# Kept because posterior is Suggests-only. Matches posterior::rhat: the larger
# of the plain and folded rank-normalized split statistics (Vehtari et al.
# 2021). `m` is draws x chains.
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
#
# Zero within-chain variance is two opposite situations and they must not be
# conflated. If the chains also agree, the quantity is constant and R-hat is
# undefined (NA). If they sit at DIFFERENT values, every chain is frozen
# somewhere else -- the loudest non-convergence there is -- and R-hat is
# infinite. Returning NA for the second hides it.
rhat_basic <- function(m) {
  n <- nrow(m)
  within <- mean(apply(m, 2L, stats::var))
  between <- n * stats::var(colMeans(m))
  if (!is.finite(within) || !is.finite(between)) {
    return(NA_real_)
  }
  if (within <= 0) {
    return(if (between > 0) Inf else NA_real_)
  }
  sqrt((between / within + n - 1) / n)
}
