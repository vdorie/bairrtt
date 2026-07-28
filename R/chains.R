# Running chains and assembling them into a fit. The diagnostics over the
# result live in R/diagnostics.R.

# Run one chain per seed. Chains cannot share a process: each owns mutable
# dbarts samplers and an XPtr-held WALNUTS state (see the IrtSampler note in
# inst/include/bairrtt_types.h), so parallelism here is fork-based, not
# threaded, and is simply unavailable on Windows. No RNG brokering is needed --
# each chain calls set.seed() on its own seed, so sequential and parallel runs
# give identical draws.
run_chains <- function(spec, seeds, n_cores) {
  n_chains <- length(seeds)
  n_cores <- min(n_cores, n_chains) # validated in irt_causal_bart()
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
    # mc.preschedule = FALSE gives each chain its own process, so a worker that
    # dies takes one chain down rather than its whole prescheduled batch. With
    # a handful of long-running chains the extra forks cost nothing.
    parallel::mclapply(
      seq_len(n_chains),
      one,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(n_chains), one)
  }

  # mclapply() never raises: a worker that ERRORS comes back as a try-error,
  # and a worker that DIES (segfault in dbarts or WALNUTS, OOM kill) comes back
  # as a plain NULL with only a warning. Both are silent partial results, and
  # NULL is the more dangerous one -- it recycles through array() downstream
  # without even a warning. Refuse both.
  died <- vapply(fits, is.null, logical(1L))
  errored <- vapply(fits, inherits, logical(1L), what = "try-error")
  failed <- which(died | errored)
  if (length(failed) > 0L) {
    detail <- if (errored[failed[1L]]) {
      conditionMessage(attr(fits[[failed[1L]]], "condition"))
    } else {
      "the worker process died without returning a result"
    }
    stop(
      "chain(s) ",
      paste(failed, collapse = ", "),
      " failed: ",
      detail
    )
  }
  fits
}

# Stack per-chain results into one object. The chain axis is always present,
# including at a single chain, so nothing downstream has to branch on chain
# count -- the cost is that `ate` is an n_sampling x n_chains matrix rather
# than a vector, which mean(), quantile(), and hist() all handle unchanged.
new_irt_causal_fit <- function(chains, seeds, n_traits, n_items, call) {
  stack <- function(nm) {
    if (is.null(chains[[1L]][[nm]])) {
      return(NULL)
    }
    bind_chain_draws(lapply(chains, `[[`, nm))
  }
  structure(
    list(
      ate = stack("ate"),
      alpha = stack("alpha"),
      beta = stack("beta"),
      sigma = stack("sigma"),
      theta = stack("theta"),
      theta_accept = stack("theta_accept"),
      theta_sd = vapply(chains, `[[`, numeric(1L), "theta_sd"),
      theta_sd_trace = stack("theta_sd_trace"),
      tuning = lapply(chains, `[[`, "tuning"),
      n_chains = length(chains),
      n_traits = n_traits,
      n_items = n_items,
      n_sampling = length(chains[[1L]]$ate),
      seeds = seeds,
      call = call
    ),
    class = "irt_causal_fit"
  )
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
