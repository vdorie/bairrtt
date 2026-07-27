# Multiple chains: the R-hat statistic itself, then the surface over a fit.
# The statistic is tested on synthetic input rather than on draws -- it has to
# be right independently of whether any particular fit converged.

rhat_matrix <- getFromNamespace("rhat_matrix", "bairrtt")

# --- the statistic -----------------------------------------------------------
set.seed(1)
iid <- matrix(rnorm(400 * 4), 400L, 4L)
expect_true(rhat_matrix(iid) < 1.05) # chains agree

shifted <- iid
shifted[, 1L] <- shifted[, 1L] + 10
expect_true(rhat_matrix(shifted) > 1.5) # one chain elsewhere

# the case this whole item exists for: half the chains on -theta, half on theta
flipped <- iid
flipped[, 1:2] <- -abs(flipped[, 1:2])
flipped[, 3:4] <- abs(flipped[, 3:4])
expect_true(rhat_matrix(flipped) > 1.5)

# chains agreeing on location but not scale: caught only by the folded half of
# the statistic, so this pins that half down
scale_only <- cbind(
  matrix(rnorm(400 * 2), 400L, 2L),
  matrix(rnorm(400 * 2, 0, 5), 400L, 2L)
)
expect_true(rhat_matrix(scale_only) > 1.1)

expect_true(is.na(rhat_matrix(matrix(3, 100L, 4L)))) # constant: undefined
expect_true(is.na(rhat_matrix(matrix(rnorm(6), 3L, 2L)))) # too few draws
expect_true(is.na(rhat_matrix(matrix(c(NA, rnorm(39)), 10L, 4L)))) # NA

# agreement with the reference implementation, to the last bit
if (requireNamespace("posterior", quietly = TRUE)) {
  set.seed(2)
  for (i in 1:5) {
    m <- matrix(rt(200 * 4, df = 3), 200L, 4L)
    expect_equal(rhat_matrix(m), posterior::rhat(m), tolerance = 1e-10)
  }
}

# --- n_chains = 1 is unchanged ----------------------------------------------
sim <- simulate_irt_causal(n_persons = 150L, n_items = 12L, seed = 11)
common <- list(
  sim$responses,
  sim$y,
  sim$z,
  n_burnin = 30L,
  n_sampling = 60L,
  warmup_start = 10L,
  seed = 11L
)

one <- do.call(irt_causal_bart, common)
expect_false(inherits(one, "irt_causal_chains"))
expect_equal(length(one$ate), 60L)
expect_error(irt_rhat(one)) # not a chains object
expect_error(irt_chain_draws(one, "ate"))

# --- multiple chains ---------------------------------------------------------
fit <- do.call(irt_causal_bart, c(common, list(n_chains = 3L)))
expect_inherits(fit, "irt_causal_chains")
expect_equal(fit$n_chains, 3L)
expect_equal(fit$seeds, c(11L, 12L, 13L)) # stride is n_traits = 1
expect_equal(length(fit$chains), 3L)

# chain 1 takes the base seed, so it must reproduce the single-chain fit
expect_identical(fit$chains[[1L]], one)

# sequential and forked runs are the same draws
fit_par <- do.call(
  irt_causal_bart,
  c(common, list(n_chains = 3L, n_cores = 2L))
)
expect_identical(fit_par$chains, fit$chains)

# --- draw extraction ---------------------------------------------------------
ate <- irt_chain_draws(fit, "ate")
expect_equal(dim(ate), c(60L, 3L))
expect_identical(ate[, 2L], fit$chains[[2L]]$ate)
expect_equal(mean(ate), mean(unlist(lapply(fit$chains, `[[`, "ate"))))

alpha <- irt_chain_draws(fit, "alpha")
expect_equal(dim(alpha), c(60L, 3L, 12L)) # draws x chains x items
expect_identical(alpha[, 3L, 5L], fit$chains[[3L]]$alpha[, 5L])

expect_equal(dim(irt_chain_draws(fit, "sigma")), c(60L, 3L))
expect_error(irt_chain_draws(fit, "theta")) # keep_theta was FALSE
expect_error(irt_chain_draws(fit, "nonesuch")) # not a known quantity

# --- the diagnostic over a fit ----------------------------------------------
rh <- irt_rhat(fit)
expect_equal(length(rh$ate), 1L)
expect_equal(length(rh$sigma), 1L)
expect_equal(length(rh$alpha), 12L) # one per item
expect_equal(length(rh$beta), 12L)
expect_null(rh$theta) # not kept, so not reported
expect_true(is.finite(rh$max))
expect_equal(
  rh$max,
  max(c(rh$ate, rh$sigma, rh$alpha, rh$beta), na.rm = TRUE)
)
expect_true(all(rh$alpha > 0.9)) # R-hat is bounded near 1 from below

expect_stdout(print(fit), "irt_causal_chains")

# --- two traits: everything per-bank ----------------------------------------
sim2 <- simulate_irt_causal(
  n_persons = 120L,
  n_items = c(10L, 14L),
  n_traits = 2L,
  seed = 12
)
fit2 <- irt_causal_bart(
  sim2$responses,
  sim2$y,
  sim2$z,
  n_burnin = 20L,
  n_sampling = 40L,
  warmup_start = 8L,
  n_chains = 2L,
  seed = 12L,
  keep_theta = TRUE
)
expect_equal(fit2$seeds, c(12L, 14L)) # stride is n_traits = 2, so no overlap
expect_equal(fit2$n_traits, 2L)

alpha2 <- irt_chain_draws(fit2, "alpha")
expect_equal(length(alpha2), 2L)
expect_equal(dim(alpha2[[1L]]), c(40L, 2L, 10L))
expect_equal(dim(alpha2[[2L]]), c(40L, 2L, 14L))

rh2 <- irt_rhat(fit2)
expect_equal(lengths(rh2$alpha), c(10L, 14L))
expect_equal(lengths(rh2$theta), c(120L, 120L))
expect_true(is.finite(rh2$max))

# --- argument checks ---------------------------------------------------------
expect_error(do.call(irt_causal_bart, c(common, list(n_chains = 0L))))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 2L, n_cores = 0L))
))
# explicit seeds must be distinct and the right length
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 2L, seeds = c(1L, 1L)))
))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 3L, seeds = c(1L, 2L)))
))
# explicit seeds override `seed`
fit_s <- do.call(
  irt_causal_bart,
  c(common, list(n_chains = 2L, seeds = c(101L, 202L)))
)
expect_equal(fit_s$seeds, c(101L, 202L))
