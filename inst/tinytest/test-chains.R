# Multiple chains: the R-hat statistic itself, then the surface over a fit.
# The statistic is tested on synthetic input rather than on draws -- it has to
# be right independently of whether any particular fit converged. It is the
# degrade path used when posterior is not installed.

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

# Every chain frozen at a DIFFERENT value is the loudest non-convergence there
# is, and must not be conflated with the constant case above. R-hat is infinite,
# not undefined. (At two chains the folded half of the statistic collapses to a
# constant and returns NA, which is what posterior::rhat does too.)
for (k in 3:5) {
  stuck <- matrix(rep(seq_len(k), each = 20L), 20L, k)
  expect_true(is.infinite(rhat_matrix(stuck)))
  if (requireNamespace("posterior", quietly = TRUE)) {
    expect_equal(rhat_matrix(stuck), posterior::rhat(stuck))
  }
}

# agreement with the reference implementation, to the last bit
if (requireNamespace("posterior", quietly = TRUE)) {
  set.seed(2)
  for (i in 1:5) {
    m <- matrix(rt(200 * 4, df = 3), 200L, 4L)
    expect_equal(rhat_matrix(m), posterior::rhat(m), tolerance = 1e-10)
  }
}

# --- the chain axis is always present ----------------------------------------
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

one <- do.call(irt_causal_bart, c(common, list(n_chains = 1L)))
expect_inherits(one, "irt_causal_fit")
expect_equal(dim(one$ate), c(60L, 1L))
expect_equal(one$n_chains, 1L)

fit <- do.call(irt_causal_bart, c(common, list(n_chains = 3L)))
expect_inherits(fit, "irt_causal_fit")
expect_equal(fit$n_chains, 3L)
expect_equal(fit$seeds, c(11L, 12L, 13L)) # stride is n_traits = 1
expect_equal(dim(fit$ate), c(60L, 3L))
expect_equal(dim(fit$alpha), c(60L, 3L, 12L))
expect_equal(length(fit$theta_sd), 3L)
expect_equal(length(fit$tuning), 3L)

# chain 1 takes the base seed, so it must reproduce the single-chain fit
expect_identical(fit$ate[, 1L], one$ate[, 1L])
expect_identical(fit$alpha[, 1L, ], one$alpha[, 1L, ])

# sequential and forked runs are the same draws
fit_par <- do.call(
  irt_causal_bart,
  c(common, list(n_chains = 3L, n_cores = 2L))
)
expect_identical(fit_par$ate, fit$ate)
expect_identical(fit_par$alpha, fit$alpha)

# --- extract -----------------------------------------------------------------
ate <- extract(fit, "ate")
expect_equal(length(ate), 180L) # 60 draws x 3 chains, pooled
expect_equal(ate[1:60], fit$ate[, 1L])
expect_equal(dim(extract(fit, "ate", combine_chains = FALSE)), c(60L, 3L))

alpha <- extract(fit, "alpha")
expect_equal(dim(alpha), c(180L, 12L)) # pooled draws x items
expect_equal(alpha[1:60, 5L], fit$alpha[, 1L, 5L])
expect_equal(
  dim(extract(fit, "alpha", combine_chains = FALSE)),
  c(60L, 3L, 12L)
)

expect_error(extract(fit, "theta")) # keep_theta was FALSE
expect_error(extract(fit, "nonesuch")) # not a known quantity

# --- summary -----------------------------------------------------------------
s <- summary(fit)
expect_inherits(s, "summary.irt_causal_fit")
expect_equal(sort(as.character(s$stats$variable)), c("ate", "sigma"))
expect_equal(s$n_chains, 3L)
expect_stdout(print(s), "ate")
expect_stdout(print(fit), "irt_causal_fit")

# item parameters are available but not the default: the max over hundreds of
# R-hats has a null distribution that is not centred at one, so thresholding it
# would flag converged fits
s_alpha <- summary(fit, vars = c("ate", "alpha"))
expect_equal(nrow(s_alpha$stats), 13L) # ate + 12 items
expect_true("alpha[1]" %in% as.character(s_alpha$stats$variable))
expect_error(summary(fit, vars = "theta")) # not kept

if (requireNamespace("posterior", quietly = TRUE)) {
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(s$stats)))
  # the draws array round-trips through posterior
  da <- posterior::as_draws_array(fit)
  expect_equal(posterior::niterations(da), 60L)
  expect_equal(posterior::nchains(da), 3L)
}

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

alpha2 <- extract(fit2, "alpha")
expect_equal(length(alpha2), 2L)
expect_equal(dim(alpha2[[1L]]), c(80L, 10L)) # 40 draws x 2 chains, pooled
expect_equal(dim(alpha2[[2L]]), c(80L, 14L))

# heterogeneous banks cannot share a rectangular axis, so they are flattened
# into named variables instead
s2 <- summary(fit2, vars = c("ate", "alpha"))
expect_equal(nrow(s2$stats), 25L) # ate + 10 + 14
expect_true("alpha[bank2,14]" %in% as.character(s2$stats$variable))

# With two banks a chain consumes seeds[i] .. seeds[i] + 1, so seeds one apart
# overlap and two chains would share a WALNUTS stream. Distinctness alone is
# not enough; the blocks have to be disjoint.
expect_error(irt_causal_bart(
  sim2$responses,
  sim2$y,
  sim2$z,
  n_burnin = 5L,
  n_sampling = 5L,
  warmup_start = 2L,
  n_chains = 3L,
  seeds = c(1L, 2L, 3L)
))
# spaced by n_traits, the same seeds are fine
expect_silent(irt_causal_bart(
  sim2$responses,
  sim2$y,
  sim2$z,
  n_burnin = 5L,
  n_sampling = 5L,
  warmup_start = 2L,
  n_chains = 3L,
  seeds = c(1L, 3L, 5L)
))

# --- argument checks ---------------------------------------------------------
expect_error(do.call(irt_causal_bart, c(common, list(n_chains = 0L))))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 2L, n_cores = 0L))
))
# n_cores is validated on every path, not only where it is used
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 1L, n_cores = -3L))
))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(n_chains = 1L, n_cores = NA_integer_))
))
# a negative theta_sd would otherwise freeze theta for the whole run and still
# return an ordinary-looking ate
expect_error(do.call(irt_causal_bart, c(common, list(theta_sd = -1))))
expect_error(do.call(irt_causal_bart, c(common, list(theta_sd = 0))))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(theta_accept_target = 0))
))
expect_error(do.call(
  irt_causal_bart,
  c(common, list(theta_accept_target = 1))
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
