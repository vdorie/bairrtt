# Two latent traits, two disjoint item banks. Covers the list-of-matrices
# interface, the per-bank return shapes, and (loosely) recovery. Kept small so
# the suite stays fast; strict recovery is left to the vignette / longer runs.

sim <- simulate_irt_causal(
  n_persons = 200L,
  n_items = c(15L, 20L),
  n_traits = 2L,
  ate = -0.2,
  prognostic = c(1.2, 0.8),
  seed = 7
)

# --- simulator shapes --------------------------------------------------------
expect_true(is.list(sim$responses))
expect_equal(length(sim$responses), 2L)
expect_equal(dim(sim$responses[[1L]]), c(200L, 15L))
expect_equal(dim(sim$responses[[2L]]), c(200L, 20L))
expect_equal(dim(sim$theta), c(200L, 2L))
expect_equal(lengths(sim$alpha), c(15L, 20L))
expect_equal(lengths(sim$beta), c(15L, 20L))
expect_equal(sim$prognostic, c(1.2, 0.8))

# n_traits = 1 keeps the bare single-trait shapes
sim1 <- simulate_irt_causal(n_persons = 50L, n_items = 10L, seed = 7)
expect_true(is.matrix(sim1$responses))
expect_true(is.null(dim(sim1$theta)))
expect_equal(length(sim1$alpha), 10L)

fit <- irt_causal_bart(
  sim$responses,
  sim$y,
  sim$z,
  n_burnin = 40L,
  n_sampling = 80L,
  warmup_start = 20L,
  seed = 7L,
  keep_theta = TRUE
)

# --- shapes: one entry per bank ----------------------------------------------
expect_equal(length(fit$ate), 80L)
expect_equal(length(fit$alpha), 2L)
expect_equal(dim(fit$alpha[[1L]]), c(80L, 15L))
expect_equal(dim(fit$alpha[[2L]]), c(80L, 20L))
expect_equal(dim(fit$beta[[1L]]), c(80L, 15L))
expect_equal(dim(fit$beta[[2L]]), c(80L, 20L))
expect_equal(length(fit$theta), 2L)
expect_equal(dim(fit$theta[[1L]]), c(80L, 200L))
expect_equal(dim(fit$theta[[2L]]), c(80L, 200L))
expect_equal(length(fit$tuning), 2L)
expect_equal(length(fit$theta_accept), 120L)

# --- invariants --------------------------------------------------------------
expect_true(all(is.finite(fit$ate)))
expect_true(all(fit$alpha[[1L]] > 0)) # discriminations positive
expect_true(all(fit$alpha[[2L]] > 0))
expect_true(all(fit$sigma > 0))
# every bank's WALNUTS sampler froze exactly once and did adapt
expect_true(all(vapply(fit$tuning, function(t) isTRUE(t$frozen), logical(1L))))
expect_true(all(vapply(fit$tuning, function(t) t$warmup_iter > 0, logical(1L))))

# recovery is loose at this size, but each bank's difficulties should track
# the truth and the ATE posterior should sit in a sane range.
expect_true(mean(fit$ate) > -1.5 && mean(fit$ate) < 0.5)
expect_true(cor(colMeans(fit$beta[[1L]]), sim$beta[[1L]]) > 0.8)
expect_true(cor(colMeans(fit$beta[[2L]]), sim$beta[[2L]]) > 0.8)

# --- reproducibility ---------------------------------------------------------
fit2 <- irt_causal_bart(
  sim$responses,
  sim$y,
  sim$z,
  n_burnin = 40L,
  n_sampling = 80L,
  warmup_start = 20L,
  seed = 7L,
  keep_theta = TRUE
)
expect_identical(fit$ate, fit2$ate)
expect_identical(fit$alpha, fit2$alpha)
expect_identical(fit$theta, fit2$theta)

# --- a one-element list is the K = 1 model, just wrapped ---------------------
# It returns the K > 1 (list) shapes, but must draw exactly what the bare
# matrix draws: same trait count, same RNG stream.
fit_bare <- irt_causal_bart(
  sim1$responses,
  sim1$y,
  sim1$z,
  n_burnin = 10L,
  n_sampling = 20L,
  warmup_start = 5L,
  seed = 8L
)
fit_list <- irt_causal_bart(
  list(sim1$responses),
  sim1$y,
  sim1$z,
  n_burnin = 10L,
  n_sampling = 20L,
  warmup_start = 5L,
  seed = 8L
)
expect_identical(fit_bare$ate, fit_list$ate)
expect_identical(fit_bare$alpha, fit_list$alpha)

# --- argument checks ---------------------------------------------------------
# banks must agree on the number of persons
expect_error(irt_causal_bart(
  list(sim$responses[[1L]], sim$responses[[2L]][-1L, ]),
  sim$y,
  sim$z
))
expect_error(irt_causal_bart(list(), sim$y, sim$z)) # no banks
expect_error(simulate_irt_causal(50L, c(10L, 12L), n_traits = 3L)) # n_items
expect_error(simulate_irt_causal(50L, 10L, n_traits = 2L, prognostic = 1:3))
