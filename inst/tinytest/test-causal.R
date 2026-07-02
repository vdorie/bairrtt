# Full causal sampler: output structure, invariants, and (loosely) recovery.
# Kept small so it runs quickly on install/check; strict recovery is left to
# the vignette / longer runs.

sim <- simulate_irt_causal(n_persons = 200L, n_items = 20L, ate = -0.2, seed = 3)

fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
                       n_burnin = 40L, n_sampling = 80L, warmup_start = 20L,
                       seed = 3L, keep_theta = TRUE)

# --- shapes ------------------------------------------------------------------
expect_equal(length(fit$ate), 80L)
expect_equal(dim(fit$alpha), c(80L, 20L))
expect_equal(dim(fit$beta), c(80L, 20L))
expect_equal(length(fit$sigma), 80L)
expect_equal(dim(fit$theta), c(80L, 200L))
expect_equal(length(fit$theta_accept), 120L)
expect_equal(length(fit$theta_sd_trace), 120L)

# --- invariants --------------------------------------------------------------
expect_true(all(is.finite(fit$ate)))
expect_true(all(fit$alpha > 0))                 # discriminations positive
expect_true(all(fit$sigma > 0))
expect_true(isTRUE(fit$tuning$frozen))          # WALNUTS frozen for sampling
expect_true(fit$tuning$warmup_iter > 0)         # adaptation actually happened
expect_true(all(fit$theta_accept >= 0 & fit$theta_accept <= 1))

# recovery is loose at this size, but the ATE posterior should be in a sane
# range and item parameters should correlate with the truth.
expect_true(mean(fit$ate) > -1.5 && mean(fit$ate) < 0.5)
expect_true(cor(colMeans(fit$beta), sim$beta) > 0.8)

# --- reproducibility: same seed -> identical draws ---------------------------
fit2 <- irt_causal_bart(sim$responses, sim$y, sim$z,
                        n_burnin = 40L, n_sampling = 80L, warmup_start = 20L,
                        seed = 3L, keep_theta = TRUE)
expect_identical(fit$ate, fit2$ate)
expect_identical(fit$alpha, fit2$alpha)
expect_identical(fit$theta, fit2$theta)

# --- keep_theta = FALSE (default) drops the theta draws ----------------------
fit_nt <- irt_causal_bart(sim$responses, sim$y, sim$z,
                          n_burnin = 20L, n_sampling = 20L, warmup_start = 10L,
                          seed = 4L)
expect_null(fit_nt$theta)
expect_equal(length(fit_nt$ate), 20L)

# --- n_burnin = 0 runs and freezes the item sampler immediately --------------
fit_nb <- irt_causal_bart(sim$responses, sim$y, sim$z,
                          n_burnin = 0L, n_sampling = 20L, seed = 5L)
expect_equal(length(fit_nb$ate), 20L)
expect_true(isTRUE(fit_nb$tuning$frozen))
expect_equal(length(fit_nb$theta_accept), 20L)

# --- argument checks ---------------------------------------------------------
expect_error(irt_causal_bart(sim$responses, sim$y[-1], sim$z))         # y length
expect_error(irt_causal_bart(sim$responses, sim$y, sim$z * 2))         # z not 0/1
y_na <- sim$y; y_na[1] <- NA
expect_error(irt_causal_bart(sim$responses, y_na, sim$z))              # NA in y
expect_error(irt_causal_bart(sim$responses, sim$y, sim$z,
                             n_theta_cutpoints = 0L))                  # bad cutpoints
# warmup_start >= n_burnin warns and clamps (does not error)
expect_warning(irt_causal_bart(sim$responses, sim$y, sim$z,
                               n_burnin = 10L, n_sampling = 10L,
                               warmup_start = 50L, seed = 6L))
