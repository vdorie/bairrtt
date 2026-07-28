# Observed covariates: the interface, the latent regression, and the thing the
# feature exists for -- that adjusting removes bias that ignoring x leaves.

# --- simulator ---------------------------------------------------------------
sim0 <- simulate_irt_causal(n_persons = 50L, n_items = 10L, seed = 2)
expect_null(sim0$x) # no covariates by default
expect_null(sim0$gamma)

sim <- simulate_irt_causal(
  n_persons = 200L,
  n_items = 20L,
  ate = -0.2,
  n_covariates = 2L,
  covariate_effect = c(1, -0.8),
  seed = 2
)
expect_equal(dim(sim$x), c(200L, 2L))
expect_equal(sim$gamma, c(1, -0.8))
expect_error(simulate_irt_causal(50L, 10L, n_covariates = -1L))
expect_error(simulate_irt_causal(
  50L,
  10L,
  n_covariates = 2L,
  covariate_effect = 1:3
))

common <- list(
  sim$responses,
  sim$y,
  sim$z,
  n_burnin = 40L,
  n_sampling = 80L,
  warmup_start = 10L,
  n_chains = 2L,
  seed = 2L
)

fit <- do.call(irt_causal_bart, c(common, list(x = sim$x)))

# --- shapes: gamma is intercept + one per covariate --------------------------
expect_equal(dim(fit$gamma), c(80L, 2L, 3L))
expect_equal(dim(extract(fit, "gamma")), c(160L, 3L))
expect_true(all(is.finite(fit$gamma)))

s <- summary(fit, vars = c("ate", "gamma"))
expect_equal(nrow(s$stats), 4L) # ate + 3 gamma
expect_true("gamma[3]" %in% as.character(s$stats$variable))

# without covariates there is no gamma at all
plain <- do.call(irt_causal_bart, common)
expect_null(plain$gamma)
expect_error(extract(plain, "gamma"))

# --- x = NULL is the previous model exactly ---------------------------------
# not merely similar: passing a zero-column x must take the same path
plain2 <- do.call(irt_causal_bart, c(common, list(x = NULL)))
expect_identical(plain$ate, plain2$ate)
expect_identical(plain$alpha, plain2$alpha)

# --- the latent regression recovers gamma ------------------------------------
# theta ~ N(x'gamma, 1), so the slopes should track the truth; the intercept is
# 0 by construction
gamma_hat <- colMeans(extract(fit, "gamma"))
expect_true(abs(gamma_hat[2L] - sim$gamma[1L]) < 0.35)
expect_true(abs(gamma_hat[3L] - sim$gamma[2L]) < 0.35)

# --- adjustment removes bias that ignoring x leaves --------------------------
# The point of the feature. If both fits covered, this test would prove nothing,
# so the unadjusted one is asserted to MISS.
big <- simulate_irt_causal(
  n_persons = 500L,
  n_items = 30L,
  ate = -0.2,
  n_covariates = 2L,
  covariate_effect = c(1, -0.8),
  seed = 15
)
big_args <- list(
  big$responses,
  big$y,
  big$z,
  n_burnin = 150L,
  n_sampling = 300L,
  warmup_start = 40L,
  n_chains = 2L,
  seed = 15L
)
unadj <- do.call(irt_causal_bart, big_args)
adj <- do.call(irt_causal_bart, c(big_args, list(x = big$x)))

covers <- function(f) {
  ci <- quantile(as.vector(f$ate), c(0.025, 0.975), names = FALSE)
  ci[1L] <= big$ate && ci[2L] >= big$ate
}
expect_false(covers(unadj)) # confounded by x, and it shows
expect_true(covers(adj)) # adjusted, and it covers
expect_true(abs(mean(adj$ate) - big$ate) < abs(mean(unadj$ate) - big$ate))

# --- argument checks ---------------------------------------------------------
expect_error(do.call(irt_causal_bart, c(common, list(x = sim$x[-1L, ]))))
x_na <- sim$x
x_na[1L, 1L] <- NA
expect_error(do.call(irt_causal_bart, c(common, list(x = x_na))))
expect_error(do.call(irt_causal_bart, c(common, list(x = "nope"))))
expect_error(do.call(irt_causal_bart, c(common, list(gamma_sd = 0))))
expect_error(do.call(irt_causal_bart, c(common, list(gamma_sd = -1))))
# a covariate named like a trait column would silently shadow it in the BART
# model frames
x_clash <- data.frame(theta = sim$x[, 1L], ok = sim$x[, 2L])
expect_error(do.call(irt_causal_bart, c(common, list(x = x_clash))))

# a data frame with a factor works: dbarts builds its own design matrix, and
# the latent regression takes model.matrix() of the same frame
x_df <- data.frame(
  num = sim$x[, 1L],
  grp = factor(rep(c("a", "b"), length.out = 200L))
)
fit_df <- do.call(irt_causal_bart, c(common, list(x = x_df)))
expect_equal(dim(fit_df$gamma)[3L], 3L) # intercept + num + grpb
expect_true(all(is.finite(fit_df$ate)))

# --- the recorded call never carries the data -------------------------------
# a fit is meant to be saved; do.call() would otherwise inline the matrices
expect_true(object.size(fit$call) < 5000)
expect_true(grepl("irt_causal_bart", paste(deparse(fit$call), collapse = "")))
