# Missing item responses: NA is accepted and dropped from the item likelihood.
# The checks that matter are the exactness ones -- that the target really is the
# observed-data likelihood, not merely something that runs. See
# docs/design/missing-responses.md.

set.seed(7)
n_persons <- 120L
n_items <- 12L
beta_sd <- 10
sim <- simulate_irt_causal(n_persons, n_items, seed = 7)
Y <- sim$responses
th <- sim$theta

par0 <- c(log(runif(n_items, 0.5, 2)), rnorm(n_items))
alpha <- exp(par0[seq_len(n_items)])
beta <- par0[n_items + seq_len(n_items)]

# ~20% missing, plus one item nobody answered
mask <- matrix(runif(n_persons * n_items) < 0.2, n_persons, n_items)
mask[, 3L] <- TRUE
Y_miss <- Y
Y_miss[mask] <- NA_real_

# --- validation: NA in, 2 still out -----------------------------------------
expect_silent(irt_item_logdensity(par0, Y_miss, th))
expect_error(irt_item_sampler(
  replace(Y, 1L, 2),
  rep(1, n_items),
  rep(0, n_items)
))
expect_inherits(
  irt_item_sampler(Y_miss, rep(1, n_items), rep(0, n_items), seed = 7L),
  "irt_item_sampler"
)
# an all-NA matrix is degenerate but not an error: the target is the prior
expect_silent(irt_item_logdensity(par0, Y * NA_real_, th))

# --- exactness: the dropped cells are exactly the missing ones ---------------
ld_full <- irt_item_logdensity(par0, Y, th, beta_sd)
ld_miss <- irt_item_logdensity(par0, Y_miss, th, beta_sd)

eta <- outer(th, beta, "-") * matrix(alpha, n_persons, n_items, byrow = TRUE)
cell_lp <- plogis((2 * Y - 1) * eta, log.p = TRUE)
expect_true(abs((ld_full$value - ld_miss$value) - sum(cell_lp[mask])) < 1e-9)

resid <- Y - plogis(eta)
alpha_mat <- matrix(alpha, n_persons, n_items, byrow = TRUE)
drop_alpha <- colSums(resid * outer(th, beta, "-") * alpha_mat * mask)
drop_beta <- colSums(-resid * alpha_mat * mask)
expect_true(
  max(abs(
    (ld_full$gradient[seq_len(n_items)] - ld_miss$gradient[seq_len(n_items)]) -
      drop_alpha
  )) <
    1e-9
)
expect_true(
  max(abs(
    (ld_full$gradient[n_items + seq_len(n_items)] -
      ld_miss$gradient[n_items + seq_len(n_items)]) -
      drop_beta
  )) <
    1e-9
)

# --- an item nobody answered is drawn from its prior -------------------------
# item 3 is all NA, so its gradient is the prior's and nothing else
expect_equal(ld_miss$gradient[3L], 1 - alpha[3L])
expect_equal(ld_miss$gradient[n_items + 3L], -beta[3L] / beta_sd^2)

# --- the finite-difference gate, through a matrix with NA --------------------
# Same tolerance as the complete-data check in test-engine.R; a sign error in
# the masked branch has nothing else standing between it and the posterior.
eps <- 1e-5
g_fd <- numeric(length(par0))
for (k in seq_along(par0)) {
  up <- dn <- par0
  up[k] <- up[k] + eps
  dn[k] <- dn[k] - eps
  g_fd[k] <- (irt_item_logdensity(up, Y_miss, th, beta_sd)$value -
    irt_item_logdensity(dn, Y_miss, th, beta_sd)$value) /
    (2 * eps)
}
expect_true(max(abs(ld_miss$gradient - g_fd)) < 1e-4)

# --- degenerate data warns, and still fits -----------------------------------
sim_w <- simulate_irt_causal(60L, 10L, seed = 3)
bad_item <- sim_w$responses
bad_item[, 2L] <- NA_real_
expect_warning(irt_causal_bart(
  bad_item,
  sim_w$y,
  sim_w$z,
  n_burnin = 10L,
  n_sampling = 20L,
  warmup_start = 5L,
  n_chains = 1L,
  seed = 3L
))

bad_person <- sim_w$responses
bad_person[4L, ] <- NA_real_
fit_bp <- suppressWarnings(irt_causal_bart(
  bad_person,
  sim_w$y,
  sim_w$z,
  n_burnin = 10L,
  n_sampling = 20L,
  warmup_start = 5L,
  n_chains = 1L,
  seed = 3L,
  keep_theta = TRUE
))
# a person with nothing observed still gets a finite theta, from the prior and
# from y and z, and still contributes to the ATE average
expect_true(all(is.finite(extract(fit_bp, "theta")[, 4L])))
expect_true(all(is.finite(fit_bp$ate)))

# --- recovery under 20% MCAR missingness -------------------------------------
sim_r <- simulate_irt_causal(500L, 30L, ate = -0.2, seed = 11)
set.seed(11)
resp_r <- sim_r$responses
resp_r[matrix(runif(500L * 30L) < 0.2, 500L, 30L)] <- NA_real_

fit_r <- irt_causal_bart(
  resp_r,
  sim_r$y,
  sim_r$z,
  n_burnin = 150L,
  n_sampling = 300L,
  warmup_start = 40L,
  n_chains = 2L,
  seed = 11L
)
ci <- quantile(fit_r$ate, c(0.025, 0.975))
expect_true(ci[1L] < sim_r$ate && sim_r$ate < ci[2L])
expect_true(cor(colMeans(extract(fit_r, "alpha")), sim_r$alpha) > 0.5)
expect_true(cor(colMeans(extract(fit_r, "beta")), sim_r$beta) > 0.9)

# --- multi-trait: one bank missing, one complete -----------------------------
# exercises the per-bank sign matrix, which is where the ratio's mask lives
sim_m <- simulate_irt_causal(
  n_persons = 200L,
  n_items = c(15L, 20L),
  n_traits = 2L,
  ate = -0.2,
  seed = 5
)
set.seed(5)
banks <- sim_m$responses
banks[[1L]][matrix(runif(200L * 15L) < 0.15, 200L, 15L)] <- NA_real_
fit_m <- irt_causal_bart(
  banks,
  sim_m$y,
  sim_m$z,
  n_burnin = 40L,
  n_sampling = 80L,
  warmup_start = 10L,
  n_chains = 1L,
  seed = 5L
)
expect_equal(length(fit_m$alpha), 2L)
expect_true(all(is.finite(fit_m$ate)))
expect_true(all(is.finite(extract(fit_m, "alpha")[[1L]])))
