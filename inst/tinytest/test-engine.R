# WALNUTS item-parameter engine: gradient, lifecycle, validation.

set.seed(1)
n_persons <- 150L
n_items <- 15L
sim <- simulate_irt_causal(n_persons, n_items, seed = 1)

# --- analytic gradient matches finite differences ---------------------------
par0 <- c(log(rep(1, n_items)), rnorm(n_items))
ld <- irt_item_logdensity(par0, sim$responses, sim$theta)
expect_true(is.numeric(ld$value) && length(ld$value) == 1L)
expect_equal(length(ld$gradient), 2L * n_items)

eps <- 1e-5
g_fd <- numeric(length(par0))
for (k in seq_along(par0)) {
  up <- dn <- par0
  up[k] <- up[k] + eps
  dn[k] <- dn[k] - eps
  g_fd[k] <- (irt_item_logdensity(up, sim$responses, sim$theta)$value -
              irt_item_logdensity(dn, sim$responses, sim$theta)$value) / (2 * eps)
}
expect_true(max(abs(ld$gradient - g_fd)) < 1e-4)

# --- lifecycle: warmup -> freeze -> draw ------------------------------------
s <- irt_item_sampler(sim$responses, rep(1, n_items), rep(0, n_items), seed = 42L)
expect_inherits(s, "irt_item_sampler")
expect_false(irt_tuning(s)$frozen)

th <- sim$theta
d_warm <- irt_warmup(s, th)
expect_equal(length(d_warm), 2L * n_items)
expect_true(all(d_warm[1:n_items] > 0))          # alpha positive on natural scale
expect_true(irt_tuning(s)$warmup_iter >= 1)

irt_freeze(s)
expect_true(irt_tuning(s)$frozen)
d_draw <- irt_draw(s, th)
expect_equal(length(d_draw), 2L * n_items)
expect_true(all(d_draw[1:n_items] > 0))

# freezing twice / warming up after freeze are errors
expect_error(irt_freeze(s))
expect_error(irt_warmup(s, th))

# --- input validation --------------------------------------------------------
expect_error(irt_item_sampler(sim$responses, rep(1, n_items + 1L), rep(0, n_items)))
expect_error(irt_item_sampler(sim$responses, rep(-1, n_items), rep(0, n_items))) # alpha <= 0
bad <- sim$responses; bad[1, 1] <- 2
expect_error(irt_item_sampler(bad, rep(1, n_items), rep(0, n_items)))             # not 0/1
expect_error(irt_draw(s, th[-1]))                                                 # wrong theta length
expect_error(irt_warmup("not a sampler", th))
