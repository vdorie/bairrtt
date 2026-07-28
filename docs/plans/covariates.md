# covariates

rng: posterior-changing (a new prior on theta and a new Gibbs block; the
no-covariate path must stay bitwise intact)
budget: ~250 lines of R, 0 lines of C++

## Goal

`irt_causal_bart()` accepts observed covariates `x`, which enter both BART
surfaces and the measurement model as a latent regression
`theta_j ~ N(x_j' gamma, 1)`. `x = NULL` keeps the current behavior and the
current draws. No C++ change.

## Context

- Design and the identification argument: `docs/design/covariates.md`. Read it
  first, in particular "Where x has to enter" (the latent regression is the half
  that gets forgotten) and "Scale identification" (residual variance is what is
  fixed at 1, not marginal).
- The assumption this removes is written up under "Assumptions" in
  `?irt_causal_bart` (`R/irt_causal_bart.R` roxygen); that text needs updating
  when this lands.
- The Gibbs scan order is documented under "Layout" in `CLAUDE.local.md`. The
  `gamma` block joins it after the item parameters.
- Multi-trait shapes: `alpha`/`beta`/`theta` are per-bank lists at K > 1
  (`docs/plans/multi-trait.md`, "Status"). `gamma` follows the same convention.

## Constraints

- Out of scope: separate covariate sets per surface (see the design doc's
  "Interface"), covariates in the item likelihood (that is DIF, TODO
  `item-models`), missing values in `x`, and any change under
  `inst/include/walnuts/` or to `IrtLogpGrad`.
- The `x = NULL` path must produce identical draws to the current code for the
  same seed. This is the acceptance gate; get it before adding any x behavior.
- `theta`'s residual SD stays fixed at 1. Do not add a free scale parameter -
  it is what identifies the 2PL.
- Keep the exported engine API (`irt_item_sampler` and friends) unchanged.

## Steps

1. Accept `x` (matrix, data frame, or `NULL`) validated in the wrapper: row
   count matches `n_persons`, no `NA`, not a response-matrix-shaped mistake.
   Build once: `x_frame` for dbarts and `x_design` = `model.matrix()` with
   intercept for the latent regression.
2. Add the `x` columns to both BART formulas and both model frames, and to
   every `predict()` frame - the two proposal predicts, the between-trait
   re-read, and the two estimand predicts. `setCutPoints` still applies only to
   the trait columns.
3. Initialize `gamma` at zero, one vector per trait, length `ncol(x_design)`.
4. Replace the `dnorm(theta, log = TRUE)` prior terms in the Metropolis ratio
   with `dnorm(theta, mean = x_design %*% gamma[[k]], sd = 1, log = TRUE)` on
   both sides. With `x = NULL` this must reduce to the current call exactly.
5. Add the conjugate `gamma` Gibbs block after the item-parameter block:
   `V = solve(crossprod(x_design) + diag(1 / gamma_sd^2))`,
   `m = V %*% crossprod(x_design, theta[, k])`, draw via the Cholesky of `V`.
   One draw per trait per scan. `gamma_sd` is a new prior argument.
6. Store `gamma` draws with the same chain axis as everything else, and add
   them to `extract()` and to `summary()`'s permitted `vars`.
7. `simulate_irt_causal()` grows `n_covariates`, generating `x`, a `gamma`
   that shifts `theta`, and coefficients putting `x` into the outcome and
   assignment models so the DGP actually needs adjustment. Returns `x` and the
   true `gamma`.
8. Update the Assumptions section (the "theta is the whole confounder" bullet
   softens to "(theta, x) is the whole confounder"), the vignette, `NEWS.md`,
   and `man/` via roxygen.

## Verification

- `R CMD INSTALL .` (no `--preclean`: no header or Makevars edit).
- `x = NULL` regression: a single-chain fit at `seed = 1` on
  `simulate_irt_causal(200, 30, ate = -0.2, seed = 1)` is `identical()`
  element-for-element to the pre-change build, as is the simulated data.
- `tinytest::test_package("bairrtt")`: 157/157 existing plus the new file.
- The finite-difference gradient check in `test-engine.R` must still pass
  untouched - if it does not, something reached the C++ that should not have.
- Recovery with covariates: simulate with a known `ate` and a `gamma` whose
  covariates genuinely confound, fit, and confirm the interval covers the ATE
  and that `colMeans` of the `gamma` draws correlate with the truth.
- Confounding bias check, the one that shows the feature works: on a DGP where
  `x` confounds and `theta` does not fully mediate it, the `x = NULL` fit should
  MISS the true ATE and the `x`-adjusted fit should cover it. If both cover, the
  test is not exercising the feature.
- `R CMD check`, lint, and `air format .` clean.

## Status

LANDED 2026-07-27, as designed. ~300 lines of R, no C++, as scoped.

Verification:
  - `x = NULL` regression: bitwise identical. Every element of a single-chain
    fit still `identical()` to the pre-multi-trait baseline, so neither the
    latent regression nor the extra Gibbs block moved the no-covariate RNG.
  - tinytest 187/187 (30 new in `test-covariates.R`), R CMD check OK, lint
    clean. The finite-difference gradient check in `test-engine.R` passes
    untouched, confirming nothing reached the C++.
  - The check that shows the feature works, on 500 persons / 30 items with two
    covariates confounding through the trait, the assignment, and the outcome:
    the unadjusted fit MISSES (95% [0.287, 1.066] against a true -0.2) and the
    adjusted one covers ([-0.358, 0.283]). Asserted in both directions in the
    test file, so a regression that silently ignored `x` would fail it.
  - `gamma` recovery: posterior means (0.99, -0.82) against a true (1, -0.8).

Two things the plan did not anticipate:
  - `match.call()` inlines the DATA when the fit is built through `do.call()`,
    which is how the tests call it - 137 KB of deparsed response matrix in one
    printed summary, and megabytes stored in a fit meant to be saved. The
    recorded call now abbreviates any inlined value to a `<matrix>` placeholder
    while keeping symbolic arguments (`sim$responses`) readable. This was
    finding 8 of the correctness review, filed as informational; the covariate
    work made it visible.
  - Covariate names are rejected if they collide with `y`, `z`, or a trait
    column, which would otherwise silently shadow a trait in the BART frames.

Left as scoped-out: separate covariate sets per surface, covariates in the item
likelihood (that is DIF, TODO `item-models`), and missing values in `x`.
