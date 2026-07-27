# multi-trait

rng: posterior-changing (a new model; the K = 1 path must stay bitwise intact)
budget: ~250 lines of R, 0 lines of C++

## Goal

`irt_causal_bart()` accepts K latent traits, each measured by its own item bank,
all confounding the same treatment and outcome. `K = 1` keeps its current
behavior and its current draws. No C++ change; no dbarts change.

## Context

- Design, scoping, and the identification argument: `docs/design/multi-trait.md`.
  Read it first - in particular "What generalizes for free" (the item blocks
  factorize exactly, so K traits means K independent samplers) and "The one real
  constraint" (the joint install is single-column).
- The motivating Stan model and its simulation script are in `scratch/`.
- The Gibbs scan order is documented under "Layout" in `CLAUDE.local.md`.
- `irt_item_sampler()` and the whole compiled engine are unchanged by this item.

## Constraints

- Out of scope: same-item multidimensional IRT (rotational indeterminacy - see
  the design doc's "Identification"), correlated trait priors, block moves on a
  person's whole trait vector (TODO `theta-block-move`), any edit under
  `inst/include/walnuts/`.
- The K = 1 path must produce identical draws to the current code for the same
  seed. This is the acceptance gate for "no regression"; get it before adding
  any K > 1 behavior.
- Keep the exported engine API (`irt_item_sampler`, `irt_warmup`, `irt_freeze`,
  `irt_draw`, `irt_tuning`, `irt_item_logdensity`) unchanged.

## Steps

1. Accept `responses` as either a matrix (K = 1, current behavior) or a list of
   matrices (K = length(list)), validated in `as_response_matrix()`'s caller.
   All banks must agree on `n_persons`; item counts may differ.
2. Make `theta` an `n_persons x K` matrix internally. Build both BART models
   with K trait columns named `theta1..thetaK` (K = 1 keeps the name `theta`, so
   the existing draws and formulas are untouched), calling `setCutPoints` per
   column.
3. Create K `irt_item_sampler` objects, one per bank. Run the
   warmup/freeze/draw lifecycle on each, per scan, in bank order.
4. Generalize the Metropolis ratio: propose all K traits for a person at once,
   sum the K bank log-likelihoods and the K `dnorm(theta_k, log = TRUE)` prior
   terms, and add the two BART terms once.
5. Install component-wise: loop k, calling
   `dbarts::updatePredictorPerObservationJointly(list(response_model,
   assignment_model), theta[, k], paste0("theta", k))` and reverting the
   persons its mask rejects. Burn-in keeps using `setPredictor(forceUpdate =
   TRUE)`, also per column.
6. Return `alpha`, `beta`, and (with `keep_theta`) `theta` as length-K lists of
   matrices when K > 1, and unchanged when K = 1.
7. `simulate_irt_causal()` grows `n_traits` and a `prognostic` vector; returns a
   list of response matrices and an `n_persons x K` `theta` when K > 1.
8. Update `man/` via roxygen and add a K = 2 section to the vignette.

## Verification

- `R CMD INSTALL .` (no `--preclean` needed - no header or Makevars edit).
- `tinytest::test_package("bairrtt")`: 46/46 existing expectations still pass.
- K = 1 regression: `irt_causal_bart(..., seed = 1)` on
  `simulate_irt_causal(200, 30, ate = -0.2, seed = 1)` returns an `ate` vector
  identical to the pre-change build's. Record the pre-change draws before
  touching anything.
- K = 2 recovery: simulate with two banks and a known `ate`, fit, and confirm
  the 95% interval covers it and each bank's `alpha`/`beta` posterior means
  correlate with truth. Add as a new `inst/tinytest/test-multi-trait.R`, sized
  so the whole suite stays under a few seconds.
- Sanity check against the colleague's model: `scratch/FINAL_multi_latent.R`
  simulates two traits at ate 0.2; the fitted interval should cover it.
