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
4. Generalize the Metropolis ratio: one move per trait, in a systematic scan.
   Trait `k`'s ratio uses bank `k`'s log-likelihood, `dnorm(theta_k)`, and the
   two BART terms read at the current value of every other trait.
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

## Status

LANDED 2026-07-27. ~330 lines of R, no C++, no dbarts change, as scoped.

Step 4 was written as a block proposal - all K traits for a person accepted or
rejected together - which the Constraints section of this same file puts out of
scope, and which "The one real constraint" in the design doc and the
`theta-block-move` TODO both rule out. Implemented as the systematic scan those
three agree on: trait `k` is proposed, accepted, and installed for every person
before trait `k + 1` is touched, and the "current" BART terms are re-read
between traits so each move conditions on the accepted values of the ones before
it. The re-read is skipped for the last trait, so K = 1 does no extra work.

Verification:
  - K = 1 regression: bitwise identical. Every element of the fit
    (`ate`, `alpha`, `beta`, `sigma`, `theta`, `theta_accept`, `theta_sd`,
    `theta_sd_trace`, `tuning`) `identical()` to the pre-change build on
    `simulate_irt_causal(200, 30, ate = -0.2, seed = 1)`, `seed = 1`, and so is
    the simulated data itself.
  - tinytest 86/86 (46 existing + 40 new in `test-multi-trait.R`), 3.3s.
  - K = 2 recovery, 400 persons / banks of 30 and 40, true ate 0.2: posterior
    mean 0.378, 95% interval [0.068, 0.640], covers. Per-bank cor(beta) 0.987
    and 0.982, cor(theta) 0.937 and 0.947.
  - Colleague's DGP (`scratch/FINAL_multi_latent.R`, 1000 persons, 100 items per
    bank, probit propensity in both traits, true gamma_z 0.2): posterior mean
    0.106, 95% interval [-0.064, 0.289], covers. Per-bank cor(beta) 0.969 and
    0.987, cor(alpha) 0.937 and 0.959, cor(theta) 0.975 both. 157s for
    500 + 1000 scans. This is the measurement `theta-block-move` was waiting on;
    component-wise mixing holds ~0.44 acceptance and costs nothing visible.

Left as scoped-out: same-item multidimensional IRT, correlated trait priors, and
a per-trait `theta_sd` (the proposal SD is still one shared, adapted scalar).

