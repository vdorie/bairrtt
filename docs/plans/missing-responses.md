# missing-responses

rng: neutral (complete-data draws must be bitwise identical; `NA` was an error
before, so there is no existing posterior to change) - but `IrtLogpGrad` is
edited, so the finite-difference gate applies and grows a missing-data case
budget: ~10 lines of C++, ~50 of R, ~130 of tests, ~120 of docs

## Goal

`NA` in a response matrix is accepted everywhere responses are taken and dropped
from the item likelihood under ignorable missingness, in the C++ gradient block
and in the `theta` Metropolis ratio alike. Complete data fits exactly as before.

## Context

- Design, the assumption being made, and why listwise deletion and
  omits-as-zero are both wrong: `docs/design/missing-responses.md`. Read
  "What is dropped, and what that assumes" and "Where the mask goes" first.
- The three edit sites are named in that doc's "Where the mask goes":
  `as_response_matrix()` (`R/engine.R`), `IrtLogpGrad::operator()`
  (`inst/include/bairrtt_types.h`), and the `sign_mat` term of the per-person
  ratio in `irt_causal_bart_chain()` (`R/irt_causal_bart.R`).
- The Assumptions section of `?irt_causal_bart` (`R/irt_causal_bart.R` roxygen)
  gains a bullet; `docs/design/covariates.md`'s Status line calls this item the
  remaining blocker and needs updating when it lands.
- Gate classes and the brevity rule: `docs/plans/README.md`.

## Constraints

- Complete-data draws must be bitwise identical to the pre-change build. Get
  that gate before adding any missing-data behavior.
- `NA` and `NaN` are one predicate: `is.na()` in R, `std::isnan()` in C++.
- Do not use `na.rm = TRUE` in the Metropolis ratio; zero the sign matrix, which
  is exact (the design doc explains the cancellation) and does not swallow an
  `NA` arriving from anywhere else.
- Never weaken the finite-difference tolerance in `test-engine.R`.
- Out of scope: models for the missingness mechanism, missing `y`/`z`/`x`, a
  `p_missing` argument to `simulate_irt_causal()` (inject `NA` in the tests
  instead), zero-variance item warnings (TODO `measurement-diagnostics`), and
  any change under `inst/include/walnuts/`.

## Steps

1. `as_response_matrix()` admits `NA`: validate only the observed cells are
   `0`/`1`, with the observed mask leading the `&` chain so `any()` stays
   logical. Error message becomes "only 0, 1, and NA".
2. `IrtLogpGrad::operator()` skips missing cells - read `Y(j, i)`, `continue` on
   `std::isnan`, before any accumulation into `lp`, `g_alpha`, or `sum_resid`.
   Item priors and the gradient's prior terms are untouched. Update the model
   comment above the struct.
3. `sign_mat` in `irt_causal_bart_chain()` is zeroed where the bank is `NA`,
   with a comment naming the cancellation.
4. Warn once per bank, from `as_response_banks()`, for items with no observed
   response and for persons with no observed response in that bank. Not from
   `as_response_matrix()`: `irt_item_logdensity()` runs it inside
   finite-difference loops.
5. Docs: `@param responses` for `irt_causal_bart()`, `irt_item_sampler()`, and
   `irt_item_logdensity()`; a new Assumptions bullet on ignorable missingness
   pointing at the design doc; the vignette; `NEWS.md`; `man/` via roxygen.
6. Retire the `missing-responses` item from `TODO`, update the Status lines in
   `docs/design/covariates.md` and `docs/design/missing-responses.md`.

## Verification

- `R CMD INSTALL . --preclean` (the header changed).
- Complete-data regression: a single-chain fit at `seed = 1` on
  `simulate_irt_causal(200, 30, ate = -0.2, seed = 1)` is `identical()`
  element-for-element to the pre-change build. This is the acceptance gate.
- `tinytest::test_package("bairrtt")`: 46 existing plus `test-missing.R`.
- The existing finite-difference check passes untouched, and a second one over a
  matrix with ~20% `NA` (including one all-`NA` item) passes at the same
  tolerance.
- Exactness, the check that proves it is the observed-data likelihood: for a
  complete matrix and a masked copy, `lp_complete - lp_masked` equals the sum of
  the masked cells' `plogis(sign * alpha * (theta - beta), log.p = TRUE)` terms
  to ~1e-10. Same for the gradient, term by term.
- Prior-only item: with item `i` all `NA`, `gradient[i] == 1 - alpha_i` and
  `gradient[n_items + i] == -beta_i / beta_sd^2` exactly.
- Recovery at 20% MCAR missing on 500 persons / 30 items: the 95% `ate`
  interval covers the truth and `cor()` of the posterior-mean item parameters
  with the truth is high. Also fit the same data with `NA` scored `0` and report
  both biases in Status - that is the comparator a referee asks about, but do
  not assert on it.
- A person with an all-`NA` bank row fits without error and yields finite
  `theta` draws under `keep_theta = TRUE`; both degenerate cases warn.
- `R CMD check`, `air format .`, and lint clean.
