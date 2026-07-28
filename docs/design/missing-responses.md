# Missing item responses

Status: PROPOSED 2026-07-27. Implementation plan in
`docs/plans/missing-responses.md`. Raised by the pre-release statistical review
as one of two blockers for applied use; the other, `covariates`, has landed.

## Problem

`as_response_matrix()` (`R/engine.R`) rejects `NA`, so every entry point that
takes responses - `irt_causal_bart()`, `irt_item_sampler()`,
`irt_item_logdensity()` - refuses any real response matrix. Test booklets rotate
forms, respondents skip items, and instruments are administered in matrix-
sampled designs on purpose. A package that cannot read an `NA` cannot be pointed
at data.

The workarounds a user would reach for are both wrong. Listwise deletion of
persons removes nearly the whole sample under item-level missingness (at 20% and
30 items, almost no one is complete) and, worse, deletes them from the *causal*
estimand as well as the measurement model - the ATE is then over a
completers-only population. Scoring omits as `0` is the achievement-testing
convention for a reason, but applied to accidental missingness it biases `theta`
downward, and here `theta` is the confounder, so the bias lands directly on the
treatment effect.

## What is dropped, and what that assumes

The 2PL likelihood factors over items given `theta`, so the observed-data
likelihood for person `j` is the product over that person's *observed* cells:

    L_j = prod_{i : y_ij observed} p_ij^{y_ij} (1 - p_ij)^{1 - y_ij}

Dropping the missing terms - which is what mirt, TAM, and ltm all do by default
- is exactly right when the missingness is **ignorable**: missing at random
given the observed data and the model, with the missingness parameters distinct
from the model's (Rubin 1976; Little and Rubin 2019, ch. 6). Nothing about the
mechanism has to be modeled under that condition, which is why the change is
small.

It is not free of assumptions, and this model has one wrinkle a plain IRT fit
does not.

* **MCAR by design is the easy case.** Planned missingness - booklet rotation,
  matrix sampling, a form a person was never shown - is missing by a mechanism
  the analyst controls, so ignorability holds by construction. This is the case
  the feature is for.

* **Missingness that depends on `theta` is not handled.** Low-ability
  respondents skipping hard items is missing *not* at random with respect to
  the latent trait, and the dropped-terms likelihood is then wrong. This is the
  usual omitted-response problem in achievement testing and it needs a model for
  the response propensity (Holman and Glas 2005; Pohl et al. 2014), which is out
  of scope. Users with omits rather than not-reached items should consider
  scoring the omits `0` deliberately, which is a substantive choice the package
  should not make for them.

* **Missingness that depends on `z` or `y` is worse here than elsewhere.**
  `theta`'s full conditional in this sampler already includes the outcome and
  assignment likelihoods (see the Metropolis ratio in `R/irt_causal_bart.R`, and
  the joint-vs-cut discussion filed as TODO `cut-diagnostic`). Missingness
  driven by treatment status is therefore not merely a measurement nuisance: it
  is differential measurement by arm, which is the exclusion-restriction
  violation the Assumptions section already warns about, arriving by a second
  route.

The honest summary for `?irt_causal_bart` is that missing responses are assumed
ignorable, and that this is an assumption about the *design*, not something the
data can check.

## Where the mask goes

Three places, and the third is the one that is easy to miss.

1. **Validation.** `as_response_matrix()` admits `NA` and keeps it. `NaN` counts
   as missing too, so the R and C++ sides agree on one predicate.

2. **The gradient block.** `IrtLogpGrad::operator()`
   (`inst/include/bairrtt_types.h`) skips a missing `(j, i)` cell entirely: no
   `lp` term, no `g_alpha` term, no `sum_resid` term. The item priors are
   untouched, so an item with no observed responses is drawn from its prior -
   well defined, and worth a warning rather than an error.

3. **The `theta` Metropolis ratio** (`irt_causal_bart_chain`,
   `R/irt_causal_bart.R`). The per-person IRT log-likelihood is a `rowSums()`
   over a sign matrix `2 * bank - 1`, which is `NA` wherever the bank is. Zeroing
   the sign at those cells is exact rather than merely convenient: a zeroed cell
   contributes `plogis(0, log.p = TRUE) = -log 2` to the current and to the
   proposed sum alike, and the ratio is a difference of the two, so the constant
   cancels. No `na.rm = TRUE`, which would also swallow an `NA` arriving from
   somewhere it should not.

## Degenerate data

Neither case is an error; both are warnings, because both are legitimate in a
matrix-sampled design and neither breaks the sampler.

* **An item nobody answered** has a flat likelihood, so its `alpha` and `beta`
  are prior draws. Related, and left to TODO `measurement-diagnostics`: an item
  answered identically by everyone has a likelihood monotone in `beta` and is
  held only by the prior. The all-missing case is the limit of that.

* **A person with no observed responses in a bank** gets a `theta` informed only
  by the prior (or latent regression) and by `y` and `z`. That is the correct
  joint-model behavior - it is what multiple imputation would do - and the
  person still contributes to the ATE average. It is also a good sign the data
  were assembled wrong, hence the warning.

## Complete data does not move

A response matrix with no `NA` must produce bitwise-identical draws: the C++
skip never fires and the zeroed sign matrix is the sign matrix. That is the
acceptance gate for the change, and it is what keeps this out of the
posterior-changing class.

## Out of scope

Models for the missingness mechanism (MNAR, response-propensity, omit-vs-
not-reached); missing `y`, `z`, or `x`; and any change under
`inst/include/walnuts/`.
