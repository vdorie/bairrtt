# bairrtt 0.1.0

* Initial version. Fits a causal model in which a latent person trait,
  measured by a two-parameter logistic item-response model, confounds a binary
  treatment and a continuous outcome.

* `irt_causal_bart()` runs the three-block Gibbs sampler: per-person Metropolis
  for the latent trait, `dbarts` tree updates for the outcome and assignment
  surfaces, and WALNUTS for the IRT item parameters.

* The item-parameter sampler is exported on its own for reuse in other models:
  `irt_item_sampler()`, `irt_warmup()`, `irt_freeze()`, `irt_draw()`,
  `irt_tuning()`, and `irt_item_logdensity()` (which exposes the analytic
  log-posterior and gradient for finite-difference checking).

* `irt_causal_bart()` accepts more than one latent trait: pass `responses` as a
  list of item banks, one per trait, all confounding the same treatment and
  outcome. The traits are updated in a systematic scan, one bank gets one
  WALNUTS sampler, and the per-bank `alpha`, `beta`, and `theta` draws come back
  as lists. A single response matrix keeps its previous interface and its
  previous draws. Only disjoint banks (each item loading on one trait) are
  supported.

* `simulate_irt_causal()` generates data from the fitted model, with `n_traits`
  banks and a per-trait `prognostic` coefficient.

* Requires dbarts 1.0-0 or newer for
  `dbarts::updatePredictorPerObservationJointly()`, and a C++20 compiler for
  the bundled WALNUTS headers.
