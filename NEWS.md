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

* `simulate_irt_causal()` generates data from the fitted model.

* Requires dbarts 1.0-0 or newer for
  `dbarts::updatePredictorPerObservationJointly()`, and a C++20 compiler for
  the bundled WALNUTS headers.
