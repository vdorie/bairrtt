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

* `irt_causal_bart()` runs several chains via `n_chains`, four by default,
  optionally in parallel via `n_cores` (defaulting to `getOption("mc.cores")`).
  A chain's seed determines it completely, so `n_cores` changes only the wall
  clock, never the draws; it is fork-based, and so a no-op on Windows. A single
  chain admits no between-chain diagnostic, and IRT models invite the sign and
  label multimodality such a diagnostic exists to catch.

* The fit is now an `"irt_causal_fit"` object with the chain axis always
  present, at any chain count, so nothing downstream branches on it. Read it
  with `extract()` --- `dbarts`'s generic, which `stan4bart` also registers on
  --- and `summary()`, which delegates to `posterior::summarise_draws()` for
  rank-normalized split R-hat with bulk and tail effective sample size, and
  degrades to quantiles when `posterior` is not installed.

* `summary()` defaults to `ate` and `sigma`. R-hat is conventionally read
  against 1.01, but the maximum of many R-hats has a null distribution that is
  not centred at one, so a max over hundreds of item and person parameters
  flags perfectly converged fits. Those parameters are available via `vars` and
  are better judged on effective sample size.

* `n_sampling` now defaults to 2000 and `warmup_start` to 125. Measured on
  300 persons and 30 items, doubling the kept draws roughly doubles every
  effective sample size at exactly twice the cost; longer burn-in, per-person
  proposal adaptation, and other `n_trees` or `theta_accept_target` settings
  were all measured and bought nothing.

* `irt_causal_bart()` accepts observed covariates via `x`. They enter both BART
  surfaces, so the adjustment set becomes `(theta, x)`, and the measurement
  model as a latent regression `theta_j ~ N(x_j' gamma, 1)` --- omitting the
  latter would shrink exactly the trait-covariate relationships the model is
  asked about (Mislevy 1991). `gamma` is drawn from its exact conjugate full
  conditional and returned alongside the other draws. `x = NULL` keeps the
  previous model and the previous draws.

* `simulate_irt_causal()` generates data from the fitted model, with `n_traits`
  banks, a per-trait `prognostic` coefficient, and `n_covariates`
  covariates that confound through the trait, the assignment, and the outcome.

* Requires dbarts 1.0-0 or newer for
  `dbarts::updatePredictorPerObservationJointly()`, and a C++20 compiler for
  the bundled WALNUTS headers.
