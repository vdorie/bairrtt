# Observed covariates

Status: LANDED 2026-07-27, as proposed. Landing notes and verification numbers
are in `docs/plans/covariates.md`, "Status". Raised by the pre-release
statistical review as one of two blockers for applied use; the other,
`missing-responses`, has since landed too.

## Problem

`irt_causal_bart()` accepts `responses`, `y`, and `z` and nothing else. The
outcome surface is `reformulate(c("z", theta_names))`, so the identifying
assumption is literally "the latent trait is the entire confounder." No applied
study satisfies that. The Assumptions section of `?irt_causal_bart` states the
limitation honestly; this item removes it.

## Where x has to enter

Two places, and the second is the one that gets forgotten.

1. **The BART surfaces.** `y ~ f(z, theta, x)` and `z ~ f(theta, x)`, so that
   the adjustment set is `(theta, x)` rather than `theta` alone.

2. **The measurement model, as a latent regression** `theta_j ~ N(x_j' gamma,
   1)`. Omitting this is Mislevy's (1991) conditioning-model bias: a
   conditioning model that leaves out variables the analyst later relates to
   `theta` shrinks those relationships toward zero. It is why NAEP conditions
   plausible values on essentially every downstream variable. Since `theta`'s
   posterior here already sees `y` and `z` (see "Traits see the outcome" in
   `multi-trait.md`'s sibling discussion, and the Metropolis ratio in
   `R/irt_causal_bart.R`), leaving `x` out of the prior while `y` and `z` are in
   the likelihood is the worst of both worlds.

## Scale identification

`theta_j ~ N(x_j' gamma, 1)` fixes the **residual** SD at 1, not the marginal
one. That is deliberate and it is what keeps the 2PL identified: with `alpha`
free, `theta -> c * theta` and `alpha -> alpha / c` is a symmetry, and pinning
the residual variance pins `c`. The marginal variance of `theta` becomes
`var(x'gamma) + 1 >= 1`, which is the point - covariates explain part of the
trait. This is the standard latent-regression IRT parameterization.

With no covariates, `gamma` is absent and the prior is exactly today's
`N(0, 1)`, so the existing path is untouched.

## The gamma block is conjugate

Given `theta_k` and `X`, and a prior `gamma_k ~ N(0, g^2 I)`, the full
conditional is normal:

```
V = (X'X + I / g^2)^{-1},   m = V X' theta_k,   gamma_k | . ~ N(m, V)
```

so this is one more Gibbs block, drawn exactly, in R. It costs a `p x p` solve
per trait per scan, with `p` the number of covariate columns - negligible beside
a BART update. One `gamma` vector per trait; the traits stay independent given
`X`.

## What does not change

- **`IrtLogpGrad` and the whole C++ side.** The item likelihood is
  `alpha_i (theta_j - beta_i)`; `x` never appears in it. `theta` remains data to
  the WALNUTS block. So this item needs no C++ and no re-derived gradient, and
  the finite-difference check in `test-engine.R` is untouched.
- **`updatePredictorPerObservationJointly()`.** It still installs one trait
  column at a time. The `x` columns are static for the whole run and are never
  updated, so the empty-leaf machinery is unaffected.
- **The estimand.** `E[f(1, theta, x) - f(0, theta, x)]` is still averaged over
  the sample, so it remains the sample ATE. It is still invariant to the latent
  metric's location and scale.

## Interface

`x` as a matrix or data frame of `n_persons` rows, defaulting to `NULL`.
A data frame is passed through to `dbarts()`, which already handles factors and
builds the design matrix; the latent regression needs a numeric design matrix,
so it uses `model.matrix()` on the same frame, with an intercept.

Deliberately NOT offered: separate covariate sets for the outcome surface, the
assignment surface, and the latent regression. Three formulas is a materially
bigger interface, and the whole point of the latent regression is that it should
contain everything downstream - splitting the sets invites exactly the omission
this item exists to prevent. One `x` for all three.

## Rejected alternative

Putting `x` only in the BART surfaces and leaving `theta ~ N(0, 1)`. Cheaper by
one Gibbs block and wrong for the reason in point 2: it reintroduces
conditioning-model bias, and it would make the package's own `x` argument a trap
- users would reasonably assume adjustment was complete.

## Open question, deferred

Whether `gamma` should be reported with a warning that it is a *conditioning*
model rather than a structural one. Latent-regression coefficients are routinely
over-interpreted as effects of `x` on ability. Leaning yes, in the docs only.
