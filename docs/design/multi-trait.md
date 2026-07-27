# Multiple latent traits

Status: LANDED 2026-07-27, as proposed. Landing notes and the verification
numbers are in `docs/plans/multi-trait.md`, "Status"; the one deviation from the
plan (its step 4 read as a block proposal, which this doc's "The one real
constraint" rules out) is recorded there. Everything below held: the item blocks
factorized, the change was pure R, and dbarts was not touched.

## Problem

`irt_causal_bart()` fits exactly one latent trait. The motivating case is a
colleague's Stan model with two: two disjoint item banks, two abilities, both
confounding the same treatment and outcome. Their model is a single joint HMC
block over everything -

```
correct1 ~ bernoulli_logit(alpha1[ii1] * (theta1[jj1] - b1[ii1]))
correct2 ~ bernoulli_logit(alpha2[ii2] * (theta2[jj2] - b2[ii2]))
z        ~ bernoulli_logit(beta0 + beta[1]*theta1 + beta[2]*theta2)
y        ~ normal(gamma0 + gamma_z*z + gamma[1]*theta1 + gamma[2]*theta2, sigma)
```

- and the estimand is `gamma_z`. bairrtt already replaces the two linear terms
(prognostic and propensity) with BART surfaces and the joint HMC block with the
three-block Gibbs sampler described in `CLAUDE.local.md`. The remaining gap is
K = 1 vs K = 2.

## What generalizes for free

**The item-parameter blocks factorize exactly.** In the motivating model
`correct1` loads only on `theta1` and `correct2` only on `theta2`. Conditional
on the traits, the two `(alpha, beta)` posteriors are independent. So K traits
means K independent `irt_item_sampler` objects, each over its own response
matrix - not one larger sampler. `IrtLogpGrad`, `IrtSampler`, and every exported
`irt_*` function carry over unmodified. This is exact, not an approximation, and
it is the reason this item is expected to need no C++.

Corollary: the per-sampler warmup/freeze/draw lifecycle just runs K times per
scan, and `irt_freeze()` is still called exactly once per sampler.

## What has to change

All in `R/irt_causal_bart.R` and `R/simulate.R`:

1. `theta` becomes an `n_persons x K` matrix and contributes K columns to both
   BART models. `setCutPoints` is called per column.
2. The per-person Metropolis acceptance ratio sums K IRT log-likelihood terms
   and K standard-normal prior terms; the `theta_mat`/`alpha_mat`/`beta_mat`
   broadcasting machinery becomes per-trait.
3. Installing the accepted trait values into both BART models has to go
   component-wise - see the constraint below.
4. The estimand is unchanged in form: `f(1, theta)` and `f(0, theta)` now
   condition on the whole trait matrix.
5. `simulate_irt_causal()` grows a K argument and returns a list of response
   matrices plus a trait matrix.

## The one real constraint

`dbarts::updatePredictorPerObservationJointly()` accepts a **single** shared
column and errors on more (its own guard: "joint updates can only be applied to
a single column"). That function is what makes an empty-leaf `theta` move get
rejected rather than collapsed, and what keeps the two BART models' `theta`
columns identical.

So with K traits the install must be a component-wise systematic scan: propose
and install trait 1 for all persons, then trait 2, and so on. That is a valid
Metropolis-within-Gibbs sweep and needs no dbarts change. What it cannot do is
move person j's whole trait vector as a single accept/reject unit; if
component-wise mixing turns out to be poor, the remedy is multi-column support
in dbarts (TODO `theta-block-move`), which should not be opened speculatively.

## Identification

The motivating model has **disjoint item banks** - each item loads on exactly
one trait. That is the identified case: within a bank, `alpha > 0` fixes the
reflection and `theta ~ N(0, 1)` fixes location and scale, exactly as at K = 1.

If both traits ever load on the *same* items - genuine multidimensional IRT -
the loading matrix is identified only up to rotation, and a constraint is
required (lower-triangular loadings, or a fixed anchor item per trait). That is
a materially larger change: it turns the data structure from "a list of response
matrices" into "one matrix plus a loading pattern," and it changes the item
likelihood from scalar `alpha_i` to a loading vector. **Confirm which case is
wanted before designing the interface.** The disjoint-bank case is the one
scoped here.

## Traits are modeled independent

The Stan model and the current code both put independent `N(0, 1)` priors on the
traits. Correlated traits would need a covariance parameter and a
prior on it, and the `theta` Metropolis proposal would want to match that
correlation. Out of scope; note it if the colleague expects correlated abilities.

## Rejected alternative

Putting `theta` into the WALNUTS block along with the item parameters. It does
not work at any K: `theta` feeds the BART surfaces, which are piecewise constant
in it, so the conditional has no gradient. See `CLAUDE.local.md` gotchas.
