# bairrtt

Causal inference with a latent, IRT-measured confounder, using BART and WALNUTS.

`bairrtt` fits a model in which a latent person trait `theta` is simultaneously

* the ability parameter of a two-parameter logistic (2PL) item-response model
  (measured by a matrix of 0/1 item responses), and
* a confounder between a binary treatment `z` and a continuous outcome `y`.

Two Bayesian Additive Regression Trees (BART) surfaces from
[dbarts](https://cran.r-project.org/package=dbarts) model the outcome and the
treatment assignment as flexible functions of `theta`; the IRT item parameters
are drawn conditional on `theta` by the gradient-based
[WALNUTS](https://github.com/bob-carpenter/walnuts) sampler (a No-U-Turn
variant), and `theta` itself by per-person Metropolis. Everything runs in a
Gibbs loop, exposed as a single entry point.

## Installation

`bairrtt` needs **dbarts (>= 1.0-0)**, which is newer than the version on CRAN
(it uses `updatePredictorPerObservationJointly()`, new in that release). Until
1.0-0 is on CRAN, install dbarts from its release branch first:

```r
# install.packages("remotes")
remotes::install_github("vdorie/dbarts@bartcore")   # dbarts 1.0-0, not yet on CRAN
remotes::install_github("vdorie/bairrtt")
```

Requires a C++20 compiler (WALNUTS uses concepts) and `RcppEigen`. On recent
toolchains the default `R CMD INSTALL` just works. If your default compiler is
too old for C++20, point R at a newer one in `~/.R/Makevars`, e.g. Homebrew
LLVM:

```make
CXX20 = /opt/homebrew/opt/llvm/bin/clang++
```

## Quick start

```r
library(bairrtt)

# simulate data from the model (treatment depends on the latent trait)
sim <- simulate_irt_causal(n_persons = 500, n_items = 60, ate = -0.2, seed = 1)

# fit
fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
                       n_burnin = 150, n_sampling = 400, seed = 1)

mean(fit$ate)                       # posterior mean treatment effect (~ -0.2)
quantile(fit$ate, c(0.025, 0.975))  # 95% credible interval
```

See `vignette("bairrtt")` for a walk-through, including recovery of the item
parameters and the WALNUTS item sampler on its own.

## Several latent traits

More than one trait can confound the same treatment and outcome, each measured
by its own item bank. Pass `responses` as a list of response matrices; every
item loads on exactly one trait, and the banks may differ in item count so long
as they cover the same persons in the same order. The per-bank quantities
(`alpha`, `beta`, and `theta`) come back as lists.

```r
sim <- simulate_irt_causal(n_persons = 500, n_items = c(60, 40), n_traits = 2,
                           ate = -0.2, prognostic = c(1.2, 0.8), seed = 1)

fit <- irt_causal_bart(sim$responses, sim$y, sim$z,
                       n_burnin = 150, n_sampling = 400, seed = 1)

mean(fit$ate)
cor(colMeans(fit$beta[[1]]), sim$beta[[1]])   # per-bank item recovery
```

Only disjoint banks are supported. Genuine multidimensional IRT, where one item
loads on several traits, is identified only up to a rotation of the loadings and
needs a constraint this interface does not express.

## Several chains

IRT models invite multimodality --- nothing separates `theta` from `-theta` but
the sign constraint on the discriminations --- and a single chain that settled in
the wrong mode looks healthy from the inside. Run several and compare them.

```r
fits <- irt_causal_bart(sim$responses, sim$y, sim$z,
                        n_burnin = 500, n_sampling = 1000,
                        n_chains = 4, n_cores = 4, seed = 1)
fits                                  # chains, pooled ate, max R-hat

irt_rhat(fits)$max                    # rank-normalized split R-hat
mean(irt_chain_draws(fits, "ate"))    # pooled over all chains
```

A chain's seed determines it completely, so `n_cores` changes only the wall
clock, never the draws. It forks, so it is a no-op on Windows.

## The WALNUTS item sampler

The item-parameter sampler is exported for reuse in other models. Hold `theta`
fixed within a draw but change it between draws:

```r
s <- irt_item_sampler(responses, alpha, beta, seed = 1)
for (i in seq_len(n_warmup)) irt_warmup(s, theta)  # adapt (theta may change)
irt_freeze(s)                                      # fix the tuning
draw <- irt_draw(s, theta)                         # c(alpha, beta)
```

## Licensing

`bairrtt` is released under the GPL (>= 2). It bundles the header-only WALNUTS
library (`inst/include/walnuts`), which is distributed under the MIT License
(Copyright (c) 2025 Bob Carpenter); see `inst/COPYRIGHTS` and
`inst/include/WALNUTS_LICENSE`.
