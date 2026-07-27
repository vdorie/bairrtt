# multiple-chains

rng: neutral (existing single-chain draws unchanged; the surface is additive)
budget: ~350 lines of R, 0 lines of C++

## Goal

Running several chains from different seeds is a supported, documented
operation, and their agreement is reportable as a rank-normalized split R-hat.
Catching the label/sign multimodality IRT invites stops requiring the user to
hand-roll an outer loop and a diagnostic.

## Context

- The backlog entry (root `TODO`, `multiple-chains`) names the blocker: each
  chain owns mutable dbarts samplers and an `Rcpp::XPtr`-held WALNUTS state, so
  parallel chains are separate processes, not threads. See the `IrtSampler`
  gotcha in `CLAUDE.local.md`.
- `irt_causal_bart()` calls `set.seed(seed)` itself and hands `seed` to the
  WALNUTS RNG, so a chain's entire stream is determined by its seed. Parallel
  chains therefore need no L'Ecuyer-CMRG streams and no seed brokering - the
  seeds are the contract. The reproducibility expectations in `test-causal.R`
  already assert this. It is why the item is cheap.
- Sibling house style puts chains on the main entry point: `chains`/`cores` in
  `stan4bart()` (`../stan4bart/R/stan4bart.R`), `n.chains`/`n.threads` in
  `dbarts::dbartsControl()`. stan4bart stacks per-chain results into arrays with
  a trailing chain dimension. See the Decision below.
- `alpha`, `beta`, and `theta` are already per-bank lists when K > 1
  (`docs/plans/multi-trait.md`, "Status"). Any combiner recurses over those.
- `.github/workflows/check-standard.yaml` runs windows-latest, where
  `parallel::mclapply()` does not fork.

## Decision

Where do chains live? VD signs off before implementation.

(a) `irt_causal_bart(..., n_chains, n_cores)`. Matches both siblings and user
    expectation; one entry point. Costs a branching return type - `n_chains = 1`
    keeps today's shape, `> 1` returns the chains object - which every caller
    then has to handle.
(b) `irt_causal_bart_chains()`, a wrapper, plus the diagnostic helper. What the
    TODO prescribes. `irt_causal_bart()` keeps an unconditional return type and
    is untouched; reversible, since it can later become the implementation
    behind (a). Costs a second entry point to discover, and diverges from the
    siblings.

Recommendation: (a), on consistency with dbarts and stan4bart, which ship in
lockstep and whose users are the same people. What would change it: a
requirement that `irt_causal_bart()`'s return type stay unconditional, which is
a real cost (a good one) that (b) buys.

RESOLVED 2026-07-27: (a). VD signed off on the sibling-consistency argument.

## Constraints

- Out of scope: effective sample size (bulk/tail ESS), any within-chain
  diagnostic, running chains on a cluster/PSOCK, and resolving the sign/label
  multimodality that R-hat will now expose - reporting it is the deliverable.
- The single-chain path must produce identical draws for the same seed. Get that
  before adding anything.
- `n_cores > 1` is fork-based (`parallel::mclapply`) and so is a no-op on
  Windows; fall back to sequential with one message, do not error.
- Keep the exported engine API (`irt_item_sampler` and friends) unchanged.

## Steps

1. Chain runner: seeds default to `seq_len(n_chains)` offset from `seed`, are
   validated distinct, and each chain is an ordinary `irt_causal_bart()` call.
   Sequential by default; `parallel::mclapply` when `n_cores > 1`. Add
   `parallel` to Imports (base R, no new dependency weight).
2. Result object of class `"irt_causal_chains"`: the per-chain fits, the seeds,
   and the call. No pooled copy of the draws - accessors extract on demand.
3. `irt_draws(x, what)` returning an `n_sampling` x `n_chains` matrix for a
   scalar quantity, or a list of `n_sampling` x `n_chains` x `n_par` arrays for
   `alpha`/`beta`/`theta`. This is what both the diagnostic and pooling need.
4. `irt_rhat(x)`: rank-normalized split R-hat, taken as the max of the plain and
   folded (`|x - median(x)|`) rank-normalized values, per Vehtari et al. (2021),
   so it matches `posterior::rhat`. Returns per-quantity values plus the max.
5. `print.irt_causal_chains`: chain count, draws per chain, pooled ATE mean and
   interval, max R-hat, and a flag when it exceeds 1.01.
6. Roxygen for the new surface, a NEWS entry, and a diagnostics section in the
   vignette showing a converged and a deliberately-not-converged fit.

## Verification

- `R CMD INSTALL .` (no `--preclean`: no header or Makevars edit).
- Single-chain regression: one chain at `seed = 1` is `identical()` to
  `irt_causal_bart(..., seed = 1)`, element for element.
- `tinytest::test_package("bairrtt")`: 86/86 existing plus the new file.
- R-hat unit checks against synthetic input, not fits: four iid N(0,1) chains
  give R-hat < 1.05; four chains with one shifted by 10 give R-hat > 1.5; and
  the value agrees with `posterior::rhat` to 1e-8 on random input (test skipped
  when `posterior` is not installed; it goes in Suggests).
- A real 4-chain fit on `simulate_irt_causal()` reports max R-hat over `ate` and
  `sigma` below 1.05, and `n_cores = 2` gives draws identical to `n_cores = 1`.

## Status

LANDED 2026-07-27, decision (a). ~400 lines of R, no C++, as scoped.

`irt_causal_bart()` split into a validating wrapper and an internal
`irt_causal_bart_chain()` worker, mirroring the R-surface/engine split in
`R/engine.R`; all argument validation moved up into the wrapper, so the worker
assumes normalized inputs. Everything else new is in `R/chains.R`.

Verification:
  - Single-chain regression: bitwise identical. The same comparison
    `docs/plans/multi-trait.md` records still passes element for element, so the
    refactor moved no RNG.
  - Chain 1 of a multi-chain run `identical()` to the single-chain fit at the
    same seed; `n_cores = 4` `identical()` to `n_cores = 1`; 3.0x on 4 cores.
  - tinytest 138/138 (52 new in `test-chains.R`), 5.9s. R CMD check OK.
  - R-hat agrees with `posterior::rhat()` to 4.4e-16 over 200 random cases
    spanning iid, shifted, t(2), and scale-mismatched chains. It reads 1.73 on
    simulated sign multimodality and 1.29 when chains differ only in scale (the
    case the folded half of the statistic exists for).
  - Convergence behaves: on 300 persons / 30 items, max R-hat falls 1.106 ->
    1.025 -> 1.004 as burn-in/sampling go 150/300 -> 500/1000 -> 2000/4000. The
    package defaults (500/1000) land at 1.025, i.e. just outside the usual 1.01
    threshold - worth knowing, and now visible.

Two things worth recording that the plan did not anticipate:
  - Chain seeds are strided by `n_traits`, not 1. Bank `k` inside a chain takes
    `seed + k - 1`, so consecutive chain seeds would have made chain `c`'s
    second bank share a WALNUTS stream with chain `c + 1`'s first - a
    between-chain dependence in exactly the quantity the diagnostic reports on.
  - The accessor is `irt_chain_draws()`, not `irt_draws()`: one letter from the
    existing `irt_draw()`, which advances a standalone item sampler and is
    unrelated.

Left as scoped-out: effective sample size, now its own TODO item
(`effective-sample-size`), and any attempt to resolve the multimodality R-hat
now exposes rather than merely report it.
