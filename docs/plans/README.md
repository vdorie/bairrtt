# Plan process

The `TODO` at the repo root is an unordered backlog; every item names a plan
file here. This document is the contract between whoever plans an item and
whoever implements it.

Adapted from dbarts's `docs/plans/README.md` and scaled to this package -
bairrtt is ~900 lines of R and C++, so there is no agent-routing table and no
worktree-per-item requirement. Everything else carries over.

## Plan format

One file per TODO item, `<item>.md`, hard cap 80 lines. Front block: `rng:`
(see the gate classes below) and `budget:` (expected diff size). Sections:

- Goal: two or three sentences; what is true after the item lands.
- Context: durable references (a code symbol, a doc section title) and design
  pointers. No narrative that a pointer can replace. Never a bare line number.
- Decision (decision-gated items only): the question, a recommendation, and
  what evidence would change it. VD signs off before implementation.
- Constraints: gates, contract freezes, explicit out-of-scope list.
- Steps: numbered; each independently verifiable.
- Verification: exact commands and expected outcomes.

## Cross-references

Cite a durable landmark, never a bare line number: line numbers shift on every
insertion above them, while the thing cited did not move. For code, name the
symbol - `IrtLogpGrad::operator()` (`bairrtt_types.h`), not
`bairrtt_types.h:114`. For a doc, name the section title. A line may trail a
symbol as a convenience, never stand as the anchor. State each fact in one home
doc; elsewhere link to it by title rather than restating it - a copied fact is a
second thing to keep in sync, and the one that rots.

## Gate classes

- **neutral**: draws unchanged (documentation, validation, R-surface renames).
  Gate: `tinytest::test_package("bairrtt")`, 46/46.
- **shifting**: draws change, the posterior does not (reordering a Gibbs block,
  changing a proposal's parameterization).
  Gates: the above, plus a recovery check - `simulate_irt_causal()` with a known
  `ate` and enough draws, confirming the 95% interval covers it and the item
  parameters correlate with truth.
- **posterior-changing**: the stationary distribution or a default changes (a
  prior, a link, the target density, the acceptance ratio).
  Gates: all of the above, plus a design note in `docs/design/`, plus - for any
  change to `IrtLogpGrad` - the finite-difference gradient check in
  `test-engine.R` passing. That check is the only barrier between a sign error
  and a silently wrong posterior; never weaken its tolerance to make it pass.

Changes to the vendored `inst/include/walnuts/` are out of scope by default:
that is upstream's code. If a change is genuinely needed there, it goes upstream
first, and the vendored copy is refreshed wholesale.

## Implementation protocol

- One item at a time. The prompt is the plan file path plus `CLAUDE.local.md`;
  the plan is the spec, so do not restate it.
- After any C++ change: `R CMD INSTALL .` before tinytest (`--preclean` after
  editing `inst/include/bairrtt_types.h` or `src/Makevars`).
- After changing an `// [[Rcpp::export]]`: rerun `Rcpp::compileAttributes()`.
- Commit each fix as it is verified rather than batching at the end.
- Stop conditions: a step fails twice; the diff exceeds 1.5x budget; a needed
  change is out of scope. Report and stop; do not improvise.

## Brevity (binding)

Final report <= 20 lines: files touched, gate results (pass/fail lines
verbatim), deviations from plan. Nothing else - no methodology, no plan recap,
no prose about what the code now does.
