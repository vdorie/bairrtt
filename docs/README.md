# docs/ orientation

Wayfinding map for `docs/`. It records no decision itself - it explains where
decisions live.

## Layout

- `docs/design/` - design proposals, rationale, and landing/decision records.
  One doc per feature or investigation: the problem, the forks considered, the
  decision, and (for landed work) the landing notes.
- `docs/plans/` - implementation plans: one file per open TODO item or landed
  feature, with Goal/Context/Steps/Verification and a closing Status note.
- `docs/plans/README.md` - the process doc: how a plan is written, which gates
  each class of change has to clear, the review and brevity rubric. Process,
  not content; it is not an index.
- repo-root `TODO` - the live, unordered backlog of OPEN work. Every entry names
  a plan file here. Landed items are trimmed out once closed; their record lives
  in git history and in the design/plan docs.

## Three surfaces, three jobs

Do not expect any one of them to do another's job:

1. **The root `TODO`** - what is still open, right now. Forward-facing only.
2. **`docs/plans/README.md`** - how work gets planned, gated, and reviewed.
3. **`docs/design/`** - why the code is the way it is. Read the relevant design
   doc before treating any behavior as accidental.

Both directories are small enough today that a listing is a directory read;
add an `INDEX.md` to either once it passes roughly twenty files.
