# V1DD functional metrics

Read [`docs/pipeline.md`](docs/pipeline.md) for how the pipeline is structured, and
`.claude/memory/` (indexed by `MEMORY.md`) for decisions and cautions that are not
derivable from the code.

Essentials, repeated because they are easy to get wrong:

- **`code/` is `/code` in the capsule.** Everything the pipeline imports must live under
  `code/`; that is why the package is `code/src/`, not repo-root `src/`.
- **Do not run `git commit`.** Leave finished work uncommitted, say what changed and
  where, and stop. The user reviews diffs and commits.
- **The pipeline refuses to start without a code version.** `$V1DD_CODE_VERSION`, else
  `code/CODE_VERSION`, else `git rev-parse`.
- **Code comments stay laconic.** A docstring saying what a function does and how its
  arguments differ. Rationale, decisions and cautions belong in `docs/`, not inline.
- **`code/src/v1dd_metrics/` must not import from `code/validation/`.** A test enforces it.
- This repo is a fresh start. It owes no output compatibility to `allen_v1dd`, the 2019
  white paper, or our own earlier assets; comparability is discussed in
  [`docs/comparability.md`](docs/comparability.md), not enforced in code.
