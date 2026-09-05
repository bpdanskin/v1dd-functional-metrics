---
name: code-version-stamp-ladder
description: "Why the pipeline refuses to start without a code version, and the three-step ladder that resolves one."
metadata: { node_type: memory, type: project, modified: 2026-09-03 }
---

A reproducible run copies `code/` **without `.git`**, so `git rev-parse` returns nothing
and an asset would ship `commit_hash: null` with nothing tying it to the code that made it.
The entry point therefore refuses to start when the version is unknown — one second lost
rather than five hours.

Ladder, in `run_pipeline.resolve_code_version`:

1. `$V1DD_CODE_VERSION`
2. `code/CODE_VERSION` (first non-comment line) — tracked, diffable, immune to image rebuilds
3. `git rev-parse HEAD`
4. otherwise exit

**The value must match `^[0-9a-fA-F]{7,40}$`.** This is not pedantry: the 2026-09-01 asset
shipped `"version": "= 17cacea..."` because Docker's legacy `ENV <key> <value>` form takes
everything after the first space as the value. That string is neither empty nor whitespace,
so a mere is-it-set guard passed it.

**Why not a Dockerfile `ENV`, as the fork used?** This repo's Dockerfile is generated from
`.codeocean/environment.json`, and Code Ocean may regenerate it — a hand-added `ENV` line
is not guaranteed to survive. The tracked file does.

`code/CODE_VERSION` ships comment-only on purpose, so it never asserts a stale commit.
Stamp it before a capture run: `git rev-parse HEAD > code/CODE_VERSION`.
