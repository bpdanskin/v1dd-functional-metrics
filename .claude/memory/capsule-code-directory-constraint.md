---
name: capsule-code-directory-constraint
description: "In a Code Ocean capsule the repo's code/ IS /code, so everything the pipeline imports must live under code/."
metadata: { node_type: memory, type: project, modified: 2026-09-03 }
---

**The repo's `code/` directory is mounted as `/code` in the capsule.** `code/run` executes
with cwd `/code`, and a module at `code/x.py` sees `Path(__file__).parents[1] == "/"`.

Consequence: **everything the pipeline imports must live under `code/`.** A repo-root
`src/` would simply not exist in the container. That is why the package is
`code/src/v1dd_metrics/`.

`code/run` puts it on the path in one line, derived from its own location so the same
command works off-capsule:

```bash
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$HERE/src${PYTHONPATH:+:$PYTHONPATH}"
```

This replaced seven separate `sys.path` insertions in the fork. **Do not reintroduce
`sys.path` surgery in modules** — if an import fails, the path setup is the thing to fix.

Also true of the capsule, and worth remembering when a test constructs a directory layout:
`/data` is read-only mounts, `/results` is captured as the asset, `/scratch` is discarded,
and **`code/` is copied without `.git`** — see [[code-version-stamp-ladder]].
