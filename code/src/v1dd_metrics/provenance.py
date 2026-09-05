"""Provenance primitives: code version, run directories, and JSON coercion.

The asset carries its own provenance record naming the seed, config and package versions
that produced it. See docs/pipeline.md.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import Any

import numpy as np

from .version import CODE_DIR, code_version, read_version_file

__all__ = ["jsonable", "code_version", "git_sha", "run_stamp", "run_dir",
           "list_runs", "latest_run", "read_version_file", "CODE_DIR"]

_STAMP_RE = re.compile(r"^(?P<name>.+)_(?P<stamp>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$")


def run_stamp(when: float | None = None) -> str:
    """A sortable, filesystem-safe timestamp for one run. Local time, second resolution.

    Call once per run and reuse it; deriving it twice can straddle a second boundary.
    """
    return time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime(when))


def run_dir(root: str, name: str, stamp: str | None = None) -> str:
    """Path for this run's outputs, ``<root>/<name>_<stamp>``. Does not create it."""
    return os.path.join(root, f"{name}_{stamp or run_stamp()}")


def list_runs(root: str, name: str) -> list[str]:
    """Existing run directories for ``name``, oldest first.

    Sorted by the stamp in the directory name, not mtime, which copying rewrites.
    """
    if not os.path.isdir(root):
        return []
    found = []
    for entry in os.listdir(root):
        match = _STAMP_RE.match(entry)
        if match and match.group("name") == name and os.path.isdir(os.path.join(root, entry)):
            found.append((match.group("stamp"), os.path.join(root, entry)))
    return [path for _, path in sorted(found)]


def latest_run(root: str, name: str) -> str | None:
    """The most recent existing run directory for ``name``, or None."""
    runs = list_runs(root, name)
    return runs[-1] if runs else None


def git_sha(repo: str | None = None) -> str | None:
    """The commit for a provenance sidecar, or None. Lenient ``code_version``."""
    return code_version(Path(repo) if repo else None, strict=False)


def jsonable(obj: Any) -> Any:
    """Recursively convert numpy scalars and containers to plain JSON types.

    Casting up front means an artifact never fails to write after an expensive run.
    """
    # float first, and ahead of any passthrough: np.float64 subclasses float, so a
    # passthrough branch would let NaN reach json.dump, which writes a bare NaN literal
    # that strict parsers reject.
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if obj is None or isinstance(obj, (bool, np.bool_)):
        return bool(obj) if obj is not None else None
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, str):
        return obj
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v) for v in obj]
    return str(obj)
