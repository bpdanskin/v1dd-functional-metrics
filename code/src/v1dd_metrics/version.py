"""Resolving the commit an asset was built from. Standard library only.

Kept dependency-free so the entry point's pre-flight check runs before anything heavy is
imported. See docs/pipeline.md for why the pipeline refuses to start without a version.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

__all__ = ["code_version", "read_version_file", "CODE_DIR"]

#: `code/` -- this package lives at code/src/v1dd_metrics/.
CODE_DIR = Path(__file__).resolve().parents[2]

_SHA_RE = re.compile(r"^[0-9a-fA-F]{7,40}$")


def read_version_file(path: Path) -> str:
    """First non-comment, non-blank line of a CODE_VERSION file, or ''."""
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return line
    return ""


def code_version(code_dir: Path | None = None, strict: bool = True) -> str | None:
    """The commit this code was built from.

    Tried in order: ``$V1DD_CODE_VERSION``, ``code/CODE_VERSION``, ``git rev-parse``.
    ``strict`` raises SystemExit when nothing answers or the value is not a git object
    name; otherwise returns None, which suits a best-effort provenance sidecar.
    """
    code_dir = code_dir or CODE_DIR
    for value in (os.environ.get("V1DD_CODE_VERSION", "").strip(),
                  read_version_file(code_dir / "CODE_VERSION")):
        if value:
            if _SHA_RE.match(value):
                return value
            if strict:
                raise SystemExit(
                    "refusing to start: code version " + repr(value)
                    + " is not a commit (expected 7-40 hex characters)."
                )
            return None
    try:
        out = subprocess.run(["git", "-C", str(code_dir), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10, check=True)
        if out.stdout.strip():
            return out.stdout.strip()
    except Exception:
        pass
    if strict:
        raise SystemExit(
            "refusing to start: the code version is unknown.\n"
            "A reproducible run copies code/ without .git, so set V1DD_CODE_VERSION\n"
            "or write the commit into code/CODE_VERSION:\n"
            "    git rev-parse HEAD > code/CODE_VERSION"
        )
    return None
