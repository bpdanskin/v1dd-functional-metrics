"""Locating the input data, on the capsule or in a local checkout.

Precedence: ``$V1DD_DATA_ROOT``, then a repo-local ``data/<probe>`` found by searching
upward, then ``/data`` (where Code Ocean mounts attached datasets).
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["find_up", "resolve_data_root", "resolve_dataset_dir"]

CAPSULE_DATA_ROOT = "/data"


def find_up(relative_path: str, max_levels: int = 4) -> str | None:
    """Search the working directory and its parents for ``relative_path``.

    Pass a nested path (``"data/409828_V1DD_Filtered"``) rather than a bare directory
    name, so walking upward cannot match something outside the repo. Returns an absolute
    path, or None.
    """
    base = Path.cwd()
    for level in range(max_levels + 1):
        candidate = base.joinpath(*([".."] * level), relative_path)
        if candidate.is_dir():
            return str(candidate.resolve())
    return None


def resolve_data_root(probe: str | None = None, max_levels: int = 4) -> str:
    """The directory the datasets sit inside.

    ``probe`` names a dataset to look for under a repo-local ``data/``; omit it to skip
    that step. The result is not guaranteed to exist -- use ``resolve_dataset_dir`` for a
    checked path.
    """
    if os.environ.get("V1DD_DATA_ROOT"):
        return os.environ["V1DD_DATA_ROOT"]
    if probe:
        local = find_up(str(Path("data") / probe), max_levels=max_levels)
        if local:
            return str(Path(local).parent)
    return CAPSULE_DATA_ROOT


def resolve_dataset_dir(*names: str, root: str | None = None,
                        required: bool = True) -> str | None:
    """The first of ``names`` that exists under ``root``, as an absolute path.

    Several names let one call accept more than one layout. Raises FileNotFoundError
    naming every path tried, so a missing dataset fails here rather than deep inside a
    later read; pass ``required=False`` to get None instead.
    """
    if not names:
        raise ValueError("give at least one dataset name")
    if root is None:
        root = resolve_data_root(names[0])

    tried = []
    for name in names:
        candidate = Path(root) / name
        if candidate.is_dir():
            return str(candidate.resolve())
        tried.append(str(candidate))

    if not required:
        return None
    raise FileNotFoundError(
        "none of these dataset directories exist:\n  " + "\n  ".join(tried)
        + f"\n(working directory: {os.getcwd()})\n"
        "Attach the dataset in Code Ocean, or set V1DD_DATA_ROOT to the directory "
        "that contains it."
    )
