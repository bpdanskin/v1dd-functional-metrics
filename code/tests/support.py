"""Shared helpers for the ported suite.

``check`` carries a name alongside the condition, which is worth keeping: these assertions
describe a property in words, and the name is what tells you which property broke.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

#: repo root -- this file is code/tests/support.py
REPO = Path(__file__).resolve().parents[2]
#: the metadata CLI, now inside the package
METADATA_CLI = REPO / "code" / "src" / "v1dd_metrics" / "metadata.py"


def check(name: str, cond, detail: str = "") -> None:
    """Assert ``cond``, reporting ``name`` and any ``detail``."""
    assert bool(cond), f"{name}" + (f" -- {detail}" if detail else "")


def synthetic_plane(n_rois: int = 12, n_frames: int = 3000, dt: float = 0.16504,
                    seed: int = 0):
    """A minimal stand-in for a loaded imaging plane.

    Enough for the array-level families: traces, timestamps, ROI identity and validity.
    Nothing here reads NWB.
    """
    rng = np.random.default_rng(seed)
    ts = np.cumsum(np.full(n_frames, dt)) + 3.0
    events = rng.gamma(2.0, 5e-4, size=(n_frames, n_rois))
    dff = rng.normal(0.0, 0.05, size=(n_frames, n_rois))

    class Plane:
        column, volume, plane, n_rois_ = 1, "3", 2, n_rois
        timestamps = ts

        def __init__(self):
            self.n_rois = n_rois
            self.dt = dt
            self.traces = {"events": events, "dff": dff}
            self.is_valid = np.ones(n_rois, dtype=bool)
            self.roi_ids = np.arange(n_rois)
            self.depth_um = 82.0
            self.roi_confidence = np.full(n_rois, 0.9)

    return Plane()
