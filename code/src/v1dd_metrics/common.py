"""Numeric helpers shared by more than one family."""

from __future__ import annotations

import numpy as np

def _nanmean(x: np.ndarray, axis) -> np.ndarray:
    """`np.nanmean` without the all-NaN-slice warning, which fires constantly here."""
    finite = np.isfinite(x)
    n = finite.sum(axis=axis)
    total = np.where(finite, x, 0.0).sum(axis=axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, total / np.maximum(n, 1), np.nan)


def _ratio(p, q, zero_to_nan: bool = True):
    """``p / q``, with a zero denominator yielding NaN or 0.

    ``zero_to_nan=False`` reproduces the historical behaviour, where an undefined ratio
    was reported as 0.
    """
    q = np.asarray(q, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    fill = np.nan if zero_to_nan else 0.0
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(q == 0, fill, p / np.where(q == 0, 1.0, q))


def _metric_index(a, b):
    """The SSI index: **NaN** when the denominator is 0, unlike `_ratio`."""
    s = np.asarray(a, dtype=np.float64) + np.asarray(b, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(s == 0, np.nan, (a - b) / np.where(s == 0, 1.0, s))


def _condition_codes(values: np.ndarray, levels: np.ndarray) -> np.ndarray:
    """Index of each value in `levels`.

    Exact lookup rather than tolerance-based: both sides come from the same column, so a
    float32 round-trip (0.04 stored as 0.039999999) matches itself.
    """
    idx = np.clip(np.searchsorted(levels, values), 0, len(levels) - 1)
    if not np.array_equal(levels[idx], values):
        raise ValueError("condition values are not exactly present in the level list")
    return idx


def _lifetime_sparseness_chunked(ta: np.ndarray, block: int = 256,
                                 over: str = "conditions") -> np.ndarray:
    """Lifetime sparseness over a (n_conditions, n_trials, n_rois) trial array.

    ``over="conditions"`` averages each condition's trials first, which is the
    Vinje & Gallant definition; ``over="trials"`` uses every individual trial response,
    which runs systematically higher because trial variance adds spread without adding
    mean. The two are not interchangeable -- see docs/comparability.md.

    Accumulated in blocks over the condition axis, because natural movie is 3,600
    conditions x 9 repeats x ~480 ROIs and flattening it would copy 125 MB per plane.
    """
    if over not in ("conditions", "trials"):
        raise ValueError("over must be 'conditions' or 'trials'")
    n_cond, n_trials, n_rois = ta.shape
    n = np.zeros(n_rois)
    s1 = np.zeros(n_rois)
    s2 = np.zeros(n_rois)
    for s in range(0, n_cond, block):
        x = ta[s : s + block]
        if over == "conditions":
            x = _nanmean(x, axis=1)                    # (block, n_rois)
        finite = np.isfinite(x)
        x0 = np.where(finite, x, 0.0)
        axes = (0,) if over == "conditions" else (0, 1)
        n += finite.sum(axis=axes)
        s1 += x0.sum(axis=axes)
        s2 += (x0 * x0).sum(axis=axes)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean, mean_sq = s1 / n, s2 / n
        return (1.0 - mean**2 / mean_sq) / (1.0 - 1.0 / n)
