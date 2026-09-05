"""Stimulus-locked responses, trial arrays, and bootstrap nulls -- the arithmetic layer.

Arrays in, arrays out: nothing here imports pandas or an NWB reader, so it can be tested
against a synthetic trace with no data mounted. See docs/pipeline.md for the prefix-sum
design and for why two different window conventions exist.
"""

from typing import Dict, Optional, Tuple

import numpy as np

__all__ = [
    "prefix_sums",
    "window_bounds",
    "window_means",
    "sweep_responses",
    "sweep_responses_frames",
    "spontaneous_null",
    "trial_array",
    "frac_trials_above_null",
    "lifetime_sparseness",
    "trial_reliability",
    "si_permutation_test",
]


# ------------------------------------------------------------------ the primitive


def prefix_sums(traces: np.ndarray) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Cumulative sums over time, for constant-cost window means.

    ``traces`` is (n_frames, n_rois), NWB-Zarr's native orientation -- do not transpose
    first. Returns ``cs``, where ``cs[k] == traces[:k].sum(axis=0)`` in float64, so the mean
    over ``traces[a:b]`` is ``(cs[b] - cs[a]) / (b - a)``; and ``counts``, the cumulative
    count of finite samples, or None when the traces contain no NaN.
    """
    traces = np.asarray(traces)
    if traces.ndim != 2:
        raise ValueError(f"expected (n_frames, n_rois), got shape {traces.shape}")

    n_frames, n_rois = traces.shape
    finite = np.isfinite(traces)
    counts = None
    if not finite.all():
        counts = np.zeros((n_frames + 1, n_rois), dtype=np.int64)
        np.cumsum(finite, axis=0, out=counts[1:])
        traces = np.where(finite, traces, 0.0)

    cs = np.zeros((n_frames + 1, n_rois), dtype=np.float64)
    np.cumsum(traces, axis=0, dtype=np.float64, out=cs[1:])
    return cs, counts


def window_bounds(
    timestamps: np.ndarray, starts: np.ndarray, t0: float, t1: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Half-open index bounds for the label-closed interval [start + t0, start + t1].

    Matches ``xarray.sel(time=slice(lo, hi))``, which includes both endpoints, so the
    returned windows have variable width.
    """
    timestamps = np.asarray(timestamps, dtype=np.float64)
    starts = np.asarray(starts, dtype=np.float64)
    a = np.searchsorted(timestamps, starts + t0, side="left")
    b = np.searchsorted(timestamps, starts + t1, side="right")
    return a, b


def window_means(
    cs: np.ndarray, counts: Optional[np.ndarray], a: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """Mean over `traces[a[i]:b[i]]` for each window i. Returns (n_windows, n_rois).

    An empty window (b <= a) yields NaN via 0/0, matching `xarray.mean` of an empty
    selection. That happens at session edges and must not raise.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    total = cs[b] - cs[a]
    n = (b - a)[:, None].astype(np.float64) if counts is None else (counts[b] - counts[a])
    with np.errstate(invalid="ignore", divide="ignore"):
        return total / n


# ------------------------------------------------------------------ sweep responses


def sweep_responses(
    traces: np.ndarray,
    timestamps: np.ndarray,
    starts: np.ndarray,
    response_window: Tuple[float, float],
    baseline_window: Optional[Tuple[float, float]] = None,
    block: int = 8192,
) -> np.ndarray:
    """Mean response per stimulus sweep. Returns (n_sweeps, n_rois).

    ``baseline_window=None`` means no baseline subtraction at all -- the shipped
    configuration for every events-based family, and not the same as subtracting zero.
    ``block`` bounds the (block, n_rois) intermediate.
    """
    cs, counts = prefix_sums(traces)
    starts = np.asarray(starts, dtype=np.float64)
    out = np.empty((starts.size, traces.shape[1]), dtype=np.float64)

    for s in range(0, starts.size, block):
        sl = slice(s, min(s + block, starts.size))
        a, b = window_bounds(timestamps, starts[sl], *response_window)
        r = window_means(cs, counts, a, b)
        if baseline_window is not None:
            a0, b0 = window_bounds(timestamps, starts[sl], *baseline_window)
            r = r - window_means(cs, counts, a0, b0)
        out[sl] = r

    return out


# ------------------------------------------------------------------ bootstrap null


def sweep_responses_frames(
    traces: np.ndarray,
    timestamps: np.ndarray,
    starts: np.ndarray,
    n_frames: int,
    offset_frames: int = 0,
) -> np.ndarray:
    """Mean of exactly ``n_frames`` samples from each onset. Returns (n_sweeps, n_rois).

    The fixed-width alternative to ``sweep_responses``, starting at the first sample at or
    after the onset. The two are not interchangeable; see docs/pipeline.md.
    """
    if n_frames < 1:
        raise ValueError("n_frames must be at least 1")
    cs, counts = prefix_sums(traces)
    timestamps = np.asarray(timestamps, dtype=np.float64)
    starts = np.asarray(starts, dtype=np.float64)

    a = np.searchsorted(timestamps, starts, side="left") + offset_frames
    a = np.clip(a, 0, len(timestamps))
    b = np.clip(a + n_frames, 0, len(timestamps))
    return window_means(cs, counts, a, b)


def _nearest_index(timestamps: np.ndarray, value: float) -> int:
    """`pandas.Index.get_loc(value, method="nearest")`, which pandas 2.0 removed.

    Ties break toward the earlier sample, matching pandas.
    """
    n = len(timestamps)
    i = int(np.clip(np.searchsorted(timestamps, value), 1, n - 1))
    return i - int(abs(value - timestamps[i - 1]) <= abs(timestamps[i] - value))


def spontaneous_null(
    traces: np.ndarray,
    timestamps: np.ndarray,
    spont_start: float,
    spont_stop: float,
    response_window: Tuple[float, float],
    baseline_window: Optional[Tuple[float, float]] = None,
    n_boot: int = 2500,
    n_means: int = 1,
    rng: Optional[np.random.Generator] = None,
    memory_budget_mb: float = 64.0,
) -> np.ndarray:
    """Bootstrap distribution of spontaneous responses. Returns (n_rois, n_boot).

    Windows are frame-indexed and fixed width, onsets are drawn uniformly with replacement
    from one spontaneous block, and ``n_means`` draws are averaged per bootstrap sample:
    1 gives the single-trial null behind ``frac_responsive_trials``, ``n_trials`` the
    multi-trial null behind ``z_score``. ``memory_budget_mb`` bounds the block size -- the
    loop is over memory blocks, not bootstraps.
    """
    rng = np.random.default_rng() if rng is None else rng
    timestamps = np.asarray(timestamps, dtype=np.float64)
    dt = float(np.median(np.diff(timestamps)))

    r0, r1 = (int(round(response_window[0] / dt)), int(round(response_window[1] / dt)))
    if baseline_window is None:
        b0 = b1 = None
        start_pad = 0
    else:
        b0, b1 = (int(round(baseline_window[0] / dt)), int(round(baseline_window[1] / dt)))
        start_pad = -b0

    lo = _nearest_index(timestamps, spont_start) + start_pad
    hi = _nearest_index(timestamps, spont_stop) - r1  # exclusive, as np.random.randint
    if hi <= lo:
        raise ValueError(
            f"spontaneous block [{spont_start:.1f}, {spont_stop:.1f}] s is too short for "
            f"response window {response_window} at dt={dt:.4f} s"
        )

    cs, counts = prefix_sums(traces)
    n_rois = traces.shape[1]
    idx = rng.integers(lo, hi, size=(n_boot, n_means))

    per_row_bytes = max(n_means * n_rois * 8, 1)
    block = max(1, int(memory_budget_mb * 1e6 // per_row_bytes))

    out = np.empty((n_boot, n_rois), dtype=np.float64)
    for s in range(0, n_boot, block):
        j = idx[s : s + block]
        f = j.ravel()
        r = window_means(cs, counts, f + r0, f + r1)
        if b0 is not None:
            r = r - window_means(cs, counts, f + b0, f + b1)
        out[s : s + block] = r.reshape(j.shape[0], n_means, n_rois).mean(axis=1)

    return out.T


# ------------------------------------------------------------------ trial arrays


def trial_array(
    sweep_resp: np.ndarray,
    condition: np.ndarray,
    n_trials: Optional[int] = None,
    n_conditions: Optional[int] = None,
) -> np.ndarray:
    """Scatter per-sweep responses into (n_conditions, n_trials, n_rois), NaN-padded.

    ``condition`` is an integer code per sweep in table order, which is chronological, so a
    sweep's trial index is its rank within its condition -- the same indexing the per-trial
    running speeds use. Conditions presented fewer than ``n_trials`` times keep NaN in the
    tail; sweeps beyond ``n_trials`` for a condition are dropped.
    """
    sweep_resp = np.asarray(sweep_resp)
    condition = np.asarray(condition)
    if condition.size != sweep_resp.shape[0]:
        raise ValueError(
            f"condition has {condition.size} entries but sweep_resp has "
            f"{sweep_resp.shape[0]} sweeps"
        )
    if condition.size and condition.min() < 0:
        raise ValueError("condition codes must be non-negative")

    # stable sort keeps time order inside each condition group
    order = np.argsort(condition, kind="stable")
    codes = condition[order]

    if n_conditions is None:
        n_conditions = int(codes.max()) + 1 if codes.size else 0
    first = np.searchsorted(codes, np.arange(n_conditions), side="left")
    trial = np.arange(codes.size) - first[codes]

    if n_trials is None:
        n_trials = int(trial.max()) + 1 if trial.size else 0

    out = np.full((n_conditions, n_trials, sweep_resp.shape[1]), np.nan, dtype=np.float64)
    keep = trial < n_trials
    out[codes[keep], trial[keep]] = sweep_resp[order[keep]]
    return out


# ------------------------------------------------------------------ reductions


def frac_trials_above_null(
    trial_resp: np.ndarray,
    null_single: np.ndarray,
    p_thresh: float = 0.05,
    block_rois: int = 64,
) -> np.ndarray:
    """Fraction of a neuron's trials that beat its spontaneous null. Returns (n_rois,).

    ``trial_resp`` is (n_rois, n_trials) at each neuron's preferred condition, NaN-padded;
    ``null_single`` is (n_rois, n_boot). A trial counts as responsive when
    ``mean(null > response) < p_thresh``, which is not the same as comparing against a
    quantile of the null. NaN-padded trials are excluded rather than counted.
    """
    trial_resp = np.asarray(trial_resp, dtype=np.float64)
    null_single = np.asarray(null_single, dtype=np.float64)
    n_rois = trial_resp.shape[0]
    out = np.full(n_rois, np.nan)

    for s in range(0, n_rois, block_rois):
        sl = slice(s, min(s + block_rois, n_rois))
        r = trial_resp[sl]                                    # (k, n_trials)
        p = (null_single[sl][:, None, :] > r[:, :, None]).mean(axis=2)
        sig = np.where(np.isnan(r), np.nan, (p < p_thresh).astype(np.float64))
        with np.errstate(invalid="ignore"):
            out[sl] = np.nanmean(sig, axis=1) if sig.shape[1] else np.nan

    return out


def lifetime_sparseness(x: np.ndarray) -> np.ndarray:
    """Olsen & Wilson (2008) lifetime sparseness. `x` is (n_rois, n_responses).

    Computed over **every individual trial response**, flattened across conditions and
    trials — not over the condition means. NaNs are dropped first, so the normaliser
    `1 - 1/n` uses each neuron's own count of real responses.
    """
    x = np.asarray(x, dtype=np.float64)
    finite = np.isfinite(x)
    n = finite.sum(axis=1).astype(np.float64)
    xs = np.where(finite, x, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = xs.sum(axis=1) / n
        mean_sq = (xs * xs).sum(axis=1) / n
        return (1.0 - mean**2 / mean_sq) / (1.0 - 1.0 / n)


def _vector_strength(tuning: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """|sum R_theta e^{i k theta}| / sum R_theta, matching `_selectivity_index`.

    A plain (NaN-propagating) sum in the denominator, deliberately: the original used
    `np.sum`, and switching to `np.nansum` here would silently change which neurons get
    a finite selectivity index.
    """
    norm = np.sum(tuning, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(tuning @ phase / np.where(norm == 0, np.nan, norm))


def si_permutation_test(
    trial_responses: np.ndarray,
    n_shuffles: int = 1000,
    rng: Optional[np.random.Generator] = None,
    block: int = 64,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Permutation test for orientation and direction selectivity.

    ``trial_responses`` is (n_rois, n_directions, n_trials) at each neuron's preferred
    spatial frequency. Returns ``{"osi": (si, p), "dsi": (si, p)}``, both reusing one set of
    shuffled tunings. This index is a circular vector strength on index-derived angles, a
    different quantity from the ratio-form ``osi``/``dsi`` in the output table.
    """
    rng = np.random.default_rng() if rng is None else rng
    trial_responses = np.asarray(trial_responses, dtype=np.float64)
    n_rois, n_dir, n_trial = trial_responses.shape

    ang = np.arange(n_dir) / n_dir * 2 * np.pi
    phase = {"dsi": np.exp(1j * ang), "osi": np.exp(2j * ang)}

    with np.errstate(invalid="ignore", divide="ignore"):
        tune_true = _nanmean_quiet(trial_responses, axis=-1)
        true = {m: _vector_strength(tune_true, p) for m, p in phase.items()}
        exceed = {m: np.zeros(n_rois) for m in phase}

        for s in range(0, n_shuffles, block):
            k = min(block, n_shuffles - s)
            xs = np.broadcast_to(trial_responses, (k, n_rois, n_dir, n_trial)).copy()
            xs = rng.permuted(xs, axis=2)          # permute directions per (roi, trial)
            ts = _nanmean_quiet(xs, axis=-1)
            for m, p in phase.items():
                exceed[m] += (true[m] < _vector_strength(ts, p)).sum(axis=0)

    return {m: (true[m], exceed[m] / n_shuffles) for m in phase}


def _nanmean_quiet(x: np.ndarray, axis: int) -> np.ndarray:
    """`np.nanmean` without the all-NaN-slice RuntimeWarning, which fires constantly on
    NaN-padded trial arrays and drowns real warnings."""
    finite = np.isfinite(x)
    n = finite.sum(axis=axis)
    total = np.where(finite, x, 0.0).sum(axis=axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, total / np.maximum(n, 1), np.nan)


def trial_reliability(ta: np.ndarray, min_conditions: int = 3) -> np.ndarray:
    """Mean pairwise between-trial correlation per ROI. Returns (n_rois,), NaN if no pair.

    ``ta`` is (n_conditions, n_trials, n_rois) from ``trial_array``. Pairs are scored
    pairwise-complete, over the conditions where both trials are finite; zero-variance
    trials are skipped rather than scored zero, and ``min_conditions`` sets the smallest
    overlap a pair may have.

    Not the same quantity as de Vries et al. 2019's "reliability", which is our
    ``frac_responsive_trials``. See docs/comparability.md.
    """
    ta = np.asarray(ta, dtype=np.float64)
    if ta.ndim != 3:
        raise ValueError(f"expected (n_conditions, n_trials, n_rois), got {ta.shape}")
    n_cond, n_trials, n_rois = ta.shape
    total = np.zeros(n_rois, dtype=np.float64)
    count = np.zeros(n_rois, dtype=np.int64)

    for i in range(n_trials):
        a_all = ta[:, i, :]
        for j in range(i + 1, n_trials):
            b_all = ta[:, j, :]
            both = np.isfinite(a_all) & np.isfinite(b_all)      # (n_cond, n_rois)
            n = both.sum(axis=0)
            ok = n >= min_conditions
            if not ok.any():
                continue
            a = np.where(both, a_all, 0.0)
            b = np.where(both, b_all, 0.0)
            n_safe = np.maximum(n, 1)
            ma, mb = a.sum(axis=0) / n_safe, b.sum(axis=0) / n_safe
            da = np.where(both, a_all - ma, 0.0)
            db = np.where(both, b_all - mb, 0.0)
            va, vb = (da * da).sum(axis=0), (db * db).sum(axis=0)
            # a flat trial has no correlation to give -- skip the pair for that ROI
            usable = ok & (va > 0) & (vb > 0)
            if not usable.any():
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                r = (da * db).sum(axis=0) / np.sqrt(np.where(usable, va * vb, 1.0))
            total += np.where(usable, r, 0.0)
            count += usable
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(count > 0, total / np.maximum(count, 1), np.nan)
