"""Natural movie: responses to a repeated clip, indexed by frame.

See docs/families/natural_movie.md."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .. import responses as tr
from ..common import _lifetime_sparseness_chunked, _nanmean
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import _reliability_on, absent_frame, roi_frame

def natural_movie_metrics(
    plane,
    trials: pd.DataFrame,
    spont: Sequence[float],
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    rng: Optional[np.random.Generator] = None,
    mouse: Optional[str] = None,
) -> pd.DataFrame:
    """Natural-movie metrics: one row per ROI.

    Each movie frame is a condition and each pass through the movie a repeat.
    ``frac_responsive_trials`` is deterministic (fraction with response > 0, no bootstrap).
    """
    if not len(trials):
        return absent_frame(plane, "natural_movie", mouse)
    rng = np.random.default_rng() if rng is None else rng
    traces = plane.traces[config.trace_type["natural_movie"]]
    window = (0.0, config.nm_response_frames * plane.dt)

    starts = trials["start_time"].to_numpy(dtype=np.float64)
    frames = trials["frame"].to_numpy()
    if np.isnan(frames).any():
        raise ValueError("natural_movie trials contain NaN frame indices")
    frames = frames.astype(int)

    frame_ids = np.arange(frames.max() + 1)
    if not np.array_equal(np.unique(frames), frame_ids):
        raise ValueError(
            "movie frames are not contiguous from 0; the original indexes them "
            "positionally and by label interchangeably, which only works if they are"
        )

    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None)
    n_repeats = int(np.bincount(frames).max())
    ta = tr.trial_array(sweeps, frames, n_trials=n_repeats, n_conditions=len(frame_ids))

    mean_resp = _nanmean(ta, axis=1)                       # (n_frames, n_rois)
    n_rois = plane.n_rois
    roi_ix = np.arange(n_rois)

    all_nan = np.all(~np.isfinite(mean_resp), axis=0)
    safe = np.where(np.isfinite(mean_resp), mean_resp, -np.inf)
    pref_idx = safe.argmax(axis=0)
    pref_response = np.where(all_nan, np.nan, mean_resp[pref_idx, roi_ix])
    pref_img = np.where(all_nan, -1, frame_ids[pref_idx])

    # fraction of repeats with any response at the preferred frame
    pref_trials = ta[pref_idx, :, roi_ix]                  # (n_rois, n_repeats)
    sig = np.where(np.isfinite(pref_trials), (pref_trials > 0).astype(float), np.nan)
    frac_responsive = _nanmean(sig, axis=1)

    null = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, None,
        n_boot=config.other_n_boot, n_means=n_repeats, rng=rng,
        memory_budget_mb=config.memory_budget_mb,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        z_score = (pref_response - null.mean(axis=1)) / null.std(axis=1)

    out = roi_frame(plane, mouse=mouse)
    out["frac_responsive_trials"] = frac_responsive
    out["lifetime_sparseness"] = _lifetime_sparseness_chunked(
        ta, over=config.lifetime_sparseness_over)
    out["reliability"] = tr.trial_reliability(ta)
    out["reliability_dff"] = _reliability_on(
        plane, "dff", starts, frames, n_trials=n_repeats,
        n_conditions=len(frame_ids), window=window)
    out["n_trials_at_pref"] = np.where(
        all_nan, np.nan, np.isfinite(pref_trials).sum(axis=1)).astype(float)
    out["pref_img"] = pref_img
    out["pref_response"] = pref_response
    out["z_score"] = z_score
    return out
