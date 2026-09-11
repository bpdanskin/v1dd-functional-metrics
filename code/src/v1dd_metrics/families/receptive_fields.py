"""Receptive fields from locally sparse noise: ON/OFF subfield maps and centres.

See docs/families/receptive_fields.md."""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .. import responses as tr
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import absent_frame, roi_frame

def _rf_pixel_to_degrees(mean_idx, centers: np.ndarray, scale_bug: bool) -> np.ndarray:
    """Map a fractional pixel index to degrees of visual angle.

    ``scale_bug=True`` reproduces the historical compressed scale (range/n instead of
    range/(n-1)); ``False`` gives the correct interpolation.
    """
    mean_idx = np.asarray(mean_idx, dtype=np.float64)
    if scale_bug:
        pitch = (centers[-1] - centers[0]) / len(centers)      # should be len - 1
        return (mean_idx + 0.5) * pitch + centers[0]
    pitch = (centers[-1] - centers[0]) / (len(centers) - 1)
    return mean_idx * pitch + centers[0]


def receptive_field_metrics(
    plane,
    trials: pd.DataFrame,
    spont: Sequence[float],
    lsn: Mapping[str, Any],
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    rng: Optional[np.random.Generator] = None,
    mouse: Optional[str] = None,
) -> pd.DataFrame:
    """ON/OFF subfield maps and centres from the locally-sparse-noise stimulus.

    Trace type and response window are set by ``config``. ``lsn`` is the dict from
    ``load_lsn_template``; its ``pixel_on``/``pixel_off`` determine the design matrix.
    Returns ``(metrics, rf_map)`` where ``rf_map`` is the continuous pre-threshold map
    ``(n_rois, 2, n_rows, n_cols)`` in float32.
    """
    images = np.asarray(lsn["images"])
    n_rows, n_cols = images.shape[1], images.shape[2]

    if not len(trials):
        empty_map = np.zeros((plane.n_rois, 2, n_rows, n_cols), dtype=np.float32)
        return absent_frame(plane, "rf_metrics", mouse), empty_map
    rng = np.random.default_rng() if rng is None else rng
    trace_key = config.trace_type["locally_sparse_noise"]
    traces = plane.traces[trace_key]
    window = (0.0, config.lsn_response_frames * plane.dt)
    # Events are non-negative and already denoised by L0's sparsity penalty;
    # baseline subtraction would make them signed, defeating the advantage.
    baseline = None if trace_key == "events" else (-1.0, 0.0)
    pixel_on, pixel_off = lsn.get("pixel_on"), lsn.get("pixel_off")
    if pixel_on is None or pixel_off is None:
        raise ValueError(
            f"could not determine ON/OFF pixel codes from the template "
            f"(values seen: {lsn.get('pixel_values')})"
        )

    starts = trials["start_time"].to_numpy(dtype=np.float64)
    frames = trials["frame"].to_numpy()
    if np.isnan(frames).any():
        raise ValueError("locally_sparse_noise trials contain NaN frame indices")
    frames = frames.astype(int)
    if frames.max() >= len(images):
        raise ValueError(
            f"frame index {frames.max()} exceeds the {len(images)}-frame template"
        )

    n_pixels = n_rows * n_cols

    # (2 * n_pixels, n_sweeps): rows 0..n_pixels-1 are ON, the rest OFF. Gray is neither.
    stim_pixels = images[frames].reshape(len(frames), n_pixels)
    design = np.concatenate([stim_pixels == pixel_on, stim_pixels == pixel_off], axis=1).T

    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, baseline)

    null = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, baseline,
        n_boot=config.other_n_boot, n_means=1, rng=rng,
        memory_budget_mb=config.memory_budget_mb,
    )
    threshold = np.quantile(null, 0.95, axis=1)                 # (n_rois,)
    significant = sweeps > threshold[None, :]                   # (n_sweeps, n_rois)

    n_pixel_trials = design.sum(axis=1)                         # (2 * n_pixels,)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = (design.astype(np.int64) @ significant).T / np.where(
            n_pixel_trials > 0, n_pixel_trials, np.nan)
    frac = np.nan_to_num(frac, nan=0.0)
    frac[~plane.is_valid] = 0.0
    rf_map = frac.reshape(plane.n_rois, 2, n_rows, n_cols).astype(np.float32).copy()
    frac[frac < config.rf_frac_thresh] = 0.0
    rf = frac.reshape(plane.n_rois, 2, n_rows, n_cols)          # dim 1: 0 = ON, 1 = OFF

    mask = rf > 0
    counts = mask.sum(axis=(2, 3))                              # (n_rois, 2)
    rows_ix = np.arange(n_rows)[None, None, :, None]
    cols_ix = np.arange(n_cols)[None, None, None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_row = (mask * rows_ix).sum(axis=(2, 3)) / counts
        mean_col = (mask * cols_ix).sum(axis=(2, 3)) / counts

    altitudes = np.asarray(lsn["altitudes"], dtype=np.float64)
    azimuths = np.asarray(lsn["azimuths"], dtype=np.float64)
    alt = _rf_pixel_to_degrees(mean_row, altitudes, config.rf_center_scale_bug)
    azi = _rf_pixel_to_degrees(mean_col, azimuths, config.rf_center_scale_bug)

    has_on = counts[:, 0] > 0
    has_off = counts[:, 1] > 0

    out = roi_frame(plane, mouse=mouse)
    out["has_rf_on"] = has_on
    out["has_rf_off"] = has_off
    out["has_rf_on_or_off"] = has_on | has_off
    out["azimuth_rf_on"] = np.where(has_on, azi[:, 0], np.nan)
    out["altitude_rf_on"] = np.where(has_on, alt[:, 0], np.nan)
    out["azimuth_rf_off"] = np.where(has_off, azi[:, 1], np.nan)
    out["altitude_rf_off"] = np.where(has_off, alt[:, 1], np.nan)
    return out, rf_map
