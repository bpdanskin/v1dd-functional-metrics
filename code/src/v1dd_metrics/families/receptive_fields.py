"""Receptive fields from locally sparse noise: ON/OFF subfield maps and centres.

Two methods, selected by ``config.rf_method``:

* ``"greedy"`` (default) — per-pixel STA of events, per-pixel bootstrap p-value
  (resample sweeps with replacement), Holm-Šidák correction across all 2*n_pixels
  tests, detect where corrected p < alpha. Two variants ship: the strict one at
  ``rf_greedy_alpha_strict`` (canonical columns) and a sensitive one at
  ``rf_greedy_alpha_sens`` (``_a05`` columns).
* ``"fraction"`` — the historical method: fraction of a pixel's presentations whose
  response beats the ROI's spontaneous 95th percentile, thresholded at
  ``rf_frac_thresh``. Used by ``REFERENCE_CONFIG``.

See docs/families/receptive_fields.md and docs/explorations/greedy_rf.md."""

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
        pitch = (centers[-1] - centers[0]) / len(centers)
        return (mean_idx + 0.5) * pitch + centers[0]
    pitch = (centers[-1] - centers[0]) / (len(centers) - 1)
    return mean_idx * pitch + centers[0]


def holm_sidak_reject(pvals: np.ndarray, alpha: float) -> np.ndarray:
    """Step-down Holm-Šidák rejection over axis 0 (tests) for each column (ROI).

    ``pvals`` is ``(n_tests, n_rois)``; returns a boolean reject mask of the same shape.
    Matches ``statsmodels.multipletests(..., method='hs')[0]`` without the dependency.
    """
    m = pvals.shape[0]
    order = np.argsort(pvals, axis=0)
    sorted_p = np.take_along_axis(pvals, order, axis=0)
    thr = 1.0 - (1.0 - alpha) ** (1.0 / (m - np.arange(m)))
    passed = sorted_p <= thr[:, None]
    rej_sorted = np.logical_and.accumulate(passed, axis=0)
    mask = np.zeros_like(pvals, dtype=bool)
    np.put_along_axis(mask, order, rej_sorted, axis=0)
    return mask


def _has_and_centres(out, suffix, mask, altitudes, azimuths, scale_bug):
    """Write the 7 RF columns for one variant from a boolean detection mask.

    ``mask`` is ``(n_rois, 2, n_rows, n_cols)`` (dim 1: 0=ON, 1=OFF). Centres are the
    centroid of the significant pixels, mapped to degrees.
    """
    n_rows, n_cols = mask.shape[2], mask.shape[3]
    counts = mask.sum(axis=(2, 3))
    rows_ix = np.arange(n_rows)[None, None, :, None]
    cols_ix = np.arange(n_cols)[None, None, None, :]
    denom = np.where(counts > 0, counts, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_row = (mask * rows_ix).sum(axis=(2, 3)) / denom
        mean_col = (mask * cols_ix).sum(axis=(2, 3)) / denom
    alt = _rf_pixel_to_degrees(mean_row, altitudes, scale_bug)
    azi = _rf_pixel_to_degrees(mean_col, azimuths, scale_bug)
    pix_deg2 = abs((altitudes[1] - altitudes[0]) * (azimuths[1] - azimuths[0]))
    area = counts * pix_deg2
    has_on, has_off = counts[:, 0] > 0, counts[:, 1] > 0
    out[f"has_rf_on{suffix}"] = has_on
    out[f"has_rf_off{suffix}"] = has_off
    out[f"has_rf_on_or_off{suffix}"] = has_on | has_off
    out[f"azimuth_rf_on{suffix}"] = np.where(has_on, azi[:, 0], np.nan)
    out[f"altitude_rf_on{suffix}"] = np.where(has_on, alt[:, 0], np.nan)
    out[f"azimuth_rf_off{suffix}"] = np.where(has_off, azi[:, 1], np.nan)
    out[f"altitude_rf_off{suffix}"] = np.where(has_off, alt[:, 1], np.nan)
    out[f"rf_on_area{suffix}"] = np.where(has_on, area[:, 0], np.nan)
    out[f"rf_off_area{suffix}"] = np.where(has_off, area[:, 1], np.nan)


def _design_matrix(images, frames, pixel_on, pixel_off):
    """(2*n_pixels, n_sweeps) design: rows 0..n_pixels-1 ON, the rest OFF; gray neither."""
    n_pixels = images.shape[1] * images.shape[2]
    stim = images[frames].reshape(len(frames), n_pixels)
    return np.concatenate([stim == pixel_on, stim == pixel_off], axis=1).T


def receptive_field_metrics(
    plane,
    trials: pd.DataFrame,
    spont: Sequence[float],
    lsn: Mapping[str, Any],
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    rng: Optional[np.random.Generator] = None,
    mouse: Optional[str] = None,
):
    """ON/OFF subfield detection and centres from the locally-sparse-noise stimulus.

    Returns ``(metrics, arrays)``. ``arrays`` is a dict keyed by ``method``:
    ``"greedy"`` carries ``sta``, ``ge`` (bootstrap shuffle counts, uint16) and
    ``strict_mask``; ``"fraction"`` carries the pre-threshold ``rf_map``. Both also
    carry ``overlap_map``/``overlap_thresh`` for ``window_containment``. All arrays are
    ``(n_rois, 2, n_rows, n_cols)``.
    """
    images = np.asarray(lsn["images"])
    n_rows, n_cols = images.shape[1], images.shape[2]
    method = config.rf_method

    if not len(trials):
        z = np.zeros((plane.n_rois, 2, n_rows, n_cols), dtype=np.float32)
        arrays = ({"method": "fraction", "rf_map": z, "overlap_map": z,
                   "overlap_thresh": config.rf_frac_thresh} if method == "fraction"
                  else {"method": "greedy", "sta": z,
                        "ge": z.astype(np.uint16), "strict_mask": z.astype(bool),
                        "overlap_map": z, "overlap_thresh": 0.5})
        return absent_frame(plane, "rf_metrics", mouse), arrays

    rng = np.random.default_rng() if rng is None else rng
    trace_key = config.trace_type["locally_sparse_noise"]
    traces = plane.traces[trace_key]
    window = (0.0, config.lsn_response_frames * plane.dt)
    pixel_on, pixel_off = lsn.get("pixel_on"), lsn.get("pixel_off")
    if pixel_on is None or pixel_off is None:
        raise ValueError(
            f"could not determine ON/OFF pixel codes from the template "
            f"(values seen: {lsn.get('pixel_values')})")

    starts = trials["start_time"].to_numpy(dtype=np.float64)
    frames = trials["frame"].to_numpy()
    if np.isnan(frames).any():
        raise ValueError("locally_sparse_noise trials contain NaN frame indices")
    frames = frames.astype(int)
    if frames.max() >= len(images):
        raise ValueError(
            f"frame index {frames.max()} exceeds the {len(images)}-frame template")

    design = _design_matrix(images, frames, pixel_on, pixel_off)
    altitudes = np.asarray(lsn["altitudes"], dtype=np.float64)
    azimuths = np.asarray(lsn["azimuths"], dtype=np.float64)
    out = roi_frame(plane, mouse=mouse)
    ctx = (out, plane, traces, design, starts, window, altitudes, azimuths,
           n_rows, n_cols, config, rng)

    if method == "fraction":
        return _fraction_rf(*ctx, spont=spont)
    if method == "greedy":
        return _greedy_rf(*ctx)
    raise ValueError(f"unknown rf_method {method!r} (expected 'greedy' or 'fraction')")


def _fraction_rf(out, plane, traces, design, starts, window, altitudes, azimuths,
                 n_rows, n_cols, config, rng, *, spont):
    """Historical method: fraction of presentations above the spontaneous 95th pct."""
    trace_key = config.trace_type["locally_sparse_noise"]
    baseline = None if trace_key == "events" else (-1.0, 0.0)
    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, baseline)
    null = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, baseline,
        n_boot=config.other_n_boot, n_means=1, rng=rng,
        memory_budget_mb=config.memory_budget_mb)
    threshold = np.quantile(null, 0.95, axis=1)
    significant = sweeps > threshold[None, :]

    n_pixel_trials = design.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = (design.astype(np.int64) @ significant).T / np.where(
            n_pixel_trials > 0, n_pixel_trials, np.nan)
    frac = np.nan_to_num(frac, nan=0.0)
    frac[~plane.is_valid] = 0.0
    rf_map = frac.reshape(plane.n_rois, 2, n_rows, n_cols).astype(np.float32).copy()

    mask = (frac.reshape(plane.n_rois, 2, n_rows, n_cols) >= config.rf_frac_thresh)
    _has_and_centres(out, "", mask, altitudes, azimuths, config.rf_center_scale_bug)
    _has_and_centres(out, "_a05", np.zeros_like(mask), altitudes, azimuths,
                     config.rf_center_scale_bug)
    arrays = {"method": "fraction", "rf_map": rf_map,
              "overlap_map": rf_map, "overlap_thresh": config.rf_frac_thresh}
    return out, arrays


def _greedy_rf(out, plane, traces, design, starts, window, altitudes, azimuths,
               n_rows, n_cols, config, rng):
    """Greedy pixelwise RF: STA + per-pixel bootstrap p-values + Holm-Šidák, two variants."""
    A = design.astype(np.float32)
    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None
                                ).astype(np.float32)
    n_sweeps = sweeps.shape[0]
    sta = A @ sweeps

    ge = np.zeros_like(sta)
    for _ in range(config.rf_greedy_n_boot):
        idx = rng.integers(0, n_sweeps, n_sweeps)
        ge += (A @ sweeps[idx]) >= sta
    pvals = ge / config.rf_greedy_n_boot

    def variant_mask(alpha):
        m = holm_sidak_reject(pvals, alpha).T.reshape(plane.n_rois, 2, n_rows, n_cols)
        m[~plane.is_valid] = False
        return m

    strict = variant_mask(config.rf_greedy_alpha_strict)
    sens = variant_mask(config.rf_greedy_alpha_sens)
    _has_and_centres(out, "", strict, altitudes, azimuths, config.rf_center_scale_bug)
    _has_and_centres(out, "_a05", sens, altitudes, azimuths, config.rf_center_scale_bug)

    sta_maps = sta.T.reshape(plane.n_rois, 2, n_rows, n_cols).astype(np.float32)
    arrays = {
        "method": "greedy",
        "sta": sta_maps,
        "ge": np.clip(ge.T.reshape(plane.n_rois, 2, n_rows, n_cols), 0, 65535).astype(np.uint16),
        "strict_mask": strict,
        "overlap_map": strict.astype(np.float32),
        "overlap_thresh": 0.5,
    }
    return out, arrays
