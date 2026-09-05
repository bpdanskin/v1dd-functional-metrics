"""The output schema: per-ROI identity, column order, and empty-family frames."""

from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from . import responses as tr

OUTPUT_COLUMNS: Dict[str, Sequence[str]] = {
    "drifting_gratings_full": [
        "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
        "pika_roi_confidence",
        "dsi", "frac_responsive_trials", "gosi", "is_responsive",
        "lifetime_sparseness", "osi", "preferred_dir", "preferred_sf", "pref_dir_mean"],
    "natural_images": [
        "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
        "pika_roi_confidence",
        "frac_responsive_trials", "lifetime_sparseness", "pref_img", "pref_response",
        "z_score", "reliability", "reliability_dff", "n_trials_at_pref"],
    "surround_suppression": [
        "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
        "pika_roi_confidence",
        "ssi", "ssi_avg", "ssi_avg_at_pref_sf", "ssi_running",
        "ssi_running_avg_at_pref_sf", "ssi_stationary",
        "ssi_stationary_avg_at_pref_sf", "ssi_tuning_fit",
        "dgw_center_azimuth", "dgw_center_elevation", "dgw_center_inferred",
        "dgw_rf_distance_on", "dgw_rf_distance_off",
        "dgw_rf_overlap_on", "dgw_rf_overlap_off"],
    "roi_summary": [
        "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
        "pika_roi_confidence",
        "snr", "signal_power", "noise_power",
        "run_frac", "spont_run_frac",
        "spont_rate", "spont_rate_run", "spont_rate_stat",
        "run_mod_dgf", "run_mod_dgw", "run_mod_spont",
        "run_corr_dff"],
    "rf_metrics": [
        "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
        "pika_roi_confidence",
        "has_rf_on", "has_rf_off", "has_rf_on_or_off",
        "azimuth_rf_on", "altitude_rf_on", "azimuth_rf_off", "altitude_rf_off"],
}


# Families that publish an identical column set share one entry.
OUTPUT_COLUMNS["drifting_gratings_windowed"] = OUTPUT_COLUMNS["drifting_gratings_full"]
OUTPUT_COLUMNS["natural_images_12"] = OUTPUT_COLUMNS["natural_images"]
OUTPUT_COLUMNS["natural_movie"] = OUTPUT_COLUMNS["natural_images"]

OUTPUT_COLUMNS["roi_position"] = [
    "roi_unique_id", "mouse", "column", "volume", "plane", "roi", "depth_um",
    "pika_roi_confidence",
    "roi_x_px", "roi_y_px", "roi_area_px", "roi_radius_px",
    "roi_x_um", "roi_y_um",
    "roi_x_um_published", "roi_y_um_published",
    "roi_x_um_retinotopic", "roi_y_um_retinotopic"]


def roi_frame(plane, mouse: Optional[str] = None) -> pd.DataFrame:
    """The identity columns every output table starts with.

    `roi_unique_id` reproduces the published format `M{mouse}_{volume}_{plane}_{roi}`,
    which **omits the column** and therefore collides across the five columns — 164,345
    published rows share only 56,449 distinct ids. It is emitted for drop-in
    compatibility; `roi_key` is the non-colliding version. **Join on
    `(column, volume, plane, roi)`, never on either string.**
    """
    n = plane.n_rois
    # The mouse comes from the file, not from a constant: this pipeline is expected to
    # run on other animals. `mouse` overrides only when a caller genuinely knows better.
    mouse_num = (mouse or "").lstrip("M") or getattr(plane, "mouse_id", "")
    if not mouse_num:
        raise ValueError(
            "no mouse id: pass mouse=, or load the plane with load_plane(), which reads "
            "it from nwb.subject via session_mouse()"
        )
    mouse = f"M{mouse_num}"
    return pd.DataFrame({
        "roi_unique_id": [f"M{mouse_num}_{plane.volume}_{plane.plane}_{r}" for r in plane.roi],
        "roi_key": [
            f"M{mouse_num}_{plane.column}_{plane.volume}_{plane.plane}_{r}" for r in plane.roi
        ],
        "mouse": [mouse] * n,
        "column": np.full(n, plane.column, dtype=int),
        "volume": [plane.volume] * n,
        "plane": np.full(n, plane.plane, dtype=int),
        "roi": plane.roi.astype(int),
        # Physical depth, which (column, volume, plane) only encodes implicitly. NaN when
        # the file does not carry it -- no metric depends on it.
        "depth_um": np.full(n, getattr(plane, "depth_um", None)
                            if getattr(plane, "depth_um", None) is not None else np.nan,
                            dtype=float),
        # Segmentation confidence, emitted so consumers can see which ROIs the pipeline
        # treated as unreliable. Without it, low-confidence ROIs are neither dropped nor
        # labelled: `preferred_dir`, the `ssi*` columns and every receptive-field column
        # are suppressed for them, but `osi`, `dsi`, `lifetime_sparseness` and the natural
        # scene metrics are populated as usual, so they enter any population average
        # unnoticed. `is_valid` is this column > 0.5.
        "pika_roi_confidence": _roi_confidence(plane),
    })


def _roi_confidence(plane) -> np.ndarray:
    """Per-ROI segmentation confidence, or NaN where the ROI table does not carry it."""
    table = getattr(plane, "roi_table", None)
    if table is None or "pika_roi_confidence" not in getattr(table, "columns", []):
        return np.full(plane.n_rois, np.nan, dtype=float)
    return pd.to_numeric(table["pika_roi_confidence"], errors="coerce").to_numpy(float)


BOOLEAN_COLUMNS = frozenset({"has_rf_on", "has_rf_off", "has_rf_on_or_off",
                             "dgw_center_inferred"})


def to_output_schema(df: pd.DataFrame, family: str) -> pd.DataFrame:
    """Reorder and dtype a metrics frame to the asset's output schema.

    Reorders to ``OUTPUT_COLUMNS[family]`` and applies the published dtypes.
    """
    cols = list(OUTPUT_COLUMNS[family])
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{family}: missing published columns {missing}")
    out = df[cols].copy()
    out["column"] = out["column"].astype(int)
    out["volume"] = out["volume"].astype(str)
    out["plane"] = out["plane"].astype(int)
    out["roi"] = out["roi"].astype(int)
    if "pref_img" in out:                        # published uses int with a -1 sentinel
        out["pref_img"] = out["pref_img"].fillna(-1).astype(int)
    if "is_responsive" in out:                   # published writes float 0.0/1.0
        out["is_responsive"] = out["is_responsive"].astype(float)
    for c in BOOLEAN_COLUMNS:
        if c in out:                             # published writes True/False
            out[c] = out[c].astype(bool)
    return out


def _reliability_on(plane, trace_key, starts, codes, *, n_trials, n_conditions,
                    window=None, frames=None):
    """`trial_reliability` recomputed on a second trace type, or NaN if it is absent.

    Reliability is reported twice — once on the trace the family's metrics use (events)
    and once on dF/F — because the two answer different questions and, on sparse events,
    disagree substantially. Events are exactly zero most of the time, so a repeat's
    response vector is mostly flat and its correlation with another repeat rests on a
    handful of frames; dF/F carries a continuous signal and is what the white paper's
    Figure 18 reports. Shipping both makes "how reproducible are the events every other
    metric is built on?" a question the asset can answer.

    Returns NaN for every ROI when `trace_key` is not loaded, so a plane loaded with a
    single trace type still produces a valid frame rather than raising.
    """
    traces = plane.traces.get(trace_key)
    if traces is None:
        return np.full(plane.n_rois, np.nan)
    if frames is not None:
        sweeps = tr.sweep_responses_frames(traces, plane.timestamps, starts, int(frames))
    else:
        sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None)
    return tr.trial_reliability(
        tr.trial_array(sweeps, codes, n_trials=n_trials, n_conditions=n_conditions))


def absent_frame(plane, family: str, mouse: Optional[str] = None) -> pd.DataFrame:
    """Identity rows with no metrics, for a session that did not run this stimulus.

    The pre-flight found all six families present in all 25 sessions of this asset, so
    this is insurance rather than a code path in daily use. It exists because
    `stimulus_trials` returns an *empty frame* for a missing stimulus instead of raising:
    without a guard, an absent stimulus would flow into the metric functions and come out
    as confident nonsense rather than as an absence.

    Booleans are set False rather than left NaN — `to_output_schema` casts them with
    `astype(bool)`, and `bool(nan)` is **True**, which would report a receptive field, or
    an imputed aperture centre, for every ROI in a session that never saw the stimulus.
    Which columns those are comes from `BOOLEAN_COLUMNS`, shared with `to_output_schema`;
    it used to be a `has_rf_` prefix test here, and adding `dgw_center_inferred`
    immediately broke it.
    """
    out = roi_frame(plane, mouse=mouse)
    for column in OUTPUT_COLUMNS[family]:
        if column in out:
            continue
        out[column] = False if column in BOOLEAN_COLUMNS else np.nan
    return out
