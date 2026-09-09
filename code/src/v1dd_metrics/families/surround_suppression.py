"""Surround suppression: the windowed grating response against the full-field one,
and whether each cell's receptive field was inside the aperture.

See docs/families/surround_suppression.md."""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from ..common import _metric_index
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import roi_frame
from .drifting_gratings import DGResult, vonmises_pref_dir, vonmises_two_peak

SSI_COLUMNS = ["ssi", "ssi_avg", "ssi_avg_at_pref_sf", "ssi_running",
               "ssi_running_avg_at_pref_sf", "ssi_stationary",
               "ssi_stationary_avg_at_pref_sf", "ssi_tuning_fit"]


def surround_suppression_metrics(
    dgw: "DGResult",
    dgf: "DGResult",
    plane,
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    mouse: Optional[str] = None,
    containment: Optional[pd.DataFrame] = None,
    center: Optional[Sequence[float]] = None,
    center_inferred: bool = False,
) -> pd.DataFrame:
    """Eight SSI variants, all ``(W - F) / (W + F)`` at the windowed preferred condition.

    ``containment`` splices in the RF-containment frame from ``window_containment``;
    omit it and those columns are NaN. ``center`` overrides the session's recorded
    aperture position (used for imputed centres). ``center_inferred`` flags the rows.
    """
    n_rois = plane.n_rois
    out = {m: np.full(n_rois, np.nan) for m in SSI_COLUMNS}

    thr = config.running_threshold_cm_s
    W, F = dgw.trial_responses, dgf.trial_responses
    W_run = np.where((dgw.trial_running_speeds > thr)[None], W, np.nan)
    W_stat = np.where((dgw.trial_running_speeds < thr)[None], W, np.nan)
    F_run = np.where((dgf.trial_running_speeds > thr)[None], F, np.nan)
    F_stat = np.where((dgf.trial_running_speeds < thr)[None], F, np.nan)

    def m(x):
        """nan-mean that returns NaN for an all-NaN slice instead of warning."""
        finite = np.isfinite(x)
        return np.nan if not finite.any() else float(np.mean(x[finite]))

    for roi in range(n_rois):
        di, si = dgw.pref_cond_index[roi]
        if di < 0 or si < 0:
            continue

        out["ssi"][roi] = _metric_index(m(W[roi, di, si]), m(F[roi, di, si]))
        out["ssi_avg"][roi] = _metric_index(m(W[roi]), m(F[roi]))
        out["ssi_avg_at_pref_sf"][roi] = _metric_index(m(W[roi, :, si]), m(F[roi, :, si]))
        out["ssi_running_avg_at_pref_sf"][roi] = _metric_index(
            m(W_run[roi, :, si]), m(F_run[roi, :, si]))
        out["ssi_stationary_avg_at_pref_sf"][roi] = _metric_index(
            m(W_stat[roi, :, si]), m(F_stat[roi, :, si]))

        for key, wa, fa in (("ssi_stationary", W_stat, F_stat),
                            ("ssi_running", W_run, F_run)):
            ws, fs = wa[roi, di, si], fa[roi, di, si]
            ws, fs = ws[np.isfinite(ws)], fs[np.isfinite(fs)]
            if len(ws) >= config.ssi_min_trials and len(fs) >= config.ssi_min_trials:
                out[key][roi] = _metric_index(ws.mean(), fs.mean())

        wp, fp = dgw.tuning_params[roi, si], dgf.tuning_params[roi, si]
        if np.isfinite(wp).all() and np.isfinite(fp).all():
            d0 = vonmises_pref_dir(wp)
            off_w = 0.0 if config.ssi_tuning_fit_includes_baseline else float(wp[-1])
            off_f = 0.0 if config.ssi_tuning_fit_includes_baseline else float(fp[-1])
            out["ssi_tuning_fit"][roi] = _metric_index(
                float(vonmises_two_peak(d0, *wp)) - off_w,
                float(vonmises_two_peak(d0, *fp)) - off_f)

    frame = roi_frame(plane, mouse=mouse)
    for k, v in out.items():
        frame[k] = v
    az, el = dgw.center if center is None else (float(center[0]), float(center[1]))
    frame["dgw_center_azimuth"] = np.full(n_rois, az)
    frame["dgw_center_elevation"] = np.full(n_rois, el)
    frame["dgw_center_inferred"] = np.full(n_rois, bool(center_inferred), dtype=bool)
    for c in CONTAINMENT_COLUMNS:
        if containment is None:
            frame[c] = np.full(n_rois, np.nan)
        else:
            if len(containment) != n_rois:
                raise ValueError(
                    f"containment has {len(containment)} rows, plane has {n_rois} ROIs")
            frame[c] = np.asarray(containment[c], dtype=np.float64)
    return frame


CONTAINMENT_COLUMNS = ["dgw_rf_distance_on", "dgw_rf_distance_off",
                       "dgw_rf_overlap_on", "dgw_rf_overlap_off"]


def _window_coverage(azimuths, altitudes, center, radius: float, sub: int = 8):
    """Fraction of each stimulus pixel's area inside the aperture disc.

    Returns ``(n_rows, n_cols)`` coverage, or ``None`` if the centre is non-finite.
    ``sub`` controls sub-pixel sampling resolution.
    """
    az = np.asarray(azimuths, dtype=np.float64)
    alt = np.asarray(altitudes, dtype=np.float64)
    caz, cel = float(center[0]), float(center[1])
    if not (np.isfinite(caz) and np.isfinite(cel)):
        return None

    def pitch(v, name):
        d = np.diff(v)
        if len(d) and not np.allclose(d, d[0]):
            raise ValueError(f"{name} are not evenly spaced: {np.unique(np.round(d, 6))}")
        return float(abs(d[0])) if len(d) else 0.0

    p_az, p_alt = pitch(az, "azimuths"), pitch(alt, "altitudes")
    offs = (np.arange(sub) + 0.5) / sub - 0.5                    # sub-cell centres
    d_az = (az[:, None] + offs[None, :] * p_az) - caz            # (n_cols, sub)
    d_alt = (alt[:, None] + offs[None, :] * p_alt) - cel         # (n_rows, sub)
    inside = ((d_alt ** 2)[:, :, None, None] + (d_az ** 2)[None, None, :, :]
              <= radius * radius)                                # (rows, sub, cols, sub)
    return inside.mean(axis=(1, 3))                              # (n_rows, n_cols)


def window_containment(
    rf_frame: pd.DataFrame,
    rf_map: np.ndarray,
    lsn: Mapping[str, Any],
    center: Sequence[float],
    *,
    config: MetricConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """RF-to-aperture distance and overlap for each ROI.

    ``rf_frame`` carries RF centres; ``rf_map`` is ``(n_rois, 2, n_rows, n_cols)``
    pre-threshold. ``lsn`` supplies the stimulus grid; ``center`` is the aperture
    ``(azimuth, elevation)``. Returns a frame of ``CONTAINMENT_COLUMNS``, all NaN where
    the centre is non-finite or the ROI has no field.
    """
    n_rois = len(rf_frame)
    out = {c: np.full(n_rois, np.nan) for c in CONTAINMENT_COLUMNS}
    cov = _window_coverage(lsn["azimuths"], lsn["altitudes"], center,
                           config.dgw_window_radius_deg)
    if cov is None:                       # no recorded aperture -> everything stays NaN
        return pd.DataFrame(out, index=rf_frame.index)

    caz, cel = float(center[0]), float(center[1])
    rf_map = np.asarray(rf_map)
    for i, sub in enumerate(("on", "off")):
        out[f"dgw_rf_distance_{sub}"] = np.hypot(
            rf_frame[f"azimuth_rf_{sub}"].to_numpy(dtype=np.float64) - caz,
            rf_frame[f"altitude_rf_{sub}"].to_numpy(dtype=np.float64) - cel)
        w = np.asarray(rf_map[:, i, :, :], dtype=np.float64).copy()
        w[w < config.rf_frac_thresh] = 0.0
        den = w.sum(axis=(1, 2))
        num = (w * cov[None, :, :]).sum(axis=(1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            out[f"dgw_rf_overlap_{sub}"] = np.where(den > 0, num / np.where(den > 0, den, 1.0),
                                                    np.nan)
    return pd.DataFrame(out, index=rf_frame.index)
