"""Drifting gratings, full field and windowed: tuning, selectivity, and the
aperture geometry the windowed condition was presented through.

See docs/families/drifting_gratings.md."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .. import responses as tr
from ..common import _condition_codes, _lifetime_sparseness_chunked, _nanmean, _ratio
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import roi_frame

def window_center(trials: pd.DataFrame) -> Tuple[float, float]:
    """The `(azimuth, elevation)` of the grating aperture for one session, or `(nan, nan)`.

    Pass **non-blank sweeps only**. Takes the first distinct non-NaN value per axis.
    """
    out = []
    for col in ("center_azimuth", "center_elevation"):
        vals = trials[col].dropna().unique() if col in trials.columns else []
        out.append(float(vals[0]) if len(vals) else np.nan)
    return out[0], out[1]


@dataclass(frozen=True)
class WindowCenters:
    """Per-session aperture centres after imputation, plus what was done and why.

    `centers` and `inferred` are keyed by `(column, volume)` with **volume as a string** —
    volumes run 1-9 and a-f across this project, and an int key would not match.
    """

    centers: Dict[Tuple[int, str], Tuple[float, float]]
    inferred: Dict[Tuple[int, str], bool]
    provenance: Dict[str, Any]


def infer_window_centers(
    observed: Mapping[Tuple[int, str], Tuple[float, float]],
    *,
    config: MetricConfig = DEFAULT_CONFIG,
) -> WindowCenters:
    """Fill missing aperture centres from the median of the same cortical column.

    ``observed`` maps ``(column, volume)`` to ``(azimuth, elevation)``; sessions with NaN
    in either axis are candidates for imputation. A session donates only if it has both
    axes, and a column with no donors is left NaN. Controlled by ``config.impute_dgw_center``.
    """
    by_column: Dict[int, list] = {}
    for key in observed:
        by_column.setdefault(int(key[0]), []).append(key)

    centers: Dict[Tuple[int, str], Tuple[float, float]] = {}
    inferred: Dict[Tuple[int, str], bool] = {}
    col_prov: Dict[str, Any] = {}
    filled: list = []
    partial: list = []

    for col in sorted(by_column):
        keys = sorted(by_column[col], key=lambda k: str(k[1]))
        donors, missing = [], []
        for k in keys:
            az, el = observed[k]
            az, el = float(az), float(el)
            if np.isfinite(az) and np.isfinite(el):
                donors.append((k, az, el))
            else:
                missing.append(k)
                # One of the two present is neither a donor nor a clean absence; say so.
                if np.isfinite(az) != np.isfinite(el):
                    partial.append({"column": col, "volume": str(k[1]),
                                    "azimuth": az if np.isfinite(az) else None,
                                    "elevation": el if np.isfinite(el) else None})

        az_vals = np.array([d[1] for d in donors], dtype=np.float64)
        el_vals = np.array([d[2] for d in donors], dtype=np.float64)
        med = ((float(np.median(az_vals)), float(np.median(el_vals)))
               if len(donors) else (np.nan, np.nan))

        for k, az, el in donors:
            centers[k] = (az, el)
            inferred[k] = False
        for k in missing:
            do_fill = bool(config.impute_dgw_center) and len(donors) > 0
            centers[k] = med if do_fill else (np.nan, np.nan)
            inferred[k] = do_fill
            if do_fill:
                filled.append({"column": col, "volume": str(k[1]),
                               "azimuth": med[0], "elevation": med[1]})

        col_prov[str(col)] = {
            "n_donors": len(donors),
            "donor_volumes": [str(d[0][1]) for d in donors],
            "median_azimuth": med[0] if len(donors) else None,
            "median_elevation": med[1] if len(donors) else None,
            "spread_azimuth": float(np.ptp(az_vals)) if len(donors) else None,
            "spread_elevation": float(np.ptp(el_vals)) if len(donors) else None,
            "n_distinct_azimuth": int(len(np.unique(az_vals))) if len(donors) else 0,
            "n_distinct_elevation": int(len(np.unique(el_vals))) if len(donors) else 0,
            "missing_volumes": [str(k[1]) for k in missing],
        }

    prov = {
        "enabled": bool(config.impute_dgw_center),
        "n_sessions": len(observed),
        "n_measured": int(sum(1 for k, v in inferred.items()
                              if not v and np.isfinite(centers[k][0]))),
        "n_inferred": int(sum(1 for v in inferred.values() if v)),
        "n_unfilled": int(sum(1 for k, v in inferred.items()
                              if not v and not np.isfinite(centers[k][0]))),
        "filled": filled,
        "partial_sessions": partial,
        "columns": col_prov,
    }
    return WindowCenters(centers=centers, inferred=inferred, provenance=prov)


@dataclass
class DGResult:
    """Drifting-gratings metrics plus the intermediates surround suppression needs."""

    metrics: pd.DataFrame                 # per-ROI, published column names
    trial_responses: np.ndarray           # (n_rois, n_dir, n_sf, n_trials), NaN-padded
    trial_running_speeds: np.ndarray      # (n_dir, n_sf, n_trials) cm/s -- no ROI axis
    pref_cond_index: np.ndarray           # (n_rois, 2) [dir_idx, sf_idx], -1 if invalid
    tuning_params: np.ndarray             # (n_rois, n_sf, 6) von Mises, NaN if no fit
    dir_list: np.ndarray
    sf_list: np.ndarray
    blank_responses: np.ndarray           # (n_rois, n_blank)
    center: Tuple[float, float] = field(default_factory=lambda: (np.nan, np.nan))


def vonmises_two_peak(x, scale_1, k_1, x0, scale_2, k_2, b):
    """Two 180-degree-opposed von Mises bumps plus an offset. x is in degrees."""
    x = np.asarray(x, dtype=np.float64)
    return (scale_1 * np.exp(k_1 * np.cos(np.deg2rad(x - x0)))
            + scale_2 * np.exp(k_2 * np.cos(np.deg2rad(x - x0 - 180)))
            + b)


_VONMISES_BOUNDS = (
    (0, 0, 0, 0, 0, 0),
    (np.inf, np.inf, 360, np.inf, np.inf, np.inf),
)


def vonmises_data_p0(x: np.ndarray, y: np.ndarray) -> tuple:
    """Initial guess taken from the curve: peak location, both peak heights, baseline.

    The two bumps are 180 degrees apart, so the second is read at the direction opposite
    the peak. Heights are divided by e because the model reaches ``scale * exp(k)`` at its
    own peak with ``k = 1``.
    """
    b = max(float(np.nanmin(y)), 0.0)
    i = int(np.nanargmax(y))
    x0 = float(x[i])
    opp = (x0 + 180.0) % 360.0
    j = int(np.argmin(np.abs(((np.asarray(x, float) - opp + 180.0) % 360.0) - 180.0)))
    amp1 = max(float(y[i]) - b, 0.0) / np.e
    amp2 = max(float(y[j]) - b, 0.0) / np.e
    return (amp1, 1.0, x0, amp2, 1.0, b)


def vonmises_two_peak_fit(x, y, p0=(0.1, 1, 180, 0.01, 1, 0.001),
                          max_fn_calls=(2000, 10000), data_p0: bool = False):
    """Least-squares fit of the two-peak von Mises, or None if it never converges.

    ``data_p0`` derives the starting guess from the curve instead of using ``p0``, which
    is faster but can settle in a different local minimum.
    """
    from scipy.optimize import curve_fit

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    good = np.isfinite(y)
    if good.sum() < 6:
        return None
    if data_p0:
        p0 = vonmises_data_p0(x[good], y[good])
    for maxfev in max_fn_calls:
        try:
            params, _ = curve_fit(vonmises_two_peak, x[good], y[good], maxfev=maxfev,
                                  bounds=_VONMISES_BOUNDS, p0=p0)
            return params
        except (RuntimeError, ValueError):
            continue
    return None


def vonmises_pref_dir(params) -> float:
    """Preferred direction of a fitted curve: whichever of the two peaks is taller.

    Peak height is baseline-subtracted, ``f(x) - b``. ``ssi_tuning_fit`` evaluates the
    curve the same way unless ``ssi_tuning_fit_includes_baseline`` is set.
    """
    x0 = float(params[2])
    x1 = (x0 + 180.0) % 360.0
    a0 = vonmises_two_peak(x0, *params) - params[-1]
    a1 = vonmises_two_peak(x1, *params) - params[-1]
    return x0 if a0 > a1 else x1


def dg_metrics_from_trials(
    ta: np.ndarray,
    dir_list: np.ndarray,
    sf_list: np.ndarray,
    is_valid: np.ndarray,
    *,
    family: str = "drifting_gratings_full",
    fit_sf_index: Optional[np.ndarray] = None,
    config: MetricConfig = DEFAULT_CONFIG,
) -> dict:
    """Grating metrics from the trial array; responsiveness excluded (needs the trace).

    ``ta`` is ``(n_rois, n_dir, n_sf, n_trials)``, NaN-padded. ``is_valid`` marks ROIs
    that passed segmentation. ``fit_sf_index`` selects one SF to fit per ROI, or ``None``
    for all.
    """
    import warnings

    n_rois, n_dir, n_sf, n_trials = ta.shape
    roi_ix = np.arange(n_rois)
    mean_tr = _nanmean(ta, axis=3)                       # (n_rois, n_dir, n_sf)

    k_fill = np.nan_to_num(mean_tr, nan=-1.0).reshape(n_rois, -1).argmax(axis=1)
    k_skip = np.where(np.isfinite(mean_tr), mean_tr, -np.inf).reshape(n_rois, -1).argmax(axis=1)
    pref_dir_idx, pref_sf_idx = np.divmod(k_skip, n_sf)

    n_diverge = int(np.sum(k_fill != k_skip))
    if n_diverge:
        warnings.warn(
            f"{family}: the two preferred-condition definitions disagree on {n_diverge} "
            f"of {n_rois} ROIs; SSI uses the fillna(-1) one, osi/dsi the other"
        )

    k_pref = k_fill if config.pref_cond_fillna else k_skip
    pref_dir_pref, pref_sf_pref = np.divmod(k_pref, n_sf)
    pref_cond_index = np.stack([pref_dir_pref, pref_sf_pref], axis=1).astype(int)
    pref_cond_index[~is_valid] = -1
    if not config.pref_cond_fillna:
        no_response = ~np.isfinite(mean_tr).any(axis=(1, 2))
        pref_cond_index[no_response] = -1

    tuning = mean_tr[roi_ix, :, pref_sf_idx]             # (n_rois, n_dir) at preferred SF
    pref = tuning[roi_ix, pref_dir_idx]
    null_r = tuning[roi_ix, (pref_dir_idx + 6) % 12]
    orth_r = 0.5 * (tuning[roi_ix, (pref_dir_idx + 3) % 12]
                    + tuning[roi_ix, (pref_dir_idx - 3) % 12])

    zn = config.zero_denominator_nan
    osi = _ratio(pref - orth_r, pref + orth_r, zero_to_nan=zn)
    dsi = _ratio(pref - null_r, pref + null_r, zero_to_nan=zn)

    theta = np.deg2rad(dir_list.astype(float))
    L_norm = tuning.sum(axis=1)
    L_ori = tuning @ np.exp(2j * theta)
    with np.errstate(invalid="ignore", divide="ignore"):
        gosi = np.abs(np.where(L_norm != 0, L_ori / np.where(L_norm != 0, L_norm, 1.0),
                               L_ori))
    vec = np.nan_to_num(tuning, nan=0.0) @ np.exp(1j * theta)
    pref_dir_mean = np.degrees(np.angle(vec)) % 360.0

    lifetime = _lifetime_sparseness_chunked(
        ta.transpose(1, 2, 3, 0).reshape(-1, n_trials, n_rois),
        over=config.lifetime_sparseness_over)

    tuning_params = np.full((n_rois, n_sf, 6), np.nan)
    if config.fit_tuning_curves:
        _fit_sf = fit_sf_index
        if _fit_sf is None and family.endswith("windowed"):
            _fit_sf = pref_cond_index[:, 1]
        if config.fit_all_sf:
            _fit_sf = None
        for roi in range(n_rois):
            if not is_valid[roi]:
                continue
            if _fit_sf is not None:
                sf_i = int(_fit_sf[roi])
                if sf_i < 0:
                    continue
                p = vonmises_two_peak_fit(dir_list, mean_tr[roi, :, sf_i],
                                          data_p0=config.vonmises_data_p0)
                if p is not None:
                    tuning_params[roi, sf_i] = p
            else:
                for sf_i in range(n_sf):
                    p = vonmises_two_peak_fit(dir_list, mean_tr[roi, :, sf_i],
                                          data_p0=config.vonmises_data_p0)
                    if p is not None:
                        tuning_params[roi, sf_i] = p

    return {
        "dsi": dsi,
        "gosi": gosi,
        "lifetime_sparseness": lifetime,
        "osi": osi,
        "pref_dir_mean": pref_dir_mean,
        "preferred_dir": np.where(pref_cond_index[:, 0] >= 0,
                                  dir_list[pref_cond_index[:, 0]], np.nan),
        "preferred_sf": np.where(pref_cond_index[:, 1] >= 0,
                                 sf_list[pref_cond_index[:, 1]], np.nan),
        "pref_cond_index": pref_cond_index,
        "tuning_params": tuning_params,
        "pref_dir_idx": pref_dir_idx,
        "pref_sf_idx": pref_sf_idx,
    }


def drifting_gratings_metrics(
    plane,
    trials: pd.DataFrame,
    is_blank: np.ndarray,
    spont: Sequence[float],
    running: Optional[Sequence[np.ndarray]] = None,
    *,
    dg_type: str = "full",
    fit_sf_index: Optional[np.ndarray] = None,
    config: MetricConfig = DEFAULT_CONFIG,
    rng: Optional[np.random.Generator] = None,
    mouse: Optional[str] = None,
) -> "DGResult":
    """Drifting-gratings metrics for one plane. Returns a ``DGResult``.

    ``dg_type`` is ``"full"`` or ``"windowed"``; the computation is identical.
    ``fit_sf_index`` selects one SF to fit per ROI (windowed uses its own preferred SF
    by default). ``running`` is ``(speed_array, running_timestamps)`` or ``None``.
    """
    if not len(trials):
        raise ValueError(
            f"no drifting_gratings_{dg_type} sweeps for column {plane.column} "
            f"volume {plane.volume} plane {plane.plane}. Surround suppression consumes "
            "this result, so there is no empty value that stays honest downstream -- "
            "skip the session instead."
        )
    rng = np.random.default_rng() if rng is None else rng
    family = f"drifting_gratings_{dg_type}"
    traces = plane.traces[config.trace_type[family]]

    if config.dg_response_seconds == "per_trial":
        raise NotImplementedError("per-trial windows need a per-sweep window API")
    window = (0.0, float(config.dg_response_seconds))

    starts = trials["start_time"].to_numpy(dtype=np.float64)
    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None)

    grat = trials.loc[~is_blank]

    # per-session grating-aperture centre (NaN for full-field, which uses (0, 0) placeholders)
    center: Tuple[float, float] = window_center(grat)

    dir_list = np.sort(grat["direction"].dropna().unique())
    sf_list = np.sort(grat["spatial_frequency"].dropna().unique())
    if len(dir_list) != 12:
        raise ValueError(
            f"{family}: found {len(dir_list)} directions, expected 12. The orthogonal and "
            "null directions are hard-coded as (i +/- 3) % 12 and (i + 6) % 12, which "
            "compute silently wrong values for any other count."
        )
    n_dir, n_sf = len(dir_list), len(sf_list)

    d = _condition_codes(grat["direction"].to_numpy(), dir_list)
    s = _condition_codes(grat["spatial_frequency"].to_numpy(), sf_list)
    code = d * n_sf + s
    n_trials = int(np.bincount(code, minlength=n_dir * n_sf).max())

    ta = tr.trial_array(sweeps[~is_blank], code, n_trials=n_trials,
                        n_conditions=n_dir * n_sf)
    ta = ta.reshape(n_dir, n_sf, n_trials, plane.n_rois).transpose(3, 0, 1, 2)

    blank = sweeps[is_blank].T if bool(is_blank.any()) else np.empty((plane.n_rois, 0))

    core = dg_metrics_from_trials(ta, dir_list, sf_list, plane.is_valid,
                                  family=family, fit_sf_index=fit_sf_index, config=config)
    n_rois = plane.n_rois
    roi_ix = np.arange(n_rois)

    null_single = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, None,
        n_boot=config.dg_n_boot, n_means=1, rng=rng,
        memory_budget_mb=config.memory_budget_mb,
    )
    pref_trials = ta[roi_ix, core["pref_dir_idx"], core["pref_sf_idx"], :]
    frac = tr.frac_trials_above_null(pref_trials, null_single, p_thresh=config.sig_p_thresh)
    is_responsive = (plane.is_valid & (frac >= config.dg_frac_thresh)).astype(float)

    # per-trial running speed, shape (n_dir, n_sf, n_trials); no ROI axis
    trs = np.full((n_dir, n_sf, n_trials), np.nan)
    if running is not None:
        speed, rts = running
        pad = config.running_pad_seconds
        gstarts = grat["start_time"].to_numpy(dtype=np.float64)
        stops = grat["stop_time"].to_numpy(dtype=np.float64)
        cs, counts = tr.prefix_sums(np.asarray(speed, dtype=np.float64)[:, None])
        a = np.searchsorted(rts, gstarts - pad, side="left")
        b = np.searchsorted(rts, stops + pad, side="right")
        per_sweep = tr.window_means(cs, counts, a, b)              # (n_sweeps, 1)
        trs = tr.trial_array(per_sweep, code, n_trials=n_trials,
                             n_conditions=n_dir * n_sf).reshape(n_dir, n_sf, n_trials)

    out = roi_frame(plane, mouse=mouse)
    out["dsi"] = core["dsi"]
    out["frac_responsive_trials"] = frac
    out["gosi"] = core["gosi"]
    out["is_responsive"] = is_responsive
    out["lifetime_sparseness"] = core["lifetime_sparseness"]
    out["osi"] = core["osi"]
    out["preferred_dir"] = core["preferred_dir"]
    out["preferred_sf"] = core["preferred_sf"]
    out["pref_dir_mean"] = core["pref_dir_mean"]

    return DGResult(metrics=out, trial_responses=ta, trial_running_speeds=trs,
                    pref_cond_index=core["pref_cond_index"],
                    tuning_params=core["tuning_params"],
                    dir_list=dir_list, sf_list=sf_list, blank_responses=blank,
                    center=center)
