"""Per-ROI signal quality, spontaneous activity, and locomotion modulation.

See docs/families/roi_quality.md."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .. import responses as tr
from ..common import _metric_index, _nanmean
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import roi_frame
from .drifting_gratings import DGResult

ROI_SUMMARY_COLUMNS = ["snr", "signal_power", "noise_power",
                       "run_frac", "spont_run_frac",
                       "spont_rate", "spont_rate_run", "spont_rate_stat",
                       "run_mod_dgf", "run_mod_dgw", "run_mod_spont",
                       "run_corr_dff"]


def _pearson_columns(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson r between each column of ``x`` (n_frames, n_rois) and the vector ``y``.

    Frames where either side is non-finite are dropped per ROI, so a trace with gaps
    is still correlated over the samples it has. Returns NaN where either side has no
    variance.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)[:, None]
    n = finite.sum(axis=0)
    xs = np.where(finite, x, 0.0)
    ys = np.where(finite, y[:, None], 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = xs.sum(axis=0) / n, ys.sum(axis=0) / n
        cov = (xs * ys).sum(axis=0) / n - mx * my
        vx = (xs * xs).sum(axis=0) / n - mx * mx
        vy = (ys * ys).sum(axis=0) / n - my * my
        r = cov / np.sqrt(vx * vy)
    return np.where((n >= 2) & (vx > 0) & (vy > 0), r, np.nan)


def spectral_snr(traces, fs: float, signal_band=(0.1, 1.5), noise_band=(2.0, 2.1),
                 demean: bool = True):
    """Per-ROI SNR from signal-band power against a white-noise reference band.

    ``traces`` is ``(n_frames, n_rois)``. Returns ``(snr, signal_power, noise_power)``,
    each ``(n_rois,)``. SNR is a power ratio, not decibels. dF/F only.
    """
    x = np.asarray(traces, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError(f"expected (n_frames, n_rois), got {x.shape}")
    x = x.T                                            # -> (n_rois, n_frames)
    n_time = x.shape[1]
    if demean:
        x = x - x.mean(axis=1, keepdims=True)          # kill the DC bin

    freqs = np.fft.rfftfreq(n_time, d=1.0 / float(fs))
    if noise_band[1] > freqs[-1]:
        raise ValueError(
            f"noise_band {noise_band} exceeds Nyquist {freqs[-1]:.3f} Hz at "
            f"fs={fs:.3f} Hz -- choose a band below it")
    power = np.abs(np.fft.rfft(x, axis=1)) ** 2

    sig_mask = (freqs >= signal_band[0]) & (freqs <= signal_band[1])
    noi_mask = (freqs >= noise_band[0]) & (freqs <= noise_band[1])
    if not sig_mask.any() or not noi_mask.any():
        raise ValueError(f"empty frequency band: signal={int(sig_mask.sum())} "
                         f"noise={int(noi_mask.sum())} bins at fs={fs:.3f} Hz")

    signal_power = power[:, sig_mask].sum(axis=1)
    noise_power = power[:, noi_mask].mean(axis=1) * int(sig_mask.sum())
    return signal_power / (noise_power + 1e-12), signal_power, noise_power


def _run_modulation(resp, speeds, thr, min_trials):
    """(R_run - R_stat) / (R_run + R_stat) per ROI, pooled over every trial.

    `resp` is (n_rois, n_dir, n_sf, n_trials), `speeds` (n_dir, n_sf, n_trials).
    """
    if resp is None or speeds is None:
        return None
    run = speeds > thr
    stat = speeds < thr                      # strict both sides: exactly thr is neither,
    if run.sum() < min_trials or stat.sum() < min_trials:   # matching the ssi convention
        return None
    r = _nanmean(np.where(run[None], resp, np.nan).reshape(resp.shape[0], -1), axis=1)
    s = _nanmean(np.where(stat[None], resp, np.nan).reshape(resp.shape[0], -1), axis=1)
    return _metric_index(r, s)


def roi_summary_metrics(
    plane,
    dgw: Optional["DGResult"],
    dgf: Optional["DGResult"],
    spont: Sequence[float],
    running: Optional[Sequence[np.ndarray]] = None,
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    mouse: Optional[str] = None,
) -> pd.DataFrame:
    """SNR, spontaneous rates, and locomotion modulation for one plane.

    ``dgw`` and ``dgf`` supply per-trial responses for ``run_mod_dgf``/``run_mod_dgw``.
    ``running`` is ``(speed_array, running_timestamps)`` or ``None``.
    ``spont`` is ``(start_time, end_time)`` of the spontaneous block.
    """
    n_rois = plane.n_rois
    out = {c: np.full(n_rois, np.nan) for c in ROI_SUMMARY_COLUMNS}
    thr = config.running_threshold_cm_s
    n_min = config.ssi_min_trials

    dff = plane.traces.get("dff")
    if dff is not None and plane.dt and np.isfinite(plane.dt):
        snr, sig, noi = spectral_snr(dff, fs=1.0 / float(plane.dt),
                                     signal_band=tuple(config.snr_signal_band),
                                     noise_band=tuple(config.snr_noise_band))
        out["snr"], out["signal_power"], out["noise_power"] = snr, sig, noi

    for name, res in (("dgf", dgf), ("dgw", dgw)):
        if res is None:
            continue
        got = _run_modulation(res.trial_responses, res.trial_running_speeds, thr, n_min)
        if got is not None:
            out[f"run_mod_{name}"] = got

    if running is not None:
        speed, rts = np.asarray(running[0], dtype=np.float64), np.asarray(running[1])
        moving = speed > thr
        out["run_frac"] = np.full(n_rois, float(moving.mean()) if moving.size else np.nan)

        t0, t1 = float(spont[0]), float(spont[1])
        in_spont = (rts >= t0) & (rts <= t1)
        out["spont_run_frac"] = np.full(
            n_rois, float((speed[in_spont] > thr).mean()) if in_spont.any() else np.nan)

        ts = np.asarray(plane.timestamps, dtype=np.float64)

        if dff is not None and plane.dt and np.isfinite(plane.dt):
            cs_all, counts_all = tr.prefix_sums(speed[:, None])
            a_all = np.searchsorted(rts, ts, side="left")
            b_all = np.searchsorted(rts, ts + plane.dt, side="right")
            speed_per_frame = tr.window_means(cs_all, counts_all, a_all, b_all)[:, 0]
            out["run_corr_dff"] = _pearson_columns(np.asarray(dff), speed_per_frame)

        frames = np.flatnonzero((ts >= t0) & (ts <= t1))
        if frames.size:
            cs, counts = tr.prefix_sums(speed[:, None])
            a = np.searchsorted(rts, ts[frames], side="left")
            b = np.searchsorted(rts, ts[frames] + plane.dt, side="right")
            per_frame = tr.window_means(cs, counts, a, b)[:, 0]
            run_f = per_frame > thr
            stat_f = per_frame < thr
            traces = plane.traces.get(config.trace_type["drifting_gratings_full"])
            if traces is not None:
                block = np.asarray(traces)[frames]           # (n_spont_frames, n_rois)
                out["spont_rate"] = _nanmean(block, axis=0)
                if run_f.sum() >= n_min:
                    out["spont_rate_run"] = _nanmean(block[run_f], axis=0)
                if stat_f.sum() >= n_min:
                    out["spont_rate_stat"] = _nanmean(block[stat_f], axis=0)
                if run_f.sum() >= n_min and stat_f.sum() >= n_min:
                    out["run_mod_spont"] = _metric_index(out["spont_rate_run"],
                                                         out["spont_rate_stat"])

    frame = roi_frame(plane, mouse=mouse)
    for k, v in out.items():
        frame[k] = v
    return frame
