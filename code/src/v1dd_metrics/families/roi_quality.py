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

    Ported from `Functional Data Cell-Cell Correlations.ipynb`
    (`estimate_snr_white_noise_model`), which ships `snr_by_cell.feather` in the CCM
    asset. Same bands and same scaling, so the two assets are directly comparable — but
    **not bit-identical**: that notebook interpolates every plane onto one reference
    timebase before the FFT, while this runs per plane on its own timestamps. Expect
    agreement in distribution, not in the last digit.

    Calcium transients live around 0.1-1.5 Hz. Above them the spectrum of a well-isolated
    ROI is close to flat, so a narrow band at 2.0-2.1 Hz estimates the white-noise floor
    per bin; scaling that by the number of signal bins gives the noise power the signal
    band would contain if it held nothing but noise.

    `traces` is `(n_frames, n_rois)` as NWB stores it -- transposed internally, unlike the
    notebook's function, which expects the transpose already done.

    **dF/F only.** On deconvolved events this measures nothing useful: deconvolution has
    already removed the noise floor the reference band is meant to sample, so the
    denominator is whatever numerical residue is left rather than a noise estimate.

    Returns `(snr, signal_power, noise_power)`, each `(n_rois,)`. SNR is a power ratio,
    not decibels -- the CCM notebook plots it on a log axis from 1 to 1e5.
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
    """How much locomotion changes each neuron's activity, across three conditions.

    `(R_run - R_stat) / (R_run + R_stat)`, the conventional index and the same form as
    every `ssi_*` column. **Not** the white paper's `C*(Rmax - Rmin)/Rmax`, which divides
    by the max rather than the sum; that would put a second convention in the same asset
    for no gain, and the two are monotonically related anyway.

    Three conditions, all on **deconvolved events only**:

    * **`run_mod_dgf` / `run_mod_dgw`** — full-field and windowed gratings. Having both
      is a cross-check on `ssi_running` / `ssi_stationary`, which split the same trials at
      the same threshold but only at each ROI's preferred condition.
    * **`run_mod_spont`** — the spontaneous block, where nothing is on the screen. This is
      the control: locomotion modulates cortex whether or not there is a stimulus, so a
      grating index should be read against this rather than against zero.

    `spont_rate` is the mean activity over the whole spontaneous block, with
    `spont_rate_run` / `spont_rate_stat` the same quantity split by state. They are here
    for two reasons. They are the magnitudes behind `run_mod_spont`, so a ratio built from
    two near-zero numbers can be gated rather than trusted. And `spont_rate` is the only
    per-ROI **baseline activity level** in the asset — how much a neuron does with nothing
    on the screen — which is useful well beyond locomotion: as a normaliser for evoked
    responses, and for spotting unusually silent or hyperactive cells.

    Note `spont_rate` is **not** recoverable from the other two. `spont_run_frac` is a
    fraction of *time* on the running trace's own ~59 Hz samples (the white paper's
    definition), while the split above classifies *imaging frames*; the two are close but
    are not the weights that would recombine the state means. So the overall level is
    computed and shipped rather than left to be derived.

    **Pooled over all trials, not computed at the preferred condition.** That is forced by
    the data, not preference. `ssi_running` needs >=3 running trials at one condition out
    of 8, and only **6.3 % of ROIs in 9 of 25 sessions** clear that bar — running here is
    close to all-or-nothing per session, with many sessions 100 % stationary and two 100 %
    running. Pooling over all 192 grating trials makes the same threshold easy wherever
    the animal ran at all.

    `run_frac` (whole session) and `spont_run_frac` (the spontaneous block alone) are
    reported beside them because they say whether any of it is interpretable — the paper
    gates its locomotion analyses at a running fraction of 0.2, and the two fractions can
    differ substantially.

    **Events only, deliberately — a ratio index is not safe on a signed trace.** With
    non-negative events `R_run + R_stat` is a sum of magnitudes and vanishes only for a
    silent cell. On signed dF/F the same expression breaks in two ways: near-cancelling
    responses of opposite sign give an unbounded index (+0.050 vs -0.049 -> 99), and when
    both responses are negative the sign inverts, so a suppressed cell that is *less*
    suppressed while running scores negative.

    For `spont_rate` the case against dF/F is stronger still, and different: dF/F is
    defined against a rolling baseline, so its mean over a long stimulus-free block is
    ~0 **by construction**. It would not be unstable, it would be uninformative.

    The white paper's `C*(Rmax - Rmin)/Rmax` does not rescue the ratio; it is worse. Its
    denominator can itself be negative -- `max(-0.01, -0.05) = -0.01` -- giving -4.0 for
    that same both-negative case. Rectifying, or dividing by `|R_run| + |R_stat|`, would
    fix the sign and the bound, but neither fixes the deeper problem below, and a dF/F
    column would then need a different formula from every other index in this asset.

    **What no denominator fixes**, and which applies to events too: when both responses
    are near zero the ratio is large and meaningless. `R_run = 1e-6, R_stat = 2e-6` gives
    -0.33 under every variant. That is a signal-to-noise problem. Gate on magnitude before
    trusting a value from a quiet cell -- for the gratings the raw per-trial responses and
    running speeds are both in `tuning_curves_*.npz`, so `R_run` and `R_stat` can be
    recomputed and thresholded however you like; for the spontaneous block, gate on
    `spont_rate_run` / `spont_rate_stat`.

    (Note the contrast with `reliability`, which *is* reported on both trace types. That
    is a correlation -- invariant to sign and scale, no denominator -- so the argument
    against a dF/F variant here does not apply there.)
    """
    n_rois = plane.n_rois
    out = {c: np.full(n_rois, np.nan) for c in ROI_SUMMARY_COLUMNS}
    thr = config.running_threshold_cm_s
    n_min = config.ssi_min_trials

    # Recording quality, from dF/F -- the one quantity in this table that is not events.
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

        # Spontaneous has no trials, so each imaging frame is the unit: classify the frame
        # by the running speed over its own interval, then average the trace within each
        # class. No padding, unlike the trial windows -- frames are contiguous, so a pad
        # would count the same running samples into both neighbours.
        t0, t1 = float(spont[0]), float(spont[1])
        in_spont = (rts >= t0) & (rts <= t1)
        out["spont_run_frac"] = np.full(
            n_rois, float((speed[in_spont] > thr).mean()) if in_spont.any() else np.nan)

        ts = np.asarray(plane.timestamps, dtype=np.float64)

        # Running against the continuous dF/F trace, whole session. On dF/F rather
        # than events, and with no state split, so it is finite even in sessions where
        # run_mod_* is all-NaN for want of trials in one state.
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
            # Same trace the gratings use, so every number in this table is in one
            # currency: mean event magnitude per imaging sample.
            traces = plane.traces.get(config.trace_type["drifting_gratings_full"])
            if traces is not None:
                block = np.asarray(traces)[frames]           # (n_spont_frames, n_rois)
                # Gated separately, on purpose. A session where the animal never ran still
                # has a baseline rate and a stationary rate; requiring both states would
                # throw those away, and 13 of 25 sessions here are one-sided.
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
