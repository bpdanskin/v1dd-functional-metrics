"""Prototype: greedy pixelwise RF (v1dd_physiology / Dan Millman) vs our fraction-threshold RF.

Both operate on the true 8x14 LSN grid and L0 events, so this isolates the
*significance method*, which is the only real difference once resolution is
controlled (the NWB's 16x28 template is a verified 2x2 upsample of an 8x14
stimulus, so it carries no finer spatial information).

  ours   : per pixel, fraction of that pixel's presentations whose response
           exceeds the ROI's bootstrapped spontaneous 95th percentile;
           detect where fraction >= rf_frac_thresh (0.25). No multiple-
           comparisons correction.
  greedy : per pixel, STA of mean events; per-pixel bootstrap p-value
           (resample sweeps with replacement); Holm-Sidak correction across
           all 2*n_pixels tests; detect where corrected p < alpha.

Reported per method (ON subfield): detection rate, single-component fraction
(contiguity), RF area (deg^2), and — in synthetic mode where truth is known —
true-positive detection, center error, and spurious-pixel / false-cell rates.

Usage
-----
  # local, no data needed — synthetic ground truth:
  python greedy_rf_comparison.py --synthetic

  # in the capsule, on real sessions:
  python greedy_rf_comparison.py --sessions 2
"""
import argparse
import pathlib
import sys
import time

import numpy as np
from scipy import ndimage

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "src"))

from v1dd_metrics import responses as tr


# ----------------------------------------------------------------- greedy method
def holm_sidak_reject(pvals, alpha):
    """Step-down Holm-Sidak rejection over axis 0 (tests) for each column (ROI).

    Matches statsmodels.multipletests(..., method='hs')[0]. Returns a boolean
    reject mask the same shape as ``pvals``.
    """
    m = pvals.shape[0]
    order = np.argsort(pvals, axis=0)
    sorted_p = np.take_along_axis(pvals, order, axis=0)
    thr = 1.0 - (1.0 - alpha) ** (1.0 / (m - np.arange(m)))     # per-rank threshold
    passed = sorted_p <= thr[:, None]
    rej_sorted = np.logical_and.accumulate(passed, axis=0)      # step-down prefix
    mask = np.zeros_like(pvals, dtype=bool)
    np.put_along_axis(mask, order, rej_sorted, axis=0)
    return mask


def greedy_sta_pvals(plane, trials, lsn, *, trace_key="events", window=(0.0, 0.5),
                     n_boot=2000, seed=1):
    """STA and per-pixel bootstrap p-values for every ROI (alpha-independent).

    Returns (STA, pvals, n_rows, n_cols), STA and pvals both (2*n_pix, n_rois)
    with ON rows then OFF rows. The bootstrap is the expensive step, so this is
    computed once and thresholded at several alphas by ``greedy_masks``.
    """
    images = np.asarray(lsn["images"])
    n_rows, n_cols = images.shape[1], images.shape[2]
    n_pix = n_rows * n_cols
    traces = plane.traces[trace_key]
    starts = trials["start_time"].to_numpy(dtype=np.float64)
    frames = trials["frame"].to_numpy().astype(int)
    pixel_on, pixel_off = lsn["pixel_on"], lsn["pixel_off"]

    stim = images[frames].reshape(len(frames), n_pix)
    A = np.concatenate([stim == pixel_on, stim == pixel_off], axis=1).T.astype(np.float32)
    sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None).astype(np.float32)
    n_sweeps = sweeps.shape[0]
    STA = A @ sweeps                               # (2*n_pix, n_rois)

    rng = np.random.default_rng(seed)
    ge = np.zeros_like(STA)                        # count of shuffles >= actual
    for _ in range(n_boot):
        idx = rng.integers(0, n_sweeps, n_sweeps)
        ge += (A @ sweeps[idx]) >= STA
    pvals = ge / n_boot
    return STA, pvals, n_rows, n_cols


def greedy_masks(STA, pvals, alpha, n_rows, n_cols):
    """Threshold precomputed greedy p-values at one alpha -> (rf_continuous, mask)."""
    m2d = holm_sidak_reject(pvals, alpha)
    n_rois = STA.shape[1]
    rf = (STA * m2d).T.reshape(n_rois, 2, n_rows, n_cols)
    mask = m2d.T.reshape(n_rois, 2, n_rows, n_cols)
    return rf.astype(np.float32), mask


# ----------------------------------------------------------------- shared metrics
def mask_metrics(mask, grid_deg):
    """Per-ROI/pol: has_rf, n_components, largest-component fraction, area (deg^2)."""
    n_rois = mask.shape[0]
    has = np.zeros((n_rois, 2), bool)
    ncomp = np.zeros((n_rois, 2), int)
    largest = np.zeros((n_rois, 2), float)
    area = np.zeros((n_rois, 2), float)
    for pol in range(2):
        for r in range(n_rois):
            m = mask[r, pol]
            tot = int(m.sum())
            if tot == 0:
                continue
            has[r, pol] = True
            area[r, pol] = tot * grid_deg * grid_deg
            lab, nc = ndimage.label(m)
            ncomp[r, pol] = nc
            sizes = ndimage.sum(m, lab, range(1, nc + 1))
            largest[r, pol] = sizes.max() / tot
    return has, ncomp, largest, area


def centroid(mask_map, alt, azi):
    """Intensity-agnostic centroid of a boolean map, in (alt, azi) degrees; nan if empty."""
    if not mask_map.any():
        return np.nan, np.nan
    rr, cc = np.nonzero(mask_map)
    return float(alt[rr].mean()), float(azi[cc].mean())


# ----------------------------------------------------------------- synthetic mode
def synthetic_plane(n_rf=40, n_noise=40, n_frames=1600, seed=0):
    """Build a synthetic plane with known ON receptive fields plus RF-less noise cells.

    Returns (plane, trials, spont, lsn, truth) where truth has per-ROI center
    (row, col) and is_rf flag. Events are non-negative and sparse (bursty), so
    both methods face the real 6 Hz detection regime.
    """
    from types import SimpleNamespace
    import pandas as pd

    rng = np.random.default_rng(seed)
    ROWS, COLS, GRID = 8, 14, 9.3
    DT = 0.16504
    PERIOD = 1.0
    N = n_rf + n_noise

    # --- LSN template: 6 ON + 6 OFF pixels per frame, sparse ---
    images = np.zeros((n_frames, ROWS, COLS), dtype=np.int8)
    for f in range(n_frames):
        on = rng.choice(ROWS * COLS, 6, replace=False)
        off = rng.choice(np.setdiff1d(np.arange(ROWS * COLS), on), 6, replace=False)
        images[f].flat[on] = 1
        images[f].flat[off] = -1

    azimuths = (np.arange(COLS) - COLS // 2 + 0.5) * GRID
    altitudes = (np.arange(ROWS) - ROWS // 2 + 0.5) * GRID
    lsn = {"images": images, "azimuths": azimuths, "altitudes": altitudes,
           "pixel_on": 1, "pixel_off": -1, "pixel_gray": 0, "pixel_values": [-1, 0, 1]}

    starts = 30.0 + np.arange(n_frames) * PERIOD
    trials = pd.DataFrame({"stim_name": "locally_sparse_noise", "start_time": starts,
                           "stop_time": starts + 0.25, "frame": np.arange(n_frames).astype(float)})
    spont_start = starts[-1] + 5.0
    spont_stop = spont_start + 300.0
    ts = np.arange(0.0, spont_stop + 20.0, DT)
    T = len(ts)

    # --- assign true RF centers to the RF cells; vary SNR ---
    centers = np.zeros((N, 2), int)
    is_rf = np.zeros(N, bool)
    snr = np.zeros(N)
    for i in range(n_rf):
        is_rf[i] = True
        centers[i] = (rng.integers(1, ROWS - 1), rng.integers(1, COLS - 1))
        snr[i] = rng.uniform(0.6, 1.6)             # event amplitude scale

    # --- spontaneous baseline events (sparse, non-negative) for every cell ---
    events = (rng.random((T, N)) < 0.02) * rng.exponential(0.01, size=(T, N))
    events = events.astype(np.float64)

    # --- inject RF-driven events: a cell fires when an ON pixel lands near its center ---
    SIGMA = 1.0
    for f, s in enumerate(starts):
        w = (ts >= s) & (ts <= s + 2 * DT)          # ~2-frame response window
        widx = np.nonzero(w)[0]
        if widx.size == 0:
            continue
        on_rc = np.argwhere(images[f] == 1)         # ON pixels this frame
        for i in range(n_rf):
            cr, cc = centers[i]
            d2 = ((on_rc[:, 0] - cr) ** 2 + (on_rc[:, 1] - cc) ** 2)
            drive = float(np.exp(-d2 / (2 * SIGMA ** 2)).max())   # nearest ON pixel only
            if drive <= 0.05:
                continue
            # bursty, sparse detection matching the 6 Hz regime (~10-25% per driving frame)
            if rng.random() < min(0.6, 0.35 * snr[i] * drive):
                amp = snr[i] * drive * rng.uniform(0.8, 1.2)
                events[widx[0], i] += amp           # event at window onset frame

    roi_table = pd.DataFrame({"column": 1, "volume": 3, "plane": 0,
                              "roi": np.arange(N), "pika_roi_confidence": 0.9})
    import v1dd_metrics.nwb as vn
    plane = vn.PlaneData(mouse_id="synthetic", depth_um=150.0, column=1, volume="3",
                         plane=0, roi=np.arange(N), is_valid=np.ones(N, bool),
                         timestamps=ts, traces={"events": events, "dff": events},
                         roi_table=roi_table, dt=DT)
    truth = SimpleNamespace(centers=centers, is_rf=is_rf, ROWS=ROWS, COLS=COLS, GRID=GRID)
    return plane, trials, (spont_start, spont_stop), lsn, truth


def _synth_report(name, mask, truth, alt, azi, n_rf, n_noise):
    """Print the ON-subfield metrics for one detection mask against ground truth."""
    has, ncomp, largest, area = mask_metrics(mask, truth.GRID)
    on_has = has[:, 0]
    tp = int((on_has & truth.is_rf).sum())
    fp_cells = int((on_has & ~truth.is_rf).sum())
    cerr, spurious, single = [], [], []
    for r in np.nonzero(on_has & truth.is_rf)[0]:
        cr, cc = truth.centers[r]
        ma, mz = centroid(mask[r, 0], alt, azi)
        cerr.append(np.hypot(ma - alt[cr], mz - azi[cc]))
        rr, ccc = np.nonzero(mask[r, 0])
        spurious.append(int((np.hypot(rr - cr, ccc - cc) > 1.5).sum()))
        single.append(ncomp[r, 0] == 1)
    cerr, spurious, single = map(np.array, (cerr, spurious, single))
    med_area = np.median(area[on_has & truth.is_rf, 0]) if tp else float("nan")
    print(f"  {name:<22} det {tp:>2}/{n_rf} ({tp/n_rf:>4.0%}) | "
          f"FP {fp_cells:>2}/{n_noise} ({fp_cells/n_noise:>4.0%}) | "
          f"single {single.mean() if tp else float('nan'):>4.0%} | "
          f"cen-err {np.median(cerr) if tp else float('nan'):>4.2f} | "
          f"area {med_area:>4.0f} | spur {spurious.mean() if tp else float('nan'):>4.2f}")


def run_synthetic(args):
    from v1dd_metrics.config import MetricConfig
    from v1dd_metrics.families.receptive_fields import receptive_field_metrics

    print("Synthetic ground-truth comparison (ON subfield)\n" + "=" * 78)
    plane, trials, spont, lsn, truth = synthetic_plane(seed=args.seed)
    n_rf = int(truth.is_rf.sum())
    n_noise = int((~truth.is_rf).sum())
    print(f"{plane.n_rois} cells: {n_rf} with a planted RF, {n_noise} RF-less (noise)")
    print("cols: detection of planted RFs | false-positive noise cells | single-component "
          "frac | median center error (deg) | median area (deg^2) | mean spurious px/cell\n")
    alt, azi = np.asarray(lsn["altitudes"]), np.asarray(lsn["azimuths"])

    cfg = MetricConfig(other_n_boot=2000)
    _, ours_map = receptive_field_metrics(plane, trials, spont, lsn, config=cfg,
                                          rng=np.random.default_rng(0))
    _synth_report("ours (frac>=0.25)", ours_map >= cfg.rf_frac_thresh,
                  truth, alt, azi, n_rf, n_noise)

    t0 = time.time()
    STA, pvals, nr, nc = greedy_sta_pvals(plane, trials, lsn, trace_key="events",
                                          window=(0.0, 2 * plane.dt), n_boot=args.n_boot, seed=1)
    for alpha in args.alphas:
        _, gmask = greedy_masks(STA, pvals, alpha, nr, nc)
        _synth_report(f"greedy (hs p<{alpha:g})", gmask, truth, alt, azi, n_rf, n_noise)
    print(f"\n(greedy: {args.n_boot} bootstraps, alpha sweep {args.alphas}, {time.time()-t0:.1f}s)")


# ----------------------------------------------------------------- real-session mode
def run_real(args):
    import pandas as pd
    from v1dd_metrics import nwb as vn
    from v1dd_metrics.config import MetricConfig
    from v1dd_metrics.families.receptive_fields import receptive_field_metrics
    from v1dd_metrics.pipeline import discover_sessions

    data_dir = args.data_dir or pathlib.Path("/data")
    sessions, _ = discover_sessions(str(data_dir))
    if args.sessions:
        sessions = sessions.iloc[:args.sessions]
    cfg = MetricConfig()
    methods = ["ours"] + [f"greedy_a{a:g}" for a in args.alphas]
    rows = []
    print(f"Real-session comparison on {len(sessions)} sessions (ON subfield)")
    print(f"greedy alpha sweep: {args.alphas}\n" + "=" * 64)

    for si, (_, sess) in enumerate(sessions.iterrows()):
        print(f"[{si+1}/{len(sessions)}] {sess['name']}")
        nwbfile, io = vn.open_session(sess["path"])
        try:
            stim = vn.load_stimulus_table(nwbfile)
            lsn_trials = vn.stimulus_trials(stim, "locally_sparse_noise")[0]
            if not len(lsn_trials):
                continue
            spont = vn.spontaneous_block(nwbfile)
            lsn = vn.load_lsn_template(nwbfile)
            grid_deg = lsn["grid_size_deg"]
            for plane_key in vn.list_planes(nwbfile):
                plane = vn.load_plane(nwbfile, plane_key, trace_types=("events", "dff"))
                if "events" not in plane.traces:
                    continue
                # our method
                _, ours_map = receptive_field_metrics(plane, lsn_trials, spont, lsn,
                                                      config=cfg, rng=np.random.default_rng(42))
                masks = {"ours": ours_map >= cfg.rf_frac_thresh}
                # greedy p-values once, thresholded at each alpha
                STA, pvals, nr, ncc = greedy_sta_pvals(
                    plane, lsn_trials, lsn, trace_key="events",
                    window=(0.0, 2 * plane.dt), n_boot=args.n_boot, seed=1)
                for a in args.alphas:
                    _, gm = greedy_masks(STA, pvals, a, nr, ncc)
                    masks[f"greedy_a{a:g}"] = gm
                mm = {m: mask_metrics(masks[m], grid_deg) for m in methods}
                for r in range(plane.n_rois):
                    row = {"session": sess["name"], "plane": plane.plane, "roi": int(plane.roi[r])}
                    for m in methods:
                        has, ncomp, largest, area = mm[m]
                        row[f"{m}_has_on"] = bool(has[r, 0])
                        row[f"{m}_single_on"] = ncomp[r, 0] == 1
                        row[f"{m}_area_on"] = area[r, 0]
                        row[f"{m}_largest_on"] = largest[r, 0]
                    rows.append(row)
                del plane
        finally:
            io.close()

    df = pd.DataFrame(rows)
    n = len(df)
    print(f"\n{n} ROIs\n" + "=" * 64)
    print(f"{'method':>12}  {'has_rf_on':>9}  {'single-comp':>11}  {'med_area':>9}  {'med_largest':>11}")
    for m in methods:
        has = df[f"{m}_has_on"]
        sub = df[has]
        sc = sub[f"{m}_single_on"].mean() if len(sub) else float("nan")
        ar = sub[f"{m}_area_on"].median() if len(sub) else float("nan")
        lf = sub[f"{m}_largest_on"].median() if len(sub) else float("nan")
        print(f"{m:>12}  {has.mean():>8.1%}  {sc:>10.1%}  {ar:>8.0f}  {lf:>11.3f}")
    print(f"\nagreement vs ours (ON detection):")
    for m in methods[1:]:
        both = (df["ours_has_on"] & df[f"{m}_has_on"]).mean()
        oo = (df["ours_has_on"] & ~df[f"{m}_has_on"]).mean()
        go = (~df["ours_has_on"] & df[f"{m}_has_on"]).mean()
        print(f"  {m:>12}: both {both:.1%} | ours-only {oo:.1%} | {m}-only {go:.1%}")
    out = REPO_ROOT / "docs" / "explorations" / "greedy_rf_comparison.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nsaved {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--synthetic", action="store_true", help="local ground-truth comparison")
    p.add_argument("--sessions", type=int, default=None, help="limit real sessions")
    p.add_argument("--data-dir", type=pathlib.Path, default=None)
    p.add_argument("--alphas", type=str, default="0.05,0.01,0.001",
                   help="comma-separated Holm-Sidak alphas to sweep (specificity knob)")
    p.add_argument("--n-boot", type=int, default=5000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    args.alphas = [float(a) for a in str(args.alphas).split(",") if a.strip()]
    if args.synthetic:
        run_synthetic(args)
    else:
        run_real(args)


if __name__ == "__main__":
    main()
