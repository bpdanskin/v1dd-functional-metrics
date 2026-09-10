"""Compare receptive field maps across trace types and response windows.

Two dimensions: trace type (dF/F vs events) and response window (1, 2, 4
imaging frames). At 6 Hz each frame is ~165 ms; the LSN stimulus presents
a new pattern every 250 ms. A 4-frame window (660 ms) spans 2.6 stimulus
presentations, smearing spatial precision. A 1-frame window (165 ms) stays
within one presentation but captures less of the calcium transient.

Compares field contiguity (connected-component count, largest-component
fraction) and has_rf rates across conditions. Saves unthresholded maps
for all conditions to an .npz archive for visual inspection.

Designed to run inside the capsule where NWB files are available.

Usage (inside capsule):
    python rf_trace_comparison.py
    python rf_trace_comparison.py --sessions 3   # fewer for speed

Outputs:
    docs/figures/evt_rf_window_comparison.png  — summary figure
    docs/rf_trace_comparison.csv               — per-ROI metrics table
    docs/rf_comparison_maps.npz                — unthresholded maps, all conditions
"""
import argparse, pathlib, sys, time
from types import MappingProxyType

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "src"))

from v1dd_metrics import nwb as vn
from v1dd_metrics.config import MetricConfig
from v1dd_metrics.families.receptive_fields import receptive_field_metrics
from v1dd_metrics.pipeline import discover_sessions


FRAME_COUNTS = [1, 2, 4]
TRACE_TYPES = ["dff", "events"]


def contiguity_metrics(rf_map, threshold):
    """Per-ROI contiguity of above-threshold RF pixels.

    Returns (n_components, largest_frac) each (n_rois, 2) for ON/OFF.
    """
    n_rois, _, n_rows, n_cols = rf_map.shape
    n_comp = np.zeros((n_rois, 2), dtype=int)
    largest_frac = np.zeros((n_rois, 2), dtype=float)

    for pol in range(2):
        for r in range(n_rois):
            mask = rf_map[r, pol] >= threshold
            total = mask.sum()
            if total == 0:
                continue
            labeled, nc = ndimage.label(mask)
            n_comp[r, pol] = nc
            sizes = ndimage.sum(mask, labeled, range(1, nc + 1))
            largest_frac[r, pol] = sizes.max() / total

    return n_comp, largest_frac


def make_config(trace_key, n_frames):
    """MetricConfig with LSN trace type and response window overridden."""
    tt = dict(MetricConfig().trace_type)
    tt["locally_sparse_noise"] = trace_key
    return MetricConfig(trace_type=MappingProxyType(tt),
                        lsn_response_frames=n_frames)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fig-dir", type=pathlib.Path,
                        default=REPO_ROOT / "docs" / "figures")
    parser.add_argument("--sessions", type=int, default=None,
                        help="Limit to first N sessions (default: all)")
    parser.add_argument("--data-dir", type=pathlib.Path, default=None,
                        help="NWB data directory (default: capsule /data)")
    args = parser.parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    data_dir = args.data_dir or pathlib.Path("/data")
    sessions, _ = discover_sessions(str(data_dir))
    if args.sessions:
        sessions = sessions.iloc[:args.sessions]

    conditions = [(tk, nf) for tk in TRACE_TYPES for nf in FRAME_COUNTS]
    configs = {(tk, nf): make_config(tk, nf) for tk, nf in conditions}
    # Accumulate unthresholded maps per condition: list of (n_rois, 2, 8, 14)
    map_parts = {cond: [] for cond in conditions}
    roi_keys = []
    lsn_grid = None

    print(f"Comparing RF maps: {len(TRACE_TYPES)} trace types x "
          f"{len(FRAME_COUNTS)} windows on {len(sessions)} sessions")
    t0 = time.time()

    # --- report LSN timing from the first session ---
    first_nwb, first_io = vn.open_session(sessions.iloc[0]["path"])
    try:
        stim = vn.load_stimulus_table(first_nwb)
        lsn_rows = stim[stim["stim_name"] == "locally_sparse_noise"].sort_values("start_time")
        if len(lsn_rows) > 1:
            isi = np.diff(lsn_rows["start_time"].values)
            durations = (lsn_rows["stop_time"] - lsn_rows["start_time"]).values
            print(f"\nLSN stimulus timing (session {sessions.iloc[0]['name']}):")
            print(f"  {len(lsn_rows)} presentations")
            print(f"  Duration per frame: median {np.median(durations):.4f} s, "
                  f"range [{np.min(durations):.4f}, {np.max(durations):.4f}]")
            print(f"  ISI (onset-to-onset): median {np.median(isi):.4f} s, "
                  f"range [{np.min(isi):.4f}, {np.max(isi):.4f}]")
    finally:
        first_io.close()

    all_rows = []
    for si, (_, sess_row) in enumerate(sessions.iterrows()):
        print(f"\n[{si+1}/{len(sessions)}] {sess_row['name']}")
        nwbfile, io = vn.open_session(sess_row["path"])
        try:
            stim = vn.load_stimulus_table(nwbfile)
            lsn_trials = vn.stimulus_trials(stim, "locally_sparse_noise")[0]
            if not len(lsn_trials):
                print("  no LSN trials, skipping")
                continue
            spont = vn.spontaneous_block(nwbfile)
            lsn = vn.load_lsn_template(nwbfile)
            if lsn_grid is None:
                lsn_grid = lsn

            for plane_key in vn.list_planes(nwbfile):
                plane = vn.load_plane(nwbfile, plane_key,
                                      trace_types=("events", "dff"))
                if "dff" not in plane.traces or "events" not in plane.traces:
                    continue

                results = {}
                for tk, nf in conditions:
                    cfg = configs[(tk, nf)]
                    rng = np.random.default_rng(42)
                    out, rfmap = receptive_field_metrics(
                        plane, lsn_trials, spont, lsn, config=cfg, rng=rng)
                    nc, lf = contiguity_metrics(rfmap, cfg.rf_frac_thresh)
                    results[(tk, nf)] = (out, nc, lf, rfmap)

                for cond in conditions:
                    map_parts[cond].append(results[cond][3])

                col = sess_row.get("column", "")
                vol = sess_row.get("volume", "")
                for r in range(plane.n_rois):
                    roi_keys.append(
                        f"c{col}_v{vol}_p{plane.plane}_{plane.roi[r]}")
                    row = {
                        "session": sess_row["name"],
                        "plane": plane.plane,
                        "roi": int(plane.roi[r]),
                        "dt": plane.dt,
                    }
                    for tk, nf in conditions:
                        tag = f"{tk}_{nf}f"
                        out, nc, lf, _ = results[(tk, nf)]
                        row[f"has_rf_on_{tag}"] = bool(out["has_rf_on"].iloc[r])
                        row[f"has_rf_off_{tag}"] = bool(out["has_rf_off"].iloc[r])
                        row[f"n_comp_on_{tag}"] = nc[r, 0]
                        row[f"n_comp_off_{tag}"] = nc[r, 1]
                        row[f"largest_frac_on_{tag}"] = lf[r, 0]
                        row[f"largest_frac_off_{tag}"] = lf[r, 1]
                    all_rows.append(row)
                del plane
        finally:
            io.close()

    if not all_rows:
        sys.exit("No planes with both traces and LSN stimulus found")

    import pandas as pd
    df = pd.DataFrame(all_rows)
    n = len(df)
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"RESULTS: {n} ROIs across {df['session'].nunique()} sessions  "
          f"({elapsed/60:.1f} min)")
    print(f"dt range: [{df['dt'].min():.5f}, {df['dt'].max():.5f}] s")
    print(f"{'='*70}")

    # --- summary table ---
    print(f"\n{'':>12}  {'has_rf_on':>10}  {'single_comp':>11}  "
          f"{'med_largest':>11}  {'has_rf_off':>10}")
    for tk in TRACE_TYPES:
        for nf in FRAME_COUNTS:
            tag = f"{tk}_{nf}f"
            n_on = df[f"has_rf_on_{tag}"].sum()
            has_both = df[f"has_rf_on_{tag}"]
            sub = df[has_both]
            sc = (sub[f"n_comp_on_{tag}"] == 1).mean() if len(sub) else float("nan")
            ml = sub[f"largest_frac_on_{tag}"].median() if len(sub) else float("nan")
            n_off = df[f"has_rf_off_{tag}"].sum()
            window_ms = nf * df["dt"].median() * 1000
            label = f"{tk:>6} {nf}f ({window_ms:.0f}ms)"
            print(f"  {label:>20}  {n_on:>5} ({n_on/n:>5.1%})  "
                  f"{sc:>10.1%}  {ml:>10.3f}  "
                  f"{n_off:>5} ({n_off/n:>5.1%})")

    # --- figure: has_rf rate and contiguity across conditions ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = {"dff": "steelblue", "events": "darkorange"}
    markers = {"dff": "o", "events": "s"}

    for col_idx, pol in enumerate(("on", "off")):
        # Panel 1: has_rf rate vs window width
        ax = axes[0, col_idx]
        for tk in TRACE_TYPES:
            rates = []
            for nf in FRAME_COUNTS:
                tag = f"{tk}_{nf}f"
                rates.append(df[f"has_rf_{pol}_{tag}"].mean())
            ax.plot(FRAME_COUNTS, rates, f"-{markers[tk]}", color=colors[tk],
                    label=tk, markersize=8, lw=2)
        ax.set_xlabel("Response window (imaging frames)")
        ax.set_ylabel(f"Fraction with RF {pol.upper()}")
        ax.set_title(f"RF {pol.upper()} detection rate")
        ax.set_xticks(FRAME_COUNTS)
        dt_med = df["dt"].median()
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(FRAME_COUNTS)
        ax2.set_xticklabels([f"{nf * dt_med * 1000:.0f} ms" for nf in FRAME_COUNTS])
        ax2.set_xlabel("Window duration")
        ax.legend()

        # Panel 2: single-component fraction vs window width
        ax = axes[1, col_idx]
        for tk in TRACE_TYPES:
            sc_rates = []
            for nf in FRAME_COUNTS:
                tag = f"{tk}_{nf}f"
                has = df[f"has_rf_{pol}_{tag}"]
                sub = df[has]
                if len(sub):
                    sc_rates.append((sub[f"n_comp_{pol}_{tag}"] == 1).mean())
                else:
                    sc_rates.append(float("nan"))
            ax.plot(FRAME_COUNTS, sc_rates, f"-{markers[tk]}", color=colors[tk],
                    label=tk, markersize=8, lw=2)
        ax.set_xlabel("Response window (imaging frames)")
        ax.set_ylabel(f"Single-component fraction (RF {pol.upper()})")
        ax.set_title(f"RF {pol.upper()} contiguity")
        ax.set_xticks(FRAME_COUNTS)
        ax.legend()

    fig.suptitle("RF mapping: trace type x response window at 6 Hz\n"
                 f"({n:,} ROIs, {df['session'].nunique()} sessions, "
                 f"LSN ISI ~250 ms, dt ~{dt_med*1000:.0f} ms)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(args.fig_dir / "evt_rf_window_comparison.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nSaved figure to {args.fig_dir / 'evt_rf_window_comparison.png'}")

    out_path = args.fig_dir.parent / "rf_trace_comparison.csv"
    df.to_csv(out_path, index=False)
    print(f"Saved comparison table to {out_path}")

    # --- save unthresholded maps for visual inspection ---
    npz_data = {"roi_key": np.array(roi_keys)}
    for (tk, nf), parts in map_parts.items():
        if parts:
            npz_data[f"rf_maps_{tk}_{nf}f"] = np.concatenate(parts, axis=0)
    # Include the stimulus grid so maps can be plotted in degrees
    npz_data["altitudes"] = np.asarray(lsn_grid["altitudes"], dtype=np.float64)
    npz_data["azimuths"] = np.asarray(lsn_grid["azimuths"], dtype=np.float64)
    npz_data["frame_counts"] = np.array(FRAME_COUNTS)
    npz_data["trace_types"] = np.array(TRACE_TYPES)
    npz_path = args.fig_dir.parent / "rf_comparison_maps.npz"
    np.savez_compressed(npz_path, **npz_data)
    print(f"Saved unthresholded maps to {npz_path} "
          f"({sum(v.nbytes for v in npz_data.values()) / 1e6:.1f} MB uncompressed)")


if __name__ == "__main__":
    main()
