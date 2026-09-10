"""Simulate L0 event detection at V1DD's 6 Hz from Huang et al. ground truth.

Downloads Emx1-s (GCaMP6s) neurons from the Allen oephys dataset, downsamples
fluorescence from ~158 Hz to 30 Hz and 6 Hz, runs OASIS deconvolution, and
compares detected events against known spike times.

Figures land in docs/figures/ by default (evt_detection_vs_ap_count.png,
evt_roc_curves.png, evt_example_traces.png).

Usage:
    python simulate_6hz.py                          # defaults
    python simulate_6hz.py --fig-dir path/to/figs   # custom output
    python simulate_6hz.py --data-dir /tmp/oephys    # custom download cache
    python simulate_6hz.py --n-neurons 5             # fewer neurons for speed
"""
import argparse, os, sys, urllib.request, pathlib
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from oasis import oasisAR1

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]

BUCKET = "https://allen-paper-supplements.s3.us-west-2.amazonaws.com"
EMXS_IDS = [
    102081, 102086, 102102, 102422, 102436, 102437, 102443, 102451,
    102511, 102516, 102519, 102525, 102529, 102858, 102860, 102862,
    102865, 102871, 102875, 102890, 102906,
    103750, 103753, 103808, 103813, 103842, 103844, 103849, 103852,
    103856, 103863, 103864,
]
TAU_GCAMPS = 1.5  # GCaMP6s decay time constant, seconds


def download_neuron(specimen_id, data_dir):
    path = data_dir / f"{specimen_id}_processed.h5"
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    url = (f"{BUCKET}/huang_published_2021/processed_data/"
           f"Emx1-s_highzoom/{specimen_id}_processed.h5")
    print(f"  Downloading {specimen_id}...", end=" ", flush=True)
    urllib.request.urlretrieve(url, path)
    print(f"{path.stat().st_size / 1e6:.0f} MB")
    return path


def load_neuron(path):
    with h5py.File(path, "r") as f:
        return {
            "f_cell": f["f_cell"][0].astype(np.float64),
            "f_np": f["f_np"][0].astype(np.float64),
            "spk": f["spk"][0].astype(np.float64),
            "iSpk": f["iSpk"][0].astype(np.int64),
            "iFrames": f["iFrames"][0].astype(np.int64),
            "dto": float(f["dto"][0]),
            "dte": float(f["dte"][0]),
            "depth": float(f["depth"][0, 0]),
        }


def compute_dff(f_cell, percentile=10):
    baseline = np.percentile(f_cell, percentile)
    if baseline <= 0:
        baseline = np.mean(f_cell)
    return (f_cell - baseline) / baseline


def downsample_trace(trace, factor):
    n = len(trace) // factor * factor
    return trace[:n].reshape(-1, factor).mean(axis=1)


def bin_spikes(spike_times_sec, bin_edges):
    counts, _ = np.histogram(spike_times_sec, bins=bin_edges)
    return counts


def run_deconv(dff, dt, tau=TAU_GCAMPS):
    gamma = np.exp(-dt / tau)
    noise_std = np.median(np.abs(np.diff(dff))) / 0.6745
    lam = noise_std * (1 - gamma)
    c, s = oasisAR1(dff.astype(np.float64), gamma, lam=lam, s_min=0.0)
    return c, s, gamma, lam


def detection_analysis(s_inferred, spike_counts_per_bin, fp_threshold=0.01):
    is_event_bin = spike_counts_per_bin > 0
    n_neg = np.sum(~is_event_bin)
    n_pos = np.sum(is_event_bin)
    if n_pos == 0 or n_neg == 0:
        return {}

    thresholds = np.sort(np.unique(s_inferred[s_inferred > 0]))
    if len(thresholds) == 0:
        return {"tp_rate": np.array([0.0]), "fp_rate": np.array([0.0]),
                "det_by_ap": {}}

    tp_rates, fp_rates = [], []
    for thr in thresholds:
        detected = s_inferred >= thr
        fp_rates.append(np.sum(detected & ~is_event_bin) / n_neg)
        tp_rates.append(np.sum(detected & is_event_bin) / n_pos)

    fp_rates = np.array(fp_rates)
    tp_rates = np.array(tp_rates)

    det_by_ap = {}
    valid = fp_rates <= fp_threshold
    if np.any(valid):
        best_idx = np.argmax(tp_rates[valid])
        thr_at_fp = thresholds[np.where(valid)[0][best_idx]]
        detected_at_thr = s_inferred >= thr_at_fp
        for n_ap in range(1, 6):
            mask = spike_counts_per_bin == n_ap
            if mask.sum() > 0:
                det_by_ap[n_ap] = np.mean(detected_at_thr[mask])
            else:
                det_by_ap[n_ap] = np.nan

    return {"tp_rate": tp_rates, "fp_rate": fp_rates,
            "det_by_ap": det_by_ap, "thresholds": thresholds}


TARGET_RATES = {"native": None, "30Hz": 30.0, "6Hz": 6.0}


def analyze_one_neuron(data, specimen_id):
    dto_native = data["dto"]
    rate_native = 1.0 / dto_native
    f_cell = data["f_cell"]
    spike_indices = data["iSpk"]
    dte = data["dte"]
    iFrames = data["iFrames"]

    spike_times_sec = spike_indices * dte
    ophys_start_sec = iFrames[0] * dte
    ophys_end_sec = iFrames[-1] * dte
    valid_spikes = spike_times_sec[
        (spike_times_sec >= ophys_start_sec) & (spike_times_sec <= ophys_end_sec)
    ]
    valid_spikes_rel = valid_spikes - ophys_start_sec
    dff_native = compute_dff(f_cell)

    results = {}
    for label, target_rate in TARGET_RATES.items():
        if target_rate is None:
            dff, dt, actual_rate = dff_native, dto_native, rate_native
        else:
            ds_factor = max(1, int(round(rate_native / target_rate)))
            dff = downsample_trace(dff_native, ds_factor)
            dt = dto_native * ds_factor
            actual_rate = 1.0 / dt

        n_frames = len(dff)
        bin_edges = np.arange(n_frames + 1) * dt
        spike_counts = bin_spikes(valid_spikes_rel, bin_edges)

        try:
            c, s, gamma, lam = run_deconv(dff, dt)
        except Exception as e:
            print(f"  {specimen_id} @ {label}: deconv failed: {e}")
            continue

        det = detection_analysis(s, spike_counts)
        results[label] = {
            "rate": actual_rate, "dt": dt, "n_frames": n_frames,
            "duration": n_frames * dt,
            "n_spikes_total": int(spike_counts.sum()),
            "n_spike_bins": int((spike_counts > 0).sum()),
            "spike_counts": spike_counts,
            "dff": dff, "c": c, "s": s,
            "gamma": gamma, "lam": lam,
            "detection": det,
        }
    return results


def make_figures(all_results, fig_dir):
    fig_dir.mkdir(parents=True, exist_ok=True)

    # --- Detection probability vs AP count ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    for ax_i, (label, color) in enumerate(
        zip(["native", "30Hz", "6Hz"], ["#2a7d9a", "#e8993a", "#c74040"])
    ):
        ax = axes[ax_i]
        ap_counts, det_means, det_sems = [], [], []
        for n_ap in range(1, 6):
            vals = [
                r["detection"]["det_by_ap"][n_ap]
                for r in all_results.values() if label in r
                and "det_by_ap" in r[label].get("detection", {})
                and n_ap in r[label]["detection"]["det_by_ap"]
                and not np.isnan(r[label]["detection"]["det_by_ap"][n_ap])
            ]
            if vals:
                ap_counts.append(n_ap)
                det_means.append(np.mean(vals))
                det_sems.append(np.std(vals) / np.sqrt(len(vals)))
        if ap_counts:
            ax.bar(ap_counts, det_means, yerr=det_sems, color=color,
                   alpha=0.8, capsize=4, edgecolor="white", linewidth=0.5)

        rate_str = next(
            (f"{r[label]['rate']:.0f} Hz" for r in all_results.values() if label in r),
            label,
        )
        ax.set_title(rate_str, fontsize=13, fontweight="bold")
        ax.set_xlabel("# APs per bin", fontsize=11)
        if ax_i == 0:
            ax.set_ylabel("Detection probability\n(at 1% FP)", fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.set_xticks(range(1, 6))
        ax.axhline(0.5, color="gray", ls="--", alpha=0.3, lw=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Event detection: GCaMP6s (Emx1-s) ground truth at 3 sampling rates",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "evt_detection_vs_ap_count.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- ROC curves ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    for ax_i, (label, color) in enumerate(
        zip(["native", "30Hz", "6Hz"], ["#2a7d9a", "#e8993a", "#c74040"])
    ):
        ax = axes[ax_i]
        fp_common = np.linspace(0, 0.05, 200)
        tp_interp = []
        for results in all_results.values():
            if label not in results:
                continue
            det = results[label].get("detection", {})
            fp = det.get("fp_rate", np.array([]))
            tp = det.get("tp_rate", np.array([]))
            if len(fp) > 0:
                ax.plot(fp, tp, color=color, alpha=0.2, lw=0.7)
            if len(fp) > 1:
                order = np.argsort(fp)
                tp_interp.append(np.interp(fp_common, fp[order], tp[order],
                                           left=0, right=1))
        if tp_interp:
            ax.plot(fp_common, np.mean(tp_interp, axis=0),
                    color=color, lw=2.5, label="Mean")

        rate_str = next(
            (f"{r[label]['rate']:.0f} Hz" for r in all_results.values() if label in r),
            label,
        )
        ax.set_title(rate_str, fontsize=13, fontweight="bold")
        ax.set_xlabel("False positive rate", fontsize=11)
        if ax_i == 0:
            ax.set_ylabel("True positive rate", fontsize=11)
        ax.set_xlim(0, 0.05)
        ax.set_ylim(0, 1.0)
        ax.axvline(0.01, color="gray", ls=":", alpha=0.4, lw=0.8)
        ax.legend(fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("ROC curves: event detection at 3 sampling rates (FP 0-5%)",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "evt_roc_curves.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Example traces ---
    best_sid = max(all_results, key=lambda s: all_results[s].get(
        "native", {}).get("n_spikes_total", 0))
    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    for ax_i, (label, color) in enumerate(
        zip(["native", "30Hz", "6Hz"], ["#2a7d9a", "#e8993a", "#c74040"])
    ):
        ax = axes[ax_i]
        r = all_results[best_sid][label]
        t = np.arange(r["n_frames"]) * r["dt"]
        mask = t <= 30

        ax.plot(t[mask], r["dff"][mask], color="#999999", alpha=0.5, lw=0.5,
                label="dF/F")
        ax.plot(t[mask], r["c"][mask], color=color, lw=1.0, label="Denoised")

        s_mask = mask & (r["s"] > 0)
        ax.scatter(t[s_mask], r["s"][s_mask] + r["c"][s_mask],
                   color=color, s=15, zorder=5, marker="v", label="Inferred")

        sc = r["spike_counts"]
        gt_mask = mask & (sc > 0)
        ax.scatter(t[gt_mask], np.full(gt_mask.sum(), r["dff"][mask].max() * 1.05),
                   color="black", s=8, marker="|", zorder=4, label="True spikes")

        ax.set_ylabel(f"{r['rate']:.0f} Hz\ndF/F", fontsize=10)
        ax.legend(fontsize=8, loc="upper right", ncol=4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Time (s)", fontsize=11)
    fig.suptitle(f"Example traces: neuron {best_sid} at 3 sampling rates (first 30s)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(fig_dir / "evt_example_traces.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Saved 3 figures to {fig_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fig-dir", type=pathlib.Path,
                        default=REPO_ROOT / "docs" / "figures",
                        help="Where to write figures (default: docs/figures/)")
    parser.add_argument("--data-dir", type=pathlib.Path,
                        default=REPO_ROOT / "scratch" / "oephys_data",
                        help="Where to cache downloaded HDF5 files")
    parser.add_argument("--n-neurons", type=int, default=10,
                        help="Number of Emx1-s neurons to analyse (max 32)")
    args = parser.parse_args()

    args.data_dir.mkdir(parents=True, exist_ok=True)
    n_neurons = min(args.n_neurons, len(EMXS_IDS))
    print(f"Downloading {n_neurons} Emx1-s (GCaMP6s) neurons...")

    all_results = {}
    for sid in EMXS_IDS[:n_neurons]:
        path = download_neuron(sid, args.data_dir)
        data = load_neuron(path)
        print(f"  {sid}: {1/data['dto']:.0f} Hz, {len(data['f_cell'])} frames, "
              f"{len(data['iSpk'])} spikes "
              f"({len(data['iSpk'])/(len(data['f_cell'])*data['dto']):.2f} Hz)")
        all_results[sid] = analyze_one_neuron(data, sid)

    # Print detection table
    print("\n" + "=" * 70)
    print("DETECTION PROBABILITY AT 1% FALSE POSITIVE RATE")
    print("=" * 70)
    for label in ["native", "30Hz", "6Hz"]:
        det_by_ap_all = {n: [] for n in range(1, 6)}
        for results in all_results.values():
            if label not in results:
                continue
            dba = results[label].get("detection", {}).get("det_by_ap", {})
            for n_ap in range(1, 6):
                if n_ap in dba and not np.isnan(dba[n_ap]):
                    det_by_ap_all[n_ap].append(dba[n_ap])

        rate_info = next(
            (f"({r[label]['rate']:.1f} Hz)" for r in all_results.values() if label in r),
            "",
        )
        print(f"\n{label} {rate_info}:")
        for n_ap in range(1, 6):
            vals = det_by_ap_all[n_ap]
            if vals:
                print(f"  {n_ap} AP: {np.mean(vals):.3f} +/- "
                      f"{np.std(vals)/np.sqrt(len(vals)):.3f} (n={len(vals)})")
            else:
                print(f"  {n_ap} AP: no data")

    make_figures(all_results, args.fig_dir)


if __name__ == "__main__":
    main()
