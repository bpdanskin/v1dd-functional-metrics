"""L0 event inference quality diagnostics for V1DD functional metrics.

Loads the stimulus_metrics.parquet asset and computes diagnostic statistics:
spontaneous event rate distribution, zero-inflation, depth dependence,
reliability (events vs dF/F), and SNR vs event rate.

Figures land in docs/figures/ by default (evt_spont_rate_distribution.png,
evt_depth_dependence.png, evt_reliability_comparison.png,
evt_snr_vs_event_rate.png).

Usage:
    python event_diagnostics.py                              # auto-find asset
    python event_diagnostics.py --asset path/to/asset_dir    # explicit asset
    python event_diagnostics.py --fig-dir path/to/figs       # custom output
"""
import argparse, pathlib, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


def find_latest_asset(results_dir):
    """Return the most recent asset directory containing stimulus_metrics.parquet."""
    candidates = sorted(
        (d for d in results_dir.iterdir()
         if d.is_dir() and (d / "stimulus_metrics.parquet").exists()),
        key=lambda d: d.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        sys.exit(f"No asset with stimulus_metrics.parquet found in {results_dir}")
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--asset", type=pathlib.Path, default=None,
                        help="Asset directory (auto-detected from results/ if omitted)")
    parser.add_argument("--fig-dir", type=pathlib.Path,
                        default=REPO_ROOT / "docs" / "figures",
                        help="Where to write figures (default: docs/figures/)")
    args = parser.parse_args()

    if args.asset is None:
        args.asset = find_latest_asset(REPO_ROOT / "results")
    parquet = args.asset / "stimulus_metrics.parquet"
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {parquet}")
    df = pd.read_parquet(parquet)
    n_rois = len(df)
    print(f"  {n_rois:,} ROIs, {len(df.columns)} columns\n")

    sr = df["spont_rate"]
    depth = df["depth_um"]
    snr = df["snr"]

    # ----- 1. Spontaneous event rate -----
    print("=" * 72)
    print("1. SPONTANEOUS EVENT RATE")
    print("=" * 72)
    print(f"  Median: {sr.median():.6f}  Mean: {sr.mean():.6f}")
    print(f"  5-95th: [{sr.quantile(0.05):.6f}, {sr.quantile(0.95):.6f}]")
    print(f"  Exactly zero: {(sr == 0).sum():,} ROIs ({(sr == 0).mean():.4f})")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(sr, bins=200, color="steelblue", edgecolor="none")
    axes[0].set_xlabel("spont_rate (mean event amplitude)")
    axes[0].set_ylabel("ROI count")
    axes[0].set_title("Spontaneous event rate distribution")
    axes[0].axvline(sr.median(), color="red", ls="--",
                    label=f"median={sr.median():.4f}")
    axes[0].axvline(sr.quantile(0.95), color="orange", ls="--",
                    label=f"95th={sr.quantile(0.95):.4f}")
    axes[0].legend(fontsize=8)

    sr_pos = sr[sr > 0]
    axes[1].hist(np.log10(sr_pos), bins=200, color="steelblue", edgecolor="none")
    axes[1].set_xlabel("log10(spont_rate)")
    axes[1].set_ylabel("ROI count")
    axes[1].set_title("Spontaneous event rate (log scale)")
    axes[1].axvline(np.log10(sr.median()), color="red", ls="--", label="median")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(args.fig_dir / "evt_spont_rate_distribution.png", dpi=150)
    plt.close(fig)

    # ----- 2. Zero-inflation -----
    print("\n" + "=" * 72)
    print("2. ZERO-INFLATION")
    print("=" * 72)
    pref_cols = {
        "ni_pref_response":   "Natural Images (118)",
        "ni12_pref_response":  "Natural Images (12)",
        "nm_pref_response":   "Natural Movie",
    }
    for col, label in pref_cols.items():
        z = (df[col] == 0).mean()
        print(f"  {label:<35s}  exact zero: {z:.4f}")

    # ----- 3. Depth dependence -----
    print("\n" + "=" * 72)
    print("3. DEPTH DEPENDENCE")
    print("=" * 72)
    r_spont, _ = stats.pearsonr(depth, sr)
    r_snr, _ = stats.pearsonr(depth, snr)
    print(f"  r(depth, spont_rate) = {r_spont:.4f}")
    print(f"  r(depth, snr)        = {r_snr:.4f}")

    depth_stats = df.groupby("depth_um").agg(
        n=("roi", "count"),
        spont_rate_median=("spont_rate", "median"),
        snr_median=("snr", "median"),
        snr_25=("snr", lambda x: x.quantile(0.25)),
        snr_75=("snr", lambda x: x.quantile(0.75)),
        spont_rate_25=("spont_rate", lambda x: x.quantile(0.25)),
        spont_rate_75=("spont_rate", lambda x: x.quantile(0.75)),
    ).reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ax = axes[0]
    ax.errorbar(depth_stats["depth_um"], depth_stats["spont_rate_median"],
                yerr=[depth_stats["spont_rate_median"] - depth_stats["spont_rate_25"],
                      depth_stats["spont_rate_75"] - depth_stats["spont_rate_median"]],
                fmt="o-", color="steelblue", markersize=4, capsize=2)
    ax.set_xlabel("Depth (um)")
    ax.set_ylabel("spont_rate (median, IQR)")
    ax.set_title(f"Spontaneous event rate vs depth\nr={r_spont:.3f}")

    ax = axes[1]
    ax.errorbar(depth_stats["depth_um"], depth_stats["snr_median"],
                yerr=[depth_stats["snr_median"] - depth_stats["snr_25"],
                      depth_stats["snr_75"] - depth_stats["snr_median"]],
                fmt="o-", color="darkorange", markersize=4, capsize=2)
    ax.set_xlabel("Depth (um)")
    ax.set_ylabel("SNR (median, IQR)")
    ax.set_title(f"SNR vs depth\nr={r_snr:.3f}")

    ax = axes[2]
    ax.bar(depth_stats["depth_um"], depth_stats["n"], width=14, color="gray", alpha=0.7)
    ax.set_xlabel("Depth (um)")
    ax.set_ylabel("ROI count")
    ax.set_title("ROI count by depth")

    plt.tight_layout()
    fig.savefig(args.fig_dir / "evt_depth_dependence.png", dpi=150)
    plt.close(fig)

    # ----- 4. Reliability: events vs dF/F -----
    print("\n" + "=" * 72)
    print("4. RELIABILITY: EVENTS vs dF/F")
    print("=" * 72)
    # reliability = dF/F (primary), reliability_events = events (companion)
    reliability_pairs = [
        ("ni_reliability_events", "ni_reliability", "Natural Images (118)"),
        ("ni12_reliability_events", "ni12_reliability", "Natural Images (12)"),
        ("nm_reliability_events", "nm_reliability", "Natural Movie"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i, (ev_col, dff_col, label) in enumerate(reliability_pairs):
        ev = df[ev_col]
        dff = df[dff_col]
        mask = np.isfinite(ev) & np.isfinite(dff)
        ev_c, dff_c = ev[mask], dff[mask]

        r_val, _ = stats.pearsonr(ev_c, dff_c)
        sign_disagree = ((ev_c > 0) != (dff_c > 0)).mean()
        ev_lower = (ev_c < dff_c).mean()
        print(f"  {label}:  r={r_val:.3f}  sign-disagree={sign_disagree:.1%}  "
              f"events<dF/F={ev_lower:.0%}")

        ax = axes[i]
        rng = np.random.default_rng(42)
        idx = rng.choice(len(ev_c), min(10000, len(ev_c)), replace=False)
        ax.scatter(dff_c.iloc[idx], ev_c.iloc[idx], s=1, alpha=0.15,
                   color="steelblue", rasterized=True)
        lims = [min(ev_c.min(), dff_c.min()), max(ev_c.max(), dff_c.max())]
        ax.plot(lims, lims, "r--", lw=1, label="identity")
        ax.set_xlabel(dff_col)
        ax.set_ylabel(ev_col)
        ax.set_title(f"{label}\nr={r_val:.3f}, sign disagree={sign_disagree:.1%}")
        ax.set_aspect("equal", adjustable="datalim")
        ax.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(args.fig_dir / "evt_reliability_comparison.png", dpi=150)
    plt.close(fig)

    # ----- 5. SNR vs event rate -----
    print("\n" + "=" * 72)
    print("5. SNR vs EVENT RATE")
    print("=" * 72)
    r_snr_sr, _ = stats.pearsonr(snr, sr)
    rho, _ = stats.spearmanr(snr, sr)
    print(f"  Pearson r = {r_snr_sr:.4f}  Spearman rho = {rho:.4f}")

    low_q = sr[snr < snr.quantile(0.25)]
    high_q = sr[snr > snr.quantile(0.75)]
    print(f"  Bottom-quartile SNR  median event rate: {low_q.median():.6f}")
    print(f"  Top-quartile SNR    median event rate:  {high_q.median():.6f}")

    snr_edges = np.unique(np.percentile(snr, np.linspace(0, 100, 21)))
    snr_mids, snr_med_sr, snr_ns = [], [], []
    for lo, hi in zip(snr_edges[:-1], snr_edges[1:]):
        m = (snr >= lo) & (snr < hi)
        if m.sum() == 0:
            continue
        snr_mids.append((lo + hi) / 2)
        snr_med_sr.append(sr[m].median())
        snr_ns.append(m.sum())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    rng = np.random.default_rng(42)
    idx = rng.choice(n_rois, min(15000, n_rois), replace=False)
    ax.scatter(snr.iloc[idx], sr.iloc[idx], s=1, alpha=0.1,
               color="steelblue", rasterized=True)
    ax.set_xlabel("SNR (spectral, from dF/F)")
    ax.set_ylabel("spont_rate (mean event amplitude)")
    ax.set_title(f"SNR vs event rate\nPearson r={r_snr_sr:.3f}, Spearman rho={rho:.3f}")

    ax = axes[1]
    ax.plot(snr_mids, snr_med_sr, "o-", color="darkorange", markersize=5)
    ax.set_xlabel("SNR (bin center)")
    ax.set_ylabel("Median spont_rate in bin")
    ax.set_title("Binned median event rate by SNR")
    for x, y, n in zip(snr_mids, snr_med_sr, snr_ns):
        ax.annotate(f"n={n}", (x, y), fontsize=6, textcoords="offset points",
                    xytext=(0, 6), ha="center")

    plt.tight_layout()
    fig.savefig(args.fig_dir / "evt_snr_vs_event_rate.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved 4 figures to {args.fig_dir}")


if __name__ == "__main__":
    main()
