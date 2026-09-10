"""Extract natural image stimuli from NWB and pair with population response.

Saves the 118 + 12 stimulus images and the population-level preferred-image
distribution. Produces a figure showing the most- and least-activating images
ranked by population response.

Usage (inside capsule):
    PYTHONPATH=/code/src python -u validation/extract_ni_images.py

Outputs:
    docs/ni_stimulus_images.npz  — pixel data + population response summary
    docs/figures/ni_preferred_images.png
"""
import argparse, pathlib, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "src"))

from v1dd_metrics import nwb as vn
from v1dd_metrics.pipeline import discover_sessions


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fig-dir", type=pathlib.Path,
                        default=REPO_ROOT / "docs" / "figures")
    parser.add_argument("--data-dir", type=pathlib.Path, default=None)
    parser.add_argument("--asset-dir", type=pathlib.Path, default=None,
                        help="Asset with condition_means.npz (default: latest in results/)")
    args = parser.parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    # --- load stimulus images from the first session's NWB ---
    data_dir = args.data_dir or pathlib.Path("/data")
    sessions, _ = discover_sessions(str(data_dir))
    nwbfile, io = vn.open_session(sessions.iloc[0]["path"])
    try:
        ni_images, ni_frames = vn._images_to_array(nwbfile.stimulus["natural_images"])
        print(f"Natural images: {ni_images.shape} (indices {ni_frames[0]}..{ni_frames[-1]})")

        nm_images, nm_frames = vn._images_to_array(nwbfile.stimulus["natural_movie"])
        print(f"Natural movie frames: {nm_images.shape}")
    finally:
        io.close()

    # --- load population responses from the asset ---
    if args.asset_dir:
        asset = args.asset_dir
    else:
        results = sorted(REPO_ROOT.glob("results/409828_V1DD_functional_metrics_*"))
        asset = results[-1] if results else None

    cm = None
    df = None
    if asset and (asset / "condition_means.npz").exists():
        cm = np.load(asset / "condition_means.npz", allow_pickle=True)
        print(f"\nCondition means from {asset.name}")
        print(f"  ni_mean: {cm['ni_mean'].shape}")

    if asset and (asset / "stimulus_metrics.parquet").exists():
        import pandas as pd
        df = pd.read_parquet(asset / "stimulus_metrics.parquet")
        print(f"  {len(df):,} ROIs")

    # --- population response per image ---
    if cm is not None:
        ni_mean = cm["ni_mean"]          # (n_rois, 118)
        ni_ids = cm["ni_images"]         # (118,)
        ni12_ids = cm["ni12_images"]     # (12,)

        pop_mean = np.nanmean(ni_mean, axis=0)  # mean across ROIs
        pop_median = np.nanmedian(ni_mean, axis=0)
        rank = np.argsort(pop_mean)[::-1]        # most activating first

        # How many ROIs prefer each image?
        if df is not None:
            pref_counts = df["ni_pref_img"].value_counts().reindex(ni_ids, fill_value=0)
        else:
            pref_counts = None

        print(f"\nPopulation mean response per image:")
        print(f"  Range: [{pop_mean.min():.6f}, {pop_mean.max():.6f}]")
        print(f"  Top 5 images (by pop mean): {ni_ids[rank[:5]]}")
        print(f"  Bottom 5 images (by pop mean): {ni_ids[rank[-5:]]}")

    # --- figure: top and bottom images ---
    if cm is not None:
        n_show = 10
        fig, axes = plt.subplots(4, n_show, figsize=(2.2 * n_show, 9.5))

        for col in range(n_show):
            # Top row: most activating
            idx = rank[col]
            img_id = ni_ids[idx]
            img = ni_images[img_id] if img_id < len(ni_images) else np.zeros((10, 10))
            ax = axes[0, col]
            ax.imshow(img, cmap="gray" if img.ndim == 2 else None)
            ax.set_title(f"#{col+1}\nimg {img_id}", fontsize=8)
            ax.axis("off")

            ax = axes[1, col]
            ax.bar(0, pop_mean[idx], color="steelblue", width=0.6)
            if pref_counts is not None:
                ax.bar(1, pref_counts.iloc[idx], color="darkorange", width=0.6)
                ax.set_xticks([0, 1])
                ax.set_xticklabels(["resp", "pref"], fontsize=7)
            ax.tick_params(labelsize=7)

            # Bottom rows: least activating
            idx_b = rank[-(col + 1)]
            img_id_b = ni_ids[idx_b]
            img_b = ni_images[img_id_b] if img_id_b < len(ni_images) else np.zeros((10, 10))
            ax = axes[2, col]
            ax.imshow(img_b, cmap="gray" if img_b.ndim == 2 else None)
            ax.set_title(f"#{118 - col}\nimg {img_id_b}", fontsize=8)
            ax.axis("off")

            ax = axes[3, col]
            ax.bar(0, pop_mean[idx_b], color="steelblue", width=0.6)
            if pref_counts is not None:
                ax.bar(1, pref_counts.iloc[idx_b], color="darkorange", width=0.6)
                ax.set_xticks([0, 1])
                ax.set_xticklabels(["resp", "pref"], fontsize=7)
            ax.tick_params(labelsize=7)

        axes[0, 0].set_ylabel("Most activating", fontsize=10)
        axes[2, 0].set_ylabel("Least activating", fontsize=10)
        fig.suptitle(
            f"Natural images ranked by population mean response\n"
            f"blue = mean response across {ni_mean.shape[0]:,} ROIs, "
            f"orange = # ROIs preferring this image",
            fontsize=11, fontweight="bold")
        plt.tight_layout()
        fig.savefig(args.fig_dir / "ni_preferred_images.png", dpi=150,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"\nSaved figure to {args.fig_dir / 'ni_preferred_images.png'}")

    # --- save npz ---
    npz_path = args.fig_dir.parent / "ni_stimulus_images.npz"
    payload = {
        "ni_images": ni_images,             # (118, H, W) or (118, H, W, 3)
        "ni_image_ids": np.array(ni_frames),
        "nm_images": nm_images,
    }
    if cm is not None:
        payload["ni_pop_mean_response"] = pop_mean
        payload["ni_pop_median_response"] = pop_median
        payload["ni_rank_by_pop_mean"] = rank
        payload["ni12_image_ids"] = cm["ni12_images"]
    np.savez_compressed(npz_path, **payload)
    sz = sum(v.nbytes for v in payload.values())
    print(f"Saved {npz_path} ({sz / 1e6:.1f} MB uncompressed)")


if __name__ == "__main__":
    main()
