"""Prototype: split-half cross-validated OSI/DSI vs our naive (single-pass) OSI/DSI.

Our naive selectivity (``families/drifting_gratings.py`` ``dg_metrics_from_trials``)
picks the preferred condition by argmax over the trial-averaged grid and then
measures OSI/DSI at that *same* condition on the *same* trials. Picking and
measuring on one draw biases selectivity **upward**: the argmax rides the noise
peak, and the orthogonal/null directions ride their noise troughs. The bias is
worst for weakly-driven cells, where noise dominates the ranking.

The **split-half cross-validation** (Jun Zhuang's ``*_shw_*`` columns) removes
that coupling: on each iteration, split the trials in two, pick the preferred
condition on one half, measure OSI/DSI/gOSI on the *other* half, and average over
many iterations. The held-out measurement cannot inherit the selection's noise,
so the estimate is (near-)unbiased. Two by-products fall out:

  * the spread of the estimate across iterations -- how noisy a 4+4-trial split
    is, i.e. how much averaging is needed to trust the number; and
  * a reliability count -- how many iterations had the held-out preferred
    response reach a z-score threshold above the blank (Jun's ``iter_num_2``),
    a candidate for folding reliability *into* OSI/DSI rather than shipping a
    separate column.

Runs fully locally on the asset -- no capsule. Inputs:
  results/<asset>/tuning_curves.npz     dgw_trials/dgf_trials (N,12,2,8), blanks,
                                        directions, spatial_frequencies
  results/<asset>/stimulus_metrics.parquet   the naive dgw_/dgf_ osi/dsi/gosi to
                                        cross-check the reproduction against

Output: docs/explorations/crossval_osi_dsi.csv (per ROI, both stimuli, naive vs
cross-validated OSI/DSI/gOSI + reliability counts + join keys).

Usage
-----
  python code/validation/crossval_osi_dsi.py                # default asset, 200 iters
  python code/validation/crossval_osi_dsi.py --n-iter 500 --stimulus dgw
"""
import argparse
import pathlib
import sys
import time

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "src"))

from v1dd_metrics.common import _ratio


def _naive_selectivity(mean_tr, dir_list):
    """Naive OSI/DSI/gOSI/pref-index from a trial-averaged grid ``(N,n_dir,n_sf)``.

    Mirrors ``dg_metrics_from_trials``: finite argmax over the (dir x sf) grid is
    the preferred condition; selectivity is read at the preferred SF. Returns a
    dict of ``(N,)`` arrays plus the preferred dir/sf indices.
    """
    n_rois, n_dir, n_sf = mean_tr.shape
    roi_ix = np.arange(n_rois)
    flat = np.where(np.isfinite(mean_tr), mean_tr, -np.inf).reshape(n_rois, -1)
    all_nan = ~np.isfinite(mean_tr).any(axis=(1, 2))
    k = flat.argmax(axis=1)
    pref_dir_idx, pref_sf_idx = np.divmod(k, n_sf)

    tuning = mean_tr[roi_ix, :, pref_sf_idx]                 # (N, n_dir) at pref SF
    pref = tuning[roi_ix, pref_dir_idx]
    null_r = tuning[roi_ix, (pref_dir_idx + 6) % n_dir]
    orth_r = 0.5 * (tuning[roi_ix, (pref_dir_idx + 3) % n_dir]
                    + tuning[roi_ix, (pref_dir_idx - 3) % n_dir])
    osi = _ratio(pref - orth_r, pref + orth_r, zero_to_nan=True)
    dsi = _ratio(pref - null_r, pref + null_r, zero_to_nan=True)

    theta = np.deg2rad(dir_list.astype(float))
    L_norm = tuning.sum(axis=1)
    L_ori = tuning @ np.exp(2j * theta)
    with np.errstate(invalid="ignore", divide="ignore"):
        gosi = np.abs(np.where(L_norm != 0, L_ori / np.where(L_norm != 0, L_norm, 1.0),
                               L_ori))
    for arr in (osi, dsi, gosi):
        arr[all_nan] = np.nan
    return {"osi": osi, "dsi": dsi, "gosi": gosi, "pref": pref,
            "pref_dir_idx": np.where(all_nan, -1, pref_dir_idx),
            "pref_sf_idx": np.where(all_nan, -1, pref_sf_idx)}


def crossval_selectivity(trials, blank, dir_list, *, n_iter=200, seed=0):
    """Split-half cross-validated OSI/DSI/gOSI over ``n_iter`` random trial splits.

    ``trials`` is ``(N, n_dir, n_sf, n_trials)`` NaN-padded; ``blank`` is
    ``(N, n_blank)``. Each iteration: split the trial slots 4/4, pick the
    preferred (dir, sf) on half A, measure selectivity on half B at that
    condition. Reliability z is the held-out preferred response minus the blank
    mean over the blank std (floored at the population median positive std, so
    the ~13% of cells with an all-silent blank stay defined). Returns per-ROI
    means, cross-iteration stds, preferred-direction stability against the naive
    pick, and the fraction of iterations with z >= 2 and z >= 3.
    """
    N, n_dir, n_sf, n_trials = trials.shape
    roi_ix = np.arange(N)
    theta = np.deg2rad(dir_list.astype(float))
    exp2 = np.exp(2j * theta)

    naive = _naive_selectivity(np.nanmean(trials, axis=3), dir_list)
    naive_dir = naive["pref_dir_idx"]

    blank_mean = np.nanmean(blank, axis=1)
    blank_std = np.nanstd(blank, axis=1)
    floor = float(np.median(blank_std[blank_std > 0]))
    blank_std = np.where(blank_std > floor, blank_std, floor)

    # streaming accumulators over iterations
    keys = ("osi", "dsi", "gosi")
    s1 = {k: np.zeros(N) for k in keys}                      # sum for the mean
    s2 = {k: np.zeros(N) for k in keys}                      # sum of squares for std
    cnt = {k: np.zeros(N) for k in keys}                     # finite-iteration count
    dir_match = np.zeros(N)
    z2 = np.zeros(N)
    z3 = np.zeros(N)

    rng = np.random.default_rng(seed)
    for _ in range(n_iter):
        perm = rng.permutation(n_trials)
        a, b = perm[:n_trials // 2], perm[n_trials // 2:]
        mA = np.nanmean(trials[:, :, :, a], axis=3)          # (N, n_dir, n_sf)
        mB = np.nanmean(trials[:, :, :, b], axis=3)

        flatA = np.where(np.isfinite(mA), mA, -np.inf).reshape(N, -1)
        valid = np.isfinite(mA).any(axis=(1, 2))
        kA = flatA.argmax(axis=1)
        pdir, psf = np.divmod(kA, n_sf)

        tunB = mB[roi_ix, :, psf]                            # (N, n_dir) held-out
        prefB = tunB[roi_ix, pdir]
        nullB = tunB[roi_ix, (pdir + 6) % n_dir]
        orthB = 0.5 * (tunB[roi_ix, (pdir + 3) % n_dir]
                       + tunB[roi_ix, (pdir - 3) % n_dir])
        vals = {
            "osi": _ratio(prefB - orthB, prefB + orthB, zero_to_nan=True),
            "dsi": _ratio(prefB - nullB, prefB + nullB, zero_to_nan=True),
        }
        L_norm = np.nansum(tunB, axis=1)
        L_ori = np.nan_to_num(tunB, nan=0.0) @ exp2
        with np.errstate(invalid="ignore", divide="ignore"):
            vals["gosi"] = np.abs(np.where(L_norm != 0, L_ori / np.where(L_norm != 0, L_norm, 1.0),
                                           L_ori))
        for k in keys:
            v = vals[k]
            fin = np.isfinite(v) & valid
            s1[k][fin] += v[fin]
            s2[k][fin] += v[fin] ** 2
            cnt[k][fin] += 1.0

        z = (prefB - blank_mean) / blank_std
        z2 += valid & np.isfinite(z) & (z >= 2.0)
        z3 += valid & np.isfinite(z) & (z >= 3.0)
        dir_match += valid & (pdir == naive_dir)

    out = {"naive_osi": naive["osi"], "naive_dsi": naive["dsi"], "naive_gosi": naive["gosi"]}
    for k in keys:
        n = np.where(cnt[k] > 0, cnt[k], np.nan)
        mean = s1[k] / n
        var = np.maximum(s2[k] / n - mean ** 2, 0.0)
        out[f"cv_{k}"] = mean
        out[f"cv_{k}_std"] = np.sqrt(var)
    out["pref_dir_stability"] = dir_match / n_iter
    out["reliab_frac_z2"] = z2 / n_iter
    out["reliab_frac_z3"] = z3 / n_iter
    return out


def _summary(tag, res, is_resp):
    """Print naive-vs-CV bias, CV noise, and reliability, overall and by responsiveness."""
    def line(name, m):
        for k in ("osi", "dsi", "gosi"):
            nv, cv = res["naive_" + k][m], res["cv_" + k][m]
            good = np.isfinite(nv) & np.isfinite(cv)
            bias = np.nanmean((nv - cv)[good])
            std = np.nanmean(res[f"cv_{k}_std"][m])
            print(f"    {name:<14} {k:>4}  naive {np.nanmean(nv[good]):.3f}  "
                  f"cv {np.nanmean(cv[good]):.3f}  bias {bias:+.3f}  "
                  f"cv-iter-std {std:.3f}")
    print(f"\n{tag}  (N={res['naive_osi'].size}, {int(is_resp.sum())} responsive)")
    line("all cells", np.ones_like(is_resp, bool))
    line("responsive", is_resp.astype(bool))
    line("non-resp", ~is_resp.astype(bool))
    print(f"    pref-dir stability (mean frac of splits matching naive pick): "
          f"{np.nanmean(res['pref_dir_stability']):.2f}  "
          f"(responsive {np.nanmean(res['pref_dir_stability'][is_resp.astype(bool)]):.2f})")
    print(f"    reliability: median frac splits z>=2 {np.nanmedian(res['reliab_frac_z2']):.2f}, "
          f"z>=3 {np.nanmedian(res['reliab_frac_z3']):.2f}  "
          f"(responsive z>=2 {np.nanmedian(res['reliab_frac_z2'][is_resp.astype(bool)]):.2f})")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--asset", type=str,
                   default="results/409828_V1DD_functional_metrics_2026-09-11_01-17-17")
    p.add_argument("--stimulus", type=str, default="both", choices=["dgw", "dgf", "both"])
    p.add_argument("--n-iter", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    import pandas as pd
    asset = REPO_ROOT / args.asset
    z = np.load(asset / "tuning_curves.npz", allow_pickle=True)
    meta = pd.read_parquet(asset / "stimulus_metrics.parquet")
    dir_list = z["directions"]
    stims = ["dgw", "dgf"] if args.stimulus == "both" else [args.stimulus]

    keep = ["mouse", "column", "volume", "plane", "roi", "depth_um", "snr"]
    df = meta[keep].copy()

    for st in stims:
        print("=" * 72)
        t0 = time.time()
        res = crossval_selectivity(z[f"{st}_trials"], z[f"{st}_blank"], dir_list,
                                   n_iter=args.n_iter, seed=args.seed)
        # cross-check the naive reproduction against the shipped parquet columns
        for k in ("osi", "dsi", "gosi"):
            shipped = meta[f"{st}_{k}"].to_numpy()
            mine = res[f"naive_{k}"]
            good = np.isfinite(shipped) & np.isfinite(mine)
            md = np.max(np.abs(shipped[good] - mine[good]))
            print(f"  reproduce {st}_{k}: max|diff| vs parquet = {md:.2e}")
        _summary(st, res, meta[f"{st}_is_responsive"].to_numpy())
        print(f"  ({args.n_iter} iterations, {time.time() - t0:.1f}s)")
        for col, vals in res.items():
            df[f"{st}_{col}"] = vals

    out = REPO_ROOT / "docs" / "explorations" / "crossval_osi_dsi.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print("\nsaved", out)


if __name__ == "__main__":
    main()
