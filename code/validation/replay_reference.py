"""Replay a shipped asset's trial arrays through this code and compare the metrics.

The NWB input is only mounted in Code Ocean, but a shipped asset carries both the
trial-level arrays and the per-ROI metrics derived from them -- so the array-to-metric
half of the pipeline can be checked with no capsule and no mounted dataset. This is the
refactor gate; see docs/pipeline.md.

    python code/validation/replay_reference.py --asset <dir> [--n-fits 300]

Metrics needing the continuous trace (responsiveness, spontaneous rates) are not
replayable and are reported as skipped.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v1dd_metrics import config as cfg                      # noqa: E402
from v1dd_metrics.families import drifting_gratings as dg   # noqa: E402
from v1dd_metrics.families import receptive_fields as rfm    # noqa: E402
from v1dd_metrics.families import surround_suppression as ssm  # noqa: E402

RTOL = 1e-6
FAILS: list[str] = []


def ulp_noise(a: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """``a`` perturbed by up to half a float32 ULP -- the precision the archive stores."""
    a = np.asarray(a, float)
    ulp = np.nextafter(np.float32(np.abs(a)), np.float32(np.inf)).astype(float) - np.abs(a)
    return a + (rng.random(a.shape) - 0.5) * ulp


def curve(params: np.ndarray, directions: np.ndarray) -> np.ndarray:
    """Evaluate fitted von Mises parameters at each direction, flattened.

    ``params`` is (n_rois, n_sf, 6); unfitted entries stay NaN and drop out of any
    comparison.
    """
    p = np.asarray(params, float).reshape(-1, 6)
    out = np.full((p.shape[0], len(directions)), np.nan)
    ok = np.isfinite(p).all(axis=1)
    for i in np.flatnonzero(ok):
        out[i] = dg.vonmises_two_peak(directions, *p[i])
    return out.ravel()


def compare(name: str, got: np.ndarray, want: np.ndarray, floor: np.ndarray | None = None,
            rtol: float = RTOL, scale: float | None = None) -> None:
    """Report agreement between a replayed array and a shipped column.

    The archives store trials as float32, so a ratio metric inherits a relative error of
    about eps32/|metric| -- unbounded as the metric approaches zero. ``floor`` is the same
    metric recomputed from ULP-perturbed input, which measures that sensitivity directly;
    a disagreement no larger than the floor is storage precision, not a defect.

    ``scale`` floors the denominator for quantities with a meaningful zero. Degrees of
    visual angle are signed and centred on the screen, so a value of exactly 0 deg makes a
    pure relative error meaningless -- pass ``scale=1.0`` and the difference is reported
    against a degree instead. Ratio metrics take no scale on purpose: there the blow-up
    near zero is the real cancellation sensitivity and should stay visible.
    """
    got, want = np.asarray(got, float), np.asarray(want, float)
    gnan, wnan = np.isnan(got), np.isnan(want)
    both = ~gnan & ~wnan
    n_mismatch_nan = int((gnan != wnan).sum())

    def spread(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        if not m.any():
            return 0, 0.0, 0.0
        den = np.maximum(np.abs(b[m]), scale if scale is not None else 1e-300)
        rel = np.abs(a[m] - b[m]) / den
        return int((rel > rtol).sum()), float(np.percentile(rel, 99)), float(rel.max())

    n_ref, p99_ref, worst_ref = spread(got, want)

    if floor is None:
        ok = n_ref == 0 and n_mismatch_nan == 0
        detail = f"n={int(both.sum()):>6}  over_rtol={n_ref:>4}  worst_rel={worst_ref:.2e}"
    else:
        n_floor, p99_floor, worst_floor = spread(np.asarray(floor, float), got)
        # Judged on the count and the tail, not the extreme. The von Mises fit is
        # ill-conditioned enough that a single ROI settling in a different local minimum
        # gives an O(1) worst case on both sides, so `worst` is reported but does not
        # decide -- it would fail or pass on which ROIs the subsample happened to draw.
        ok = (n_mismatch_nan == 0 and n_ref <= 1.5 * n_floor + 5
              and p99_ref <= max(5 * p99_floor, rtol))
        detail = (f"n={int(both.sum()):>6}  over_rtol={n_ref:>4}|{n_floor:<4}"
                  f"  p99={p99_ref:.1e}|{p99_floor:.1e}"
                  f"  worst={worst_ref:.1e}|{worst_floor:.1e}")

    print(f"  {'PASS' if ok else 'FAIL'}  {name:<26} {detail}"
          + (f"  nan_mismatch={n_mismatch_nan}" if n_mismatch_nan else ""))
    if not ok:
        FAILS.append(name)


def check(name: str, cond: bool, detail: str = "") -> None:
    status = "PASS" if cond else "FAIL"
    print(f"  {status}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        FAILS.append(name)


def as_shipped(**kw):
    """The config that produced the reference asset.

    Today's defaults deliberately differ on three settings, so reproducing a shipped
    asset means reverting exactly those. Anything else that diverges is a regression.
    """
    return dataclasses.replace(cfg.DEFAULT_CONFIG,
                               lifetime_sparseness_over="trials",
                               zero_denominator_nan=False,
                               ssi_tuning_fit_includes_baseline=True, **kw)


def report_metric_changes(df: pd.DataFrame, tc: dict, is_valid: np.ndarray) -> None:
    """How far the corrected metrics move from what the reference asset shipped.

    Reported rather than asserted: these changes are intended, and the point is to make
    their size visible so no later diff is ambiguous about why a number moved.
    """
    print("\nintended metric changes (shipped convention -> current default)")
    dirs, sfs = tc["directions"], tc["spatial_frequencies"]
    old_cfg = as_shipped(fit_tuning_curves=False)
    new_cfg = dataclasses.replace(cfg.DEFAULT_CONFIG, fit_tuning_curves=False)

    for kind in ("dgw", "dgf"):
        ta = tc[kind + "_trials"].astype(np.float64)
        family = "drifting_gratings_" + ("windowed" if kind == "dgw" else "full")
        old = dg.dg_metrics_from_trials(ta, dirs, sfs, is_valid, family=family,
                                        config=old_cfg)
        new = dg.dg_metrics_from_trials(ta, dirs, sfs, is_valid, family=family,
                                        config=new_cfg)
        for metric in ("lifetime_sparseness", "osi", "dsi"):
            a, b = np.asarray(old[metric], float), np.asarray(new[metric], float)
            both = np.isfinite(a) & np.isfinite(b)
            moved = int((both & ~np.isclose(a, b, rtol=1e-9)).sum())
            newly_nan = int((np.isfinite(a) & ~np.isfinite(b)).sum())
            r = (np.corrcoef(a[both], b[both])[0, 1] if both.sum() > 2 else float("nan"))
            print(f"        {kind}_{metric:<22} median {np.nanmedian(a):.4f} -> "
                  f"{np.nanmedian(b):.4f}   moved {moved:>6}   "
                  f"newly NaN {newly_nan:>5}   r={r:.4f}")


def replay_gratings(df: pd.DataFrame, tc: dict, is_valid: np.ndarray, n_fits: int) -> None:
    """Drifting gratings, both types, from the shipped trial arrays."""
    dirs, sfs = tc["directions"], tc["spatial_frequencies"]
    no_fit = as_shipped(fit_tuning_curves=False)

    rng = np.random.default_rng(0)
    windowed = None
    for kind, prefix in (("dgw", "dgw_"), ("dgf", "dgf_")):
        ta = tc[kind + "_trials"].astype(np.float64)
        family = "drifting_gratings_" + ("windowed" if kind == "dgw" else "full")
        print(f"\n{kind}: metrics from {kind}_trials {ta.shape}")
        core = dg.dg_metrics_from_trials(ta, dirs, sfs, is_valid,
                                         family=family, config=no_fit)
        noise = dg.dg_metrics_from_trials(ulp_noise(ta, rng), dirs, sfs, is_valid,
                                          family=family, config=no_fit)
        if kind == "dgw":
            windowed = core
        for metric in ("dsi", "osi", "gosi", "pref_dir_mean", "lifetime_sparseness",
                       "preferred_dir", "preferred_sf"):
            compare(prefix + metric, core[metric], df[prefix + metric].to_numpy(),
                    floor=noise[metric])
        print(f"  SKIP  {prefix}frac_responsive_trials, {prefix}is_responsive "
              "(need the continuous trace)")

    if n_fits <= 0:
        return

    # Fits are ~0.1 s each, so a subsample rather than all 39,407.
    print(f"\nvon Mises fits on a {n_fits}-ROI subsample")
    rng = np.random.default_rng(0)
    for kind in ("dgw", "dgf"):
        shipped = tc[kind + "_params"]
        fitted_any = np.isfinite(shipped).any(axis=(1, 2)) & is_valid
        pool = np.flatnonzero(fitted_any)
        if pool.size == 0:
            check(f"{kind}_params subsample", False, "no fitted ROIs in the asset")
            continue
        idx = rng.choice(pool, size=min(n_fits, pool.size), replace=False)
        ta = tc[kind + "_trials"][idx].astype(np.float64)
        family = "drifting_gratings_" + ("windowed" if kind == "dgw" else "full")
        # Windowed self-selects its preferred SF; full field is told which one to fit.
        fit_sf = None if kind == "dgw" else windowed["pref_cond_index"][idx, 1]
        kw = dict(fit_sf_index=fit_sf, config=as_shipped(), family=family)
        core = dg.dg_metrics_from_trials(ta, dirs, sfs, is_valid[idx], **kw)
        # A bounded 6-parameter least squares amplifies input noise nonlinearly, so the
        # fit needs the same float32 floor the ratio metrics get.
        noise = dg.dg_metrics_from_trials(ulp_noise(ta, rng), dirs, sfs, is_valid[idx], **kw)
        # Compare the fitted CURVE, not the parameters. Several parameters are bounded at
        # zero and park there, so their relative error is against ~1e-30 and meaningless;
        # the curve is also the only thing ssi_tuning_fit reads.
        compare(f"{kind}_tuning_curve", curve(core["tuning_params"], dirs),
                curve(shipped[idx], dirs), floor=curve(noise["tuning_params"], dirs))


class _StandInPlane:
    """The identity a metrics frame needs, taken from the shipped table.

    Surround suppression is computed per plane in production, so the replay walks planes
    too -- running speeds have no ROI axis and key on the plane.
    """

    def __init__(self, rows: pd.DataFrame):
        first = rows.iloc[0]
        self.n_rois = len(rows)
        self.column, self.volume, self.plane = int(first["column"]), str(first["volume"]), int(first["plane"])
        self.roi = rows["roi"].to_numpy()
        self.mouse_id = str(first["mouse"]).lstrip("M")
        self.depth_um = float(first["depth_um"])
        self.roi_table = pd.DataFrame(
            {"pika_roi_confidence": rows["pika_roi_confidence"].to_numpy()})
        self.is_valid = rows["pika_roi_confidence"].to_numpy() > 0.5


def _ssi_frame(df: pd.DataFrame, tc: dict, config, trials=None, params=None) -> pd.DataFrame:
    """Recompute surround suppression for every plane from the shipped arrays.

    ``trials`` and ``params`` override the archives, which is how the noise floor is
    measured; both default to what the asset shipped.
    """
    dirs, sfs = tc["directions"], tc["spatial_frequencies"]
    trials = trials or {k: tc[k + "_trials"].astype(np.float64) for k in ("dgw", "dgf")}
    params = params or {k: tc[k + "_params"].astype(np.float64) for k in ("dgw", "dgf")}
    roi_plane = np.array([k.rsplit("_", 1)[0] for k in tc["roi_key"]])

    out = []
    for pi, pkey in enumerate(tc["plane_key"]):
        sel = np.flatnonzero(roi_plane == pkey)
        if not sel.size:
            continue
        rows = df.iloc[sel]
        plane = _StandInPlane(rows)
        results = {}
        for kind in ("dgw", "dgf"):
            core = dg.dg_metrics_from_trials(
                trials[kind][sel], dirs, sfs, plane.is_valid,
                family="drifting_gratings_" + ("windowed" if kind == "dgw" else "full"),
                config=dataclasses.replace(config, fit_tuning_curves=False))
            results[kind] = dg.DGResult(
                metrics=pd.DataFrame(), trial_responses=trials[kind][sel],
                trial_running_speeds=tc[kind + "_running"][pi].astype(np.float64),
                pref_cond_index=core["pref_cond_index"], tuning_params=params[kind][sel],
                dir_list=dirs, sf_list=sfs, blank_responses=tc[kind + "_blank"][sel])
        out.append(ssm.surround_suppression_metrics(
            results["dgw"], results["dgf"], plane, config=config,
            center=(float(rows["dgw_center_azimuth"].iloc[0]),
                    float(rows["dgw_center_elevation"].iloc[0]))))
    return pd.concat(out, ignore_index=True)


def replay_surround_suppression(df: pd.DataFrame, tc: dict) -> None:
    """The eight ssi_* columns, from the same arrays that produced them."""
    print("\nsurround suppression: from dg{w,f}_trials, _running and _params")
    shipped_cfg = as_shipped(fit_tuning_curves=False)
    got = _ssi_frame(df, tc, shipped_cfg)

    order = df["roi_key"].tolist()
    got = got.set_index("roi_key").reindex(order)
    check("recomputed rows cover the shipped ROIs", got.index.equals(pd.Index(order)),
          f"{len(got)} vs {len(order)}")

    # Perturb BOTH archives. ssi_tuning_fit reads tuning_params and never touches the
    # trials, so a floor built from trials alone measures nothing for it -- it came out
    # exactly 0, which is the tell that the wrong input was being varied.
    rng = np.random.default_rng(0)
    noise_trials = {k: ulp_noise(tc[k + "_trials"].astype(np.float64), rng)
                    for k in ("dgw", "dgf")}
    noise_params = {k: ulp_noise(tc[k + "_params"].astype(np.float64), rng)
                    for k in ("dgw", "dgf")}
    floor = _ssi_frame(df, tc, shipped_cfg, trials=noise_trials, params=noise_params
                       ).set_index("roi_key").reindex(order)

    for col in ssm.SSI_COLUMNS:
        compare(col, got[col].to_numpy(), df[col].to_numpy(),
                floor=floor[col].to_numpy())
    for col in ("dgw_center_azimuth", "dgw_center_elevation"):
        compare(col, got[col].to_numpy(), df[col].to_numpy())
    print("  ....  dgw_rf_* are replayed from rf_maps, in the next section")


def report_ssi_baseline_change(df: pd.DataFrame, tc: dict) -> None:
    """How far ssi_tuning_fit moves once the baseline is handled consistently.

    Driven from the *shipped* fit parameters, so the only difference between the two runs
    is the baseline convention -- the fit's own instability cannot contaminate it.
    """
    print("\nintended change: ssi_tuning_fit baseline (shipped fit parameters, so this is"
          " exact arithmetic)")
    old = _ssi_frame(df, tc, as_shipped(fit_tuning_curves=False))
    new = _ssi_frame(df, tc, dataclasses.replace(
        cfg.DEFAULT_CONFIG, fit_tuning_curves=False))
    a = old.set_index("roi_key")["ssi_tuning_fit"].reindex(df["roi_key"]).to_numpy()
    b = new.set_index("roi_key")["ssi_tuning_fit"].reindex(df["roi_key"]).to_numpy()
    both = np.isfinite(a) & np.isfinite(b)
    moved = int((both & ~np.isclose(a, b, rtol=1e-9)).sum())
    r = np.corrcoef(a[both], b[both])[0, 1] if both.sum() > 2 else float("nan")
    print(f"        ssi_tuning_fit          median {np.nanmedian(a):.4f} -> "
          f"{np.nanmedian(b):.4f}   moved {moved:>6} of {int(both.sum())}   "
          f"newly NaN {int((np.isfinite(a) & ~np.isfinite(b)).sum()):>5}   r={r:.4f}")


def _rf_centres(maps: np.ndarray, rf: dict, config) -> dict:
    """ON/OFF centres from the pre-threshold maps: threshold, then unweighted centroid.

    Invalid ROIs were already zeroed before the maps were saved, so nothing here needs
    ``is_valid`` -- a blank map means excluded, and falls out as no field.
    """
    m = np.asarray(maps, dtype=np.float64).copy()
    m[m < config.rf_frac_thresh] = 0.0
    mask = m > 0
    n_rows, n_cols = m.shape[2], m.shape[3]
    counts = mask.sum(axis=(2, 3))
    rows_ix = np.arange(n_rows)[None, None, :, None]
    cols_ix = np.arange(n_cols)[None, None, None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_row = (mask * rows_ix).sum(axis=(2, 3)) / counts
        mean_col = (mask * cols_ix).sum(axis=(2, 3)) / counts

    alt = rfm._rf_pixel_to_degrees(mean_row, np.asarray(rf["altitudes"], float),
                                   config.rf_center_scale_bug)
    azi = rfm._rf_pixel_to_degrees(mean_col, np.asarray(rf["azimuths"], float),
                                   config.rf_center_scale_bug)
    has = {"on": counts[:, 0] > 0, "off": counts[:, 1] > 0}
    out = {}
    for i, side in enumerate(("on", "off")):
        out[f"has_rf_{side}"] = has[side]
        out[f"azimuth_rf_{side}"] = np.where(has[side], azi[:, i], np.nan)
        out[f"altitude_rf_{side}"] = np.where(has[side], alt[:, i], np.nan)
    return out


def _containment_frame(df: pd.DataFrame, maps: np.ndarray, rf: dict, config
                       ) -> pd.DataFrame:
    """Window containment for every ROI, grouped by aperture centre.

    Coverage depends only on the window position, so it is computed once per distinct
    centre rather than once per ROI. The centre is constant within a column, so this is a
    handful of groups -- and the two sessions that record none fall into their own group,
    where ``_window_coverage`` returns None and every column stays NaN.
    """
    lsn = {"azimuths": rf["azimuths"], "altitudes": rf["altitudes"]}
    parts = []
    for (az, el), rows in df.groupby(["dgw_center_azimuth", "dgw_center_elevation"],
                                     dropna=False):
        idx = rows.index.to_numpy()
        parts.append(ssm.window_containment(rows, maps[idx], lsn, (az, el),
                                            config=config))
    return pd.concat(parts).reindex(df.index)


def replay_receptive_fields(df: pd.DataFrame, rf: dict) -> None:
    """RF centres and window containment, both from the shipped pre-threshold maps."""
    print("\nreceptive fields: centres and window containment, from rf_maps")
    cfg_ship = as_shipped()
    maps = np.asarray(rf["rf_maps"], dtype=np.float64)
    rng = np.random.default_rng(0)
    noisy = ulp_noise(maps, rng)

    got, floor = _rf_centres(maps, rf, cfg_ship), _rf_centres(noisy, rf, cfg_ship)
    for side in ("on", "off"):
        agree = float((got[f"has_rf_{side}"] == (df[f"has_rf_{side}"].to_numpy() > 0)).mean())
        check(f"has_rf_{side} reproduced", agree == 1.0, f"{agree:.4%}")
        for axis in ("azimuth", "altitude"):
            col = f"{axis}_rf_{side}"
            # degrees of visual angle: signed and centred on zero
            compare(col, got[col], df[col].to_numpy(), floor=floor[col],
                    scale=1.0)

    cont = _containment_frame(df, maps, rf, cfg_ship)
    cont_floor = _containment_frame(df, noisy, rf, cfg_ship)
    for col in ssm.CONTAINMENT_COLUMNS:
        # distances are degrees; overlaps are fractions in [0, 1]
        compare(col, cont[col].to_numpy(), df[col].to_numpy(),
                floor=cont_floor[col].to_numpy(), scale=1.0)

    # The two sessions with no recorded aperture: every containment column must be NaN,
    # and nothing else may be.
    no_centre = ~np.isfinite(df["dgw_center_azimuth"].to_numpy())
    if no_centre.any():
        blank = cont.loc[no_centre, ssm.CONTAINMENT_COLUMNS]
        check("sessions with no aperture centre get NaN containment",
              blank.isna().all().all(), f"{int(no_centre.sum())} ROIs")


def check_exact_relationships(df: pd.DataFrame, rf: dict, tc: dict,
                              is_valid: np.ndarray) -> None:
    """Three relationships that hold exactly, so any drift is unambiguous."""
    print("\nexact relationships")

    thresh = cfg.DEFAULT_CONFIG.rf_frac_thresh
    maps = rf["rf_maps"]
    for i, side in enumerate(("on", "off")):
        got = (maps[:, i] >= thresh).any(axis=(1, 2)) & is_valid
        want = df[f"has_rf_{side}"].to_numpy() > 0
        agree = float((got == want).mean())
        check(f"rf_maps at {thresh} reproduces has_rf_{side}", agree == 1.0,
              f"{agree:.4%} of ROIs")

    no_fit = as_shipped(fit_tuning_curves=False)
    for kind in ("dgw", "dgf"):
        ta = tc[kind + "_trials"].astype(np.float64)
        core = dg.dg_metrics_from_trials(
            ta, tc["directions"], tc["spatial_frequencies"], is_valid,
            family="drifting_gratings_" + ("windowed" if kind == "dgw" else "full"),
            config=no_fit)
        frac = df[f"{kind}_frac_responsive_trials"].to_numpy()
        ri = np.arange(len(ta))
        n = np.isfinite(ta[ri, core["pref_dir_idx"], core["pref_sf_idx"], :]).sum(axis=1)
        finite = np.isfinite(frac) & (n > 0)
        prod = frac[finite] * n[finite]
        resid = np.abs(prod - np.round(prod))
        check(f"{kind} frac x n is an integer", float(resid.max()) < 1e-9,
              f"n={int(finite.sum())}, max residual {float(resid.max()):.2e}, "
              f"n in {sorted(set(n[finite].tolist()))}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--asset", required=True, help="a shipped asset directory")
    ap.add_argument("--n-fits", type=int, default=300,
                    help="ROIs to re-fit for the von Mises comparison; 0 to skip")
    args = ap.parse_args(argv)

    asset = Path(args.asset)
    (feather,) = list(asset.glob("*.feather"))
    df = pd.read_feather(feather)
    tc = dict(np.load(next(asset.glob("tuning_curves*.npz")), allow_pickle=True))
    rf = dict(np.load(next(asset.glob("rf_maps*.npz")), allow_pickle=True))

    print(f"asset   : {asset.name}")
    print(f"table   : {df.shape[0]} ROIs x {df.shape[1]} columns")
    print("paired values are  ours-vs-reference | float32-noise-floor")

    for name, arr in (("tuning_curves", tc), ("rf_maps", rf)):
        same = np.array_equal(np.asarray(arr["roi_key"], dtype=object),
                              df["roi_key"].to_numpy(dtype=object))
        check(f"{name} roi_key aligns with the table", same)
    if FAILS:
        return 1

    is_valid = df["pika_roi_confidence"].to_numpy() > 0.5
    check("is_valid from pika_roi_confidence", int((~is_valid).sum()) == 1038,
          f"{int((~is_valid).sum())} invalid")

    replay_gratings(df, tc, is_valid, args.n_fits)
    replay_surround_suppression(df, tc)
    replay_receptive_fields(df, rf)
    check_exact_relationships(df, rf, tc, is_valid)
    report_metric_changes(df, tc, is_valid)
    report_ssi_baseline_change(df, tc)

    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)} check(s): " + ", ".join(FAILS))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
