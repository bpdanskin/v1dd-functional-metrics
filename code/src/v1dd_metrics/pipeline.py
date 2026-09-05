"""Session and plane orchestration: load, run each family, assemble the asset.

One run walks every session in the mounted asset, computes each family for every imaging
plane, and writes the per-ROI table plus the arrays that do not fit a table. See
docs/pipeline.md for the stage order and docs/outputs.md for the asset layout.
"""

from __future__ import annotations

import dataclasses
import json
import platform
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np
import pandas as pd

from . import nwb as vn
from . import provenance as prov
from .config import DEFAULT_CONFIG, REFERENCE_CONFIG, MetricConfig
from .families import drifting_gratings as dgm
from .families import natural_images as nim
from .families import natural_movie as nmm
from .families import receptive_fields as rfm
from .families import roi_quality as rqm
from .families import surround_suppression as ssm
from .schema import to_output_schema

FAMILIES = ["roi_summary", "drifting_gratings_full", "drifting_gratings_windowed",
            "surround_suppression", "natural_images", "natural_images_12",
            "natural_movie", "rf_metrics"]

#: Column prefix each family contributes to the wide table. Families whose published
#: names already carry the stimulus keep an empty prefix.
PREFIX = {"drifting_gratings_full": "dgf_", "drifting_gratings_windowed": "dgw_",
          "surround_suppression": "", "roi_summary": "", "natural_images": "ni_",
          "natural_images_12": "ni12_", "natural_movie": "nm_", "rf_metrics": ""}

ID_COLS = ["roi_unique_id", "roi_key", "mouse", "column", "volume", "plane", "roi",
           "depth_um", "pika_roi_confidence"]
KEYS = ["column", "volume", "plane", "roi"]

_TUNING_PARTS = ("trials", "blank", "params", "running", "plane_key")


class Accumulator:
    """Per-plane results, gathered across the whole run."""

    def __init__(self) -> None:
        self.parts: dict[str, list] = {f: [] for f in FAMILIES}
        self.rf_maps: list[np.ndarray] = []
        self.lsn_grid: Optional[dict] = None
        self.tuning = {k: {p: [] for p in _TUNING_PARTS} for k in ("dgw", "dgf")}
        self.tuning_axes: Optional[tuple] = None
        self.cond_means: dict[str, list] = {"natural_images": [], "natural_images_12": []}
        self.cond_ids: dict[str, Any] = {}
        self.plane_log: list[dict] = []
        self.failures: list[dict] = []

    def tables(self) -> dict[str, pd.DataFrame]:
        return {f: pd.concat(v, ignore_index=True) for f, v in self.parts.items() if v}


def discover_sessions(functional_dir: str,
                      session_filter: Optional[Sequence[tuple]] = None
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Readable sessions in the asset, and the full inventory including failures.

    ``session_filter`` is a list of ``(column, volume)`` pairs restricting the run; the
    subset is recorded in provenance so a partial asset cannot pass for a complete one.
    """
    paths = vn.find_sessions(functional_dir)
    inventory = pd.DataFrame([vn.peek_session(p) for p in paths])
    inventory["format"] = [vn.nwb_format(p) for p in paths]
    sessions = inventory[inventory["error"].isna()].reset_index(drop=True)
    if session_filter is not None:
        want = {(int(c), str(v)) for c, v in session_filter}
        sessions = sessions[[(int(c), str(v)) in want
                             for c, v in zip(sessions["column"], sessions["volume"])]
                            ].reset_index(drop=True)
    return sessions, inventory


def aperture_centers(sessions: pd.DataFrame, config: MetricConfig):
    """Per-session windowed-grating aperture centre, filled from the column median.

    A pre-pass rather than a fix-up: filling needs every session in a column, and the
    per-plane loop cannot know them when it reaches the first one. Returns
    ``(WindowCenters, read_failures)``; an unreadable session yields an *absent* centre,
    which the imputation already handles, never a wrong one.
    """
    observed, failures = {}, []
    for _, row in sessions.iterrows():
        key = (int(row["column"]), str(row["volume"]))
        try:
            with vn.session(row["path"]) as nwbfile:
                trials, blank = vn.stimulus_trials(
                    vn.load_stimulus_table(nwbfile), "drifting_gratings_windowed",
                    vn.DG_PARAM_COLUMNS)
            observed[key] = dgm.window_center(trials.loc[~blank])
        except Exception as exc:                                        # noqa: BLE001
            failures.append({"session": row["name"], "column": key[0], "volume": key[1],
                             "error": repr(exc)[:200]})
            observed[key] = (np.nan, np.nan)
    return dgm.infer_window_centers(observed, config=config), failures


def process_plane(plane, ctx: dict, acc: Accumulator, config: MetricConfig,
                  seed: int) -> None:
    """Run every family for one imaging plane and append to ``acc``.

    ``rng`` is a factory returning a freshly seeded generator per family, so the call
    order below does not change any number.
    """
    def rng() -> np.random.Generator:
        return np.random.default_rng(seed)

    # Receptive fields first: surround suppression reports how much of each field the
    # grating aperture covered, so the maps must exist before it runs.
    rf_df, rf_map = rfm.receptive_field_metrics(
        plane, ctx["lsn_trials"], ctx["spont"], ctx["lsn"], config=config, rng=rng())
    acc.parts["rf_metrics"].append(rf_df)
    acc.rf_maps.append(rf_map)

    trials, blank = ctx["dg_trials"]["windowed"]
    dgw = dgm.drifting_gratings_metrics(
        plane, trials, blank, ctx["spont"], ctx["running"], dg_type="windowed",
        config=config, rng=rng())
    acc.parts["drifting_gratings_windowed"].append(dgw.metrics)

    # Full field fits only the spatial frequency surround suppression reads, which is
    # windowed's preferred one per ROI.
    trials, blank = ctx["dg_trials"]["full"]
    dgf = dgm.drifting_gratings_metrics(
        plane, trials, blank, ctx["spont"], ctx["running"], dg_type="full",
        fit_sf_index=dgw.pref_cond_index[:, 1], config=config, rng=rng())
    acc.parts["drifting_gratings_full"].append(dgf.metrics)

    containment = ssm.window_containment(rf_df, rf_map, ctx["lsn"], ctx["center"],
                                         config=config)
    acc.parts["surround_suppression"].append(ssm.surround_suppression_metrics(
        dgw, dgf, plane, config=config, containment=containment,
        center=ctx["center"], center_inferred=ctx["center_inferred"]))

    # Locomotion spans both grating types and the spontaneous block, so it is its own
    # family rather than columns bolted onto one of them.
    acc.parts["roi_summary"].append(rqm.roi_summary_metrics(
        plane, dgw, dgf, ctx["spont"], ctx["running"], config=config))

    for key, res in (("dgw", dgw), ("dgf", dgf)):
        part = acc.tuning[key]
        part["trials"].append(res.trial_responses.astype(np.float32))
        part["blank"].append(res.blank_responses.astype(np.float32))
        part["params"].append(res.tuning_params.astype(np.float32))
        part["running"].append(res.trial_running_speeds.astype(np.float32))
        # Running speeds have no ROI axis, so they key on the plane. Taken as roi_key's
        # prefix rather than rebuilt, so the two cannot drift.
        part["plane_key"].append(res.metrics["roi_key"].iloc[0].rsplit("_", 1)[0]
                                 if len(res.metrics) else "")
    if acc.tuning_axes is None:
        acc.tuning_axes = (dgw.dir_list, dgw.sf_list)

    for fam in ("natural_images", "natural_images_12"):
        df, means = nim.natural_images_metrics(
            plane, ctx["ni_trials"][fam], ctx["spont"], ns_type=fam, config=config,
            rng=rng())
        acc.parts[fam].append(df)
        if means is not None:
            acc.cond_means[fam].append(means[0])
            acc.cond_ids.setdefault(fam, means[1])

    acc.parts["natural_movie"].append(nmm.natural_movie_metrics(
        plane, ctx["nm_trials"], ctx["spont"], config=config, rng=rng()))

    acc.plane_log.append({"session": ctx["session_name"], "column": plane.column,
                          "volume": plane.volume, "plane": plane.plane,
                          "n_rois": plane.n_rois, "depth_um": plane.depth_um,
                          "dt": round(plane.dt, 6)})


def process_session(row: pd.Series, centers, acc: Accumulator, config: MetricConfig,
                    seed: int) -> None:
    """Open one session and run every plane in it. Raises; the caller records failures."""
    nwbfile, io = vn.open_session(row["path"])
    try:
        stim = vn.load_stimulus_table(nwbfile)
        key = (int(row["column"]), str(row["volume"]))
        ctx = {
            "session_name": row["name"],
            "spont": vn.spontaneous_block(nwbfile),
            "running": vn.load_running_speed(nwbfile),
            "dg_trials": {t: vn.stimulus_trials(stim, f"drifting_gratings_{t}",
                                                vn.DG_PARAM_COLUMNS)
                          for t in ("full", "windowed")},
            "ni_trials": {f: vn.stimulus_trials(stim, f)[0]
                          for f in ("natural_images", "natural_images_12")},
            "nm_trials": vn.stimulus_trials(stim, "natural_movie")[0],
            "lsn_trials": vn.stimulus_trials(stim, "locally_sparse_noise")[0],
            "lsn": vn.load_lsn_template(nwbfile),
            "center": centers.centers[key],
            "center_inferred": centers.inferred[key],
        }
        if acc.lsn_grid is None:
            acc.lsn_grid = {
                "altitudes": np.asarray(ctx["lsn"]["altitudes"], dtype=np.float64),
                "azimuths": np.asarray(ctx["lsn"]["azimuths"], dtype=np.float64)}

        for plane_key in vn.list_planes(nwbfile):
            plane = vn.load_plane(nwbfile, plane_key, trace_types=("events", "dff"))
            process_plane(plane, ctx, acc, config, seed)
            del plane
    finally:
        io.close()


def check_roi_coverage(tables: dict[str, pd.DataFrame]) -> int:
    """Every family must cover exactly the same ROIs; returns how many."""
    ref = set(map(tuple, tables["natural_movie"][KEYS].astype({"volume": str})
                  .to_numpy().tolist()))
    for fam, table in tables.items():
        got = set(map(tuple, table[KEYS].astype({"volume": str}).to_numpy().tolist()))
        if got != ref:
            raise AssertionError(f"{fam}: ROI set differs -- {len(got - ref)} extra, "
                                 f"{len(ref - got)} missing")
        if len(table) != len(ref):
            raise AssertionError(f"{fam}: {len(table)} rows for {len(ref)} distinct ROIs")
    return len(ref)


def build_wide(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict]:
    """One per-ROI table: identity once, then every family's metrics, prefixed."""
    def keyed(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["volume"] = out["volume"].astype(str)
        for k in ("column", "plane", "roi"):
            out[k] = out[k].astype("int64")
        return out

    wide = keyed(tables["natural_movie"][ID_COLS])
    manifest = {}
    for fam in FAMILIES:
        pub = keyed(to_output_schema(tables[fam], fam))
        metric_cols = [c for c in pub.columns if c not in ID_COLS]
        prefix = PREFIX[fam]
        part = pub[KEYS + metric_cols].rename(
            columns={c: prefix + c for c in metric_cols})
        # one_to_one fails loudly rather than fanning the table out on a duplicate key.
        wide = wide.merge(part, on=KEYS, how="left", validate="one_to_one")
        manifest[fam] = {"prefix": prefix, "columns": [prefix + c for c in metric_cols]}

    if not wide.columns.is_unique:
        dupes = wide.columns[wide.columns.duplicated()].tolist()
        raise AssertionError(f"duplicate columns: {dupes}")
    return wide, manifest


def _guarded(label: str, outputs: list, errors: list) -> Callable:
    """Run an array writer without letting a failure forfeit the run's provenance.

    Provenance and the manifest are written after the arrays, so an unguarded raise here
    costs a multi-hour run its record of what it did -- which has happened. Record the
    failure, say so loudly, and carry on.
    """
    def wrap(fn: Callable[[], Optional[str]]) -> None:
        try:
            name = fn()
            if name:
                outputs.append(name)
        except Exception as exc:                                        # noqa: BLE001
            errors.append({"output": label, "error": f"{type(exc).__name__}: {exc}",
                           "traceback": traceback.format_exc(limit=3)})
            print(f"!! {label} was not written: {exc}", flush=True)
    return wrap


def write_rf_maps(acc: Accumulator, tables: dict, save_dir: Path, seed: int
                  ) -> Optional[str]:
    """Pre-threshold ON/OFF subfield maps, (n_rois, 2, n_rows, n_cols).

    Altitudes, azimuths and the seed travel with the maps or pixel indices cannot be
    turned into degrees and per-pixel significance cannot be reproduced.
    """
    if not acc.rf_maps or acc.lsn_grid is None:
        return None
    path = save_dir / "receptive_field_maps.npz"
    np.savez_compressed(
        path, rf_maps=np.concatenate(acc.rf_maps, axis=0).astype(np.float32),
        roi_key=tables["rf_metrics"]["roi_key"].to_numpy(),
        altitudes=acc.lsn_grid["altitudes"], azimuths=acc.lsn_grid["azimuths"],
        seed=np.array(seed))
    return path.name


def write_tuning_curves(acc: Accumulator, tables: dict, save_dir: Path,
                        config: MetricConfig) -> Optional[str]:
    """Per-trial grating responses, blank sweeps, von Mises fits and running speeds.

    Trial-level rather than trial means: means are one nanmean away and the reverse is
    not. Blank sweeps are ragged across sessions (5-8), so they are NaN-padded to the
    widest plane and ``dg{w,f}_n_blank`` records each plane's true width.
    """
    if not all(acc.tuning[k]["trials"] for k in ("dgw", "dgf")) or acc.tuning_axes is None:
        return None

    payload, shapes = {}, {}
    for key in ("dgw", "dgf"):
        part = acc.tuning[key]
        for name, axes in (("trials", (1, 2, 3)), ("params", (1, 2))):
            seen = {tuple(a.shape[i] for i in axes) for a in part[name]}
            if len(seen) != 1:
                raise AssertionError(f"{key}_{name}: planes disagree on shape "
                                     f"{sorted(seen)}; cannot stack along the ROI axis")
        if len({a.shape for a in part["running"]}) != 1:
            raise AssertionError(f"{key}_running: planes disagree on shape")

        n_blank = np.array([a.shape[1] for a in part["blank"]], dtype=np.int32)
        pad_to = int(n_blank.max())
        blanks = [a if a.shape[1] == pad_to
                  else np.pad(a, ((0, 0), (0, pad_to - a.shape[1])),
                              constant_values=np.nan)
                  for a in part["blank"]]

        payload[f"{key}_trials"] = np.concatenate(part["trials"], axis=0)
        payload[f"{key}_blank"] = np.concatenate(blanks, axis=0)
        payload[f"{key}_n_blank"] = n_blank
        payload[f"{key}_params"] = np.concatenate(part["params"], axis=0)
        payload[f"{key}_running"] = np.stack(part["running"], axis=0)
        shapes[key] = payload[f"{key}_trials"].shape

    roi_key = tables["drifting_gratings_windowed"]["roi_key"].to_numpy()
    plane_key = np.asarray(acc.tuning["dgw"]["plane_key"])
    if payload["dgw_trials"].shape[0] != len(roi_key):
        raise AssertionError(f"tuning curves cover {payload['dgw_trials'].shape[0]} ROIs "
                             f"but the windowed table has {len(roi_key)}")
    if shapes["dgw"] != shapes["dgf"]:
        raise AssertionError(f"windowed {shapes['dgw']} and full field {shapes['dgf']} "
                             f"differ in shape")
    if payload["dgw_running"].shape[0] != len(plane_key):
        raise AssertionError("running speeds and plane keys are out of step")
    if not np.isin(np.array([k.rsplit("_", 1)[0] for k in roi_key]), plane_key).all():
        raise AssertionError("roi_key prefixes do not all appear in plane_key")

    path = save_dir / "tuning_curves.npz"
    np.savez_compressed(
        path, roi_key=roi_key, plane_key=plane_key,
        directions=np.asarray(acc.tuning_axes[0], dtype=np.float64),
        spatial_frequencies=np.asarray(acc.tuning_axes[1], dtype=np.float64),
        trace_type=np.array(config.trace_type["drifting_gratings_windowed"]),
        **payload)
    return path.name


def write_condition_means(acc: Accumulator, tables: dict, save_dir: Path,
                          config: MetricConfig) -> Optional[str]:
    """Neuron-by-image trial-mean responses for the two natural-image sets.

    Every published natural-image column is a reduction of this matrix, and a reduction
    cannot be un-taken. Natural movie is excluded: (n_rois, 3600) is ~541 MiB even after
    averaging over repeats.
    """
    if not any(acc.cond_means[f] for f in acc.cond_means):
        return None
    payload = {"roi_key": tables["natural_images"]["roi_key"].to_numpy(),
               "trace_type": np.array(config.trace_type["natural_images"])}
    for fam, short in (("natural_images", "ni"), ("natural_images_12", "ni12")):
        blocks = acc.cond_means[fam]
        if not blocks:
            continue
        widths = {b.shape[1] for b in blocks}
        if len(widths) != 1:
            raise AssertionError(f"{fam}: planes disagree on image count {sorted(widths)}")
        arr = np.concatenate(blocks, axis=0)
        if arr.shape[0] != len(payload["roi_key"]):
            raise AssertionError(f"{fam}: {arr.shape[0]} ROIs in the matrix vs "
                                 f"{len(payload['roi_key'])} in the table")
        if arr.shape[1] != len(acc.cond_ids[fam]):
            raise AssertionError(f"{fam}: {arr.shape[1]} columns vs "
                                 f"{len(acc.cond_ids[fam])} image ids")
        payload[f"{short}_mean"] = arr
        payload[f"{short}_images"] = np.asarray(acc.cond_ids[fam])

    path = save_dir / "condition_means.npz"
    np.savez_compressed(path, **payload)
    return path.name


def config_dict(cfg: MetricConfig) -> dict:
    """``MetricConfig`` as a plain dict; ``asdict`` cannot copy the MappingProxyType."""
    d = {f.name: getattr(cfg, f.name) for f in dataclasses.fields(MetricConfig)}
    d["trace_type"] = dict(d["trace_type"])
    return d


def _package_version(name: str) -> Optional[str]:
    try:
        return __import__(name).__version__
    except Exception:                                                   # noqa: BLE001
        return None


def build_provenance(*, asset_name: str, stamp: str, mouse_label: str, config: MetricConfig, seed: int,
                     input_asset: Path, sessions: pd.DataFrame, inventory: pd.DataFrame,
                     planes: pd.DataFrame, wide: pd.DataFrame, wall_seconds: float,
                     session_filter, failures: list, write_errors: list,
                     window_centers: dict, center_read_failures: list,
                     wide_name: str, manifest: dict, arrays: list) -> dict:
    """The asset's own provenance record.

    Built separately from writing it so its contents can be checked without a run --
    ``metadata.py`` reads several of these keys, and a rename there once cost an asset
    its record of whether validation had happened.
    """
    defaults, reference = config_dict(config), config_dict(REFERENCE_CONFIG)
    return {
        "asset": asset_name, "run_stamp": stamp,
        "generated_utc": pd.Timestamp.utcnow().isoformat(timespec="seconds"),
        "git_sha": prov.git_sha(), "mouse": mouse_label, "seed": seed,
        "input_asset": str(input_asset),
        "n_sessions": int(len(sessions)), "n_planes": int(len(planes)),
        "n_rois": int(len(wide)), "wall_seconds": round(wall_seconds, 1),
        "session_filter": list(session_filter) if session_filter else None,
        "complete_asset": bool(session_filter is None and not failures
                               and not write_errors),
        "sessions": [{"name": r["name"], "format": r["format"], "column": int(r["column"]),
                      "volume": str(r["volume"]), "n_planes": int(r["n_planes"])}
                     for _, r in sessions.iterrows()],
        "unreadable_sessions": int(inventory["error"].notna().sum()),
        "failed_sessions": failures,
        "failed_outputs": write_errors,
        "config": defaults,
        "window_centers": {**window_centers,
                           "read_failures": center_read_failures},
        # The delta rather than prose, so the claim stays checkable.
        "differs_from_reference_config": {
            k: {"used": defaults[k], "historical": reference[k]}
            for k in defaults if defaults[k] != reference[k]},
        "outputs": {"table": wide_name, "columns": manifest,
                    "arrays": sorted(arrays)},
        "environment": {
            "python": platform.python_version(), "platform": platform.platform(),
            "packages": {m: _package_version(m) for m in
                         ("numpy", "pandas", "scipy", "pynwb", "hdmf", "hdmf_zarr",
                          "h5py", "zarr", "pyarrow")}},
    }


def run(input_asset: Path, results_dir: Path, asset_prefix: str = "V1DD_functional_metrics",
        *, config: MetricConfig = DEFAULT_CONFIG, seed: int = 0,
        session_filter: Optional[Sequence[tuple]] = None) -> Path:
    """Build one asset directory under ``results_dir`` and return its path.

    ``input_asset`` is the mounted NWB dataset. ``session_filter`` restricts the run to a
    list of ``(column, volume)`` pairs, which is recorded so a partial asset cannot pass
    for a complete one.
    """
    t_start = time.time()
    sessions, inventory = discover_sessions(str(input_asset), session_filter)
    if not len(sessions):
        raise FileNotFoundError(f"no readable sessions under {input_asset}")

    with vn.session(sessions["path"].iloc[0]) as first:
        mouse_id, mouse_label = vn.session_mouse(first, sessions["path"].iloc[0])

    asset_name = f"{mouse_id}_{asset_prefix}"
    stamp = prov.run_stamp()
    save_dir = Path(prov.run_dir(str(results_dir), asset_name, stamp))
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(sessions)} session(s), {int(sessions['n_planes'].sum())} planes, "
          f"mouse {mouse_label}\nsave_dir: {save_dir}", flush=True)

    centers, center_read_failures = aperture_centers(sessions, config)
    wc = centers.provenance
    print(f"aperture centres: {wc['n_measured']} measured, {wc['n_inferred']} inferred, "
          f"{wc['n_unfilled']} unknown", flush=True)

    acc = Accumulator()
    for n, (_, row) in enumerate(sessions.iterrows(), 1):
        t_sess = time.time()
        try:
            process_session(row, centers, acc, config, seed)
            print(f"  [{n:>2}/{len(sessions)}] col{row['column']} vol{row['volume']} "
                  f"{row['name'][:34]:<34} {time.time() - t_sess:>6.1f}s", flush=True)
        except Exception as exc:                                        # noqa: BLE001
            acc.failures.append({"name": row["name"],
                                 "error": f"{type(exc).__name__}: {exc}",
                                 "traceback": traceback.format_exc(limit=3)})
            print(f"  [{n:>2}/{len(sessions)}] {row['name'][:34]:<34} FAILED: {exc}",
                  flush=True)

    wall_seconds = time.time() - t_start
    tables = acc.tables()
    if not tables:
        raise RuntimeError("no session produced any metrics")
    planes = pd.DataFrame(acc.plane_log)
    n_rois = check_roi_coverage(tables)
    print(f"\n{len(planes)} planes, {n_rois} ROIs, {wall_seconds / 60:.1f} min",
          flush=True)

    wide, manifest = build_wide(tables)
    wide_path = save_dir / "stimulus_metrics.parquet"
    wide.to_parquet(wide_path, index=False)
    print(f"wrote {wide_path.name}: {len(wide)} ROIs x {len(wide.columns)} columns",
          flush=True)

    # Array writers are guarded individually: they run before provenance, and an
    # unguarded raise in one of them costs the run its record of everything else.
    arrays: list[str] = []
    write_errors: list[dict] = []
    _guarded("receptive_field_maps.npz", arrays, write_errors)(
        lambda: write_rf_maps(acc, tables, save_dir, seed))
    _guarded("tuning_curves.npz", arrays, write_errors)(
        lambda: write_tuning_curves(acc, tables, save_dir, config))
    _guarded("condition_means.npz", arrays, write_errors)(
        lambda: write_condition_means(acc, tables, save_dir, config))

    record = build_provenance(
        asset_name=asset_name, stamp=stamp, mouse_label=mouse_label,
        config=config, seed=seed,
        input_asset=input_asset, sessions=sessions, inventory=inventory,
        planes=planes, wide=wide, wall_seconds=wall_seconds,
        session_filter=session_filter, failures=acc.failures,
        write_errors=write_errors, window_centers=wc,
        center_read_failures=center_read_failures,
        wide_name=wide_path.name, manifest=manifest, arrays=arrays)
    with open(save_dir / "provenance.json", "w", encoding="utf-8") as fh:
        json.dump(prov.jsonable(record), fh, indent=2, sort_keys=True, allow_nan=False)
        fh.write("\n")

    print(f"\n{save_dir}", flush=True)
    for f in sorted(save_dir.iterdir()):
        print(f"  {f.name:<40} {f.stat().st_size / 1024:>10.1f} KB", flush=True)
    if acc.failures:
        print(f"!! {len(acc.failures)} session(s) failed -- the asset is incomplete")
    if write_errors:
        print(f"!! {len(write_errors)} array output(s) missing: "
              f"{[e['output'] for e in write_errors]}")
    return save_dir
