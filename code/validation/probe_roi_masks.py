"""Read-only probe: how ROI masks are stored, and whether positions can be derived.

Answers the questions that decide the design of an ROI-position step, none of which can
be settled off-capsule because they are about the NWB files themselves:

1. Which mask column each session carries -- ``pixel_mask``, ``image_mask``, or both --
   and whether that tracks the storage format.
2. What is actually in them: field names, dtypes, and **whether the weights are all 1**.
   A weighted centroid over a graded mask and an unweighted one over a binary mask are
   different definitions, so mixing them across sessions would be silent.
3. Which axis is which. ``pixel_mask`` is a list of coordinates whose field order is not
   guaranteed to be (x, y); ``image_mask`` is (row, col). Getting this backwards
   transposes the field of view and nothing downstream would notice.
4. Whether a centroid is already stored by the segmentation, making all of this moot.
5. Whether the imaging plane carries ``grid_spacing`` / ``origin_coords``, which is what
   turns pixels into micrometres and, if origins are in a common frame, what would let
   ROIs from different columns share one anatomical space.
6. What it costs to read -- an ``image_mask`` is dense, so a plane of ~300 ROIs at
   512x512 is far larger than the ragged form.

    python code/validation/probe_roi_masks.py [output-dir] [--sessions N] [--rois N]

Writes ``roi_mask_probe.json`` and prints a summary. Never raises on one session: a file
that cannot be read is recorded and the walk continues.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v1dd_metrics import nwb as vn                              # noqa: E402
from v1dd_metrics import provenance as prov                     # noqa: E402

MASK_COLUMNS = ("pixel_mask", "image_mask", "voxel_mask")


def _jsonable(obj):
    return prov.jsonable(obj)


def describe_pixel_mask(ps, n_rois: int, sample: int) -> dict:
    """Field names, weight values and per-ROI pixel counts from the ragged form."""
    info: dict = {"kind": "pixel_mask"}
    col = ps["pixel_mask"]
    rows = []
    for i in range(min(sample, n_rois)):
        entry = np.asarray(col[i])
        rows.append(entry)
    if not rows:
        return info

    first = rows[0]
    info["entry_dtype"] = str(first.dtype)
    info["fields"] = list(first.dtype.names) if first.dtype.names else None
    info["entry_shape"] = list(np.shape(first))
    counts = [len(r) for r in rows]
    info["pixels_per_roi"] = {"min": int(min(counts)), "max": int(max(counts)),
                              "median": float(np.median(counts))}

    flat = np.concatenate([np.asarray(r).reshape(len(r), -1) for r in rows]) \
        if first.dtype.names is None else None
    if first.dtype.names:
        cols = {name: np.concatenate([np.asarray(r[name]) for r in rows])
                for name in first.dtype.names}
    else:
        # unstructured (n, 3): assume the NWB convention (x, y, weight) but REPORT the
        # ranges so the assumption is checkable rather than believed
        cols = {f"col{j}": flat[:, j] for j in range(flat.shape[1])}
    info["column_ranges"] = {k: [float(v.min()), float(v.max())] for k, v in cols.items()}
    weights = cols.get("weight", cols.get(f"col{len(cols) - 1}"))
    if weights is not None:
        uniq = np.unique(weights)
        info["weights_all_one"] = bool(np.allclose(weights, 1.0))
        info["distinct_weights"] = int(uniq.size)
        info["weight_sample"] = [float(x) for x in uniq[:8]]
    return info


def describe_image_mask(ps, n_rois: int, sample: int) -> dict:
    """dtype, shape and value set of the dense form."""
    info: dict = {"kind": "image_mask"}
    col = ps["image_mask"]
    data = getattr(col, "data", col)
    info["stored_shape"] = list(getattr(data, "shape", []))
    info["stored_dtype"] = str(getattr(data, "dtype", ""))

    counts, uniq_all = [], set()
    for i in range(min(sample, n_rois)):
        m = np.asarray(data[i])
        counts.append(int((m > 0).sum()))
        uniq_all.update(np.unique(m).tolist()[:8])
    if counts:
        info["nonzero_per_roi"] = {"min": min(counts), "max": max(counts),
                                   "median": float(np.median(counts))}
    info["value_sample"] = sorted(uniq_all)[:8]
    info["is_binary"] = bool(uniq_all <= {0, 1, 0.0, 1.0, True, False})
    return info


def centroid_from_pixel_mask(entry) -> tuple:
    """Unweighted centroid of a ragged entry, reported per field so axes stay explicit."""
    e = np.asarray(entry)
    if e.dtype.names:
        names = list(e.dtype.names)
        return tuple(float(np.mean(e[n])) for n in names[:2])
    arr = e.reshape(len(e), -1)
    return float(arr[:, 0].mean()), float(arr[:, 1].mean())


def probe_plane(nwbfile, plane_key: str, sample: int) -> dict:
    """Everything about one plane's masks and imaging geometry."""
    module = nwbfile.processing[plane_key]
    seg = module["image_segmentation"]
    ps = seg[list(seg.plane_segmentations)[0]] if hasattr(seg, "plane_segmentations") \
        else list(seg.children)[0]

    colnames = list(getattr(ps, "colnames", []))
    out: dict = {"plane": plane_key, "n_rois": len(ps.id),
                 "roi_table_columns": colnames,
                 "mask_columns": [c for c in MASK_COLUMNS if c in colnames],
                 # a centroid already stored by the segmentation would make this moot
                 "centroidish_columns": [c for c in colnames if any(
                     k in c.lower() for k in ("centroid", "center", "_x", "_y", "pos"))]}

    n_rois = len(ps.id)
    t0 = time.time()
    if "pixel_mask" in colnames:
        out["pixel_mask"] = describe_pixel_mask(ps, n_rois, sample)
    if "image_mask" in colnames:
        out["image_mask"] = describe_image_mask(ps, n_rois, sample)
    out["read_seconds_for_sample"] = round(time.time() - t0, 3)

    # Both present is the one case that can settle the axis convention outright.
    if "pixel_mask" in colnames and "image_mask" in colnames:
        agree = []
        for i in range(min(sample, n_rois)):
            a = centroid_from_pixel_mask(ps["pixel_mask"][i])
            m = np.asarray(getattr(ps["image_mask"], "data", ps["image_mask"])[i])
            yx = np.argwhere(m > 0)
            if not len(yx):
                continue
            b = (float(yx[:, 1].mean()), float(yx[:, 0].mean()))   # (col=x, row=y)
            agree.append([a, b])
        out["both_formats_centroids"] = agree[:5]

    plane_obj = getattr(ps, "imaging_plane", None)
    geom = {}
    for attr in ("grid_spacing", "grid_spacing_unit", "origin_coords",
                 "origin_coords_unit", "imaging_rate", "location", "description",
                 "unit"):
        val = getattr(plane_obj, attr, None)
        if val is not None:
            geom[attr] = np.asarray(val).tolist() if hasattr(val, "__len__") \
                and not isinstance(val, str) else val
    out["imaging_plane"] = geom
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", nargs="?", default=".", help="where to write the JSON")
    ap.add_argument("--input-asset", default="/data/409828_V1DD_Filtered")
    ap.add_argument("--sessions", type=int, default=0,
                    help="probe only the first N sessions (0 = all)")
    ap.add_argument("--planes", type=int, default=1,
                    help="planes to probe per session")
    ap.add_argument("--rois", type=int, default=12,
                    help="ROIs to sample per plane")
    args = ap.parse_args(argv)

    paths = vn.find_sessions(args.input_asset)
    if args.sessions:
        paths = paths[:args.sessions]
    print(f"{len(paths)} session(s) under {args.input_asset}\n")

    report = {"input_asset": args.input_asset, "sessions": []}
    for path in paths:
        entry = {"path": str(path), "name": Path(path).name,
                 "format": vn.nwb_format(path)}
        try:
            with vn.session(path) as nwbfile:
                keys = vn.list_planes(nwbfile)
                entry["n_planes"] = len(keys)
                entry["planes"] = [probe_plane(nwbfile, k, args.rois)
                                   for k in keys[:args.planes]]
        except Exception as exc:                                    # noqa: BLE001
            entry["error"] = f"{type(exc).__name__}: {exc}"
            entry["traceback"] = traceback.format_exc(limit=3)
        report["sessions"].append(entry)

        masks = sorted({m for p in entry.get("planes", []) for m in p["mask_columns"]})
        print(f"  {entry['name'][:44]:<44} {entry['format']:<5} "
              f"{'/'.join(masks) or entry.get('error', 'no mask column')}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "roi_mask_probe.json").write_text(
        json.dumps(_jsonable(report), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # The headline questions, answered across everything that was readable
    planes = [p for s in report["sessions"] for p in s.get("planes", [])]
    by_format: dict = {}
    for s in report["sessions"]:
        for p in s.get("planes", []):
            by_format.setdefault(s["format"], set()).update(p["mask_columns"])
    print("\nmask column by storage format:")
    for fmt, cols in sorted(by_format.items()):
        print(f"  {fmt:<6} {sorted(cols)}")

    weights = [p["pixel_mask"].get("weights_all_one") for p in planes if "pixel_mask" in p]
    if weights:
        print(f"\npixel_mask weights all 1: {weights.count(True)}/{len(weights)} planes")
        fields = {tuple(p['pixel_mask'].get('fields') or []) for p in planes
                  if 'pixel_mask' in p}
        print(f"pixel_mask field names seen: {sorted(fields)}")
    binaries = [p["image_mask"].get("is_binary") for p in planes if "image_mask" in p]
    if binaries:
        print(f"image_mask binary: {binaries.count(True)}/{len(binaries)} planes")

    centroidish = sorted({c for p in planes for c in p["centroidish_columns"]})
    print(f"\ncolumns that might already hold a position: {centroidish or 'none'}")
    geoms = [p["imaging_plane"] for p in planes]
    have_spacing = sum(1 for g in geoms if g.get("grid_spacing") is not None)
    have_origin = sum(1 for g in geoms if g.get("origin_coords") is not None)
    print(f"grid_spacing present : {have_spacing}/{len(geoms)} planes"
          f"   -> pixels to micrometres")
    print(f"origin_coords present: {have_origin}/{len(geoms)} planes"
          f"   -> a shared frame across columns")
    if geoms:
        print(f"first plane geometry : {geoms[0]}")
    print(f"\nwrote {out_dir / 'roi_mask_probe.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
