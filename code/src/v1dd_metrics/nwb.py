"""Reading V1DD two-photon sessions out of NWB.

The only module that touches ``hdmf_zarr`` or ``pynwb.NWBHDF5IO``, so a storage-layout
change moves one file. The asset mixes NWB-Zarr directories and plain HDF5, traces keep
NWB's (n_frames, n_rois) orientation, and every stimulus shares one concatenated table.
See docs/data_access.md.
"""

import os
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

__all__ = [
    "PlaneData",
    "open_session",
    "session",
    "nwb_format",
    "find_sessions",
    "peek_session",
    "list_planes",
    "load_plane",
    "session_mouse",
    "plane_depth_um",
    "load_stimulus_table",
    "stimulus_trials",
    "epoch_table",
    "spontaneous_block",
    "load_running_speed",
    "load_lsn_template",
    "RoiMasks",
    "load_roi_masks",
    "STIMULUS_TABLE_COLUMNS",
    "DG_PARAM_COLUMNS",
]

#: Columns the NWB writer is expected to emit. Verified against a real file in M1;
#: `schema_report` reports what is actually present rather than assuming this.
STIMULUS_TABLE_COLUMNS = [
    "stim_name", "start_time", "stop_time", "temporal_frequency", "spatial_frequency",
    "center_azimuth", "center_elevation", "direction", "frame", "image_order",
    "image_index", "stimulus_condition_id",
]

#: The columns that define a drifting-gratings condition. Blank (gray) sweeps are rows
#: that are NaN in these — NOT rows that are NaN anywhere, which in the concatenated
#: table is every drifting-gratings row.
DG_PARAM_COLUMNS = ["temporal_frequency", "spatial_frequency", "direction"]

#: `ImagingPlane.location` when it carries a depth: a number and an optional unit, and
#: nothing else. Anchored on purpose -- see `plane_depth_um`.
_DEPTH_RE = re.compile(r"^(-?\d+(?:\.\d+)?)\s*(?:um|µm|μm|micron|micrometers?)?$",
                       re.IGNORECASE)


@dataclass
class PlaneData:
    """One imaging plane's traces and ROI identity."""

    column: int
    volume: str                       # str: volumes run 1..9 and a..f across the project
    plane: int
    roi: np.ndarray                   # (n_rois,) ROI ids from the file, not 0..n-1
    is_valid: np.ndarray              # (n_rois,) bool, pika_roi_confidence > 0.5
    timestamps: np.ndarray            # (n_frames,) seconds
    traces: Dict[str, np.ndarray]     # trace_type -> (n_frames, n_rois)
    roi_table: pd.DataFrame
    dt: float                         # median inter-frame interval, seconds
    mouse_id: str = ""                # bare number, e.g. "409828"; see session_mouse
    depth_um: Optional[float] = None  # imaging depth, None if the file does not say

    @property
    def n_rois(self) -> int:
        return len(self.roi)

    @property
    def fs(self) -> float:
        return 1.0 / self.dt


def nwb_format(path) -> str:
    """`"zarr"` for a `.nwb.zarr` directory, `"hdf5"` for a plain `.nwb` file."""
    return "zarr" if str(path).endswith(".nwb.zarr") else "hdf5"


def open_session(path):
    """Open a session read-only, in whichever format it was written.

    Returns ``(nwbfile, io)``; the caller owns closing ``io``, which is not optional for
    HDF5, where the handle lazily backs every array. Prefer ``session()`` in a loop.
    """
    path = str(path)
    if nwb_format(path) == "zarr":
        from hdmf_zarr import NWBZarrIO

        io = NWBZarrIO(path, mode="r")
    else:
        from pynwb import NWBHDF5IO

        io = NWBHDF5IO(path, mode="r")
    return io.read(), io


@contextmanager
def session(path):
    """Context manager yielding an open ``NWBFile``, closing the handle on exit.

    Prefer this over ``open_session`` in loops: an unclosed ``NWBHDF5IO`` keeps a file
    handle per session and a long loop reaches the open-file limit.
    """
    nwb, io = open_session(path)
    try:
        yield nwb
    finally:
        io.close()


def _walk_sessions(root: Path) -> Tuple[List[Path], List[Path]]:
    """(zarr directories, hdf5 files) under ``root``, from one walk that prunes Zarr stores.

    ``os.walk`` rather than ``rglob``, which descends into a Zarr store's thousands of chunk
    directories and effectively never finishes. Nothing inside a store is ever a session.
    """
    zarr_paths: List[Path] = []
    hdf5_paths: List[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        here = Path(dirpath)
        # Collect the stores before pruning them -- a `.nwb.zarr` *is* a session, it is
        # only its interior that is uninteresting.
        zarr_paths.extend(here / d for d in dirnames if d.endswith(".nwb.zarr"))
        dirnames[:] = [d for d in dirnames if not d.endswith(".nwb.zarr")]
        hdf5_paths.extend(here / f for f in filenames if f.endswith(".nwb"))
    return sorted(zarr_paths), sorted(hdf5_paths)


def find_sessions(functional_dir, prefer: str = "zarr") -> List[Path]:
    """One NWB path per session directory, across both storage formats.

    A session present in both formats yields one path, ``prefer`` deciding which. The
    crawling fallback fires only when the expected layout yields nothing at all -- see
    docs/data_access.md for why a per-format fallback hangs on a single-format asset.
    """
    root = Path(functional_dir)
    if not root.is_dir():
        raise FileNotFoundError(
            f"{functional_dir} does not exist. On CodeOcean, confirm the functional "
            "data asset is attached; locally, set V1DD_DATA_ROOT."
        )

    zarr_paths = sorted(root.glob("*/*.nwb.zarr"))
    hdf5_paths = sorted(root.glob("*/*.nwb"))
    if not zarr_paths and not hdf5_paths:
        zarr_paths, hdf5_paths = _walk_sessions(root)

    if prefer not in ("zarr", "hdf5"):
        raise ValueError("prefer must be 'zarr' or 'hdf5'")
    first, second = (zarr_paths, hdf5_paths) if prefer == "zarr" else (hdf5_paths, zarr_paths)

    by_session: Dict[Path, Path] = {}
    for p in first:
        by_session.setdefault(p.parent, p)
    for p in second:
        by_session.setdefault(p.parent, p)

    if not by_session:
        raise FileNotFoundError(f"no *.nwb.zarr or *.nwb found under {functional_dir}")
    return [by_session[k] for k in sorted(by_session)]


def peek_session(path) -> Dict[str, Any]:
    """Identify a session without loading any traces.

    Reads `(column, volume)` from the ROI table *inside* the file rather than from the
    directory name, and records which storage format it came from. A file that cannot be
    opened yields a row with NaN identifiers and an `error` string, so one bad session
    shortens nothing and is visible in the index.
    """
    path = Path(path)
    info: Dict[str, Any] = {
        "path": str(path), "name": path.parent.name, "format": nwb_format(path),
        "column": np.nan, "volume": np.nan, "n_planes": np.nan, "session_id": "",
        "error": None,
    }
    io = None
    try:
        nwb, io = open_session(path)
        planes = list_planes(nwb)
        info["n_planes"] = len(planes)
        info["session_id"] = str(nwb.session_id)
        if planes:
            rois = nwb.processing[planes[0]]["dff"].rois.to_dataframe()
            info["column"] = int(pd.to_numeric(rois["column"]).iloc[0])
            info["volume"] = _as_volume_str(rois["volume"].iloc[0])
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if io is not None:
            try:
                io.close()
            except Exception:
                pass
    return info


def list_planes(nwb) -> List[str]:
    """Imaging-plane processing module names, ordered by plane number."""
    keys = [k for k in nwb.processing if k.startswith("plane-")]
    return sorted(keys, key=lambda k: int(k.split("-", 1)[1]))


def _as_volume_str(value) -> str:
    """Volumes are 1..9 and a..f; normalise to str once, here, and never re-parse."""
    try:
        return str(int(value))
    except (TypeError, ValueError):
        return str(value)


def session_mouse(nwb, path=None) -> Tuple[str, str]:
    """``(mouse_id, mouse_label)`` for a session, e.g. ``("409828", "M409828")``.

    The bare number goes inside ``roi_unique_id``; the M-prefixed label is what the ``mouse``
    column and filenames carry. Resolved from ``nwb.subject.subject_id``, else the leading
    token of the session directory name when ``path`` is given, else raises -- never
    silently defaulted, since a wrong id corrupts every ROI identifier in the output.
    """
    subject = getattr(nwb, "subject", None)
    raw = getattr(subject, "subject_id", None) if subject is not None else None
    if raw is None or not str(raw).strip():
        token = Path(path).parent.name.split("_", 1)[0] if path else ""
        if not token.isdigit():
            raise ValueError(
                "cannot determine the mouse: nwb.subject.subject_id is empty and the "
                f"session name {token!r} does not start with a numeric id"
            )
        raw = token
    mouse_id = str(raw).strip().lstrip("Mm")
    return mouse_id, f"M{mouse_id}"


def plane_depth_um(nwb, plane) -> Optional[float]:
    """Imaging depth of one plane in micrometres, or None if the file does not say.

    Read from ``ImagingPlane.location``, which this asset writes as free text ("50 um").
    The pattern is anchored -- a number and an optional unit and nothing else -- because a
    loose match on free text returns a plausible wrong depth, which is worse than none.
    """
    key = f"plane-{int(plane)}" if not isinstance(plane, str) else plane
    try:
        series = nwb.processing[key][
            "dff" if "dff" in nwb.processing[key].data_interfaces else "events"]
        location = getattr(series.rois.table.imaging_plane, "location", None)
    except Exception:                                        # noqa: BLE001
        return None
    if location is None:
        return None
    match = _DEPTH_RE.match(str(location).strip())
    return float(match.group(1)) if match else None


def load_plane(
    nwb,
    plane,
    trace_types: Sequence[str] = ("dff", "events"),
    valid_only: bool = False,
) -> PlaneData:
    """Load one plane's traces, timestamps and ROI table.

    ``trace_types`` are the data interfaces to pull; loading several at once is cheaper than
    reopening the plane. The shipped metrics use ``events`` for everything except receptive
    fields, which use ``dff``.
    """
    key = f"plane-{int(plane)}" if not isinstance(plane, str) else plane
    module = nwb.processing[key]

    missing = [t for t in trace_types if t not in module.data_interfaces]
    if missing:
        raise KeyError(
            f"{key} has no {missing}; available: {sorted(module.data_interfaces)}"
        )

    first = module[trace_types[0]]
    timestamps = np.asarray(first.timestamps[:], dtype=np.float64)

    roi_table = first.rois.to_dataframe().reset_index(drop=False)
    roi_table = roi_table.rename(columns={"id": "nwb_roi_row_id"})
    roi_table = roi_table.drop(columns=[c for c in ("pixel_mask",) if c in roi_table])

    traces: Dict[str, np.ndarray] = {}
    for t in trace_types:
        series = module[t]
        data = np.asarray(series.data[:])
        if data.shape[0] != timestamps.shape[0]:
            if data.ndim == 2 and data.shape[1] == timestamps.shape[0]:
                data = data.T
            else:
                raise ValueError(
                    f"{key}/{t}: data {data.shape} incompatible with timestamps "
                    f"{timestamps.shape}"
                )
        if data.shape[1] != len(roi_table):
            raise ValueError(
                f"{key}/{t}: {data.shape[1]} trace columns vs {len(roi_table)} roi rows"
            )
        traces[t] = data

    if "pika_roi_confidence" in roi_table:
        is_valid = roi_table["pika_roi_confidence"].to_numpy() > 0.5
    elif "is_soma" in roi_table:
        is_valid = roi_table["is_soma"].to_numpy().astype(bool)
    else:
        is_valid = np.ones(len(roi_table), dtype=bool)

    if valid_only:
        traces = {t: d[:, is_valid] for t, d in traces.items()}
        roi_table = roi_table[is_valid].reset_index(drop=True)
        is_valid = np.ones(len(roi_table), dtype=bool)

    return PlaneData(
        column=int(pd.to_numeric(roi_table["column"]).iloc[0]),
        volume=_as_volume_str(roi_table["volume"].iloc[0]),
        plane=int(pd.to_numeric(roi_table["plane"]).iloc[0]),
        roi=pd.to_numeric(roi_table["roi"]).to_numpy(dtype=int),
        is_valid=is_valid,
        timestamps=timestamps,
        traces=traces,
        roi_table=roi_table,
        dt=float(np.median(np.diff(timestamps))),
        mouse_id=session_mouse(nwb)[0],
        depth_um=plane_depth_um(nwb, key),
    )


def load_stimulus_table(nwb) -> pd.DataFrame:
    """The flat per-trial stimulus table, chronological, index reset."""
    df = nwb.intervals["stimulus_table"].to_dataframe().reset_index(drop=False)
    return df.sort_values("start_time").reset_index(drop=True)


def epoch_table(nwb) -> pd.DataFrame:
    """Stimulus-block boundaries, sorted by start time."""
    df = nwb.intervals["epochs"].to_dataframe().reset_index(drop=False)
    if "duration" not in df:
        df["duration"] = df["stop_time"] - df["start_time"]
    return df.sort_values("start_time").reset_index(drop=True)


def stimulus_trials(
    stim_table: pd.DataFrame,
    stim_name: str,
    param_columns: Sequence[str] = (),
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Rows for one stimulus, chronological, index reset, plus a blank-sweep mask.

    ``param_columns`` are the columns defining a condition for this stimulus; a sweep is
    blank when it is NaN in *those* columns. Restricting the test is essential in the
    concatenated table, where every row is NaN in some other stimulus's columns. Pass ``()``
    for stimuli with no blank sweeps.
    """
    trials = (
        stim_table.loc[stim_table["stim_name"] == stim_name]
        .sort_values("start_time")
        .reset_index(drop=True)
    )
    if len(param_columns):
        present = [c for c in param_columns if c in trials.columns]
        is_blank = trials[present].isna().any(axis=1).to_numpy() if present else np.zeros(len(trials), bool)
    else:
        is_blank = np.zeros(len(trials), dtype=bool)
    return trials, is_blank


def spontaneous_block(nwb, which: int = 0) -> Tuple[float, float]:
    """Start/stop of a spontaneous epoch, in seconds.

    The original drew its null distribution from `spont_stim_table.iloc[0]` — the first
    spontaneous block only. These sessions have exactly one, so first and only coincide,
    but `which` is explicit so that a session with more than one cannot silently change
    the answer.
    """
    ep = epoch_table(nwb)
    spont = ep.loc[ep["stim_name"] == "spontaneous"].reset_index(drop=True)
    if not len(spont):
        raise LookupError("no spontaneous epoch in this session")
    if which >= len(spont):
        raise IndexError(f"requested spontaneous block {which} of {len(spont)}")
    row = spont.iloc[which]
    return float(row["start_time"]), float(row["stop_time"])


def load_running_speed(nwb) -> Tuple[np.ndarray, np.ndarray]:
    """`(speed_cm_s, timestamps)`, on the behaviour clock — not the imaging clock.

    Already differentiated by the NWB writer. The old client computed a central
    difference of cumulative distance itself, so small disagreements here propagate into
    the running/stationary surround-suppression variants, which threshold at exactly
    1 cm/s.
    """
    series = nwb.processing["behavior"]["running_speed"]
    return (
        np.asarray(series.data[:], dtype=np.float64),
        np.asarray(series.timestamps[:], dtype=np.float64),
    )


def _images_to_array(images) -> Tuple[np.ndarray, List[int]]:
    """Stack an NWB `Images` container into (n, rows, cols), ordered by integer name."""
    named = getattr(images, "images", None) or {}
    keys = sorted(named, key=lambda k: int(k))
    stack = np.stack([np.asarray(named[k].data[:]) for k in keys])
    return stack, [int(k) for k in keys]


def load_lsn_template(
    nwb,
    downsample_to: Optional[Tuple[int, int]] = (8, 14),
    grid_size_deg: float = 9.3,
) -> Dict[str, Any]:
    """The locally-sparse-noise template, reduced to the grid the RF code expects.

    NWB embeds a 16x28 grid where the original analysis used 8x14, so this block-reduces and
    asserts every 2x2 block is uniform; a failure means the two are not the same stimulus.
    Returns a dict (``images``, ``frames``, ``azimuths``, ``altitudes``, ``native_shape``,
    ``reduced``, ``blocks_uniform``) so a caller can inspect that rather than only catch it.
    See docs/data_access.md.
    """
    images, frames = _images_to_array(nwb.stimulus["locally_sparse_noise"])
    native_shape = tuple(images.shape[1:])
    out: Dict[str, Any] = {
        "frames": frames,
        "native_shape": native_shape,
        "reduced": False,
        "blocks_uniform": None,
        "grid_size_deg": grid_size_deg,
    }

    if downsample_to is not None and native_shape != tuple(downsample_to):
        want_r, want_c = downsample_to
        have_r, have_c = native_shape
        if have_r % want_r or have_c % want_c:
            out["images"] = images
            out["error"] = (
                f"template {native_shape} is not an integer multiple of {downsample_to}; "
                "cannot reduce to the grid the RF code was written against"
            )
            out["azimuths"], out["altitudes"] = _grid_degrees(have_c, have_r, grid_size_deg)
            out.update(_pixel_codes(images))
            return out

        fr, fc = have_r // want_r, have_c // want_c
        blocks = images.reshape(len(images), want_r, fr, want_c, fc)
        uniform = bool(np.all(blocks == blocks[:, :, :1, :, :1]))
        out["blocks_uniform"] = uniform
        out["block_factor"] = (fr, fc)
        if not uniform:
            out["images"] = images
            out["error"] = (
                f"{fr}x{fc} blocks are not uniform, so the {native_shape} template is not "
                f"a {fr}x{fc} upsample of a {downsample_to} grid — receptive fields "
                "cannot be reproduced from this asset"
            )
            out["azimuths"], out["altitudes"] = _grid_degrees(have_c, have_r, grid_size_deg)
            out.update(_pixel_codes(images))
            return out

        images = blocks[:, :, 0, :, 0]
        out["reduced"] = True

    rows, cols = images.shape[1:]
    out["images"] = images
    out["azimuths"], out["altitudes"] = _grid_degrees(cols, rows, grid_size_deg)
    out.update(_pixel_codes(images))
    return out


def _pixel_codes(images: np.ndarray) -> Dict[str, Any]:
    """Which pixel value means bright, dark, and background.

    `locally_sparse_noise.py` hard-codes `pixel_on=255`, `pixel_off=0`, `pixel_gray=127`,
    but this asset encodes the template as **-1 / 0 / 1**. Porting those constants
    literally makes the ON and OFF design matrices all-False, so every ROI comes out with
    no receptive field — a silent wrong answer that reads as a biological result rather
    than a bug. So derive the mapping from the data: brightest is ON, darkest is OFF,
    the middle value is background.
    """
    values = sorted(int(v) for v in np.unique(images))
    out: Dict[str, Any] = {"pixel_values": values}
    if len(values) == 3:
        out["pixel_off"], out["pixel_gray"], out["pixel_on"] = values
    elif len(values) == 2:                      # no background pixels in this template
        out["pixel_off"], out["pixel_on"] = values
        out["pixel_gray"] = None
    else:
        out["pixel_off"] = out["pixel_gray"] = out["pixel_on"] = None
        out["pixel_error"] = (
            f"expected 2 or 3 distinct pixel values, found {len(values)}: {values[:10]}"
        )
    return out


def _grid_degrees(n_cols: int, n_rows: int, grid: float) -> Tuple[np.ndarray, np.ndarray]:
    """Screen-centred pixel-centre coordinates, matching `LocallySparseNoise._load_frames`."""
    azimuths = (np.arange(n_cols) - n_cols // 2 + 0.5) * grid
    altitudes = (np.arange(n_rows) - n_rows // 2 + 0.5) * grid
    return azimuths, altitudes


@dataclass
class RoiMasks:
    """One plane's ROI footprints, flattened to pixel coordinates.

    ``row``/``col`` hold every in-mask pixel for every ROI concatenated, and ``offsets``
    slices them per ROI, so a reduction over all ROIs is one pass. ``shape`` is the field
    of view in pixels where the file says, else None; ``source`` names the column the
    coordinates came from.
    """

    row: np.ndarray                   # (n_pixels,) pixel row, i.e. y
    col: np.ndarray                   # (n_pixels,) pixel column, i.e. x
    offsets: np.ndarray               # (n_rois + 1,) slice bounds into row/col
    source: str                       # "pixel_mask" or "image_mask"
    shape: Optional[Tuple[int, int]] = None
    weights_all_one: bool = True

    @property
    def n_rois(self) -> int:
        return len(self.offsets) - 1


def _plane_segmentation(nwbfile, plane):
    """The PlaneSegmentation for one plane."""
    key = f"plane-{int(plane)}" if not isinstance(plane, str) else plane
    seg = nwbfile.processing[key]["image_segmentation"]
    names = list(getattr(seg, "plane_segmentations", {}) or {})
    return seg[names[0]] if names else list(seg.children)[0]


def load_roi_masks(nwbfile, plane) -> RoiMasks:
    """One plane's ROI footprints, from whichever mask column the file carries.

    ``pixel_mask`` is a ragged list of ``(x, y, weight)`` per ROI; ``image_mask`` is a
    dense per-ROI boolean image. Both reduce to the same pixel coordinates, so nothing
    downstream needs to know which was stored -- see docs/families/roi_position.md.

    Read per ROI rather than whole-column: a dense mask is (n_rois, 512, 512), which is
    far larger than the footprints it describes.
    """
    ps = _plane_segmentation(nwbfile, plane)
    colnames = list(getattr(ps, "colnames", []))
    n_rois = len(ps.id)
    rows, cols, counts = [], [], []
    weights_one = True

    if "pixel_mask" in colnames:
        source, shape = "pixel_mask", None
        column = ps["pixel_mask"]
        for i in range(n_rois):
            entry = np.asarray(column[i])
            if entry.dtype.names:                      # named fields: x, y, weight
                x, y = entry["x"], entry["y"]
                if "weight" in entry.dtype.names and not np.allclose(entry["weight"], 1.0):
                    weights_one = False
            else:                                      # unstructured (n, 3) as (x, y, w)
                flat = entry.reshape(len(entry), -1)
                x, y = flat[:, 0], flat[:, 1]
                if flat.shape[1] > 2 and not np.allclose(flat[:, 2], 1.0):
                    weights_one = False
            cols.append(np.asarray(x, dtype=np.int32))
            rows.append(np.asarray(y, dtype=np.int32))
            counts.append(len(x))

    elif "image_mask" in colnames:
        source = "image_mask"
        column = ps["image_mask"]
        data = getattr(column, "data", column)
        shape = tuple(int(v) for v in np.shape(data)[1:3]) or None
        for i in range(n_rois):
            yx = np.argwhere(np.asarray(data[i]) > 0)
            rows.append(yx[:, 0].astype(np.int32))
            cols.append(yx[:, 1].astype(np.int32))
            counts.append(len(yx))

    else:
        raise KeyError(
            f"plane {plane} has no ROI mask column; available: {sorted(colnames)}"
        )

    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    return RoiMasks(
        row=np.concatenate(rows) if rows else np.zeros(0, np.int32),
        col=np.concatenate(cols) if cols else np.zeros(0, np.int32),
        offsets=offsets, source=source, shape=shape, weights_all_one=weights_one)
