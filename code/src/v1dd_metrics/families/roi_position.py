"""Where each ROI sits: centroid and size in the imaging plane, plus two candidate
placements of the imaging columns in a shared anatomical frame.

Plane-local position is solid -- it comes straight from the segmentation footprint. The
column offsets are not: the NWB records no usable geometry, so both offsets here are
reconstructions and they are shipped side by side rather than reconciled. See
docs/families/roi_position.md.
"""

from __future__ import annotations

import itertools
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import roi_frame

ROI_POSITION_COLUMNS = ["roi_x_px", "roi_y_px", "roi_area_px", "roi_radius_px",
                        "roi_x_um", "roi_y_um",
                        "roi_x_um_published", "roi_y_um_published",
                        "roi_x_um_retinotopic", "roi_y_um_retinotopic"]

#: Columns that tile the grid, and the one at its centre. From the white paper's imaging
#: strategy: "Four columns create a grid covering 800 X 800 um, while the 5th column
#: captured the center." Which data column index is the centre is not stated in the text
#: -- it is inferred here and reported, see ``assign_columns``.
GRID_SPAN_UM = 800.0
N_GRID_COLUMNS = 4

#: Quadrant offsets of a 2x2 grid covering ``GRID_SPAN_UM``, in (x, y) micrometres.
QUADRANTS = ((-0.25, -0.25), (-0.25, 0.25), (0.25, -0.25), (0.25, 0.25))


def centroids(masks) -> Dict[str, np.ndarray]:
    """Centroid, pixel area and equivalent radius per ROI, from a ``RoiMasks``.

    The centroid is the unweighted mean over in-mask pixels. That is not a choice this
    asset forces us to make: ``pixel_mask`` weights are all exactly 1 and ``image_mask``
    is binary, so weighted and unweighted agree.
    """
    off = np.asarray(masks.offsets)
    counts = np.diff(off).astype(np.int64)
    n_rois = len(counts)
    empty = counts == 0

    # bincount over an ROI index rather than reduceat over offsets: reduceat returns the
    # *element* at a zero-width segment instead of 0, so an ROI with no surviving pixels
    # both gets a bogus centroid and corrupts its neighbours.
    roi_index = np.repeat(np.arange(n_rois), counts)
    sum_col = np.bincount(roi_index, weights=np.asarray(masks.col, float),
                          minlength=n_rois)
    sum_row = np.bincount(roi_index, weights=np.asarray(masks.row, float),
                          minlength=n_rois)
    with np.errstate(invalid="ignore", divide="ignore"):
        safe_n = np.where(empty, 1.0, counts.astype(float))
        x = np.where(empty, np.nan, sum_col / safe_n)
        y = np.where(empty, np.nan, sum_row / safe_n)
    area = np.where(empty, np.nan, counts.astype(float))
    return {"roi_x_px": x, "roi_y_px": y, "roi_area_px": area,
            "roi_radius_px": np.sqrt(area / np.pi)}


def published_offsets(assignment: Mapping[int, Optional[int]]) -> Dict[int, Tuple[float, float]]:
    """Per-column (x, y) offsets in micrometres from the paper's nominal layout.

    ``assignment`` maps a column index to a quadrant index 0-3, or to None for the column
    at the centre.
    """
    out: Dict[int, Tuple[float, float]] = {}
    for col, quad in assignment.items():
        if quad is None:
            out[col] = (0.0, 0.0)
        else:
            fx, fy = QUADRANTS[quad]
            out[col] = (fx * GRID_SPAN_UM, fy * GRID_SPAN_UM)
    return out


def assign_columns(centers: Mapping[int, Sequence[float]]) -> dict:
    """Work out which column is the centre and which quadrant each other one occupies.

    The paper puts the layout in a figure, not in words, so it is recovered here from the
    only spatial information the asset carries: the windowed-grating aperture position,
    which the paper says was "determined separately to align with the population receptive
    fields of imaged neurons" for each column. Retinotopy in V1 is monotonic, so relative
    aperture position stands in for relative cortical position.

    The centre column is taken to be the one closest to the centroid of the others, and
    the remaining four are matched to quadrants by exhaustive search over all 24
    permutations. Returns the assignment, the fitted micrometres-per-degree, and the
    residual -- **read the residual before trusting the assignment.**
    """
    cols = sorted(k for k, v in centers.items()
                  if v is not None and np.isfinite(v[0]) and np.isfinite(v[1]))
    if len(cols) < N_GRID_COLUMNS + 1:
        return {"assignment": {c: None for c in cols}, "um_per_degree": None,
                "residual_um": None,
                "note": f"need {N_GRID_COLUMNS + 1} columns with a centre, have {len(cols)}"}

    pts = {c: np.asarray(centers[c], dtype=float) for c in cols}
    # the centre column is the one nearest the centroid of the rest
    def offness(c):
        others = np.stack([pts[k] for k in cols if k != c])
        return float(np.linalg.norm(pts[c] - others.mean(axis=0)))
    centre = min(cols, key=offness)
    rest = [c for c in cols if c != centre]

    rel = {c: pts[c] - pts[centre] for c in rest}
    best = None
    for perm in itertools.permutations(range(N_GRID_COLUMNS)):
        target = np.stack([np.array(QUADRANTS[q]) * GRID_SPAN_UM for q in perm])
        source = np.stack([rel[c] for c in rest])
        # one isotropic scale, no rotation: degrees to micrometres
        denom = float((source * source).sum())
        scale = float((source * target).sum() / denom) if denom else 0.0
        resid = float(np.sqrt(((source * scale - target) ** 2).sum(axis=1).mean()))
        if best is None or resid < best[0]:
            best = (resid, perm, scale)

    resid, perm, scale = best
    assignment: Dict[int, Optional[int]] = {centre: None}
    assignment.update({c: q for c, q in zip(rest, perm)})

    # Per-axis scales, reported because a single isotropic one hides the thing that
    # matters: the apertures separate the columns far better in azimuth than in
    # elevation, so this layout is close to one-dimensional. A large ratio here means the
    # second grid axis is not actually constrained by the data.
    target = np.stack([np.array(QUADRANTS[q]) * GRID_SPAN_UM for q in perm])
    source = np.stack([rel[c] for c in rest])
    per_axis = []
    for j in (0, 1):
        den = float((source[:, j] ** 2).sum())
        per_axis.append(float((source[:, j] * target[:, j]).sum() / den) if den else None)
    spread = [float(np.ptp(source[:, j])) for j in (0, 1)]

    return {"assignment": assignment, "centre_column": centre,
            "um_per_degree": scale, "residual_um": resid,
            "um_per_degree_per_axis": {"azimuth": per_axis[0], "elevation": per_axis[1]},
            "aperture_spread_deg": {"azimuth": spread[0], "elevation": spread[1]},
            "anisotropy": (per_axis[1] / per_axis[0]
                           if per_axis[0] and per_axis[1] else None),
            "centre_offness_deg": {c: round(offness(c), 3) for c in cols}}


def retinotopic_offsets(centers: Mapping[int, Sequence[float]], centre_column: int,
                        um_per_degree: float) -> Dict[int, Tuple[float, float]]:
    """Per-column offsets implied by the aperture retinotopy alone, in micrometres.

    Unsnapped: this is the measured layout scaled into micrometres, not fitted to the
    paper's grid. Comparing it against ``published_offsets`` is the point.
    """
    if um_per_degree is None or centre_column not in centers:
        return {}
    base = np.asarray(centers[centre_column], dtype=float)
    out = {}
    for col, pos in centers.items():
        p = np.asarray(pos, dtype=float)
        if not np.isfinite(p).all():
            continue
        d = (p - base) * um_per_degree
        out[col] = (float(d[0]), float(d[1]))
    return out


def roi_position_metrics(
    plane,
    masks,
    *,
    config: MetricConfig = DEFAULT_CONFIG,
    mouse: Optional[str] = None,
    published: Optional[Mapping[int, Sequence[float]]] = None,
    retinotopic: Optional[Mapping[int, Sequence[float]]] = None,
) -> pd.DataFrame:
    """Per-ROI position and size for one plane.

    ``masks`` is a ``nwb.RoiMasks``. ``published`` and ``retinotopic`` are per-column
    ``(x, y)`` offsets in micrometres; both are optional and their columns come out NaN
    when absent, because neither is recorded in the file and both are reconstructions.
    """
    out = {c: np.full(plane.n_rois, np.nan) for c in ROI_POSITION_COLUMNS}
    if masks is not None:
        if masks.n_rois != plane.n_rois:
            raise ValueError(
                f"masks cover {masks.n_rois} ROIs, plane has {plane.n_rois}")
        out.update(centroids(masks))

        scale = config.um_per_pixel
        if scale:
            out["roi_x_um"] = out["roi_x_px"] * scale
            out["roi_y_um"] = out["roi_y_px"] * scale
            for name, table in (("published", published), ("retinotopic", retinotopic)):
                off = (table or {}).get(int(plane.column))
                if off is not None:
                    out[f"roi_x_um_{name}"] = out["roi_x_um"] + float(off[0])
                    out[f"roi_y_um_{name}"] = out["roi_y_um"] + float(off[1])

    frame = roi_frame(plane, mouse=mouse)
    for key, value in out.items():
        frame[key] = value
    return frame
