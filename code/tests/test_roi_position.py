"""ROI centroids and size, and the two reconstructions of the column layout."""

import numpy as np
import pytest

from support import check
from v1dd_metrics import config as cfg
from v1dd_metrics.families import roi_position as rpm
from v1dd_metrics.nwb import RoiMasks
from v1dd_metrics.schema import OUTPUT_COLUMNS


def masks_from_footprints(footprints, source="pixel_mask"):
    """A RoiMasks from a list of (rows, cols) arrays."""
    rows = np.concatenate([f[0] for f in footprints]) if footprints else np.zeros(0, int)
    cols = np.concatenate([f[1] for f in footprints]) if footprints else np.zeros(0, int)
    offs = np.concatenate([[0], np.cumsum([len(f[0]) for f in footprints])])
    return RoiMasks(row=rows.astype(np.int32), col=cols.astype(np.int32),
                    offsets=offs.astype(np.int64), source=source)


def square(y0, x0, side):
    yy, xx = np.meshgrid(np.arange(y0, y0 + side), np.arange(x0, x0 + side), indexing="ij")
    return yy.ravel(), xx.ravel()


def test_centroid_area_and_radius_on_a_known_square():
    # a 4x4 block at rows 10-13, cols 20-23: centroid (11.5, 21.5), area 16
    got = rpm.centroids(masks_from_footprints([square(10, 20, 4)]))
    check("y centroid", got["roi_y_px"][0] == pytest.approx(11.5))
    check("x centroid", got["roi_x_px"][0] == pytest.approx(21.5))
    check("area is the pixel count", got["roi_area_px"][0] == 16)
    check("radius is the equivalent disc",
          got["roi_radius_px"][0] == pytest.approx(np.sqrt(16 / np.pi)))


def test_centroids_match_a_naive_per_roi_loop():
    """The reduceat pass must agree with the obvious implementation."""
    rng = np.random.default_rng(0)
    foot = [(rng.integers(0, 512, k), rng.integers(0, 512, k))
            for k in (171, 196, 228, 1, 40)]
    got = rpm.centroids(masks_from_footprints(foot))
    for i, (r, c) in enumerate(foot):
        check(f"roi {i} y", got["roi_y_px"][i] == pytest.approx(r.mean()))
        check(f"roi {i} x", got["roi_x_px"][i] == pytest.approx(c.mean()))
        check(f"roi {i} area", got["roi_area_px"][i] == len(r))


def test_an_empty_footprint_is_nan_not_zero():
    """A segmented ROI with no surviving pixels has no position, which is not (0, 0)."""
    got = rpm.centroids(masks_from_footprints(
        [square(0, 0, 3), (np.zeros(0, int), np.zeros(0, int)), square(5, 5, 2)]))
    check("empty ROI x is NaN", np.isnan(got["roi_x_px"][1]))
    check("empty ROI y is NaN", np.isnan(got["roi_y_px"][1]))
    check("empty ROI area is NaN", np.isnan(got["roi_area_px"][1]))
    check("its neighbours are unaffected",
          got["roi_x_px"][0] == pytest.approx(1.0)
          and got["roi_x_px"][2] == pytest.approx(5.5))


def test_the_two_mask_formats_give_the_same_centroid():
    """23 sessions store pixel_mask and 2 store image_mask; they must not disagree.

    Both reduce to the same pixel coordinates, so this is the property that lets the rest
    of the pipeline ignore which was stored.
    """
    foot = [square(10, 20, 4), square(100, 300, 6)]
    a = rpm.centroids(masks_from_footprints(foot, source="pixel_mask"))
    b = rpm.centroids(masks_from_footprints(foot, source="image_mask"))
    for key in a:
        check(f"{key} agrees across formats", np.allclose(a[key], b[key], equal_nan=True))


def test_published_offsets_are_the_papers_grid():
    """Four columns tiling 800 um, one at the centre."""
    assignment = {1: None, 2: 0, 3: 1, 4: 2, 5: 3}
    off = rpm.published_offsets(assignment)
    check("the centre column is the origin", off[1] == (0.0, 0.0))
    xs = sorted({o[0] for c, o in off.items() if c != 1})
    ys = sorted({o[1] for c, o in off.items() if c != 1})
    check("x offsets are +/- a quarter span", xs == [-200.0, 200.0], str(xs))
    check("y offsets are +/- a quarter span", ys == [-200.0, 200.0], str(ys))
    check("the four span the published grid",
          max(xs) - min(xs) == rpm.GRID_SPAN_UM / 2)


def test_assign_columns_recovers_a_perfect_grid():
    """On a synthetic layout the fit must be exact, with no anisotropy."""
    scale = 20.0                                   # um per degree
    centers = {7: (0.0, 0.0)}                       # the centre, deliberately not index 1
    for col, (fx, fy) in zip((2, 3, 4, 5), rpm.QUADRANTS):
        centers[col] = (fx * rpm.GRID_SPAN_UM / scale, fy * rpm.GRID_SPAN_UM / scale)
    fit = rpm.assign_columns(centers)
    check("finds the centre column", fit["centre_column"] == 7, str(fit["centre_column"]))
    check("recovers the scale", fit["um_per_degree"] == pytest.approx(scale))
    check("residual is ~0", fit["residual_um"] < 1e-6, str(fit["residual_um"]))
    check("no anisotropy", fit["anisotropy"] == pytest.approx(1.0))
    check("every non-centre column gets a distinct quadrant",
          sorted(v for v in fit["assignment"].values() if v is not None) == [0, 1, 2, 3])


def test_assign_columns_reports_anisotropy_when_one_axis_is_flat():
    """The real asset's apertures barely separate the columns in elevation.

    The fit must make that visible rather than absorbing it into a residual, because a
    flat axis means the second grid dimension is not constrained by the data.
    """
    centers = {1: (0.0, 0.0), 2: (-10.0, -1.0), 3: (-10.0, 1.0),
               4: (10.0, -1.0), 5: (10.0, 1.0)}
    fit = rpm.assign_columns(centers)
    check("anisotropy is large", fit["anisotropy"] > 5, str(fit["anisotropy"]))
    check("the flat axis has the inflated scale",
          fit["um_per_degree_per_axis"]["elevation"]
          > fit["um_per_degree_per_axis"]["azimuth"])
    check("spread is reported per axis",
          fit["aperture_spread_deg"]["azimuth"] > fit["aperture_spread_deg"]["elevation"])


def test_too_few_columns_declines_rather_than_guessing():
    fit = rpm.assign_columns({1: (0.0, 0.0), 2: (1.0, 1.0)})
    check("no scale is invented", fit["um_per_degree"] is None)
    check("and it says why", "need" in fit.get("note", ""), str(fit))


def test_metrics_frame_matches_the_published_schema():
    from support import synthetic_plane
    plane = synthetic_plane(n_rois=3)
    plane.column, plane.volume, plane.plane = 1, "3", 2
    plane.roi = np.arange(3)
    plane.roi_table = None
    masks = masks_from_footprints([square(10, 20, 4), square(30, 40, 5),
                                   square(50, 60, 6)])
    frame = rpm.roi_position_metrics(
        plane, masks, mouse="M409828",
        published={1: (0.0, 0.0)}, retinotopic={1: (0.0, 0.0)})
    for col in OUTPUT_COLUMNS["roi_position"]:
        check(f"{col} present", col in frame.columns, str(list(frame.columns)))
    check("micrometres scale from pixels",
          np.allclose(frame["roi_x_um"], frame["roi_x_px"] * cfg.DEFAULT_CONFIG.um_per_pixel))


def test_no_pixel_scale_means_no_micrometres():
    """um_per_pixel is inferred, not recorded; None must leave the um columns empty."""
    import dataclasses
    from support import synthetic_plane
    plane = synthetic_plane(n_rois=2)
    plane.column, plane.volume, plane.plane, plane.roi = 1, "3", 2, np.arange(2)
    plane.roi_table = None
    masks = masks_from_footprints([square(0, 0, 3), square(9, 9, 3)])
    frame = rpm.roi_position_metrics(
        plane, masks, mouse="M409828",
        config=dataclasses.replace(cfg.DEFAULT_CONFIG, um_per_pixel=None))
    check("pixels are still reported", np.isfinite(frame["roi_x_px"]).all())
    check("micrometres are NaN", frame["roi_x_um"].isna().all())
    check("and so are both anatomical frames",
          frame["roi_x_um_published"].isna().all()
          and frame["roi_x_um_retinotopic"].isna().all())
