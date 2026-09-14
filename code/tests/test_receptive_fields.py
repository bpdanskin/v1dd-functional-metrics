"""receptive_field_metrics against synthetic locally-sparse-noise data.

Covers both methods: the default greedy pixelwise RF (bootstrap + Holm-Šidák, two
variants) and the historical fraction-threshold method (rf_method="fraction", what
REFERENCE_CONFIG uses), plus the method-agnostic window-containment geometry.
"""

import numpy as np
import pandas as pd

from support import check
from v1dd_metrics import nwb as vn
from v1dd_metrics.config import MetricConfig
from v1dd_metrics.families.receptive_fields import (
    _rf_pixel_to_degrees, holm_sidak_reject, receptive_field_metrics)
from v1dd_metrics.families.surround_suppression import (
    CONTAINMENT_COLUMNS, POP_RF_COLUMNS, _window_coverage, window_containment)
from v1dd_metrics.schema import OUTPUT_COLUMNS, to_output_schema

DT = 0.16504
ROWS, COLS, GRID = 8, 14, 9.3
N_FRAMES, PERIOD = 1600, 1.8
N_ROIS = 4
TARGET_R, TARGET_C = 3, 5


def _scenario(seed=19):
    """Synthetic LSN plane: ROI0 selective to one bright pixel, ROI1/3 silent, ROI2 responds
    to everything. Returns (plane, trials, spont, lsn, altitudes, azimuths)."""
    RNG = np.random.default_rng(seed)
    images = np.zeros((N_FRAMES, ROWS, COLS), dtype=np.int8)
    for f in range(N_FRAMES):
        on = RNG.choice(ROWS * COLS, 6, replace=False)
        off = RNG.choice(np.setdiff1d(np.arange(ROWS * COLS), on), 6, replace=False)
        images[f].flat[on] = 1
        images[f].flat[off] = -1
    # the driving pixel is bright on a SMALL fraction of frames — sparsity keeps non-RF
    # pixels below threshold (that is what "locally sparse noise" buys).
    bright_frames = np.arange(0, N_FRAMES, 20)
    images[bright_frames, TARGET_R, TARGET_C] = 1
    dark_only = np.setdiff1d(np.arange(N_FRAMES), bright_frames)
    images[dark_only, TARGET_R, TARGET_C] = np.where(
        images[dark_only, TARGET_R, TARGET_C] == 1, 0, images[dark_only, TARGET_R, TARGET_C])

    azimuths = (np.arange(COLS) - COLS // 2 + 0.5) * GRID
    altitudes = (np.arange(ROWS) - ROWS // 2 + 0.5) * GRID
    lsn = {"images": images, "azimuths": azimuths, "altitudes": altitudes,
           "pixel_on": 1, "pixel_off": -1, "pixel_gray": 0, "pixel_values": [-1, 0, 1]}

    starts = 30.0 + np.arange(N_FRAMES) * PERIOD
    trials = pd.DataFrame({"stim_name": "locally_sparse_noise", "start_time": starts,
                           "stop_time": starts + 0.30, "frame": np.arange(N_FRAMES).astype(float)})
    spont_start = starts[-1] + 5.0
    spont_stop = spont_start + 300.0
    ts = np.arange(0.0, spont_stop + 20.0, DT)
    dff = RNG.normal(0.0, 0.02, size=(len(ts), N_ROIS))
    for f, s in enumerate(starts):
        w = (ts >= s) & (ts <= s + 4 * DT)
        if images[f, TARGET_R, TARGET_C] == 1:
            dff[w, 0] += 2.0
        dff[w, 2] += 2.0
    evt = np.maximum(dff, 0.0)
    roi_table = pd.DataFrame({"column": 1, "volume": 3, "plane": 0,
                              "roi": np.arange(N_ROIS), "pika_roi_confidence": 0.9})
    plane = vn.PlaneData(mouse_id="409828", depth_um=150.0, column=1, volume="3", plane=0,
                         roi=np.arange(N_ROIS), is_valid=np.ones(N_ROIS, bool), timestamps=ts,
                         traces={"dff": dff, "events": evt}, roi_table=roi_table, dt=DT)
    return plane, trials, (spont_start, spont_stop), lsn, altitudes, azimuths


def test_holm_sidak_reject():
    # a clearly-significant test (p=0) is rejected; a null one (p large) is not
    p = np.array([[0.0], [0.5], [0.9]])
    rej = holm_sidak_reject(p, 0.05)
    check("holm-sidak rejects p=0", bool(rej[0, 0]))
    check("holm-sidak keeps p=0.5/0.9", not rej[1, 0] and not rej[2, 0])
    # step-down: the largest p (0.06 > alpha) fails its rank, so it alone is not rejected
    p2 = np.array([[1e-6], [0.06], [1e-7]])
    rej2 = holm_sidak_reject(p2, 0.05)
    check("step-down rejects the two tiny p, not the 0.06 one",
          bool(rej2[0, 0]) and bool(rej2[2, 0]) and not bool(rej2[1, 0]),
          str(rej2.ravel()))


def test_receptive_fields_greedy():
    plane, trials, spont, lsn, altitudes, azimuths = _scenario()
    cfg = MetricConfig(rf_greedy_n_boot=1000)          # default rf_method="greedy"

    print("[greedy] receptive_field_metrics")
    out, arr = receptive_field_metrics(plane, trials, spont, lsn, config=cfg,
                                       rng=np.random.default_rng(0))
    check("method is greedy", arr["method"] == "greedy")
    check("arrays carry sta/ge/strict_mask",
          {"sta", "ge", "strict_mask", "overlap_map", "overlap_thresh"} <= set(arr))
    check("sta shape", arr["sta"].shape == (N_ROIS, 2, ROWS, COLS), str(arr["sta"].shape))
    check("ge dtype uint16", arr["ge"].dtype == np.uint16, str(arr["ge"].dtype))
    check("all 14 RF columns present",
          {"has_rf_on", "has_rf_on_a05", "azimuth_rf_on", "azimuth_rf_on_a05"} <= set(out.columns))

    check("selective ROI0 has an ON receptive field", bool(out.has_rf_on[0]))
    check("silent ROIs 1 and 3 have none",
          not out.has_rf_on_or_off[1] and not out.has_rf_on_or_off[3])
    check("no-RF ROI gets NaN centres", bool(np.isnan(out.azimuth_rf_on[1])))

    print("\n[greedy] centre lands on the driving pixel")
    check("ON altitude within one pixel of the target",
          abs(out.altitude_rf_on[0] - altitudes[TARGET_R]) <= GRID + 1e-6,
          f"{out.altitude_rf_on[0]:.2f} vs {altitudes[TARGET_R]:.2f}")
    check("ON azimuth within one pixel of the target",
          abs(out.azimuth_rf_on[0] - azimuths[TARGET_C]) <= GRID + 1e-6,
          f"{out.azimuth_rf_on[0]:.2f} vs {azimuths[TARGET_C]:.2f}")

    print("\n[greedy] sensitive variant detects at least as much as strict")
    check("sens has_rf_on_or_off >= strict, per ROI",
          bool((out.has_rf_on_or_off_a05.to_numpy() >=
                out.has_rf_on_or_off.to_numpy()).all()))
    check("selective ROI0 detected by both variants",
          bool(out.has_rf_on[0]) and bool(out.has_rf_on_a05[0]))

    print("\n[greedy] schema and guardrails")
    pub = to_output_schema(out, "rf_metrics")
    check("column order matches published (14 RF cols)",
          list(pub.columns) == list(OUTPUT_COLUMNS["rf_metrics"]))
    check("has_rf_* (both variants) written as bool",
          all(pub[c].dtype == bool for c in
              ("has_rf_on", "has_rf_off", "has_rf_on_or_off",
               "has_rf_on_a05", "has_rf_off_a05", "has_rf_on_or_off_a05")))
    # empty trials -> absent frame + zero arrays of the right shape
    empty_out, empty_arr = receptive_field_metrics(plane, trials.iloc[:0], spont, lsn, config=cfg)
    check("empty trials -> absent frame, no RF", not empty_out.has_rf_on.any())
    check("empty trials -> greedy arrays right shape",
          empty_arr["sta"].shape == (N_ROIS, 2, ROWS, COLS))
    # frame index past the template raises
    bad = trials.copy(); bad.loc[0, "frame"] = 9999.0
    try:
        receptive_field_metrics(plane, bad, spont, lsn, config=cfg, rng=np.random.default_rng(0))
        check("raises on a frame index past the template", False)
    except ValueError as e:
        check("raises on a frame index past the template", "exceeds" in str(e))


def test_receptive_fields_fraction():
    plane, trials, spont, lsn, altitudes, azimuths = _scenario()
    cfg = MetricConfig(rf_method="fraction", other_n_boot=2000)

    print("[fraction] receptive_field_metrics")
    out, arr = receptive_field_metrics(plane, trials, spont, lsn, config=cfg,
                                       rng=np.random.default_rng(0))
    rf_map = arr["rf_map"]
    check("method is fraction", arr["method"] == "fraction")
    check("selective ROI0 has an ON receptive field", bool(out.has_rf_on[0]))
    check("silent ROIs have none",
          not out.has_rf_on[1] and not out.has_rf_off[1] and not out.has_rf_on_or_off[1])
    check("no-RF ROIs get NaN centres", bool(np.isnan(out.azimuth_rf_on[1])))
    check("_a05 columns are absent under fraction (bools False)",
          not out.has_rf_on_a05.any() and bool(np.isnan(out.azimuth_rf_on_a05[0])))

    print("\n[fraction] centre lands on the driving pixel (corrected interpolation)")
    check("altitude matches the target pixel",
          abs(out.altitude_rf_on[0] - altitudes[TARGET_R]) < 1e-9,
          f"{out.altitude_rf_on[0]:.4f} vs {altitudes[TARGET_R]:.4f}")
    check("azimuth matches the target pixel",
          abs(out.azimuth_rf_on[0] - azimuths[TARGET_C]) < 1e-9)

    # the historical scale-bug config compresses centres by (n-1)/n
    hist, _ = receptive_field_metrics(
        plane, trials, spont, lsn,
        config=MetricConfig(rf_method="fraction", other_n_boot=2000, rf_center_scale_bug=True),
        rng=np.random.default_rng(0))
    exp_alt = (TARGET_R + 0.5) * ((altitudes[-1] - altitudes[0]) / ROWS) + altitudes[0]
    exp_azi = (TARGET_C + 0.5) * ((azimuths[-1] - azimuths[0]) / COLS) + azimuths[0]
    check("historical config reproduces the compressed centre",
          abs(hist.altitude_rf_on[0] - exp_alt) < 1e-9 and abs(hist.azimuth_rf_on[0] - exp_azi) < 1e-9)
    check("corrected == historical * n/(n-1)",
          abs(out.altitude_rf_on[0] - hist.altitude_rf_on[0] * (ROWS / (ROWS - 1))) < 1e-9
          and abs(out.azimuth_rf_on[0] - hist.azimuth_rf_on[0] * (COLS / (COLS - 1))) < 1e-9)

    print("\n[fraction] the point_to_alt_azi scale bug")
    last_alt_bug = _rf_pixel_to_degrees(ROWS - 1, altitudes, True)
    last_azi_bug = _rf_pixel_to_degrees(COLS - 1, azimuths, True)
    check("buggy altitude at the last row is 28.481", abs(last_alt_bug - 28.4812) < 1e-3)
    check("buggy azimuth at the last column is 56.132", abs(last_azi_bug - 56.1321) < 1e-3)
    check("corrected altitude is 32.55",
          abs(_rf_pixel_to_degrees(ROWS - 1, altitudes, False) - 32.55) < 1e-9)
    check("corrected azimuth is 60.45",
          abs(_rf_pixel_to_degrees(COLS - 1, azimuths, False) - 60.45) < 1e-9)

    print("\n[fraction] pixel encoding — where a literal port breaks")
    wrong = dict(lsn, pixel_on=255, pixel_off=0)
    wrong_out, _ = receptive_field_metrics(plane, trials, spont, wrong, config=cfg,
                                           rng=np.random.default_rng(0))
    check("hard-coded pixel_on=255 erases every ON receptive field", not wrong_out.has_rf_on.any())
    missing = dict(lsn, pixel_on=None)
    try:
        receptive_field_metrics(plane, trials, spont, missing, config=cfg,
                                rng=np.random.default_rng(0))
        check("raises when the codes cannot be determined", False)
    except ValueError as e:
        check("raises when the codes cannot be determined", "pixel codes" in str(e))

    print("\n[fraction] rf_map shape and schema")
    check("rf_map shape", rf_map.shape == (N_ROIS, 2, ROWS, COLS), str(rf_map.shape))
    check("rf_map in [0,1]", float(rf_map.min()) >= 0.0 and float(rf_map.max()) <= 1.0)
    check("rf_map non-zero at the target pixel", rf_map[0, 0, TARGET_R, TARGET_C] > 0)
    pub = to_output_schema(out, "rf_metrics")
    check("column order matches published", list(pub.columns) == list(OUTPUT_COLUMNS["rf_metrics"]))


def test_window_containment():
    RADIUS = 15.0
    cfg_c = MetricConfig(dgw_window_radius_deg=RADIUS)
    alt_c = (np.arange(ROWS) - (ROWS - 1) / 2) * GRID
    azi_c = (np.arange(COLS) - (COLS - 1) / 2) * GRID
    lsn_c = {"altitudes": alt_c, "azimuths": azi_c}

    cov = _window_coverage(azi_c, alt_c, (0.0, 0.0), RADIUS)
    area = cov.sum() * GRID * GRID
    check("coverage recovers the disc area to <1%",
          abs(area - np.pi * RADIUS ** 2) / (np.pi * RADIUS ** 2) < 0.01)
    check("coverage is a fraction everywhere", bool(((cov >= 0) & (cov <= 1)).all()))
    check("centre pixel fully covered", cov[ROWS // 2, COLS // 2] == 1.0)
    check("far corner untouched", cov[0, 0] == 0.0)

    def containment(rf_map, centres, center=(0.0, 0.0), thresh=None):
        frame = pd.DataFrame({
            "azimuth_rf_on": [centres[0]], "altitude_rf_on": [centres[1]],
            "azimuth_rf_off": [np.nan], "altitude_rf_off": [np.nan]})
        return window_containment(frame, rf_map, lsn_c, center, config=cfg_c, thresh=thresh)

    m = np.zeros((1, 2, ROWS, COLS), dtype=np.float32)
    m[0, 0, ROWS // 2, COLS // 2] = 1.0
    near = (azi_c[COLS // 2], alt_c[ROWS // 2])
    c = containment(m, near)
    check("overlap = 1.0 for the pixel nearest the aperture centre",
          abs(float(c["dgw_rf_overlap_on"][0]) - 1.0) < 1e-12)
    check("distance = analytic hypot(pitch/2, pitch/2)",
          abs(float(c["dgw_rf_distance_on"][0]) - np.hypot(GRID / 2, GRID / 2)) < 1e-12)
    # population RF: a lone ON field -> pop RF is that field's centre; per-axis distance to
    # the aperture (here centred at 0,0) is |centre|
    check("pop RF azimuth equals the lone ON field's centre",
          abs(float(c["pop_rf_azimuth"][0]) - near[0]) < 1e-9, str(c["pop_rf_azimuth"][0]))
    check("pop RF per-axis distance = |pop - aperture|",
          abs(float(c["pop_rf_dis_azimuth"][0]) - abs(near[0])) < 1e-9
          and abs(float(c["pop_rf_dis_altitude"][0]) - abs(near[1])) < 1e-9)

    m = np.ones((1, 2, ROWS, COLS), dtype=np.float32)
    expect = (np.pi * RADIUS ** 2) / ((ROWS * GRID) * (COLS * GRID))
    c = containment(m, (0.0, 0.0))
    check("uniform map -> overlap = disc/screen",
          abs(float(c["dgw_rf_overlap_on"][0]) - expect) < 0.01 * expect)

    # greedy passes a binary strict mask with thresh=0.5; a lone significant pixel at the
    # centre must give overlap 1.0 exactly like the fraction case
    gm = np.zeros((1, 2, ROWS, COLS), dtype=np.float32)
    gm[0, 0, ROWS // 2, COLS // 2] = 1.0
    c = containment(gm, near, thresh=0.5)
    check("greedy binary mask (thresh=0.5): overlap 1.0 for the centre pixel",
          abs(float(c["dgw_rf_overlap_on"][0]) - 1.0) < 1e-12)

    # sub-threshold pixels do not dilute (fraction default thresh)
    m = np.zeros((1, 2, ROWS, COLS), dtype=np.float32)
    m[0, 0, ROWS // 2, COLS // 2] = 1.0
    m[0, 0, 0, 0] = float(cfg_c.rf_frac_thresh) - 0.01
    c = containment(m, (azi_c[COLS // 2], alt_c[ROWS // 2]))
    check("a sub-threshold pixel does not dilute the overlap",
          abs(float(c["dgw_rf_overlap_on"][0]) - 1.0) < 1e-12)
    m[0, 0, 0, 0] = float(cfg_c.rf_frac_thresh) + 0.01
    c = containment(m, (azi_c[COLS // 2], alt_c[ROWS // 2]))
    check("the same pixel above threshold dilutes it", float(c["dgw_rf_overlap_on"][0]) < 1.0)

    m = np.zeros((1, 2, ROWS, COLS), dtype=np.float32)
    c = containment(m, (np.nan, np.nan))
    check("no suprathreshold pixels -> overlap NaN", np.isnan(float(c["dgw_rf_overlap_on"][0])))
    c = containment(np.ones((1, 2, ROWS, COLS), dtype=np.float32), (0.0, 0.0), center=(np.nan, np.nan))
    check("unknown aperture -> all containment columns NaN",
          bool(c[CONTAINMENT_COLUMNS].isna().all().all()))

    try:
        _window_coverage(np.array([0.0, 9.3, 30.0]), alt_c, (0.0, 0.0), RADIUS)
        check("raises on an unevenly spaced stimulus grid", False)
    except ValueError as e:
        check("raises on an unevenly spaced stimulus grid", "evenly spaced" in str(e))

    check("containment then pop-RF columns close OUTPUT_COLUMNS",
          list(OUTPUT_COLUMNS["surround_suppression"][-8:]) == CONTAINMENT_COLUMNS + POP_RF_COLUMNS)
