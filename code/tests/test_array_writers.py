"""The three array archives: what they must contain, and what must make them refuse.

Replaces the fork's `test_tuning_export`, which read the processing notebook's JSON. The
writers are functions now, so they can be called directly.
"""

import numpy as np
import pandas as pd
import pytest

from support import check
from v1dd_metrics import pipeline as pl

N_DIR, N_SF, N_TRIALS, N_BLANK = 12, 2, 8, (8, 6, 7)


def accumulator(n_planes=3, n_rois=5, blanks=N_BLANK):
    """An accumulator filled as the per-plane loop would fill it."""
    rng = np.random.default_rng(0)
    acc = pl.Accumulator()
    acc.tuning_axes = (np.arange(0.0, 360.0, 30.0), np.array([0.02, 0.04]))
    acc.lsn_grid = {"altitudes": np.arange(8.0), "azimuths": np.arange(14.0)}
    keys = []
    for p in range(n_planes):
        for kind in ("dgw", "dgf"):
            part = acc.tuning[kind]
            part["trials"].append(
                rng.gamma(1.0, 1e-3, (n_rois, N_DIR, N_SF, N_TRIALS)).astype(np.float32))
            part["blank"].append(rng.gamma(1.0, 1e-3, (n_rois, blanks[p])).astype(np.float32))
            part["params"].append(rng.random((n_rois, N_SF, 6)).astype(np.float32))
            part["running"].append(rng.random((N_DIR, N_SF, N_TRIALS)).astype(np.float32))
            part["plane_key"].append(f"M409828_1_3_{p}")
        keys += [f"M409828_1_3_{p}_{r}" for r in range(n_rois)]
        acc.rf_maps.append(rng.random((n_rois, 2, 8, 14)).astype(np.float32))
    tables = {"drifting_gratings_windowed": pd.DataFrame({"roi_key": keys}),
              "rf_metrics": pd.DataFrame({"roi_key": keys})}
    return acc, tables


def test_tuning_curves_round_trips_with_everything_needed_to_read_it(tmp_path):
    acc, tables = accumulator()
    name = pl.write_tuning_curves(acc, tables, tmp_path, pl.DEFAULT_CONFIG)
    check("a file was written", name == "tuning_curves.npz")

    z = dict(np.load(tmp_path / name, allow_pickle=True))
    for key in ("roi_key", "plane_key", "directions", "spatial_frequencies", "trace_type",
                "dgw_trials", "dgw_blank", "dgw_n_blank", "dgw_params", "dgw_running",
                "dgf_trials", "dgf_blank", "dgf_n_blank", "dgf_params", "dgf_running"):
        check(f"{key} present", key in z, str(sorted(z)))

    check("trials are (rois, dir, sf, trials)",
          z["dgw_trials"].shape == (15, N_DIR, N_SF, N_TRIALS), str(z["dgw_trials"].shape))
    check("running has no ROI axis and keys on the plane",
          z["dgw_running"].shape == (3, N_DIR, N_SF, N_TRIALS) == (len(z["plane_key"]),) +
          z["dgw_running"].shape[1:], str(z["dgw_running"].shape))
    check("every roi_key sits under a plane_key",
          np.isin([k.rsplit("_", 1)[0] for k in z["roi_key"]], z["plane_key"]).all())


def test_ragged_blank_sweeps_are_padded_and_their_true_width_recorded(tmp_path):
    """Sessions ran 5-8 grey sweeps. Padding to the widest is the answer, not a surprise:
    `nanmean(blank, axis=1)` is then the baseline whatever the session ran, and
    `*_n_blank` tells a reader which columns are real."""
    acc, tables = accumulator(blanks=(8, 6, 7))
    pl.write_tuning_curves(acc, tables, tmp_path, pl.DEFAULT_CONFIG)
    z = dict(np.load(tmp_path / "tuning_curves.npz", allow_pickle=True))

    check("padded to the widest plane", z["dgw_blank"].shape == (15, 8),
          str(z["dgw_blank"].shape))
    check("true per-plane widths recorded", list(z["dgw_n_blank"]) == [8, 6, 7],
          str(z["dgw_n_blank"]))
    narrow = z["dgw_blank"][5:10]                     # the 6-sweep plane
    check("the pad is NaN, not zero", np.isnan(narrow[:, 6:]).all())
    check("and the real columns are finite", np.isfinite(narrow[:, :6]).all())


def test_it_refuses_when_planes_disagree_on_shape(tmp_path):
    acc, tables = accumulator()
    acc.tuning["dgw"]["trials"][1] = acc.tuning["dgw"]["trials"][1][:, :, :1, :]
    with pytest.raises(AssertionError, match="disagree on shape"):
        pl.write_tuning_curves(acc, tables, tmp_path, pl.DEFAULT_CONFIG)


def test_it_refuses_when_the_roi_axis_does_not_match_the_table(tmp_path):
    acc, tables = accumulator()
    tables["drifting_gratings_windowed"] = tables["drifting_gratings_windowed"].iloc[:-1]
    with pytest.raises(AssertionError, match="ROIs"):
        pl.write_tuning_curves(acc, tables, tmp_path, pl.DEFAULT_CONFIG)


def test_receptive_field_maps_carry_their_degree_axes_and_seed(tmp_path):
    """Without altitudes, azimuths and the seed the maps are uninterpretable."""
    acc, tables = accumulator()
    name = pl.write_rf_maps(acc, tables, tmp_path, seed=7)
    z = dict(np.load(tmp_path / name, allow_pickle=True))
    check("maps are (rois, 2, rows, cols)", z["rf_maps"].shape == (15, 2, 8, 14),
          str(z["rf_maps"].shape))
    check("ON/OFF axis is second", z["rf_maps"].shape[1] == 2)
    for key in ("roi_key", "altitudes", "azimuths", "seed"):
        check(f"{key} travels with the maps", key in z, str(sorted(z)))
    check("the seed is the one used", int(z["seed"]) == 7)


def test_a_writer_with_nothing_to_write_returns_none(tmp_path):
    """Empty is not a failure -- the guard must be able to tell them apart."""
    check("no tuning data", pl.write_tuning_curves(pl.Accumulator(), {}, tmp_path,
                                                   pl.DEFAULT_CONFIG) is None)
    check("no rf maps", pl.write_rf_maps(pl.Accumulator(), {}, tmp_path, 0) is None)
    check("no condition means", pl.write_condition_means(pl.Accumulator(), {}, tmp_path,
                                                         pl.DEFAULT_CONFIG) is None)
