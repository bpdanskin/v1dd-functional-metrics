"""Assembling the asset: the wide table, ROI coverage, and array-writer isolation.

These cover the orchestration around the metrics rather than the metrics themselves, so
they need no NWB input and run anywhere.
"""

import numpy as np
import pandas as pd
import pytest

from v1dd_metrics import pipeline as pl
from v1dd_metrics.schema import OUTPUT_COLUMNS


def family_frame(family: str, n_rois: int = 6, plane: int = 2) -> pd.DataFrame:
    """A minimal frame carrying every column ``family`` is expected to publish."""
    roi = np.arange(n_rois)
    df = pd.DataFrame({
        "roi_unique_id": [f"M409828_1{plane}_{r}" for r in roi],
        "roi_key": [f"M409828_1_3_{plane}_{r}" for r in roi],
        "mouse": "M409828", "column": 1, "volume": "3", "plane": plane, "roi": roi,
        "depth_um": 50.0 + 16 * plane, "pika_roi_confidence": 0.9,
    })
    for col in OUTPUT_COLUMNS[family]:
        if col not in df.columns:
            df[col] = np.linspace(0.0, 1.0, n_rois)
    return df


@pytest.fixture
def tables() -> dict:
    return {fam: family_frame(fam) for fam in pl.FAMILIES}


def test_build_wide_has_identity_once_and_every_metric(tables):
    wide, manifest = pl.build_wide(tables)

    assert len(wide) == 6, "one row per ROI"
    assert list(wide.columns[:len(pl.ID_COLS)]) == pl.ID_COLS
    assert wide.columns.is_unique

    # every family's metrics arrive, under that family's prefix
    for fam in pl.FAMILIES:
        expected = [c for c in OUTPUT_COLUMNS[fam] if c not in pl.ID_COLS]
        assert manifest[fam]["columns"] == [pl.PREFIX[fam] + c for c in expected]
        for col in manifest[fam]["columns"]:
            assert col in wide.columns

    # the identity block appears once, not once per family
    assert sum(c == "depth_um" for c in wide.columns) == 1


def test_build_wide_rejects_a_duplicate_key(tables):
    tables["natural_images"] = pd.concat(
        [tables["natural_images"], tables["natural_images"].iloc[:1]], ignore_index=True)
    with pytest.raises(Exception):
        pl.build_wide(tables)


def test_check_roi_coverage_counts_and_detects_a_gap(tables):
    assert pl.check_roi_coverage(tables) == 6

    tables["rf_metrics"] = tables["rf_metrics"].iloc[:-1]
    with pytest.raises(AssertionError, match="ROI set differs"):
        pl.check_roi_coverage(tables)


def test_guarded_writer_records_a_failure_and_lets_the_others_run():
    """A raising array writer must not stop the ones after it.

    Provenance is written after the array writers, so an unguarded raise costs a
    multi-hour run its record of everything else -- which is how a previous run lost
    three files.
    """
    written, errors = [], []
    guard = pl._guarded

    guard("first.npz", written, errors)(lambda: "first.npz")
    guard("second.npz", written, errors)(
        lambda: (_ for _ in ()).throw(ValueError("planes disagree on shape")))
    guard("third.npz", written, errors)(lambda: "third.npz")

    assert written == ["first.npz", "third.npz"], "the writer after the failure still ran"
    assert [e["output"] for e in errors] == ["second.npz"]
    assert "planes disagree on shape" in errors[0]["error"]
    assert "traceback" in errors[0]


def test_guarded_writer_treats_none_as_nothing_written():
    """A writer returning None had no data, which is not a failure."""
    written, errors = [], []
    pl._guarded("empty.npz", written, errors)(lambda: None)
    assert written == [] and errors == []


def test_config_dict_is_json_safe():
    from v1dd_metrics import provenance as prov
    d = pl.config_dict(pl.DEFAULT_CONFIG)
    assert isinstance(d["trace_type"], dict), "MappingProxyType must be copied out"
    import json
    json.dumps(prov.jsonable(d), allow_nan=False)


def test_reference_config_divergence_is_reported():
    """Every deliberate departure from the historical settings must be visible."""
    used = pl.config_dict(pl.DEFAULT_CONFIG)
    historical = pl.config_dict(pl.REFERENCE_CONFIG)
    differs = {k for k in used if used[k] != historical[k]}
    assert differs == {
        "rf_center_scale_bug", "pref_cond_fillna", "ni_response_frames",
        "fit_all_sf", "impute_dgw_center", "ssi_tuning_fit_includes_baseline",
        "lifetime_sparseness_over", "zero_denominator_nan",
    }, f"unexpected divergence set: {sorted(differs)}"


def test_a_family_that_never_ran_is_named(tables):
    """The failure the first capsule run hit, turned into a readable one.

    `roi_position` was wired into FAMILIES and the schema but its call was missing from
    the per-plane loop. Accumulator.tables() drops empty families, so build_wide raised a
    bare KeyError after 15 minutes of work -- naming the family but not the reason.
    """
    del tables["roi_position"]
    with pytest.raises(RuntimeError, match="no rows were produced"):
        pl.check_families_ran(tables)
    with pytest.raises(RuntimeError, match="roi_position"):
        pl.check_families_ran(tables)


def test_every_family_has_a_schema_and_a_prefix():
    """FAMILIES, OUTPUT_COLUMNS and PREFIX must agree, or the wide table cannot be built."""
    for fam in pl.FAMILIES:
        check_in = fam in OUTPUT_COLUMNS
        assert check_in, f"{fam} has no OUTPUT_COLUMNS entry"
        assert fam in pl.PREFIX, f"{fam} has no PREFIX entry"


def test_process_plane_appends_to_every_family():
    """Static guard: each family must be appended somewhere in the per-plane loop.

    Cheap, and it is the check that was missing -- the wiring was in FAMILIES, in the
    schema and in PREFIX, so everything except an actual run looked complete.
    """
    import inspect
    src = inspect.getsource(pl.process_plane)
    for fam in pl.FAMILIES:
        if fam in ("natural_images", "natural_images_12"):
            continue                      # appended through a loop variable
        assert f'"{fam}"' in src, f"process_plane never appends {fam}"
    assert 'acc.parts[fam].append' in src, "the natural-image loop should still be here"
