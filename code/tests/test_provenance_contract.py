"""The provenance record is metadata.py's input, so its keys are a contract.

A rename on either side once cost a published asset its record of whether validation had
run, and the tests at the time passed because they built a layout production never has.
"""

import json

import pandas as pd
import pytest

from v1dd_metrics import pipeline as pl
from v1dd_metrics import provenance as prov


@pytest.fixture
def record():
    sessions = pd.DataFrame([{"name": "s1", "format": "zarr", "column": 1,
                              "volume": "3", "n_planes": 6, "error": None}])
    inventory = sessions.copy()
    planes = pd.DataFrame([{"plane": p, "n_rois": 10} for p in range(6)])
    wide = pd.DataFrame({"roi_key": [f"r{i}" for i in range(60)]})
    return pl.build_provenance(
        asset_name="409828_V1DD_functional_metrics", stamp="2026-09-03_12-00-00",
        mouse_label="M409828", config=pl.DEFAULT_CONFIG, seed=0,
        input_asset=pl.Path("/data/409828_V1DD_Filtered"), sessions=sessions,
        inventory=inventory, planes=planes, wide=wide, wall_seconds=123.4,
        session_filter=None, failures=[], write_errors=[],
        window_centers={"n_measured": 1, "n_inferred": 0, "n_unfilled": 0},
        center_read_failures=[], wide_name="stimulus_metrics.parquet",
        manifest={}, arrays=["tuning_curves.npz"])


#: Every key ``metadata.py`` reads out of the provenance file.
METADATA_READS = ["generated_utc", "wall_seconds", "seed", "n_sessions", "n_planes",
                  "n_rois", "complete_asset", "session_filter", "config",
                  "differs_from_reference_config", "environment"]


def test_metadata_can_read_every_key_it_needs(record):
    missing = [k for k in METADATA_READS if k not in record]
    assert not missing, f"metadata.py reads keys the pipeline does not write: {missing}"
    assert "packages" in record["environment"]


def test_record_is_strict_json(record):
    """allow_nan=False, because a bare NaN literal is not valid JSON."""
    json.dumps(prov.jsonable(record), allow_nan=False, sort_keys=True)


def test_complete_asset_is_false_when_anything_went_wrong(record):
    assert record["complete_asset"] is True

    def rebuild(**kw):
        args = dict(
            asset_name="a", stamp="s", mouse_label="M409828", config=pl.DEFAULT_CONFIG,
            seed=0, input_asset=pl.Path("/data/x"),
            sessions=pd.DataFrame([{"name": "s1", "format": "zarr", "column": 1,
                                    "volume": "3", "n_planes": 6, "error": None}]),
            inventory=pd.DataFrame([{"error": None}]),
            planes=pd.DataFrame([{"plane": 0}]), wide=pd.DataFrame({"roi_key": ["r"]}),
            wall_seconds=1.0, session_filter=None, failures=[], write_errors=[],
            window_centers={}, center_read_failures=[], wide_name="t.parquet",
            manifest={}, arrays=[])
        args.update(kw)
        return pl.build_provenance(**args)

    assert rebuild(failures=[{"name": "s1"}])["complete_asset"] is False
    assert rebuild(write_errors=[{"output": "x.npz"}])["complete_asset"] is False
    assert rebuild(session_filter=[(1, "3")])["complete_asset"] is False


def test_every_array_that_landed_is_named(record):
    """An earlier asset shipped a 6 MB archive provenance never mentioned."""
    assert record["outputs"]["arrays"] == ["tuning_curves.npz"]
    assert record["outputs"]["table"] == "stimulus_metrics.parquet"


def test_metadata_reads_what_run_tests_writes(tmp_path):
    """The other half of the contract: validation's output must reach processing.json.

    `_validation_summary` returning None means the asset records one process where there
    should be two, silently. That has happened twice -- once from a path mismatch, once
    when the validation notebook was replaced by pytest and the lookup still expected the
    notebook's files. So drive the real writer's output through the real reader.
    """
    import sys
    sys.path.insert(0, str(pl.Path(__file__).resolve().parents[1] / "validation"))
    import run_tests
    from v1dd_metrics.metadata import _validation_summary

    # a JUnit report shaped like pytest's, with one failure and one skip
    (tmp_path / "tests.xml").write_text(
        '<?xml version="1.0" encoding="utf-8"?><testsuites><testsuite name="pytest" '
        'errors="0" failures="1" skipped="1" tests="3">'
        '<testcase classname="t" name="ok"/>'
        '<testcase classname="t" name="bad"><failure message="boom"/></testcase>'
        '<testcase classname="t" name="skipped"><skipped/></testcase>'
        '</testsuite></testsuites>', encoding="utf-8")
    summary = run_tests.summarise(tmp_path / "tests.xml")
    summary["seconds"] = 1.0
    (tmp_path / "tests.json").write_text(json.dumps(summary), encoding="utf-8")

    got = _validation_summary(tmp_path)
    assert got is not None, "metadata could not read what run_tests wrote"
    assert got["unit_tests_passed"] == 1, got
    assert got["unit_tests_failed"] == 1, got
    assert got["unit_tests_skipped"] == 1, got
    assert got["unit_tests_total"] == 3, got
    assert got["failed"] == ["t::bad"], got


def test_no_validation_artifacts_is_reported_not_guessed(tmp_path):
    from v1dd_metrics.metadata import _validation_summary
    assert _validation_summary(tmp_path) is None
    assert _validation_summary(None) is None
