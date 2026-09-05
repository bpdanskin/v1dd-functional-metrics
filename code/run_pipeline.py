"""Capsule entry point: build the V1DD functional-metrics asset and describe it.

Stages run in a fixed order -- version gate, processing, validation, metadata --
because metadata records what validation found. See docs/pipeline.md.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from v1dd_metrics.version import code_version

CODE = Path(__file__).resolve().parent
RESULTS = Path(os.environ.get("V1DD_RESULTS_DIR", "/results"))
INPUT_ASSET = Path(os.environ.get("V1DD_INPUT_ASSET", "/data/409828_V1DD_Filtered"))
VALIDATION_DIR = Path(os.environ.get("V1DD_VALIDATION_DIR",
                                     "/scratch/v1dd_metrics_validation"))
ASSET_PREFIX = "V1DD_functional_metrics"


def describe_input_asset(path: Path) -> dict:
    """Identify the mounted asset from the sidecars its sessions carry.

    Which version of the input is mounted is not obvious from the path, and it has
    changed under this project before -- session directory names lost a ``_filtered_``
    suffix between one run and the next. Each session carries a ``data_description.json``
    whose ``name`` and ``creation_time`` say which build it came from, so read those
    rather than inferring from directory names.
    """
    sessions = sorted(d for d in path.iterdir()
                      if d.is_dir() and (d / "data_description.json").is_file())
    names, created = [], {}
    for session in sessions:
        try:
            raw = json.loads((session / "data_description.json").read_text(encoding="utf-8"))
        except Exception:                                           # noqa: BLE001
            continue
        names.append(raw.get("name") or session.name)
        stamp = str(raw.get("creation_time") or "")[:10]
        created[stamp] = created.get(stamp, 0) + 1
    return {"n_sessions": len(sessions), "names": names,
            "creation_dates": dict(sorted(created.items())),
            "filtered_suffix": sum("_filtered_" in n for n in names)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-env", action="store_true",
                    help="resolve the version and paths, report, and exit without running")
    ap.add_argument("--session", action="append", metavar="COLUMN:VOLUME", default=None,
                    help="restrict the run to these sessions, e.g. --session 1:3. "
                         "Repeatable. Recorded in provenance, so a partial asset cannot "
                         "pass for a complete one.")
    args = ap.parse_args(argv)

    session_filter = None
    if args.session:
        session_filter = []
        for spec in args.session:
            if ":" not in spec:
                raise SystemExit(f"--session expects COLUMN:VOLUME, got {spec!r}")
            col, vol = spec.split(":", 1)
            session_filter.append((int(col), str(vol)))

    version = code_version(CODE)
    print("code version : " + version)
    print("input asset  : " + str(INPUT_ASSET) + ("" if INPUT_ASSET.is_dir() else "   [MISSING]"))
    print("results      : " + str(RESULTS))
    print("validation   : " + str(VALIDATION_DIR))
    if session_filter:
        print(f"session filter: {session_filter}   [PARTIAL ASSET]")

    if args.check_env:
        if not INPUT_ASSET.is_dir():
            return 1
        info = describe_input_asset(INPUT_ASSET)
        print(f"\nmounted asset: {info['n_sessions']} session(s) with a sidecar")
        print(f"  names carrying _filtered_ : {info['filtered_suffix']}/{info['n_sessions']}")
        print(f"  sidecar creation dates    : {info['creation_dates']}")
        if info["names"]:
            print(f"  first session name        : {info['names'][0]}")
        # Off-schedule builds are the thread behind the two anomalous sessions; a single
        # date means one build, several means the asset was assembled in pieces.
        if len(info["creation_dates"]) > 1:
            print("  !! more than one build date -- sessions were not all made together")
        return 0

    RESULTS.mkdir(parents=True, exist_ok=True)
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    os.environ.update(V1DD_OUTPUT_TARGET="results", V1DD_CODE_VERSION=version)

    from v1dd_metrics import pipeline

    print("\n=== processing", flush=True)
    asset_dir = pipeline.run(input_asset=INPUT_ASSET, results_dir=RESULTS,
                             asset_prefix=ASSET_PREFIX, session_filter=session_filter)

    # Non-fatal: a failed check should surface in its artifacts rather than destroy an
    # asset that took hours to build. Before metadata, because processing.json records
    # that validation ran -- written first it could only describe one that had not.
    print("\n=== validation", flush=True)
    rc = subprocess.run([sys.executable, str(CODE / "validation" / "run_tests.py"),
                         str(VALIDATION_DIR)], env=os.environ).returncode
    if rc:
        print(f"!! validation reported failures ({rc}); the asset is still in {asset_dir}",
              flush=True)

    print("\n=== metadata", flush=True)
    subprocess.run(
        [sys.executable, "-u", str(CODE / "src" / "v1dd_metrics" / "metadata.py"),
         "--asset-dir", str(asset_dir), "--input-asset", str(INPUT_ASSET),
         "--results-dir", str(RESULTS), "--validation-dir", str(VALIDATION_DIR)],
        check=True, env=os.environ)

    print("\nasset: " + str(asset_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
