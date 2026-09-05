"""Capsule entry point: build the V1DD functional-metrics asset and describe it.

Stages run in a fixed order -- version gate, processing, validation, metadata --
because metadata records what validation found. See docs/pipeline.md.
"""

from __future__ import annotations

import argparse
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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-env", action="store_true",
                    help="resolve the version and paths, report, and exit without running")
    args = ap.parse_args(argv)

    version = code_version(CODE)
    print("code version : " + version)
    print("input asset  : " + str(INPUT_ASSET) + ("" if INPUT_ASSET.is_dir() else "   [MISSING]"))
    print("results      : " + str(RESULTS))
    print("validation   : " + str(VALIDATION_DIR))

    if args.check_env:
        return 0 if INPUT_ASSET.is_dir() else 1

    RESULTS.mkdir(parents=True, exist_ok=True)
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    os.environ.update(V1DD_OUTPUT_TARGET="results", V1DD_CODE_VERSION=version)

    from v1dd_metrics import pipeline

    print("\n=== processing", flush=True)
    asset_dir = pipeline.run(input_asset=INPUT_ASSET, results_dir=RESULTS,
                             asset_prefix=ASSET_PREFIX)

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
