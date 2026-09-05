"""Run the unit suite and write a machine-readable summary.

The run happens where the data is and the result gets read somewhere else, so the summary
is a file rather than a console tally: ``metadata.py`` reads it to record that checking
happened, and the capsule log is thousands of lines long by the time anyone looks.

    python code/validation/run_tests.py [output-dir]
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def summarise(xml_path: Path) -> dict:
    """Counts and failing test ids from a JUnit XML report."""
    root = ET.parse(xml_path).getroot()
    suite = root.find("testsuite") if root.tag == "testsuites" else root
    failed = []
    for case in suite.iter("testcase"):
        if case.find("failure") is not None or case.find("error") is not None:
            failed.append(f"{case.get('classname')}::{case.get('name')}")
    total = int(suite.get("tests", 0))
    n_err = int(suite.get("errors", 0)) + int(suite.get("failures", 0))
    n_skip = int(suite.get("skipped", 0))
    return {"n_tests": total, "n_failed": n_err, "n_skipped": n_skip,
            "n_passed": total - n_err - n_skip, "failed": failed}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    out_dir = Path(argv[0]) if argv else REPO
    out_dir.mkdir(parents=True, exist_ok=True)
    xml_path = out_dir / "tests.xml"

    started = time.time()
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--tb=short",
         f"--junit-xml={xml_path}", str(REPO / "code" / "tests")],
        cwd=str(REPO), text=True)

    summary = {"exit_code": proc.returncode, "seconds": round(time.time() - started, 1)}
    if xml_path.is_file():
        summary.update(summarise(xml_path))
    else:
        summary["error"] = "pytest produced no report"

    (out_dir / "tests.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out_dir / 'tests.json'}: "
          + ", ".join(f"{k}={summary[k]}" for k in
                      ("n_passed", "n_failed", "n_skipped") if k in summary))
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
