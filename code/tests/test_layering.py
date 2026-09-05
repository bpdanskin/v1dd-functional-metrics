"""The package must not depend on anything outside itself.

Replaces the fork's `test_import_boundary`, which grepped `code/utils` for imports of
`code/validation`. The shape changed -- the pipeline is a package now -- but the property
is the same: validation and tests may read the pipeline, never the reverse. That boundary
erodes on its own unless something checks it.
"""

import ast
import subprocess
import sys
from pathlib import Path

from support import REPO, check

PKG = REPO / "code" / "src" / "v1dd_metrics"
FORBIDDEN = ("validation", "tests", "run_pipeline", "support", "pytest")


def package_modules():
    return sorted(PKG.rglob("*.py"))


def test_package_imports_nothing_from_validation_or_tests():
    offenders = []
    for path in package_modules():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            for name in names:
                head = name.split(".")[0]
                if head in FORBIDDEN:
                    offenders.append(f"{path.name}: {name}")
    check("no package module imports validation, tests or the entry point",
          not offenders, str(offenders))


def test_package_imports_with_nothing_else_on_the_path():
    """The real check: import it in a subprocess whose path holds only code/src.

    A static scan misses a runtime `sys.path` insertion; this does not.
    """
    code = ("import importlib, sys; "
            "mods = ['config','common','schema','paths','provenance','version',"
            "'responses','nwb','metadata','pipeline']; "
            "[importlib.import_module('v1dd_metrics.' + m) for m in mods]; "
            "print('ok')")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(REPO), env={"PYTHONPATH": str(PKG.parent),
                                             "PATH": "", "SYSTEMROOT": "C:\\Windows"})
    check("every package module imports with only code/src on the path",
          out.returncode == 0 and "ok" in out.stdout,
          (out.stderr or out.stdout)[-400:])


def test_families_do_not_import_each_other_except_the_documented_one():
    """Surround suppression and roi_quality read drifting gratings; nothing else pairs up.

    Family-to-family imports are how a session missing one stimulus takes an unrelated
    family down with it, so the allowed set is written down rather than discovered.
    """
    allowed = {("surround_suppression", "drifting_gratings"),
               ("roi_quality", "drifting_gratings")}
    found = set()
    for path in sorted((PKG / "families").glob("*.py")):
        if path.stem == "__init__":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                found.add((path.stem, node.module.split(".")[-1]))
    check("no undocumented family-to-family import", found <= allowed,
          str(sorted(found - allowed)))
