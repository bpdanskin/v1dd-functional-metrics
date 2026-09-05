"""The version gate, exercised the way the capsule runs it.

Replaces the fork's `test_entrypoint`. The ladder changed -- there is no Dockerfile ENV
line any more -- but the property has not: a run that cannot name its own commit must
refuse to start, because a reproducible run copies `code/` without `.git` and an asset
with no version cannot be tied to the code that made it.
"""

import shutil
import subprocess
import sys

import pytest

from support import REPO, check

ENTRY = REPO / "code" / "run_pipeline.py"
SRC = REPO / "code" / "src"


def run_entry(env_extra=None, cwd=None, entry=None, src=None):
    env = {"PYTHONPATH": str(src or SRC), "PATH": "", "SYSTEMROOT": "C:\\Windows"}
    env.update(env_extra or {})
    return subprocess.run([sys.executable, str(entry or ENTRY), "--check-env"],
                          capture_output=True, text=True, cwd=str(cwd or REPO), env=env)


def test_an_explicit_version_is_used_verbatim():
    out = run_entry({"V1DD_CODE_VERSION": "abc1234"})
    check("the environment value wins", "abc1234" in out.stdout, out.stdout + out.stderr)


def test_a_blank_version_is_treated_as_unset_not_as_malformed(tmp_path, monkeypatch):
    """Whitespace means "nobody set this", so the ladder continues to the next rung.

    Distinct from a malformed value, which stops the run: an empty string carries no
    claim about which commit built the asset, whereas `"= <sha>"` carries a false one.
    """
    monkeypatch.setenv("V1DD_CODE_VERSION", "   ")
    (tmp_path / "CODE_VERSION").write_text("feedface\n", encoding="utf-8")
    from v1dd_metrics.version import code_version
    check("falls through to the file", code_version(tmp_path, strict=True) == "feedface")


@pytest.mark.parametrize("bad", ["= 17cacea", "not-a-sha", "abc", "zzzzzzz"])
def test_a_value_that_is_not_a_commit_is_refused(bad):
    """Being set is not the same as being usable.

    A malformed Docker `ENV <key> <value>` line once made this `"= <sha>"`, which is
    neither empty nor whitespace, so a guard that only asked "is it set?" shipped it.
    """
    out = run_entry({"V1DD_CODE_VERSION": bad})
    check(f"refuses {bad!r}", out.returncode != 0, out.stdout + out.stderr)
    check("and says why", "not a commit" in (out.stdout + out.stderr)
          or "unknown" in (out.stdout + out.stderr), out.stdout + out.stderr)


def test_the_shipped_version_file_asserts_nothing():
    """`code/CODE_VERSION` ships comment-only, so it can never claim a stale commit.

    Skipped where it has been stamped: writing the commit into it before a capture run is
    exactly what the file is for, so a stamped file is a deployment state rather than a
    defect. This is a claim about what the repository ships, not about any runtime.
    """
    from v1dd_metrics.version import read_version_file
    stamped = read_version_file(REPO / "code" / "CODE_VERSION")
    if stamped:
        pytest.skip(f"CODE_VERSION has been stamped ({stamped[:12]}) -- "
                    "a deployment state, not something to assert against")
    check("no version is asserted until someone stamps it", stamped == "")


def test_a_stamped_version_file_is_read(tmp_path, monkeypatch):
    """The file is only consulted when the environment does not answer first.

    The environment must be cleared explicitly: the entry point exports
    ``V1DD_CODE_VERSION`` before launching validation, so in a capsule this test runs with
    it already set and would otherwise be reading the environment, not the file.
    """
    monkeypatch.delenv("V1DD_CODE_VERSION", raising=False)
    (tmp_path / "CODE_VERSION").write_text("# a comment\n\nfeedface\n", encoding="utf-8")
    from v1dd_metrics.version import code_version
    check("first non-comment line is the version",
          code_version(tmp_path, strict=True) == "feedface")


def test_it_refuses_in_a_git_less_copy_with_no_version(tmp_path):
    """The shape a reproducible run actually has.

    Three provenance defects in earlier runs passed their tests and failed in the capsule
    because the tests ran in a shape production never has. This is that shape.
    """
    code = tmp_path / "code"
    shutil.copytree(REPO / "code", code,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    # Control the fixture rather than inherit it. `code/CODE_VERSION` ships comment-only,
    # but a capsule stamps it before a run -- and a copy carrying a stamped version is not
    # the situation under test.
    (code / "CODE_VERSION").write_text(
        "# blanked by the test: this copy must have no version\n", encoding="utf-8")
    check("the copy really has no .git above it", not (tmp_path / ".git").exists())
    check("and no version stamped into it",
          not [ln for ln in (code / "CODE_VERSION").read_text(encoding="utf-8").splitlines()
               if ln.strip() and not ln.startswith("#")])

    out = run_entry({}, cwd=tmp_path, entry=code / "run_pipeline.py", src=code / "src")
    check("refuses rather than shipping a null version", out.returncode != 0,
          out.stdout + out.stderr)
    check("and points at the two ways to set one",
          "V1DD_CODE_VERSION" in (out.stdout + out.stderr)
          and "CODE_VERSION" in (out.stdout + out.stderr), out.stdout + out.stderr)

    out = run_entry({"V1DD_CODE_VERSION": "0ab9142"}, cwd=tmp_path,
                    entry=code / "run_pipeline.py", src=code / "src")
    check("but starts once told", "0ab9142" in out.stdout, out.stdout + out.stderr)


def _in_git_checkout() -> bool:
    out = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True)
    return out.returncode == 0


@pytest.mark.skipif(not _in_git_checkout(),
                    reason="no git repository -- the fallback cannot apply here")
def test_the_git_fallback_answers_in_a_checkout():
    """Skipped rather than failed where it cannot apply.

    The fork's version of this failed in the capsule, and the validation notebook printed
    "unit tests failed" over a clean asset -- twice, because the first fix could only skip
    whole files. Marking the single inapplicable check is what pytest gives for free.
    """
    from v1dd_metrics.version import code_version
    sha = code_version(REPO / "code", strict=True)
    check("resolves from git here", sha and len(sha) >= 7, str(sha))
    check("looks like a full sha", len(sha) == 40, str(sha))
