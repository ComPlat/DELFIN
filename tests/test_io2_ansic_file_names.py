"""IO2 round 15: ANSI-C quoting must not hide a real script file.

`sbatch $'/x/walk.sh'` was ALLOWED while `sbatch /x/walk.sh` refused
(measured at b1071469c). Cause: `shlex` does not expand ANSI-C quoting, so
the token arrives as `$/x/walk.sh`, the existing-file test finds nothing,
and the "a name that does not exist goes to sbatch, which fails by itself"
branch allows. But bash DOES expand `$'...'` before sbatch sees it, so
sbatch receives a real, walking file -- the gate's allow is a fail-open,
not the deliberate missing-name pass.

Both directions pinned:
  * `$'<real walking file>'`  -> refuse (the file's content is judged)
  * `$'<real clean file>'`    -> allow  (no over-refusal)
  * an ordinary path          -> unchanged in both directions
  * a `$'...'` naming nothing -> allow (sbatch fails on it, as before)
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def _scripts(tmp_path):
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    clean = tmp_path / "clean.sh"
    clean.write_text(CLEAN, encoding="utf-8")
    return walk, clean


def test_ansi_c_quoted_walking_file_is_refused(lustre, tmp_path):
    walk, _clean = _scripts(tmp_path)
    verdict, text, may_run = _v(f"sbatch $'{walk}'")
    assert verdict == "refuse", (
        f"FAIL-OPEN: `sbatch $'<walking file>'` was allowed "
        f"(verdict={verdict!r}); bash expands the quoting and sbatch runs "
        f"the file")
    assert may_run is False
    assert "tree_walk_on_lustre" in text


def test_ansi_c_quoted_walking_file_with_options_is_refused(lustre, tmp_path):
    walk, _clean = _scripts(tmp_path)
    verdict, text, may_run = _v(f"sbatch --test-only -p dev_cpu $'{walk}'")
    assert verdict == "refuse", (
        f"FAIL-OPEN with options (verdict={verdict!r})")
    assert may_run is False


def test_ansi_c_quoted_clean_file_still_runs(lustre, tmp_path):
    _walk, clean = _scripts(tmp_path)
    verdict, text, may_run = _v(f"sbatch $'{clean}'")
    assert verdict == "allow", (
        f"OVER-REFUSAL: `sbatch $'<clean file>'` was refused: {text[:160]}")
    assert may_run is True


def test_ansi_c_quoted_clean_file_with_options_still_runs(lustre, tmp_path):
    _walk, clean = _scripts(tmp_path)
    verdict, text, may_run = _v(f"sbatch -J probe -p dev_cpu $'{clean}'")
    assert verdict == "allow", (
        f"OVER-REFUSAL with options: {text[:160]}")
    assert may_run is True


def test_plain_paths_are_unchanged(lustre, tmp_path):
    walk, clean = _scripts(tmp_path)
    assert _v(f"sbatch {walk}")[0] == "refuse"
    assert _v(f"sbatch {clean}")[0] == "allow"


def test_ansi_c_quoting_a_missing_file_is_still_allowed(lustre, tmp_path):
    """The deliberate missing-name pass is preserved, not replaced."""
    verdict, text, may_run = _v(f"sbatch $'{tmp_path}/nosuch.sh'")
    assert verdict == "allow", (
        f"a name that exists nowhere should go to sbatch, which fails on "
        f"it; got {verdict!r}")
    assert may_run is True
