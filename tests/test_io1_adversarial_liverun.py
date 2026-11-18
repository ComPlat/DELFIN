"""Reviewer adversarial tests for the operator live-run findings L1-L5
(operator-A compute-node run @ 7a219348a, 2026-11-18). Each case is red
against 7a219348a; the builder's fix must turn them green. Uses only
injected mounts text (MountsSource fixture) or a monkeypatched default
path -- never the host's real /proc/mounts.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from delfin import hpc_io_check

LUSTRE = "nodefs /pfs/data6 lustre rw 0 0\n"
MOUNTS = hpc_io_check.MountsSource(text=LUSTRE)
ENV = {"HOME": "/pfs/data6/home/user"}


def _rules(script, env=ENV):
    return {f.rule for f in hpc_io_check.check(script, env=env, mounts_path=MOUNTS)}


def test_L1_default_mounts_source_reads_the_file(monkeypatch, tmp_path):
    """check() with NO mounts_path must read /proc/mounts via the file
    path, not treat the path string as mounts text (L1: the production
    wiring returned [] for every script)."""
    fake = tmp_path / "mounts"
    fake.write_text("nodefs /pfs/data6 lustre rw 0 0\n", encoding="utf-8")
    src = hpc_io_check.MountsSource(str(fake))
    assert "nodefs" in src.read(), "path-string source must be read as a FILE"
    findings = hpc_io_check.check(
        "for f in a b; do du -sh /pfs/data6/w; done\n",
        env=ENV, mounts_path=src)
    assert findings, "with a real mounts FILE the check must produce findings"


def test_L1_no_mounts_argument_produces_findings(monkeypatch, tmp_path):
    """The production call shape: check(script, env) with no mounts_path.
    A job script walking /pfs/data6 must be caught, not silently passed."""
    fake = tmp_path / "mounts"
    fake.write_text("nodefs /pfs/data6 lustre rw 0 0\n", encoding="utf-8")
    monkeypatch.setattr(hpc_io_check, "DEFAULT_MOUNTS", str(fake), raising=False)
    findings = hpc_io_check.check(
        "du -sh /pfs/data6/workspace\n", env=ENV)
    assert "tree_walk_on_lustre" in {f.rule for f in findings}


def test_L2_module_main_runs_and_exits_1_on_refusal(tmp_path):
    """python -m delfin.hpc_io_check <script> must work (L2: __main__
    block above the class definitions -> NameError) and exit 1 on a
    refusal. The mounts file is redirected so the test never reads the
    host's real /proc/mounts."""
    script = tmp_path / "walk.sh"
    script.write_text("du -sh /pfs/data6/workspace\n", encoding="utf-8")
    mounts = tmp_path / "mounts"
    mounts.write_text("nodefs /pfs/data6 lustre rw 0 0\n", encoding="utf-8")
    code = (
        "import sys, delfin.hpc_io_check as m\n"
        f"m.DEFAULT_MOUNTS = {str(mounts)!r}\n"
        "sys.argv = ['hpc_io_check', " + repr(str(script)) + "]\n"
        "raise SystemExit(m.main())\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, timeout=60)
    assert proc.returncode == 1, (
        f"expected exit 1 (refusal), got {proc.returncode}: "
        f"stdout={proc.stdout!r} stderr={proc.stderr!r}")


def test_L3_absolute_path_venv_python_in_loop_is_warned():
    """The 14.11. shape: a venv python called by ABSOLUTE path inside a
    loop (L3: m.group(0) = '/python ' defeats the re.match head)."""
    script = ("for f in a b c; do\n"
              "  /pfs/data6/home/user/env/bin/python job.py\n"
              "done\n")
    rules = _rules(script)
    assert "lustre_interpreter_in_loop" in rules, (
        "absolute-path venv python in a loop on Lustre must warn")


def test_L4_python_dash_m_pytest_from_home_is_refused():
    """The 16.11. shape: `cd <HOME checkout>; .venv/bin/python -m pytest`
    (L4: pytest only matched as the interpreter word)."""
    script = ("cd /pfs/data6/home/user/delfin\n"
              ".venv/bin/python -m pytest tests/\n")
    rules = _rules(script)
    assert "pytest_from_lustre_checkout" in rules, (
        "python -m pytest from a HOME checkout must refuse")


def test_L5_stagein_copy_to_tmpdir_is_not_refused():
    """cp -r of ONE job folder from the workspace to $TMPDIR is the FIX
    the rule asks for, not a tree walk (L5 false refusal)."""
    script = "cp -r /pfs/data6/workspace/runs/water $TMPDIR/\n"
    rules = _rules(script)
    assert "tree_walk_on_lustre" not in rules
    script2 = "rsync -a /pfs/data6/workspace/runs/water/ $TMPDIR/water/\n"
    assert "tree_walk_on_lustre" not in _rules(script2)


def test_L5_tree_walks_over_home_still_refused():
    """Guard against over-correction: du/find over HOME and a workspace
    root stay refused (bad_walk.sh must stay red)."""
    assert "tree_walk_on_lustre" in _rules("du -sh /pfs/data6/workspace\n")
    assert "tree_walk_on_lustre" in _rules(
        "find /pfs/data6/home/user -name '*.out'\n")
