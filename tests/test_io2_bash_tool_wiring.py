"""IO2: the bash tool's sbatch route goes through the Lustre check.

This is the wiring contract for the api_client.py hook (the px pair in
.gate/px/, applied by the operator -- that file is protected, so the
hook is not part of this branch's diff). The test drives the REAL tool
path, ``_DocToolExecutor.execute("bash", ...)`` -- the function the
hook sits in -- not the checker function. A stub ``sbatch`` first on
PATH writes a marker file and exits 0, so nothing real is submitted:

- a Lustre-walking script: the tool result carries code
  ``io_check_refused`` and the marker does NOT exist (the submit never
  reaches the sbatch binary);
- a clean script: the stub IS reached (no over-refusal), marker exists;
- a warn finding: the command runs and the warning is in
  ``payload["notes"]``;
- a checker crash: an ordinary command still runs (a crash in the
  checker must not refuse every bash command of every session), while
  an ``sbatch`` command is refused (fail closed where the command can
  submit).

Measured control, operator node run (2026-11-21, px/io2 264162ed1):
without the patch 3 failed / 153 passed (refused-walk reaches the
stub, warn has no notes, crash-on-sbatch not refused); with the patch
156 passed, twice. The three offline-stand-in points (calc_policy's
sbatch dialog stubbed, MountsSource read patched, shell isolation
resolved to the plain path because nested bwrap exits 126 inside the
test cage) are in the fixtures below, each with its reason.

WHICH LAYER EACH CASE EXERCISES (reviewer nacht-s2's point: do not
assume a green here proves more wiring than it does). The stubs make
the path ABOVE the hook permissive, not invisible:
- calc_policy's resource dialog and its sbatch script scan are
  stubbed to pass, so an sbatch command reaches _execute_bash instead
  of being diverted to the dialog. In production calc_policy may
  refuse first; this test measures the hook layer BELOW that.
- _DocToolExecutor.execute is the real entry (api_client.py
  _execute_bash is the function the px hunk sits in); hpc_io_check
  itself is real except where a case monkeypatches it to raise.
- The bash run resolves to the plain in-process path (cage off): the
  command really executes, which is what makes the marker file a
  witness. The production cage layer is exercised by the operator's
  node build, not by this file.
"""
from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest

from delfin.agent import calc_policy as _calc_policy
from delfin import hpc_io_check

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho done > /dev/null\n"
WARN = "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n"
ENV = {"HOME": "/home/x"}


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


@pytest.fixture
def stub_sbatch(tmp_path, monkeypatch):
    """A fake sbatch first on PATH: writes a marker, exits 0."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    stub = bindir / "sbatch"
    stub.write_text(
        "#!/bin/sh\n"
        f"echo submitted >> {tmp_path / 'marker.txt'}\n"
        "exit 0\n", encoding="utf-8")
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP
               | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ.get('PATH', '')}")
    # calc_policy judges every sbatch segment before the hook runs; make
    # it permissive OFFLINE so the IO2 refusal is the only thing measured
    # (a real scontrol/squeue call is refused in the gate cage).
    monkeypatch.setattr(_calc_policy, "partition_limits",
                        lambda name=None: {"MaxTime": "08:00:00",
                                           "MaxNodes": "4",
                                           "MaxCPUsPerNode": "64",
                                           "Default": "YES"})
    monkeypatch.setattr(_calc_policy, "queued_jobs", lambda: 0)
    # calc_policy.check_sbatch reads the script for orca/xtb lines and
    # refuses anything else; the IO2 gate is what is measured here, so
    # the resource half is stubbed offline ("stubbed", "") -- the walk
    # check still sees the script text unchanged.
    monkeypatch.setattr(_calc_policy, "check_sbatch",
                        lambda run, cwd, limits_of=None, jobs_now=None:
                        ("stubbed sbatch", ""))
    return bindir


def _run_bash(cmd: str, tmp_path: Path, monkeypatch=None) -> dict:
    """The real tool path the hook lives in. The cage is switched off
    because nested bwrap exits 126 inside the test cage and would hide
    both halves of the control -- with the cage off the command really
    runs, so the marker tells the refused from the submitted case.
    calc_policy's sbatch dialog is made permissive by the stub_sbatch
    fixture; the classifier second opinion needs no model in a test
    (None runner -> None)."""
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
    if monkeypatch is not None:
        # Resolve the shell isolation to the plain in-process path
        # instead of the real host probes (operator's tested pattern,
        # test_a2w_adversarial_hook).
        monkeypatch.setattr(A, "_bash_isolation_decision",
                            lambda perms, mode=None: ("none", "off"))
        monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
        monkeypatch.setattr(A, "_landlock_functional", lambda: False)
        # _guard_without_cage still wraps in the socket-guard seccomp
        # helper (nested cage -> exit 126 in the gate cage); the
        # documented escape hatch _in_process_cage itself reads is
        # DELFIN_PROCESS_CAGE=off, which makes the wrapper a no-op.
        monkeypatch.setenv("DELFIN_PROCESS_CAGE", "off")
    perms = KitToolPermissions(mode="default",
                               workspace=str(tmp_path))
    out = _DocToolExecutor().execute(
        "bash", {"command": cmd, "description": "io2 wiring check"}, perms)
    return json.loads(out)


def _marker(tmp_path: Path) -> bool:
    return (tmp_path / "marker.txt").exists()


def _script(tmp_path: Path, text: str, name: str = "run.sh") -> str:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_refused_walk_script_never_reaches_sbatch(lustre, stub_sbatch,
                                                  tmp_path, monkeypatch):
    script = _script(tmp_path, WALK)
    out = _run_bash(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert not _marker(tmp_path), (
        "a refused sbatch command must never reach the sbatch binary")


def test_clean_script_still_reaches_sbatch(lustre, stub_sbatch, tmp_path,
                                           monkeypatch):
    script = _script(tmp_path, CLEAN)
    out = _run_bash(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("exit_code") == 0, out
    assert _marker(tmp_path), "a clean submit must run"


def test_warn_finding_runs_and_rides_in_notes(lustre, stub_sbatch,
                                              tmp_path, monkeypatch):
    script = _script(tmp_path, WARN)
    out = _run_bash(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("exit_code") == 0, out
    assert _marker(tmp_path)
    notes = out.get("notes") or []
    assert any("warn" in n or "scratch" in n.lower() for n in notes), out


def test_checker_crash_does_not_block_ordinary_commands(lustre,
                                                        stub_sbatch,
                                                        tmp_path,
                                                        monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("io2 drill")

    monkeypatch.setattr(hpc_io_check, "sbatch_gate", _boom)
    out = _run_bash("echo ordinary", tmp_path, monkeypatch)
    assert out.get("exit_code") == 0, out
    assert out.get("code") != "io_check_refused", out


def test_checker_crash_still_refuses_sbatch(lustre, stub_sbatch,
                                            tmp_path, monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("io2 drill")

    monkeypatch.setattr(hpc_io_check, "sbatch_gate", _boom)
    script = _script(tmp_path, CLEAN)
    out = _run_bash(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert not _marker(tmp_path)


def test_cd_relative_walk_script_is_refused(lustre, stub_sbatch,
                                            tmp_path, monkeypatch):
    """The operator's live-run gap (2026-11-21): `cd <dir> &&
    sbatch walk.sh` resolves the RELATIVE script against the cd
    target. Unfixed, the gate resolved it against the process cwd,
    found nothing and passed the submit through to sbatch."""
    sub = tmp_path / "sub"
    sub.mkdir()
    _script(sub, WALK, name="walk.sh")
    out = _run_bash(f"cd {sub} && sbatch walk.sh", tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert not _marker(tmp_path), (
        "a refused cd-relative submit must never reach sbatch")
