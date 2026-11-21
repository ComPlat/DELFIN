"""IO2 round 10, Finding 4: the BACKGROUND bash route must be gated too.

The agent has two bash tools. The px hook was added to `_execute_bash`
(foreground) only, so the identical command typed into
`_execute_bash_background` started a job that ran sbatch with no Lustre
check at all: the foreground route refused, the background route
submitted. Found by QS lead nacht-s8 and confirmed by the operator with a
stub sbatch on a node.

This file drives the REAL background route -- `_DocToolExecutor.
_execute_bash_background` directly, not the foreground one and not the
helper -- so it measures the wiring, not the check. The stub sbatch
writes a marker on the filesystem, so the marker is the witness that the
command actually reached the binary:

- a walking script: `code == "io_check_refused"`, no marker, and no job;
- a clean script: the job starts and the stub writes the marker;
- a warn finding: the job starts and the warning is in `notes`;
- a checker crash: an ordinary command still starts a job (a bug in the
  checker is not a permission decision), an sbatch command is refused.

Red on the previous commit: `_execute_bash_background` had no IO2 call at
all, so every refusal case returned `{"status": "started"}` and the
marker appeared.
"""
from __future__ import annotations

import json
import os
import stat
import time
from pathlib import Path

import pytest

from delfin import hpc_io_check
from delfin.agent import calc_policy as _calc_policy

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho done > /dev/null\n"
WARN = "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n"


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
    monkeypatch.setattr(_calc_policy, "partition_limits",
                        lambda name=None: {"MaxTime": "08:00:00",
                                           "MaxNodes": "4",
                                           "MaxCPUsPerNode": "64",
                                           "Default": "YES"})
    monkeypatch.setattr(_calc_policy, "queued_jobs", lambda: 0)
    monkeypatch.setattr(_calc_policy, "check_sbatch",
                        lambda run, cwd, limits_of=None, jobs_now=None:
                        ("stubbed sbatch", ""))
    return bindir


def _run_background(cmd: str, tmp_path: Path, monkeypatch) -> dict:
    """The real BACKGROUND tool path -- the function the hook now lives in."""
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
    # Same offline stand-ins as the foreground wiring test: the cage cannot
    # nest inside the test cage (exit 126) and the host probes are not what
    # this file measures.
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("none", "off"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setenv("DELFIN_PROCESS_CAGE", "off")
    perms = KitToolPermissions(mode="default", workspace=str(tmp_path))
    out = _DocToolExecutor()._execute_bash_background(
        {"command": cmd, "description": "io2 background wiring check"},
        perms)
    return json.loads(out)


def _marker(tmp_path: Path, wait_s: float = 5.0) -> bool:
    """The stub's marker. A background job needs a moment to reach exec."""
    deadline = time.monotonic() + wait_s
    while time.monotonic() < deadline:
        if (tmp_path / "marker.txt").exists():
            return True
        time.sleep(0.05)
    return (tmp_path / "marker.txt").exists()


def _script(tmp_path: Path, text: str, name: str = "run.sh") -> str:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_background_walking_script_is_refused(lustre, stub_sbatch, tmp_path,
                                              monkeypatch):
    script = _script(tmp_path, WALK)
    out = _run_background(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert "status" not in out, out
    assert not _marker(tmp_path, wait_s=1.0), (
        "a refused background sbatch must never reach the sbatch binary")


def test_background_clean_script_starts_and_reaches_sbatch(lustre,
                                                          stub_sbatch,
                                                          tmp_path,
                                                          monkeypatch):
    script = _script(tmp_path, CLEAN)
    out = _run_background(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("status") == "started", out
    assert out.get("job_id"), out
    assert _marker(tmp_path), "a clean background submit must run"


def test_background_warn_finding_rides_in_notes(lustre, stub_sbatch, tmp_path,
                                                monkeypatch):
    script = _script(tmp_path, WARN)
    out = _run_background(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("status") == "started", out
    notes = out.get("notes") or []
    assert any("warn" in n or "scratch" in n.lower() for n in notes), out


def test_background_checker_crash_refuses_sbatch(lustre, stub_sbatch,
                                                tmp_path, monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("io2 background drill")

    monkeypatch.setattr(hpc_io_check, "sbatch_gate", _boom)
    script = _script(tmp_path, CLEAN)
    out = _run_background(f"sbatch {script}", tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert not _marker(tmp_path, wait_s=1.0)


def test_background_checker_crash_does_not_block_ordinary_commands(
        lustre, stub_sbatch, tmp_path, monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("io2 background drill")

    monkeypatch.setattr(hpc_io_check, "sbatch_gate", _boom)
    out = _run_background("echo ordinary", tmp_path, monkeypatch)
    assert out.get("status") == "started", out
    assert out.get("code") != "io_check_refused", out


def test_background_heredoc_submit_is_refused(lustre, stub_sbatch, tmp_path,
                                             monkeypatch):
    """The round-10 shapes reach this route too, not only the foreground."""
    cmd = ("sbatch <<'EOF'\n"
           "#!/bin/bash\n"
           "du -sh /pfs/data6/workspace/\n"
           "EOF")
    out = _run_background(cmd, tmp_path, monkeypatch)
    assert out.get("code") == "io_check_refused", out
    assert not _marker(tmp_path, wait_s=1.0)
