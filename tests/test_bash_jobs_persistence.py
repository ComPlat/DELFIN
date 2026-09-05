"""Persistent bash-job registry: crash-safe re-attach + completion events.

A dashboard/kernel restart used to orphan every running calculation: the
setsid child survived, but the new process had no record of job ids,
output-file paths, or start times, and the only completion signal was a
blocking bash_status poll. These tests cover the additive persistence
layer in delfin.agent.bash_jobs:

- <workspace>/.delfin/bash_jobs.json is written at job start and updated
  by the watchdog at exit (atomic temp-file + os.replace writes),
- drain_finished_events() reports each finished job exactly once, with
  the acknowledged flag surviving a simulated restart,
- unknown job ids re-attach from the registry file: a live pid reports
  running (pid-reuse guarded via /proc start ticks), a dead pid reports
  finished with exit_code=None (unrecoverable once init reaped it),
- records older than ~7 days are pruned opportunistically.

All jobs are short-lived real subprocesses — no SLURM, no LLM.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ


@pytest.fixture(autouse=True)
def _isolated_registry(tmp_path, monkeypatch):
    # Keep the per-user locator index out of the real ~/.delfin, and make
    # sure no in-memory job leaks between tests (process-wide singleton).
    monkeypatch.setattr(BJ, "_INDEX_PATH", tmp_path / "bash_jobs_index.json")
    BJ._REGISTRY._jobs.clear()
    yield
    BJ._REGISTRY._jobs.clear()


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _read_record(ws: Path, job_id: str) -> dict:
    data = json.loads((ws / ".delfin" / "bash_jobs.json").read_text())
    return data["jobs"][job_id]


def _wait_for_exit_record(ws: Path, job_id: str, timeout: float = 10.0) -> dict:
    # The watchdog writes exit_code + finished_at at the moment of exit.
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        rec = _read_record(ws, job_id)
        if rec.get("finished_at") is not None:
            return rec
        time.sleep(0.05)
    raise AssertionError(f"watchdog never recorded exit for {job_id}")


def _write_registry(ws: Path, records: dict) -> None:
    reg = ws / ".delfin" / "bash_jobs.json"
    reg.parent.mkdir(parents=True, exist_ok=True)
    reg.write_text(json.dumps({"jobs": records}))


def _record(ws: Path, job_id: str, **overrides) -> dict:
    rec = {
        "job_id": job_id,
        "pid": 999999999,                 # definitely dead
        "proc_start_ticks": None,
        "command": "orca big.inp",
        "description": "crashed session",
        "cwd": str(ws),
        "workspace": str(ws),
        "stdout_path": str(ws / f"{job_id}.stdout"),
        "stderr_path": str(ws / f"{job_id}.stderr"),
        "started_at": time.time() - 120,
        "timeout_s": 3600,
        "exit_code": None,
        "finished_at": None,
        "acknowledged": False,
    }
    rec.update(overrides)
    return rec


# ---------------------------------------------------------------------------
# Registry file lifecycle
# ---------------------------------------------------------------------------


def test_registry_file_written_on_start_and_updated_on_exit(workspace):
    job = BJ.get_registry().start(
        "sleep 0.5; echo done", cwd=str(workspace), workspace=workspace)
    rec = _read_record(workspace, job.job_id)
    assert rec["pid"] == job.proc.pid
    assert rec["command"].startswith("sleep 0.5")
    assert rec["cwd"] == str(workspace)
    assert rec["stdout_path"] == str(job.stdout_path)
    assert rec["stderr_path"] == str(job.stderr_path)
    assert rec["started_at"] == pytest.approx(time.time(), abs=30)
    assert rec["exit_code"] is None and rec["finished_at"] is None

    rec = _wait_for_exit_record(workspace, job.job_id)
    assert rec["exit_code"] == 0
    assert rec["finished_at"] >= rec["started_at"]


# ---------------------------------------------------------------------------
# Completion events: exactly once, restart-safe
# ---------------------------------------------------------------------------


def test_drain_returns_finished_job_exactly_once_and_survives_restart(workspace):
    job = BJ.get_registry().start(
        "echo out-line; echo err-line >&2; exit 3",
        cwd=str(workspace), workspace=workspace)
    _wait_for_exit_record(workspace, job.job_id)

    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == [job.job_id]
    ev = events[0]
    assert ev["exit_code"] == 3
    assert "out-line" in ev["stdout_tail"]
    assert "err-line" in ev["stderr_tail"]
    assert ev["runtime_s"] >= 0
    assert "echo out-line" in ev["command"]

    # Exactly once within this process...
    assert BJ.drain_finished_events(workspace) == []
    # ...and across a simulated restart: the in-memory dict is gone, the
    # acknowledged flag lives in the registry file.
    BJ._REGISTRY._jobs.clear()
    assert BJ.drain_finished_events(workspace) == []

    # New jobs after the "restart" still produce events.
    job2 = BJ.get_registry().start(
        "echo two", cwd=str(workspace), workspace=workspace)
    _wait_for_exit_record(workspace, job2.job_id)
    events2 = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events2] == [job2.job_id]
    assert events2[0]["exit_code"] == 0


def test_drain_reports_job_that_died_while_unattached(workspace):
    # Crash scenario: the watchdog never wrote the exit (agent process died
    # mid-run) and the child's pid is long gone -> the drain itself detects
    # the death and reports exit_code=None exactly once.
    Path(_record(workspace, "deadbeef")["stdout_path"]).write_text(
        "tail of a finished run\n")
    _write_registry(workspace, {"deadbeef": _record(workspace, "deadbeef")})

    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == ["deadbeef"]
    assert events[0]["exit_code"] is None            # unrecoverable
    assert "finished run" in events[0]["stdout_tail"]
    assert BJ.drain_finished_events(workspace) == []


def test_old_records_are_pruned_on_drain(workspace):
    old = _record(workspace, "aaaa0000",
                  started_at=time.time() - 8 * 24 * 3600)
    fresh = _record(workspace, "bbbb1111",
                    finished_at=time.time() - 5, exit_code=0)
    _write_registry(workspace, {"aaaa0000": old, "bbbb1111": fresh})

    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == ["bbbb1111"]   # old one: no event
    data = json.loads((workspace / ".delfin" / "bash_jobs.json").read_text())
    assert "aaaa0000" not in data["jobs"]                  # pruned from file


# ---------------------------------------------------------------------------
# Re-attach after a restart
# ---------------------------------------------------------------------------


def test_reattach_live_job_after_restart(workspace):
    job = BJ.get_registry().start(
        "sleep 30", cwd=str(workspace), workspace=workspace)
    jid = job.job_id
    time.sleep(0.1)
    BJ._REGISTRY._jobs.clear()                       # simulated restart

    # Lookup WITHOUT a workspace goes through the per-user locator index
    # (bash_status/bash_output/bash_kill carry no workspace argument).
    re_job = BJ.get_registry().get(jid)
    assert re_job is not None and re_job is not job
    status = re_job.status_dict()
    assert status["running"] is True
    assert status["reattached"] is True
    assert status["stdout_path"] == str(job.stdout_path)
    assert 0 <= status["elapsed_s"] < 60             # epoch-based, sane

    # Explicit workspace works too, and kill() reaches the re-attached pid.
    assert BJ.get_registry().get(jid, workspace) is not None
    ok, _msg = BJ.get_registry().kill(jid)
    assert ok is True
    deadline = time.monotonic() + 5
    while (BJ.get_registry().get(jid).poll() is None
           and time.monotonic() < deadline):
        time.sleep(0.05)
    assert BJ.get_registry().get(jid).status_dict()["running"] is False


def test_dead_pid_reported_finished_with_unknown_exit(workspace):
    rec = _record(workspace, "deadbeef")
    Path(rec["stdout_path"]).write_text("tail of a finished ORCA run\n")
    _write_registry(workspace, {"deadbeef": rec})

    job = BJ.get_registry().get("deadbeef", workspace)
    assert job is not None
    status = job.status_dict()
    assert status["running"] is False
    assert status["exit_code"] is None               # unrecoverable
    assert "unknown" in status.get("note", "")
    # The observed death was persisted for later status calls / drains.
    assert _read_record(workspace, "deadbeef")["finished_at"] is not None


def test_pid_reuse_guard_via_proc_start_ticks(workspace):
    # A recycled pid must not masquerade as the old job: our own (alive) pid
    # with a deliberately wrong recorded start time reads as finished.
    if BJ._proc_start_ticks(os.getpid()) is None:
        pytest.skip("/proc/<pid>/stat not available on this platform")
    rec = _record(workspace, "feedface",
                  pid=os.getpid(), proc_start_ticks=1)
    _write_registry(workspace, {"feedface": rec})

    job = BJ.get_registry().get("feedface", workspace)
    assert job.status_dict()["running"] is False


def test_unknown_job_id_stays_unknown(workspace):
    assert BJ.get_registry().get("ffffffff") is None
    assert BJ.get_registry().get("ffffffff", workspace) is None
