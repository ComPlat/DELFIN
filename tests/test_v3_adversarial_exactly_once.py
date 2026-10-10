"""Adversarial tests: the two-phase completion-delivery contract.

Package V3 (background shells/watchers) promises: one completion record
per finished job/watch, persisted, delivered exactly once to the session
that started it (also after a restart), carrying exit code, tail of
output and the original purpose.

This file attacks the two-phase delivery at the heart of that promise --
``drain_finished_events`` CLAIMS an event, ``confirm_finished_events``
ACKNOWLEDGES it, and an unconfirmed claim expires after
``_EVENT_CLAIM_GRACE_S`` so a turn that dies between drain and confirm
re-delivers instead of losing the only notice of a finished run.

The claim/grace path is deliberately exercised with fabricated dead-pid
records (no real subprocess, no sleeps beyond the grace monkeypatch), so
the tests are deterministic and run inside the write cage.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ


@pytest.fixture(autouse=True)
def _isolated_registry(tmp_path, monkeypatch):
    monkeypatch.setattr(BJ, "_INDEX_PATH", tmp_path / "bash_jobs_index.json")
    BJ._REGISTRY._jobs.clear()
    yield
    BJ._REGISTRY._jobs.clear()


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _write_registry(ws: Path, records: dict) -> None:
    reg = ws / ".delfin" / "bash_jobs.json"
    reg.parent.mkdir(parents=True, exist_ok=True)
    reg.write_text(json.dumps({"jobs": records}))


def _read_record(ws: Path, job_id: str) -> dict:
    data = json.loads((ws / ".delfin" / "bash_jobs.json").read_text())
    return data["jobs"][job_id]


def _finished_record(ws: Path, job_id: str, *, machine: dict | None = None,
                     session_id: str = "") -> dict:
    """A finished job whose pid is long gone (init reaped it).

    ``machine`` mirrors what ``this_machine()`` records; when given, the
    drain's ``record_runs_elsewhere`` must still see it as OURS, so the
    machine goes unpinned (the per-user index is monkeypatched away)."""
    rec = {
        "job_id": job_id,
        "pid": 999999999,                 # definitely dead
        "proc_start_ticks": None,
        "command": "bench --heavy",
        "description": "the original purpose",
        "cwd": str(ws),
        "workspace": str(ws),
        "stdout_path": str(ws / f"kit_bg_{job_id}.stdout"),
        "stderr_path": str(ws / f"kit_bg_{job_id}.stderr"),
        "started_at": time.time() - 60,
        "timeout_s": 3600,
        "exit_code": 7,
        "finished_at": time.time() - 30,
        "acknowledged": False,
    }
    if session_id:
        rec["session_id"] = session_id
    (ws / f"kit_bg_{job_id}.stdout").write_text("tail of done output\n")
    (ws / f"kit_bg_{job_id}.stderr").write_text("tail of stderr\n")
    return rec
def test_drain_returns_finished_job_then_confirm_retires_it(workspace):
    """A finished job is reported; after confirm it never returns, even
    across a simulated restart (the acknowledged flag lives in the file)."""
    _write_registry(workspace, {"aaaa1111": _finished_record(workspace, "aaaa1111")})

    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == ["aaaa1111"]
    ev = events[0]
    assert ev["exit_code"] == 7
    assert ev["description"] == "the original purpose"
    assert "tail of done output" in ev["stdout_tail"]
    assert "tail of stderr" in ev["stderr_tail"]

    BJ.confirm_finished_events(events)

    # Acknowledged: no second report this process, or after a restart.
    assert BJ.drain_finished_events(workspace) == []
    BJ._REGISTRY._jobs.clear()
    assert BJ.drain_finished_events(workspace) == []


def test_drained_but_unconfirmed_event_returns_after_grace(workspace, monkeypatch):
    """The claim is a hold, not a retirement: a turn that drains but never
    confirms (dies before delivering) loses nothing -- after the grace the
    event comes back. This is the exactly-once-across-restart guarantee."""
    _write_registry(workspace, {"bbbb2222": _finished_record(workspace, "bbbb2222")})

    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == ["bbbb2222"]
    # Unconfirmed. While the claim is fresh the event is NOT re-reported...
    assert BJ.drain_finished_events(workspace) == []

    # ...but once the grace expires it is offered again -- also to a fresh
    # process that knows nothing about the first drain (a restart).
    monkeypatch.setattr(BJ, "_EVENT_CLAIM_GRACE_S", -1.0)
    BJ._REGISTRY._jobs.clear()
    again = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in again] == ["bbbb2222"]


def test_confirm_of_an_already_acknowledged_event_is_safe(workspace):
    """confirm is idempotent: acknowledging twice, or an id that is already
    gone from the file, must not raise and must not resurrect an event."""
    _write_registry(workspace, {"cccc3333": _finished_record(workspace, "cccc3333")})
    events = BJ.drain_finished_events(workspace)
    assert BJ.confirm_finished_events(events) == 1
    assert BJ.confirm_finished_events(events) == 0
    assert BJ.confirm_finished_events(events) == 0    # already acked
    assert BJ.drain_finished_events(workspace) == []


def test_confirm_of_a_junk_event_cannot_raise(workspace):
    """Confirm handed garbage (a non-list, wrong shapes) is a no-op, not a
    crash -- the deliver path is best-effort and never raises."""
    assert BJ.confirm_finished_events(None) == 0
    assert BJ.confirm_finished_events([None, "junk", 42]) == 0
    assert BJ.confirm_finished_events([{"job_id": "x", "workspace": ""}]) == 0


def test_a_claimed_event_still_reports_once_even_across_restart(workspace):
    """The 'exactly once' contract holds across a restart in the happy path:
    drain -> confirm in one living process, a second process sees nothing."""
    _write_registry(workspace, {"dddd4444": _finished_record(workspace, "dddd4444")})
    first = BJ.drain_finished_events(workspace)
    BJ.confirm_finished_events(first)
    BJ._REGISTRY._jobs.clear()
    assert BJ.drain_finished_events(workspace) == []


def test_orphan_pid_anonymous_completion_claims_exit_none(workspace):
    """A job whose pid vanished while unattached reports exit_code=None
    (unrecoverable) but still reaches the owning session once."""
    rec = _finished_record(workspace, "eeee5555")
    rec["exit_code"] = None
    rec["finished_at"] = None          # watchdog never wrote the exit
    _write_registry(workspace, {"eeee5555": rec})
    events = BJ.drain_finished_events(workspace)
    assert [e["job_id"] for e in events] == ["eeee5555"]
    assert events[0]["exit_code"] is None
