"""Package V3 -- phase 1 RED control (gap V3-1, reviewer-confirmed).

A watched job's completion must be delivered exactly once to its owning
session, however the session ends up reading it. Today it is delivered
TWICE when the same session reads it once idle and once busy:

  * idle peek: ``job_wake.finished_watched_jobs`` -> ``check_agent_jobs(
    consume=False, marker="wake_notified")`` keeps the entry and sets
    ``wake_notified`` (job_monitor.py:768);
  * busy drain: the engine's background block calls ``check_agent_jobs(
    consume=True, marker="daemon_notified")`` which appends the SAME
    completion again (job_monitor.py:752-766) and pops the entry -- the
    per-reader markers are not honoured by the consuming call.

This violates phase-2 "delivered exactly once". All three cases below are
RED on the unchanged tree (main 56b2eac2): double delivery happens.
"""
from __future__ import annotations

import pytest

from delfin.agent import job_monitor as jm

_SCHEDULE = {}


def _scheduler(cmd):
    if cmd[0] == "squeue":
        return "".join(f"{jid} RUNNING"
                       for jid, st in _SCHEDULE.items() if st == "RUNNING")
    if cmd[0] == "sacct":
        return "".join(f"{jid}   {st}" for jid, st in _SCHEDULE.items()
                       if st and st != "RUNNING")
    raise AssertionError(f"unexpected command: {cmd!r}")


@pytest.fixture(autouse=True)
def _reset(monkeypatch, tmp_path):
    _SCHEDULE.clear()
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")


@pytest.fixture
def ws(tmp_path):
    return tmp_path


def _register(ws, jid, *, session_id, state):
    jm.register_agent_job(ws, jid, description=f"job {jid}",
                          extra={"session_id": session_id})
    _SCHEDULE[jid] = state


def _peek_then_drain(ws, session_id="session-a"):
    """The reviewer's realistic trigger: a job finishes while the owning
    session is idle (peek), then the user types and the next turn's busy
    block drains. Both reads are for the same session."""
    peek = jm.check_agent_jobs(ws, run_fn=_scheduler, consume=False,
                               marker="wake_notified", session_id=session_id)
    drain = jm.check_agent_jobs(ws, run_fn=_scheduler, consume=True,
                                session_id=session_id)
    return peek, drain


def test_peek_then_drain_delivers_a_completion_once(ws):
    """One finished wat is announced once, not once per surface."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    peek, drain = _peek_then_drain(ws)
    assert [e["job_id"] for e in peek] == ["12345"]
    # The busy drain must NOT re-announce what the idle peek already
    # told this session.
    assert [e["job_id"] for e in drain] == []


def test_peek_then_drain_does_not_duplicate_a_failure(ws):
    """The failure case is the one that must never be whispered twice."""
    _register(ws, "67890", session_id="session-a", state="FAILED")
    peek, drain = _peek_then_drain(ws)
    assert [e["job_id"] for e in peek] == ["67890"]
    assert [e["job_id"] for e in drain] == []


def test_drain_alone_still_delivers(ws):
    """A session that skipped the idle peek still gets the completion
    from its busy drain -- the fix must not silence the lone reader."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    drain = jm.check_agent_jobs(ws, run_fn=_scheduler, consume=True,
                                session_id="session-a")
    assert [e["job_id"] for e in drain] == ["12345"]
