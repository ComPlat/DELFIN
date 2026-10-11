"""Package V3 -- daemon/owner completion fence (reviewer s36 finding).

A finished session-owned watch must reach its OWNING session exactly once,
however the entries are read. This pins the failure reviewer `nacht-s36`
found: a daemon-side read that CONSUMES (``check_agent_jobs`` default is
``consume=True``; a bare call with ``session_id=None`` pops every entry it
sees) returns a session-a-owned completion to the daemon AND pops it, so
the owning session's later ``consume=True`` gets ``[]`` -- the only notice
of the finished calculation is silently eaten by another reader.

The fence: a consuming terminal branch pops an entry only when the caller
is a real session (``session_id is not None``) OR the entry has no owner
(dangling/daemon-owned). A session-owned entry is never popped by a
``session_id=None`` caller -- the owner is a SEPARATE reader and is owed
the completion. ``check_all_agent_jobs`` already defaults to
``consume=False`` for the real daemon; this closes the consuming path too.

RED on 717eb73b (owner's first consume returned []); GREEN with the
session-scoped pop fence on all three terminal branches (slurm/ci/bash).
Fakes only: an injected run_fn stands in for squeue/sacct, no sleeps.
"""
from __future__ import annotations

import pytest

from delfin.agent import job_monitor as jm

_SCHEDULE = {}
_RUNNING = "RUNNING"


def _scheduler(cmd):
    if cmd[0] in ("squeue", "sacct"):
        lines = [f"{jid} {state}" for jid, state in _SCHEDULE.items()
                 if state and state != _RUNNING]
        return "\n".join(lines)
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
                          extra={"kind": "slurm", "session_id": session_id})
    _SCHEDULE[jid] = state


def test_daemon_consume_peek_then_owner_consume_delivers_once(ws):
    """The daemon's CONSUMING read (bare check_agent_jobs, session_id=None)
    must not pop the entry the owning session-a is owed -- the owner's
    first consume delivers it, the second is silent."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    assert jm.check_agent_jobs(ws, run_fn=_scheduler) != []      # daemon peek
    first = jm.check_agent_jobs(ws, run_fn=_scheduler,
                                session_id="session-a", consume=True)
    assert first != [], "owner must get the completion once after the daemon peek"
    assert [e["job_id"] for e in first] == ["12345"]
    assert jm.check_agent_jobs(ws, run_fn=_scheduler,
                               session_id="session-a", consume=True) == []


def test_another_session_never_gets_the_owned_completion(ws):
    """The completion is scoped: session-b must never consume session-a's."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    assert jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a",
                               consume=True) != []
    assert jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-b",
                               consume=True) == []


def test_the_daemon_alone_consumes_an_unowned_entry(ws):
    """A dangling/unowned entry is nobody's: a session_id=None consuming
    read still pops it -- the ownerless completion is the daemon's to
    report and retire (no owner is waiting)."""
    jm.register_agent_job(ws, "99999", description="unowned job")
    _SCHEDULE["99999"] = "COMPLETED"
    assert jm.check_agent_jobs(ws, run_fn=_scheduler) != []
    assert jm.check_agent_jobs(ws, run_fn=_scheduler) == []
