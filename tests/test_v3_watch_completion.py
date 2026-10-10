"""Package V3 -- phase 1 RED control (gap G2): a finished WATCHED job
leaves no persisted completion record for its owning session.

The shell half of the completion system persists an exactly-once, restart-safe
receipt (``bash_jobs.drain_finished_events`` / ``confirm_finished_events``,
claim/confirm into the registry file). A watched SLURM/CI job has no such
receipt: ``check_agent_jobs(consume=True)`` POPS the watch entry at the same
call that returns the completion (job_monitor.py:765-767), so a turn that
dies between the check and the delivery loses the only notice of the finished
calculation -- precisely what the shell half protects against.

This asserts the phase-2 requirement: after a finished watched job is
consumed, a persisted completion record carrying description + terminal
state and scoped to the OWNING session must be recoverable after a simulated
restart, exactly once -- through the same public path the shell half exposes
(``drain_watch_completions`` / ``confirm_watch_completions``). RED on the
unchanged tree: no such API exists.

Fakes only: an injected run_fn stands in for squeue/sacct; nothing real,
no sleeps.
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
    """The per-test workspace: the agent watch file lives under tmp_path."""
    return tmp_path


def _register(ws, jid, *, session_id, state, description="job desc"):
    jm.register_agent_job(ws, jid, description=description,
                          extra={"session_id": session_id})
    _SCHEDULE[jid] = state


def test_a_finished_watched_job_has_a_persisted_completion_record(ws):
    """Phase-2 requirement: consuming the live event must not destroy the
    only copy. A persisted receipt for the owning session survives restart;
    a neighbouring session's read must not consume it."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    done = jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a")
    assert [d["job_id"] for d in done] == ["12345"]

    own = jm.drain_watch_completions(ws, session_id="session-a")
    assert [r["job_id"] for r in own] == ["12345"]
    assert own[0]["description"] == "job desc"
    assert own[0]["state"] == "COMPLETED"

    # Exactly once for this session...
    assert jm.drain_watch_completions(ws, session_id="session-a") == []


def test_a_watched_completion_is_scoped_to_its_owning_session(ws):
    """The receipt belongs to the conversation that was waiting for it."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a")

    assert jm.drain_watch_completions(ws, session_id="session-b") == []
    own = jm.drain_watch_completions(ws, session_id="session-a")
    assert [r["job_id"] for r in own] == ["12345"]


def test_the_consume_pop_does_not_lose_a_failed_watch(ws):
    """A FAILED terminal state is the one that must never be lost."""
    _register(ws, "67890", session_id="session-a", state="FAILED")
    done = jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a")
    assert [d["job_id"] for d in done] == ["67890"]
    own = jm.drain_watch_completions(ws, session_id="session-a")
    assert [r["job_id"] for r in own] == ["67890"]
    assert own[0]["state"] == "FAILED"


def test_drain_reports_a_watched_completion_from_an_earlier_process(ws):
    """A restart is a new process: the drain re-reads the persisted file,
    never an in-memory copy. The store is file-backed, so a fresh read IS
    the restart; the only exclusivity is the confirm flag in the file."""
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a")

    # no in-process state: re-read the file as a new process would
    own = jm.drain_watch_completions(ws, session_id="session-a")
    assert [r["job_id"] for r in own] == ["12345"]
    jm.confirm_watch_completions(ws, own, session_id="session-a")
    # ...and the acknowledgement lives in the file, so a second process
    # (another drain) does not see it twice.
    assert jm.drain_watch_completions(ws, session_id="session-a") == []
