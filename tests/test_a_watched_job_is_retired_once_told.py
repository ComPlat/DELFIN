"""A finished watch is reported once per reader and then leaves the file.

Two readers of ``check_agent_jobs`` kept a finished entry they never
retired:

* An unscoped consuming read (``consume=True``, no ``session_id``: the
  engine of a session without an id) does not pop an entry another session
  owns, but it did not mark it either, so it reported the same finished job
  on every turn.
* The owner's consuming read skips an entry the idle wake already told it
  (``wake_notified``) but left it in the file, where every following tick
  queried the scheduler for it until the seven-day prune.

Input: a fake scheduler and a watch file under tmp_path. Output: the
reports each reader receives, and the entries left in the file.
"""
from __future__ import annotations

import pytest

from delfin.agent import job_monitor as jm


def _scheduler(cmd):
    if cmd[0] == "squeue":
        return ""
    if cmd[0] == "sacct":
        return "123   COMPLETED\n"
    raise AssertionError(f"unexpected command: {cmd!r}")


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    d = tmp_path / "ws"
    d.mkdir()
    jm.register_agent_job(d, "123", description="opt",
                          extra={"kind": "slurm", "session_id": "S1"})
    return d


def _left(ws) -> list[str]:
    return list(jm.load_watched(jm._agent_watch_path(ws)).get("jobs", {}))


def test_an_unscoped_consuming_read_reports_an_owned_job_once(ws):
    first = jm.check_agent_jobs(ws, run_fn=_scheduler)
    second = jm.check_agent_jobs(ws, run_fn=_scheduler)
    assert [e["job_id"] for e in first] == ["123"]
    assert second == [], "the same finished job was reported again"
    # The owner is still owed it.
    owner = jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="S1")
    assert [e["job_id"] for e in owner] == ["123"]
    assert _left(ws) == []


def test_the_owner_retires_what_the_wake_already_told_it(ws):
    woke = jm.check_agent_jobs(ws, run_fn=_scheduler, consume=False,
                               marker="wake_notified", session_id="S1")
    drain = jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="S1")
    assert [e["job_id"] for e in woke] == ["123"]
    assert drain == []
    assert _left(ws) == [], "a told and finished watch stays polled"
