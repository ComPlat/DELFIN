"""A finished background job in another folder was never reported.

``drain_finished_events(workspace)`` reads ONE workspace's registry, and
the engine calls it with the workspace the session happens to be in right
now. A job is registered under the workspace it was STARTED in.

Those are not always the same folder:

  - a subagent running in an isolated worktree starts the job there;
  - the user switches the office folder, or the workspace, mid-session;
  - a resumed session comes back with a different root.

In each case the job runs to completion, its registry record is updated
with the exit code, and nothing ever drains it. The one signal the engine
has that a multi-hour calculation finished never fires -- so the agent
either babysits it with blocking status polls or, more often, forgets it.

The locator index that makes this fixable already exists: every start
writes ``job_id -> workspace`` into ``~/.delfin/bash_jobs_index.json``,
pruned to the same seven days as the registries. Draining the workspaces
named there, in addition to the current one, costs one file read.

The exactly-once contract is per record and unchanged: the acknowledged
flag lives in the registry file, so a job drained from its own workspace
is not reported a second time from the sweep.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import bash_jobs as bj


@pytest.fixture
def two_workspaces(tmp_path):
    a = tmp_path / "here"
    b = tmp_path / "elsewhere"
    for p in (a, b):
        (p / ".delfin").mkdir(parents=True)
    return a, b


def _finished_job(ws, job_id: str, *, rc: int = 0) -> None:
    """Register an already-finished job under ``ws``.

    The timestamps have to be recent: records older than the seven-day
    horizon are pruned on read, before anything can be drained."""
    import time
    now = time.time()
    rec = {
        "job_id": job_id,
        "command": f"echo {job_id}",
        "description": "",
        "cwd": str(ws),
        "pid": 999_999_999,
        "started_at": now - 2.0,
        "finished_at": now,
        "exit_code": rc,
        "stdout_path": "",
        "stderr_path": "",
    }
    bj._persist_job_start(str(ws), rec)
    bj._note_job_workspace(job_id, str(ws))


def _ids(events) -> list[str]:
    return sorted(e["job_id"] for e in events)


# ---------------------------------------------------------------------------
# The sweep reaches the other folder
# ---------------------------------------------------------------------------

def test_a_job_from_another_workspace_is_reported(two_workspaces):
    here, elsewhere = two_workspaces
    _finished_job(elsewhere, "far")
    assert _ids(bj.drain_all_finished_events(here)) == ["far"]


def test_the_current_workspace_is_still_drained(two_workspaces):
    here, _ = two_workspaces
    _finished_job(here, "near")
    assert _ids(bj.drain_all_finished_events(here)) == ["near"]


def test_both_are_reported_together(two_workspaces):
    here, elsewhere = two_workspaces
    _finished_job(here, "near")
    _finished_job(elsewhere, "far")
    assert _ids(bj.drain_all_finished_events(here)) == ["far", "near"]


def test_the_exit_code_survives_the_sweep(two_workspaces):
    here, elsewhere = two_workspaces
    _finished_job(elsewhere, "failed", rc=3)
    assert bj.drain_all_finished_events(here)[0]["exit_code"] == 3


# ---------------------------------------------------------------------------
# ...exactly once, as before
# ---------------------------------------------------------------------------

def test_a_job_is_not_reported_twice(two_workspaces):
    here, elsewhere = two_workspaces
    _finished_job(elsewhere, "far")
    assert _ids(bj.drain_all_finished_events(here)) == ["far"]
    assert bj.drain_all_finished_events(here) == []


def test_draining_the_owner_first_silences_the_sweep(two_workspaces):
    """The acknowledged flag lives in the record, not in the caller."""
    here, elsewhere = two_workspaces
    _finished_job(elsewhere, "far")
    assert _ids(bj.drain_finished_events(elsewhere)) == ["far"]
    assert bj.drain_all_finished_events(here) == []


def test_the_same_workspace_named_twice_reports_once(two_workspaces):
    here, _ = two_workspaces
    _finished_job(here, "near")
    bj._note_job_workspace("near", str(here) + "/")
    assert _ids(bj.drain_all_finished_events(here)) == ["near"]


# ---------------------------------------------------------------------------
# ...and never at the cost of the turn
# ---------------------------------------------------------------------------

def test_a_workspace_that_no_longer_exists_is_skipped(two_workspaces):
    here, elsewhere = two_workspaces
    _finished_job(elsewhere, "far")
    bj._note_job_workspace("ghost", "/nonexistent/gone")
    assert _ids(bj.drain_all_finished_events(here)) == ["far"]


def test_an_unreadable_index_does_not_take_the_turn_down(two_workspaces):
    here, _ = two_workspaces
    _finished_job(here, "near")
    bj._INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    bj._INDEX_PATH.write_text("{not json", encoding="utf-8")
    assert _ids(bj.drain_all_finished_events(here)) == ["near"]


def test_an_empty_index_reports_the_current_workspace(two_workspaces):
    here, _ = two_workspaces
    _finished_job(here, "near")
    bj._INDEX_PATH.write_text(json.dumps({"jobs": {}}), encoding="utf-8")
    assert _ids(bj.drain_all_finished_events(here)) == ["near"]


def test_nothing_finished_is_no_events(two_workspaces):
    here, _ = two_workspaces
    assert bj.drain_all_finished_events(here) == []


# ---------------------------------------------------------------------------
# The engine uses the sweep, not the single-folder drain
# ---------------------------------------------------------------------------

def test_the_engine_calls_the_sweep():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "engine.py").read_text(encoding="utf-8")
    assert "drain_all_finished_events" in src, (
        "the engine still drains only the folder it happens to be in")
