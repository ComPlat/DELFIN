"""A finished watched job's state says which code it ended with.

WHY THIS FILE EXISTS. A real run (probe ``failed-job-not-called-done``, on
this branch) started ``echo hi; exit 3`` with ``bash_background``, watched
it, and was handed two notices in one turn:

    [background] job 89dc4f48 (Run exact test command that exits 3)
    finished with exit code 3. Read its output with bash_output ...

    # Background jobs finished since your last turn (act on these results now)
    - bash 89dc4f48 [FINISHED] Watch test job that exits 3

The second line is the one a hurried reader takes at a glance, and it says
"finished" and nothing else about a job that FAILED. The base tree's single
line said ``exit code 3``; so on this one point the branch was weaker than
main, which is why this file exists.

The information is not missing at the source. ``check_agent_jobs`` builds
its entry with ``exit_code`` and ``ok`` in the same dict
(``job_monitor.py``, the ``job.status_dict()`` branch), and then sets

    "state": ("FINISHED (children still running)"
              if status.get("children_running") else "FINISHED")

-- the state, which is the only field the renderer shows besides the
description, drops the code the dict is carrying. So the state must say
which code it ended with, in the same reading order the shell row uses:
``FINISHED`` for a clean exit, ``FINISHED, exit N`` for a non-zero one.

This is measured through ``check_agent_jobs`` itself, on a real job started
through the dispatch, not by constructing an entry: an assertion about a
dict we build ourselves would pass whatever the code did.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import bash_jobs as bj
from delfin.agent import job_monitor
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def job_registry():
    bj._REGISTRY._jobs.clear()
    yield bj
    bj._REGISTRY._jobs.clear()


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _finished_watched_job(ws, command: str) -> tuple[str, int]:
    """Start, wait for, and watch a job through the model's own calls."""
    perms = KitToolPermissions(workspace=ws, mode="default")
    started = json.loads(_doc_executor.execute(
        "bash_background", {"command": command, "description": "the work"},
        perms))
    assert started.get("status") == "started", started
    job_id = started["job_id"]

    status = {}
    for _ in range(400):
        status = json.loads(_doc_executor.execute(
            "bash_status", {"job_id": job_id}, perms))
        if not status.get("running"):
            break
        time.sleep(0.02)
    assert status.get("running") is False, status

    watched = json.loads(_doc_executor.execute(
        "watch_job", {"job_id": job_id}, perms))
    assert watched.get("status") == "watching", watched
    return job_id, status["exit_code"]


def _entry(ws, job_id):
    """The watch entry ``check_agent_jobs`` reports for this job."""
    for _ in range(200):
        for ev in job_monitor.check_agent_jobs(
                ws, consume=False, marker="wake_notified") or []:
            if ev.get("job_id") == job_id:
                return ev
        time.sleep(0.02)
    return None


def test_a_failing_watched_job_state_names_the_code(job_registry, workspace):
    """The input that makes this fail: a job that ends non-zero.

    ``exit 3`` cannot end 0, so the assertion below runs on a failing code
    wherever it runs -- the literal is read from the job, never assumed.
    """
    job_id, exit_code = _finished_watched_job(
        workspace, "echo nothing-is-run-here; exit 3")

    ev = _entry(workspace, job_id)
    assert ev is not None, "the watched job was never reported by job_monitor"
    assert exit_code not in (0, None), (
        f"the job did not fail as the test needs it to: exit_code={exit_code}")
    assert str(exit_code) in str(ev.get("state") or ""), (
        "the state of a job that ended " + str(exit_code) + " does not say "
        "so, so the line built from it reads as a clean finish: "
        + repr(ev.get("state")))


def test_a_clean_watched_job_state_stays_quiet(job_registry, workspace):
    """The control: a job that ended well must not gain a code to explain.

    Host-independent by construction: whatever code the job ended with
    (0 on a host that runs it, 126 in a sandbox that refuses it), the state
    must still read as ``FINISHED`` first -- the fix adds the code for a
    failure without replacing the word -- and must never spell out
    ``exit 0``. Both halves pass on the base and on this branch.
    """
    job_id, exit_code = _finished_watched_job(
        workspace, "echo nothing-is-run-here")

    ev = _entry(workspace, job_id)
    assert ev is not None, "the watched job was never reported by job_monitor"
    state = str(ev.get("state") or "")
    assert state.startswith("FINISHED"), (
        "the state no longer leads with FINISHED, which is what every "
        "caller that shows it has always meant: " + repr(state))
    assert "exit 0" not in state, (
        "a clean exit should read as FINISHED, not as a code to explain: "
        + repr(state))
