"""A SLURM calculation the session itself submitted wakes the terminal idle prompt.

Wave-10 finding (s15/s17): two terminal sessions sat idle for 90 minutes
although their SLURM jobs had long finished — only an operator message
woke them. The dashboard tick drains the per-workspace agent watch file
(tab_agent._job_wake_tick -> check_agent_jobs, consume=False,
marker="wake_notified"); the terminal's only wake pull (_wake_text)
reports session messages and finished shells and nothing else.

The producer for the missing half already exists and is session-scoped:
every SLURM submission path registers the job with the submitting
session's id (api_client.py:8752 sbatch-stdout scan, api_client.py:11809
watch tool, bash_jobs.py:1023 background shell) into
<workspace>/.delfin/agent_watched_jobs.json, and
job_monitor.check_agent_jobs(session_id=...) reports exactly that
session's watches. Only the terminal pull is missing.

All fakes, no scheduler: the fake run_fn answers for every job id without
squeue; no real time, nothing waits.
"""

from __future__ import annotations

import pytest

from delfin.agent import job_monitor as JM
from delfin.agent import job_wake
from delfin.agent import repl as R


@pytest.fixture()
def agent(tmp_path, monkeypatch):
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    a.opts = type("O", (), {"cwd": str(tmp_path)})()
    a.engine = type("E", (), {"session_id": "sess-me"})()
    monkeypatch.setattr(job_wake, "wake_enabled", lambda *a, **k: True)
    return a


def _fake_run(states):
    """A run_fn that answers squeue/sacct line shapes from a dict.

    ``squeue -h -o %i %T`` and ``sacct -n -X -o JobID,State`` both print
    whitespace-separated ``JOBID STATE`` lines — exactly what
    ``_parse_state_lines`` reads.
    """

    def run(cmd):
        for jid, st in states.items():
            if jid in cmd:
                return f"{jid} {st}"
        return ""

    return run


def _register(tmp_path, job_id, session_id="sess-me", **extra):
    JM.register_agent_job(tmp_path, job_id, description="calc my molecule",
                          extra={"session_id": session_id, **extra})


# -- the producer: finished_watched_jobs ------------------------------------

def test_a_finished_slurm_job_wakes_the_terminal(agent, tmp_path):
    """The reported case: sbatch finished, the terminal sat idle."""
    run = _fake_run({"12345": "COMPLETED"})
    _register(tmp_path, "12345")
    seen: set = set()
    events = job_wake.finished_watched_jobs(
        str(tmp_path), seen, session_id="sess-me", run_fn=run)
    assert events, "a finished own SLURM job produced no wake event"
    ev = events[0]
    assert ev["kind"] == "slurm"
    assert ev["job_id"] == "12345"
    assert ev["ok"] is True
    assert ev["state"] == "COMPLETED"


def test_the_wake_prompt_names_the_job(agent, tmp_path):
    run = _fake_run({"12345": "COMPLETED"})
    _register(tmp_path, "12345")
    seen: set = set()
    events = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                            session_id="sess-me", run_fn=run)
    text = job_wake.wake_prompt(events)
    assert "12345" in text
    assert "not the user" in text.lower()
    assert "COMPLETED" in text


def test_another_sessions_job_does_not_wake_me(agent, tmp_path):
    """The 2026-09-28 scoping rule, now for the terminal too."""
    run = _fake_run({"99999": "COMPLETED"})
    _register(tmp_path, "99999", session_id="sess-other")
    seen: set = set()
    events = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                            session_id="sess-me", run_fn=run)
    assert events == []


def test_a_running_job_wakes_nobody(agent, tmp_path):
    run = _fake_run({"12345": "RUNNING"})
    _register(tmp_path, "12345")
    seen: set = set()
    events = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                            session_id="sess-me", run_fn=run)
    assert events == []


def test_each_job_wakes_exactly_once(agent, tmp_path):
    run = _fake_run({"12345": "COMPLETED"})
    _register(tmp_path, "12345")
    seen: set = set()
    first = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                           session_id="sess-me", run_fn=run)
    second = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                            session_id="sess-me", run_fn=run)
    assert first and second == [], "the same job woke the prompt twice"


def test_a_failed_job_is_a_failed_event(agent, tmp_path):
    run = _fake_run({"12345": "FAILED"})
    _register(tmp_path, "12345", folder=str(tmp_path))
    seen: set = set()
    events = job_wake.finished_watched_jobs(str(tmp_path), seen,
                                            session_id="sess-me", run_fn=run)
    assert events and events[0]["ok"] is False


# -- the wiring: repl._wake_text pulls the watch file too -------------------

def test_wired_into_the_terminal_wake_text(agent, tmp_path, monkeypatch):
    """The end-to-end path: idle prompt, empty chunk, own SLURM job done."""
    run = _fake_run({"12345": "COMPLETED"})
    _register(tmp_path, "12345")
    real = job_wake.finished_watched_jobs

    def _wired(workspace, seen, session_id="", run_fn=None):
        return real(workspace, seen, session_id=session_id, run_fn=run)

    monkeypatch.setattr(job_wake, "finished_watched_jobs", _wired)
    agent._wake_last_look = 0.0
    text = agent._wake_text("")
    assert "12345" in text, "the terminal wake never looked at the watch file"
    assert "not the user" in text.lower()


def test_the_terminal_look_is_not_a_poll(agent, monkeypatch):
    """One look per throttle window even when the box ticks fast."""
    calls = {"n": 0}

    def _count(workspace, seen, session_id="", run_fn=None):
        calls["n"] += 1
        return []

    monkeypatch.setattr(job_wake, "finished_watched_jobs", _count)
    for _ in range(50):
        agent._wake_text("")
    assert calls["n"] == 1, f"asked {calls['n']} times in one interval"


def test_the_terminal_wake_never_raises(agent, monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("watch file unreadable")
    monkeypatch.setattr(job_wake, "finished_watched_jobs", _boom)
    assert agent._wake_text("") == ""
