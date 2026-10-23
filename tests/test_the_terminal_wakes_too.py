"""Nothing in the terminal woke an agent. The dashboard did.

A background job that finishes DURING a turn is reported at the end of
it. One that finishes afterwards sat there: the drain is a pull, the only
thing that pulls is a turn, and the only thing that starts a turn is the
user typing. Reported from the field on 2026-09-18 — a session started a
run with ``bash_background`` + ``watch_job``, ended its turn, and was
never woken; the user had to nudge it by hand. The dashboard has had a
timer for this all along (``tab_agent._job_wake_tick``); the engine even
says so in a comment: "nothing in the codebase wakes an idle agent".

The idle prompt ticks now — the framed box reads keys on a timeout — so
the look costs nothing extra. What it must not do is take a half-written
message away, or fire for a user who turned it off.

The words are the dashboard's own, from delfin.agent.job_wake, which both
surfaces read. A second copy of that pair is what produced six live
wake-ups reading "shell None [?]".
"""

from __future__ import annotations

import pytest

from delfin.agent import job_wake
from delfin.agent import repl as R


class _Job:
    def __init__(self, job_id, code, command="python -m pytest -q"):
        self.job_id = job_id
        self.command = command
        self._code = code

    def poll(self):
        return self._code


class _Registry:
    def __init__(self, jobs):
        self._jobs = jobs

    def list_jobs(self, include_finished=False):
        return list(self._jobs)


@pytest.fixture()
def agent(monkeypatch):
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    monkeypatch.setattr(job_wake, "wake_enabled", lambda *a, **k: True)
    return a


#: The real producer, captured before anything patches the name — a
#: lambda that reaches for job_wake.finished_shells would call the patch.
_REAL_FINISHED = job_wake.finished_shells


def _finished(monkeypatch, jobs):
    monkeypatch.setattr(
        job_wake, "finished_shells",
        lambda seen, **k: _REAL_FINISHED(seen, registry=_Registry(jobs)))


# -- the reported case ------------------------------------------------------

def test_a_job_that_finished_while_idle_wakes_the_prompt(agent, monkeypatch):
    _finished(monkeypatch, [_Job("01e5b151", 0)])
    text = agent._wake_text("")
    assert "01e5b151" in text
    assert "pytest" in text
    assert "not the user" in text.lower()


def test_a_failed_job_says_how_it_failed(agent, monkeypatch):
    _finished(monkeypatch, [_Job("deadbeef", -9)])
    assert "killed by SIGKILL" in agent._wake_text("")


def test_nothing_finished_wakes_nobody(agent, monkeypatch):
    _finished(monkeypatch, [_Job("running1", None)])
    assert agent._wake_text("") == ""


# -- the two guards ---------------------------------------------------------

def test_a_half_written_message_is_never_taken_away(agent, monkeypatch):
    _finished(monkeypatch, [_Job("01e5b151", 0)])
    assert agent._wake_text("what I was typin") == "", (
        "a wake-up replaced the user's unsent message")


def test_a_user_who_turned_it_off_is_not_woken(agent, monkeypatch):
    _finished(monkeypatch, [_Job("01e5b151", 0)])
    monkeypatch.setattr(job_wake, "wake_enabled", lambda *a, **k: False)
    assert agent._wake_text("") == ""


def test_each_job_wakes_once(agent, monkeypatch):
    _finished(monkeypatch, [_Job("once1", 0)])
    assert agent._wake_text("")
    agent._wake_last_look = 0.0          # allow the next look
    assert agent._wake_text("") == "", "the same job woke the prompt twice"


def test_the_look_is_not_a_poll(agent, monkeypatch):
    """The read loop ticks ten times a second; the registry is asked far
    less often than that."""
    looks = {"n": 0}

    def _count(seen, **k):
        looks["n"] += 1
        return []
    monkeypatch.setattr(job_wake, "finished_shells", _count)
    for _ in range(50):
        agent._wake_text("")
    assert looks["n"] == 1, f"asked {looks['n']} times in one interval"


def test_it_never_raises(agent, monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("registry is gone")
    monkeypatch.setattr(job_wake, "finished_shells", _boom)
    assert agent._wake_text("") == ""


# -- one pair of words, not two ---------------------------------------------

def test_the_dashboard_and_the_terminal_say_the_same_thing():
    from delfin.dashboard import tab_agent as T
    ev = [{"kind": "shell", "job_id": "01e5b151", "state": "ok",
           "description": "python -m pytest -q"}]
    assert T._job_wake_prompt(ev) == job_wake.wake_prompt(ev)
