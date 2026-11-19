"""A one-shot run waits for its session's background job and wakes on it.

Input: ``delfin-agent run`` whose turn started a background shell (or
watches a SLURM job) and ended -- ``watch_job`` tells the model "you will be
notified in a future turn, end your turn".

Before: the run ended with that answer, its shells were stopped at exit and
nothing woke anyone. Seen on 2026-11-14: an agent in a batch job started a
calculation chain with bash_background, ended its turn saying it would be
woken -- and the chain died seconds later with the run.
Now the run stays up while a shell of its session runs or a watched job is
pending, and runs a wake turn (the terminal's system-event message) when one
finished -- as a background command in Claude Code reports back on its own.
Bounded by ``--job-wait`` seconds (0: the old behaviour) and a wake count.
"""
from __future__ import annotations

import io

import pytest

from delfin.agent import cli, job_monitor, job_wake


class _Engine:
    session_id = "s-1"


class _Clock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t

    def sleep(self, s):
        self.t += s


@pytest.fixture
def world(monkeypatch):
    """running: shell ids still running; finished: events the next look sees."""
    w = {"running": [], "watched": False, "finished": [], "turns": []}
    monkeypatch.setattr(job_wake, "wake_enabled", lambda *a, **k: True)
    monkeypatch.setattr(cli, "_own_running_shells", lambda sid, registry=None: list(w["running"]))
    monkeypatch.setattr(job_monitor, "session_is_watched", lambda ws, sid: w["watched"])

    def shells(seen, *, registry=None, session_id=None):
        out, w["finished"] = list(w["finished"]), []
        if out:
            w["running"], w["watched"] = [], False
        return out

    monkeypatch.setattr(job_wake, "finished_shells", shells)
    monkeypatch.setattr(job_wake, "finished_watched_jobs", lambda ws, seen, **k: [])

    def run_once(engine, prompt, **kw):
        w["turns"].append(prompt)
        return {"text": "the chain finished", "error": None, "stop_reason": "done"}

    monkeypatch.setattr(cli, "_run_once", run_once)
    import delfin.agent.run_limits as rl
    monkeypatch.setattr(rl, "ceiling_reason", lambda engine: "")
    return w


def _wait(tmp_path, wait_s, clock, err):
    first = {"text": "Started it; I'll be told when it is done.", "error": None,
             "stop_reason": "done"}
    return cli._wait_for_watched_jobs(_Engine(), tmp_path, first, wait_s=wait_s,
                                      max_tokens=0, sleep=clock.sleep, clock=clock,
                                      err=err)


def test_a_running_shell_keeps_the_run_up_and_its_end_wakes_the_agent(tmp_path, world):
    world["running"] = ["fd1d30dd"]
    clock, err = _Clock(), io.StringIO()
    world["finished"] = [{"kind": "shell", "job_id": "fd1d30dd", "state": "exit 0",
                          "description": "chain1"}]
    out = _wait(tmp_path, 3600, clock, err)
    assert len(world["turns"]) == 1
    assert "fd1d30dd" in world["turns"][0]
    assert "not the user" in world["turns"][0]
    assert out["text"].startswith("Started it")
    assert out["text"].endswith("the chain finished")
    assert "waiting for it" in err.getvalue()


def test_a_watched_slurm_job_is_waited_for_too(tmp_path, world, monkeypatch):
    world["watched"] = True
    events = [[], [{"kind": "slurm", "job_id": "4711", "state": "COMPLETED"}]]

    def watched(ws, seen, **k):
        ev = events.pop(0) if events else []
        if ev:
            world["watched"] = False
        return ev

    monkeypatch.setattr(job_wake, "finished_watched_jobs", watched)
    _wait(tmp_path, 3600, _Clock(), io.StringIO())
    assert len(world["turns"]) == 1 and "4711" in world["turns"][0]


def test_nothing_running_ends_the_run_at_once(tmp_path, world):
    clock = _Clock()
    out = _wait(tmp_path, 3600, clock, io.StringIO())
    assert world["turns"] == [] and clock.t == 0.0
    assert out["text"].startswith("Started it")


def test_job_wait_zero_is_the_old_behaviour(tmp_path, world):
    world["running"] = ["fd1d30dd"]
    clock = _Clock()
    _wait(tmp_path, 0, clock, io.StringIO())
    assert world["turns"] == [] and clock.t == 0.0


def test_the_wait_is_bounded_and_says_so(tmp_path, world):
    world["running"] = ["forever"]
    clock, err = _Clock(), io.StringIO()
    _wait(tmp_path, 120, clock, err)
    assert world["turns"] == []
    assert 120 <= clock.t <= 120 + cli._JOB_WAIT_POLL_S
    assert "the run ends now" in err.getvalue()


def test_the_run_command_takes_job_wait():
    parser = cli.build_parser()
    assert parser.parse_args(["run", "hi"]).job_wait == cli._JOB_WAIT_DEFAULT_S
    assert parser.parse_args(["run", "--job-wait", "0", "hi"]).job_wait == 0


def test_own_running_shells_are_this_sessions_only():
    class J:
        def __init__(self, jid, sid, code):
            self.job_id, self.session_id, self._c = jid, sid, code

        def poll(self):
            return self._c

    class R:
        def list_jobs(self, include_finished=True):
            return [J("a", "s-1", None), J("b", "s-2", None), J("c", "s-1", 0)]

    assert cli._own_running_shells("s-1", registry=R()) == ["a"]
