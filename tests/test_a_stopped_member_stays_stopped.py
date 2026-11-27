"""A member a person stopped stays stopped.

A stop that the tick undoes is not a stop. The record: on 2026-11-26
``swarm wave stop w59 t1-a2`` wrote ``stopped``, yet the next tick resumed
t1-a2 under a new pid as if it had ended on an error, and it started the
package's integrator because the stopped arm counted as an ended one. The
wider evidence (2026-11-27) is the same defect after an emergency stop: once
the marker was cleared, the members the stop had ended were classed "silent"
and resumed, and a lander was started although its ``after`` named three
markers, two of which did not exist.

The stop is the person's decision, so the decision has to outlive the pid.
``wave.stop`` records it in the run's state, and both rules follow from that
one record:

  stopped  no path restarts a member a person stopped -- not the stage's
           error resume, not its "was ended from outside" resume, not the
           "ask for a real run" resume, and not ``wave resume`` by hand.
  ended    an ``after`` that names a stopped member is NOT met: the stop is
           not the member finishing. The gate holds and ``status`` says why.

Every case here drives the verbs a person and the service drive --
``cli_swarm.main(["wave", ...])``, ``guard.tick()`` and ``wave.status`` --
and judges the run by what it observes: the pids in ``state.json`` and the
markers on disk. Nothing here reads the live waves root or the real slice:
the ``waves`` fixture points the root at a tmp dir, as the wave suite does.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import cli_swarm
from delfin.agent.swarm import guard, stages, wave

_MEMBER = r'''
import json, os, subprocess, sys, time
args = sys.argv[1:]
brief = args[args.index("-p") + 1]
wt = args[args.index("--cwd") + 1]
session = args[args.index("--session") + 1] if "--session" in args else ""
os.chdir(wt)
for word in brief.split():
    if word.startswith("sleep:"):
        time.sleep(float(word[6:]))
print(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                  "session_id": "sess-" + str(os.getpid()),
                  "tool_calls": [], "text": brief[:100]}), flush=True)
'''


@pytest.fixture
def waves(tmp_path, monkeypatch):
    """A waves root, a git repository with a main branch, and members that
    are a stand-in for a model session -- no model and no network."""
    root = tmp_path / "waves"
    repo = tmp_path / "repo"
    repo.mkdir()
    for args in (["init", "-q", "-b", "main"], ["config", "user.email", "t@example.org"],
                 ["config", "user.name", "t"]):
        subprocess.run(["git", *args], cwd=repo, check=True)
    (repo / "README").write_text("r\n")
    (repo / ".gitignore").write_text("scratch/\n")
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "commit", "-qm", "init"], cwd=repo, check=True)
    member = tmp_path / "member.py"
    member.write_text(_MEMBER)
    # A slice and a root of the test's own: the shipped names belong to the
    # live swarm, and an emergency stop here must never reach it.
    cfg = dict(wave._DEFAULTS, waves_root=str(root), min_available_gb=0, max_load_per_core=1e9,
               slice=f"delfin-test-{os.getpid()}.slice")
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    monkeypatch.setattr(guard, "_systemctl", lambda *a: subprocess.CompletedProcess(a, 0, "", ""))
    monkeypatch.setattr(wave, "_scoped", lambda cmd, cfg: cmd)
    monkeypatch.setattr(wave, "_refuse_on_a_login_node", lambda: None)
    monkeypatch.setattr(wave, "member_argv", lambda run, task, session="", cfg=None, python=None:
                        [sys.executable, str(member), "-p", task["brief"],
                         "--cwd", str(wave.run_dir(run) / task["id"])]
                        + (["--session", session] if session else []))
    monkeypatch.setattr(guard, "_attention", lambda *a, **k: None)
    # The host seams the tick reads, so no case measures this machine.
    monkeypatch.setattr(guard, "outside_envelope", lambda: "")
    monkeypatch.setattr(guard, "_children", lambda pid: {int(pid)})
    monkeypatch.setattr(guard, "listeners_beyond_loopback", lambda pids: [])
    # The error resume is what a person's stop is up against: without this
    # the pause would hide the resume behind 15 minutes of silence.
    monkeypatch.setattr(stages, "ERROR_RESUME_AFTER_S", 0)
    return {"root": root, "repo": repo}


def _wait_ended(run, tid, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        a = wave.load_state(run)["agents"].get(tid)
        if a and not wave.alive(a["pid"]):
            return a
        time.sleep(0.05)
    raise AssertionError(f"{run}/{tid} did not end")


def _task(tid, brief="work", **kw):
    return {"id": tid, "brief": brief, "base": "main", **kw}


def _cli(*argv) -> dict:
    """Run one ``delfin-agent swarm wave VERB`` and read what it printed."""
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = cli_swarm.main(["wave", *argv])
    assert rc == 0, f"wave {argv[0]} exited {rc}"
    text = buf.getvalue().strip()
    # The verb prints its JSON answer; the whole of it is the last block.
    start = text.find("{")
    assert start >= 0, f"wave {argv[0]} printed no answer: {text!r}"
    return json.loads(text[start:])


def _start_like_the_service(run, tasks, repo) -> dict:
    """A run the tick drives, started through the verb the service uses."""
    path = Path(repo) / f"{run}-tasks.json"
    path.write_text(json.dumps(tasks))
    return _cli("start", run, str(path), "--repo", str(repo))


def _stopped_member(waves, run, tid="arm", brief="work", **kw):
    """A member that ran, ended, and then a person stopped it.

    A member that is already ended is exactly the case the record names: the
    stop must hold on the member's NAME, because there is no pid left to
    hold it on.
    """
    _start_like_the_service(run, [_task(tid, brief, **kw)], waves["repo"])
    _wait_ended(run, tid)
    before = wave.load_state(run)["agents"][tid]["pid"]
    said = _cli("stop", run, tid)
    return before, said


# -- 1. the stage does not resume a stopped member ------------------------------------

def test_a_stopped_member_is_not_resumed_by_the_tick(waves):
    """The record's first half: three ticks over a stopped member resume it
    zero times, and write no resume counter.

    Red on the base for the reason the case names: the tick saw the member as
    "ended from outside" (its log has a result line, but the stage's other
    branch reads the same silence) and called ``wave.resume`` once per tick up
    to ``MAX_ERROR_RESUMES``. The control below is the same run without the
    stop, where the resume is what the stage is for.
    """
    before, said = _stopped_member(waves, "s1", "arm")
    assert said["stopped"] == ["arm"], said
    d = wave.run_dir("s1")
    for _ in range(3):
        guard.tick()
    # A resumed member would show as a new pid in the state, or as a counter
    # the stage increments before it resumes.
    after = wave.load_state("s1")["agents"]["arm"]["pid"]
    assert after == before, f"the stopped member was resumed under pid {after}"
    assert not wave.alive(after)
    assert not (d / "arm.error_resumes").exists(), \
        "the stage counted an error resume for a member a person stopped"
    assert not (d / "arm.asked").exists(), "the stage asked a stopped member for real runs"


def test_without_the_stop_the_same_member_is_final(waves):
    """Control for the case above: the same run, no stop, and the stage does
    what it is for -- a member that ended without asking is asked once. This
    passes on the base and after the change; it says the assertion above is
    about the stop and not about a stage that never acts.
    """
    _start_like_the_service("s2", [_task("arm", "work")], waves["repo"])
    _wait_ended("s2", "arm")
    d = wave.run_dir("s2")
    guard.tick()
    assert (d / "arm.asked").exists(), "the stage did not ask a member that simply ended"


def test_a_stopped_reviewer_is_not_marked_final(waves):
    """The second half: a reviewer (``no_realrun``) that a person stopped is
    not marked final. ``.final`` is what a later stage's ``after`` waits for,
    so writing it turns the person's stop into a finished arm and lets the
    stage depend on it.

    Red on the base: ``stages.advance`` marked every stopped ``no_realrun``
    member final on the first tick.
    """
    _before, said = _stopped_member(waves, "s3", "rev", no_realrun=True)
    assert said["stopped"] == ["rev"], said
    d = wave.run_dir("s3")
    for _ in range(2):
        guard.tick()
    assert not (d / "rev.final").exists(), \
        "a stopped reviewer was marked final, so a stage may now wait on it as done"


# -- 2. an `after` that names a stopped member is not met -----------------------------

def test_a_stage_waits_for_a_stopped_arm_and_says_why(waves):
    """The evidence's second half: ``zz-land`` was started although its
    ``after`` named markers that did not exist, because the stopped arms
    counted as ended. Here the arm is stopped and the dependent task must not
    start -- and ``status`` must say it is waiting on a stop, not on a run.
    """
    _start_like_the_service("s4", [
        _task("arm", "work"),
        _task("land", "work", after=["s4/arm"]),
    ], waves["repo"])
    _wait_ended("s4", "arm")
    assert wave.ended("s4/arm"), "the arm must have ended before it is stopped"
    _cli("stop", "s4", "arm")
    d = wave.run_dir("s4")
    for _ in range(2):
        guard.tick()
    assert "land" not in wave.load_state("s4")["agents"], \
        "a stage whose `after` names a stopped arm was started anyway"
    assert not (d / "land.log").exists()
    rows = {m["id"]: m for m in wave.status("s4")["members"]}
    assert rows["land"]["state"] == "waiting", rows["land"]
    assert "stop" in json.dumps(rows["land"]).lower(), \
        f"status does not say the stage waits on a stop: {rows['land']}"


def test_a_stage_does_start_once_the_arm_really_ended(waves):
    """Control for the case above: the same task list where the arm is left
    to end on its own. The stage starts, so the gate above is about the stop
    and not about an ``after`` that never resolves.
    """
    _start_like_the_service("s5", [
        _task("arm", "work"),
        _task("land", "work", after=["s5/arm"]),
    ], waves["repo"])
    _wait_ended("s5", "arm")
    guard.tick()
    assert "land" in wave.load_state("s5")["agents"], "a ready stage was not started"


# -- 3. the stop record and the status word ------------------------------------------

def test_status_marks_a_stopped_member_as_stopped(waves):
    """A person reading ``wave status`` must see the stop, not ``ended``:
    the two are the same pid state and opposite meanings, and the record's
    defect was that nothing told them apart.
    """
    _before, said = _stopped_member(waves, "s6", "arm")
    assert said["stopped"] == ["arm"], said
    row = {m["id"]: m for m in wave.status("s6")["members"]}["arm"]
    assert row["state"] == "stopped", row
    # The marks a person needs beside it: the stop is why it is not running.
    assert "stopped" in row["marks"]


def test_a_member_that_ended_on_its_own_is_not_marked_stopped(waves):
    """Control for the status word: an ordinary ending is still ``ended``.
    Without this the case above would pass for a status that calls everything
    stopped.
    """
    _start_like_the_service("s7", [_task("arm", "work")], waves["repo"])
    _wait_ended("s7", "arm")
    row = {m["id"]: m for m in wave.status("s7")["members"]}["arm"]
    assert row["state"] == "ended", row
    assert "stopped" not in row["marks"]


# -- 4. the stop is a person's decision, also by hand --------------------------------

def test_resume_by_hand_does_not_restart_a_stopped_member(waves):
    """``wave resume`` is a verb a person runs; it must not undo a stop by
    another person's command. Red on the base: the CLI resumed the stopped
    member under a new pid.
    """
    before, _said = _stopped_member(waves, "s8", "arm")
    out = _cli("resume", "s8", "arm", "--note", "carry on")
    assert out.get("resumed") is False, out
    assert out.get("why") == "stopped by a person", out
    after = wave.load_state("s8")["agents"]["arm"]["pid"]
    assert after == before, f"resume started the stopped member as pid {after}"


def test_resume_after_the_stop_is_lifted_starts_the_member(waves):
    """Control for the case above: once a person starts it again (``wave
    start`` for a task the run does not have yet is the shipped way), the
    resume works. Here the lift is the record's own: a member with no stop
    record resumes as it always did.
    """
    _start_like_the_service("s9", [_task("arm", "work")], waves["repo"])
    _wait_ended("s9", "arm")
    out = _cli("resume", "s9", "arm", "--note", "carry on")
    assert out.get("resumed") is True, out
    a = wave.load_state("s9")["agents"]["arm"]
    assert wave.alive(a["pid"])


def test_a_person_lifts_a_stop_with_after_stop_and_then_the_stage_goes_on(waves):
    """The stop is a decision, not a tombstone. A person who stopped a member
    by mistake must be able to start it again -- ``wave resume --after-stop``
    -- and after that the member is an ordinary member again: the stop record
    is older than the launch, so the stage treats it as running and, when it
    ends, as ended. Without this the rule above would be a member nothing in
    DELFIN could ever start again.
    """
    before, _said = _stopped_member(waves, "s10", "arm")
    out = _cli("resume", "s10", "arm", "--note", "carry on", "--after-stop")
    assert out.get("resumed") is True, out
    a = wave.load_state("s10")["agents"]["arm"]
    assert a["pid"] != before and wave.alive(a["pid"])
    assert not wave.stopped_by_a_person("s10", "arm"), \
        "the stop record still outranks the newer launch, so the member is still held"


def test_a_failed_real_run_verdict_does_not_resume_a_stopped_member(waves):
    """``verdict_fail`` is the fourth restarter and does not go through the
    stage at all: it calls ``wave.resume`` from ``merge.py``. A stop must hold
    there too, so the assertion is on the verb a person runs, and it names the
    reason rather than a silent no-op.
    """
    from delfin.agent.swarm import merge
    before, said = _stopped_member(waves, "s11", "arm")
    assert said["stopped"] == ["arm"], said
    out = merge.record_verdict("s11", "arm", "fail", by="integrator",
                               detail="the probe measured nothing")
    assert out.get("resumed") is False, out
    after = wave.load_state("s11")["agents"]["arm"]["pid"]
    assert after == before, f"the verdict started the stopped member as pid {after}"

