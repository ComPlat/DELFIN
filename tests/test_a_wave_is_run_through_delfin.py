"""A wave is run, watched and steered through DELFIN, not through scripts
around it.

Until 2026-11-24 every wave was driven by scripts outside the repository
(a starter, a stage watcher, a wave watcher, a host guard, a friction
scan) that the orchestrator rewrote as each defect showed: a zombie that
counted as alive and stalled a stage for an hour, a starter whose lock
deadlocked the stage, real runs that imported DELFIN from their own
workspace. ``delfin/agent/swarm/`` now holds the verbs (``wave``), the
stage (``stages``), the real runs (``realrun``), the tick and the brake
(``guard``), the friction count (``friction``) and the MCP surface
(``mcp_server``). These tests run real git worktrees and real member
processes; the member is a stand-in for ``delfin-agent chat`` that
writes the same stream-json, so they need no model and no network.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent.swarm import friction, guard, mcp_server as M, realrun, stages, wave

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
    if word == "commit":
        open("work.txt", "a").write("x\n")
        subprocess.run(["git", "add", "work.txt"], check=True)
        subprocess.run(["git", "-c", "user.name=m", "-c", "user.email=m@example.org",
                        "commit", "-qm", "member work"], check=True)
    if word == "prmd":
        os.makedirs("scratch", exist_ok=True)
        open("scratch/PR.md", "w").write("# done\n")
    if word == "askrun":
        os.makedirs("scratch", exist_ok=True)
        json.dump({"probes": [{"name": "p1", "prompt": "show it"}]},
                  open("scratch/realrun.json", "w"))
print(json.dumps({"type": "permission_denied", "tool": "bash", "reason_code": "outside_ws",
                  "text": "/home/someone/x"}), flush=True)
print(json.dumps({"type": "notice", "text": "⏳ Transient API error (RateLimitError); retrying 1/7 in 2s"}))
print(json.dumps({"type": "result", "stop_reason": "done", "session_id": "sess-" + str(os.getpid()),
                  "tool_calls": [{"name": "bash"}, {"name": "read_file"}],
                  "text": "brief was: " + brief[:200], "resumed_from": session}))
'''


@pytest.fixture
def waves(tmp_path, monkeypatch):
    """A waves root, a git repository with a main branch, and members that
    are the stand-in instead of a model session."""
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
    # A slice of the test's own: the shipped name is the live swarm's, and an
    # emergency stop in a test stopped it, ending a running wave (2026-11-24).
    cfg = dict(wave._DEFAULTS, waves_root=str(root), min_available_gb=0, max_load_per_core=1e9,
               slice=f"delfin-test-{os.getpid()}.slice")
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    systemctl = []
    monkeypatch.setattr(guard, "_systemctl", lambda *a: systemctl.append(list(a))
                        or subprocess.CompletedProcess(a, 0, "", ""))
    monkeypatch.setattr(wave, "_scoped", lambda cmd, cfg: cmd)
    monkeypatch.setattr(wave, "_refuse_on_a_login_node", lambda: None)
    monkeypatch.setattr(wave, "member_argv", lambda run, task, session="", cfg=None, python=None:
                        [sys.executable, str(member), "-p", task["brief"],
                         "--cwd", str(wave.run_dir(run) / task["id"])]
                        + (["--session", session] if session else []))
    monkeypatch.setattr(guard, "_attention", lambda *a, **k: None)
    return {"root": root, "repo": repo, "systemctl": systemctl}


def _wait_ended(run, tid, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        a = wave.load_state(run)["agents"].get(tid)
        if a and not wave.alive(a["pid"]):
            return a
        time.sleep(0.05)
    raise AssertionError(f"{run}/{tid} did not end")


def _task(tid, brief, **kw):
    return {"id": tid, "brief": brief, "base": "main", **kw}


# -- the task list and the places ---------------------------------------------------

def test_a_bad_task_list_is_refused_with_every_problem_named():
    with pytest.raises(wave.WaveError) as exc:
        wave.validate_tasks([{"id": "a b", "brief": ""}, {"id": "x", "brief": "b", "after": ["nope"]},
                             {"id": "x", "brief": "b", "collect": ["r/t.final"], "turns": -1}])
    text = str(exc.value)
    for part in ("id must be", "brief is empty", "after must be", "id used twice",
                 "collect names a member", "turns must be"):
        assert part in text, (part, text)


def test_the_waves_root_is_never_under_the_hidden_home(monkeypatch):
    monkeypatch.setattr(wave, "settings", lambda: dict(wave._DEFAULTS, waves_root="~/.delfin/waves"))
    with pytest.raises(wave.WaveError, match="hides"):
        wave.root()


# -- start, dependencies, collect -------------------------------------------------------

def test_a_member_works_in_its_own_worktree_and_cannot_push(waves):
    out = wave.start("r1", [_task("a", "commit prmd")], waves["repo"], wait=False)
    assert out["started"] == ["a"]
    a = _wait_ended("r1", "a")
    wt = Path(a["worktree"])
    assert a["branch"] == "wave/r1/a"
    assert subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=wt,
                          capture_output=True, text=True).stdout.strip() == "wave/r1/a"
    push = subprocess.run(["git", "config", "--get", "remote.origin.pushurl"], cwd=wt,
                          capture_output=True, text=True).stdout.strip()
    assert push == wave._PUSH_NOWHERE
    row = wave.status("r1")["members"][0]
    assert row["state"] == "ended" and row["commits"] == 1 and row["pr_md"] is True
    assert row["stop_reason"] == "done"


def test_a_task_waits_for_what_it_depends_on_and_collects_its_work(waves):
    tasks = [_task("arm", "prmd sleep:0.5"),
             _task("int", "commit", after=["r2/arm"], collect=["r2/arm"])]
    first = wave.start("r2", tasks, waves["repo"], wait=False)
    assert first["started"] == ["arm"] and first["waiting"] == ["int"]
    _wait_ended("r2", "arm")
    second = wave.start("r2", tasks, waves["repo"], wait=False)
    assert second["started"] == ["int"]
    collected = wave.run_dir("r2") / "int" / "scratch" / "arms" / "arm" / "PR.md"
    assert collected.read_text() == "# done\n"


def test_a_quorum_starts_after_its_wait(waves, monkeypatch):
    task = _task("q", "x", after=["r3/a", "r3/b"], after_min=1, after_wait_s=100)
    ended = {"r3/a": True, "r3/b": False}
    monkeypatch.setattr(wave, "ended", lambda ref: ended[ref])
    monkeypatch.setattr(wave, "_ended_at", lambda ref: 1000.0)
    assert wave.ready(task, now=1050.0) is False
    assert wave.ready(task, now=1100.0) is True
    ended["r3/b"] = True
    assert wave.ready(task, now=1001.0) is True


def test_a_zombie_counts_as_ended(tmp_path):
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    deadline = time.time() + 10
    while time.time() < deadline:                  # exited, not yet reaped
        try:
            state = Path(f"/proc/{proc.pid}/stat").read_text().rsplit(")", 1)[1].split()[0]
        except OSError:
            break
        if state == "Z":
            break
        time.sleep(0.05)
    assert wave.alive(proc.pid) is False
    proc.wait()


def test_a_started_task_is_never_started_twice(waves):
    tasks = [_task("a", "sleep:0.3")]
    wave.start("r4", tasks, waves["repo"], wait=False)
    pid = wave.load_state("r4")["agents"]["a"]["pid"]
    again = wave.start("r4", tasks, waves["repo"], wait=False)
    assert again["started"] == [] and wave.load_state("r4")["agents"]["a"]["pid"] == pid


def test_a_stale_read_cannot_start_a_task_twice(waves, monkeypatch):
    """Two starters at once: each read the queue before the other launched.
    The launch re-reads the state under the lock and starts nothing the
    other already started."""
    tasks = [_task("a", "sleep:0.3")]
    wave.start("r4b", tasks, waves["repo"], wait=False)
    pid = wave.load_state("r4b")["agents"]["a"]["pid"]
    real = wave.load_state
    reads = {"n": 0}

    def stale_first(run):
        reads["n"] += 1
        st = real(run)
        if reads["n"] == 2:     # the queue's read, outside the lock (the 1st is under it)
            st = dict(st, agents={})
        return st

    monkeypatch.setattr(wave, "load_state", stale_first)
    again = wave.start("r4b", tasks, waves["repo"], wait=False)
    monkeypatch.setattr(wave, "load_state", real)
    assert again["started"] == [] and real("r4b")["agents"]["a"]["pid"] == pid


def test_the_harness_moves_to_main_and_back_when_main_does_not_import(waves, tmp_path, monkeypatch):
    h = tmp_path / "harness"
    subprocess.run(["git", "clone", "-q", str(waves["repo"]), str(h)], check=True)
    (h / "delfin").mkdir()
    (h / "delfin" / "__init__.py").write_text("")
    subprocess.run(["git", "-C", str(h), "add", "."], check=True)
    subprocess.run(["git", "-C", str(h), "-c", "user.name=t", "-c", "user.email=t@example.org",
                    "commit", "-qm", "pkg"], check=True)
    subprocess.run(["git", "-C", str(h), "update-ref", "refs/remotes/origin/main", "HEAD"], check=True)
    cfg = dict(wave.settings(), harness=str(h))
    monkeypatch.setattr(subprocess, "run", _no_fetch(subprocess.run))
    assert wave.refresh_harness(cfg)["moved"] is False            # already on main
    head = subprocess.check_output(["git", "-C", str(h), "rev-parse", "HEAD"], text=True).strip()
    (waves["repo"] / "x").write_text("x")                        # a new main that does not import
    subprocess.run(["git", "-C", str(waves["repo"]), "add", "x"], check=True)
    subprocess.run(["git", "-C", str(waves["repo"]), "commit", "-qm", "x"], check=True)
    subprocess.run(["git", "-C", str(h), "fetch", "-q", str(waves["repo"]), "+main:refs/remotes/origin/main"],
                   check=True)
    with pytest.raises(wave.WaveError, match="does not import"):
        wave.refresh_harness(cfg)
    back = subprocess.check_output(["git", "-C", str(h), "rev-parse", "HEAD"], text=True).strip()
    assert back == head


def _no_fetch(real):
    """git fetch from origin is the network: the test sets origin/main itself."""
    def run(cmd, *a, **k):
        if isinstance(cmd, list) and "fetch" in cmd and "origin" in cmd and cmd[-1] == "origin":
            return subprocess.CompletedProcess(cmd, 0, "", "")
        return real(cmd, *a, **k)
    return run


# -- events and watch --------------------------------------------------------------------

def test_watch_returns_when_a_member_ends(waves):
    wave.start("r5", [_task("a", "prmd sleep:0.3")], waves["repo"], wait=False)
    got = wave.watch(["r5"], timeout_s=30, poll_s=0.1)
    kinds = [e["kind"] for e in got["events"]]
    if kinds == ["started"]:
        got = wave.watch(["r5"], timeout_s=30, poll_s=0.1, cursors=got["cursors"])
        kinds = [e["kind"] for e in got["events"]]
    ended = [e for e in got["events"] if e["kind"] == "ended"]
    assert ended and ended[0]["stop_reason"] == "done" and ended[0]["pr_md"] is True, kinds
    quiet = wave.watch(["r5"], timeout_s=0.3, poll_s=0.1, cursors=got["cursors"])
    assert quiet["events"] == []                       # nothing twice


# -- resume and stop ----------------------------------------------------------------------

def test_a_resume_continues_the_session_with_the_note(waves):
    wave.start("r6", [_task("a", "prmd")], waves["repo"], wait=False)
    first = _wait_ended("r6", "a")
    sid = wave._last_result(Path(first["log"]))["session_id"]
    out = wave.resume("r6", "a", "fix the test", wait_for_room=False)
    assert out["resumed"] is True
    a = _wait_ended("r6", "a")
    assert a["resumes"] == 1
    res = wave._last_result(Path(a["log"]))
    assert res["resumed_from"] == sid
    assert "fix the test" in a["task"]["brief"] and a["task"]["orig_brief"] == "prmd"


def test_a_running_member_is_not_resumed_and_stop_ends_it(waves):
    wave.start("r7", [_task("a", "sleep:30")], waves["repo"], wait=False)
    assert wave.resume("r7", "a", "x", wait_for_room=False)["resumed"] is False
    assert wave.stop("r7")["stopped"] == ["a"]
    assert not wave.alive(wave.load_state("r7")["agents"]["a"]["pid"])


def test_a_later_stage_is_appended_to_a_running_wave(waves):
    M.h_start("r16", [_task("land", "prmd")], str(waves["repo"]))
    _wait_ended("r16", "land")
    out = wave.add_tasks("r16", [_task("rev", "prmd", after=["r16/land"], no_realrun=True)])
    assert out["added"] == ["rev"] and out["tasks"] == 2
    with pytest.raises(wave.WaveError, match="already has"):
        wave.add_tasks("r16", [_task("land", "x")])
    assert guard.tick()["runs"]["r16"].get("started") == ["rev"]


# -- the stage ------------------------------------------------------------------------

def test_the_stage_asks_once_runs_what_is_asked_and_marks_final(waves):
    tasks = [_task("rev", "prmd", no_realrun=True), _task("quiet", "prmd"),
             _task("asks", "askrun")]
    wave.start("r8", tasks, waves["repo"], wait=False)
    for t in tasks:
        _wait_ended("r8", t["id"])
    started = []
    did = stages.advance("r8", start_real_run=lambda run, tid, busy: started.append(tid))
    by = {d["id"]: d["did"] for d in did}
    assert by == {"rev": "final", "quiet": "asked for a real run", "asks": "real run"}
    assert started == ["asks"]
    d = wave.run_dir("r8")
    assert (d / "rev.final").exists() and (d / "quiet.asked").exists()
    assert (d / "asks.realrun_running").exists()       # no second tick starts it again
    _wait_ended("r8", "quiet")                         # its resume ended without asking
    again = {x["id"]: x.get("why") for x in stages.advance("r8", start_real_run=lambda *a: None)}
    assert again.get("quiet") == "ended twice without asking"


def test_a_member_ended_by_the_endpoint_is_resumed_not_finished(waves, monkeypatch):
    wave.start("r17", [_task("a", "prmd")], waves["repo"], wait=False)
    a = _wait_ended("r17", "a")
    (Path(a["worktree"]) / "scratch" / "PR.md").unlink()
    log = Path(a["log"])
    log.write_text(log.read_text() + json.dumps({"type": "result", "stop_reason": "error",
                                                 "error": "HTTP 503 temporarily unavailable"}) + "\n")
    resumed = []
    monkeypatch.setattr(wave, "resume", lambda run, tid, note="", **k: resumed.append(note))
    monkeypatch.setattr(stages, "ERROR_RESUME_AFTER_S", 0)
    for _ in range(stages.MAX_ERROR_RESUMES):
        assert stages.advance("r17")[0]["did"] == "resumed after an endpoint error"
    assert len(resumed) == 3 and "503" in resumed[0]
    after = stages.advance("r17", start_real_run=lambda *x: None)
    assert after[0]["did"] == "asked for a real run"      # the stage goes on as usual


def test_a_member_killed_from_outside_is_resumed_not_finished(waves, monkeypatch):
    """A session with no result line was killed (a stopped slice, an OOM
    kill), not finished: asking it for real runs and then marking it final
    turned ten killed integrators into ten empty finals (2026-11-24)."""
    import signal
    wave.start("r31", [_task("a", "sleep:60")], waves["repo"], wait=False)
    a = wave.load_state("r31")["agents"]["a"]
    # begun before measuring, as the brief asks: it does not make a kill an ending
    (Path(a["worktree"]) / "scratch").mkdir(exist_ok=True)
    (Path(a["worktree"]) / "scratch" / "PR.md").write_text("# PENDING\n")
    os.kill(int(a["pid"]), signal.SIGKILL)
    _wait_ended("r31", "a")
    assert wave._last_result(Path(a["log"])) == {}
    resumed = []
    monkeypatch.setattr(wave, "resume", lambda run, tid, note="", **k: resumed.append(note))
    monkeypatch.setattr(stages, "ERROR_RESUME_AFTER_S", 0)
    for _ in range(stages.MAX_ERROR_RESUMES):
        assert stages.advance("r31")[0]["did"] == "resumed after being ended from outside"
    assert "ended from outside" in resumed[0]
    assert not (wave.run_dir("r31") / "a.final").exists()
    after = stages.advance("r31", start_real_run=lambda *x: None)
    assert after[0]["did"] == "asked for a real run"      # then the stage goes on as usual


def test_a_real_run_request_is_checked():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "scratch").mkdir()
        for bad, why in (({"probes": []}, "1 to 3"),
                         ({"probes": [{"name": "Bad Name", "prompt": "x"}]}, "name"),
                         ({"probes": [{"name": "ok", "prompt": "x", "workspace": "/"}]}, "workspace")):
            (Path(d) / "scratch" / "realrun.json").write_text(json.dumps(bad))
            with pytest.raises(wave.WaveError, match=why):
                realrun.read_request(Path(d))


def test_the_real_runs_hand_each_provider_the_key_its_client_reads():
    import inspect
    from delfin.agent import api_client as A
    src = inspect.getsource(A.create_client)
    for provider, key in realrun.PROVIDER_KEYS.items():
        assert key in src or key in ("OPENAI_API_KEY", "ANTHROPIC_API_KEY") and key in inspect.getsource(A), (provider, key)
    assert 'key_env_var="KIT_TOOLBOX_API_KEY"' in src and 'key_env_var="BLABLADOR_API_KEY"' in src


# -- the tick and the brake -----------------------------------------------------------

def test_the_tick_starts_what_became_ready(waves):
    tasks = [_task("a", "prmd"), _task("b", "prmd", after=["r9/a"])]
    M.h_start("r9", tasks, str(waves["repo"]))
    _wait_ended("r9", "a")
    did = guard.tick()
    assert did["runs"]["r9"].get("started") == ["b"]


def test_a_stop_holds_until_a_person_clears_it(waves):
    tasks = [_task("a", "sleep:30"), _task("b", "prmd", after=["r10/a"])]
    M.h_start("r10", tasks, str(waves["repo"]))
    out = guard.emergency_stop("test")
    assert out["stopped"] == {"r10": ["a"]}
    assert "error" in json.loads(M.h_start("r11", [_task("c", "x")], str(waves["repo"])))
    assert guard.tick().get("held")                    # records, starts nothing
    assert "b" not in wave.load_state("r10")["agents"]
    assert guard.clear_stop() is True and guard.clear_stop() is False


def test_two_ticks_outside_the_envelope_stop_everything(waves, monkeypatch):
    M.h_start("r12", [_task("a", "sleep:30")], str(waves["repo"]))
    monkeypatch.setattr(guard, "outside_envelope", lambda: "load 999 on 1 cores")
    assert guard.tick()["stopped"] is None
    assert guard.tick()["stopped"]["stopped"] == {"r12": ["a"]}
    guard.clear_stop()


def test_a_member_listening_beyond_the_loopback_stops_everything(waves, monkeypatch):
    M.h_start("r13", [_task("a", "sleep:30")], str(waves["repo"]))
    monkeypatch.setattr(guard, "listeners_beyond_loopback", lambda pids: ["t 0.0.0.0:1025"])
    out = guard.tick()["stopped"]
    assert "0.0.0.0:1025" in out["why"] and out["stopped"] == {"r13": ["a"]}
    guard.clear_stop()


def test_the_tick_loop_outlives_a_tick_that_raises(waves, monkeypatch):
    calls, lines, slept = [], [], []
    def flaky():
        calls.append(1)
        if len(calls) == 2:
            raise RuntimeError("disk gone")
        return {"runs": {"r": {}, "q": {"started": ["a"], "events": [{"kind": "ended"}]}}}
    monkeypatch.setattr(guard, "tick", flaky)
    assert guard.loop(5, ticks=3, sleep=slept.append, log=lines.append) == 3
    assert len(calls) == 3 and slept == [5, 5]
    assert "tick failed: RuntimeError: disk gone" in lines[1]
    last = json.loads(lines[2])
    assert last["runs"] == 2 and last["moved"] == {"q": ["ended", "started"]}   # a quiet run is not listed


def test_an_emergency_stop_in_a_test_never_reaches_the_live_slice(waves):
    guard.emergency_stop("test")
    guard.clear_stop()
    assert waves["systemctl"] == [["stop", f"delfin-test-{os.getpid()}.slice"]]
    assert wave._DEFAULTS["slice"] not in json.dumps(waves["systemctl"])


def test_a_second_tick_while_one_runs_does_nothing(waves):
    with wave._flock(wave.root() / ".tick.lock") as held:
        assert held
        assert guard.tick() == {"runs": {}, "stopped": None, "busy": True}
    assert "busy" not in guard.tick()


def test_two_looks_at_once_record_each_event_once(waves):
    """The service's tick and an orchestrator's watch observe the same run;
    run side by side, an ended member is still one event."""
    import threading
    wave.start("r30", [_task("a", "write x")], str(waves["repo"]), wait=True)
    _wait_ended("r30", "a")
    got, errors = [], []
    def look():
        try:
            got.extend(wave.observe("r30"))
        except Exception as exc:
            errors.append(exc)
    threads = [threading.Thread(target=look) for _ in range(8)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    assert not errors
    assert [e["kind"] for e in wave.events("r30")["events"]].count("ended") == 1
    assert not list(wave.run_dir("r30").glob("*.tmp"))


def _fake_systemctl(monkeypatch, manager=True):
    calls = []
    def run(args):
        calls.append(list(args))
        ok = manager or args[0] != "show-environment"
        return subprocess.CompletedProcess(args, 0 if ok else 1, "active\n" if args[0] == "is-active" else "",
                                           "" if ok else "Failed to connect to bus")
    monkeypatch.setattr(guard, "_systemctl", lambda *a: run(a))
    return calls


def test_the_tick_service_runs_the_loop_from_the_harness(waves, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    cfg = dict(wave.settings(), harness=str(tmp_path / "har%ness"))
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    calls = _fake_systemctl(monkeypatch)
    out = guard.install_unit(60)
    unit = Path(out["unit"])
    assert unit == tmp_path / "cfg" / "systemd" / "user" / guard.UNIT and out["active"] == "active"
    text = unit.read_text()
    assert 'Environment="PYTHONPATH=' + str(tmp_path / "har%%ness") + '"' in text   # % escaped
    assert '"-P"' in text and '"swarm" "wave" "tick" "--every" "60"' in text
    assert "WorkingDirectory=" + str(waves["root"]) in text and "Restart=on-failure" in text
    assert ["enable", "--now", guard.UNIT] in calls
    # the command the unit runs is one the command line takes
    from delfin.agent import cli_swarm
    import argparse
    ns = cli_swarm.configure_parser(argparse.ArgumentParser()).parse_args(["wave", "tick", "--every", "60"])
    assert ns.wave_cmd == "tick" and ns.every == 60
    assert guard.remove_unit() is True and not unit.exists() and guard.remove_unit() is False


def test_a_host_without_systemd_is_told_how_to_tick_instead(waves, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    _fake_systemctl(monkeypatch, manager=False)
    with pytest.raises(wave.WaveError, match="wave tick --every 60"):
        guard.install_unit(60)
    assert not (tmp_path / "cfg").exists()
    assert '"error"' in M.h_guard_install(1, systemd=True)[:40]


def test_the_scheduler_runs_the_tick_as_a_pass_not_a_model_turn(waves, monkeypatch):
    from delfin.agent import scheduler_daemon as SD
    called = []
    monkeypatch.setattr(guard, "tick", lambda: called.append(1) or {"runs": {}})
    monkeypatch.setattr(SD, "run_entry", lambda *a, **k: pytest.fail("a model turn was started"))

    class Entry:
        id, kind, reason = "e1", "interval", ""
        prompt = f"{guard.PROMPT_MARKER} tick"
        workspace = str(waves["repo"])

    SD.make_fire_callback(log=lambda *_: None)(Entry())
    assert called == [1]
    assert SD.SWARM_GUARD_ENTRY_PREFIX == guard.PROMPT_MARKER


# -- friction and the MCP surface -------------------------------------------------------

def test_friction_is_counted_per_reason_and_per_model(waves):
    wave.start("r14", [_task("a", "prmd", model="m1")], waves["repo"], wait=False)
    _wait_ended("r14", "a")
    rep = friction.report(["r14"])
    reasons = {(r["kind"], r["code"]) for r in rep["by_reason"]}
    assert ("denial", "bash:outside_ws") in reasons and ("api_retry", "RateLimitError") in reasons
    m1 = next(m for m in rep["per_model"] if m["model"] == "m1")
    assert m1["tool_calls"] == 2 and m1["denial_per_100"] == 50.0
    assert all("/home/someone" not in r["example"] for r in rep["by_reason"])


def test_review_and_log_show_what_a_member_did(waves):
    wave.start("r15", [_task("a", "commit prmd")], waves["repo"], wait=False)
    _wait_ended("r15", "a")
    rev = json.loads(M.h_review("r15", "a"))
    assert rev["pr_md"] == "# done\n" and len(rev["commits"]) == 1
    assert rev["last_session"]["stop_reason"] == "done"
    log = json.loads(M.h_log("r15", "a"))
    assert any(s["type"] == "permission_denied" for s in log["steps"])
    assert "/home/someone" not in json.dumps(log)


def test_every_wave_verb_is_on_the_command_line():
    from delfin.agent import cli_swarm
    import argparse
    parser = cli_swarm.configure_parser(argparse.ArgumentParser(prog="swarm"))
    for argv in (["wave", "runs"], ["wave", "status", "r"], ["wave", "watch", "r", "--timeout", "1"],
                 ["wave", "review", "r", "a"], ["wave", "resume", "r", "a", "--note", "n"],
                 ["wave", "stop", "r"], ["wave", "tick"], ["wave", "emergency-stop", "why"],
                 ["wave", "clear-stop"], ["wave", "friction"], ["wave", "guard", "install"],
                 ["wave", "tick", "--every", "60"], ["wave", "guard", "install", "--systemd"],
                 ["wave", "mcp"], ["wave", "start", "r", "tasks.json"], ["wave", "add", "r", "more.json"],
                 ["wave", "land", "r1", "--push"], ["wave", "land", "--calibrate-pair", "a", "b"],
                 ["wave", "compose", "r", "spec.json", "--out", "t.json"]):
        ns = parser.parse_args(argv)
        assert ns.swarm_cmd == "wave" and ns.wave_cmd == argv[1]


def test_the_mcp_tools_name_every_handler():
    src = Path(M.__file__).read_text()
    tools = set(re.findall(r"@mcp\.tool\(\)\n    def (\w+)", src))
    assert tools == {"wave_compose", "wave_runs", "wave_status", "wave_start", "wave_add_tasks", "wave_watch", "wave_review",
                     "wave_log", "wave_resume", "wave_stop", "wave_tick", "wave_emergency_stop",
                     "wave_clear_stop", "wave_friction", "wave_guard_install", "wave_guard_remove",
                     "wave_land"}
    assert os.path.basename(M.__file__) == "mcp_server.py"


# -- the wave's shape -----------------------------------------------------------------

_M = {"glm": {"provider": "blablador", "model": "alias-glm-5.3-flash"},
      "ds": {"provider": "kit", "model": "kit.deepseek-v4-flash"}}


def test_a_wave_is_composed_in_the_shape_every_wave_had():
    from delfin.agent.swarm import compose
    pk = [{"id": "p1", "title": "one", "body": "do one"}, {"id": "p2", "title": "two", "body": "do two"}]
    tasks = compose.compose("w", pk, arms=[_M["glm"], _M["ds"]], integrator=_M["ds"],
                            lander=_M["ds"], reviewer=_M["glm"], fixer=_M["ds"])
    ids = [t["id"] for t in tasks]
    assert ids == ["p1-a1", "p1-a2", "p1-int", "p2-a1", "p2-a2", "p2-int", "zz-land",
                   "zz-rev-behaviour", "zz-rev-tests", "zz-rev-safety", "zz-rev-claims", "zz-fix"]
    by = {t["id"]: t for t in tasks}
    assert by["p1-a1"]["model"] != by["p1-a2"]["model"] and by["p1-a1"]["no_realrun"]
    assert by["p1-int"]["after"] == ["w/p1-a1", "w/p1-a2"] == by["p1-int"]["collect"]
    assert by["zz-land"]["after"] == ["w/p1-int.final", "w/p2-int.final"]
    assert by["zz-rev-safety"]["after"] == ["w/zz-land.final"] and by["zz-rev-safety"]["land"] == "none"
    assert by["zz-fix"]["after"] == [f"w/zz-rev-{k}" for k in compose.REVIEW_FOCUS]
    assert all(compose.DELIVERY_RULES.strip() in by[i]["brief"] for i in ("p1-a1", "p1-int", "zz-land"))
    real = compose.REAL_RUNS.strip()
    assert real in by["p1-int"]["brief"] and real not in by["p1-a1"]["brief"]


def test_reviewers_of_the_landers_family_are_refused():
    from delfin.agent.swarm import compose
    with pytest.raises(wave.WaveError, match="not independent"):
        compose.compose("w", [{"id": "p", "title": "t", "body": "b"}], arms=[_M["ds"]],
                        integrator=_M["ds"], lander=_M["ds"], reviewer=_M["ds"], fixer=_M["ds"])


def test_a_composed_wave_runs_through_the_tick(waves):
    """Arms first; the integrator once its arms ended; nothing of the later
    stages before their .final marks exist."""
    from delfin.agent.swarm import compose
    tasks = compose.compose("r20", [{"id": "p", "title": "t", "body": "prmd"}],
                            arms=[_M["glm"], _M["ds"]], integrator=_M["ds"], lander=_M["ds"],
                            reviewer=_M["glm"], fixer=_M["ds"], base="main")
    for t in tasks:
        t["brief"] = "prmd"
    M.h_start("r20", tasks, str(waves["repo"]))
    assert sorted(wave.load_state("r20")["agents"]) == ["p-a1", "p-a2"]
    _wait_ended("r20", "p-a1"), _wait_ended("r20", "p-a2")
    guard.tick()
    assert "p-int" in wave.load_state("r20")["agents"]
    assert "zz-land" not in wave.load_state("r20")["agents"]


def test_a_real_run_may_name_the_empty_workspace():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "scratch").mkdir()
        (Path(d) / "scratch" / "realrun.json").write_text(json.dumps(
            {"probes": [{"name": "p", "prompt": "x", "workspace": "empty"}]}))
        assert realrun.read_request(Path(d))["probes"][0]["name"] == "p"


def test_a_stopped_integration_stage_is_continued_with_the_arms_that_ran(waves):
    """Arms that already ran in another run are read and collected; no new
    arm is started."""
    from delfin.agent.swarm import compose
    wave.start("old", [_task("p-x", "prmd"), _task("p-y", "prmd")], waves["repo"], wait=False)
    _wait_ended("old", "p-x"), _wait_ended("old", "p-y")
    tasks = compose.compose("new", [{"id": "p", "title": "t", "body": "b"}], arms=[],
                            integrator=_M["ds"], lander=_M["ds"], reviewer=_M["glm"], fixer=_M["ds"],
                            stages=("integrate", "land"), arm_refs={"p": ["old/p-x", "old/p-y"]})
    by = {t["id"]: t for t in tasks}
    assert list(by) == ["p-int", "zz-land"]
    assert by["p-int"]["after"] == ["old/p-x", "old/p-y"] == by["p-int"]["collect"]
    assert "branch wave/old/p-x" in by["p-int"]["brief"]


def test_the_tick_drives_only_runs_started_through_the_interface(waves):
    """A run some other tool left in the root keeps its queue: the tick
    starts nothing in it."""
    wave.start("other", [_task("a", "prmd")], waves["repo"], wait=False)
    _wait_ended("other", "a")
    wave.add_tasks("other", [_task("b", "prmd")])
    assert guard.tick()["runs"]["other"].get("started") is None
    assert "b" not in wave.load_state("other")["agents"]
