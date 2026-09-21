"""An emergency stop ends every agent of the user, on every login node.

On 2026-09-16 three agent sessions were kept on one login node; the next
login landed on another, and nothing there could stop them -- a lifeline
reaches only its own machine. What the login nodes share is the home
directory, so the stop is a file there, and every DELFIN process that
started before it ends itself. Afterwards nothing starts on its own until
somebody sends a session something by hand.

Every test here points the stop at its own directory: written into the real
home it would end the user's dashboards.
"""

from __future__ import annotations

import ast
import inspect
import os
import pathlib
import subprocess
import sys
import threading
import time
import uuid

import pytest

from delfin.agent import lifeline as L
from delfin.agent import scheduler as sched_mod
from delfin.agent import stop_all as S

_REPO = pathlib.Path(__file__).resolve().parents[1]
_REAL_STOP = pathlib.Path.home() / ".delfin" / "stop_all.json"

_WATCHER = r"""
import os, sys, time
from pathlib import Path
from delfin.agent import lifeline, stop_all
stop_all._PATH = Path(sys.argv[1])
lifeline.watch(lambda: os._exit(7), poll_s=0.05)
print("ready", flush=True)
time.sleep(60)
"""


@pytest.fixture
def stop_file(tmp_path, monkeypatch):
    path = tmp_path / "home" / ".delfin" / "stop_all.json"
    monkeypatch.setattr(S, "_PATH", path)
    monkeypatch.setattr(sched_mod, "_DEFAULT_PATH",
                        tmp_path / "home" / ".delfin" / "cron.json")
    assert S._PATH != _REAL_STOP and _REAL_STOP.parent not in path.parents
    return path


def _child_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items()
           if k not in (L.ENV_PID, L.ENV_TICKS)}
    env["PYTHONPATH"] = str(_REPO)
    env.update(extra)
    return env


def _watching_child(stop_path) -> subprocess.Popen:
    child = subprocess.Popen(
        [sys.executable, "-c", _WATCHER, str(stop_path)], cwd=str(_REPO),
        env=_child_env(), stdout=subprocess.PIPE, text=True)
    assert child.stdout.readline().strip() == "ready"
    return child


def _sleeper(*argv: str, **env) -> subprocess.Popen:
    return subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(60)", *argv],
        env=_child_env(**env))


def test_without_a_stop_nothing_is_stopped(stop_file):
    assert S.last_stop() == {}
    assert not S.stopped_since_start()
    assert S.wakes_allowed(0.0)


def test_a_stop_ends_a_process_that_started_before_it(stop_file):
    child = _watching_child(stop_file)
    try:
        S.request("test", end_here=False)
        assert child.wait(timeout=15) == 7
    finally:
        child.kill()


def test_a_process_started_after_the_stop_carries_on(stop_file):
    S.request("test", end_here=False)
    # A process's start is known to about a second (the boot time is kept
    # in whole seconds); one started well after the stop is not older.
    time.sleep(1.5)
    child = _watching_child(stop_file)
    try:
        with pytest.raises(subprocess.TimeoutExpired):
            child.wait(timeout=1.0)
    finally:
        child.kill()
        child.wait()


def test_a_process_under_no_lifeline_still_obeys_the_stop(stop_file, monkeypatch,
                                                         real_lifeline_watch):
    monkeypatch.delenv(L.ENV_PID, raising=False)
    monkeypatch.delenv(L.ENV_TICKS, raising=False)
    monkeypatch.setattr(S, "_STARTED_AT", time.time() - 60)
    gone = threading.Event()
    assert L.watch(gone.set, poll_s=0.05) is not None
    S.request("test", end_here=False)
    assert gone.wait(10)


def test_the_stop_disables_every_schedule(stop_file, tmp_path):
    sch = sched_mod.Scheduler()
    first = sch.schedule_once(delay_seconds=600, prompt="look again",
                              workspace=str(tmp_path))
    second = sch.schedule_once(delay_seconds=900, prompt="and again",
                               workspace=str(tmp_path))
    summary = S.request("test", end_here=False)
    assert summary["schedules_disabled"] == 2
    entries = {e.id: e for e in sched_mod.Scheduler().list_entries()}
    for entry in (first, second):
        assert entries[entry.id].disabled
        assert "emergency stop" in entries[entry.id].disabled_reason


def test_after_a_stop_a_session_waits_for_somebody_to_send(stop_file):
    before = time.time()
    S.request("test", end_here=False)
    assert not S.wakes_allowed(0.0), "a session nobody has written to"
    assert not S.wakes_allowed(before), "a message sent before the stop"
    assert S.wakes_allowed(time.time() + 0.01), "a message sent after it"
    note = S.held_note()
    assert "emergency stop" in note and S._hostname() in note


def test_the_stop_ends_what_a_dashboard_started_on_this_machine(stop_file, monkeypatch):
    mark = f"DELFIN_TEST_DASHBOARD_{uuid.uuid4().hex}"
    monkeypatch.setattr(S, "_DASHBOARD_ENV", (mark + "=",))
    started_by_a_dashboard = _sleeper(**{mark: "8866"})
    the_server_itself = _sleeper("jupyter-server", **{mark: "8866"})
    unrelated = _sleeper()
    try:
        summary = S.request("test", settle_s=0.2, grace_s=5)
        assert summary["found_here"] == [started_by_a_dashboard.pid]
        assert started_by_a_dashboard.wait(timeout=10) is not None
        assert the_server_itself.poll() is None
        assert unrelated.poll() is None
    finally:
        for proc in (started_by_a_dashboard, the_server_itself, unrelated):
            proc.kill()
            proc.wait()


def test_a_daemon_or_a_terminal_agent_is_found_by_its_command_line(stop_file, monkeypatch):
    mark = f"delfin-test-daemon-{uuid.uuid4().hex}"
    monkeypatch.setattr(S, "_COMMAND_MARKERS", (mark,))
    daemon = _sleeper(mark)
    try:
        summary = S.request("test", settle_s=0.2, grace_s=5)
        assert summary["found_here"] == [daemon.pid]
        assert daemon.wait(timeout=10) is not None
    finally:
        daemon.kill()
        daemon.wait()


def test_a_process_that_ends_on_its_own_is_not_signalled(stop_file, monkeypatch):
    mark = f"delfin-test-daemon-{uuid.uuid4().hex}"
    monkeypatch.setattr(S, "_COMMAND_MARKERS", (mark,))
    quick = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(0.5)", mark],
        env=_child_env())
    try:
        summary = S.request("test", settle_s=10, grace_s=5)
        assert summary["found_here"] == [quick.pid]
        assert summary["ended_here"] == []
    finally:
        quick.kill()
        quick.wait()


def test_the_command_asks_first_and_stops_nothing_on_no(stop_file, monkeypatch, capsys):
    from delfin.agent import cli
    called = []
    monkeypatch.setattr(S, "request", lambda *a, **k: called.append(1))
    monkeypatch.setattr("builtins.input", lambda prompt="": "n")
    assert cli.main(["stop-all"]) == 1
    assert called == []
    assert "Nothing stopped" in capsys.readouterr().out


def test_the_command_says_what_it_did(stop_file, monkeypatch, capsys):
    from delfin.agent import cli
    monkeypatch.setattr(S, "request", lambda reason="", **k: {
        "at": time.time(), "host": "uc3n990", "pid": 1, "reason": reason,
        "count": 1, "schedules_disabled": 2, "found_here": [11, 12],
        "ended_here": [12]})
    assert cli.main(["stop-all", "--yes"]) == 0
    out = capsys.readouterr().out
    assert "Emergency stop given" in out
    assert "2 agent process(es) found; 1 ended on their own, 1 had to be ended" in out
    assert "Schedules disabled: 2" in out
    assert "Nothing starts on its own again" in out


# -- the dashboard ------------------------------------------------------------

def _tab_agent_source() -> ast.Module:
    from delfin.dashboard import tab_agent
    return ast.parse(pathlib.Path(inspect.getfile(tab_agent)).read_text(
        encoding="utf-8"))


def _fn(tree: ast.Module, name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found")


def _calls(fn: ast.FunctionDef) -> set[str]:
    return {node.func.id for node in ast.walk(fn)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}


@pytest.mark.parametrize("starter", [
    "_deliver_session_messages",    # a message from another session
    "_on_wake",                     # a scheduled wake-up
    "_job_wake_tick",               # a finished job or background agent
])
def test_a_turn_nobody_typed_goes_through_the_stop_check(starter):
    calls = _calls(_fn(_tab_agent_source(), starter))
    assert "_send_on_its_own" in calls
    assert "_on_send" not in calls


def test_a_message_sent_by_hand_is_what_re_arms_a_session():
    tree = _tab_agent_source()
    door = ast.unparse(_fn(tree, "_send_on_its_own"))
    # ast.unparse writes every string with single quotes
    assert "wakes_allowed" in door and "'_on_its_own'" in door
    send = ast.unparse(_fn(tree, "_on_send"))
    assert "state['_armed_at'] = time.time()" in send
    assert "if not state.get('_on_its_own')" in send


def test_a_kernel_ending_with_its_lifeline_lets_go_of_its_kept_session(monkeypatch):
    from delfin.dashboard import session as kept
    from delfin.dashboard import tab_agent
    dropped = []
    monkeypatch.setattr(tab_agent, "_stop_what_this_process_started", lambda: None)
    monkeypatch.setattr(kept, "drop_record", lambda *a, **k: dropped.append(k))
    monkeypatch.setattr(os, "_exit", lambda code: (_ for _ in ()).throw(SystemExit(code)))
    with pytest.raises(SystemExit):
        tab_agent._kernel_lost_its_lifeline()
    assert len(dropped) == 1


def test_the_stop_file_is_redirected_for_every_test():
    from delfin.agent.state_paths import USER_STATE_SINKS
    assert ("delfin.agent.stop_all", "_PATH", "stop_all.json") in USER_STATE_SINKS
    assert S._PATH != _REAL_STOP


# -- the check ----------------------------------------------------------------

def test_the_check_changes_nothing_and_says_when_nothing_runs(stop_file, capsys):
    from delfin.agent import cli
    assert cli.main(["stop-all", "--check"]) == 0
    assert not stop_file.exists(), "a check is not a stop"
    assert "no agent is running or scheduled" in capsys.readouterr().out


def test_the_check_names_a_schedule_that_would_start_an_agent(stop_file, tmp_path, capsys):
    from delfin.agent import cli
    sched_mod.Scheduler().schedule_once(delay_seconds=600, prompt="look again",
                                        workspace=str(tmp_path))
    assert cli.main(["stop-all", "--check"]) == 1
    out = capsys.readouterr().out
    assert "Active schedules:        1" in out and "look again" in out


def test_the_check_finds_a_process_whose_dashboard_is_gone(stop_file):
    gone = subprocess.Popen([sys.executable, "-c", "pass"])
    gone_ticks = L._start_ticks(gone.pid)
    gone.wait()
    orphan = _sleeper(**{L.ENV_PID: str(gone.pid), L.ENV_TICKS: str(gone_ticks)})
    try:
        deadline = time.time() + 5
        rows = []
        while time.time() < deadline and not rows:
            rows = [r for r in S._outlived_their_start() if r["pid"] == orphan.pid]
            time.sleep(0.1)
        assert rows and "lifeline" in rows[0]["why"]
    finally:
        orphan.kill()
        orphan.wait()


def test_the_shell_a_stop_is_typed_into_is_not_one_of_the_agents():
    assert os.getppid() in S._ancestors(os.getpid())


def _report(**overrides):
    report = {"host": "uc3n990", "at": time.time(), "last_stop": {},
              "open_sessions": [], "kept_sessions": [], "active_schedules": [],
              "daemon_pid_files": {}, "agent_processes_here": [],
              "outlived_their_start_here": []}
    report.update(overrides)
    return report


def test_a_session_silent_for_minutes_is_not_reported_as_running(capsys):
    from delfin.agent import cli
    code = cli._print_stop_all_check(_report(open_sessions=[
        {"host": "uc3n991", "title": "Session C", "seconds_since_heartbeat": 600}]))
    out = capsys.readouterr().out
    assert code == 0
    assert "ended" in out and "left over" in out


def test_a_session_with_a_fresh_heartbeat_is_running(capsys):
    from delfin.agent import cli
    code = cli._print_stop_all_check(_report(open_sessions=[
        {"host": "uc3n991", "title": "Session C", "seconds_since_heartbeat": 30}]))
    assert code == 1
    assert "something is running" in capsys.readouterr().out


def test_left_over_records_are_named_but_not_counted_as_running(capsys):
    from delfin.agent import cli
    code = cli._print_stop_all_check(_report(
        kept_sessions=[{"session": "uc3n991-e2d3", "host": "uc3n991",
                        "since": time.time(), "here": False, "alive_here": False}],
        daemon_pid_files={"scheduler_daemon": {"pid": 4242, "alive_here": False}}))
    out = capsys.readouterr().out
    assert code == 0
    assert "run this check on uc3n991" in out
    assert "left over unless it runs on another login node" in out


def test_a_daemon_running_here_is_running(capsys):
    from delfin.agent import cli
    code = cli._print_stop_all_check(_report(
        daemon_pid_files={"job_monitor": {"pid": 4242, "alive_here": True}}))
    assert code == 1
