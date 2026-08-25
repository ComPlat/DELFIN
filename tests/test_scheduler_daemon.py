"""Tests for delfin.agent.scheduler_daemon (headless schedule executor).

No real LLM calls, no network: the engine factory is always mocked and
the session store is redirected into tmp_path.
"""

from __future__ import annotations

import argparse
import os
import signal
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import job_monitor as JM
from delfin.agent import scheduler as S
from delfin.agent import scheduler_daemon as D
from delfin.agent import session_store as SS


@pytest.fixture
def cron_path(tmp_path):
    return tmp_path / "cron.json"


@pytest.fixture
def pid_path(tmp_path):
    return tmp_path / "scheduler_daemon.pid"


@pytest.fixture
def sessions_dir(tmp_path, monkeypatch):
    d = tmp_path / "sessions"
    monkeypatch.setattr(SS, "_SESSIONS_DIR", d)
    return d


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


class FakeEngine:
    """Just enough surface for cli._run_once / cli._save_session."""

    def __init__(self):
        self.prompts: list[str] = []
        self.token_usage = {"input": 3, "output": 5}
        self.session_id = ""
        self.outcomes: list[tuple] = []

    def stream_response(self, *, user_message, on_token=None,
                        on_tool_use=None, max_tokens=0):
        self.prompts.append(user_message)
        if on_token:
            on_token("done")
        return "done"

    def export_state(self):
        return {"mode": "solo", "engine_messages": [],
                "token_usage": dict(self.token_usage), "cost_usd": 0.0}

    def record_cycle_outcome(self, verdict, user_task, **_kw):
        self.outcomes.append((verdict, user_task))


def _backdate(sch, ent):
    ent.next_fire_at = time.time() - 1
    sch._save()


# ---------------------------------------------------------------------------
# PID single-instance guard
# ---------------------------------------------------------------------------

def test_pid_guard_blocks_second_instance(pid_path):
    assert D.acquire_pid_lock(pid_path) is True
    # Our own (live) pid is in the file — a second instance must refuse.
    assert D.acquire_pid_lock(pid_path) is False
    D.release_pid_lock(pid_path)
    assert not pid_path.exists()
    assert D.acquire_pid_lock(pid_path) is True
    D.release_pid_lock(pid_path)


def test_stale_pid_lock_is_taken_over(pid_path, monkeypatch):
    pid_path.write_text("99999")
    monkeypatch.setattr(JM, "_pid_alive", lambda _pid: False)
    assert D.acquire_pid_lock(pid_path) is True
    assert pid_path.read_text().strip() == str(os.getpid())
    D.release_pid_lock(pid_path)


def test_daemon_status_reports_running_and_entries(
        pid_path, cron_path, workspace):
    st = D.daemon_status(path=pid_path, cron_path=cron_path)
    assert st == {"running": False, "pid": 0, "entries": 0, "disabled": 0}

    sch = S.Scheduler(path=cron_path)
    sch.schedule_once(delay_seconds=60, prompt="x", workspace=str(workspace))
    assert D.acquire_pid_lock(pid_path) is True
    try:
        st = D.daemon_status(path=pid_path, cron_path=cron_path)
        assert st["running"] is True
        assert st["pid"] == os.getpid()
        assert st["entries"] == 1
        assert st["disabled"] == 0
    finally:
        D.release_pid_lock(pid_path)


# ---------------------------------------------------------------------------
# Firing semantics (engine factory mocked)
# ---------------------------------------------------------------------------

def test_once_entry_fires_exactly_once_and_is_consumed(
        cron_path, workspace, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="do the thing", reason="nightly",
        workspace=str(workspace))
    _backdate(sch, ent)

    engines: list[tuple[str, FakeEngine]] = []

    def factory(ws, settings=None):
        e = FakeEngine()
        engines.append((ws, e))
        return e

    fire = D.make_fire_callback(engine_factory=factory, log=lambda _s: None)
    assert sch.tick(fire_callback=fire) == 1

    # Engine built for the entry's workspace, prompt ran once.
    assert len(engines) == 1
    ws_arg, engine = engines[0]
    assert ws_arg == str(workspace)
    assert len(engine.prompts) == 1
    assert "do the thing" in engine.prompts[0]
    assert f"[scheduled:{ent.id}]" in engine.prompts[0]
    assert "nightly" in engine.prompts[0]
    # Cycle outcome recorded as PASS.
    assert engine.outcomes and engine.outcomes[0][0] == "PASS"
    # Session was saved.
    assert list(sessions_dir.glob("*.json"))
    # "once" entries are consumed — a second tick must not re-fire.
    assert sch.list_entries() == []
    assert sch.tick(fire_callback=fire) == 0
    assert len(engines) == 1


def test_interval_entry_is_rescheduled(cron_path, workspace, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_interval(
        every_seconds=600, prompt="poll", workspace=str(workspace))
    _backdate(sch, ent)

    fire = D.make_fire_callback(
        engine_factory=lambda _ws, _s=None: FakeEngine(),
        log=lambda _s: None)
    assert sch.tick(fire_callback=fire) == 1
    entries = sch.list_entries()
    assert len(entries) == 1
    assert entries[0].fire_count == 1
    assert entries[0].fail_count == 0
    assert entries[0].next_fire_at > time.time() + 500
    # Not due any more — no second fire.
    assert sch.tick(fire_callback=fire) == 0


def test_failing_entry_does_not_kill_loop(cron_path, tmp_path, sessions_dir):
    ws_bad = tmp_path / "bad"
    ws_bad.mkdir()
    ws_ok = tmp_path / "ok"
    ws_ok.mkdir()

    sch = S.Scheduler(path=cron_path)
    bad = sch.schedule_once(
        delay_seconds=60, prompt="explode", workspace=str(ws_bad))
    good = sch.schedule_once(
        delay_seconds=60, prompt="fine", workspace=str(ws_ok))
    _backdate(sch, bad)
    _backdate(sch, good)

    def factory(ws, settings=None):
        if ws == str(ws_bad):
            raise RuntimeError("engine boom")
        return FakeEngine()

    fire = D.make_fire_callback(engine_factory=factory, log=lambda _s: None)
    # The raising entry is isolated; the good one still fires.
    assert sch.tick(fire_callback=fire) == 1
    remaining = {e.id: e for e in sch.list_entries()}
    assert good.id not in remaining          # consumed
    assert remaining[bad.id].fail_count == 1
    assert remaining[bad.id].disabled is False
    # Failure count is persisted.
    reloaded = {e.id: e for e in S.Scheduler(path=cron_path).list_entries()}
    assert reloaded[bad.id].fail_count == 1


def test_three_consecutive_failures_disable_entry(
        cron_path, workspace, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_interval(
        every_seconds=600, prompt="poll", workspace=str(workspace))
    _backdate(sch, ent)

    calls = []

    def factory(_ws, _settings=None):
        calls.append(1)
        raise RuntimeError("always broken")

    fire = D.make_fire_callback(engine_factory=factory, log=lambda _s: None)
    for _ in range(3):
        assert sch.tick(fire_callback=fire) == 0
    assert len(calls) == 3
    ent2 = sch.list_entries()[0]
    assert ent2.disabled is True
    assert "consecutive failures" in ent2.disabled_reason
    # Persisted, and never fired again.
    reloaded = S.Scheduler(path=cron_path).list_entries()[0]
    assert reloaded.disabled is True
    assert sch.tick(fire_callback=fire) == 0
    assert len(calls) == 3


def test_missing_workspace_skips_and_disables(cron_path, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="x",
        workspace="/nonexistent/delfin/workspace/zz9")
    _backdate(sch, ent)

    calls = []
    fire = D.make_fire_callback(
        engine_factory=lambda *_a, **_k: calls.append(1) or FakeEngine(),
        log=lambda _s: None)
    assert sch.tick(fire_callback=fire) == 0
    assert calls == []                       # no engine, zero tokens
    ent2 = sch.list_entries()[0]
    assert ent2.disabled is True
    assert "no longer exists" in ent2.disabled_reason
    # Persisted flag.
    assert S.Scheduler(path=cron_path).list_entries()[0].disabled is True


def test_entry_without_recorded_workspace_is_disabled(
        cron_path, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="legacy")
    ent.workspace = ""                       # simulate a pre-daemon entry
    _backdate(sch, ent)

    fire = D.make_fire_callback(
        engine_factory=lambda *_a, **_k: FakeEngine(), log=lambda _s: None)
    assert sch.tick(fire_callback=fire) == 0
    ent2 = sch.list_entries()[0]
    assert ent2.disabled is True
    assert "no workspace recorded" in ent2.disabled_reason


def test_budget_threads_into_engine_cost_cap(
        cron_path, workspace, sessions_dir):
    class CapEngine(FakeEngine):
        def _cost_hard_cap(self):
            return 50.0

    engine = CapEngine()
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="capped", workspace=str(workspace),
        budget_usd=1.5)
    result = D.run_entry(ent, engine_factory=lambda *_a, **_k: engine)
    assert engine._cost_hard_cap() == 1.5    # per-entry budget wins
    assert engine.prompts and "capped" in engine.prompts[0]
    assert result["session_id"]


def test_run_entry_raises_on_turn_error(cron_path, workspace, sessions_dir):
    class ErrorEngine(FakeEngine):
        def stream_response(self, **_kw):
            raise RuntimeError("api down")

    engine = ErrorEngine()
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="x", workspace=str(workspace))
    with pytest.raises(RuntimeError, match="api down"):
        D.run_entry(ent, engine_factory=lambda *_a, **_k: engine)
    # Outcome still recorded as FAIL before the raise.
    assert engine.outcomes and engine.outcomes[0][0] == "FAIL"


# ---------------------------------------------------------------------------
# Loop + signals
# ---------------------------------------------------------------------------

def test_run_loop_executes_due_entry_from_disk(
        cron_path, workspace, sessions_dir):
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="from disk", workspace=str(workspace))
    _backdate(sch, ent)

    engines: list[FakeEngine] = []

    def factory(_ws, _settings=None):
        e = FakeEngine()
        engines.append(e)
        return e

    sleeps: list[float] = []
    rc = D.run_loop(
        cron_path=cron_path,
        engine_factory=factory,
        max_iterations=2,
        poll_s=0.05,
        sleep_fn=sleeps.append,
        log=lambda _s: None,
    )
    assert rc == 0
    assert len(engines) == 1                 # fired exactly once
    assert sleeps == [0.05]                  # slept only between passes
    # The consumed entry is gone on disk too.
    assert S.Scheduler(path=cron_path).list_entries() == []


def test_sigterm_handler_sets_stop_and_loop_exits(cron_path):
    stop = threading.Event()
    old_term = signal.getsignal(signal.SIGTERM)
    old_int = signal.getsignal(signal.SIGINT)
    try:
        handler = D.install_signal_handlers(stop)
        assert not stop.is_set()
        handler(signal.SIGTERM, None)        # invoke directly (no real kill)
        assert stop.is_set()
    finally:
        signal.signal(signal.SIGTERM, old_term)
        signal.signal(signal.SIGINT, old_int)

    # With the stop flag set the loop exits immediately: no tick, no engine.
    calls = []
    rc = D.run_loop(
        cron_path=cron_path,
        engine_factory=lambda *_a, **_k: calls.append(1),
        max_iterations=10,
        sleep_fn=lambda _s: pytest.fail("must not sleep after stop"),
        stop_event=stop,
        log=lambda _s: None,
    )
    assert rc == 0
    assert calls == []


# ---------------------------------------------------------------------------
# CLI surface: delfin-agent scheduler start|status|stop
# ---------------------------------------------------------------------------

def test_cli_scheduler_subcommand_registered():
    from delfin.agent import cli as C
    p = C.build_parser()
    for argv, action in (
        (["scheduler"], "status"),           # default action
        (["scheduler", "status"], "status"),
        (["scheduler", "start"], "start"),
        (["scheduler", "stop"], "stop"),
    ):
        args = p.parse_args(argv)
        assert args.func is C.cmd_scheduler
        assert args.scheduler_action == action


def test_cli_scheduler_start_spawns_detached_daemon(
        tmp_path, monkeypatch, capsys, pid_path, cron_path):
    from delfin.agent import cli as C
    monkeypatch.setattr(D, "_PID_PATH", pid_path)
    monkeypatch.setattr(S, "_DEFAULT_PATH", cron_path)
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    spawned = {}

    class FakePopen:
        def __init__(self, cmd, **kw):
            spawned["cmd"] = cmd
            spawned["kw"] = kw

    monkeypatch.setattr("subprocess.Popen", FakePopen)
    rc = C.cmd_scheduler(argparse.Namespace(scheduler_action="start"))
    assert rc == 0
    assert spawned["cmd"][1:] == ["-m", "delfin.agent.scheduler_daemon"]
    assert spawned["kw"]["start_new_session"] is True
    assert "started" in capsys.readouterr().out


def test_cli_scheduler_stop_sends_sigterm(
        monkeypatch, capsys, pid_path, cron_path):
    from delfin.agent import cli as C
    monkeypatch.setattr(D, "_PID_PATH", pid_path)
    monkeypatch.setattr(S, "_DEFAULT_PATH", cron_path)
    assert D.acquire_pid_lock(pid_path) is True

    killed: list[int] = []

    def fake_kill(pid, sig):
        if sig == signal.SIGTERM:
            killed.append(pid)
        # sig 0 (liveness probe in _pid_alive) succeeds silently

    monkeypatch.setattr("os.kill", fake_kill)
    try:
        rc = C.cmd_scheduler(argparse.Namespace(scheduler_action="stop"))
        assert rc == 0
        assert killed == [os.getpid()]
        assert "SIGTERM" in capsys.readouterr().out
    finally:
        D.release_pid_lock(pid_path)   # unlink only — no os.kill involved


def test_cli_scheduler_status_lists_entries(
        monkeypatch, capsys, pid_path, cron_path, workspace):
    from delfin.agent import cli as C
    monkeypatch.setattr(D, "_PID_PATH", pid_path)
    monkeypatch.setattr(S, "_DEFAULT_PATH", cron_path)
    sch = S.Scheduler(path=cron_path)
    ent = sch.schedule_once(
        delay_seconds=60, prompt="hello world", workspace=str(workspace))
    rc = C.cmd_scheduler(argparse.Namespace(scheduler_action="status"))
    out = capsys.readouterr().out
    assert rc == 0
    assert "not running" in out
    assert ent.id in out
    assert "hello world" in out
