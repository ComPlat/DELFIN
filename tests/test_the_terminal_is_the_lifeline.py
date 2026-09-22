"""Everything an agent starts ends with the terminal that started DELFIN.

Ctrl+C in the delfin-voila terminal -- or that terminal closing -- is the
one stop for everything: server, kernels, kept sessions, the agent's
background shells, MCP servers and the daemons. Before this, a server ran
five days after its launcher was gone, and the daemons were started to
survive the dashboard on purpose.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from delfin.agent import lifeline as L


@pytest.fixture
def root(tmp_path, monkeypatch):
    """This test process as the lifeline, with its own ledger."""
    monkeypatch.setattr(L, "_DIR", tmp_path / "lifeline")
    monkeypatch.setenv(L.ENV_PID, "0")
    monkeypatch.setenv(L.ENV_TICKS, "")
    L.claim_root()
    return tmp_path


def _detached(cmd: str = "sleep 60") -> subprocess.Popen:
    return subprocess.Popen(["bash", "-c", cmd], start_new_session=True)


def _ended(proc, seconds: float = 10.0) -> bool:
    try:
        proc.wait(timeout=seconds)
        return True
    except subprocess.TimeoutExpired:
        return False


def test_the_root_names_itself_by_pid_and_start_time(root):
    assert L.current() == (os.getpid(), L._start_ticks(os.getpid()))


def test_the_root_stopping_ends_every_detached_child(root):
    child = _detached()
    try:
        L.record_child(child.pid, "shell")
        assert L.end_children() == [child.pid]
        assert _ended(child)
    finally:
        child.kill()


def test_a_pid_that_now_belongs_to_another_process_is_left_alone(root):
    other = _detached()
    try:
        ledger = L._ledger(L.current())
        ledger.parent.mkdir(parents=True, exist_ok=True)
        ledger.write_text(json.dumps({
            "pid": other.pid, "ticks": (L._start_ticks(other.pid) or 0) + 1,
        }) + "\n", encoding="utf-8")
        assert L.end_children() == []
        assert other.poll() is None
    finally:
        other.kill()


def test_a_watcher_notices_that_its_lifeline_is_gone(tmp_path, monkeypatch,
                                                    real_lifeline_watch):
    launcher = subprocess.Popen([sys.executable, "-c",
                                 "import time; time.sleep(0.5)"])
    monkeypatch.setenv(L.ENV_PID, str(launcher.pid))
    monkeypatch.setenv(L.ENV_TICKS, str(L._start_ticks(launcher.pid) or ""))
    gone = threading.Event()
    L.watch(gone.set, poll_s=0.05)
    launcher.wait()
    assert gone.wait(10)


def test_a_daemon_started_from_a_terminal_follows_that_terminal(monkeypatch):
    monkeypatch.delenv(L.ENV_PID, raising=False)
    monkeypatch.delenv(L.ENV_TICKS, raising=False)
    assert L.child_env({})[L.ENV_PID] == str(os.getsid(0))


def test_an_agent_shell_is_on_the_ledger(root):
    from delfin.agent import bash_jobs as bj
    job = bj.get_registry().start(command="sleep 60", cwd=str(root),
                                  workspace=str(root))
    try:
        assert L.end_children() == [job.proc.pid]
        deadline = time.time() + 10
        while job.poll() is None and time.time() < deadline:
            time.sleep(0.1)
        assert job.poll() is not None
    finally:
        bj.get_registry().kill(job.job_id)


def test_the_daemons_end_with_their_lifeline(monkeypatch):
    from delfin.agent import bug_watcher, job_monitor, scheduler_daemon
    guarded = []
    monkeypatch.setattr(L, "guard_daemon", lambda: guarded.append(1))
    monkeypatch.setattr(job_monitor, "monitor_settings", lambda *a, **k: {
        "enabled": False, "interval_s": 60, "auto_diagnose": False})
    monkeypatch.setattr(scheduler_daemon, "acquire_pid_lock", lambda: False)
    monkeypatch.setattr(bug_watcher, "watcher_settings", lambda *a, **k: {
        "enabled": False, "interval_s": 60, "auto_analyze": False,
        "propose_fix": False})

    assert job_monitor.main() == 2
    assert scheduler_daemon.main() == 3
    assert bug_watcher.main() == 2
    assert guarded == [1, 1, 1], "every daemon watches its lifeline first"


def test_the_launcher_ends_what_the_dashboard_started(root, capsys):
    from delfin import cli_voila
    child = _detached()
    try:
        L.record_child(child.pid, "mcp")
        cli_voila._end_agent_processes()
        assert _ended(child)
        assert "Stopped 1 background process" in capsys.readouterr().out
    finally:
        child.kill()


def test_a_closed_terminal_stops_the_launcher_like_ctrl_c():
    from delfin import cli_voila
    before = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGHUP)}
    try:
        cli_voila._install_stop_signals()
        assert signal.getsignal(signal.SIGHUP) is cli_voila._terminal_gone
        assert signal.getsignal(signal.SIGTERM) is cli_voila._terminal_gone
        with pytest.raises(KeyboardInterrupt):
            cli_voila._terminal_gone(signal.SIGHUP, None)
    finally:
        for sig, handler in before.items():
            signal.signal(sig, handler)


def test_a_dashboard_kernel_watches_its_lifeline(monkeypatch):
    import atexit

    from delfin.dashboard import tab_agent
    watched = []
    monkeypatch.setattr(tab_agent, "_PROCESS_CLEANUP_REGISTERED", False)
    monkeypatch.setattr(atexit, "register", lambda *a, **kw: None)
    monkeypatch.setattr(L, "watch", lambda fn, **kw: watched.append(fn))
    tab_agent._register_process_exit_cleanup()
    assert watched == [tab_agent._kernel_lost_its_lifeline]
