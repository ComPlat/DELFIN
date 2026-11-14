"""``delfin-voila status`` and ``delfin-voila stop``.

A server that no terminal holds has to be found and ended from any
terminal. ``status`` says whether it runs, where, for how long, which
dashboards are kept (each with the address back into it and its last
activity) and how many agent sessions are kept. ``stop`` ends it the way
Ctrl+C did: SIGTERM to the launcher, which takes the Ctrl+C path
(cli_voila._terminal_gone) -- the server shuts its kernels down, and the
agent sessions with Keep session on are handed to the agent service as on
any end. It asks nothing and says what it ended.

The launcher is named by pid AND start time: a record whose pid was
handed to another process ends nothing.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone

import pytest

from delfin import cli_voila
from delfin.agent import proc_identity
from delfin.dashboard import kept_agent_sessions as KAS
from delfin.dashboard import server_record as SR
from delfin.dashboard import session as S
from conftest import child_env

RESUME = "/voila/render/delfin_voila_runtime/delfin_resume.ipynb"


@pytest.fixture(autouse=True)
def _private_stores(monkeypatch, tmp_path):
    monkeypatch.setenv(KAS.DIR_ENV, str(tmp_path / "kept_agent_sessions"))
    monkeypatch.setattr(S, "_server_stdout", lambda: None)


def _record(pid=None, proc_start=None, **extra):
    rec = SR.own_record(pid or os.getpid())
    if proc_start is not None:
        rec["proc_start"] = proc_start
    rec.update({"port": 8899, "ip": "127.0.0.1", "scheme": "http",
                "url": "/voila/render/delfin_voila_runtime/delfin_dashboard.ipynb",
                "resume_path": RESUME, "token": "k" * 43,
                "started_at": time.time() - 3725, "detached": True,
                "log": str(SR.log_file())})
    rec.update(extra)
    return rec


def test_status_says_when_nothing_runs(capsys):
    assert cli_voila.main(["status"]) == 1
    assert "No dashboard server runs on this host" in capsys.readouterr().out


def test_status_lists_the_server_and_what_it_keeps(monkeypatch, capsys):
    SR.write(_record())
    S.write_record("node-ab12", kid="kid-kept")
    KAS.arm("sess-1", "/work", pid=os.getpid())
    three_min = (datetime.now(timezone.utc) - timedelta(minutes=3)).isoformat()
    monkeypatch.setattr(SR, "kernels", lambda record, timeout=5.0: [
        {"id": "kid-kept", "connections": 0, "last_activity": three_min},
        {"id": "kid-open", "connections": 1, "last_activity": three_min},
    ])
    assert cli_voila.main(["status"]) == 0
    out = capsys.readouterr().out
    assert "running" in out and "port 8899" in out
    assert "up 1 h 2 min" in out
    assert "node-ab12" in out
    assert f"http://localhost:8899{RESUME}?session=node-ab12" in out
    assert "3 min ago" in out and "no window" in out
    assert "2 dashboards open" in out
    assert "1 kept agent session" in out
    assert str(SR.log_file()) in out
    assert "k" * 43 not in out, "the token is not written into a non-terminal"
    assert "delfin-voila stop" in out


_STAND_IN = r"""
# a stand-in for the delfin cli_voila launcher: SIGTERM is its Ctrl+C
import pathlib, signal, sys, time
marker, ready = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
def bye(*_):
    marker.write_text("TERM")
    sys.exit(0)
signal.signal(signal.SIGTERM, bye)
ready.write_text("ready")
while True:
    time.sleep(0.05)
"""


@pytest.fixture
def stand_in(tmp_path):
    marker, ready = tmp_path / "marker", tmp_path / "ready"
    proc = subprocess.Popen([sys.executable, "-c", _STAND_IN, str(marker),
                             str(ready)], env=child_env(tmp_path))
    deadline = time.monotonic() + 20
    while not ready.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    assert ready.exists()
    yield proc, marker
    if proc.poll() is None:
        proc.kill()
    proc.wait()


def test_stop_ends_the_launcher_the_way_ctrl_c_does(stand_in, monkeypatch,
                                                     capsys):
    proc, marker = stand_in
    SR.write(_record(pid=proc.pid))
    S.write_record("node-cd34", kid="kid-1")
    monkeypatch.setattr(SR, "kernels", lambda record, timeout=5.0: [
        {"id": "kid-1", "connections": 0, "last_activity": ""},
        {"id": "kid-2", "connections": 1, "last_activity": ""}])
    assert cli_voila.main(["stop"]) == 0
    assert marker.read_text() == "TERM", "SIGTERM: the launcher's Ctrl+C path"
    out = capsys.readouterr().out
    assert "Stopped the dashboard server" in out
    assert f"pid {proc.pid}" in out and "port 8899" in out
    assert "2 dashboards" in out and "node-cd34" in out
    assert proc_identity.alive(proc.pid) is False


def test_stop_does_not_touch_a_pid_that_is_another_process_now(
        stand_in, capsys):
    proc, marker = stand_in
    SR.write(_record(pid=proc.pid, proc_start="1"))
    assert cli_voila.main(["stop"]) == 0
    assert not marker.exists(), "a recycled pid is not the launcher"
    assert proc.poll() is None
    assert "No dashboard server runs on this host" in capsys.readouterr().out


def test_stop_says_which_kept_agent_sessions_went_to_the_service(
        stand_in, monkeypatch, capsys):
    proc, marker = stand_in
    SR.write(_record(pid=proc.pid))
    monkeypatch.setattr(SR, "kernels", lambda record, timeout=5.0: [])
    # Held by the stand-in's process: when it ends, the hand-over resumes
    # it in the agent service and removes the record (here: by hand).
    KAS.arm("sess-handed", "/work", pid=proc.pid, title="Analysis")
    monkeypatch.setattr(cli_voila, "_HANDOVER_WAIT_S", 3.0)

    def _hand_over_once():
        for rec in KAS.records():
            os.unlink(rec["_path"])

    monkeypatch.setattr(cli_voila, "_await_hand_over",
                        lambda ids, timeout: (_hand_over_once(), (ids, []))[1])
    assert cli_voila.main(["stop"]) == 0
    out = capsys.readouterr().out
    assert "1 kept agent session" in out
    assert "agent service" in out
