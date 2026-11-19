"""One dashboard server per user and host, found again by its record.

A second ``delfin-voila`` used to start a second server on the next free
port. Now the launcher leaves a record (server_record) once its server
answers, and the next launch opens that one.

Pinned here:

  the record          owner-only, in an owner-only directory, with the
                      launcher's pid AND start time
  believed only       while that pid runs with that start time: a dead
                      pid or a recycled one is a stale record, removed
  a second launch     finds the running server and starts nothing; with
                      --strict-port on another port it refuses (an SSH
                      tunnel points at the port it asked for)
  the lock            a second launcher on the host is refused while the
                      first holds it
  the foreground      a foreground launcher writes the record when its
                      server answers, and removes only its own
"""

from __future__ import annotations

import os
import stat
import subprocess
import sys

import pytest

from delfin import cli_voila
from delfin.dashboard import server_record as SR
from conftest import child_env


def _record(**extra):
    rec = SR.own_record()
    rec.update({"port": 8899, "ip": "127.0.0.1", "scheme": "http",
                "url": "/voila/render/delfin_voila_runtime/delfin_dashboard.ipynb",
                "resume_path": "/voila/render/delfin_voila_runtime/delfin_resume.ipynb",
                "token": "k" * 43, "started_at": 0.0, "detached": True})
    rec.update(extra)
    return rec


def test_the_record_is_the_owners_alone():
    path = SR.write(_record())
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700
    assert SR.read()["token"] == "k" * 43


def test_a_running_launcher_is_believed():
    SR.write(_record())
    assert SR.running() is not None


def test_a_dead_launcher_is_a_stale_record_and_is_removed(tmp_path):
    proc = subprocess.Popen([sys.executable, "-c", "pass"],
                            env=child_env(tmp_path))
    proc.wait()
    SR.write(_record(pid=proc.pid, proc_start="12345"))
    assert SR.running() is None
    assert not SR.record_file().exists(), "a stale record is replaced"


def test_a_recycled_pid_is_not_the_launcher():
    """Same pid, another start time: another process took the number."""
    SR.write(_record(proc_start="1"))
    assert SR.alive(SR.read()) is False
    assert SR.running() is None


def test_a_record_of_another_host_says_nothing_here():
    SR.write(_record(host="another-node"))
    assert SR.alive(SR.read("another-node")) is False
    assert SR.read() is None


def test_withdraw_removes_only_its_own_record(tmp_path):
    path = SR.write(_record(pid=1, proc_start="1"))
    SR.withdraw(path)
    assert path.exists(), "a record of another launcher is left alone"
    path = SR.write(_record())
    SR.withdraw(path)
    assert not path.exists()


def test_the_lock_admits_one_launcher_per_host():
    first, second = SR.Lock(), SR.Lock()
    assert first.acquire() is True
    try:
        assert second.acquire() is False
    finally:
        first.release()
    assert second.acquire() is True
    second.release()


@pytest.fixture
def quiet_launch(monkeypatch, tmp_path):
    from delfin.agent import process_guard

    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("off")))
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(tmp_path / "runtime"))
    for name in ("TMUX", "STY", "SLURM_JOB_ID", "DELFIN_PREVIEW_BRANCH",
                 "DELFIN_PREVIEW_WORKTREE", "DELFIN_LIFELINE_PID",
                 "DELFIN_LIFELINE_TICKS", "DELFIN_LAUNCH_CWD"):
        monkeypatch.setenv(name, "-")
        monkeypatch.delenv(name)
    monkeypatch.setattr(cli_voila.subprocess, "Popen",
                        lambda *a, **k: pytest.fail("no second server"))
    monkeypatch.setattr(SR, "port_answers", lambda record, timeout=0.5: True)
    monkeypatch.setattr(SR, "kernels", lambda record, timeout=5.0: [])


def test_a_second_launch_opens_the_running_server(quiet_launch, capsys):
    SR.write(_record())
    rc = cli_voila.main(["--no-browser"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "already runs" in out
    assert "localhost:8899" in out
    assert "k" * 43 not in out, "no token in output that is not a terminal"


def test_strict_port_refuses_a_server_on_another_port(quiet_launch, capsys):
    SR.write(_record())
    rc = cli_voila.main(["--no-browser", "--strict-port", "--port", "9001"])
    assert rc == 1
    assert "8899" in capsys.readouterr().err


def test_a_branch_preview_starts_its_own_server(monkeypatch):
    """A preview runs its branch's server beside the real one: neither
    the record nor the lock is its."""
    monkeypatch.setenv("DELFIN_PREVIEW_BRANCH", "b")
    monkeypatch.setenv("DELFIN_PREVIEW_WORKTREE", "/w")
    assert cli_voila._single_server_applies(os.environ) is False
    monkeypatch.delenv("DELFIN_PREVIEW_BRANCH")
    monkeypatch.delenv("DELFIN_PREVIEW_WORKTREE")
    monkeypatch.setenv("SLURM_JOB_ID", "1")
    assert cli_voila._single_server_applies(os.environ) is False
    monkeypatch.delenv("SLURM_JOB_ID")
    assert cli_voila._single_server_applies(os.environ) is True


def test_a_foreground_launcher_records_its_server_once_it_answers(
        monkeypatch, tmp_path):
    from delfin.agent import process_guard

    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("off")))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(tmp_path / "runtime"))
    for name in ("TMUX", "STY", "SLURM_JOB_ID", "DELFIN_LIFELINE_PID",
                 "DELFIN_LIFELINE_TICKS", "DELFIN_LAUNCH_CWD",
                 "DELFIN_VOILA_TOKEN"):
        monkeypatch.setenv(name, "-")
        monkeypatch.delenv(name)
    root = tmp_path / "root"
    root.mkdir()
    nb = tmp_path / "delfin_dashboard.ipynb"
    nb.write_text('{"cells": []}', encoding="utf-8")
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(nb))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_trust_notebook", lambda n: None)
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port",
                        lambda host, port, timeout=10.0: True)
    seen = {}
    withdrawn = []
    monkeypatch.setattr(SR, "withdraw", lambda path=None: withdrawn.append(path))

    class _Proc:
        pid = 4242

        def wait(self, timeout=None):
            # The record is written from a thread once the port answers;
            # the launcher waits for its server meanwhile.
            import time
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline and SR.read() is None:
                time.sleep(0.02)
            seen["record"] = SR.read()
            return 0

    def _popen(cmd, env=None, **k):
        seen["env"] = dict(env or {})
        return _Proc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
    with pytest.raises(SystemExit):
        cli_voila.main(["--no-browser", "--port", "9123", "--foreground"])
    rec = seen["record"]
    assert rec["pid"] == os.getpid() and rec["server_pid"] == 4242
    assert rec["port"] == 9123 and len(rec["token"]) >= 32
    assert rec["proc_start"] == SR.own_record()["proc_start"]
    assert rec["detached"] is False
    assert rec["resume_path"].endswith("delfin_resume.ipynb")
    assert withdrawn, "the launcher takes its record back when it ends"
    # A kernel finds the dashboard's own address (the header's "New
    # dashboard") in its environment, as it finds the resume address.
    assert seen["env"]["DELFIN_DASHBOARD_URL_PATH"].startswith("/voila/render/")
