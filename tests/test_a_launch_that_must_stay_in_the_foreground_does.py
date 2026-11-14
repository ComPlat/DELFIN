"""What ``delfin-voila`` does today in the launches that must not change.

Characterization, written before the launcher learned to run the server
detached from the terminal. Each case is a caller that owns the
dashboard's lifetime itself and needs the launcher in the foreground:

  no terminal                a script, a test, a service manager (systemd
                             on the Pi, launchd): the server runs in
                             this process, its exit is the launcher's
  inside tmux                the Windows launcher runs it in a tmux
                             session and finds it by that session's name
                             in the dashboard note
  inside a SLURM job         ``delfin-voila --job`` execs it in the job;
                             the job ends with it
  the Pi's unit              loopback, its own token, stays up without a
                             window (DELFIN_DASHBOARD_STAY_UP)
  --keep without a terminal  arms Keep dashboard for every window
  Ctrl+C                     stops the server and says so

Input: ``cli_voila.main`` with the server process, the port probe and the
notebook trust replaced by stand-ins. Output: whether the server command
was started from this process, with which arguments and environment.
"""

from __future__ import annotations

import sys

import pytest

from delfin import cli_voila
from delfin.dashboard import resume_server as R


class _Stdin:
    """A stand-in for sys.stdin that says whether it is a terminal."""

    def __init__(self, tty: bool):
        self._tty = tty

    def isatty(self) -> bool:
        return self._tty

    def fileno(self) -> int:
        raise OSError("no descriptor in a test")


class _Proc:
    pid = 0

    def __init__(self, interrupt: bool = False):
        self._interrupt = interrupt
        self.terminated = False

    def wait(self, timeout=None):
        if self._interrupt:
            self._interrupt = False
            raise KeyboardInterrupt
        return 0

    def terminate(self):
        self.terminated = True

    def kill(self):
        pass

    def poll(self):
        return 0


@pytest.fixture
def launch(monkeypatch, tmp_path):
    """Run main() with a stand-in server; return what it started."""
    from delfin.agent import process_guard

    def _no_preflight(*_a):
        raise RuntimeError("credential maintenance isolated in launcher tests")

    monkeypatch.setattr(process_guard, "protect", _no_preflight)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(tmp_path / "runtime"))
    root = tmp_path / "root"
    root.mkdir()
    notebook = tmp_path / "pkg" / "delfin_dashboard.ipynb"
    notebook.parent.mkdir()
    notebook.write_text('{"cells": []}', encoding="utf-8")
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    # Set, then removed: monkeypatch then restores each one, also those
    # main() sets for its server (the lifeline, the launch folder).
    for name in ("TMUX", "STY", "SLURM_JOB_ID", "DELFIN_KEEP_SESSIONS",
                 R.STAY_UP_ENV, "DELFIN_VOILA_TOKEN", "INVOCATION_ID",
                 "DELFIN_LIFELINE_PID", "DELFIN_LIFELINE_TICKS",
                 "DELFIN_LAUNCH_CWD", "DELFIN_RESUME_SESSION",
                 "JUPYTER_RUNTIME_DIR", "DELFIN_VOILA_DETACHED"):
        monkeypatch.setenv(name, "-")
        monkeypatch.delenv(name)
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env",
                        lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_trust_notebook", lambda nb: None)
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port",
                        lambda host, port, timeout=10.0: False)
    seen: dict = {"started": []}

    def run(argv, *, tty=False, env=None, interrupt=False):
        for key, value in (env or {}).items():
            monkeypatch.setenv(key, value)
        monkeypatch.setattr(sys, "stdin", _Stdin(tty))
        proc = _Proc(interrupt)

        def _popen(cmd, env=None, **kwargs):
            seen["started"].append({"cmd": list(cmd), "env": dict(env or {}),
                                    "kwargs": kwargs})
            return proc

        monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
        import signal

        sigint = signal.getsignal(signal.SIGINT)
        try:
            seen["rc"] = cli_voila.main(list(argv))
        except SystemExit as exc:
            seen["rc"] = exc.code
        finally:
            # The Ctrl+C path ignores a second Ctrl+C; not in the suite.
            signal.signal(signal.SIGINT, sigint)
        seen["proc"] = proc
        return seen

    return run


def _server_started_here(seen) -> dict:
    """The one jupyter-server command main() started, from this process."""
    servers = [s for s in seen["started"]
               if s["cmd"][1:4] == ["-m", "jupyter", "server"]]
    assert len(servers) == 1, seen["started"]
    return servers[0]


def test_without_a_terminal_the_server_runs_in_this_process(launch):
    seen = launch(["--no-browser", "--port", "9011"], tty=False)
    server = _server_started_here(seen)
    assert "--port=9011" in server["cmd"]
    assert seen["rc"] == 0


def test_inside_tmux_the_server_runs_in_this_process(launch):
    seen = launch(["--no-browser", "--port", "9012"], tty=True,
                  env={"TMUX": "/tmp/tmux-1/default,1,0"})
    _server_started_here(seen)


def test_inside_a_slurm_job_the_server_runs_in_this_process(launch):
    seen = launch(["--no-browser", "--port", "9013"], tty=True,
                  env={"SLURM_JOB_ID": "4711"})
    _server_started_here(seen)


def test_the_pi_unit_runs_in_the_foreground_and_stays_up(launch):
    """scripts/pi/pi_bootstrap.sh: exec delfin-voila --ip 127.0.0.1
    --port 8866 --no-browser --token ..., under systemd (no terminal),
    with DELFIN_DASHBOARD_STAY_UP=1."""
    token = "t" * 43
    seen = launch(["--ip", "127.0.0.1", "--port", "8866", "--no-browser",
                   "--token", token], tty=False, env={R.STAY_UP_ENV: "1"})
    server = _server_started_here(seen)
    assert server["env"]["JUPYTER_TOKEN"] == token
    assert token not in " ".join(server["cmd"]), "the token is in no command line"
    assert not any(a.startswith("--ServerApp.shutdown_no_activity_timeout")
                   for a in server["cmd"]), "a service stays up without a window"


def test_a_foreground_server_ends_after_ten_minutes_without_a_kernel(launch):
    seen = launch(["--no-browser", "--port", "9014"], tty=False)
    server = _server_started_here(seen)
    assert "--ServerApp.shutdown_no_activity_timeout=600" in server["cmd"]


def test_keep_without_a_terminal_arms_every_window(launch):
    import os

    seen = launch(["--no-browser", "--port", "9015", "--keep"], tty=False)
    _server_started_here(seen)
    assert os.environ.get("DELFIN_KEEP_SESSIONS") == "1"


def test_the_windows_launcher_finds_its_tmux_dashboard_in_the_note(
        launch, monkeypatch):
    """tools/windows-launcher/remote_launcher.py runs, inside a tmux
    session it named, ``python -c 'from delfin.cli_voila import main ...'
    --no-browser --ip 127.0.0.1 --strict-port --port N --keep --resume
    latest`` and waits for the dashboard note to name that session."""
    from delfin.agent import where

    monkeypatch.setattr(where, "_tmux_session", lambda: "delfin-win-8870-x")
    seen = launch(["--no-browser", "--ip", "127.0.0.1", "--strict-port",
                   "--port", "8870", "--keep", "--resume", "latest"],
                  tty=True, env={"TMUX": "/tmp/tmux-1/default,1,0"})
    server = _server_started_here(seen)
    assert "--ServerApp.port_retries=0" in server["cmd"]
    note = where.dashboard()
    assert note.get("tmux") == "delfin-win-8870-x"
    assert int(note.get("port") or 0) == 8870
    assert len(str(note.get("token") or "")) >= 32


def test_ctrl_c_stops_the_server_and_says_so(launch, capsys):
    seen = launch(["--no-browser", "--port", "9016"], tty=False,
                  interrupt=True)
    _server_started_here(seen)
    assert seen["proc"].terminated, "Ctrl+C ends the server it started"
    out = capsys.readouterr().out
    assert "Stopping the dashboard and everything it started" in out
    assert "Stopped." in out
