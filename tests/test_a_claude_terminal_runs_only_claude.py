"""The dashboard's Claude terminal runs the configured CLI and nothing else.

The stock Jupyter terminal manager lets the POST body override the command,
the environment and the directory; the control below shows a request for
``bash`` getting bash. The DELFIN manager must give the same request the
configured CLI, in the configured directory, without the client's
environment and without the dashboard token.
"""

from __future__ import annotations

import json
import sys

import pytest

pytest.importorskip("jupyter_server_terminals")

import terminado.management as TM  # noqa: E402
from jupyter_server_terminals.terminalmanager import TerminalManager  # noqa: E402

from delfin.dashboard import claude_terminal as CT  # noqa: E402

HOSTILE = {"shell_command": ["/bin/bash", "-c", "id"],
           "extra_env": {"LD_PRELOAD": "/tmp/evil.so"},
           "cwd": "/", "name": "../../x", "height": 30, "width": 100}


@pytest.fixture
def spawned(monkeypatch):
    seen = []

    class _Pty:
        def __init__(self, argv, env=None, cwd=None):
            seen.append({"argv": list(argv), "env": dict(env or {}),
                         "cwd": cwd})
            self.clients = []
            self.ptyproc = type("P", (), {"pid": 1, "fd": -1})()
            self.read_buffer = []

        def resize_to_smallest(self):
            pass

    monkeypatch.setattr(TM, "PtyWithClients", _Pty)
    monkeypatch.setattr(TM.TermManagerBase, "start_reading",
                        lambda self, p: None)
    return seen


@pytest.fixture
def configured(monkeypatch, tmp_path):
    argv = [sys.executable, "-c", "print('claude')"]
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps(argv))
    monkeypatch.setenv(CT.CWD_ENV, str(tmp_path))
    monkeypatch.setenv("JUPYTER_TOKEN", "secret-dashboard-token")
    return argv, str(tmp_path)


def test_control_the_stock_manager_obeys_the_request(spawned, configured):
    TerminalManager(shell_command=["/bin/sh"]).new_named_terminal(**HOSTILE)
    assert spawned[0]["argv"] == ["/bin/bash", "-c", "id"]
    assert spawned[0]["env"].get("LD_PRELOAD") == "/tmp/evil.so"


def test_the_request_cannot_choose_command_env_or_directory(
        spawned, configured):
    argv, cwd = configured
    CT.ClaudeOnlyTerminalManager().new_named_terminal(**HOSTILE)
    got = spawned[0]
    assert got["argv"] == argv
    assert got["cwd"] == cwd
    assert "LD_PRELOAD" not in got["env"]
    assert "JUPYTER_TOKEN" not in got["env"]
    assert got["env"]["COLUMNS"] == "100" and got["env"]["LINES"] == "30"


def test_a_name_that_is_not_a_plain_word_is_replaced(spawned, configured):
    mgr = CT.ClaudeOnlyTerminalManager()
    name, _ = mgr.new_named_terminal(**HOSTILE)
    assert name != "../../x" and name.isalnum()


def test_without_a_configured_command_there_is_no_terminal(monkeypatch):
    monkeypatch.delenv(CT.ARGV_ENV, raising=False)
    with pytest.raises(RuntimeError):
        CT.ClaudeOnlyTerminalManager()
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps(["claude"]))   # not absolute
    with pytest.raises(RuntimeError):
        CT.ClaudeOnlyTerminalManager()


def test_at_most_two_terminals(spawned, configured):
    mgr = CT.ClaudeOnlyTerminalManager()
    mgr.new_named_terminal()
    mgr.new_named_terminal()
    with pytest.raises(Exception):
        mgr.new_named_terminal()


def test_a_taken_name_is_not_reused(spawned, configured):
    mgr = CT.ClaudeOnlyTerminalManager()
    first, _ = mgr.new_named_terminal(name="abc")
    second, _ = mgr.new_named_terminal(name="abc")
    assert first == "abc" and second != "abc"


# -- the launcher --------------------------------------------------------------

def _launch(monkeypatch, tmp_path, extra):
    from delfin import cli_voila
    from delfin.agent import process_guard

    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("x")))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(tmp_path / "runtime"))
    root = tmp_path / "home"
    root.mkdir(exist_ok=True)
    nb = tmp_path / "delfin_dashboard.ipynb"
    nb.write_text('{"cells":[]}', encoding="utf-8")
    captured = {}

    class _Proc:
        def wait(self):
            return 0

    def _popen(cmd, env=None, **kw):
        captured["cmd"], captured["env"] = cmd, env
        return _Proc()

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(nb))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/v")
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port",
                        lambda host, port, timeout=10.0: False)
    monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
    monkeypatch.setattr(cli_voila.subprocess, "run",
                        lambda cmd, env, check: (captured.update(
                            cmd=cmd, env=env), type("R", (), {
                                "returncode": 0})())[1])
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    code = 0
    try:
        cli_voila.main(["--no-browser", "--port", "9001"] + extra)
    except SystemExit as exc:
        code = exc.code
    return code, captured


def _fake_claude(monkeypatch, tmp_path):
    exe = tmp_path / "bin" / "claude"
    exe.parent.mkdir(exist_ok=True)
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:/usr/bin:/bin")
    return exe


def test_without_the_flag_there_are_no_terminals(monkeypatch, tmp_path):
    _fake_claude(monkeypatch, tmp_path)
    code, cap = _launch(monkeypatch, tmp_path, [])
    assert code == 0
    assert "--ServerApp.terminals_enabled=False" in cap["cmd"]
    ext = next(a for a in cap["cmd"] if "jpserver_extensions" in a)
    assert "'jupyter_server_terminals': False" in ext
    assert CT.ARGV_ENV not in cap["env"]


def test_the_flag_wires_the_restricted_manager(monkeypatch, tmp_path):
    exe = _fake_claude(monkeypatch, tmp_path)
    code, cap = _launch(monkeypatch, tmp_path, ["--claude-terminal"])
    assert code == 0
    assert "--ServerApp.terminals_enabled=True" in cap["cmd"]
    assert ("--TerminalsExtensionApp.terminal_manager_class="
            "delfin.dashboard.claude_terminal.ClaudeOnlyTerminalManager"
            ) in cap["cmd"]
    assert json.loads(cap["env"][CT.ARGV_ENV]) == [str(exe)]
    assert cap["env"][CT.ENABLED_ENV] == "1"


def test_the_flag_is_refused_off_loopback(monkeypatch, tmp_path):
    _fake_claude(monkeypatch, tmp_path)
    code, _ = _launch(monkeypatch, tmp_path, [
        "--claude-terminal", "--ip", "0.0.0.0", "--allow-remote-bind"])
    assert code == 2


def test_the_flag_is_refused_without_the_cli(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    code, _ = _launch(monkeypatch, tmp_path, ["--claude-terminal"])
    assert code == 2
