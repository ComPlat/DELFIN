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
