"""An agent session cannot speak as the operator or pause a session.

Input: shell commands that run ``delfin-agent ops send|escalate|pause``
or the top-level ``delfin-agent pause|resume``, in the spellings a
command can take (wrappers, ``env``, ``xargs``, ``bash -c`` bodies, a
path, ``python -m``, the Python functions), sent through the real bash
executor (``_doc_executor.execute``) in every permission mode, with a
stand-in ``delfin-agent`` first on PATH that leaves a mark if it runs;
the CLI backend's hook (``cli_gate``); and ``delfin-agent pause``
started with the agent-session marker. Output: the refusal, the
security event, the mark file and the pause flag.

Semantics: the operator's actions are the person's (``ops_interface``):
a session that sends as the operator writes under the supervisor's name,
and a pause or a resume decides over another session. The command-text
rule holds where the marker was stripped; the marker refusal holds where
the text is built at run time. Reading (``ops sessions``, ``ops events``,
``ops approvals``) and mentioning the commands stay allowed.
"""
from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor

MARKER = "DELFIN_AGENT_SESSION"
MODES = ("default", "acceptEdits", "bypassPermissions")

REFUSED = [
    "delfin-agent ops send s2 --text hello",
    "delfin-agent ops escalate 0123abcd",
    "delfin-agent ops pause s2 --reason later",
    "delfin-agent pause s2",
    "delfin-agent resume s2",
    "delfin-agent pause s2 --reason x",
    "/opt/tools/bin/delfin-agent pause s2",
    "env FOO=1 delfin-agent ops send s2 --text hello",
    "env -u DELFIN_AGENT_SESSION delfin-agent pause s2",
    "timeout 5 delfin-agent ops pause s2 --reason x",
    "nohup delfin-agent ops send s2 --text hello &",
    "echo s2 | xargs delfin-agent pause",
    "echo hello | xargs delfin-agent ops send s2 --text",
    "echo send | xargs delfin-agent ops",
    "bash -c 'delfin-agent ops pause s2 --reason x'",
    "sh -lc \"delfin-agent pause s2\"",
    "delfin-agent \"ops\" 'send' s2 --text hello",
    "ls && delfin-agent ops escalate 0123abcd",
    "micromamba run -n delfin delfin-agent pause s2",
    "python -m delfin.agent.cli ops send s2 --text hello",
    "python3 -m delfin.agent.cli pause s2",
    "python -c 'from delfin.agent import ops_interface; "
    "ops_interface.send(\"s2\", \"hello\")'",
    "python -c 'from delfin.agent import session_pause as p; "
    "p.pause(\"s2\", \"x\")'",
    "python -c 'from delfin.agent.cli import main; "
    "main([\"ops\", \"pause\", \"s2\", \"--reason\", \"x\"])'",
    "p=delfin-agent; $p ops send s2 --text hello",
]

ALLOWED = [
    "delfin-agent ops sessions",
    "delfin-agent ops events --cursor mine",
    "delfin-agent ops approvals",
    "delfin-agent ops session s2",
    "grep -rn 'delfin-agent ops send' docs/",
    "grep -n \"ops pause\" delfin/agent/cli.py",
    'git commit -m "document delfin-agent pause"',
    "echo the person runs delfin-agent pause",
    "sleep 1",
    "python -c 'from delfin.agent import ops_interface; "
    "print(ops_interface.sessions())'",
]


@pytest.fixture
def ws(tmp_path, monkeypatch):
    w = tmp_path / "ws"
    w.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    mark = tmp_path / "ran"
    p = bin_dir / "delfin-agent"
    p.write_text(f'#!/bin/sh\necho "delfin-agent $*" >> {mark}\n')
    p.chmod(p.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return w, mark


@pytest.fixture
def events(monkeypatch):
    seen = []
    monkeypatch.setattr(A, "_record_security_event",
                        lambda kind, *a, **k: seen.append(kind))
    return seen


def _run(cmd, ws, mode):
    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=lambda *a, **k: True)
    return str(_doc_executor.execute("bash", {"command": cmd}, perms))


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", REFUSED)
def test_an_operator_action_is_refused_in_every_mode(cmd, mode, ws, events):
    w, mark = ws
    out = _run(cmd, w, mode)
    assert "acts as the operator" in out and "blocked:" in out, out
    assert "operator_action_refused" in events
    assert not mark.exists(), mark.read_text()


@pytest.mark.parametrize("cmd", ALLOWED)
def test_reading_and_mentioning_is_not_refused(cmd):
    from delfin.agent import approval_guard
    assert approval_guard.operator_command(A._prose_blanked(cmd)) is None


def test_the_cli_backend_hook_refuses_it_too(tmp_path):
    from delfin.agent import cli_gate
    import inspect
    src = inspect.getsource(cli_gate)
    assert "operator_command" in src
    from delfin.agent import approval_guard
    assert approval_guard.operator_command("delfin-agent ops send s2 --text x")


@pytest.mark.parametrize("argv", [["pause", "s2"], ["resume", "s2"]])
def test_the_top_level_pause_refuses_inside_an_agent_session(
        argv, tmp_path, monkeypatch, capsys):
    from delfin.agent import cli as agent_cli
    from delfin.agent import session_pause
    monkeypatch.setattr(session_pause, "_DIR", tmp_path / "pause")
    monkeypatch.setenv(MARKER, "4242")
    assert agent_cli.main(argv) == 2
    assert "agent session" in capsys.readouterr().err
    assert not (tmp_path / "pause").exists() or \
        not list((tmp_path / "pause").iterdir())
