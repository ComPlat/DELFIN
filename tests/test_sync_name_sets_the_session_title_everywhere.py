"""``/sync NAME`` names the session with its one title field on every surface.

Input: ``/sync NAME`` typed in the terminal REPL and in the dashboard. In a
service session it already set the session title (worker session.setTitle);
in the REPL and the dashboard it switched SimpleX on and the name went
nowhere, so the phone chat, the session list and a resume showed three
different names.

Output: the REPL sets the agent's ``session_title`` (the field its saves
and the worker's setTitle use) and its turn-end save carries it; the
dashboard sets its session title (``state["_messenger_titles"][sid]``, what
its save writes as ``title`` and its messenger adapter answers) and saves.
Both still switch the session on; ``/sync``, ``on``, ``off`` and ``status``
name nothing.
"""
from __future__ import annotations

import ast
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import repl, repl_commands
from delfin.messenger import activation


@pytest.fixture(autouse=True)
def _records(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")


def _agent(tmp_path, sid="s1"):
    agent = repl.TerminalAgent.__new__(repl.TerminalAgent)
    agent.engine = SimpleNamespace(session_id=sid)
    agent.opts = SimpleNamespace(cwd=tmp_path)
    shown = []
    agent.transcript = SimpleNamespace(chrome=shown.append,
                                       theme=SimpleNamespace(dim=str, yellow=str, red=str))
    agent._say_chrome = shown.append
    agent._auto_controller = lambda: None
    return agent, shown


def test_sync_name_in_the_repl_sets_the_session_title(tmp_path, monkeypatch):
    agent, shown = _agent(tmp_path)
    assert agent._handle_line("/sync Butan Konformere") == ""
    assert activation.enabled("s1")
    assert agent.session_title == "Butan Konformere"
    saved = {}
    from delfin.agent import cli
    monkeypatch.setattr(cli, "_save_session", lambda engine, root, title="": saved.update(title=title))
    agent._checkpoint_session()
    assert saved == {"title": "Butan Konformere"}


@pytest.mark.parametrize("args", ["", "on", "status"])
def test_plain_sync_in_the_repl_names_nothing(tmp_path, args):
    agent, _ = _agent(tmp_path)
    agent._handle_line(("/sync " + args).strip())
    assert not getattr(agent, "session_title", "")


def _nested(name):
    from delfin.dashboard import tab_agent
    tree = ast.parse(inspect.getsource(tab_agent))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(name)


def test_sync_name_in_the_dashboard_sets_the_session_title(monkeypatch):
    class Local:
        def __init__(self, sid, callback):
            self.path = Path("/run/fake-endpoint.sock")

    monkeypatch.setattr(activation, "LocalSession", Local)
    saves = []
    ns = {"state": {}, "_ensure_engine": lambda: SimpleNamespace(session_id="d1"),
          "_auto_save_session": lambda: saves.append(True)}
    for name in ("_answer_ask_from_chat", "_on_stop", "_on_send", "_on_approve", "_on_deny",
                 "_handle_slash_command", "_awaiting_user_answer"):
        ns[name] = lambda *a, **k: None
    ns["input_textarea"] = SimpleNamespace(value="")
    exec(_nested("_messenger_command"), ns)
    out = ns["_messenger_command"]("Butan Konformere")
    assert "aktiviert" in out and activation.enabled("d1")
    assert ns["state"]["_messenger_titles"]["d1"] == "Butan Konformere"
    assert saves
    ns["_messenger_command"]("off")
    assert ns["state"]["_messenger_titles"]["d1"] == "Butan Konformere"
