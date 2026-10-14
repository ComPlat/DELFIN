"""The terminal session, the one-shot run and the dashboard give the model
the same project rules and the same per-turn hints.

Input: one workspace holding an ``AGENTS.md``, one prompt, driven through
``cli.main`` twice -- interactive (a terminal on stdin, the real
``TerminalAgent``) and one-shot (``-p``) -- with the model's HTTP stream
scripted. Output: the first request each path sends (system prompt, the
non-system messages, the tool schemas). Asserted:

* the two requests are byte-identical once the per-run temporary home is
  replaced by a placeholder;
* the project rules are IN that request -- before this, only the
  dashboard loaded ``DELFIN.MD`` / ``AGENTS.md`` (through
  ``memory_context``), so ``delfin-agent init`` scaffolded a file the
  terminal never read;
* a factual ORCA keyword question carries the verify hint the dashboard
  appends, and a bare greeting the greeting hint, on both CLI paths;
* the dashboard composes its turn from the same helpers
  (``turn_context``), checked on its source because building the widget
  tree needs a notebook kernel.

Every store is under ``tmp_path``; no model endpoint is contacted.
"""

from __future__ import annotations

import inspect
import io
import json
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import mcp_client as M

_RULE = "RULE-7f3a: every energy is reported in Hartree."


class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _answer() -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta("done."), "stop")],
                           usage=_Usage())])


@pytest.fixture
def drive(monkeypatch, tmp_path):
    """``drive(interactive, prompt) -> first request kwargs``."""
    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    (ws / "AGENTS.md").write_text(f"# Project rules\n\n{_RULE}\n")
    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    from delfin.agent import provider_profile as PP
    baseline = tmp_path / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    monkeypatch.setattr(PP, "_DEFAULT_PATH", baseline)
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    real_build = agent_cli._build_engine

    def _fresh_home() -> None:
        # Each run starts from the same empty home at the same path: what
        # one run's session end writes (outcome history, episodes) is read
        # by the next run's briefing, and two runs compared must not see
        # each other.
        import shutil
        shutil.rmtree(home)
        ws.mkdir(parents=True)
        (ws / "AGENTS.md").write_text(f"# Project rules\n\n{_RULE}\n")
        monkeypatch.chdir(ws)

    def run(*, interactive: bool, prompt: str = "do the task") -> dict:
        _fresh_home()
        requests: list[dict] = []

        def _create(**kwargs):
            requests.append(kwargs)
            return _answer()

        def _build(args):
            engine = real_build(args)
            engine.client.client.chat.completions.create = _create
            return engine

        monkeypatch.setattr(agent_cli, "_build_engine", _build)
        common = ["--cwd", str(ws), "--provider", "ollama", "--model", "m"]
        if interactive:
            stdin = io.StringIO(prompt + "\n")
            stdin.isatty = lambda: True        # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            rc = agent_cli.main(common)
        else:
            stdin = io.StringIO("")
            stdin.isatty = lambda: False       # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            rc = agent_cli.main(["-p", prompt, *common])
        assert rc == 0
        assert requests, "no request reached the model"
        return requests[0]

    run.home = str(home)
    return run


def _view(request: dict, home: str) -> dict:
    msgs = request.get("messages") or []
    text = json.dumps({
        "system": [m.get("content") for m in msgs if m.get("role") == "system"],
        "messages": [m for m in msgs if m.get("role") != "system"],
        "tools": request.get("tools") or [],
    }, sort_keys=True, default=repr)
    return json.loads(text.replace(home, "<HOME>"))


def _system(request: dict) -> str:
    return "\n".join(str(m.get("content") or "")
                     for m in request.get("messages") or []
                     if m.get("role") == "system")


def _last_user(request: dict) -> str:
    users = [m for m in request.get("messages") or []
             if m.get("role") == "user"]
    content = users[-1].get("content") if users else ""
    if isinstance(content, list):
        content = " ".join(str(p.get("text", "")) for p in content
                           if isinstance(p, dict))
    return str(content or "")


def test_the_terminal_and_the_one_shot_run_send_one_request(drive):
    tui = drive(interactive=True)
    one_shot = drive(interactive=False)
    assert _view(tui, drive.home) == _view(one_shot, drive.home)


@pytest.mark.parametrize("interactive", [True, False],
                         ids=["terminal", "one-shot"])
def test_the_project_rules_reach_the_model(drive, interactive):
    request = drive(interactive=interactive)
    assert _RULE in _system(request)


@pytest.mark.parametrize("interactive", [True, False],
                         ids=["terminal", "one-shot"])
def test_a_keyword_question_carries_the_verify_hint(drive, interactive):
    from delfin.agent import turn_context
    question = "Welches Keyword nutzt ORCA für die Dispersionskorrektur?"
    hint = turn_context.verify_hint(question)
    assert hint, "the question is the shape the hint exists for"
    request = drive(interactive=interactive, prompt=question)
    assert hint in _last_user(request)


@pytest.mark.parametrize("interactive", [True, False],
                         ids=["terminal", "one-shot"])
def test_a_bare_greeting_carries_the_greeting_hint(drive, interactive):
    from delfin.agent import turn_context
    hint = turn_context.greeting_hint("hallo")
    assert hint
    request = drive(interactive=interactive, prompt="hallo")
    assert hint in _system(request)


def test_the_dashboard_composes_its_turn_from_the_same_helpers():
    from delfin.agent import turn_context
    from delfin.dashboard import tab_agent
    source = inspect.getsource(tab_agent)
    for name in ("project_instructions", "greeting_hint"):
        assert f"turn_context.{name}(" in source, name
    # The verify hint is the shared one, not a second copy of the matcher.
    assert tab_agent._build_verify_hint is turn_context.verify_hint
    # And no second loader of the instruction files.
    assert "load_project_memory(" not in source
