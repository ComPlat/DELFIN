"""The terminal session and the one-shot run bound a turn the same way.

Input: one workspace, one prompt, driven through ``cli.main`` as an
interactive session (a terminal on stdin, the real ``TerminalAgent``) and as
``chat -p``, with the model's HTTP stream scripted. Output: the first
request each path sends. Asserted:

* the output ceiling (``max_tokens``) is the same: the role's budget. ``chat``
  declares ``--max-tokens`` default 0 ("the role decides"), and ``cmd_run``
  turned that 0 into 4096 with ``args.max_tokens or 4096`` -- an eighth of the
  solo budget, on the -p path only;
* ``run --max-tokens 0`` means the role's budget, as its help says; ``run``
  without the flag keeps its historical 4096 (the scheduler and the
  benchmark depend on it);
* a session over its context budget starts fresh on both paths, as it does
  in the dashboard (``prepare_fresh_turn``): the terminal's
  ``ReplOptions.context_budget`` defaulted to 0, which disables it, and the
  one-shot path never called it.

Every store is under ``tmp_path``; no model endpoint is contacted.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import mcp_client as M

_OLD = "OLD-HISTORY-4c1d "


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
    """``drive(door, argv=(), over_budget=False) -> (request, engine)``.

    ``door`` is "terminal", "chat-p" or "run".
    """
    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
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

    def run(door: str, argv=(), *, over_budget: bool = False):
        requests: list[dict] = []
        built: dict = {}

        def _create(**kwargs):
            requests.append(kwargs)
            return _answer()

        def _build(args):
            engine = real_build(args)
            engine.client.client.chat.completions.create = _create
            if over_budget:
                # A restored conversation well past a small budget.
                engine.context_budget = 200
                engine.messages = [
                    {"role": "user", "content": _OLD + "x" * 4000},
                    {"role": "assistant", "content": _OLD + "y" * 4000},
                ]
            built["engine"] = engine
            return engine

        monkeypatch.setattr(agent_cli, "_build_engine", _build)
        common = ["--cwd", str(ws), "--provider", "ollama", "--model", "m"]
        if door == "terminal":
            stdin = io.StringIO("do the task\n")
            stdin.isatty = lambda: True        # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            rc = agent_cli.main([*common, *argv])
        else:
            stdin = io.StringIO("")
            stdin.isatty = lambda: False       # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            if door == "run":
                rc = agent_cli.main(["run", "do the task", *common, *argv])
            else:
                rc = agent_cli.main(["-p", "do the task", *common, *argv])
        assert rc == 0
        assert requests, "no request reached the model"
        return requests[0], built["engine"]

    run.home = str(home)
    return run


def _role_budget(engine) -> int:
    return engine.max_tokens_for_role(engine.current_role or "solo_agent")


def test_the_terminal_and_the_one_shot_run_have_one_output_ceiling(drive):
    tui, engine = drive("terminal")
    one_shot, _ = drive("chat-p")
    assert one_shot["max_tokens"] == tui["max_tokens"] == _role_budget(engine)


def test_run_with_max_tokens_zero_takes_the_role_budget(drive):
    request, engine = drive("run", ["--max-tokens", "0"])
    assert request["max_tokens"] == _role_budget(engine)


def test_run_without_the_flag_keeps_its_ceiling(drive):
    request, _ = drive("run")
    assert request["max_tokens"] == 4096


def _sent(request: dict, home: str) -> list:
    msgs = [m for m in request.get("messages") or []
            if m.get("role") != "system"]
    return json.loads(json.dumps(msgs, default=repr).replace(home, "<HOME>"))


@pytest.mark.parametrize("door", ["terminal", "chat-p"])
def test_a_session_over_its_budget_starts_fresh(drive, door):
    request, engine = drive(door, over_budget=True)
    sent = json.dumps(_sent(request, drive.home))
    assert _OLD not in sent, "the history over the budget was re-sent"
    assert engine.cut_history(), "the dropped history was not archived"


def test_both_paths_send_the_same_fresh_turn(drive):
    tui, _ = drive("terminal", over_budget=True)
    one_shot, _ = drive("chat-p", over_budget=True)
    assert _sent(tui, drive.home) == _sent(one_shot, drive.home)
