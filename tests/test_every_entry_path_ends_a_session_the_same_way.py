"""The terminal session and the one-shot run end a session the same way.

Input: one workspace, one prompt, driven through ``cli.main`` as an
interactive session (a terminal on stdin, the real ``TerminalAgent``) and as
``chat -p``, with the model's HTTP stream scripted. Output: which
session-end stages ran, in which order, and what the outcome stage
recorded. Asserted: the same stages in the same order on both paths --
session report, skill learning, indexing, cycle outcome, episodic memory,
eval loop -- and the outcome is recorded against the task the user gave.

Before, the terminal ran the first three only: an interactive session
added nothing to the outcome history the briefing reads, nothing to the
provider profile, and no episode a later session could recall. The
dashboard records an outcome per cycle and an episode per save.

Every stage is stubbed to a recorder; every store is under ``tmp_path``;
no model endpoint is contacted.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import mcp_client as M


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
    """``drive(interactive) -> list of (stage, detail)`` in call order."""
    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")

    from delfin.agent import engine as engine_mod
    from delfin.agent import episodes, eval_loop, session_end, session_report

    calls: list[tuple[str, str]] = []
    monkeypatch.setattr(session_report, "write_session_report",
                        lambda *a, **k: calls.append(("report", "")))
    monkeypatch.setattr(session_end, "learn_at_session_end",
                        lambda *a, **k: calls.append(("learning", "")))
    monkeypatch.setattr(session_end, "index_at_session_end",
                        lambda *a, **k: calls.append(("indexing", "")))
    monkeypatch.setattr(
        engine_mod.AgentEngine, "record_cycle_outcome",
        lambda self, verdict, user_task, **k: calls.append(
            ("cycle outcome", f"{verdict}:{user_task}")) or {})
    monkeypatch.setattr(episodes, "save_episode",
                        lambda *a, **k: calls.append(
                            ("episode", str(k.get("verdict", "")))))
    monkeypatch.setattr(eval_loop, "maybe_run_scheduled",
                        lambda *a, **k: calls.append(("eval loop", "")))
    real_build = agent_cli._build_engine

    def run(*, interactive: bool) -> list[tuple[str, str]]:
        calls.clear()

        def _build(args):
            engine = real_build(args)
            engine.client.client.chat.completions.create = \
                lambda **kw: _answer()
            return engine

        monkeypatch.setattr(agent_cli, "_build_engine", _build)
        common = ["--cwd", str(ws), "--provider", "ollama", "--model", "m"]
        if interactive:
            stdin = io.StringIO("do the task\n")
            stdin.isatty = lambda: True        # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            rc = agent_cli.main(common)
        else:
            stdin = io.StringIO("")
            stdin.isatty = lambda: False       # type: ignore[method-assign]
            monkeypatch.setattr("sys.stdin", stdin)
            rc = agent_cli.main(["-p", "do the task", *common])
        assert rc == 0
        return list(calls)

    return run


_STAGES = ["report", "learning", "indexing", "cycle outcome", "episode",
           "eval loop"]


def test_the_one_shot_run_runs_every_stage(drive):
    assert [name for name, _ in drive(interactive=False)] == _STAGES


def test_the_terminal_session_runs_the_same_stages(drive):
    one_shot = drive(interactive=False)
    terminal = drive(interactive=True)
    assert terminal == one_shot


def test_the_outcome_is_the_users_task(drive):
    terminal = dict(drive(interactive=True))
    assert terminal.get("cycle outcome") == "PASS:do the task"
    assert terminal.get("episode") == "PASS"
