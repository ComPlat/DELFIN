"""A headless run keeps going while its answer announces undone work.

``cmd_run`` (the one-shot door both ``chat -p`` and ``run`` take) used to
stop for good after one ``_run_once``: a run whose answer ends by
announcing a step it did not take ("I'll continue the review now") then
stopped with that intention as its whole output.

Contract pinned here: continuations happen only on cmd_run's one-shot
path (never the interactive session), are bounded by
``agent.announce_continuations`` / ``DELFIN_ANNOUNCE_CONTINUATIONS``
(default 5), count each continuation turn against the run's ceilings,
and stop on a real final answer, a ceiling, an engine stop, an error,
or no progress (a continuation that called no tool).
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


_ANNOUNCING = ("The test file is clear. I'll continue the review now.")
_FINAL = "Done: the review is complete, all checks pass."


class _AnnouncingEngine:
    """Scripted replies; each turn's answer and tool calls are scripted."""

    session_id = "sess-ann-1"
    token_usage = {"input": 1, "output": 1}

    def __init__(self, turns):
        # turns: list of (text, tool_calls) per call, last repeated.
        self.turns = list(turns)
        self.calls = 0
        self.messages = [{"role": "user", "content": "do it"}]
        self.run_budget_turns = 0

    def stream_response(self, **kwargs):
        turn = (self.turns[self.calls]
                if self.calls < len(self.turns) else self.turns[-1])
        self.calls += 1
        self.messages.append({"role": "assistant", "content": turn[0]})
        return turn[0], list(turn[1])

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None


def _result(engine, engine_call) -> dict:
    """A ``_run_once``-shaped result from one scripted turn.

    ``engine`` is the fake engine whose ``_tools_after_last_text`` attribute
    is set from the turn, mirroring how the real ``_run_once`` carries that
    flag (on the engine, not in the returned dict -- whose keys are pinned
    to ``repl.TURN_KEYS``).

    The tuple is (text, tool_calls, tools_after_last_text?). A missing third
    element leaves the flag at the real convention -- True when the turn
    ended ON a tool round (tools non-empty), False when it closed on text.
    A caller that models the measured failure (many tool rounds, then a
    closing announcing text) passes ``tools_after_last_text=False``
    explicitly so the turn reads as closed on text despite non-empty
    tool_calls.
    """
    text, tools = engine_call[:2]
    tools_after = (engine_call[2] if len(engine_call) > 2 else None)
    if tools_after is None:
        tools_after = bool(tools)
    try:
        engine._tools_after_last_text = bool(tools_after)
    except AttributeError:
        pass
    return {"text": text, "tool_calls": tools,
            "input_tokens": 1, "output_tokens": 1,
            "error": "", "stop_reason": "done"}


def _wire(monkeypatch, tmp_path, engine, run_once):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: engine)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda *a, **k: "sess-ann-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda *a, **k: None)
    monkeypatch.setattr("sys.stdin",
                        io.TextIOWrapper(io.BytesIO(b"")))
    monkeypatch.setattr(agent_cli, "_run_once", run_once)
    monkeypatch.chdir(Path.cwd())


@pytest.fixture
def _script(monkeypatch):
    """Scripted client: what _run_once returns per call."""
    results: list = []

    def run_once(engine, prompt, **kw):
        idx = min(engine.calls - 1 if engine.calls else 0,
                  len(engine.turns) - 1)
        results.append(_result(engine, engine.turns[idx]))
        if engine.calls < len(engine.turns):
            engine.calls += 1
        return dict(results[-1])

    return run_once


# -- the behaviours -----------------------------------------------------------

def test_many_tools_then_announcing_close_gets_one_continuation(
        monkeypatch, tmp_path):
    """The measured failure: many tool rounds, then a closing announcing
    answer. The turn closed on text, so one continuation happens and the
    nudge's answer ends the run."""
    engine = _AnnouncingEngine([
        (_ANNOUNCING, [{"name": "read_file", "input": {}}], False),
        (_FINAL, []),
    ])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        idx = min(engine_.calls, len(engine_.turns) - 1)
        engine_.calls += 1
        return _result(engine_, engine_.turns[idx])

    _wire(monkeypatch, tmp_path, engine, run_once)
    rc = agent_cli.main(["-p", "do it", "--output-format", "json"])
    assert rc == 0
    # original turn that closed on an announcement, then ONE continuation
    assert len(calls) == 2
    assert engine.calls == 2


def test_continuations_are_bounded_by_the_configurable_cap(
        monkeypatch, tmp_path):
    """An agent that announces on every turn is stopped at the cap."""
    monkeypatch.setenv("DELFIN_ANNOUNCE_CONTINUATIONS", "2")
    # The original turn closes on an announcement (no tool after). Every
    # continuation makes a tool call (real progress) but still announces,
    # so only the cap -- 2 -- limits the run, not the no-progress rule.
    engine = _AnnouncingEngine([(_ANNOUNCING, [])])
    progress = (_ANNOUNCING, [{"name": "bash", "input": {}}], True)
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return _result(engine_,
                       progress if len(calls) > 1 else (_ANNOUNCING, []))

    _wire(monkeypatch, tmp_path, engine, run_once)
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    # one original turn + 2 continuations, not 10
    assert len(calls) == 3


def test_no_progress_stops_the_continuations(monkeypatch, tmp_path):
    """A continuation that calls no tool is not re-nudged."""
    engine = _AnnouncingEngine([(_ANNOUNCING, [])])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return _result(engine_, (_ANNOUNCING, []))

    _wire(monkeypatch, tmp_path, engine, run_once)
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    assert len(calls) == 2   # original + one no-progress continuation


def test_a_ceiling_ends_the_continuations(monkeypatch, tmp_path):
    """A reached run ceiling stops continuations before any nudge runs."""
    from delfin.agent import run_limits as rl
    engine = _AnnouncingEngine([(_ANNOUNCING, [])])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return _result(engine_, (_ANNOUNCING, []))

    _wire(monkeypatch, tmp_path, engine, run_once)
    monkeypatch.setattr(rl, "ceiling_reason", lambda eng: "run_budget_s")
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    # the original turn only: the ceiling is checked before the first
    # continuation, so no nudge is spent against a spent ceiling
    assert len(calls) == 1


def test_zero_disables_the_continuations(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_ANNOUNCE_CONTINUATIONS", "0")
    engine = _AnnouncingEngine([(_ANNOUNCING, [])])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return _result(engine_, (_ANNOUNCING, []))

    _wire(monkeypatch, tmp_path, engine, run_once)
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    assert len(calls) == 1


def test_a_real_final_answer_never_gets_a_continuation(monkeypatch, tmp_path):
    engine = _AnnouncingEngine([(_FINAL, [])])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return _result(engine_, (_FINAL, []))

    _wire(monkeypatch, tmp_path, engine, run_once)
    rc = agent_cli.main(["-p", "do it", "--output-format", "json"])
    assert rc == 0
    assert len(calls) == 1


def test_an_error_or_engine_stop_is_left_alone(monkeypatch, tmp_path):
    out = {"text": _ANNOUNCING, "tool_calls": [], "input_tokens": 1,
           "output_tokens": 1, "error": "boom", "stop_reason": "error"}
    engine = _AnnouncingEngine([(_ANNOUNCING, [])])
    calls = []

    def run_once(engine_, prompt, **kw):
        calls.append(1)
        return dict(out)

    _wire(monkeypatch, tmp_path, engine, run_once)
    rc = agent_cli.main(["-p", "do it", "--output-format", "json"])
    assert rc != 0
    assert len(calls) == 1


def test_resolution_default_env_and_settings(monkeypatch, tmp_path):
    monkeypatch.delenv("DELFIN_ANNOUNCE_CONTINUATIONS", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    assert agent_cli._resolve_announce_continuations() == 5
    monkeypatch.setenv("DELFIN_ANNOUNCE_CONTINUATIONS", "3")
    assert agent_cli._resolve_announce_continuations() == 3
    monkeypatch.setenv("DELFIN_ANNOUNCE_CONTINUATIONS", "junk")
    assert agent_cli._resolve_announce_continuations() == 5
