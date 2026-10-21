"""Phase 1 (KG4): a chat turn writes its agent_metrics line.

The metrics contract agent_metrics.py states ("one JSONL line per
finished agent turn") was kept by exactly one caller — the dashboard's
turn footer (delfin/dashboard/tab_agent.py:21585). The interactive chat
path (`delfin-agent chat`, and the window's worker turn) never recorded
a line: worker.turn (delfin/agent/service/worker.py:657) computes every
number turn.ended needs (durationMs, inputTokens, outputTokens, costUsd,
toolCalls at worker.py:718-726) and sends them to the UI, and nothing
calls record_turn.

This test drives the real turn path through the public surface —
repl.run_turn via the worker's turn — with a redirected home
(state_paths.RedirectedUserState, which moves agent_metrics._LOG_PATH)
and a fake engine client, and asserts one line appears whose values
match what the turn actually measured.

Red on the unpatched code (0 lines); green after run_turn records.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import agent_metrics
from delfin.agent import repl
from delfin.agent import state_paths


class _FakeClient:
    """Delivers what a real chat client delivers: a streamed answer."""

    model = "fake-model"

    def stream(self, messages=None, **kwargs):
        def _gen():
            yield {"type": "text", "text": "hello"}
        return _gen()


class _FakeEngine:
    """Duck-typed enough for repl.run_turn."""

    def __init__(self, tmp_home: Path):
        self.session_id = "sess-test-1"
        self.client = _FakeClient()
        self.messages = []
        self.token_usage = {"input": 0, "output": 0, "cached": 0}
        self.cost_usd = 0.0
        self.mode = "solo"
        self._model_name = "fake-model"
        self._tmp_home = tmp_home

    @property
    def model(self) -> str:
        # The real engine exposes its model name as a property
        # (delfin/agent/engine.py:1281); run_turn's metrics hook reads
        # that shape, so the fake must match it.
        return self._model_name

    def stream_response(self, **kwargs):
        on_token = kwargs.get("on_token")
        if on_token:
            on_token("hello")
        self.token_usage["input"] += 120
        self.token_usage["output"] += 45
        return "hello"

    def estimate_context_tokens(self):
        return 0

    def delegate_spend(self):
        return agent_metrics.DelegateLedger()

    def get_status(self):
        return {"input_tokens": self.token_usage["input"],
                "output_tokens": self.token_usage["output"],
                "total_cost_usd": self.cost_usd}


@pytest.fixture()
def redirected_home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    with state_paths.RedirectedUserState(home):
        yield home


def _metrics_lines(home: Path) -> list[dict]:
    p = home / ".delfin" / "agent_metrics.jsonl"
    if not p.exists():
        return []
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            out.append(json.loads(line))
    return out


def test_one_chat_turn_writes_one_metrics_line(redirected_home, monkeypatch):
    """One finished chat turn -> exactly one agent_metrics line, carrying
    what the turn measured (not zeros)."""
    engine = _FakeEngine(redirected_home)
    result = repl.run_turn(engine, "say hello",
                           sink=lambda item: None)
    assert result.error == ""
    assert result.output_tokens == 45

    lines = _metrics_lines(redirected_home)
    assert len(lines) == 1, (
        f"expected exactly 1 metrics line after one chat turn, "
        f"got {len(lines)}")
    rec = lines[0]
    # The record carries what the turn measured, not zeros.
    assert rec["session_id"] == "sess-test-1"
    assert rec["model"] == "fake-model"
    assert rec["input_tokens"] == 120
    assert rec["output_tokens"] == 45
    assert rec["duration_s"] > 0.0
    assert rec["ts"] > 0.0


def test_a_failing_metrics_write_does_not_break_the_turn(redirected_home):
    """Best-effort contract (agent_metrics.py:274): a metrics write that
    fails never breaks the turn it describes."""
    engine = _FakeEngine(redirected_home)

    def _boom(metrics, *, path=None):
        raise RuntimeError("disk on fire")

    monkey = pytest.MonkeyPatch()
    monkey.setattr(agent_metrics, "record_turn", _boom)
    try:
        result = repl.run_turn(engine, "say hello",
                               sink=lambda item: None)
    finally:
        monkey.undo()
    assert result.error == ""
    assert result.text == "hello"
