"""Compaction is tried before a fresh start, which is the fallback.

Input: an engine whose history exceeds the fresh-start budget (the model's
window; after a switch to a model with a smaller window), driven through
``prepare_fresh_turn`` (dashboard) and ``run_turn`` (terminal) with a
scripted client that answers the summariser's request.

Output: when a forced compaction brings the context under the budget, the
turn is NOT a fresh start: the history begins with the compaction summary
(the task stays in context) and the prompt goes out unchanged. When the
compaction cannot fit (the kept recent messages alone exceed the budget) or
raises, the fresh start runs as before.

Cost: one summariser call, only on a turn whose context is over the
fresh-start budget; a turn under it makes no extra call.
"""
from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import repl as R
from delfin.agent.engine import AgentEngine

from tests.test_t2_engine_fresh_start import _pack_tree


@pytest.fixture(autouse=True)
def _no_background_window_probe(monkeypatch, request):
    """These tests set context_window_tokens by hand after building an
    engine. The construction-time window probe runs on a daemon thread and,
    landing after that assignment, replaced the test's window with the
    resolved one (a compaction then never became due). Only the probe's
    background start is switched off; a test that checks the probe keeps it.
    """
    if "probe_runs_in_background" in request.node.name:
        return
    from delfin.agent.engine import AgentEngine
    original = AgentEngine._refresh_context_window

    def _refresh(self, *, background=False):
        if background:
            return None
        return original(self, background=background)

    monkeypatch.setattr(AgentEngine, "_refresh_context_window", _refresh)


WINDOW = 4_000


class _Client:
    def __init__(self):
        self.summary_calls = 0
        self.sent = None

    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        if "summarisation assistant" in (system or ""):
            self.summary_calls += 1
            yield SimpleNamespace(type="text_delta",
                                  text="## Goal\nRefine the cobalt complex.")
            yield SimpleNamespace(type="message_delta", input_tokens=0,
                                  output_tokens=0, cost_usd=0.0,
                                  stop_reason="end_turn", text="")
            return
        self.sent = list(messages or [])
        yield SimpleNamespace(type="text_delta", text="continuing")
        yield SimpleNamespace(type="message_delta", input_tokens=40,
                              output_tokens=2, cost_usd=0.0,
                              stop_reason="end_turn", text="")


@pytest.fixture(autouse=True)
def _no_archive(monkeypatch):
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "archive_pre_compaction_transcript",
                        lambda *a, **k: None)


def _engine(tmp_path):
    pack = _pack_tree(tmp_path / "pack_root")
    client = _Client()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=pack, backend="cli", mode="quick",
                             pack_dir=pack)
    engine.client = client
    # An API-style backend: the summariser is a stateless side call there.
    engine.backend = "api"
    engine.context_window_tokens = WINDOW
    return engine, client


def _history(n_pairs, size=400):
    out = []
    for i in range(n_pairs):
        out.append({"role": "user", "content": f"step {i} " + "x" * size})
        out.append({"role": "assistant", "content": f"done {i} " + "y" * size})
    return out


def test_a_compaction_that_fits_replaces_the_fresh_start(tmp_path):
    engine, client = _engine(tmp_path)
    engine.messages = _history(30)
    assert R._context_tokens(engine) > WINDOW

    msg = R.prepare_fresh_turn(engine, "go on")

    assert msg == "go on", "not a fresh start"
    assert not getattr(engine, "start_fresh", False)
    assert client.summary_calls == 1
    assert "[Conversation summary" in engine.messages[0]["content"]
    assert "cobalt" in engine.messages[0]["content"]
    assert R._context_tokens(engine) <= WINDOW
    engine.stream_response(user_message=msg)
    assert "cobalt" in client.sent[0]["content"]
    assert client.sent[-1]["content"].startswith("go on")


def test_the_terminal_turn_compacts_first_too(tmp_path):
    engine, client = _engine(tmp_path)
    engine.messages = _history(30)
    R.run_turn(engine, "go on", sink=lambda item: None)
    assert client.summary_calls == 1
    assert R._FRESH_CONTEXT_MARKER not in str(client.sent)
    assert "cobalt" in client.sent[0]["content"]


def test_a_compaction_that_cannot_fit_falls_back_to_the_fresh_start(tmp_path):
    engine, client = _engine(tmp_path)
    # The four kept recent messages alone exceed the window.
    engine.messages = _history(10) + _history(2, size=10_000)
    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg.startswith(R._FRESH_CONTEXT_MARKER)
    assert engine.start_fresh


def test_a_compaction_that_raises_falls_back_to_the_fresh_start(tmp_path):
    engine, client = _engine(tmp_path)
    engine.messages = _history(30)

    def _boom(**kw):
        raise RuntimeError("summariser down")
    engine._compact_history = _boom
    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg.startswith(R._FRESH_CONTEXT_MARKER)


def test_a_turn_under_the_budget_makes_no_extra_call(tmp_path):
    engine, client = _engine(tmp_path)
    engine.messages = _history(2)
    assert R.prepare_fresh_turn(engine, "go on") == "go on"
    assert client.summary_calls == 0
