"""Tests for the LLM-backed conversation compaction path.

The engine calls a cheap model to summarise the older half of the
conversation; on any failure it falls back to extractive summarisation.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pytest

from delfin.agent.engine import AgentEngine


@dataclass
class _FakeEvent:
    type: str
    text: str = ""


def _bare_engine_with_client(text_chunks: list[str] | None, *, raise_exc: bool = False):
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None

    class _FakeClient:
        def stream_message(self, messages=None, system=None, max_tokens=None):
            if raise_exc:
                raise RuntimeError("client unavailable")
            for chunk in (text_chunks or []):
                yield _FakeEvent(type="text_delta", text=chunk)
            yield _FakeEvent(type="message_delta")

    eng.client = _FakeClient()
    return eng


def test_llm_summarize_collects_text_deltas():
    eng = _bare_engine_with_client(["## Goal\nFinish ", "the task.\n"])
    out = eng._llm_summarize_old_messages([
        {"role": "user", "content": "user goal"},
        {"role": "assistant", "content": "ack"},
    ])
    assert "Finish the task" in out


def test_llm_summarize_returns_empty_on_exception():
    eng = _bare_engine_with_client(None, raise_exc=True)
    out = eng._llm_summarize_old_messages([
        {"role": "user", "content": "anything"},
    ])
    assert out == ""


def test_llm_summarize_with_no_messages_returns_empty():
    eng = _bare_engine_with_client(["ignored"])
    assert eng._llm_summarize_old_messages([]) == ""


def test_llm_summarize_with_no_client_returns_empty():
    eng = _bare_engine_with_client(["ignored"])
    eng.client = None
    assert eng._llm_summarize_old_messages([{"role": "user", "content": "x"}]) == ""


def test_llm_summarize_skipped_on_cli_backend():
    """The CLI client uses a persistent subprocess scoped to the user's
    session id; pushing the summariser through it would graft the side-
    conversation onto the live transcript. The compactor must therefore
    decline LLM mode for cli backends and let the caller fall back to
    extractive summarisation."""
    eng = _bare_engine_with_client(["should-not-be-reached"])
    eng.backend = "cli"
    out = eng._llm_summarize_old_messages([
        {"role": "user", "content": "anything"},
        {"role": "assistant", "content": "ack"},
    ])
    assert out == ""


def test_llm_summarize_runs_on_api_backend():
    eng = _bare_engine_with_client(["API-summary"])
    eng.backend = "api"
    out = eng._llm_summarize_old_messages([
        {"role": "user", "content": "go"},
    ])
    assert out == "API-summary"


def test_llm_summarize_caps_giant_messages_individually():
    """A single 200kB tool_result must not flood the summariser's input;
    the per-message cap (4000 chars) keeps it tractable."""
    eng = _bare_engine_with_client(["summary OK"])
    huge = "X" * 200_000
    out = eng._llm_summarize_old_messages([
        {"role": "user", "content": "goal"},
        {"role": "tool", "content": huge},
    ])
    # Reaches the fake client (returns "summary OK") which means the
    # summariser did NOT explode on the huge input
    assert out == "summary OK"


def test_llm_summarize_flattens_structured_content():
    """Some backends pass content as a list of part-dicts instead of a
    plain string. The summariser should flatten those into text."""
    eng = _bare_engine_with_client(["ok"])
    out = eng._llm_summarize_old_messages([
        {"role": "assistant", "content": [
            {"type": "text", "text": "first part"},
            {"type": "text", "text": "second part"},
        ]},
    ])
    assert out == "ok"


def test_compact_history_uses_llm_path_when_available(monkeypatch):
    """End-to-end smoke: _compact_history should swap in the LLM summary
    when the LLM path succeeds, and the resulting summary message should
    contain the LLM's text rather than extracted chars."""
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.backend = "api"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0

    # Compaction is token-driven now: a small window puts this 13-message
    # history over the auto-compact threshold (the 12-message count is only
    # a floor, no longer a trigger).
    eng.context_window_tokens = 50
    for i in range(13):
        role = "user" if i % 2 == 0 else "assistant"
        eng.messages.append({"role": role, "content": f"msg #{i} body content"})

    class _FakeClient:
        def stream_message(self, **kw):
            yield _FakeEvent(type="text_delta", text="LLM_SUMMARY_HERE")

        def kill(self):
            pass

    eng.client = _FakeClient()
    # Skip the user_settings probe — keep llm path on
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"llm_compaction": True}},
        raising=False,
    )

    eng._compact_history()
    # The first message after compaction should contain the LLM summary
    assert any("LLM_SUMMARY_HERE" in str(m.get("content", "")) for m in eng.messages)
    assert eng.last_compaction_info["messages_compacted"] > 0


def test_compact_history_falls_back_to_extractive_on_llm_failure(monkeypatch):
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.backend = "api"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0

    eng.context_window_tokens = 50  # token pressure drives compaction
    for i in range(13):
        role = "user" if i % 2 == 0 else "assistant"
        eng.messages.append({"role": role, "content": f"important content {i}"})

    class _BrokenClient:
        def stream_message(self, **kw):
            raise RuntimeError("API down")

        def kill(self):
            pass

    eng.client = _BrokenClient()
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"llm_compaction": True}},
        raising=False,
    )

    eng._compact_history()
    # Extractive path still produced *something* in the summary message
    summary_msg = eng.messages[0]
    assert "Conversation summary" in summary_msg["content"]
    assert "important content" in summary_msg["content"]


def test_compact_history_skipped_when_setting_disables_llm(monkeypatch):
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.backend = "api"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0

    eng.context_window_tokens = 50  # ensure compaction actually fires
    for i in range(13):
        role = "user" if i % 2 == 0 else "assistant"
        eng.messages.append({"role": role, "content": f"item {i}"})

    called = {"n": 0}

    class _Spy:
        def stream_message(self, **kw):
            called["n"] += 1
            yield _FakeEvent(type="text_delta", text="should-not-appear")

        def kill(self):
            pass

    eng.client = _Spy()
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"llm_compaction": False}},
        raising=False,
    )
    eng._compact_history()
    assert called["n"] == 0  # LLM was NOT called
