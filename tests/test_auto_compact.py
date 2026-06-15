"""Tests for token-budget auto-compaction in AgentEngine."""

from __future__ import annotations

from unittest.mock import MagicMock

from delfin.agent.engine import AgentEngine


def _bare_engine() -> AgentEngine:
    """Build an engine instance without instantiating a real backend."""
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.80
    eng.last_compaction_info = None
    eng.session_id = ""
    eng.backend = "api"
    eng.client = MagicMock()
    eng.current_role_index = 0
    eng.route = ["solo_agent"]
    return eng


def test_estimate_returns_zero_for_empty_history():
    eng = _bare_engine()
    assert eng._estimate_context_tokens() == 0


def test_estimate_grows_with_message_size():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "x" * 4000}]
    assert eng._estimate_context_tokens() == 1000


def test_should_auto_compact_false_below_threshold():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "tiny"}]
    assert not eng._should_auto_compact()


def test_should_auto_compact_true_above_threshold():
    eng = _bare_engine()
    # 100k * 0.80 = 80k tokens budget = 320k chars
    eng.messages = [{"role": "user", "content": "x" * 400_000}]
    assert eng._should_auto_compact()


def test_disabled_when_pct_zero():
    eng = _bare_engine()
    eng.auto_compact_pct = 0
    eng.messages = [{"role": "user", "content": "x" * 1_000_000}]
    assert not eng._should_auto_compact()


def test_compaction_records_info():
    """Full-compaction path (irreducible token pressure) records last_compaction_info.

    The gentle sliding-window trim (stage 1) only shrinks large *assistant*
    payloads — user messages are preserved verbatim as the GOALS. To force the
    full summarisation path we use a history dominated by large user messages
    that the slide cannot relieve, so usage stays above the budget and full
    compaction fires. (We also need > _COMPACTION_THRESHOLD=12 messages — that
    floor guards against summarising a too-short conversation.)
    """
    eng = _bare_engine()
    big_chunk = "y" * 25_000  # ~6.25k tokens each
    eng.messages = [
        {"role": "user", "content": f"msg-{i} " + big_chunk}
        for i in range(20)
    ]
    # 20 user messages * 25k chars = 500k chars = 125k tokens > 80k budget,
    # and the slide can't trim user content -> full compaction must fire.
    assert eng._should_auto_compact()
    eng._compact_history()
    info = eng.last_compaction_info
    assert info is not None
    assert info["messages_compacted"] >= 1
    assert info["tokens_before"] > info["tokens_after"]
    # Recent messages preserved
    assert len(eng.messages) <= 6  # 2 summary + last 4 recent


def test_no_compaction_for_short_history():
    eng = _bare_engine()
    eng.messages = [
        {"role": "user", "content": "small"},
        {"role": "assistant", "content": "ok"},
    ]
    eng._compact_history()
    assert eng.last_compaction_info is None
    assert len(eng.messages) == 2
