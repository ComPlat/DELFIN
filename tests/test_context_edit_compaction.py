"""Deterministic context-editing before the LLM summary (robustness #3).

DELFIN already trims the oldest tool output at 0.70 (sliding window). This adds
a structure-preserving HARD-CLEAR of old bulky tool output at the 0.95 cliff
BEFORE any LLM summary, so the lossy + non-deterministic + paid summary is a
genuine last resort. User goals and assistant reasoning are never cleared.
"""

from __future__ import annotations

from unittest.mock import MagicMock

from delfin.agent.engine import AgentEngine


def _bare_engine() -> AgentEngine:
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


# --- unit: the deterministic clear itself -----------------------------------

def test_hard_clear_stubs_large_nonuser_messages_only():
    eng = _bare_engine()
    eng.context_window_tokens = 100      # budget = 80 tokens = 320 chars
    eng.auto_compact_pct = 0.80
    eng.messages = [
        {"role": "user", "content": "GOAL keep me verbatim"},
        {"role": "assistant", "content": "A" * 4000},
        {"role": "tool", "content": "B" * 4000},
        {"role": "assistant", "content": "short reply"},   # < 800 → untouched
    ]
    n = eng._hard_clear_old_tool_results(eng.messages)
    assert n == 2
    assert eng.messages[0]["content"] == "GOAL keep me verbatim"   # user goal safe
    assert eng.messages[1]["content"].startswith("[cleared:")
    assert "AAAA" in eng.messages[1]["content"]                    # keeps a head
    assert eng.messages[2]["content"].startswith("[cleared:")
    assert eng.messages[3]["content"] == "short reply"             # small untouched


def test_hard_clear_is_idempotent():
    eng = _bare_engine()
    eng.context_window_tokens = 100
    eng.messages = [{"role": "assistant", "content": "A" * 4000}]
    eng._hard_clear_old_tool_results(eng.messages)
    stub = eng.messages[0]["content"]
    # Second pass must not re-clear an already-cleared stub.
    assert eng._hard_clear_old_tool_results(eng.messages) == 0
    assert eng.messages[0]["content"] == stub


# --- integration: summary is a LAST RESORT ----------------------------------

def test_compact_uses_context_edit_and_skips_summary():
    eng = _bare_engine()
    eng.context_window_tokens = 2500     # budget 2000 tokens; slide floor keeps
    eng.auto_compact_pct = 0.80          # us above 0.8 so hard-clear must run
    called = {"summary": False}
    eng._llm_summarize_old_messages = lambda old: called.__setitem__("summary", True) or "SUMMARY"
    # 20 large OLD messages + 4 small RECENT (kept intact) — clearing old frees
    # enough that the summary is not needed.
    msgs = [{"role": "assistant", "content": "x" * 4000} for _ in range(20)]
    msgs += [{"role": "assistant", "content": "recent"} for _ in range(4)]
    eng.messages = msgs
    assert eng._should_auto_compact()
    eng._compact_history()
    assert called["summary"] is False
    assert eng.last_compaction_info.get("kind") == "context_edit"


def test_compact_falls_back_to_summary_when_clear_insufficient():
    eng = _bare_engine()
    eng.context_window_tokens = 2500
    eng.auto_compact_pct = 0.80
    called = {"summary": False}
    eng._llm_summarize_old_messages = lambda old: called.__setitem__("summary", True) or "SUMMARY"
    # RECENT four are themselves huge → clearing OLD cannot get under budget.
    msgs = [{"role": "assistant", "content": "x" * 4000} for _ in range(20)]
    msgs += [{"role": "assistant", "content": "y" * 8000} for _ in range(4)]
    eng.messages = msgs
    assert eng._should_auto_compact()
    eng._compact_history()
    assert called["summary"] is True
