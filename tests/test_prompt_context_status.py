"""Tests for the live context-status block injected into the solo prompt.

The block tells the model how close it is to the engine's auto-compaction
trigger (12 msgs OR 95% of 100k tokens) so it can proactively delegate to
subagents before compaction fires.
"""

from __future__ import annotations

from delfin.agent.engine import AgentEngine


def _bare_engine() -> AgentEngine:
    """Build a minimal engine for unit-testing _build_context_status_block
    without touching any backend, provider, or session state."""
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    return eng


def test_status_block_zero_usage_when_no_messages():
    eng = _bare_engine()
    block = eng._build_context_status_block()
    assert "Compaction trigger: 12 msgs OR 95% of window" in block
    assert "0 msgs" in block
    assert "(none this session)" in block
    # Empty conversation must NOT trigger the >=80% warning
    assert "WARNING" not in block


def test_status_block_warns_near_threshold():
    eng = _bare_engine()
    # The estimator counts ~chars/4. Fill enough to clear 80% of the 100k
    # token window: 320 000 chars = ~80 000 tokens.
    eng.messages = [
        {"role": "user", "content": "x" * 40_000} for _ in range(8)
    ]
    block = eng._build_context_status_block()
    assert "WARNING: nearing auto-compaction" in block, block


def test_status_block_reflects_last_compaction():
    eng = _bare_engine()
    eng.last_compaction_info = {
        "messages_compacted": 9,
        "tokens_before": 95_000,
        "tokens_after": 8_000,
        "tokens_saved": 87_000,
    }
    block = eng._build_context_status_block()
    assert "compacted 9" in block
    assert "87,000" in block or "87000" in block


def test_status_block_only_injected_for_solo_role(monkeypatch):
    """Building the system prompt for solo_agent must include the block;
    other roles must not be polluted with it."""
    eng = _bare_engine()
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    eng.mode = "solo"
    eng.mode_description = ""
    eng.role_outputs = {}
    eng._live_state = ""
    eng._prompt_session_serial = 0
    # Avoid touching session error context in this isolated test
    monkeypatch.setattr(eng, "format_error_context", lambda: "")

    captured = {}

    class _FakeLoader:
        def build_system_prompt(self, **kwargs):
            captured.update(kwargs)
            return "OK"

    eng.loader = _FakeLoader()
    eng._build_current_system_prompt()
    assert "# Context status" in (captured.get("live_state") or ""), captured

    # Same engine, different role — no status block expected
    captured.clear()
    eng.route = ["builder_agent"]
    eng._build_current_system_prompt()
    assert "# Context status" not in (captured.get("live_state") or "")
