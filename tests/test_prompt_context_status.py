"""Tests for the live context-status block injected into the solo prompt.

The block tells the model how close it is to the engine's auto-compaction
trigger (95% of the token window, with a gentle trim from 70%) so it can
proactively delegate to subagents before compaction fires. Compaction is now
purely token-driven — the legacy 12-message trigger was removed because it
fired at ~15% window full and threw away live context (Jerome 2026-06-13).
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
    assert "Compaction trigger: 95% of window" in block
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


# ---------------------------------------------------------------------------
# Open-tasks reminder (Jerome 2026-06-13: agent forgets to check off tasks)
# ---------------------------------------------------------------------------

def _perms_for(tmp_path, sid):
    class _Perms:
        workspace = tmp_path
        task_session_id = sid
    return _Perms()


def test_open_tasks_block_lists_open_and_nudges_completion(tmp_path, monkeypatch):
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    t1 = store.create("Wire parser into ops_server", session_id="s1")
    store.create("Add regression test", session_id="s1")
    store.update(t1["id"], status="in_progress")

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "s1")),
    )
    block = eng._build_open_tasks_block()
    assert "Open tasks" in block
    # Session-relative number for display, global id kept for task_update
    # (bug 172400): "[in_progress] task 1 (id 1) Wire parser…".
    assert "[in_progress] task 1 (id 1)" in block and "Wire parser" in block
    assert "[pending]" in block and "Add regression test" in block
    # An in_progress task present → an explicit completion nudge fires.
    assert "status='completed'" in block


def test_open_tasks_block_empty_when_no_open_tasks(tmp_path, monkeypatch):
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    done = store.create("already finished", session_id="s2")
    store.update(done["id"], status="completed")

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "s2")),
    )
    assert eng._build_open_tasks_block() == ""


def test_open_tasks_block_caps_long_backlog(tmp_path, monkeypatch):
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    for i in range(20):
        store.create(f"task number {i}", session_id="s3")

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "s3")),
    )
    block = eng._build_open_tasks_block()
    # The list is capped; the remainder is surfaced as a count, not dumped.
    assert "more open task(s)" in block
    assert block.count("[pending]") <= 8
