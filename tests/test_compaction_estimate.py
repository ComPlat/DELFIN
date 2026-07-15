"""Compaction budget accounts for the whole request, and the rebuilt history
stays role-alternating.

The estimate used to sum only self.messages, ignoring the system prompt (repo
map, memory, playbook — routinely 10-40k tokens) and the tool schemas, so a
request could pass the auto-compact budget yet overflow the provider window
(400). And the post-compaction rebuild could emit two consecutive assistant
turns, which strict local chat templates reject.
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


# --- B1: the estimate reflects the whole request ---------------------------

def test_estimate_folds_in_system_prompt():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "hi"}]
    base = eng._estimate_context_tokens()
    eng.last_system_prompt = "S" * 8000          # ~2000 tokens
    assert eng._estimate_context_tokens() == base + 2000


def test_estimate_never_below_real_input_count():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "x" * 4000}]   # ~1000 tokens
    # The provider reported a much larger real count (system + schemas + cache);
    # the estimate must not fall below it.
    eng._last_input_tokens = 50_000
    assert eng._estimate_context_tokens() == 50_000


def test_estimate_counts_image_parts():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": [
        {"type": "image_url", "image_url": {"url": "data:...."}}]}]
    # An image must not count as ~0 tokens.
    assert eng._estimate_context_tokens() >= 1000


# --- B2: compaction leaves an alternating, sendable history -----------------

def test_compaction_yields_alternating_roles():
    eng = _bare_engine()
    eng.context_window_tokens = 2000
    eng.auto_compact_pct = 0.5
    eng._llm_summarize_old_messages = lambda old: "SUMMARY"
    msgs = []
    for _ in range(20):
        msgs.append({"role": "user", "content": "u" * 2000})
        msgs.append({"role": "assistant", "content": "a" * 2000})
    # Mimic _compact_history's call site: it runs right after a user turn is
    # appended, so the list ends on a user message and `recent` starts assistant.
    msgs.append({"role": "user", "content": "latest question"})
    eng.messages = msgs

    assert eng._should_auto_compact()
    eng._compact_history()

    roles = [m["role"] for m in eng.messages]
    for a, b in zip(roles, roles[1:]):
        assert a != b, f"consecutive {a!r} messages after compaction: {roles}"
    # The summary block is present and the recent tail survived.
    assert eng.messages[0]["content"].startswith("[Conversation summary")
    assert eng.messages[-1]["content"] == "latest question"


def test_compaction_resets_stale_input_floor():
    eng = _bare_engine()
    eng.context_window_tokens = 2000
    eng.auto_compact_pct = 0.5
    eng._last_input_tokens = 999_999          # stale, pre-compaction
    eng._llm_summarize_old_messages = lambda old: "SUMMARY"
    msgs = []
    for _ in range(20):
        msgs.append({"role": "user", "content": "u" * 2000})
        msgs.append({"role": "assistant", "content": "a" * 2000})
    msgs.append({"role": "user", "content": "latest"})
    eng.messages = msgs
    eng._compact_history()
    # The stale floor must be cleared so the estimate reflects the small
    # compacted size (else the next turn over-compacts).
    assert eng._last_input_tokens == 0
