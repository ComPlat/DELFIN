"""A long autonomous turn compacts between its rounds.

The engine compacts only when a user turn starts; an autonomous session
is one turn of hundreds of tool rounds, so it never compacted (four long
sessions on 2026-09-26, zero compactions), and a session with a capped
window ran into context_length_exceeded. The loop now (1) budgets its
elision against the capped window and (2) replaces older rounds by a
deterministic summary when the elision alone does not reach the budget.
"""

from __future__ import annotations

from types import SimpleNamespace

from delfin.agent import api_client


def _msgs(n_rounds, result_chars=4000):
    msgs = [{"role": "user", "content": "the user's goal for this turn"}]
    for i in range(n_rounds):
        msgs.append({"role": "assistant", "content": f"round {i}",
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": "{}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}",
                     "content": "x" * result_chars})
    return msgs


def test_over_budget_the_turn_is_compacted_in_place(monkeypatch):
    logged = []
    monkeypatch.setattr("delfin.agent.compaction_log.record_compaction",
                        lambda sid, **kw: logged.append((sid, kw)))
    msgs = _msgs(30)
    before = len(msgs)
    note = api_client._compact_between_rounds(msgs, 20_000, 0, "sess")
    assert note and "compacted within the turn" in note
    assert len(msgs) < before
    assert msgs[0]["content"] == "the user's goal for this turn"
    # every surviving tool result still answers a surviving tool call
    calls = {tc["id"] for m in msgs for tc in (m.get("tool_calls") or [])}
    answers = {m["tool_call_id"] for m in msgs if m.get("role") == "tool"}
    assert answers <= calls
    assert logged and logged[0][1]["kind"] == "in_turn"


def test_under_budget_nothing_changes():
    msgs = _msgs(3, result_chars=100)
    snapshot = [dict(m) for m in msgs]
    assert api_client._compact_between_rounds(msgs, 100_000, 0, "s") == ""
    assert msgs == snapshot


def test_the_elision_budget_follows_the_capped_window(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    caps = SimpleNamespace(context_window=131_072)
    assert api_client._tool_context_char_budget(caps) == int(40_000 * 0.45 * 4)
