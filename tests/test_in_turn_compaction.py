"""Deterministic compaction BETWEEN two rounds of one tool-call turn.

Control (red before delfin/agent/in_turn_compaction.py existed): nothing
compacted inside a turn — the loop's context grew past the capped window
with no record and no summary. The public shape is the one the wiring
will call:

    compact(messages, budget_chars, *, keep_last_rounds=2,
            turn_rows_base=0) -> (new_messages, record)
"""

from __future__ import annotations

from delfin.agent.in_turn_compaction import compact


def _msgs(n_rounds: int, result_chars: int = 4000) -> list[dict]:
    """A plausible in-turn history: prefix (earlier turns), then rounds
    of assistant tool_calls + tool results, all with distinct ids."""
    msgs = [
        {"role": "user", "content": "goal from an earlier turn"},
        {"role": "assistant", "content": "ack"},
    ]
    for i in range(n_rounds):
        msgs.append({
            "role": "assistant",
            "content": f"round {i} reasoning",
            "tool_calls": [{
                "id": f"call-{i}",
                "type": "function",
                "function": {"name": "bash", "arguments":
                             "{\"command\": \"ls\"}"},
            }],
        })
        msgs.append({
            "role": "tool",
            "tool_call_id": f"call-{i}",
            "content": f"[Command results]\nr{i} " + "x" * result_chars,
        })
    return msgs


class TestPairIntegrity:
    def test_no_tool_result_without_its_call(self):
        msgs = _msgs(6)
        new, _rec = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        calls = set()
        for m in new:
            for tc in m.get("tool_calls") or []:
                calls.add(tc.get("id"))
        for m in new:
            if m.get("role") == "tool":
                assert m.get("tool_call_id") in calls, (
                    "a tool result survived whose call was compacted away")

    def test_kept_rounds_stay_verbatim(self):
        # Intent: the last ``keep_last_rounds`` rounds survive untouched
        # WHEN THE BUDGET ALLOWS them (the original 3000-char budget was
        # smaller than the two protected rounds themselves — an
        # impossible demand the budget rule must win, see
        # test_budget_beats_the_keep_preference).
        msgs = _msgs(6)
        new, _rec = compact(msgs, budget_chars=9_000, keep_last_rounds=2)
        kept = [m["content"] for m in new if m.get("role") == "tool"][-2:]
        assert any("r5" in c for c in kept)
        assert any("r4" in c for c in kept)

    def test_budget_beats_the_keep_preference(self):
        # When keeping the preferred rounds cannot reach the budget,
        # older "protected" rounds go too — the turn survives; only the
        # LAST round (the open tool-call pair) is untouchable.
        msgs = _msgs(6)  # 4 rounds ~16k tool chars, budget 5k
        new, rec = compact(msgs, budget_chars=5_000, keep_last_rounds=2)
        tool_chars = sum(len(str(m.get("content", "")))
                         for m in new if m.get("role") == "tool")
        assert tool_chars <= 5_000
        # the very last round's pair is intact
        calls = {tc.get("id") for m in new
                 for tc in (m.get("tool_calls") or [])}
        last_tools = [m for m in new if m.get("role") == "tool"]
        assert last_tools[-1].get("tool_call_id") in calls


class TestBudget:
    def test_result_is_under_budget(self):
        msgs = _msgs(8, result_chars=5000)
        new, rec = compact(msgs, budget_chars=10_000, keep_last_rounds=2)
        tool_chars = sum(len(str(m.get("content", "")))
                         for m in new if m.get("role") == "tool")
        assert tool_chars <= 10_000
        assert rec["messages_compacted"] >= 1

    def test_under_budget_is_a_noop(self):
        msgs = _msgs(3, result_chars=500)
        new, rec = compact(msgs, budget_chars=1_000_000,
                           keep_last_rounds=2)
        assert new == msgs
        assert rec is None


class TestProtectedContent:
    def test_user_goals_survive(self):
        msgs = _msgs(6)
        new, _rec = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        assert any(m.get("role") == "user"
                   and "goal from an earlier turn" in str(m.get("content"))
                   for m in new)

    def test_prefix_before_turn_rows_base_is_untouched(self):
        msgs = _msgs(6)
        new, _rec = compact(msgs, budget_chars=3000, keep_last_rounds=2,
                            turn_rows_base=2)
        assert new[0] == msgs[0]
        assert new[1] == msgs[1]


class TestWhatSurvives:
    def test_summary_names_the_working_state(self):
        msgs = _msgs(6)
        # a fact the working-state extractors find: a test outcome
        msgs.append({"role": "user", "content":
                     "[Command results]\ngate tests/test_q.py -> 3 passed"})
        new, rec = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        joined = "\n".join(str(m.get("content", "")) for m in new)
        assert "3 passed" in joined, (
            "the summary must carry the facts the working state extracts")

    def test_record_counts_and_is_deterministic(self):
        msgs = _msgs(6)
        a_new, a_rec = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        b_new, b_rec = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        assert a_new == b_new
        # the record's stable parts (no timestamps)
        assert a_rec["messages_compacted"] == b_rec["messages_compacted"]
        assert a_rec["kind"] == "in_turn"

    def test_idempotent(self):
        msgs = _msgs(6)
        once, _ = compact(msgs, budget_chars=3000, keep_last_rounds=2)
        twice, rec2 = compact(once, budget_chars=3000, keep_last_rounds=2)
        assert twice == once
        assert rec2 is None
