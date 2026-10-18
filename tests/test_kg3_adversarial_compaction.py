"""Adversarial tests for the in-turn compaction contract (KG3 review).

These pin the invariants that any keep-list / quality change (KG3 phase 2)
must not break:
1. Compaction still happens when the budget is hopeless — the
   context_length cliff outranks keeping anything.
2. The last round's open tool-call pair is never torn.
3. Same input, same output (determinism contract).
4. A person's goal outside the compacted rounds stays verbatim.
5. The summary never grows past its cap.
"""

import pytest

from delfin.agent.in_turn_compaction import compact


def _turn(n_rounds, chars_per_result=3000):
    """A synthetic single turn: `n_rounds` complete tool-call rounds."""
    msgs = [{"role": "user", "content": "Goal: measure the elision rate."}]
    for i in range(n_rounds):
        cid = f"call_{i}"
        msgs.append({
            "role": "assistant",
            "tool_calls": [{"id": cid, "type": "function",
                            "function": {"name": "bash", "arguments": "{}"}}],
        })
        msgs.append({
            "role": "tool", "tool_call_id": cid,
            "content": "x" * chars_per_result,
        })
    return msgs


def test_hopeless_budget_still_compacts():
    """Even when no budget can be met, the oldest rounds must go —
    the cliff outranks keeping."""
    msgs = _turn(8, chars_per_result=4000)
    out, record = compact(msgs, budget_chars=1, keep_last_rounds=2)
    assert record is not None
    assert record["rounds_compacted"] >= 1
    assert out != msgs


def test_last_round_is_never_torn():
    """The newest round survives whole: its call has its result."""
    msgs = _turn(6)
    out, _record = compact(msgs, budget_chars=1000, keep_last_rounds=2)
    last_call = out[-2]
    last_result = out[-1]
    assert last_result["role"] == "tool"
    assert (last_result["tool_call_id"]
            == last_call["tool_calls"][0]["id"])


def test_determinism_same_input_same_output():
    msgs = _turn(7)
    a, ra = compact(msgs, budget_chars=5000, keep_last_rounds=2)
    b, rb = compact(msgs, budget_chars=5000, keep_last_rounds=2)
    assert a == b
    assert ra == rb


def test_goal_outside_rounds_stays_verbatim():
    """The user's goal row sits before the compacted rounds and must not
    be rewritten or dropped."""
    goal = "Goal: measure the elision rate over three sessions."
    msgs = _turn(5)
    msgs[0]["content"] = goal
    out, _record = compact(msgs, budget_chars=1000, keep_last_rounds=2)
    assert any(m.get("role") == "user" and m.get("content") == goal
               for m in out)


def test_summary_is_capped():
    """The summary text never exceeds the module's own cap, even with a
    huge number of dropped results."""
    msgs = _turn(30, chars_per_result=9000)
    out, record = compact(msgs, budget_chars=500, keep_last_rounds=2)
    assert record is not None
    summaries = [m for m in out if m.get("role") == "user"
                 and str(m.get("content", "")).startswith("[Conversation summary")]
    assert len(summaries) == 1
    text = summaries[0]["content"]
    assert len(text) <= 100000  # module cap; far above the cap is a bug
