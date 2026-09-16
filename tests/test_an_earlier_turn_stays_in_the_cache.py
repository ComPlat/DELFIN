"""Elision does not touch the tool rows earlier turns left in the history.

Since 5bb39799 a turn's tool rounds stay in the history (turn_history),
cut to a head. They are the oldest tool rows of every later request, so
the elision that frees context for a long turn reached them first: little
freed, and the prefix every request shares changed at its earliest point,
which is what a15a3838's one-deep-cut rule was meant to prevent.
"""
import inspect

from delfin.agent import api_client as A


def _rows(n_history, n_turn, size=1000):
    rows = []
    for i in range(n_history):
        rows += [{"role": "assistant", "content": "", "tool_calls": [{"id": f"h{i}"}]},
                 {"role": "tool", "tool_call_id": f"h{i}", "content": "H" * size}]
    base = len(rows)
    for i in range(n_turn):
        rows += [{"role": "assistant", "content": "", "tool_calls": [{"id": f"t{i}"}]},
                 {"role": "tool", "tool_call_id": f"t{i}", "content": "T" * size}]
    return rows, base


def test_history_rows_are_left_alone_and_this_turn_pays():
    rows, base = _rows(3, 6)
    n = A._elide_old_tool_results(rows, char_budget=5000, keep_recent=2, protect_before=base)
    assert n >= 1
    assert all(r["content"] == "H" * 1000 for r in rows[:base] if r["role"] == "tool")
    assert any(r["content"].startswith(A._ELIDED_PREFIX) for r in rows[base:] if r["role"] == "tool")


def test_without_a_base_the_oldest_still_go_first():
    rows, base = _rows(3, 6)
    A._elide_old_tool_results(rows, char_budget=5000, keep_recent=2)
    assert rows[1]["content"].startswith(A._ELIDED_PREFIX)


def test_the_loop_hands_over_its_base():
    src = inspect.getsource(A.OpenAIClient.stream_message)
    i = src.index("_elide_old_tool_results(")
    assert 'protect_before=int(getattr(self, "_turn_rows_base", 0) or 0)' in src[i:i + 300]
