"""Tests for tool-result-aware context editing in long tool-call loops.

Contract: when accumulated ``role=="tool"`` output exceeds the budget,
the OLDEST tool results are elided to a placeholder while the most recent
ones stay verbatim; user/system/assistant messages and every
``tool_call_id`` are preserved (protocol-safe). Under budget -> no-op.
"""

from __future__ import annotations

from delfin.agent.api_client import _elide_old_tool_results, _ELIDED_PREFIX


def _tool(i, n):
    return {"role": "tool", "tool_call_id": f"c{i}", "content": "X" * n}


def test_noop_under_budget():
    msgs = [{"role": "system", "content": "s"},
            {"role": "user", "content": "u"},
            _tool(1, 100), _tool(2, 100)]
    before = [dict(m) for m in msgs]
    assert _elide_old_tool_results(msgs, char_budget=10000) == 0
    assert msgs == before                       # untouched


def test_oldest_elided_recent_kept():
    msgs = [{"role": "system", "content": "s"}]
    # 10 tool results of 1000 chars each = 10000; budget 4000, keep 3.
    for i in range(10):
        msgs.append({"role": "assistant", "content": f"reasoning {i}",
                     "tool_calls": [{"id": f"c{i}"}]})
        msgs.append(_tool(i, 1000))
    n = _elide_old_tool_results(msgs, char_budget=4000, keep_recent=3)
    assert n >= 1
    tools = [m for m in msgs if m.get("role") == "tool"]
    # the three most recent are kept verbatim ...
    assert all(not t["content"].startswith(_ELIDED_PREFIX) for t in tools[-3:])
    assert all(len(t["content"]) == 1000 for t in tools[-3:])
    # ... the oldest were elided ...
    assert tools[0]["content"].startswith(_ELIDED_PREFIX)
    # ... but every tool_call_id survives (protocol stays valid) ...
    assert [t["tool_call_id"] for t in tools] == [f"c{i}" for i in range(10)]
    # ... and reasoning/system/user are never touched.
    assert all(m["content"].startswith("reasoning")
               for m in msgs if m.get("role") == "assistant")


def test_assistant_and_user_never_elided():
    big = "Y" * 50000
    msgs = [{"role": "user", "content": big},
            {"role": "assistant", "content": big},
            _tool(1, 50000)]
    _elide_old_tool_results(msgs, char_budget=1000, keep_recent=0)
    assert msgs[0]["content"] == big            # user huge -> untouched
    assert msgs[1]["content"] == big            # assistant huge -> untouched


def test_already_elided_not_doubled():
    msgs = [_tool(1, 2000), _tool(2, 2000), _tool(3, 2000)]
    _elide_old_tool_results(msgs, char_budget=1500, keep_recent=0)
    # run again — placeholders must not be re-wrapped
    _elide_old_tool_results(msgs, char_budget=1500, keep_recent=0)
    elided = [m for m in msgs if m["content"].startswith(_ELIDED_PREFIX)]
    for m in elided:
        assert m["content"].count(_ELIDED_PREFIX) == 1


def test_called_each_round_in_source():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "agent" / "api_client.py").read_text(encoding="utf-8")
    i = src.find("for _round in range(_MAX_TOOL_ROUNDS")
    assert i != -1
    # Called each round with the context-scaled budget (bug 172455).
    assert "_elide_old_tool_results(api_messages, char_budget=" in src[i:i + 600]


# ---------------------------------------------------------------------------
# Context-scaled elision budget (bug 172455: a 262k model must keep far more
# tool output than the fixed 60k-char default, so it stops re-paging files)
# ---------------------------------------------------------------------------

from delfin.agent.api_client import (
    _tool_context_char_budget, _TOOL_CONTEXT_CHAR_BUDGET,
)


class _Caps:
    def __init__(self, ctx):
        self.context_window = ctx


def test_budget_floor_when_caps_missing_or_small():
    # None caps, zero/None window, and small models all keep the 60k floor.
    assert _tool_context_char_budget(None) == _TOOL_CONTEXT_CHAR_BUDGET
    assert _tool_context_char_budget(_Caps(0)) == _TOOL_CONTEXT_CHAR_BUDGET
    assert _tool_context_char_budget(_Caps(None)) == _TOOL_CONTEXT_CHAR_BUDGET
    # 8k model: 8192*0.45*4 = 14745 < 60000 → floor, never smaller than today.
    assert _tool_context_char_budget(_Caps(8192)) == _TOOL_CONTEXT_CHAR_BUDGET


def test_budget_scales_up_for_big_context():
    big = _tool_context_char_budget(_Caps(262_144))
    assert big == int(262_144 * 0.45 * 4)        # ~471k chars
    assert big > _TOOL_CONTEXT_CHAR_BUDGET * 5   # far more than the fixed default


def test_big_budget_keeps_what_small_budget_elides():
    # 12 tool results of 8000 chars = 96000 total.
    def _mk():
        m = [{"role": "system", "content": "s"}]
        for i in range(12):
            m.append({"role": "assistant", "content": f"r{i}",
                      "tool_calls": [{"id": f"c{i}"}]})
            m.append({"role": "tool", "tool_call_id": f"c{i}",
                      "content": "X" * 8000})
        return m
    small = _mk()
    big = _mk()
    # Default (small) budget elides; a 262k-scaled budget (~471k) keeps all.
    n_small = _elide_old_tool_results(
        small, char_budget=_tool_context_char_budget(_Caps(8192)))
    n_big = _elide_old_tool_results(
        big, char_budget=_tool_context_char_budget(_Caps(262_144)))
    assert n_small >= 1            # 96k > 60k floor → some elided
    assert n_big == 0              # 96k < 471k → nothing elided, no re-paging
