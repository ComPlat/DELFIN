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
    assert "_elide_old_tool_results(api_messages)" in src[i:i + 600]
