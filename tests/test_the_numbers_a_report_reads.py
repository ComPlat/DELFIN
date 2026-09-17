"""What a diagnostic report can read, and what the numbers mean.

Three sessions built reports on 2026-09-17 and each hit the same three
walls, so the walls are what this pins:

  a trace is readable elsewhere   read(session, root=...) instead of
                                  walking the directory and parsing
                                  the JSONL by hand, which all three did
  one bad field is not fatal      a record with cost_usd "x" used to
                                  take the whole aggregate down
  input tokens are reported       recorded on every turn, summed by
                                  nothing, so a report had to do it
  zero cost means something       an unpriced turn's 0.00 read as free;
                                  now the turns are counted as what they
                                  are
"""

from __future__ import annotations

import json

from delfin.agent import tool_trace
from delfin.agent.agent_metrics import aggregate_by_model


# -- the trace, from anywhere ----------------------------------------------

def _write_trace(root, session: str, entries: list[dict]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{session}.jsonl"
    path.write_text("".join(json.dumps(e) + "\n" for e in entries),
                    encoding="utf-8")


def test_a_trace_can_be_read_where_it_lies(tmp_path):
    root = tmp_path / "traces"
    _write_trace(root, "s1", [{"tool": "bash", "ok": True},
                              {"tool": "read_file", "ok": False}])
    entries = tool_trace.read("s1", root=root)
    assert [e["tool"] for e in entries] == ["bash", "read_file"]


def test_the_sessions_in_a_directory_can_be_listed(tmp_path):
    root = tmp_path / "traces"
    _write_trace(root, "older", [{"tool": "bash"}])
    _write_trace(root, "newer", [{"tool": "bash"}])
    import os
    import time
    os.utime(root / "older.jsonl", (time.time() - 60, time.time() - 60))
    assert tool_trace.sessions(root=root) == ["newer", "older"]


def test_a_directory_that_is_not_there_lists_nothing(tmp_path):
    assert tool_trace.sessions(root=tmp_path / "never-made") == []
    assert tool_trace.read("s1", root=tmp_path / "never-made") == []


def test_the_last_n_still_works_with_a_root(tmp_path):
    root = tmp_path / "traces"
    _write_trace(root, "s1", [{"tool": f"t{i}"} for i in range(5)])
    assert [e["tool"] for e in tool_trace.read("s1", last_n=2, root=root)] \
        == ["t3", "t4"]


# -- the aggregate ----------------------------------------------------------

def _turn(**kw) -> dict:
    row = {"model": "kit.glm-5.3", "ts": 1_700_000_000.0}
    row.update(kw)
    return row


def test_one_unreadable_field_does_not_take_the_aggregate_down():
    rows = [_turn(cost_usd="x", input_tokens=100, output_tokens=10),
            _turn(cost_usd=0.5, input_tokens=200, output_tokens=20)]
    out = aggregate_by_model(rows)["kit.glm-5.3"]
    assert out["n_turns"] == 2
    assert out["total_cost_usd"] == 0.5          # the readable half counts


def test_input_tokens_are_reported_like_output_tokens():
    rows = [_turn(input_tokens=100, output_tokens=10),
            _turn(input_tokens=300, output_tokens=30)]
    out = aggregate_by_model(rows)["kit.glm-5.3"]
    assert out["total_tokens_in"] == 400
    assert out["avg_tokens_in"] == 200
    assert out["total_tokens_out"] == 40


def test_a_zero_that_was_never_priced_is_counted_as_such():
    rows = [_turn(cost_usd=0.0, price_state="unknown"),
            _turn(cost_usd=0.0, price_state="unknown"),
            _turn(cost_usd=0.0, price_state="non_billing"),
            _turn(cost_usd=0.25, price_state="measured")]
    out = aggregate_by_model(rows)["kit.glm-5.3"]
    assert out["total_cost_usd"] == 0.25
    assert out["unpriced_turns"] == 2
    assert out["non_billing_turns"] == 1


def test_an_old_record_without_the_field_counts_as_neither():
    out = aggregate_by_model([_turn(cost_usd=0.0)])["kit.glm-5.3"]
    assert out["unpriced_turns"] == 0
    assert out["non_billing_turns"] == 0


def test_the_turn_record_carries_what_its_cost_means():
    from dataclasses import asdict

    from delfin.agent.agent_metrics import TurnMetrics

    assert "price_state" in asdict(TurnMetrics())
