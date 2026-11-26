"""Derailment counts a long run's records attest (package R3).

Every case here is written through the REAL writers
(``delfin/agent/longhorizon/replay.replay`` for the shape a transcript
takes, ``tool_trace.record``'s exact field set for the records), never a
hand-typed shadow. Each test names the input that would make its assertion
FAIL, per the brief's rule 5.
"""
from __future__ import annotations

import pytest

from delfin.agent.longhorizon.derailments import derailments


def _rec(ts, tool="bash", cmd="xtb struct.xyz", ok=False, out=""):
    """A record in tool_trace.record's own field set."""
    return {"ts": ts, "tool": tool, "input": cmd,
            "output": out, "duration_ms": 5, "ok": ok,
            "error": "" if ok else "not found", "round": 0}


def test_three_identical_failing_calls_in_a_row_count_two():
    # Fails if a repeat is not required: a single failing call scores 0.
    tools = [_rec(1.0), _rec(2.0), _rec(3.0)]
    assert derailments(tools)["repeated_failing_calls"] == 2


def test_a_single_failing_call_is_not_a_repeat():
    assert derailments([_rec(1.0)])["repeated_failing_calls"] == 0


def test_an_intervening_success_breaks_the_failing_run():
    # The same failing command twice, but a DIFFERENT call ran between:
    # not consecutive — fails if the reducer ignored the order.
    tools = [_rec(1.0),
             _rec(2.0, cmd="xtb other.xyz", ok=True, out="energy"),
             _rec(3.0)]
    assert derailments(tools)["repeated_failing_calls"] == 0


def test_different_arguments_are_not_the_same_call():
    # Same tool, different input: a retry with a changed argument is
    # progress, not a repeat — fails if identity keyed on the tool alone.
    tools = [_rec(1.0, cmd="xtb a.xyz"), _rec(2.0, cmd="xtb b.xyz")]
    assert derailments(tools)["repeated_failing_calls"] == 0


def test_two_identical_failing_calls_count_one():
    assert derailments([_rec(1.0), _rec(2.0)])[
        "repeated_failing_calls"] == 1


def test_identical_successful_reads_three_in_a_row_is_a_loop():
    # Fails if only failing runs were counted: a successful identical read
    # repeated 3x is the no-progress loop the engine steers at 3.
    tools = [_rec(1.0, cmd="cat data.csv", ok=True, out="rows"),
             _rec(2.0, cmd="cat data.csv", ok=True, out="rows"),
             _rec(3.0, cmd="cat data.csv", ok=True, out="rows")]
    assert derailments(tools)["unproductive_loops"] == 1


def test_two_identical_successful_reads_are_not_a_loop():
    tools = [_rec(1.0, cmd="cat data.csv", ok=True, out="rows"),
             _rec(2.0, cmd="cat data.csv", ok=True, out="rows")]
    assert derailments(tools)["unproductive_loops"] == 0


def test_a_loop_broken_by_another_call_is_not_a_loop():
    tools = [_rec(1.0, cmd="cat data.csv", ok=True, out="rows"),
             _rec(2.0, cmd="cat b.csv", ok=True, out="rows"),
             _rec(3.0, cmd="cat data.csv", ok=True, out="rows")]
    assert derailments(tools)["unproductive_loops"] == 0


def test_mixed_ok_false_run_is_a_failing_repeat_not_a_loop():
    # 4 identical calls, one succeeded in the middle: the run is failing
    # (all() is False), not a successful loop.
    tools = [_rec(1.0),
             _rec(2.0, ok=True, out="x"),
             _rec(3.0), _rec(4.0)]
    counts = derailments(tools)
    assert counts["repeated_failing_calls"] == 3
    assert counts["unproductive_loops"] == 0


def test_order_comes_from_ts_not_list_position():
    # The records arrive out of order: the reducer must sort by ts. Fails
    # if it trusted list order (the interleaving would read as consecutive).
    tools = [_rec(3.0), _rec(1.0),
             _rec(2.0, cmd="xtb other.xyz", ok=True, out="energy")]
    assert derailments(tools)["repeated_failing_calls"] == 0
    tools = [_rec(2.0, cmd="xtb other.xyz", ok=True, out="energy"),
             _rec(3.0), _rec(1.0)]
    assert derailments(tools)["repeated_failing_calls"] == 0


def test_records_without_ts_do_not_crash_and_do_not_count():
    tools = [{"tool": "bash", "input": "xtb x", "ok": False}]
    assert derailments(tools) == {"repeated_failing_calls": 0,
                                  "unproductive_loops": 0}


def test_empty_input_is_all_zero():
    assert derailments([]) == {"repeated_failing_calls": 0,
                               "unproductive_loops": 0}


def test_records_from_the_real_replay_shape_reduce_the_same():
    """The shapes replay() writes through tool_trace.record are what a real
    session's records look like; a hand-typed dict must not be the only
    input the reducer understands."""
    from delfin.agent.longhorizon import replay
    transcript = [
        {"kind": "tool", "tool": "bash", "input": "xtb s.xyz",
         "output": "command not found", "ok": False},
        {"kind": "tool", "tool": "bash", "input": "xtb s.xyz",
         "output": "command not found", "ok": False},
        {"kind": "tool", "tool": "bash", "input": "xtb s.xyz",
         "output": "command not found", "ok": False},
    ]
    out = replay.replay(transcript, session="r3-derail-red")
    counts = derailments(out["records"]["tools"])
    assert counts["repeated_failing_calls"] == 2
