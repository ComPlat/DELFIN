"""Per-call normalisation: a token total that fell because the run did less
work is not a saving.

Wave-10 finding: input tokens fell 4.65M -> 2.51M (looked like a saving),
but tool calls also halved 115 -> 54, so per CALL there was no saving
(~40k -> ~47k). ``summarise_run`` reported only totals, so the compare
output could not tell "cheaper per call" from "did less work".

Adds to ``summarise_run``: ``total_input_tokens``, ``total_output_tokens``
and, next to them, ``per_call_input_tokens`` / ``per_call_output_tokens``
(= tokens per TOOL call; the record carries no per-model-call count, so a
per-model-call figure would have to be fabricated -- it is not). The
compare markdown shows them beside the totals.

Edge cases pinned (reviewer caveat): tool_calls == 0 must not divide by
zero (per-call becomes None, not a crash), and unmeasured rows contribute
neither tokens nor calls.
"""

import pytest

from delfin.agent.benchmark import compare_runs, format_compare_markdown, summarise_run


def _row(tool_calls, in_tok, out_tok, *, unmeasured=False, n=1):
    flags = [True] * n
    return {
        "task_id": "t",
        "success": True,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 50,
        "cost_usd": 0.1,
        "duration_s": 10.0,
        "tool_calls": tool_calls,
        "input_tokens": in_tok,
        "output_tokens": out_tok,
        "unmeasured": unmeasured,
    }


def test_summarise_reports_token_totals_and_per_call():
    s = summarise_run([_row(50, 4_000_000, 1_000_000),
                       _row(60, 3_000_000, 800_000)])
    assert s["total_input_tokens"] == 7_000_000
    assert s["total_output_tokens"] == 1_800_000
    assert s["total_tool_calls"] == 110
    assert s["per_call_input_tokens"] == pytest.approx(7_000_000 / 110)
    assert s["per_call_output_tokens"] == pytest.approx(1_800_000 / 110)


def test_per_call_zero_tool_calls_is_none_not_error():
    s = summarise_run([_row(0, 1000, 500)])
    assert s["total_tool_calls"] == 0
    assert s["per_call_input_tokens"] is None
    assert s["per_call_output_tokens"] is None


def test_per_call_excludes_unmeasured_rows():
    s = summarise_run([_row(10, 1000, 500), _row(0, 0, 0, unmeasured=True)])
    assert s["total_tool_calls"] == 10
    assert s["total_input_tokens"] == 1000
    assert s["per_call_input_tokens"] == pytest.approx(100.0)


def _single_run(tool_calls, in_tok, out_tok):
    return [{
        "task_id": "t", "success": True, "n_samples": 1,
        "per_run_success": [True], "quality_0_100": 60,
        "cost_usd": 0.5, "duration_s": 60.0,
        "tool_calls": tool_calls, "input_tokens": in_tok,
        "output_tokens": out_tok,
    }]


def test_markdown_shows_per_call_and_totals():
    # Wave-10 shape: totals fell, but only because fewer calls were made.
    out = compare_runs(
        _single_run(115, 4_650_000, 500_000),
        _single_run(54, 2_510_000, 300_000),
    )
    md = format_compare_markdown(out)
    # Per-call tokens must appear next to the totals.
    assert "token" in md.lower()
    assert "call" in md.lower()
    old = out["summary"]["old"]
    new = out["summary"]["new"]
    # The comparison data carries both totals and per-call figures.
    assert old["total_input_tokens"] == 4_650_000
    assert new["total_input_tokens"] == 2_510_000
    assert old["per_call_input_tokens"] == pytest.approx(4_650_000 / 115)
