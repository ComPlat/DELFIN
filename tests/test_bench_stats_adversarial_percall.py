"""Adversarial tests for per-call token normalisation (package D, reviewer s18).

Target: d7727763 — summarise_run's total_input_tokens / total_output_tokens
and per_call_input_tokens / per_call_output_tokens, and their rendering in
format_compare_markdown.  The goal is to break the normalisation: totals
over the WRONG row set, per-call over the wrong denominator, aggregate
rows vs raw rows, None/zero edge shapes, and the markdown's n/a path.
"""

import pytest

from delfin.agent.benchmark import compare_runs, format_compare_markdown, summarise_run


def _row(task_id, tool_calls, in_tok, out_tok, **over):
    n = over.pop("n", 1)
    flags = over.pop("per_run_success", None) or [True] * n
    d = {
        "task_id": task_id,
        "success": True,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 50,
        "cost_usd": 0.1,
        "duration_s": 10.0,
        "tool_calls": tool_calls,
        "input_tokens": in_tok,
        "output_tokens": out_tok,
    }
    d.update(over)
    return d


# ── A. totals over the wrong row set ──────────────────────────────────────

def test_token_totals_exclude_unmeasured_but_cost_includes_all():
    # The empty-rows branch and the scored branch must agree on one rule.
    # Unmeasured rows never ran: their tokens AND calls must be excluded,
    # while cost (real spend, even of a failed attempt) stays summed over
    # all rows.  The current code sums tokens over `scored` and cost over
    # `rows` — verify the split is deliberate and not accidental.
    s = summarise_run([
        _row("t1", 10, 1000, 500),
        _row("t2", 0, 0, 0, unmeasured=True, cost_usd=0.3),
    ])
    assert s["total_input_tokens"] == 1000
    assert s["total_output_tokens"] == 500
    assert s["total_cost_usd"] == pytest.approx(0.4)
    # The tokens must not include the unmeasured row's (zero) tokens only
    # because they HAPPEN to be zero: use a non-zero token count on it.
    s2 = summarise_run([
        _row("t1", 10, 1000, 500),
        _row("t2", 4, 777, 333, unmeasured=True),
    ])
    assert s2["total_input_tokens"] == 1000
    assert s2["per_call_input_tokens"] == pytest.approx(100.0)


def test_total_tool_calls_counts_unmeasured_rows_but_tokens_do_not():
    # Inconsistency probe: total_tool_calls sums over ALL rows (pre-existing
    # `for r in rows`) while the token sums now run over `scored`.  An
    # unmeasured row with recorded tool calls would be counted in the
    # denominator of nothing but present in the tool-call total -- the two
    # numbers on the same line then disagree about what a "call" is.
    s = summarise_run([
        _row("t1", 10, 1000, 500),
        _row("t2", 4, 777, 333, unmeasured=True),
    ])
    assert s["total_tool_calls"] == 10, (
        "total_tool_calls must cover the same rows the token totals cover")


# ── B. per-call over the wrong denominator ────────────────────────────────

def test_per_call_uses_scored_calls_not_total_calls():
    # If total_tool_calls counts unmeasured rows but the token sum does
    # not, a per-call division over ALL calls would understate tokens/call.
    s = summarise_run([
        _row("t1", 10, 1000, 500),
        _row("t2", 4, 777, 333, unmeasured=True),
    ])
    assert s["per_call_input_tokens"] == pytest.approx(100.0)


def test_per_call_none_when_every_row_is_unmeasured():
    s = summarise_run([_row("t1", 5, 999, 999, unmeasured=True)])
    assert s["total_input_tokens"] == 0
    assert s["per_call_input_tokens"] is None
    assert s["per_call_output_tokens"] is None


def test_per_call_exact_with_mixed_zero_call_rows():
    # A scored row with zero tool calls still contributes tokens; the
    # denominator is the SUM of calls, not the count of tasks.
    s = summarise_run([
        _row("t1", 3, 600, 300),
        _row("t2", 0, 400, 100),
    ])
    assert s["per_call_input_tokens"] == pytest.approx(1000 / 3)
    assert s["per_call_output_tokens"] == pytest.approx(400 / 3)


def test_per_call_survives_negative_token_rows():
    # input_tokens is clamped >= 0 upstream (cli.py max(0, ...)), but a
    # corrupted run file with a negative value must not produce a negative
    # per-call figure silently.
    s = summarise_run([_row("t1", 2, -500, 100)])
    assert s["per_call_input_tokens"] >= 0 or s["per_call_input_tokens"] is None


# ── C. markdown rendering of the per-call table ───────────────────────────

def _single_run(tool_calls, in_tok, out_tok):
    return [_row("t", tool_calls, in_tok, out_tok)]


def test_markdown_shows_n_a_when_no_calls_either_side():
    out = compare_runs(_single_run(0, 500, 200), _single_run(0, 600, 250))
    md = format_compare_markdown(out)
    assert "n/a" in md


def test_markdown_delta_is_n_a_when_one_side_has_no_calls():
    out = compare_runs(_single_run(0, 500, 200), _single_run(54, 2_500, 900))
    md = format_compare_markdown(out)
    assert "n/a" in md
    # The side WITH calls must still show its per-call figure.
    assert f"{2_500 / 54:,.1f}" in md or f"{2500/54:.1f}" in md


def test_markdown_token_delta_signed_correctly_for_a_real_saving():
    out = compare_runs(_single_run(115, 4_650_000, 500_000),
                       _single_run(54, 2_510_000, 300_000))
    md = format_compare_markdown(out)
    total_line = [ln for ln in md.splitlines() if "In tokens" in ln][0]
    assert "-2,140,000" in total_line.replace(" ", "") or "-2140000" in total_line


def test_summarise_empty_run_has_none_per_call_and_zero_totals():
    s = summarise_run([])
    assert s["total_input_tokens"] == 0
    assert s["total_output_tokens"] == 0
    assert s["per_call_input_tokens"] is None
    assert s["per_call_output_tokens"] is None


def test_markdown_per_call_values_are_exact_not_median_of_medians():
    # Aggregated rows carry MEDIAN tokens and MEDIAN tool calls, so a
    # per-call figure computed from them is a ratio of medians, not a
    # mean per-call cost.  Pin what it IS so nobody reads it as a mean.
    rows = [
        _row("t1", 10, 1000, 500, n=3, per_run_success=[True, True, True]),
        _row("t2", 2, 100, 50),
    ]
    s = summarise_run(rows)
    assert s["total_input_tokens"] == 1100
    assert s["total_tool_calls"] == 12
    assert s["per_call_input_tokens"] == pytest.approx(1100 / 12)
