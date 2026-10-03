"""A compare must say when N is too small to decide AND what N would be
needed -- a verdict alone, however honest, is not a next step.

At N=3 even a 3-of-3 vs 0-of-3 sweep has two-sided Fisher p = 0.10, so the
label is always "within noise" and a reader is left guessing whether more
repeats could ever settle it. This adds ``min_n_for_delta`` (the smallest
per-arm N that resolves a given two-sided gap at a stated alpha) and the
compare markdown prints it whenever the pooled result is not significant.

Known values (two-sided Fisher on the extreme 100%-vs-0% table
[[n,0],[0,n]] -> p = 2/C(2n,n)):
  alpha .05: n=3 -> 2/20 = 0.10 (not significant); n=4 -> 2/70 ~ 0.0286
             (significant)  => min_n_for_delta(1.0, 0.0) = 4
  alpha .01: n=4 -> 0.0286 (not); n=5 -> 2/252 ~ 0.0079 (significant)
             => min_n_for_delta(1.0, 0.0, alpha=0.01) = 5
  identical rates cannot be resolved  => None
"""

import pytest

from delfin.agent.benchmark import compare_runs, format_compare_markdown, min_n_for_delta


# ── pure helper: smallest per-arm N that resolves a gap ────────────────────

def test_min_n_full_swing_alpha_05():
    # 3/3 vs 0/3 cannot be called significant; 4/4 vs 0/4 can.
    assert min_n_for_delta(1.0, 0.0) == 4


def test_min_n_stricter_alpha_needs_more_repeats():
    assert min_n_for_delta(1.0, 0.0, alpha=0.01) == 5
    assert min_n_for_delta(1.0, 0.0, alpha=0.01) > min_n_for_delta(1.0, 0.0)


def test_min_n_identical_rates_are_impossible():
    # A zero gap is zero at every N -- no repeat count separates them.
    assert min_n_for_delta(0.5, 0.5) is None
    # Regressive arguments (low > high) are the same non-answer, not a crash.
    assert min_n_for_delta(0.0, 1.0) is None


def test_min_n_unresolvable_at_cap_returns_none():
    # A tiny gap may never clear by max_n; must return None, not loop forever.
    assert min_n_for_delta(0.51, 0.50, max_n=50) is None


# ── the markdown says when N is too small and what N is needed ─────────────

def _row(task_id: str, n_pass: int, n: int) -> dict:
    flags = [True] * n_pass + [False] * (n - n_pass)
    return {
        "task_id": task_id, "success": n_pass * 2 >= n,
        "success_rate": n_pass / n, "n_samples": n,
        "per_run_success": flags, "quality_0_100": 50,
        "cost_usd": 0.1, "duration_s": 10.0, "tool_calls": 5,
        "input_tokens": 400, "output_tokens": 200,
    }


def test_markdown_says_N_too_small_and_needed_repeats():
    # 3/3 vs 0/3: within noise at N=3; the output must say so and name the
    # repeat count that would settle it (4 per arm for this gap).
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    md = format_compare_markdown(out)
    assert "too small" in md.lower()
    assert "4 repeats per arm" in md


def test_markdown_significant_result_is_not_called_too_small():
    out = compare_runs([_row("t", 9, 10)], [_row("t", 1, 10)])
    md = format_compare_markdown(out)
    assert "beyond noise" in md
    assert "too small" not in md.lower()


def test_markdown_identical_pooled_is_too_small_not_a_number():
    # A zero pooled gap has no resolving N; the message says so instead of
    # printing a bogus figure.
    out = compare_runs([_row("t", 2, 3)], [_row("t", 2, 3)])
    md = format_compare_markdown(out)
    assert "too small" in md.lower()
    assert "zero" in md.lower()
