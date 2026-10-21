"""Adversarial tests for the N-guidance (package D, reviewer s18).

Target: aeebff92 + e2d039db — min_n_for_delta and the "N too small to
decide" guidance in format_compare_markdown.  Round-1 file
(test_bench_stats_adversarial_significance.py) covers the pure helper's
rounding artefact and validation.  This file attacks the MARKDOWN side:
when the guidance fires, what N it names, and whether the message can
mislead a reader about what to run next.
"""

import pytest

from delfin.agent.benchmark import compare_runs, format_compare_markdown


def _row(task_id: str, n_pass: int, n: int) -> dict:
    flags = [True] * n_pass + [False] * (n - n_pass)
    return {
        "task_id": task_id,
        "success": n_pass * 2 >= n,
        "success_rate": n_pass / n,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 50,
        "cost_usd": 0.1,
        "duration_s": 10.0,
        "tool_calls": 5,
        "input_tokens": 400,
        "output_tokens": 200,
    }


# ── A. the guidance must fire exactly when the pooled result is unclear ───

def test_guidance_absent_for_significant_pooled_result():
    out = compare_runs([_row("t", 80, 80)], [_row("t", 40, 80)])
    md = format_compare_markdown(out)
    assert "too small" not in md.lower()


def test_guidance_present_for_thin_verdict_with_overlap():
    # verdict "thin" (< 3 overlap) still computes significance; if the
    # pooled result is not significant the reader must get N guidance,
    # not just a bare verdict.
    out = compare_runs([_row("t", 1, 1)], [_row("t", 0, 1)])
    md = format_compare_markdown(out)
    assert "too small" in md.lower()


# ── B. the named N must actually resolve the observed pooled gap ─────────

def test_named_n_resolves_the_observed_gap():
    # The message names a repeat count; a reader who runs THAT count with
    # the observed rates must land significant.  3/3 vs 0/3 -> promised 4.
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    md = format_compare_markdown(out)
    import re
    m = re.search(r"~(\d+) repeats per arm", md)
    assert m, md
    named = int(m.group(1))
    from delfin.agent.benchmark import _fisher_exact_2x2_pvalue
    p = _fisher_exact_2x2_pvalue(named, 0, 0, named)
    assert p < 0.05, f"named N={named} does not resolve 1.0 vs 0.0 (p={p})"


def test_named_n_is_not_absurdly_larger_than_needed():
    # The helper searches round(p*n) snaps; if it overshoots the honest
    # minimum by more than a couple repeats the message overstates the
    # cost of settling the question.  1.0-vs-0.0 resolves at 4; the
    # named N must be exactly 4, not 5 or 10.
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    md = format_compare_markdown(out)
    import re
    m = re.search(r"~(\d+) repeats per arm", md)
    assert m and int(m.group(1)) == 4


# ── C. the zero-gap message must not fire for unequal rates ──────────────

def test_zero_gap_message_only_for_equal_rates():
    # 2/4 vs 1/2 both have rate 0.5 -> truly equal -> zero message.
    out = compare_runs([_row("t", 2, 4)], [_row("t", 1, 2)])
    md = format_compare_markdown(out)
    assert "zero" in md.lower()

    # 3/4 vs 2/4 rates 0.75 vs 0.5 -> NOT zero; must not say zero.
    out2 = compare_runs([_row("t", 3, 4)], [_row("t", 2, 4)])
    md2 = format_compare_markdown(out2)
    assert "zero" not in md2.lower()
    assert "too small" in md2.lower() or "beyond noise" in md2.lower()


# ── D. mixed directions: pooled guidance must survive a reversed gap ─────

def test_guidance_for_candidate_better_reversed_gap():
    # old 1/10, new 4/10: the candidate is BETTER; the guidance must name
    # the needed N for the same gap and never claim the gap is zero.
    out = compare_runs([_row("t", 1, 10)], [_row("t", 4, 10)])
    md = format_compare_markdown(out)
    assert "zero" not in md.lower()
    import re
    m = re.search(r"~(\d+) repeats per arm", md)
    assert m, "reversed gap lost the N guidance entirely"


# ── E. huge pooled counts: guidance must not spin or lie ─────────────────

def test_guidance_with_massive_pooled_counts_within_noise():
    # 200/200 vs 199/200: nearly equal, not significant; the helper caps
    # at max_n=200.  The message must appear promptly (no long spin) and
    # must not claim zero gap (rates 1.0 vs 0.995 differ).
    out = compare_runs([_row("t", 200, 200)], [_row("t", 199, 200)])
    md = format_compare_markdown(out)
    assert "too small" in md.lower()
    assert "zero" not in md.lower()


def test_guidance_resolves_at_named_n_for_massive_counts():
    # 60/60 vs 30/60: p=1.4e-11 -> SIGNIFICANT, so no guidance at all.
    # This pins that a large-N significant result is not drowned in
    # "too small" noise.
    out = compare_runs([_row("t", 60, 60)], [_row("t", 30, 60)])
    md = format_compare_markdown(out)
    assert "beyond noise" in md
    assert "too small" not in md.lower()
