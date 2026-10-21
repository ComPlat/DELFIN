"""Adversarial tests for the significance/noise labels (package D, reviewer s18).

Target: the "bench compare" statistical statement the builder added in
commits 9de593bf..e2d039db.  The goal is to BREAK the fix, not to repeat
its happy path: contradictions between the printed p and the printed
label, zero-N rows, duplicate task ids, missing noise fields, and the
N-guidance promise.
"""

import pytest

from delfin.agent.benchmark import (
    _fisher_exact_2x2_pvalue,
    classify_pass_delta,
    compare_runs,
    format_compare_markdown,
    min_n_for_delta,
    _SIG_ALPHA,
)


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


# ── A. the printed p must never contradict the printed label ──────────────

def test_markdown_p_not_inverted_for_a_very_small_p():
    # compare_runs stores significance_p ROUNDED to 4 decimals; a pooled
    # sweep of 80/80 -> 40/80 has Fisher p = 2.5e-15, which rounds to 0.0.
    # The markdown reader `float(summary.get("significance_p") or 1.0)`
    # then treats the falsy 0.0 as "absent" and prints p=1 -- the exact
    # opposite of the truth, beside a "beyond noise" label.
    out = compare_runs([_row("t", 80, 80)], [_row("t", 40, 80)])
    s = out["summary"]
    assert s["significant"] is True
    md = format_compare_markdown(out)
    # The p that is printed must be compatible with the label: if the
    # result is beyond noise at alpha .05 the printed p must be < 0.05.
    sig_line = [ln for ln in md.splitlines() if "Fisher exact" in ln][0]
    assert "p=1 " not in sig_line and "p=1." not in sig_line
    assert "beyond noise" in md


def test_significance_p_zero_rounds_without_becoming_one():
    # The rounded zero must reach the markdown AS a zero (p < 0.0001 is
    # how strong significance reads), not as 1.0 through an `or` fallback.
    out = compare_runs([_row("t", 80, 80)], [_row("t", 40, 80)])
    assert out["summary"]["significance_p"] == 0.0


def test_markdown_p_label_agreement_mid_range():
    # Sanity companion at a realistic p: label beyond, p printed < alpha.
    out = compare_runs([_row("t", 9, 10)], [_row("t", 1, 10)])
    md = format_compare_markdown(out)
    line = [ln for ln in md.splitlines() if "Fisher exact" in ln][0]
    assert "p=0.0011" in line
    assert "beyond noise" in md


# ── B. degenerate inputs must not lie ─────────────────────────────────────

def test_per_task_counts_zero_n_samples_is_not_one():
    # A row with n_samples=0 and no per_run_success must contribute ZERO
    # samples, not one phantom sample: `int(row.get("n_samples") or 1)`
    # reads 0 as absent and invents a trial that never ran, which then
    # lands in the pooled counts and in every delta label.
    from delfin.agent.benchmark import _per_task_counts
    row = {"task_id": "t", "success": False, "n_samples": 0,
           "per_run_success": []}
    assert _per_task_counts(row) == (0, 0)


def test_pooled_counts_exclude_zero_sample_rows():
    # The pooled Fisher table must not gain phantom trials from rows that
    # never ran (n_samples=0).  Two zero-sample arms give p=1 and no
    # significance, and the markdown's zero-gap branch fires -- but the
    # pooled counts must be [0, 0], not [0, 1].
    out = compare_runs(
        [{"task_id": "t", "success": False, "n_samples": 0,
          "per_run_success": []}],
        [{"task_id": "t", "success": False, "n_samples": 0,
          "per_run_success": []}],
    )
    assert out["summary"]["pooled_success"] == {"old": [0, 0], "new": [0, 0]}


def test_markdown_skips_significance_when_no_trials_either_side():
    # With zero trials pooled, the section must not print a p for a test
    # that was never performed.
    out = compare_runs(
        [{"task_id": "t", "success": False, "n_samples": 0,
          "per_run_success": []}],
        [{"task_id": "t", "success": False, "n_samples": 0,
          "per_run_success": []}],
    )
    md = format_compare_markdown(out)
    assert "Fisher exact" not in md


def test_fisher_with_zero_total_does_not_crash():
    assert _fisher_exact_2x2_pvalue(0, 0, 0, 0) == 1.0


# ── C. duplicate task ids silently shrink the pooled counts ───────────────

def test_compare_runs_duplicate_ids_collide():
    # A run file with the same task twice (repeats written un-aggregated)
    # collapses through `by_id_old = {r.get("task_id"): r for r in ...}`:
    # the SECOND row silently replaces the first and the pooled counts
    # halve with no notice.  Pre-existing behaviour, but the new pooled
    # significance inherits it, so it is now a lie about certainty.
    out = compare_runs([_row("t", 3, 3), _row("t", 3, 3)],
                       [_row("t", 0, 3)])
    assert out["summary"]["pooled_success"]["old"][1] == 6


def test_compare_runs_duplicate_ids_reported():
    # If duplicates exist the compare must SAY so rather than silently
    # collapsing: a `duplicate_task_ids` list in the summary, or the
    # pooled counts visibly smaller than the sum of per-task counts.
    out = compare_runs([_row("t", 3, 3), _row("t", 3, 3)],
                       [_row("t", 0, 3)])
    n = out["summary"]["pooled_success"]["old"]
    per_task_n = sum(r["success_counts_old"][1] for r in out["per_task"])
    assert n[1] == per_task_n  # pooled must equal the sum it claims to pool


# ── D. missing delta_noise field must not read as "within" ────────────────

def test_markdown_noise_column_marks_missing_label():
    # `row.get("delta_noise") or "within"` prints a MISSING label as
    # "within noise" -- a default that asserts certainty the data never
    # gave.  A row built by an older caller (or a hand-built cmp_result)
    # without delta_noise must render as unknown, not as within.
    cmp_result = {
        "summary": {"n_overlap": 1, "n_better": 0, "n_worse": 0,
                    "n_neutral": 1,
                    "old": {"pass_rate": 1.0, "avg_quality": 50.0,
                            "total_cost_usd": 0.1, "total_duration_s": 10.0},
                    "new": {"pass_rate": 0.0, "avg_quality": 50.0,
                            "total_cost_usd": 0.1, "total_duration_s": 10.0},
                    "significant": False, "significance_p": 1.0,
                    "significance_alpha": _SIG_ALPHA,
                    "pooled_success": {"old": [3, 3], "new": [0, 3]}},
        "per_task": [{"task_id": "t", "class": "worse",
                      "old_quality": 50, "new_quality": 50,
                      "d_quality": 0, "d_cost_usd": 0.0,
                      "d_duration_s": 0.0}],
        "verdict": "worse",
    }
    md = format_compare_markdown(cmp_result)
    row_line = [ln for ln in md.splitlines() if ln.startswith("| `t")][0]
    # Either the field is optional and printed as unknown, or it is
    # required and the reader recomputes it -- "within" as a silent
    # default is the one answer that must not appear.
    assert "within" not in row_line


# ── E. min_n_for_delta: the promise must survive its own snapping ─────────

def test_min_n_promise_holds_at_returned_n():
    # The docstring promise: at the RETURNED n, running the observed rates
    # snaps to integer counts whose Fisher p is below alpha.  With
    # banker's rounding at p*n = x.5 the snap can be one SHORT of the
    # extreme count, and the p at the returned n can sit above alpha
    # (observed: 0.5-vs-0.0 returns 10, but n=9 with a 5/9 snap is
    # already p=0.029; the returned 10 is a rounding artefact, not the
    # honest minimum).  The returned value must be the true minimum.
    assert min_n_for_delta(0.5, 0.0) == 9


def test_min_n_uses_floor_never_ceiling_for_low_arm():
    # Symmetric check: the LOW arm must snap to its floor so the gap is
    # never understated.
    from delfin.agent.benchmark import _fisher_exact_2x2_pvalue
    n = min_n_for_delta(0.51, 0.0)
    a = round(0.51 * n)
    assert _fisher_exact_2x2_pvalue(a, n - a, 0, n) < _SIG_ALPHA


def test_min_n_alpha_zero_or_negative_returns_none_not_hang():
    # alpha <= 0 can never be beaten (p >= 0); must return None, not
    # spin 200 Fisher evaluations pretending to search.
    assert min_n_for_delta(1.0, 0.0, alpha=0.0) is None


def test_min_n_rates_outside_unit_interval_rejected():
    # Rates above 1 or below 0 are corrupted data, not a bigger effect.
    with pytest.raises(ValueError):
        min_n_for_delta(1.5, 0.5)
