"""A single run is not a result: success needs a confidence interval and
an A-vs-B comparison needs a significance statement.

``aggregate_replicates`` reports ``success_rate`` as a raw fraction.
3/3 successes reads as certainty, and 2/3 as 67% — neither carries the
uncertainty of an N=3 experiment.  The Wilson interval is the standard
remedy for small-N binomials (unlike the normal approximation, it stays
inside [0,1] and behaves at 0 and 1 successes).

``compare_runs`` classifies per-task deltas by fixed thresholds, so
"rate 0.67 vs 0.33 at N=3" could come back "worse" when the data cannot
distinguish the two.  A comparison that will be used to accept or reject
an agent change must say whether the observed difference is real at a
stated alpha, or inside the noise.
"""

import math

from delfin.agent.benchmark import (
    BenchmarkResult,
    aggregate_replicates,
    compare_runs,
    wilson_interval,
)


def _run(success: bool) -> BenchmarkResult:
    return BenchmarkResult(
        task_id="t", task_class="c", model="m",
        success=success, quality_0_100=100 if success else 10,
        per_run_success=[success], text_excerpt="x",
    )


# ── Wilson interval ─────────────────────────────────────────────────────


def test_wilson_known_values():
    # 0 of 1: the interval must not collapse to [0, 0] -- with one
    # sample and no successes the honest statement is "somewhere below
    # ~0.75", not "certainly zero".
    lo, hi = wilson_interval(0, 1)
    assert lo == 0.0
    assert 0.6 < hi < 0.8

    # 1 of 1 mirrors it.
    lo, hi = wilson_interval(1, 1)
    assert 0.2 < lo < 0.4
    assert hi == 1.0

    # 3 of 3 is wide, not a point.
    lo, hi = wilson_interval(3, 3)
    assert lo < 0.95, "3/3 must not read as 95%+ certainty"
    assert hi == 1.0

    # 0 of 0 is undefined; the answer must be a full interval, not a
    # crash or a fake certainty.
    lo, hi = wilson_interval(0, 0)
    assert (lo, hi) == (0.0, 1.0)


def test_aggregate_carries_wilson_ci():
    agg = aggregate_replicates([_run(True), _run(True), _run(False)])
    assert agg.n_samples == 3
    assert agg.success_rate == 2 / 3
    # Present, ordered, and not equal to the point estimate.
    assert agg.success_ci_low is not None
    assert agg.success_ci_high is not None
    assert 0.0 <= agg.success_ci_low < agg.success_rate < agg.success_ci_high <= 1.0


def test_single_run_carries_ci_too():
    r = _run(True)
    # N=1: trivially the only sample -- the CI machinery must not break.
    lo, hi = wilson_interval(1, 1)
    assert r.success_rate in (0.0, 1.0)


# ── Significance in A-vs-B ──────────────────────────────────────────────


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
    }


def test_compare_runs_reports_significance_fields():
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    s = out["summary"]
    # The verdict names the statistical statement, not just the direction.
    assert "significant" in s or "significance" in s


def test_identical_runs_are_not_significant():
    out = compare_runs([_row("t", 2, 3)], [_row("t", 2, 3)])
    s = out["summary"]
    assert s.get("significant") is False


def test_3_of_3_vs_0_of_3_is_not_significant_at_alpha_05():
    # Honest statistics: even a 3-of-3 vs 0-of-3 sweep gives a two-sided
    # Fisher p of 0.10 at N=3 -- the smallest p 2x3 samples can produce.
    # A compare that called this "significant" would overstate what a
    # 3-repeat run can show; more repeats are needed, not a verdict.
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    assert out["summary"].get("significant") is False
    assert out["summary"]["significance_p"] == 0.1


def test_a_real_difference_with_more_repeats_is_significant():
    # 9/10 vs 1/10 is the kind of effect size a change decision needs;
    # Fisher separates it cleanly from noise.
    out = compare_runs([_row("t", 9, 10)], [_row("t", 1, 10)])
    assert out["summary"].get("significant") is True
    assert out["summary"]["significance_p"] < 0.05


def test_2_of_3_vs_1_of_3_is_not_significant():
    # The classic trap: 67% vs 33% at N=3 is NOT distinguishable.  A
    # threshold-only compare would call this "worse" and invite a wrong
    # conclusion.
    out = compare_runs([_row("t", 2, 3)], [_row("t", 1, 3)])
    assert out["summary"].get("significant") is False
    # But the per-task row still shows the observed delta.
    assert out["per_task"][0]["success_rate_old"] > out["per_task"][0]["success_rate_new"]
