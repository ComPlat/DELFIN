"""A suite can hold 11/11 while every answer gets harder to read.

Pass rate and quality are hit-side numbers. Nearly every fix in this
framework has one shape -- a path that used to be silent now says
something -- and each one spends a caveat, a check or a refusal. None of
that can move a pass/fail number. So a wave of guard work can leave the
score untouched and the product worse, and nothing would report it.

These are the numbers that can move: caveats per answer, output length,
and how often the gates refused. The rule they encode is that a rise is a
regression even when the rate holds, and that a quantity nobody could
observe is reported as unobserved rather than as zero.
"""

from __future__ import annotations

import types

from delfin.agent.benchmark import (
    BenchmarkResult, Trajectory, caveat_count, score_outcome, summarise_run,
)
from delfin.agent import benchmark_runner as br


# ---------------------------------------------------------------------------
# Counting what the framework actually emits
# ---------------------------------------------------------------------------

def test_the_markers_counted_are_the_ones_the_framework_writes():
    """Not keyword-guessing: these two strings are what the caveat
    builders prepend, so this counts what was really appended."""
    from delfin.agent import verify_guard as vg
    from delfin.agent import office

    emitted = [
        vg.count_vs_enumeration_caveat([(5, 3)]),
        vg.truncated_output_caveat(["12"], ["read_document"]),
        vg.ambiguous_column_caveat(["Betrag"]),
        office.figure_caveat([]),
    ]
    for text in emitted:
        if text:
            assert caveat_count(text) >= 1, repr(text[:60])


def test_three_caveats_count_as_three():
    answer = ("Die Summe ist 12,50 EUR."
              "\n\n> ⚠️ eins"
              "\n\n> ⚠️ zwei"
              "\n\n[verify] Caveat: drei")
    assert caveat_count(answer) == 3


def test_a_clean_answer_costs_nothing():
    assert caveat_count("Die Summe ist 12,50 EUR.") == 0
    assert caveat_count("") == 0
    assert caveat_count(None) == 0


# ---------------------------------------------------------------------------
# ... and reaches the score card
# ---------------------------------------------------------------------------

def _task(**kw):
    from delfin.agent.benchmark import Task
    return Task(id="t", prompt="p", task_class="office", mode="office", **kw)


def test_the_score_card_carries_the_cost_side():
    traj = Trajectory(text="Antwort.\n\n> ⚠️ eins\n\n> ⚠️ zwei",
                      output_tokens=140)
    out = score_outcome(_task(), traj, model="m")
    assert out.caveats == 2
    assert out.answer_chars == len(traj.text)


def test_the_summary_reports_hedging_beside_the_rate():
    rows = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, quality_0_100=93, caveats=0,
                        output_tokens=100, answer_chars=200),
        BenchmarkResult(task_id="b", task_class="office", model="m",
                        success=True, quality_0_100=93, caveats=4,
                        output_tokens=300, answer_chars=900),
    ]
    s = summarise_run(rows)
    # The thing this whole file exists for: a perfect rate, and the cost
    # side still has something to say.
    assert s["pass_rate"] == 1.0
    assert s["avg_caveats"] == 2.0
    assert s["max_caveats"] == 4
    assert s["avg_output_tokens"] == 200.0
    assert s["n_fail"] == 0


# ---------------------------------------------------------------------------
# Unobserved is not zero
# ---------------------------------------------------------------------------

def test_denials_are_unobserved_when_nobody_could_look():
    """A mocked engine has no workspace and no audit log. Reporting 0
    there says "nothing was refused", which is a different sentence from
    "nobody looked" -- and this framework has shipped that confusion
    before."""
    assert br._denials_during(types.SimpleNamespace(), since_ts="") is None
    assert br._denials_during(object(), since_ts="") is None


def test_one_unobserved_task_makes_the_suite_total_unobserved():
    rows = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, denials=2),
        BenchmarkResult(task_id="b", task_class="office", model="m",
                        success=True, denials=None),
    ]
    assert summarise_run(rows)["total_denials"] is None


def test_a_suite_that_did_look_reports_the_number():
    rows = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, denials=2),
        BenchmarkResult(task_id="b", task_class="office", model="m",
                        success=True, denials=0),
    ]
    assert summarise_run(rows)["total_denials"] == 2


def test_an_unmeasured_task_does_not_contribute_to_the_cost_side_either():
    """The same rule the rate already follows: a request that never
    reached the model says nothing about hedging."""
    rows = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, quality_0_100=90, caveats=0),
        BenchmarkResult(task_id="b", task_class="office", model="m",
                        success=False, quality_0_100=0, caveats=9,
                        unmeasured=True),
    ]
    s = summarise_run(rows)
    assert s["avg_caveats"] == 0.0 and s["max_caveats"] == 0
    assert s["n_unmeasured"] == 1


# ---------------------------------------------------------------------------
# Aggregation across replicates
# ---------------------------------------------------------------------------

def test_replicates_report_the_median_hedging(monkeypatch):
    from delfin.agent.benchmark import aggregate_replicates
    reps = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, quality_0_100=90, caveats=c,
                        answer_chars=100 * (i + 1))
        for i, c in enumerate((1, 3, 3))
    ]
    agg = aggregate_replicates(reps)
    assert agg.caveats == 3
    assert agg.answer_chars == 200


def test_three_unmeasured_replicates_are_not_a_failing_model():
    """The single-run path had this guard and the replicate path dropped
    the flag, so --repeats bypassed the whole mechanism. Found live: an
    engine that never started -- a missing provider argument -- came back
    0/11 at quality 35 with no NOT MEASURED notice, on its way into the
    file baselines are compared against."""
    from delfin.agent.benchmark import aggregate_replicates
    reps = [BenchmarkResult(task_id="a", task_class="office", model="m",
                            success=False, quality_0_100=35, unmeasured=True,
                            error="engine init failed: no api key")
            for _ in range(3)]
    agg = aggregate_replicates(reps)
    assert agg.unmeasured is True
    s = summarise_run([agg])
    assert s["n_unmeasured"] == 1
    assert s["n_tasks"] == 0, "an outage must not be counted as a scored task"


def test_one_outage_among_three_does_not_drag_the_median_down():
    """The commoner shape: two replicates ran, the third hit an endpoint
    with no capacity. Medianing its zero in punishes the model for the
    network."""
    from delfin.agent.benchmark import aggregate_replicates
    reps = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, quality_0_100=93),
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, quality_0_100=91),
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=False, quality_0_100=35, unmeasured=True,
                        error="no deployments available"),
    ]
    agg = aggregate_replicates(reps)
    assert agg.unmeasured is False
    assert agg.quality_0_100 == 92
    assert agg.n_samples == 2, "the outage is not a sample of the model"
    assert agg.success is True


def test_one_unobserved_replicate_leaves_the_aggregate_unobserved():
    from delfin.agent.benchmark import aggregate_replicates
    reps = [
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, denials=1),
        BenchmarkResult(task_id="a", task_class="office", model="m",
                        success=True, denials=None),
    ]
    assert aggregate_replicates(reps).denials is None
