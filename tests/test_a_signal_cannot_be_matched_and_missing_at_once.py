"""The aggregate said four signals matched AND were missing.

``aggregate_replicates`` unions ``matched_signals`` and
``missing_signals`` alike across replicates. The union is right for
violations — a rule broken once is broken — and it is a contradiction for
the other two: a signal that matched in two samples of three came back in
BOTH lists.

Observed live, 2026-08-13, workflow_plan_before_act at ``--repeats 3``::

    matched_signals = [expected[1], expected[2], expected[3], expected[4]]
    missing_signals = [expected[0], expected[1], expected[2], expected[3],
                       expected[4]]

The audit then printed five missing patterns where exactly one was
genuinely absent, and the four it added were the four the model had
matched in every sample it was scored on. It sent me looking for a
harness defect in the action extractor that does not exist, and it would
send anyone else the same way.

``missing`` keeps the union: a signal the run dropped even once is not
something the run can vouch for. ``matched`` becomes the complement, so
the two are disjoint by construction — and the difference between "never
matched" and "matched twice of three" gets its own name instead of being
flattened into a verdict.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import BenchmarkResult, aggregate_replicates


def _sample(*, matched, missing, quality=50, success=False):
    return BenchmarkResult(
        task_id="t", task_class="c", model="m",
        quality_0_100=quality, success=success,
        matched_signals=list(matched), missing_signals=list(missing),
    )


# ---------------------------------------------------------------------------
# The contradiction
# ---------------------------------------------------------------------------

def test_the_two_lists_never_overlap():
    agg = aggregate_replicates([
        _sample(matched=["a", "b"], missing=["c"]),
        _sample(matched=["a"], missing=["b", "c"]),
        _sample(matched=["a", "b"], missing=["c"]),
    ])
    assert set(agg.matched_signals) & set(agg.missing_signals) == set()


def test_a_signal_that_dropped_once_is_not_reported_as_matched():
    agg = aggregate_replicates([
        _sample(matched=["a", "b"], missing=[]),
        _sample(matched=["a"], missing=["b"]),
    ])
    assert "b" not in agg.matched_signals
    assert "b" in agg.missing_signals


def test_a_signal_that_matched_every_time_is_reported_as_matched():
    agg = aggregate_replicates([
        _sample(matched=["a"], missing=["z"]),
        _sample(matched=["a"], missing=["z"]),
    ])
    assert agg.matched_signals == ["a"]


def test_a_signal_that_never_matched_is_only_missing():
    agg = aggregate_replicates([
        _sample(matched=[], missing=["z"]),
        _sample(matched=[], missing=["z"]),
    ])
    assert agg.missing_signals == ["z"]
    assert agg.flaky_signals == []


# ---------------------------------------------------------------------------
# Flaky is its own fact
# ---------------------------------------------------------------------------

def test_a_signal_that_came_and_went_is_named_flaky():
    agg = aggregate_replicates([
        _sample(matched=["a"], missing=[]),
        _sample(matched=[], missing=["a"]),
    ])
    assert agg.flaky_signals == ["a"]
    assert "a" in agg.missing_signals          # still not vouched for
    assert "a" not in agg.matched_signals


def test_the_live_shape_comes_out_readable():
    """The exact case from the run file: one genuine miss, four that
    matched every scored sample."""
    agg = aggregate_replicates([
        _sample(matched=["e1", "e2", "e3", "e4"], missing=["e0"]),
        _sample(matched=["e1", "e2", "e3", "e4"], missing=["e0"]),
        _sample(matched=["e1", "e2", "e3", "e4"], missing=["e0"]),
    ])
    assert agg.missing_signals == ["e0"]
    assert agg.matched_signals == ["e1", "e2", "e3", "e4"]
    assert agg.flaky_signals == []


def test_a_single_run_has_nothing_to_be_flaky_about():
    agg = aggregate_replicates([_sample(matched=["a"], missing=["b"])])
    assert agg.flaky_signals == []
    assert agg.matched_signals == ["a"]
    assert agg.missing_signals == ["b"]


# ---------------------------------------------------------------------------
# What the union is still right for
# ---------------------------------------------------------------------------

def test_a_violation_seen_once_survives_the_aggregate():
    """The reason the union exists. Unchanged."""
    a = _sample(matched=[], missing=[])
    b = _sample(matched=[], missing=[])
    b.violated_signals = ["forbidden"]
    agg = aggregate_replicates([a, b])
    assert agg.violated_signals == ["forbidden"]


def test_a_budget_violation_seen_once_survives_too():
    a = _sample(matched=[], missing=[])
    b = _sample(matched=[], missing=[])
    b.budget_violations = ["cost"]
    agg = aggregate_replicates([a, b])
    assert agg.budget_violations == ["cost"]


# ---------------------------------------------------------------------------
# The field is optional, so an older run file still loads
# ---------------------------------------------------------------------------

def test_a_result_without_the_field_is_still_a_result():
    r = BenchmarkResult(task_id="t", task_class="c", model="m")
    assert r.flaky_signals == []


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_the_lists_stay_disjoint_for_any_replicate_count(n):
    samples = [_sample(matched=["a", "b"][: (i % 2) + 1], missing=["c"])
               for i in range(n)]
    agg = aggregate_replicates(samples)
    assert set(agg.matched_signals) & set(agg.missing_signals) == set()
