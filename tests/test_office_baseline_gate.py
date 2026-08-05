"""The office benchmark's reference standard and the gate that reads it.

The hard part is not comparing two numbers, it is not crying wolf. The
measured reference has one task scoring 45/93/93 across three runs of
unchanged code -- a gate that fires on any drop fires on that, and a gate
that is always red is a gate nobody reads.

These tests pin both directions: the gate must stay quiet for the noise
the reference itself recorded, and it must still catch a task that dies
and a suite that sinks.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from delfin.agent.benchmark_baseline import (
    SCHEMA_VERSION,
    Baseline,
    TaskReference,
    baseline_from_results,
    compare_to_baseline,
    format_baseline_report,
    load_baseline,
    save_baseline,
)

BASELINE_PATH = pathlib.Path(__file__).parent / "fixtures" / "office_baseline.json"


def _row(task_id, n_samples, per_run_success, quality=90):
    return {
        "task_id": task_id,
        "n_samples": n_samples,
        "per_run_success": list(per_run_success),
        "success_rate": sum(1 for x in per_run_success if x) / n_samples,
        "quality_0_100": quality,
        "quality_stdev": 0.5,
        "model": "test-model",
        "mode": "office",
    }


def _reference(**tasks):
    """Build a reference from {task_id: (n_passed, n_samples)}."""
    return Baseline(
        model="test-model", mode="office", measured_at="2026-07-28",
        commit="deadbeef", repeats=3,
        tasks={
            tid: TaskReference(task_id=tid, n_samples=n, n_passed=p,
                               quality_mean=90.0, quality_stdev=1.0)
            for tid, (p, n) in tasks.items()
        },
    )


_SIX_SOLID = dict(t1=(3, 3), t2=(3, 3), t3=(3, 3), t4=(3, 3), t5=(3, 3), t6=(3, 3))


# ---------------------------------------------------------------------------
# The committed reference
# ---------------------------------------------------------------------------

def test_the_committed_reference_loads():
    b = load_baseline(BASELINE_PATH)
    assert b is not None, f"no reference at {BASELINE_PATH}"
    assert b.mode == "office"
    assert b.tasks, "a reference with no tasks gates nothing"
    assert b.commit, "the reference must record the code it was measured on"
    assert b.measured_at, "an undated reference cannot be judged for staleness"


def test_the_committed_reference_matches_the_run_it_came_from():
    """Guards against hand-edited numbers: every task is whole samples."""
    b = load_baseline(BASELINE_PATH)
    for tid, t in b.tasks.items():
        assert 0 <= t.n_passed <= t.n_samples, f"{tid}: {t.n_passed}/{t.n_samples}"
        if t.per_run_success:
            assert len(t.per_run_success) == t.n_samples, tid
            assert sum(1 for x in t.per_run_success if x) == t.n_passed, (
                f"{tid}: the recorded samples do not add up to n_passed")


def test_a_reference_survives_its_own_round_trip_exactly(tmp_path):
    """Storing a rounded rate instead of counts made this fail before."""
    b = load_baseline(BASELINE_PATH)
    p = tmp_path / "again.json"
    save_baseline(b, p)
    assert load_baseline(p).to_dict() == b.to_dict()
    assert load_baseline(p).suite_pass_rate == b.suite_pass_rate


def test_an_unknown_schema_is_refused(tmp_path):
    """Numbers whose meaning may have changed must not be reinterpreted."""
    p = tmp_path / "future.json"
    p.write_text(json.dumps({"schema": SCHEMA_VERSION + 1, "tasks": {}}))
    with pytest.raises(ValueError, match="schema"):
        load_baseline(p)


def test_a_missing_reference_is_not_an_error(tmp_path):
    assert load_baseline(tmp_path / "nothing.json") is None


# ---------------------------------------------------------------------------
# The gate stays quiet for noise
# ---------------------------------------------------------------------------

def test_an_identical_run_passes():
    ref = _reference(**_SIX_SOLID)
    rows = [_row(t, 3, [True] * 3) for t in _SIX_SOLID]
    assert compare_to_baseline(rows, ref)["verdict"] == "pass"


def test_one_flaky_sample_does_not_fail_the_gate():
    """A single sample flipping is what the reference itself recorded."""
    ref = _reference(**_SIX_SOLID)
    rows = [_row("t1", 3, [True, False, True])] + \
           [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[1:]]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "pass", result["regressions"]
    assert result["warnings"], "it should still be reported, just not fatal"


def test_the_known_unstable_task_may_fail_entirely_without_failing_the_gate():
    """The reference records this task at 2/3, and a task that unstable
    returning 0/3 happens about once in 27 runs. The gate must sit through
    it -- both the per-task rule, which does not apply to a task the
    reference already knows is shaky, and the suite rule, whose tolerance
    is fed by that very instability."""
    ref = _reference(t1=(2, 3), t2=(3, 3), t3=(3, 3), t4=(3, 3), t5=(3, 3), t6=(3, 3))
    rows = [_row("t1", 3, [False] * 3)] + \
           [_row(t, 3, [True] * 3) for t in ("t2", "t3", "t4", "t5", "t6")]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "pass", result["regressions"]
    assert "t1" not in " ".join(result["regressions"])
    assert result["warnings"], "it must still be reported"


def test_an_unstable_task_widens_the_tolerance_only_for_itself():
    """Pooling the suite's noise would smear one shaky task's variance over
    five solid ones, and hide a real break in the solid part."""
    from delfin.agent.benchmark_baseline import _suite_tolerance

    with_flaky = _reference(t1=(2, 3), t2=(3, 3), t3=(3, 3),
                            t4=(3, 3), t5=(3, 3), t6=(3, 3))
    all_solid = _reference(**_SIX_SOLID)
    assert _suite_tolerance(all_solid) < _suite_tolerance(with_flaky)

    # ...and the solid reference still catches a task that breaks.
    rows = [_row("t1", 3, [False] * 3)] + \
           [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[1:]]
    assert compare_to_baseline(rows, all_solid)["verdict"] == "regressed"


def test_a_run_that_beats_the_reference_passes():
    ref = _reference(t1=(2, 3), t2=(3, 3), t3=(3, 3), t4=(3, 3), t5=(3, 3), t6=(3, 3))
    rows = [_row(t, 3, [True] * 3) for t in _SIX_SOLID]
    assert compare_to_baseline(rows, ref)["verdict"] == "pass"


# ---------------------------------------------------------------------------
# The gate still catches real regressions
# ---------------------------------------------------------------------------

def test_a_solid_task_dying_outright_is_a_regression():
    ref = _reference(**_SIX_SOLID)
    rows = [_row("t1", 3, [False] * 3)] + \
           [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[1:]]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "regressed"
    assert any("t1" in r for r in result["regressions"])


def test_a_sinking_suite_is_a_regression_even_with_no_task_dead():
    """Every task half-failing: no single task is dead, the suite is."""
    ref = _reference(**_SIX_SOLID)
    rows = [_row(t, 3, [True, False, False]) for t in _SIX_SOLID]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "regressed"
    assert any("suite pass rate" in r for r in result["regressions"])


def test_a_task_that_stops_being_measured_is_a_regression():
    """Deleting a failing task is not a way to make the gate green."""
    ref = _reference(**_SIX_SOLID)
    rows = [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[:-1]]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "regressed"
    assert any("t6" in r for r in result["regressions"])


def test_too_little_overlap_is_reported_as_thin_not_as_pass():
    """A run sharing almost nothing with the reference proves nothing."""
    ref = _reference(**_SIX_SOLID)
    rows = [_row("t1", 3, [True] * 3)]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "thin"


def test_a_new_task_is_reported_but_not_gated():
    ref = _reference(**_SIX_SOLID)
    rows = [_row(t, 3, [True] * 3) for t in _SIX_SOLID] + \
           [_row("t_new", 3, [False] * 3)]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "pass"
    assert any("t_new" in n for n in result["notes"])


# ---------------------------------------------------------------------------
# Tolerance behaviour
# ---------------------------------------------------------------------------

def test_a_perfect_reference_still_tolerates_one_sample():
    """With a reference at 100% the standard error is zero; without a floor
    any single failure anywhere would be a regression."""
    ref = _reference(**_SIX_SOLID)
    rows = [_row("t1", 3, [True, True, False])] + \
           [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[1:]]
    assert compare_to_baseline(rows, ref)["verdict"] == "pass"


def test_more_reference_samples_give_a_tighter_gate():
    """The tolerance is derived from the reference, not fixed, so measuring
    more makes the ratchet stricter without editing a threshold."""
    from delfin.agent.benchmark_baseline import _suite_tolerance

    small = _reference(**{f"t{i}": (3, 3) for i in range(1, 7)})
    large = _reference(**{f"t{i}": (30, 30) for i in range(1, 7)})
    assert _suite_tolerance(large) < _suite_tolerance(small)


# ---------------------------------------------------------------------------
# Reading a run row
# ---------------------------------------------------------------------------

def _single_sample_row(task_id, success):
    """Exactly what a repeats=1 run writes.

    ``success_rate`` is present and 0.0 -- the dataclass default, only ever
    filled when repeats > 1 -- and ``per_run_success`` is empty. Synthetic
    rows that set success_rate honestly cannot expose the trap.
    """
    return {
        "task_id": task_id, "n_samples": 1, "per_run_success": [],
        "success": success, "success_rate": 0.0, "quality_0_100": 91,
        "quality_stdev": 0.0, "model": "test-model", "mode": "office",
    }


def test_a_passing_single_sample_run_is_not_read_as_a_total_failure():
    """The field exists but means nothing at N=1; trusting it scored a
    perfect run as zero and called every task dead."""
    from delfin.agent.benchmark_baseline import observed

    assert observed(_single_sample_row("t1", True)) == (1, 1)
    assert observed(_single_sample_row("t1", False)) == (1, 0)


def test_the_gate_agrees_with_a_single_sample_run():
    ref = _reference(**_SIX_SOLID)
    rows = [_single_sample_row(t, True) for t in _SIX_SOLID]
    result = compare_to_baseline(rows, ref)
    assert result["verdict"] == "pass", result["regressions"]
    assert result["summary"]["candidate_pass_rate"] == 1.0


def test_recorded_samples_win_over_a_summary_that_disagrees():
    from delfin.agent.benchmark_baseline import observed

    row = _row("t1", 3, [True, False, True])
    row["success_rate"] = 0.9999          # a summary the samples contradict
    assert observed(row) == (3, 2)


def test_a_multi_sample_row_without_recorded_samples_uses_the_rate():
    from delfin.agent.benchmark_baseline import observed

    assert observed({"task_id": "t1", "n_samples": 4,
                     "per_run_success": [], "success_rate": 0.75}) == (4, 3)


# ---------------------------------------------------------------------------
# Building a reference
# ---------------------------------------------------------------------------

def test_a_reference_is_built_from_the_recorded_samples_not_the_rate():
    rows = [_row("t1", 3, [True, False, True])]
    rows[0]["success_rate"] = 0.9999      # a lie the samples contradict
    b = baseline_from_results(rows, measured_at="2026-07-28")
    assert b.tasks["t1"].n_passed == 2


def test_an_empty_run_cannot_become_a_reference():
    with pytest.raises(ValueError):
        baseline_from_results([], measured_at="2026-07-28")


def test_the_report_names_what_went_wrong():
    ref = _reference(**_SIX_SOLID)
    rows = [_row("t1", 3, [False] * 3)] + \
           [_row(t, 3, [True] * 3) for t in list(_SIX_SOLID)[1:]]
    text = format_baseline_report(compare_to_baseline(rows, ref), ref)
    assert "REGRESSION" in text and "t1" in text
    assert "re-measure" in text, "a red gate must say what to do about it"
