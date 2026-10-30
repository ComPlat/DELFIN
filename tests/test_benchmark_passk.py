"""pass@k (>=1 of k repeats) and pass^k (ALL k repeat) metrics.

These are reported ALONGSIDE the existing per-task pass criterion
(``aggregate_replicates``'s majority ``success``) and never replace a
task's own pass criterion. Both are derived from the recorded
``per_run_success`` flags of one task's k repeats and are therefore
computed offline from recorded results -- no live model call is
involved. They quantify the two ends a single pass rate flattens:

- ``pass_at_k``  (pass@k): at least ONE of the k repeats passed -- the
  optimistic user experience ("it can do it").
- ``pass_pow_k`` (pass^k): ALL k repeats passed -- the reliability a
  user actually experiences day to day.

For a run with one repeat (k = 1) both collapse to the task's
``success``, which is the correct degenerate reading.
"""

import pytest

import delfin.agent.benchmark as bm


def _result(task_id, *, success, quality=70):
    """A single-attempt BenchmarkResult, as score_outcome would emit one."""
    return bm.BenchmarkResult(
        task_id=task_id, task_class="misc", model="model_a",
        mode="solo", ts=1700000000.0,
        success=success, quality_0_100=quality,
        components={"success_pts": 40 if success else 0,
                    "routing_pts": 15, "speed_pts": 10,
                    "cost_pts": 8, "clean_pts": 12},
        duration_s=2.0, cost_usd=0.01, tool_calls=1,
    )


def _aggregate(successes):
    """Aggregate k repeats of one task carrying the given per-run flags."""
    reps = [_result("t", success=s) for s in successes]
    return bm.aggregate_replicates(reps)


# ---------------------------------------------------------------------------
# per-task: pass@k and pass^k over the k repeats
# ---------------------------------------------------------------------------


def test_pass_at_k_true_when_one_of_three_passes():
    """pass@k = at least one repeat passes (the optimistic reading)."""
    agg = _aggregate([True, True, True])
    # Existing criterion unchanged: majority of 3 → success.
    assert agg.success is True
    assert agg.pass_at_k is True
    assert agg.pass_pow_k is True


def test_pass_at_k_true_when_exactly_a_minority_passes():
    """One pass out of three still satisfies pass@k but not pass^k."""
    agg = _aggregate([True, False, False])
    assert agg.success is False          # minority → existing criterion fails
    assert agg.pass_at_k is True         # ... but pass@k holds
    assert agg.pass_pow_k is False       # and pass^k does not
    assert agg.success_rate == pytest.approx(1 / 3)


def test_pass_pow_k_requires_every_repeat():
    """pass^k needs all k: a single flake drops it, exactly the reliability
    a user feels when the same task fails on the second try."""
    agg = _aggregate([True, False, True])
    assert agg.success is True           # majority still passes today
    assert agg.pass_at_k is True
    assert agg.pass_pow_k is False       # ... but it is not reliable


def test_single_repeat_collapses_both_metrics():
    """k = 1: pass@k == pass^k == the task's success."""
    for flag in (True, False):
        agg = _aggregate([flag])
        assert agg.pass_at_k is flag
        assert agg.pass_pow_k is flag
        assert agg.success is flag


# ---------------------------------------------------------------------------
# per model: summary rolls the two new readings of a run
# ---------------------------------------------------------------------------


def test_summarise_run_reports_pass_at_k_and_pass_pow_k():
    """One model's run: counts and rates of tasks reaching each bar."""
    rows = [
        _aggregate([True, True, True]),    # @k and ^k   -- reliable
        _aggregate([True, False, False]),  # @k only     -- flaky minority
        _aggregate([False, False, False]), # neither
    ]
    s = bm.summarise_run(rows)
    assert s["n_tasks"] == 3
    # Existing criterion (majority) unchanged: only row 0 passes today.
    assert s["n_pass"] == 1
    assert s["pass_rate"] == pytest.approx(1 / 3)
    # pass@k: rows 0 and 1 each had >= 1 repeat pass.
    assert s["n_pass_at_k"] == 2
    assert s["pass_at_k_rate"] == pytest.approx(2 / 3)
    # pass^k: only row 0 had every repeat pass.
    assert s["n_pass_pow_k"] == 1
    assert s["pass_pow_k_rate"] == pytest.approx(1 / 3)


def test_summarise_run_metrics_exclude_unmeasured():
    """A task whose request never reached the model stays out of the two
    new rates exactly as it stays out of pass_rate."""
    measured = _aggregate([True, True, True])
    unmeasured = _result("t2", success=False)
    unmeasured.unmeasured = True
    s = bm.summarise_run([measured, unmeasured])
    assert s["n_tasks"] == 1
    assert s["pass_at_k_rate"] == pytest.approx(1.0)
    assert s["pass_pow_k_rate"] == pytest.approx(1.0)


def test_summarise_run_empty_defaults_include_new_keys():
    s = bm.summarise_run([])
    assert s["n_pass_at_k"] == 0
    assert s["pass_at_k_rate"] == 0.0
    assert s["n_pass_pow_k"] == 0
    assert s["pass_pow_k_rate"] == 0.0


# ---------------------------------------------------------------------------
# JSON round-trip: the new fields land in the recorded row
# ---------------------------------------------------------------------------


def test_new_fields_serialize_into_run_row(tmp_path):
    """The JSON result is ``asdict(BenchmarkResult)``; the two metrics must
    survive write_run/read_run so downstream baselines can see them."""
    agg = _aggregate([True, False])
    path = bm.write_run([agg], model="model_a", runs_dir=tmp_path)
    rows = bm.read_run(path)
    assert len(rows) == 1
    row = rows[0]
    assert row["pass_at_k"] is True
    assert row["pass_pow_k"] is False
