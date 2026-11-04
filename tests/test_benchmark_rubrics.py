"""Optional weighted rubrics per task, scored deterministically.

A ``rubric`` is an OPTIONAL list of weighted criteria on a task (planning,
geometry, input, execution, analysis ...), each scored by a deterministic
check against the recorded trajectory. The rubric score is reported
ALONGSIDE ``success`` and ``quality_0_100`` and never rewrites either:
the existing per-task pass criterion stays exactly as it is.

Numerical criteria are judged with a STATED TOLERANCE through
``ExpectedValue`` (relative tolerance + absolute floor) -- never by string
matching the digits an answer happened to write.

All scores here are computed offline from recorded results; no live model
call is involved.
"""

import pytest

import delfin.agent.benchmark as bm
from delfin.agent.benchmark import (
    ExpectedValue,
    RubricCriterion,
    Task,
    Trajectory,
)


def _task(rubric=None, *, prompt="do it"):
    return Task(
        id="t_rubric", task_class="misc", mode="solo", prompt=prompt,
        rubric=rubric,
    )


def _traj(*, text="done", error="", verify_ok=None):
    return Trajectory(text=text, error=error, verify_ok=verify_ok)


# ---------------------------------------------------------------------------
# weight combination
# ---------------------------------------------------------------------------


def test_rubric_combines_weights_into_0_100():
    """score = round(100 * sum(w*criterion) / sum(w))."""
    rubric = [
        RubricCriterion(name="planning", weight=2, kind="signals"),
        RubricCriterion(name="execution", weight=1, kind="verify"),
    ]
    # planning hits all required signals (none -> 1.0), verify fails.
    task = _task(rubric)
    res = bm.score_outcome(task, _traj(text="ok", verify_ok=False))
    assert res.rubric_components["planning"] == pytest.approx(
        100 * 2 / 3)
    assert res.rubric_components["execution"] == 0.0
    assert res.rubric_0_100 == int(round(100 * (2 * 1.0 + 1 * 0.0) / 3))


def test_rubric_half_credit_midway():
    rubric = [
        RubricCriterion(name="planning", weight=1, kind="signals"),
        RubricCriterion(name="execution", weight=1, kind="verify"),
    ]
    task = _task(rubric)
    res = bm.score_outcome(task, _traj(text="ok", verify_ok=True))
    assert res.rubric_0_100 == 100


# ---------------------------------------------------------------------------
# numerical values are judged with a stated tolerance, never string-match
# ---------------------------------------------------------------------------


def test_rubric_value_criterion_uses_tolerance_not_string_match():
    """A value criterion passes within ExpectedValue's tolerance and fails
    outside it -- an answer's DIALECT of the number must not matter."""
    expected = ExpectedValue(value=1.50, tolerance=0.01)
    rubric = [
        RubricCriterion(name="analysis", weight=1, kind="values"),
    ]
    task = Task(
        id="t", task_class="misc", mode="solo", prompt="p",
        expected_values=(expected,),
        rubric=rubric,
    )
    # 1,50 (German decimal) within +/-1% of 1.50 -> matched by tolerance.
    res_ok = bm.score_outcome(task, _traj(text="1,50"))
    assert expected.judge("1,50") == "matched"
    assert res_ok.rubric_0_100 == 100
    # 1.60 is 6.7% away -> outside tolerance, judged wrong (not absent).
    res_wrong = bm.score_outcome(task, _traj(text="1,60"))
    assert expected.judge("1,60") == "wrong"
    assert res_wrong.rubric_0_100 == 0


# ---------------------------------------------------------------------------
# unknown rubric kinds are rejected loudly, not scored 0 silently
# ---------------------------------------------------------------------------


def test_unknown_rubric_kind_raises():
    rubric = [RubricCriterion(name="x", weight=1, kind="no_such_check")]
    with pytest.raises(ValueError):
        bm.score_outcome(_task(rubric), _traj())


# ---------------------------------------------------------------------------
# the rubric never rewrites the existing pass criterion
# ---------------------------------------------------------------------------


def test_rubric_does_not_alter_success():
    """A task whose existing required signal fails still fails; the rubric is
    an alongside metric, never a way back to a pass."""
    bad = Task(
        id="t", task_class="misc", mode="solo", prompt="p",
        expected_signals=(bm.Signal(pattern="NEVER_SAID"),),
        rubric=[RubricCriterion(name="planning", weight=1, kind="signals")],
    )
    res = bm.score_outcome(bad, _traj(text=""))
    assert res.success is False          # existing criterion unaffected
    assert res.rubric_0_100 == 0         # ... and the rubric still scored


# ---------------------------------------------------------------------------
# aggregation and JSON persistence
# ---------------------------------------------------------------------------


def test_rubric_aggregates_as_median_across_replicates():
    rubric = [RubricCriterion(name="execution", weight=1, kind="verify")]
    task = _task(rubric)
    r1 = bm.score_outcome(task, _traj(text="a", verify_ok=True))
    r2 = bm.score_outcome(task, _traj(text="b", verify_ok=False))
    agg = bm.aggregate_replicates([r1, r2])
    # Median of [100, 0] is 50 for the rubric score and its components.
    assert agg.rubric_0_100 == 50
    assert agg.rubric_components["execution"] == pytest.approx(50.0)


def test_rubric_fields_serialize_into_run_row(tmp_path):
    rubric = [RubricCriterion(name="analysis", weight=1, kind="verify")]
    res = bm.score_outcome(_task(rubric), _traj(text="a", verify_ok=True))
    path = bm.write_run([res], model="m", runs_dir=tmp_path)
    row = bm.read_run(path)[0]
    assert row["rubric_0_100"] == 100
    assert row["rubric_components"]["analysis"] == pytest.approx(100.0)


# Task loaders accept a rubric from YAML-like dicts.
def test_coerce_task_reads_rubric():
    raw = {
        "id": "t", "task_class": "misc", "mode": "solo", "prompt": "p",
        "rubric": [
            {"name": "geometry", "weight": 3, "kind": "verify"},
            {"name": "analysis", "weight": 1, "kind": "values"},
        ],
    }
    t = bm._coerce_task(raw)
    assert [c.name for c in t.rubric] == ["geometry", "analysis"]
    assert t.rubric[0].weight == 3
