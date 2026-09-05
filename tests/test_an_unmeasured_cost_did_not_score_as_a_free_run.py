"""A cost nobody measured was scored as a cost of zero.

Two places treated a missing number as a good number.

``benchmark.score_outcome`` awards up to 10 "cost_pts" for staying under
``task.max_cost_usd``. The guard was::

    if task.max_cost_usd > 0 and traj.cost_usd >= 0:

No cost can fail ``>= 0``. So a trajectory whose ``cost_usd`` was 0.0 --
which, given the pricing bugs this change fixes, is exactly what a model
with no known rate produced -- scored ``10/10``: a perfect mark for thrift
it never demonstrated. Measured: ``Trajectory(cost_usd=0.0)`` scored
``cost_pts == 10``, indistinguishable from a run that genuinely cost
nothing.

The distinction that was missing is between a *measured* zero and an
*absent* measurement. A KIT-Toolbox or Ollama run really does cost $0.00
and has earned its 10 points. A run on a model with no rate on record has
not been measured at all, and now scores 0 with the reason recorded in
``budget_violations`` so the missing points are explainable rather than
mysterious.

``briefing._analyse_outcomes`` had the mirror-image of the same habit::

    costs = [o.cost_usd for o in class_outcomes if o.cost_usd > 0]
    avg_cost = sum(costs) / len(costs) if costs else 0.0

...printed as ``Avg cost: $X/task``. Outcomes with no recorded cost are
dropped from the denominator, so with three $0.30 tasks and seven that
reported nothing, the briefing told the agent a task costs $0.30 when the
recorded evidence covers three tasks out of ten. It is a conditional mean
labelled as an average. It still is a conditional mean -- there is no
honest way to invent the missing seven -- but it now says so.
"""

from __future__ import annotations

from delfin.agent import benchmark as bm
from delfin.agent import briefing
from delfin.agent.outcome_tracker import CycleOutcome


def _task(**kw):
    base = dict(
        id="t", task_class="c", mode="solo", prompt="p",
        expected_signals=(bm.Signal(pattern="ok", against="text"),),
        max_duration_s=10.0, max_cost_usd=0.10, max_tool_calls=3,
    )
    base.update(kw)
    return bm.Task(**base)


# ---------------------------------------------------------------------------
# Benchmark: an unmeasured run is not a cheap run
# ---------------------------------------------------------------------------


def test_a_run_with_no_known_rate_does_not_collect_full_cost_marks():
    task = _task()
    traj = bm.Trajectory(text="ok", duration_s=1.0, cost_usd=0.0)
    r = bm.score_outcome(task, traj, model="opus")
    assert r.components["cost_pts"] == 0


def test_the_missing_cost_marks_come_with_a_reason():
    task = _task()
    traj = bm.Trajectory(text="ok", duration_s=1.0, cost_usd=0.0)
    r = bm.score_outcome(task, traj, model="opus")
    assert any(v.startswith("cost_usd unmeasured") for v in r.budget_violations)


def test_a_genuinely_free_provider_still_earns_its_cost_marks():
    """A quota-funded or local model's zero is a real zero."""
    task = _task()
    traj = bm.Trajectory(text="ok", duration_s=1.0, cost_usd=0.0)
    for model in ("kit.qwen3.5-397b-A17b", "qwen3-coder:32b"):
        r = bm.score_outcome(task, traj, model=model)
        assert r.components["cost_pts"] == 10, model
        assert r.budget_violations == [], model


def test_observed_spend_is_scored_whatever_the_table_knows():
    """The runner's own measurement outranks the price table.

    The CLI backend reports the provider's ``total_cost_usd`` per turn, so
    a positive cost is measured even for an id with no rate on record.
    """
    task = _task(max_cost_usd=0.10)
    cheap = bm.score_outcome(
        task, bm.Trajectory(text="ok", cost_usd=0.01), model="opus")
    dear = bm.score_outcome(
        task, bm.Trajectory(text="ok", cost_usd=0.09), model="opus")
    assert cheap.components["cost_pts"] > dear.components["cost_pts"]
    assert not any(v.startswith("cost_usd unmeasured")
                   for v in cheap.budget_violations)


def test_an_unmeasured_run_cannot_outscore_a_measured_cheap_one():
    """The whole point: unknown must not beat known-and-thrifty."""
    task = _task()
    unmeasured = bm.score_outcome(
        task, bm.Trajectory(text="ok", cost_usd=0.0), model="opus")
    measured = bm.score_outcome(
        task, bm.Trajectory(text="ok", cost_usd=0.001),
        model="claude-opus-4-20250514")
    assert measured.quality_0_100 > unmeasured.quality_0_100


def test_quality_is_still_the_sum_of_its_components():
    task = _task()
    r = bm.score_outcome(task, bm.Trajectory(text="ok"), model="opus")
    assert r.quality_0_100 == sum(r.components.values())


# ---------------------------------------------------------------------------
# Briefing: a conditional mean says which condition
# ---------------------------------------------------------------------------


def _outcomes(n_priced: int, n_unpriced: int) -> list[CycleOutcome]:
    out = [
        CycleOutcome(provider="claude", task_class="code", mode="solo",
                     verdict="PASS", cost_usd=0.30)
        for _ in range(n_priced)
    ]
    out += [
        CycleOutcome(provider="claude", task_class="code", mode="solo",
                     verdict="PASS", cost_usd=0.0)
        for _ in range(n_unpriced)
    ]
    return out


def test_the_briefing_names_the_denominator_when_it_is_not_every_task():
    stats = briefing._analyse_outcomes(_outcomes(3, 7), "code")
    assert stats["cost_n"] == 3
    assert stats["total"] == 10
    line = _cost_line(stats)
    assert "3/10" in line, line


def test_a_fully_priced_history_reads_as_a_plain_average():
    stats = briefing._analyse_outcomes(_outcomes(4, 0), "code")
    assert stats["cost_n"] == stats["total"] == 4
    assert _cost_line(stats) == "Avg cost: $0.30/task"


def _cost_line(stats: dict) -> str:
    """Rebuild the one briefing line under test from the same stats."""
    n, total = stats.get("cost_n", 0), stats.get("total", 0)
    over = "" if n >= total else f" (of {n}/{total} with a recorded cost)"
    return f"Avg cost: ${stats['avg_cost']:.2f}/task{over}"


def test_the_generated_briefing_carries_the_qualifier(tmp_path):
    path = tmp_path / "outcomes.jsonl"
    from delfin.agent import outcome_tracker as ot
    for o in _outcomes(3, 7):
        ot.append_outcome(o, path=path)
    text = briefing.generate_briefing("claude", "fix the parser",
                                      outcome_path=path)
    assert "Avg cost" in text, text
    assert "with a recorded cost" in text, text
