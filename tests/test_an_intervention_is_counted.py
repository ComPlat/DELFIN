"""An intervention count: what a run cost beyond tokens.

A headless benchmark run has no user to intervene, but it can be
BLOCKED: every permission denial is the harness stepping in and steering
the agent away from what it tried.  Two agent versions can tie on
success rate and quality while one of them needed three denials to get
there -- a difference the score-card must carry, not flatten.

Counted via the engine's ``on_permission_denied`` callback and collected
by the runner through the same pass-a-callback pattern the Wave-4 tool
results use (the ``_run_once`` return contract stays five keys).  A stub
engine that never passes the callback reads 0, not None: in a benchmark
there is no "unobserved" -- the run was watched start to finish.
"""

from delfin.agent.benchmark import (
    Signal,
    Task,
    Trajectory,
    aggregate_replicates,
    score_outcome,
)


def _task() -> Task:
    return Task(id="t", task_class="c", mode="solo",
                prompt="do it",
                expected_signals=[Signal(pattern=r"done")])


def _traj(**kw) -> Trajectory:
    base = dict(text="done", duration_s=1.0, cost_usd=0.01)
    base.update(kw)
    return Trajectory(**base)


def test_score_outcome_carries_intervention_count():
    r = score_outcome(_task(), _traj(user_interventions=3))
    assert r.user_interventions == 3


def test_default_intervention_count_is_zero_not_none():
    r = score_outcome(_task(), _traj())
    # Zero, not None: the run was observed end-to-end; "unobserved" is
    # the audit-log denials question, a different axis.
    assert r.user_interventions == 0


def test_aggregate_takes_median_interventions():
    runs = [
        score_outcome(_task(), _traj(user_interventions=n))
        for n in (1, 3, 8)
    ]
    agg = aggregate_replicates(runs)
    assert agg.user_interventions == 3


def test_summarise_run_reports_total_interventions():
    from delfin.agent.benchmark import summarise_run
    rows = [score_outcome(_task(), _traj(user_interventions=n))
            for n in (2, 3)]
    rows = [r.__dict__ if hasattr(r, "__dict__") else r for r in rows]
    s = summarise_run(rows)
    assert s["total_interventions"] == 5
