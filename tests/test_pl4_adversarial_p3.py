"""PL4 reviewer adversarial tests (s19) — phase 3 submit-path.

Attacks the boundary between the nightly loop and the injected submitter.
Red here = a finding to nacht-s38 with the test id.
"""
import tempfile
from pathlib import Path

from delfin.agent.bench_nightly import run as run_nightly
from delfin.agent.benchmark import Task

from tests.test_pl4_budget_and_seed import (_FixtureAccounting,
                                            _FixtureGrader,
                                            _FixtureSampler,
                                            _runner_stub)


def _tasks(n=2, cls="family-a_consistency"):
    return [Task(id=f"t{i}", task_class=cls, mode="solo", prompt="do it")
            for i in range(n)]


def _config(tmpdir, tasks, **over):
    cfg = {
        "model": "test-model",
        "sampler": _FixtureSampler(tasks),
        "grader": _FixtureGrader({}),
        "accounting": _FixtureAccounting(estimate=0.5),
        "daily_budget_core_hours": 50.0,
        "instances": len(tasks),
        "holdout_seed_path": str(Path(tmpdir) / "seed"),
        "mode": "holdout",
        "runner": _runner_stub(tasks, {}),
    }
    with open(Path(tmpdir) / "seed", "w", encoding="utf-8") as fh:
        fh.write("nightly-seed-123\n")
    cfg.update(over)
    return cfg


def test_broken_submitter_does_not_kill_the_round():
    """A submitter that raises must be recorded on the instance and the
    round must complete with all instances accounted for."""
    def _boom(task, lines, name):
        raise RuntimeError("slurm down")

    with tempfile.TemporaryDirectory() as tmp:
        rec = run_nightly(_config(tmp, _tasks(3), submitter=_boom,
                                  tool_python="/opt/venv/bin/python"))
    assert rec.ended_at > 0
    assert len(rec.instances) == 3
    for inst in rec.instances:
        # The instance never executed: submit_path slurm-failed + error,
        # but ran=True is recorded (bench_nightly.py:234 sets ran before
        # the submitter is tried, and the except path leaves it True).
        assert inst.submit_path == "slurm-failed"
        assert "slurm down" in inst.error
        assert inst.ran is False  # RED: currently True


def test_broken_submitter_keeps_budget_accounting_intact():
    """A failed submission still consumed the budget: spent must include
    every instance the loop attempted to submit, and the EF1 charge must
    have fired once per attempt (no free retries)."""
    class _Acct:
        def __init__(self):
            self.n = 0
        def estimate(self, task):
            return 0.5
        def charge(self, ch):
            self.n += 1

    def _boom(task, lines, name):
        raise RuntimeError("down")

    acct = _Acct()
    with tempfile.TemporaryDirectory() as tmp:
        rec = run_nightly(_config(tmp, _tasks(3), submitter=_boom,
                                  accounting=acct,
                                  tool_python="/opt/venv/bin/python"))
    assert rec.spent_core_hours == 1.5  # 3 x 0.5, no discount for failures
    assert acct.n == 3


def test_local_runner_runs_after_failed_submission():
    """A slurm-failed instance must NOT silently fall back to the local
    runner, and its record must not claim it ran: DELFIN's submit path is
    not optional decoration."""
    def _boom(task, lines, name):
        raise RuntimeError("down")

    with tempfile.TemporaryDirectory() as tmp:
        rec = run_nightly(_config(tmp, _tasks(1), submitter=_boom,
                                  tool_python="/opt/venv/bin/python"))
    inst = rec.instances[0]
    assert inst.submit_path == "slurm-failed"
    assert inst.ran is False  # RED: currently True (record integrity)


def test_practice_mode_uses_public_seed_with_submitter():
    """Practice mode + submitter: the public seed must drive the sampler
    (holdout refusal logic must not leak into practice runs)."""
    seen = {}

    def submitter(task, lines, name):
        return "job-%s" % task.id

    with tempfile.TemporaryDirectory() as tmp:
        rec = run_nightly(_config(tmp, _tasks(2), mode="practice",
                                  public_seed="open-seed",
                                  submitter=submitter))
    assert rec.holdout_ran is False
    assert rec.mode == "practice"
    assert [i.job_id for i in rec.instances] == ["job-t0", "job-t1"]
    assert rec.seed_log == "open-seed"
