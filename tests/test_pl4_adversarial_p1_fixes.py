"""PL4 P1 fixes: red tests for reviewer s19's adversarial findings F1-F4.

F1/F2: a broken accounting call must not kill the round.
F3: an empty seed file must refuse the holdout.
F4: a negative estimate must not let later instances fit for free.
"""
import tempfile
from pathlib import Path

from delfin.agent.bench_nightly import RunRecord, run as run_nightly
from delfin.agent.benchmark import Task

from test_pl4_budget_and_seed import (
    _FixtureGrader, _FixtureSampler, _config, _tasks)


class _BrokenEstimate:
    def estimate(self, task):
        raise ValueError("accounting backend is down")

    def charge(self, core_hours):
        pass


class _BrokenCharge:
    def estimate(self, task):
        return 0.1

    def charge(self, core_hours):
        raise ValueError("ledger write failed")


def test_broken_estimate_does_not_kill_the_round():
    """F1: estimate() raising is recorded on the instance; the round
    completes and the record survives."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks())
        cfg["accounting"] = _BrokenEstimate()
        record = run_nightly(cfg)
        assert isinstance(record, RunRecord)
        assert len(record.instances) == 3
        for inst in record.instances:
            assert inst.ran is False
            assert "ValueError" in inst.skip_reason


def test_broken_charge_does_not_kill_the_round():
    """F2: charge() raising is recorded on the instance; the round
    completes."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks())
        cfg["accounting"] = _BrokenCharge()
        record = run_nightly(cfg)
        assert isinstance(record, RunRecord)
        assert len(record.instances) == 3
        for inst in record.instances:
            assert "ValueError" in inst.skip_reason


def test_empty_seed_file_refuses_holdout():
    """F3: a seed file that exists but is EMPTY refuses the holdout —
    an empty string is not a seed."""
    with tempfile.TemporaryDirectory() as tmp:
        seed_path = Path(tmp) / "seed.txt"
        seed_path.write_text("   \n")
        cfg = _config(tmp, _tasks())
        cfg["holdout_seed_path"] = str(seed_path)
        record = run_nightly(cfg)
        assert record.holdout_ran is False
        assert record.holdout_refusal_reason != ""


def test_negative_estimate_never_grows_the_budget():
    """F4: a negative estimate is refused (stated reason) and never
    reduces the spent counter, so later instances do not fit for free."""
    with tempfile.TemporaryDirectory() as tmp:
        tasks = _tasks()
        cfg = _config(tmp, tasks, daily_budget=0.5, estimate=0.0)
        ests = {"t0": -1.0, "t1": 0.3, "t2": 0.3}

        class _PerTaskEst:
            def estimate(self, task):
                return ests[task.id]

            def charge(self, core_hours):
                pass

        cfg["accounting"] = _PerTaskEst()
        record = run_nightly(cfg)
        by_id = {i.task_id: i for i in record.instances}
        assert by_id["t0"].ran is False
        assert "negative" in by_id["t0"].skip_reason
        # the budget available for t1/t2 is NOT enlarged by t0's -1.0
        assert by_id["t1"].ran is True
        assert by_id["t2"].ran is False  # 0.3 + 0.3 > 0.5
        assert record.spent_core_hours >= 0.0
