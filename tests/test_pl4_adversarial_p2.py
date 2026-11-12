"""PL4 reviewer adversarial tests (s19) — P1 fix @ 31944f81 + P2 @ 14979198.

Attacks the fix commit and the report writer. Red here = a finding to
nacht-s38 with the test id.
"""
import json
import tempfile
from pathlib import Path

import pytest

from delfin.agent.bench_nightly import (InstanceRecord, RunRecord,
                                        run as run_nightly)
from delfin.agent.bench_report import write as write_report
from delfin.agent.benchmark import Task

from tests.test_pl4_budget_and_seed import (_FixtureAccounting,
                                            _FixtureGrader,
                                            _FixtureSampler,
                                            _runner_stub)


def _tasks(n=3, cls="family-a_consistency"):
    return [Task(id=f"t{i}", task_class=cls, mode="solo", prompt="do it")
            for i in range(n)]


def _config(tmpdir, tasks, *, accounting=None, seed_file=True,
            seed_text="nightly-seed-123\n", scores=None, **over):
    with open(Path(tmpdir) / "seed", "w", encoding="utf-8") as fh:
        if seed_file:
            fh.write(seed_text)
    if accounting is None:
        accounting = _FixtureAccounting(estimate=0.0)
    cfg = {
        "model": "test-model",
        "runner": _runner_stub(tasks, scores or {}),
        "sampler": _FixtureSampler(tasks),
        "grader": _FixtureGrader(scores or {}),
        "accounting": accounting,
        "daily_budget_core_hours": 50.0,
        "instances": len(tasks),
        "holdout_seed_path": str(Path(tmpdir) / "seed"),
    }
    cfg.update(over)
    return cfg


# --- attacks on the F1-F4 fix ------------------------------------------------


class _CountingAccounting:
    """Charges ledger-style so a double charge is measurable."""

    def __init__(self, estimate=0.6):
        self._est = estimate
        self.ledger = 0.0
        self.charges = 0

    def estimate(self, task):
        return self._est

    def charge(self, core_hours):
        self.ledger += core_hours
        self.charges += 1


def test_double_charge_of_started_instances():
    """charge() must fire once per started instance. On 31944f81 it fired
    twice (bench_nightly.py:210 inside the try, again at :219). Assert the
    CALL COUNT, not a float sum: 3 x 0.6 accumulates 1.7999999999999998."""
    acct = _CountingAccounting(estimate=0.6)
    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(3), accounting=acct))
    assert record.instances[0].ran
    # Three started instances, exactly three charges — on 31944f81: 6.
    assert acct.charges == 3
    assert acct.ledger == pytest.approx(1.8)


def test_spent_equals_ledger_after_round():
    """record.spent_core_hours and EF1's ledger must agree after a round."""
    acct = _CountingAccounting(estimate=0.6)
    with tempfile.TemporaryDirectory() as second:
        record = run_nightly(_config(second, _tasks(3), accounting=acct))
    assert record.spent_core_hours == pytest.approx(acct.ledger)


def test_seed_log_never_carries_the_seed_value_in_fix_path():
    """The fix path (present, non-empty seed file) must log presence only,
    never the seed value ('holdout seed file present')."""
    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks()))
    assert record.seed_log == "holdout seed file present"


# --- attacks on bench_report.write -------------------------------------------


def _record_with_instances():
    from delfin.agent.bench_nightly import RunRecord
    rec = RunRecord(model="m", mode="holdout", spent_core_hours=1.8,
                    daily_budget_core_hours=50.0, holdout_ran=True,
                    seed_log="holdout seed file present")
    rec.instances = [
        InstanceRecord(task_id="t0", task_class="fam_a_consistency",
                       ran=True, success=True, quality=80, core_hours=0.6),
        InstanceRecord(task_id="t1", task_class="fam_a_consistency",
                       ran=True, success=False, quality=40, core_hours=0.6),
        InstanceRecord(task_id="t2", task_class="fam_b_trap",
                       ran=False, skip_reason="over daily budget"),
    ]
    return rec


def test_report_reproducible_from_the_record():
    """Same record + history -> byte-identical report (target 5)."""
    rec = _record_with_instances()
    assert write_report(rec) == write_report(rec)


def test_report_never_contains_seed_value():
    """The seed value must not enter the report even via seed_log."""
    rec = _record_with_instances()
    rec.seed_log = "holdout seed file present"
    text = write_report(rec, history=[])
    assert "nightly-seed-123" not in text
    assert "seed-path" not in text


def test_stop_verdict_after_two_flat_rounds():
    """Two consecutive rounds without gain -> should_stop=True."""
    from delfin.agent.bench_nightly import RunRecord
    prev = RunRecord(previous_rounds_gain=[0.0])
    rec = RunRecord(previous_rounds_gain=[0.0])
    verdict = json.loads(write_report(rec, history=[prev]))
    assert verdict["stop"]["should_stop"] is True
    assert "gain" in verdict["stop"]["reason"]


def test_stop_not_triggered_by_one_flat_round():
    """One flat round alone must not stop the loop."""
    rec = RunRecord(previous_rounds_gain=[0.0])
    verdict = json.loads(write_report(rec, history=[]))
    assert verdict["stop"]["should_stop"] is False


def test_stop_on_budget_exhausted():
    rec = RunRecord(spent_core_hours=50.0, daily_budget_core_hours=50.0)
    verdict = json.loads(write_report(rec, history=[]))
    assert verdict["stop"]["should_stop"] is True
    assert "budget" in verdict["stop"]["reason"]


def test_weakest_queue_ordered_weakest_first():
    """t0 quality 80 is >= the weakest threshold (60): excluded from the
    queue — only t1 (40) may appear."""
    rec = _record_with_instances()
    report = json.loads(write_report(rec))
    q = report["weakest_queue"]
    assert [e["task_id"] for e in q] == ["t1"]
    assert q[0]["quality"] == 40
