"""PL4 phase 2 red test: bench_report.write(record, history).

Pins the PL4 brief target 4: score per tool family and per correctness
level (consistency, reproduction, reference, trap), cost next to it, the
weakest tasks as a queue for the next round, and the stop verdict
(two rounds without gain, budget exhausted).
"""
import json
import tempfile

from delfin.agent.bench_nightly import InstanceRecord, RunRecord
from delfin.agent.bench_report import write as write_report


def _record(spent=3.0, budget=50.0, gain=0.0, quality_by_class=None):
    quality_by_class = quality_by_class or {"family-a": (False, 40)}
    rec = RunRecord(model="m", mode="holdout",
                    daily_budget_core_hours=budget,
                    spent_core_hours=spent,
                    previous_rounds_gain=[gain])
    for i, (fam, (ok, q)) in enumerate(sorted(quality_by_class.items())):
        rec.instances.append(InstanceRecord(
            task_id=f"t{i}", task_class=fam, ran=True,
            success=ok, quality=q, core_hours=1.0))
    return rec


def test_report_reproducible_from_record():
    """write(record, history) -> text that parses as JSON (reproducible
    from the record: the same record yields the same bytes)."""
    rec = _record()
    text = write_report(rec, history=[])
    again = write_report(rec, history=[])
    assert json.loads(text) == json.loads(again)
    report = json.loads(text)
    assert report["families"]["family-a"]["quality"] == 40
    assert report["families"]["family-a"]["success"] is False


def test_report_scores_per_correctness_level():
    """Correctness levels come from the task_class suffix: _consistency,
    _reproduction, _reference, _trap — each scored separately."""
    rec = _record(quality_by_class={
        "family-a_consistency": (True, 80),
        "family-a_reproduction": (False, 30),
        "family-b_reference": (True, 70),
        "family-b_trap": (False, 10),
    })
    report = json.loads(write_report(rec, history=[]))
    levels = report["correctness_levels"]
    assert levels["consistency"]["quality"] == 80
    assert levels["reproduction"]["quality"] == 30
    assert levels["reference"]["quality"] == 70
    assert levels["trap"]["quality"] == 10


def test_report_cost_next_to_score():
    """Cost sits beside each score: core-hours spent vs budget."""
    rec = _record(spent=3.0, budget=50.0)
    report = json.loads(write_report(rec, history=[]))
    assert report["cost"]["spent_core_hours"] == 3.0
    assert report["cost"]["daily_budget_core_hours"] == 50.0
    assert report["families"]["family-a"]["core_hours"] == 1.0


def test_report_weakest_tasks_queue():
    """The lowest-quality tasks come back as the next round's queue,
    weakest first."""
    rec = _record(quality_by_class={
        "good_consistency": (True, 90),
        "mid_consistency": (True, 60),
        "bad_trap": (False, 10),
        "worse_reproduction": (False, 5),
    })
    report = json.loads(write_report(rec, history=[]))
    queue = report["weakest_queue"]
    ids = [t["task_id"] for t in queue]
    assert ids[0] == "worse_reproduction"
    assert ids[1] == "bad_trap"
    # strong tasks are not queued
    assert "good_consistency" not in ids


def test_report_stop_verdict_two_flat_rounds():
    """Two rounds without gain -> stop verdict (quality not improved)."""
    rec = _record(gain=0.0)
    history = [_record(gain=0.0)]  # one earlier round, also flat
    report = json.loads(write_report(rec, history=history))
    assert report["stop"]["should_stop"] is True
    assert "no gain" in report["stop"]["reason"]


def test_report_no_stop_when_gaining_or_mid_run():
    """One flat round only, or rounds WITH gain -> no stop. Budget
    exhausted stops for a different reason."""
    rec = _record(gain=0.0)
    report = json.loads(write_report(rec, history=[]))
    assert report["stop"]["should_stop"] is False  # only one flat round

    rec2 = _record(gain=5.0)
    report2 = json.loads(write_report(rec2, history=[_record(gain=0.0)]))
    assert report2["stop"]["should_stop"] is False

    rec3 = _record(gain=5.0, spent=50.0, budget=50.0)
    report3 = json.loads(write_report(rec3, history=[]))
    assert report3["stop"]["should_stop"] is True
    assert "budget" in report3["stop"]["reason"]
