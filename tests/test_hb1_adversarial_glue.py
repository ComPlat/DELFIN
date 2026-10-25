"""HB1 adversarial review tests (nacht-s25).

End-to-end wiring the package claims: a sampled ``bench_holdout``
instance goes through ``bench_nightly.run`` and reaches a runner without
crashing; the budget gate refuses above the limit; the holdout seed
value never appears in any output. Fakes for the runner (no model call
in tests); the real ``bench_holdout.sample`` sampler.

Red on the base (the glue does not exist yet): ``bench_holdout.Instance``
carries ``instance_id`` and no ``task_class``, so ``bench_nightly.run``
crashes with AttributeError before recording anything.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "delfin" / "agent"))

import bench_holdout as bh  # noqa: E402
import bench_nightly as bn  # noqa: E402


class FakeAccounting:
    """Minimal Accounting: fixed estimate, records charges."""

    def __init__(self, est: float = 1.0):
        self.est = est
        self.charged: list[float] = []

    def estimate(self, task) -> float:
        return self.est

    def charge(self, core_hours: float) -> None:
        self.charged.append(core_hours)


class FakeRunner:
    """Records the task it received; returns a marker result."""

    def __init__(self):
        self.tasks: list = []

    def __call__(self, task, **kwargs):
        self.tasks.append(task)
        return {"ok": True}


def _fake_grader():
    def grade(task, result):
        return (True, 3)

    return type("G", (), {"grade": staticmethod(grade)})()


def test_sampled_instance_reaches_runner_without_crashing():
    """A real sampled instance must run end to end: no AttributeError on
    task.id / task.task_class, task handed to the runner."""
    runner = FakeRunner()
    cfg = {
        "sampler": bh,
        "grader": _fake_grader(),
        "accounting": FakeAccounting(),
        "runner": runner,
        "model": "fake",
        "instances": 2,
        "mode": "practice",
        "public_seed": str(bh.PUBLIC_SEED),
    }
    record = bn.run(cfg)
    assert record.instances, "no instance recorded at all"
    for inst in record.instances:
        assert inst.error == "", "instance failed: %s" % inst.error
        assert inst.ran, "instance did not run: %s" % inst.skip_reason
    assert len(runner.tasks) == 2


def test_budget_gate_refuses_above_limit():
    """With a small daily budget, later instances are refused, not run."""
    runner = FakeRunner()
    cfg = {
        "sampler": bh,
        "grader": _fake_grader(),
        "accounting": FakeAccounting(est=2.0),
        "runner": runner,
        "model": "fake",
        "instances": 5,
        "daily_budget_core_hours": 3.0,
        "mode": "practice",
        "public_seed": str(bh.PUBLIC_SEED),
    }
    record = bn.run(cfg)
    assert record.spent_core_hours <= 3.0
    assert record.refused_over_budget > 0
    assert len(runner.tasks) < 5


def test_holdout_without_seed_file_refuses_and_hides_seed():
    """Holdout mode, no seed file: refusal, zero instances, and the
    PUBLIC_SEED value must not leak into the record."""
    cfg = {
        "sampler": bh,
        "grader": _fake_grader(),
        "accounting": FakeAccounting(),
        "runner": FakeRunner(),
        "model": "fake",
        "instances": 3,
        "mode": "holdout",
        "holdout_seed_path": "",  # nothing configured
    }
    record = bn.run(cfg)
    assert not record.holdout_ran
    assert record.holdout_refusal_reason != ""
    assert record.instances == []
    assert str(bh.PUBLIC_SEED) not in record.to_dict().__str__()


def test_seed_value_never_in_record_when_practice():
    """Even in practice mode the record must never carry a secret-like
    seed: only the PUBLIC seed is logged there, by contract. A holdout
    seed file read from tmp_path must NOT appear in the record."""
    import tempfile
    import os

    fd, path = tempfile.mkstemp()
    os.close(fd)
    marker = "4242-holdout-only"
    with open(path, "w") as fh:
        fh.write(marker)
    runner = FakeRunner()
    cfg = {
        "sampler": bh,
        "grader": _fake_grader(),
        "accounting": FakeAccounting(),
        "runner": runner,
        "model": "fake",
        "instances": 2,
        "mode": "holdout",
        "holdout_seed_path": path,
    }
    try:
        record = bn.run(cfg)
        blob = record.to_dict().__str__()
        assert marker not in blob, "seed value leaked into the record"
        assert record.holdout_ran
        assert record.seed_log == "holdout seed file present"
    finally:
        os.unlink(path)
