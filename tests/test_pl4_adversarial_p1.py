"""PL4 reviewer adversarial tests (s19) — phase 1 bench_nightly.run.

Findings-first: each test below attacks one claim the builder's own tests
do not cover. Red here = a finding to nacht-s38 with the test id.
"""
import tempfile
from pathlib import Path

from delfin.agent.bench_nightly import run as run_nightly
from delfin.agent.benchmark import Task

from tests.test_pl4_budget_and_seed import (_FixtureAccounting,
                                            _FixtureGrader,
                                            _FixtureSampler)


def _tasks(n=3, cls="family-a"):
    return [Task(id=f"t{i}", task_class=cls, mode="solo", prompt="do it")
            for i in range(n)]


def _config(tmpdir, tasks, *, accounting=None, seed_file=True,
            seed_text="nightly-seed-123\n"):
    from tests.test_pl4_budget_and_seed import _runner_stub

    cfg = {
        "model": "test-model",
        "sampler": _FixtureSampler(tasks),
        "grader": _FixtureGrader({}),
        "accounting": accounting or _FixtureAccounting(estimate=0.1),
        "daily_budget_core_hours": 50.0,
        "runner": _runner_stub(tasks, {}),
        "instances": len(tasks),
    }
    if seed_file:
        seed_path = Path(tmpdir) / "holdout_seed.txt"
        seed_path.write_text(seed_text)
        cfg["holdout_seed_path"] = str(seed_path)
    return cfg


def test_broken_estimate_does_not_kill_the_round():
    """A raising accounting.estimate must be recorded on the instance,
    not crash the whole nightly run (bench_nightly.py:175 is outside the
    per-instance try)."""
    class _BrokenEstimate:
        def estimate(self, task):
            raise ValueError("accounting exploded")
        def charge(self, core_hours):
            pass

    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(),
                                     accounting=_BrokenEstimate()))
    assert record.ended_at > 0
    assert len(record.instances) == 3
    # Recorded on the instance (skip_reason per the fix's contract; 'error'
    # stays reserved for a failed run itself), visible to the report.
    for inst in record.instances:
        assert inst.ran is False
        assert "accounting" in inst.skip_reason.lower()
        assert "exploded" in inst.skip_reason.lower()


def test_broken_charge_does_not_kill_the_round():
    """A raising accounting.charge must be recorded, not crash the round
    (bench_nightly.py:186 is outside the per-instance try)."""
    class _BrokenCharge:
        def estimate(self, task):
            return 0.1
        def charge(self, core_hours):
            raise ValueError("charge exploded")

    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(),
                                     accounting=_BrokenCharge()))
    assert record.ended_at > 0
    assert len(record.instances) == 3
    for inst in record.instances:
        assert inst.ran is False
        assert "charge" in inst.skip_reason.lower()
        assert "exploded" in inst.skip_reason.lower()


def test_empty_seed_file_refuses_holdout():
    """A seed file that exists but is EMPTY (or whitespace) is not a
    usable holdout seed: the gate is on file presence, not content
    (bench_nightly.py:156 checks is_file, :167 checks 'seed is None' —
    an empty string passes)."""
    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(), seed_text="  \n"))
        assert record.holdout_ran is False
        assert record.holdout_refusal_reason != ""


def test_no_holdout_seed_path_key_at_all():
    """A config without the holdout_seed_path key must refuse the holdout
    with a stated reason — not crash and not fall back to a seed."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks())
        cfg.pop("holdout_seed_path")
        record = run_nightly(cfg)
        assert record.holdout_ran is False
        assert record.holdout_refusal_reason != ""
        assert all(not i.ran for i in record.instances)


def test_negative_estimate_never_grows_the_budget():
    """A negative estimate must not INCREASE the day's allowance: with a
    negative first estimate the spent counter drops below zero and every
    later instance fits for free (bench_nightly.py:178-185)."""
    class _NegativeAccounting:
        def __init__(self):
            self.calls = 0
        def estimate(self, task):
            self.calls += 1
            return -1.0 if self.calls == 1 else 0.6
        def charge(self, core_hours):
            pass

    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(),
                                     accounting=_NegativeAccounting()))
        ran = [i.task_id for i in record.instances if i.ran]
        # budget 50: only t1 (0.6) and t2 (1.2 total) fit; t0's negative
        # estimate must not let all three run for free.
        assert ran == ["t1", "t2"], ran
        assert record.spent_core_hours >= 0.0


def test_practice_mode_wins_even_with_seed_file_present():
    """mode='practice' with a seed file PRESENT stays a practice run:
    holdout_ran stays False and the seed_log names the public seed."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks())
        cfg["mode"] = "practice"
        cfg["public_seed"] = "public-seed"
        record = run_nightly(cfg)
        assert record.holdout_ran is False
        assert record.seed_log == "public-seed"
        assert record.holdout_refusal_reason == ""


def test_skip_reason_carries_no_seed():
    """The budget skip_reason is composed of numbers only — no seed
    value, no seed path content ever appears in it."""
    seed_value = "s3cr3t-seed-xyz"
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), seed_text=seed_value + "\n")
        cfg["daily_budget_core_hours"] = 0.0
        record = run_nightly(cfg)
        blob = str(record.to_dict())
        assert seed_value not in blob
        assert all("s3cr3t" not in i.skip_reason
                   for i in record.instances)


def test_huge_estimate_is_refused_not_started():
    """An absurd estimate (1e9 core-hours) is refused, the run continues
    with the next instance, and the record reports it."""
    from tests.test_pl4_budget_and_seed import _FixtureAccounting
    class _Huge(_FixtureAccounting):
        def estimate(self, task):
            return 1e9
    with tempfile.TemporaryDirectory() as tmp:
        record = run_nightly(_config(tmp, _tasks(), accounting=_Huge()))
        assert record.refused_over_budget == 3
        assert record.spent_core_hours == 0.0
        assert all(not i.ran for i in record.instances)
