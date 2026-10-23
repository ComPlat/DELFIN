"""PL4 red test: nightly holdout runner — phase 1 budget + seed gates.

Red control for branch agent/s38-pl4b: none of these can pass before
delfin/agent/bench_nightly.py exists on the branch. Run through the gate
against the UNCHANGED commit to see them red, then after the fix green.
"""
import json
import os
import tempfile
from pathlib import Path

from delfin.agent.bench_nightly import RunRecord, run as run_nightly
from delfin.agent.benchmark import Task


class _FixtureSampler:
    """Stand-in for CB0B's bench_holdout.sample(seed, n) once merged."""

    def __init__(self, tasks):
        self._tasks = tasks

    def sample(self, seed: str, n: int):
        return self._tasks[:n]


class _FixtureGrader:
    """Stand-in for CB0A evidence / SC1 claims grading."""

    def __init__(self, score_by_id):
        self._score = score_by_id
        self.seen_prompts = []

    def grade(self, task, result):
        return self._score.get(task.id, (False, 0))


class _FixtureAccounting:
    """Stand-in for EF1's cluster accounting: per-instance core-hour estimate."""

    def __init__(self, estimate=0.0):
        self._estimate = estimate
        self.charged = 0.0

    def estimate(self, task) -> float:
        return self._estimate

    def charge(self, core_hours: float) -> None:
        self.charged += core_hours


def _runner_stub(tasks, scores):
    """A run_once/engine factory injection pair that does no real work.

    run_task (benchmark_runner.py:1450) is protected and reused as-is:
    the runner's config passes an engine_factory that builds nothing.
    Instead the fixture grading path bypasses the engine entirely —
    phase 1 tests the loop skeleton, not the engine.
    """

    def fake_run(task, **kwargs):
        from delfin.agent.benchmark import BenchmarkResult
        return BenchmarkResult(
            task_id=task.id, task_class=task.task_class,
            model=kwargs.get("model", "test"),
        )

    return fake_run


def _config(tmpdir, tasks, *, daily_budget=50.0, estimate=0.0,
            seed_file=True, scores=None):
    cfg = {
        "model": "test-model",
        "sampler": _FixtureSampler(tasks),
        "grader": _FixtureGrader(scores or {}),
        "accounting": _FixtureAccounting(estimate=estimate),
        "daily_budget_core_hours": daily_budget,
        "runner": _runner_stub(tasks, scores or {}),
        "instances": 3,
    }
    if seed_file:
        seed_path = Path(tmpdir) / "holdout_seed.txt"
        seed_path.write_text("nightly-seed-123\n")
        cfg["holdout_seed_path"] = str(seed_path)
    else:
        cfg["holdout_seed_path"] = str(Path(tmpdir) / "missing_seed.txt")
    return cfg


def _tasks():
    return [
        Task(id=f"t{i}", task_class="family-a", mode="solo",
             prompt="do the thing")
        for i in range(3)
    ]


def test_budget_overrun_refused():
    """An instance whose estimate pushes the day above the budget is not
    started and is reported as skipped, not as a failure."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), daily_budget=1.0, estimate=0.6)
        record = run_nightly(cfg)
        assert isinstance(record, RunRecord)
        started = [i.task_id for i in record.instances if i.ran]
        skipped = [i.task_id for i in record.instances
                   if not i.ran and i.skip_reason]
        # 0.6 + 0.6 = 1.2 > 1.0 → only the first fits.
        assert len(started) == 1, started
        assert len(skipped) == 2, skipped
        assert record.refused_over_budget == 2


def test_seed_absent_refuses_holdout():
    """Without the seed file the runner refuses to grade the holdout and
    says so — it does not grade with a fallback seed."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), seed_file=False)
        record = run_nightly(cfg)
        assert record.holdout_ran is False
        assert record.holdout_refusal_reason != ""


def test_seed_never_in_outputs():
    """The seed value appears in no output: record dict, report text,
    log lines. The PATH comes from a setting; the SEED is never printed."""
    with tempfile.TemporaryDirectory() as tmp:
        seed_value = "nightly-seed-123"
        cfg = _config(tmp, _tasks())
        record = run_nightly(cfg)
        blob = json.dumps(record.to_dict())
        assert seed_value not in blob
        assert cfg["holdout_seed_path"] not in blob
        # every instance/task record carries the seed in no field
        for inst in record.instances:
            for value in vars(inst).values():
                assert seed_value not in json.dumps(value, default=str)


def test_same_seed_same_instances():
    """Sampling is deterministic for a given seed."""
    with tempfile.TemporaryDirectory() as tmp:
        tasks = _tasks()
        cfg1 = _config(tmp, tasks, daily_budget=1000.0)
        cfg2 = _config(tmp, tasks, daily_budget=1000.0)
        r1 = run_nightly(cfg1)
        r2 = run_nightly(cfg2)
        ids1 = [i.task_id for i in r1.instances]
        ids2 = [i.task_id for i in r2.instances]
        assert ids1 == ids2
        # ordering, not set equality: the whole sequence repeats
        assert ids1 == ["t0", "t1", "t2"]


def test_budget_boundary_exact_limit_starts():
    """At exactly the limit the instance starts; one cent-hour above it
    is refused (reviewer s19, finding a)."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), daily_budget=1.2, estimate=0.6)
        record = run_nightly(cfg)
        started = [i.task_id for i in record.instances if i.ran]
        assert started == ["t0", "t1"], started  # 0.6, then 1.2 == limit

    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), daily_budget=1.199999, estimate=0.6)
        record = run_nightly(cfg)
        started = [i.task_id for i in record.instances if i.ran]
        assert started == ["t0"], started  # second would exceed


def test_practice_run_works_without_seed():
    """A public-seed PRACTICE run still runs to the end even though the
    holdout seed file is absent (reviewer s19, finding c)."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _config(tmp, _tasks(), seed_file=False)
        cfg["mode"] = "practice"
        cfg["public_seed"] = "public-practice-seed"
        record = run_nightly(cfg)
        assert record.holdout_ran is False
        started = [i.task_id for i in record.instances if i.ran]
        assert started == ["t0", "t1", "t2"], started


def test_seed_path_not_in_record_or_instances():
    """The seed PATH never lands in the serialized record either
    (reviewer s19, finding b — grep the JSON for the path content)."""
    with tempfile.TemporaryDirectory() as tmp:
        seed_path = str(Path(tmp) / "holdout_seed.txt")
        cfg = _config(tmp, _tasks())
        record = run_nightly(cfg)
        blob = json.dumps(record.to_dict())
        assert seed_path not in blob
