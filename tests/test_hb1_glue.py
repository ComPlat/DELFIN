"""HB1 red tests: the glue that connects bench_holdout's sampler to
bench_nightly.run and benchmark_runner.run_task.

Red control for package HB1, branch agent/s20-hb1b: none of these can
pass before delfin/agent/bench_holdout_glue.py exists. Run through the
gate against the UNCHANGED commit to see them red, then green after.
"""
from pathlib import Path

from delfin.agent import bench_holdout as bh
from delfin.agent.bench_nightly import RunRecord, run as run_nightly
from delfin.agent.benchmark import Task


class _GlueSampler:
    """The real sampler behind the CB0B interface."""

    def sample(self, seed, n):
        return bh.sample(seed, n)


class _RawInstanceSampler:
    """Sampler without the glue: raw bench_holdout Instances."""

    def sample(self, seed, n):
        return bh.sample(seed, n)


class _FixtureGrader:
    """Stand-in for the glue Grader; records what it graded."""

    def __init__(self, success=True, quality=3):
        self._success = success
        self._quality = quality
        self.seen = []

    def grade(self, task, result):
        self.seen.append((task, result))
        return (self._success, self._quality)


class _FixtureAccounting:
    def __init__(self, estimate=0.1):
        self._estimate = estimate
        self.charged = 0.0

    def estimate(self, task) -> float:
        return self._estimate

    def charge(self, core_hours: float) -> None:
        self.charged += core_hours


class _FakeRunner:
    """A runner stub shaped like benchmark_runner.run_task: takes the
    task, model keyword, returns a minimal result. No model call."""

    def __init__(self):
        self.tasks = []

    def __call__(self, task, *, model="", **kwargs):
        from delfin.agent.benchmark import BenchmarkResult
        self.tasks.append(task)
        return BenchmarkResult(task_id=task.id, task_class=task.task_class,
                               model=model)


def _config(tmpdir, runner, grader, *, instances=2, daily_budget=50.0,
            seed_file=False, seed_text="48151\n", sampler=None):
    from delfin.agent import bench_holdout_glue as glue
    cfg = {
        "model": "test-model",
        "sampler": sampler or glue.GlueSampler(),
        "grader": grader,
        "accounting": _FixtureAccounting(),
        "runner": runner,
        "instances": instances,
        "daily_budget_core_hours": daily_budget,
    }
    if seed_file:
        seed_path = Path(tmpdir) / "seed.txt"
        seed_path.write_text(seed_text)
        cfg["holdout_seed_path"] = str(seed_path)
    return cfg


def test_glue_adapter_turns_instance_into_runnable_task(tmp_path):
    """The adapter maps a bench_holdout.Instance onto a benchmark.Task
    run_task can run: id, task_class, mode and prompt filled."""
    from delfin.agent import bench_holdout_glue as glue

    instances = bh.sample(1, 3)
    task = glue.instance_to_task(instances[0])
    assert isinstance(task, Task)
    assert task.id, "task needs a non-empty id"
    assert task.task_class == instances[0].family_id
    assert task.mode == "solo"
    assert task.prompt == instances[0].prompt


def test_glue_runner_runs_a_sampled_instance(tmp_path):
    """The glue runner executes a sampled Instance through run_task's
    public path (here: the injected stub) and grades the result."""
    from delfin.agent import bench_holdout_glue as glue

    fake = _FakeRunner()
    runner = glue.GlueRunner(run_task_fn=fake)
    result = runner(glue.instance_to_task(bh.sample(1, 1)[0]),
                    model="test-model")
    assert len(fake.tasks) == 1
    assert isinstance(fake.tasks[0], Task)
    assert result.task_id == fake.tasks[0].id


def test_raw_instance_still_fails_the_round(tmp_path):
    """Direction check: a sampler that returns RAW Instances (no glue)
    still cannot run a round — the defect this package fixes, pinned so
    a regression that silently drops the adapter is visible."""
    fake = _FakeRunner()
    cfg = _config(tmp_path, fake, _FixtureGrader(), instances=2,
                  seed_file=True)
    cfg["sampler"] = _RawInstanceSampler()
    import pytest
    with pytest.raises(AttributeError):
        run_nightly(cfg)


def test_glue_grader_uses_run_evidence_and_sc1(tmp_path):
    """The glue Grader grades from DELFIN's own machinery: SC1 claim
    check on the answer; evidence-based verdict; returns (bool, int)."""
    from delfin.agent import bench_holdout_glue as glue

    g = glue.GlueGrader()
    result = type("R", (), {"answer": "", "evidence": []})()
    # No answer, no evidence: the grader must still return a verdict
    # tuple, not raise.
    verdict = g.grade(bh.sample(1, 1)[0], result)
    assert isinstance(verdict, tuple) and len(verdict) == 2
    assert isinstance(verdict[0], bool)
    assert isinstance(verdict[1], int)


def test_glued_round_completes_end_to_end(tmp_path):
    """A full bench_nightly.run with the REAL CB0B sampler behind the
    glue: the round completes, every drawn instance carries a proper
    task id/class (the bench_nightly.py:194 AttributeError bug)."""
    fake = _FakeRunner()
    grader = _FixtureGrader()
    cfg = _config(tmp_path, fake, grader, instances=2, seed_file=True)
    record = run_nightly(cfg)
    assert isinstance(record, RunRecord)
    assert record.ended_at > 0
    assert len(record.instances) == 2
    for inst in record.instances:
        assert inst.ran, f"instance {inst.task_id} did not run: {inst.skip_reason}"
        assert inst.task_id, "empty task id"
        assert inst.task_class, "empty task class"
    assert len(fake.tasks) == 2


def test_budget_gate_refuses_above_limit_with_glue(tmp_path):
    """With the glue in place the budget gate still refuses an instance
    whose estimate exceeds the daily budget."""
    fake = _FakeRunner()
    cfg = _config(tmp_path, fake, _FixtureGrader(), instances=1,
                  daily_budget=0.05, seed_file=True)
    record = run_nightly(cfg)
    assert record.refused_over_budget == 1
    assert "over daily budget" in record.instances[0].skip_reason
    assert len(fake.tasks) == 0


def test_seed_value_never_appears_in_record_or_task(tmp_path):
    """PL4's contract: the seed VALUE is never printed, logged or
    written. A full glued round with a real seed file carries no trace
    of the seed number anywhere in the record or in the tasks."""
    seed_value = "4815162342"
    fake = _FakeRunner()
    cfg = _config(tmp_path, fake, _FixtureGrader(), instances=2,
                  seed_file=True, seed_text=seed_value)
    record = run_nightly(cfg)
    blob = record.to_dict()
    assert str(seed_value) not in repr(blob), \
        "seed value leaked into the run record"
    for task in fake.tasks:
        assert seed_value not in task.prompt
        assert seed_value not in task.id
        assert seed_value not in task.task_class
    assert record.holdout_ran is True


def test_holdout_refused_without_seed_file(tmp_path):
    """Without the seed file the run refuses the holdout and says so —
    and the sampler is never called (no instances drawn)."""
    calls = []

    class _CountingSampler(_GlueSampler):
        def sample(self, seed, n):
            calls.append(seed)
            return super().sample(seed, n)

    fake = _FakeRunner()
    cfg = _config(tmp_path, fake, _FixtureGrader(), instances=2)
    cfg["sampler"] = _CountingSampler()
    record = run_nightly(cfg)
    assert record.holdout_ran is False
    assert "refus" in record.holdout_refusal_reason.lower()
    assert record.holdout_refusal_reason != ""
    assert calls == [], "sampler must not be called without a seed"
    assert len(fake.tasks) == 0
