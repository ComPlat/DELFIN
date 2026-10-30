"""HB2 red tests (package HB2, builder nacht-s9, branch agent/s9-hb2b).

The 31.10. holdout baseline ran 4 instances in 4 s and recorded every
one as ``success: false, quality: 0, error: ""`` — a dead run (the
engine never started, "engine init failed: No Anthropic API key
found") was recorded as a MEASURED failure and scored. These tests pin
the fix BEFORE it exists: run them against the unchanged commit and
watch them fail.
"""
from delfin.agent.benchmark import BenchmarkResult
from delfin.agent.bench_nightly import InstanceRecord, run as run_nightly


class _DeadRunner:
    """A run_task stub shaped like benchmark_runner.run_task that
    returns what the real one returns on engine-init failure: an
    unmeasured result carrying the error text."""

    def __call__(self, task, *, model="", **kwargs):
        self.tasks.append(task)
        self.kwargs = kwargs
        return BenchmarkResult(
            task_id=task.id, task_class=task.task_class, model=model,
            success=False, quality_0_100=0,
            error="engine init failed: No Anthropic API key found.",
            unmeasured=True,
        )

    def __init__(self):
        self.tasks = []
        self.kwargs = {}


class _HealthyRunner:
    """The control twin of _DeadRunner: a result from a run that
    actually happened."""

    def __call__(self, task, *, model="", **kwargs):
        return BenchmarkResult(
            task_id=task.id, task_class=task.task_class, model=model,
            success=True, quality_0_100=100,
        )


class _GraderSpy:
    """Grades nothing on its own; records every call so a test can
    assert an unmeasured result is never graded."""

    def __init__(self, verdict=(True, 5)):
        self._verdict = verdict
        self.seen = []

    def grade(self, task, result):
        self.seen.append(result)
        return self._verdict


class _Accounting:
    def estimate(self, task) -> float:
        return 0.1

    def charge(self, core_hours: float) -> None:
        pass


def _config(tmp_path, runner, grader, *, instances=1):
    """A holdout round config with a real seed file, a real CB0B
    sampler behind the glue, and the given runner/grader."""
    from delfin.agent import bench_holdout_glue as glue

    seed = tmp_path / "seed.txt"
    seed.write_text("48151\n")
    return {
        "model": "kit.glm-5.3-flash",
        "sampler": glue.GlueSampler(),
        "grader": grader,
        "accounting": _Accounting(),
        "runner": runner,
        "instances": instances,
        "holdout_seed_path": str(seed),
    }


def test_unmeasured_result_carries_error_text(tmp_path):
    """A result that says 'never ran' is recorded AS never having run:
    the instance record carries the error text and the unmeasured mark
    — not success=False/quality=0/error='' (the 31.10. defect)."""
    runner = _DeadRunner()
    record = run_nightly(_config(tmp_path, runner, _GraderSpy()))
    assert record.instances, "round drew no instances"
    inst = record.instances[0]
    assert inst.unmeasured is True, "unmeasured mark lost on the record"
    assert "No Anthropic API key" in inst.error, \
        f"error text lost on the record: {inst.error!r}"


def test_unmeasured_result_is_not_graded(tmp_path):
    """The grader never sees a result that says the run did not
    happen — grading a dead run is a fabricated measurement."""
    runner = _DeadRunner()
    grader = _GraderSpy()
    record = run_nightly(_config(tmp_path, runner, grader))
    assert grader.seen == [], \
        f"unmeasured result was graded: {grader.seen!r}"
    inst = record.instances[0]
    assert inst.quality is None, "unmeasured instance was scored"
    assert inst.success is None, "unmeasured instance was scored"


def test_unmeasured_not_counted_in_report_score(tmp_path):
    """A round whose single instance never ran has NO family score:
    quality stays None and n stays 0 — a dead run must not drag the
    mean down as a 0-quality failure."""
    from delfin.agent import bench_report

    runner = _DeadRunner()
    record = run_nightly(_config(tmp_path, runner, _GraderSpy()))
    report = bench_report.write(record)
    import json
    data = json.loads(report)
    families = data["families"]
    assert families, "family entries vanished entirely"
    for fam, score in families.items():
        assert score["n"] == 0, f"{fam} scored unmeasured instances"
        assert score["quality"] is None, f"{fam} has a fabricated mean"


def test_all_unmeasured_round_names_them_and_exit_nonzero(tmp_path):
    """The report names the unmeasured instances, and the CLI exits
    non-zero when EVERY drawn instance is unmeasured (a silent 0
    baseline is exactly what hid the 31.10. defect)."""
    runner = _DeadRunner()
    record = run_nightly(_config(tmp_path, runner, _GraderSpy()))
    from delfin.agent import bench_report
    report = bench_report.write(record)
    assert "engine init failed" in report, "report hides the dead runs"

    from delfin.agent import bench_nightly
    argv = ["--model", "kit.glm-5.3-flash", "--instances", "1"]
    fake_cfg = _config(tmp_path, runner, _GraderSpy())
    import unittest.mock as mock
    from delfin.agent import bench_holdout_glue
    with mock.patch.object(
            bench_holdout_glue, "build_default_config",
            return_value=fake_cfg):
        rc = bench_nightly.main(argv)
    assert rc != 0, "all-unmeasured round exited 0"


def test_measured_round_still_exits_zero(tmp_path):
    """Counter-case (QS ruling): an ordinary, measured round still
    grades normally and exits 0 — the unmeasured rule never refuses a
    healthy run."""
    runner = _HealthyRunner()
    grader = _GraderSpy()
    record = run_nightly(_config(tmp_path, runner, grader))
    assert len(grader.seen) == 1, "measured result was not graded"
    inst = record.instances[0]
    assert inst.success is True
    assert inst.quality == 5

    from delfin.agent import bench_nightly
    argv = ["--model", "kit.glm-5.3-flash", "--instances", "1"]
    import unittest.mock as mock
    from delfin.agent import bench_holdout_glue
    with mock.patch.object(
            bench_holdout_glue, "build_default_config",
            return_value=_config(tmp_path, runner, grader)):
        rc = bench_nightly.main(argv)
    assert rc == 0, f"healthy round exited {rc}"


def test_provider_reaches_run_task_for_kit_model(tmp_path):
    """D1: the provider derived for a kit.<model> reaches run_task —
    the engine must be built as kit, not as Anthropic."""
    from delfin.agent import bench_holdout_glue as glue

    runner = _DeadRunner()
    from delfin.agent import bench_holdout as bh
    task = glue.instance_to_task(bh.sample("48151", 1)[0])
    glue.GlueRunner(run_task_fn=runner)(task, model="kit.glm-5.3-flash")
    assert runner.kwargs.get("provider") == "kit", \
        f"provider never reached run_task: kwargs={runner.kwargs!r}"
