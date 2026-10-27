"""HB1 phase-3 red tests: the CLI entry point
``python -m delfin.agent.bench_nightly``.

Red control: no __main__/argparse exists in bench_nightly.py (verified
this session), so every case here fails with AttributeError /
SystemExit-nonzero on the phase-2 commit. All CLI behavior is exercised
through the module's public main() (fakes injected via monkeypatched
build_default_config / run) — no model call, no real run_task.
"""
import json

import pytest

from delfin.agent import bench_nightly as bn
from delfin.agent import bench_holdout_glue as glue


class _FakeAccounting:
    def __init__(self, estimate=0.1):
        self._estimate = estimate
        self.charged = 0.0

    def estimate(self, task):
        return self._estimate

    def charge(self, core_hours):
        self.charged += core_hours


@pytest.fixture()
def fake_pipeline(monkeypatch):
    """Replace the glue config with a fully faked pipeline: sampler
    returns two recorded tasks, the runner is a stub, the grader passes
    everything. run() itself stays REAL — the CLI drives it."""
    from delfin.agent.benchmark import Task

    tasks = [Task(id=f"t{i}", task_class="classic_opt_freq", mode="solo",
                  prompt=f"p{i}") for i in range(2)]
    seen = {}

    class _Sampler:
        def sample(self, seed, n):
            seen["seed"] = seed
            seen["n"] = n
            return tasks[:n]

    class _Grader:
        def grade(self, task, result):
            return (True, 4)

    def fake_runner(task, *, model="", **kw):
        from delfin.agent.benchmark import BenchmarkResult
        return BenchmarkResult(task_id=task.id,
                               task_class=task.task_class, model=model)

    def fake_build(**kwargs):
        cfg = {
            "model": kwargs["model"],
            "sampler": _Sampler(),
            "grader": _Grader(),
            "accounting": _FakeAccounting(),
            "runner": fake_runner,
            "instances": kwargs["instances"],
            "mode": kwargs["mode"],
            "holdout_seed_path": kwargs.get("holdout_seed_path", ""),
            "public_seed": kwargs.get("public_seed", 7),
        }
        seen["config"] = cfg
        return cfg

    monkeypatch.setattr(glue, "build_default_config", fake_build)
    return seen


def _run_main(capsys, argv):
    rc = bn.main(argv)
    out = capsys.readouterr().out
    return rc, out


def test_main_without_seed_file_refuses(capsys, fake_pipeline, tmp_path,
                                        monkeypatch):
    """Holdout without the seed file: refuse, say so on stderr, exit
    non-zero; the sampler is never called."""
    rc, out = _run_main(capsys, ["--model", "M", "--instances", "2",
                                 "--seed-path", str(tmp_path / "none.txt")])
    assert rc != 0
    assert "refus" in out.lower()
    assert "seed" not in fake_pipeline, "sampler must not run"
    # Nothing was graded: the record shows zero instances.
    assert _record(out)["instances"] == []


def _record(out):
    """The CLI's stdout IS the record (pretty-printed JSON)."""
    return json.loads(out)


def test_main_holdout_with_seed_file_runs(capsys, fake_pipeline, tmp_path):
    seed_file = tmp_path / "seed.txt"
    seed_file.write_text("424242\n")
    rc, out = _run_main(capsys, ["--model", "M", "--instances", "2",
                                 "--seed-path", str(seed_file)])
    assert rc == 0
    record = _record(out)
    assert record["holdout_ran"] is True
    assert record["mode"] == "holdout"
    # Two instances were drawn and ran through the (fake) pipeline.
    assert len(record["instances"]) == 2
    assert all(i["ran"] for i in record["instances"])
    # The seed VALUE never reaches stdout.
    assert "424242" not in out


def test_main_practice_mode_needs_no_seed_file(capsys, fake_pipeline):
    rc, out = _run_main(capsys, ["--model", "M", "--instances", "2",
                                 "--practice"])
    assert rc == 0
    record = _record(out)
    assert record["mode"] == "practice"
    assert record["holdout_ran"] is False
    assert len(record["instances"]) == 2


def test_main_practice_does_not_read_seed_path(capsys, fake_pipeline,
                                               tmp_path):
    """--practice ignores --seed-path entirely (PL4 contract: practice
    runs use the public seed, never the holdout seed file)."""
    seed_file = tmp_path / "seed.txt"
    seed_file.write_text("424242\n")
    rc, out = _run_main(capsys, ["--model", "M", "--instances", "1",
                                 "--seed-path", str(seed_file),
                                 "--practice"])
    assert rc == 0
    assert "424242" not in out


def test_main_requires_model(capsys, fake_pipeline):
    rc = bn.main(["--instances", "2", "--practice"])
    assert rc != 0
