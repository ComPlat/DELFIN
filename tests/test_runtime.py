"""Tests for the runtime layer: async runs, run store, events, cancel."""

from __future__ import annotations

import time

from delfin.tools import Application, OutputSpec, Pipeline
from delfin.tools._application import register_application
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._runtime import RunStatus, RunStore, Runtime
from delfin.tools._types import StepStatus


class _FastEnergy(StepAdapter):
    name = "rt_fast_energy"
    description = "fast fake energy for runtime tests"
    produces_geometry = False

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, time.monotonic(),
            data={"energy_Eh": -2.0},
        )


class _SlowStep(StepAdapter):
    name = "rt_slow_step"
    description = "slow fake step for runtime cancel tests"
    produces_geometry = True

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        time.sleep(0.2)
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, time.monotonic())


register(_FastEnergy())
register(_SlowStep())


def _energy_app():
    pipe = Pipeline("rt_app_pipe")
    pipe.add("rt_fast_energy", label="calc")
    return Application.from_pipeline(
        pipe, name="rt_app",
        outputs=(OutputSpec("energy_Eh", step="calc", key="energy_Eh", unit="Eh"),),
    )


def test_submit_wait_success(tmp_path):
    rt = Runtime(RunStore(tmp_path / "store"))
    register_application(_energy_app())
    h = rt.submit_application("rt_app", work_dir=tmp_path / "wd")
    rec = h.wait(timeout=10)
    assert rec.status == RunStatus.SUCCESS.value
    assert rec.outputs == {"energy_Eh": -2.0}
    assert any(e["event"] == "run_finished" for e in rec.events)
    assert rec.metrics.get("steps") == 1


def test_submit_defaults_work_dir_under_store(tmp_path):
    """With no work_dir, output files land in a predictable per-run dir."""
    rt = Runtime(RunStore(tmp_path / "store"))
    register_application(_energy_app())
    handle = rt.submit_application("rt_app")     # no work_dir given
    rec = handle.wait(timeout=10)
    assert rec.status == RunStatus.SUCCESS.value
    assert rec.work_dir and str(tmp_path / "store") in rec.work_dir
    assert "work" in rec.work_dir


def test_submit_unknown_application_fails(tmp_path):
    rt = Runtime(RunStore(tmp_path / "store"))
    h = rt.submit_application("does_not_exist")
    rec = h.wait(timeout=5)
    assert rec.status == RunStatus.FAILED.value
    assert "unknown application" in rec.error


def test_run_store_persists_and_lists(tmp_path):
    base = tmp_path / "store"
    rt = Runtime(RunStore(base))
    register_application(_energy_app())
    h = rt.submit_application("rt_app", work_dir=tmp_path / "wd")
    h.wait(timeout=10)
    # a fresh store reads the persisted record from disk
    assert RunStore(base).get(h.id).status == RunStatus.SUCCESS.value
    assert any(r.id == h.id for r in rt.list_runs())


def test_cancel_between_steps(tmp_path):
    rt = Runtime(RunStore(tmp_path / "store"))
    pipe = Pipeline("rt_slow_pipe")
    pipe.add("rt_slow_step", label="s1")
    pipe.add("rt_slow_step", label="s2")
    register_application(Application.from_pipeline(pipe, name="rt_slow_app"))

    h = rt.submit_application("rt_slow_app", work_dir=tmp_path / "wd")
    h.cancel()  # set before step 1 (0.2 s) finishes → checked in on_step after step 1
    rec = h.wait(timeout=10)
    assert rec.status == RunStatus.CANCELLED.value


def test_execute_run_runs_and_updates_record(tmp_path):
    """The compute-node executor runs a submitted record and writes its result."""
    from delfin.tools._runtime import RunRecord, RunStatus, RunStore, execute_run

    register_application(_energy_app())
    store_base = tmp_path / "store"
    store = RunStore(store_base)
    store.save(RunRecord(
        id="x1", kind="application", name="rt_app", inputs={},
        created_at="now", work_dir=str(tmp_path / "wd"), metrics={"cores": 1},
    ))

    execute_run("x1", str(store_base))

    done = store.get("x1")
    assert done.status == RunStatus.SUCCESS.value
    assert done.outputs == {"energy_Eh": -2.0}


def test_submit_slurm_without_sbatch_fails_cleanly(tmp_path, monkeypatch):
    """backend='slurm' on a host without SLURM marks the run failed, no crash."""
    from delfin.tools._runtime import RunStatus, RunStore, Runtime

    monkeypatch.setattr("shutil.which", lambda name=None: None)   # no sbatch
    rt = Runtime(RunStore(tmp_path / "store"))
    register_application(_energy_app())

    handle = rt.submit_application("rt_app", inputs={}, backend="slurm")
    rec = handle.record()
    assert rec.status == RunStatus.FAILED.value
    assert "sbatch not found" in rec.error
    assert rec.metrics.get("backend") == "slurm"


def test_facade_async_run(tmp_path, monkeypatch):
    import delfin.tools._runtime as rtmod

    monkeypatch.setattr(rtmod, "_RUNTIME", Runtime(RunStore(tmp_path / "store")))
    from delfin.tools import platform

    register_application(_energy_app())
    run_id = platform.submit_application("rt_app", work_dir=tmp_path / "wd")
    rec = platform.wait_run(run_id, timeout=10)
    assert rec.status == RunStatus.SUCCESS.value
    assert platform.run_status(run_id) == RunStatus.SUCCESS.value
    assert any(r.id == run_id for r in platform.list_runs())
