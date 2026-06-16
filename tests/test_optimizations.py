"""Tests for cross-cutting optimizations: step cache, recovery, run metrics."""

from __future__ import annotations

import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepStatus


# --- step cache -----------------------------------------------------------

_calls = {"count": 0}


class _CountingStep(StepAdapter):
    name = "cache_count_step"
    description = "counts executions and writes a geometry"
    produces_geometry = True

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        _calls["count"] += 1
        g = work_dir / "out.xyz"
        g.write_text("1\n\nH 0.0 0.0 0.0\n")
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, time.monotonic(),
            geometry=g, data={"val": kwargs.get("x", 1)},
        )


register(_CountingStep())


def test_cache_hit_skips_execution(tmp_path):
    from delfin.tools._cache import StepCache, run_step_cached

    cache = StepCache(tmp_path / "cache")
    _calls["count"] = 0

    r1 = run_step_cached("cache_count_step", work_dir=tmp_path / "w1", cache=cache, x=5)
    assert r1.ok and _calls["count"] == 1
    assert r1.data.get("_cache") == "miss"

    r2 = run_step_cached("cache_count_step", work_dir=tmp_path / "w2", cache=cache, x=5)
    assert _calls["count"] == 1                      # not re-executed
    assert r2.data.get("_cache") == "hit"
    assert r2.data.get("val") == 5
    assert r2.geometry is not None and r2.geometry.is_file()   # materialised into w2

    r3 = run_step_cached("cache_count_step", work_dir=tmp_path / "w3", cache=cache, x=7)
    assert _calls["count"] == 2 and r3.data.get("_cache") == "miss"  # different param


def test_cache_key_is_stable_and_param_sensitive():
    from delfin.tools._cache import cache_key

    assert cache_key("s", {"a": 1}) == cache_key("s", {"a": 1})
    assert cache_key("s", {"a": 1}) != cache_key("s", {"a": 2})
    # internal keys are ignored
    assert cache_key("s", {"a": 1, "_prev": 9}) == cache_key("s", {"a": 1})


# --- recovery policies ----------------------------------------------------


def test_classify_error():
    from delfin.tools._recovery import classify_error

    assert classify_error("SCF NOT CONVERGED") is ErrorKind.CONVERGENCE
    assert classify_error("the process timed out") is ErrorKind.TIMEOUT
    assert classify_error("") is ErrorKind.NONE
    assert classify_error("weird boom") is ErrorKind.TOOL_FAILED


def test_suggest_recovery_convergence():
    from delfin.tools._recovery import suggest_recovery

    a = suggest_recovery(ErrorKind.CONVERGENCE, {})
    assert a.action == "modify"
    assert a.kwargs_overrides["scf_maxiter"] >= 300
    assert "SOSCF" in a.kwargs_overrides["scf_extra"]
    assert suggest_recovery(ErrorKind.BINARY_NOT_FOUND).action == "give_up"


class _FlakyStep(StepAdapter):
    name = "recover_flaky"
    description = "fails with a convergence error until scf_maxiter is raised"
    produces_geometry = False

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        if int(kwargs.get("scf_maxiter", 0) or 0) >= 300:
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir,
                                     time.monotonic(), data={"ok": True})
        return self._make_result(self.name, StepStatus.FAILED, work_dir,
                                 time.monotonic(), error="SCF NOT CONVERGED")


register(_FlakyStep())


def test_run_step_with_recovery_recovers(tmp_path):
    from delfin.tools._recovery import run_step_with_recovery

    r = run_step_with_recovery("recover_flaky", work_dir=tmp_path / "rec", max_attempts=3)
    assert r.ok
    log = r.data["_recovery"]
    assert len(log) >= 2                       # failed, then recovered
    assert log[0]["error_kind"] == "convergence"
    assert log[0]["action"] == "modify"


# --- run metrics ----------------------------------------------------------


def test_aggregate_runs(tmp_path):
    from delfin.tools._metrics import aggregate_runs
    from delfin.tools._runtime import RunRecord, RunStatus, RunStore

    store = RunStore(tmp_path / "store")
    store.save(RunRecord(id="a", kind="application", name="app1",
                         status=RunStatus.SUCCESS.value, metrics={"elapsed_s": 2.0}))
    store.save(RunRecord(id="b", kind="application", name="app1",
                         status=RunStatus.FAILED.value, metrics={"elapsed_s": 1.0}))

    m = aggregate_runs(store)
    assert m["total_runs"] == 2
    assert m["by_status"]["success"] == 1 and m["by_status"]["failed"] == 1
    assert m["total_elapsed_s"] == 3.0
    app = m["by_application"]["app1"]
    assert app["count"] == 2 and app["success_rate"] == 0.5
