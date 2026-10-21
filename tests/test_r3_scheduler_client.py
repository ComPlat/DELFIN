"""R3 phase 2: job_monitor's scheduler queries are throttled process-wide.

Finding 1: a scheduler query would have called ``squeue`` every 5 s per
idle terminal; cluster operators had objected to per-widget ``squeue``
load before. The fix is ``delfin.scheduler_client`` — the one place that
asks the scheduler — a process-wide client that lets a distinct query
kind reach the scheduler at most once per >=25 s and serves the shared
cache while the gate is closed.

These tests use fake runners and a fake clock: never the real scheduler,
never a real sleep, so nothing depends on the machine or the wall clock.
"""
from __future__ import annotations

import time

import pytest

from delfin import scheduler_client as sc
import delfin.agent.job_monitor as jm


def _fake_runner(answer: str):
    """A runner that records every call and returns ``answer``."""
    calls: list[list[str]] = []

    def run(cmd):
        calls.append(list(cmd))
        return answer

    run.calls = calls  # type: ignore[attr-defined]
    return run


class _FakeClock:
    """Monotonic-ish clock the test advances by hand."""

    def __init__(self, start: float = 1000.0):
        self.t = start

    def __call__(self) -> float:
        return self.t

    def advance(self, seconds: float) -> None:
        self.t += seconds


@pytest.fixture(autouse=True)
def _reset_global_client():
    """Keep the process-wide client clean between tests."""
    sc.reset()
    yield
    sc.reset()


@pytest.fixture()
def clock() -> _FakeClock:
    return _FakeClock()


def _client(clock: _FakeClock, runner) -> sc.SchedulerClient:
    return sc.SchedulerClient(min_interval_s=sc.MIN_INTERVAL_S,
                              clock=clock, runner=runner)


# ---------------------------------------------------------------------------
# The client's own contract (new module).
# ---------------------------------------------------------------------------

def test_min_interval_is_at_least_25_seconds():
    assert sc.MIN_INTERVAL_S >= 25.0


def test_same_query_served_from_cache_within_interval(clock):
    runner = _fake_runner("JOB1 RUNNING")
    client = _client(clock, runner)
    assert client.query("squeue", ["squeue", "-j", "1"]) == "JOB1 RUNNING"
    # second ask within the interval must NOT reach the runner again
    assert client.query("squeue", ["squeue", "-j", "1"]) == "JOB1 RUNNING"
    assert len(runner.calls) == 1
    assert client.queries_run == 1
    assert client.served_from_cache == 1


def test_after_interval_query_runs_again(clock):
    runner = _fake_runner("JOB1 DONE")
    client = _client(clock, runner)
    client.query("squeue", ["squeue", "-j", "1"])
    clock.advance(sc.MIN_INTERVAL_S + 0.1)
    assert client.query("squeue", ["squeue", "-j", "1"]) == "JOB1 DONE"
    assert len(runner.calls) == 2


def test_different_kinds_each_get_an_ask(clock):
    runner = _fake_runner("ok")
    client = _client(clock, runner)
    client.query("squeue", ["squeue", "-j", "1"])
    # sacct is a different kind: the gate for it is open
    client.query("sacct", ["sacct", "-j", "2"])
    assert len(runner.calls) == 2


def test_other_argv_same_kind_is_throttled_not_stale(clock):
    runner = _fake_runner("ok")
    client = _client(clock, runner)
    client.query("squeue", ["squeue", "-j", "1"])
    # Different argv, same kind, gate still closed: must NOT run the runner
    # and must NOT invent an answer for job 2.
    assert client.query("squeue", ["squeue", "-j", "2"]) is sc.THROTTLED
    assert len(runner.calls) == 1


def test_runner_returning_none_is_cached_like_an_answer(clock):
    runner = _fake_runner(None)  # scheduler said "no such job"
    client = _client(clock, runner)
    assert client.query("squeue", ["squeue", "-j", "9"]) is None
    assert client.query("squeue", ["squeue", "-j", "9"]) is None
    assert len(runner.calls) == 1


def test_kind_of_argv_uses_basename():
    assert sc.kind_of_argv(["/usr/bin/squeue", "-j", "1"]) == "squeue"
    assert sc.kind_of_argv(["sacct", "-j", "1"]) == "sacct"


def test_module_level_query_routes_through_default_client(clock):
    # ``sc.query`` is the module-level helper; it must be throttled by the
    # SAME process-wide default_client, so any component that calls it shares
    # the gate with job_monitor instead of opening a fresh one.
    runner = _fake_runner("RUNNING")
    sc.default_client._runner = runner
    sc.default_client._clock = clock
    assert sc.query("squeue", ["squeue", "-j", "1"]) == "RUNNING"
    assert sc.query("squeue", ["squeue", "-j", "1"]) == "RUNNING"
    assert len(runner.calls) == 1
    assert sc.default_client.queries_run == 1


# ---------------------------------------------------------------------------
# The actual fix: job_monitor's default runner goes through the client.
# RED before, green after. This is the control.
# ---------------------------------------------------------------------------

def test_job_monitor_default_run_serves_cached_below_interval(clock):
    """Two consecutive ``_default_run`` calls of the same squeue argv must
    reach the runner exactly once — the throttle bites right at the seam
    job_monitor uses, so a watch loop cannot hammer the scheduler."""
    runner = _fake_runner("JOB1 RUNNING")
    sc.default_client._runner = runner  # inject fake into the process-wide gate
    sc.default_client._clock = clock

    first = jm._default_run(["squeue", "-j", "1"])
    second = jm._default_run(["squeue", "-j", "1"])
    assert first == "JOB1 RUNNING"
    assert second == "JOB1 RUNNING"
    assert len(runner.calls) == 1
    assert sc.default_client.queries_run == 1


def test_job_monitor_default_run_asks_again_after_interval(clock):
    runner = _fake_runner("JOB1 DONE")
    sc.default_client._runner = runner
    sc.default_client._clock = clock
    jm._default_run(["squeue", "-j", "1"])
    clock.advance(sc.MIN_INTERVAL_S + 0.1)
    assert jm._default_run(["squeue", "-j", "1"]) == "JOB1 DONE"
    assert len(runner.calls) == 2
