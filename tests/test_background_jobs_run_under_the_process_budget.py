"""A background job runs under the process budget (phase 3).

bash_background jobs outlive the tool call that started them, so the
foreground guard does not fit: it would die with the call. The budget
must be enforced by something that lives exactly as long as the job —
the job's own watchdog thread — and a breach must be visible in
bash_status and bash_output, including after a reattach from the
persistent registry.
"""
import json
import os
import signal
import time

import pytest

from delfin.agent import bash_jobs


@pytest.fixture
def registry(tmp_path):
    reg = bash_jobs._Registry()
    yield reg
    for job in reg.list_jobs(include_finished=False):
        reg.kill(job.job_id)


def _script_fanout(n: int, sleep_s: int = 20) -> str:
    return f"for i in $(seq 1 {n}); do sleep {sleep_s} & done; sleep {sleep_s}"


def _wait_for(condition, timeout_s=15.0, what="condition"):
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        value = condition()
        if value:
            return value
        time.sleep(0.2)
    raise AssertionError(f"timed out waiting for {what}")


def test_registry_start_passes_a_budget_profile(registry, monkeypatch,
                                                tmp_path):
    """The public start() arms a guard with the background profile."""
    guards = []
    real_guard = bash_jobs.process_budget.BudgetGuard

    class SpyGuard:
        def __init__(self, root_pid, **kw):
            self.kw = kw
            guards.append(self)
            self._real = real_guard(root_pid, **kw)

        def poll_once(self):
            return self._real.poll_once()

    monkeypatch.setattr(bash_jobs.process_budget, "BudgetGuard", SpyGuard)
    job = registry.start("echo hello", cwd=str(tmp_path),
                         timeout_s=30, workspace=str(tmp_path))
    try:
        assert guards, "no guard was armed for the background job"
        assert guards[0].kw.get("profile") == "background"
    finally:
        registry.kill(job.job_id)


def test_breach_ends_the_job_and_names_the_limit(registry, monkeypatch,
                                                 tmp_path):
    """A fan-out past the background profile's process limit is killed;
    bash_status reports the breach; the stderr log carries the message."""
    # Shrink the profile for the test: 4 processes max.
    monkeypatch.setattr(
        bash_jobs.process_budget, "limits_for_profile",
        lambda name: bash_jobs.process_budget.BudgetLimits(
            max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9))
    job = registry.start(_script_fanout(12), cwd=str(tmp_path),
                         timeout_s=60, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to be killed")
    rc = job.poll()
    assert rc is not None
    # Killed by our guard (SIGTERM/SIGKILL group) -> negative or nonzero.
    assert rc != 0
    # The record lands in the watchdog's finally, which may lag the kill.
    _wait_for(lambda: job.budget_breach is not None, 10,
              "the breach to be recorded")
    status = job.status_dict()
    assert status.get("budget_breach"), status
    breach = status["budget_breach"]
    assert breach["limit"] == "max_processes"
    assert breach["profile"] == "background"
    assert breach["value"] > breach["ceiling"]
    # The message reached the job's stderr log, so bash_output sees it.
    err = job.stderr_path.read_text()
    assert "process budget" in err.lower()
    assert "max_processes" in err


def test_breach_survives_a_reattach(registry, monkeypatch, tmp_path):
    """The breach is persisted in the registry, so a reattached job
    (restart orphaned it) still reports it in status."""
    monkeypatch.setattr(
        bash_jobs.process_budget, "limits_for_profile",
        lambda name: bash_jobs.process_budget.BudgetLimits(
            max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9))
    job = registry.start(_script_fanout(12), cwd=str(tmp_path),
                         timeout_s=60, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to be killed")
    # Simulate the restart: forget the in-memory job, reattach from disk.
    _wait_for(lambda: job.budget_breach is not None, 10,
              "the breach to be recorded")
    reattached = bash_jobs._reattach(job.job_id, workspace=str(tmp_path))
    assert reattached is not None
    _wait_for(lambda: reattached.poll() is not None, 10, "reattach poll")
    status = reattached.status_dict()
    assert status.get("budget_breach"), status
    assert status["budget_breach"]["limit"] == "max_processes"


def test_an_in_budget_job_is_not_killed(registry, monkeypatch, tmp_path):
    """A quiet job under the (shrunken) limit runs to its natural end."""
    monkeypatch.setattr(
        bash_jobs.process_budget, "limits_for_profile",
        lambda name: bash_jobs.process_budget.BudgetLimits(
            max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9))
    job = registry.start("echo hello", cwd=str(tmp_path),
                         timeout_s=30, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 15, "the job to finish")
    assert job.poll() == 0
    status = job.status_dict()
    assert "budget_breach" not in status
