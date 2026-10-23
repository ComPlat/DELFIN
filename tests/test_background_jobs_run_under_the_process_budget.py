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


def test_the_record_is_written_before_the_flag_is_set(registry, monkeypatch,
                                                      tmp_path):
    """The durable copy is the commit point, and the order says so.

    A job killed for a budget breach has to be able to say why, and only
    the registry survives a restart. So the in-memory flag -- which is
    what `bash_status` and every caller observes -- must never be set
    before the record is on disk.

    It was: the field was assigned, stderr was flushed, and the record
    written last. A reattach inside that window read a job that had been
    killed for no stated reason, which is how CI saw this as a flake:
    `status.get("budget_breach")` was None on a job reattached 2.2 s
    after the kill (2026-10-08).

    Asserted from inside the write, because the window is the defect:
    checking afterwards cannot tell which happened first.
    """
    seen: list[object] = []
    real = bash_jobs._update_job_record

    def _watch_order(ws, jid, **fields):
        if "budget_breach" in fields:
            # Whatever the job's observable flag is AT THIS MOMENT.
            seen.append(_the_jobs_flag(jid))
        return real(ws, jid, **fields)

    def _the_jobs_flag(jid):
        job = registry._jobs.get(jid) if hasattr(registry, "_jobs") else None
        return getattr(job, "budget_breach", "no such job")

    monkeypatch.setattr(bash_jobs, "_update_job_record", _watch_order)
    monkeypatch.setattr(
        bash_jobs.process_budget, "limits_for_profile",
        lambda name: bash_jobs.process_budget.BudgetLimits(
            max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9))

    job = registry.start(_script_fanout(12), cwd=str(tmp_path),
                         timeout_s=60, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to be killed")
    _wait_for(lambda: job.budget_breach is not None, 10,
              "the breach to be recorded")

    assert seen, "the breach was never written to the registry"
    assert seen[0] in (None, "no such job"), (
        "the flag was already set when the record was being written, so a "
        f"reattach in that window reads no reason: {seen[0]!r}")


def test_a_registry_that_cannot_be_written_still_explains_the_kill(
        registry, monkeypatch, tmp_path):
    """An explanation that reached nothing is worse than one that reached
    only this process, so a failed write does not swallow the flag."""
    real = bash_jobs._update_job_record
    refused: list[int] = []

    def _refuse_the_breach_write(ws, jid, **fields):
        # Only the breach write, and only the first one. The watchdog's
        # completion write carries the same field and is not guarded, so a
        # stub that refused every call would end the watcher thread and
        # report that instead of what this test is about.
        if "budget_breach" in fields and not refused:
            refused.append(1)
            raise OSError("read-only registry")
        return real(ws, jid, **fields)

    monkeypatch.setattr(bash_jobs, "_update_job_record",
                        _refuse_the_breach_write)
    monkeypatch.setattr(
        bash_jobs.process_budget, "limits_for_profile",
        lambda name: bash_jobs.process_budget.BudgetLimits(
            max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9))

    job = registry.start(_script_fanout(12), cwd=str(tmp_path),
                         timeout_s=60, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to be killed")
    _wait_for(lambda: job.budget_breach is not None, 10,
              "the breach to be recorded in memory anyway")
    assert job.budget_breach["limit"] == "max_processes"
    assert refused, "the write this test is about never happened"


def test_a_registry_that_cannot_be_written_still_ends_the_job(registry,
                                                              monkeypatch,
                                                              tmp_path):
    """A job that has ended must be able to say so.

    The watchdog's completion write was the last statement in its
    `finally`, unguarded. Anything it raised left the watcher thread dead
    with an unhandled exception -- and the record then kept
    `exit_code=None` with no `finished_at`, so every reattach and every
    restart read a finished job as still RUNNING. That is the one state
    the registry exists to rule out.

    Asserted through the thread rather than by calling the writer: the
    defect was where the exception went, and a direct call cannot see
    that.
    """
    import threading

    thread_errors: list = []
    monkeypatch.setattr(threading, "excepthook",
                        lambda args: thread_errors.append(args))

    real = bash_jobs._update_job_record

    def _refuse_the_completion(ws, jid, **fields):
        if "exit_code" in fields:
            raise OSError("read-only registry")
        return real(ws, jid, **fields)

    monkeypatch.setattr(bash_jobs, "_update_job_record",
                        _refuse_the_completion)

    job = registry.start("echo done", cwd=str(tmp_path),
                         timeout_s=30, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to end")
    _wait_for(lambda: job.record_error, 10, "the failure to be recorded")

    assert not thread_errors, (
        "the watcher thread died with an unhandled exception: "
        f"{thread_errors}")
    status = job.status_dict()
    assert status["record_error"], status
    assert "still running" in status["note"], status
    # And where a reader actually looks for the job's output.
    assert "registry could not be written" in job.stderr_path.read_text()


def test_a_failed_submitted_job_scan_does_not_cost_the_exit_code(
        registry, monkeypatch, tmp_path):
    """Two facts, and the exit code is the one that matters.

    The scan for queued cluster jobs ran before the completion write and
    was unguarded, so its failure took the exit code with it.
    """
    def _refuse(ws, rec):
        raise OSError("cannot scan")

    monkeypatch.setattr(bash_jobs, "_watch_submitted_jobs", _refuse)
    job = registry.start("echo done", cwd=str(tmp_path),
                         timeout_s=30, workspace=str(tmp_path))
    _wait_for(lambda: job.poll() is not None, 25, "the job to end")

    # The REGISTRY RECORD, read straight from the file -- not
    # job.finished_at, which is set before the write (waiting on one copy
    # and asserting on another is the mistake the breach ordering was
    # fixed for), and not a reattached object, which adds a second thing
    # that can be absent: this timed out on CI while passing here,
    # because what it was waiting for was `_reattach` and not the record.
    def _record():
        jobs = (bash_jobs._load_registry_file(tmp_path) or {}).get("jobs", {})
        rec = jobs.get(job.job_id) or {}
        return rec if rec.get("exit_code") is not None else None

    rec = _wait_for(_record, 20, "the exit code to reach the registry")
    assert rec.get("exit_code") == 0, rec
    assert rec.get("finished_at"), rec
    assert rec.get("watched_slurm_jobs") == [], (
        "the scan failed, so there is no watched list -- but the exit "
        f"code is written anyway: {rec}")
    assert "submitted-job scan failed" in job.stderr_path.read_text()
