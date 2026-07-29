"""Guard against jobs that can never fit into the pool.

Regression cover for the scheduler deadlock behind job 4646143: a job whose
memory request exceeded the pool budget was rejected on every scheduling pass,
re-queued forever, and the run idled at 0% CPU until the walltime kill.
"""

import logging
import queue
import threading
from itertools import count

import pytest

from delfin.workflows.scheduling.pool import DynamicCorePool, JobPriority, PoolJob


def _make_pool(total_cores=40, total_memory_mb=136000, max_jobs=4):
    """Build a pool without starting its monitor and scheduler threads."""
    pool = DynamicCorePool.__new__(DynamicCorePool)
    pool.total_cores = total_cores
    pool.total_memory_mb = total_memory_mb
    pool.max_concurrent_jobs = max_jobs
    pool.config = {}
    pool._allocated_cores = 0
    pool._allocated_memory = 0
    pool._lock = threading.RLock()
    pool._condition = threading.Condition(pool._lock)
    pool._running_jobs = {}
    pool._job_queue = queue.PriorityQueue()
    pool._completed_jobs = []
    pool._failed_jobs = []
    pool._job_counter = count()
    pool._job_metrics = None
    pool._shutdown = False
    pool._resource_event = threading.Event()
    pool._capacity_clamp_warnings = {}
    pool._idle_stall_warnings = {}
    return pool


def _make_job(job_id="initial_fob_1", *, cores_min=4, cores_optimal=13,
              cores_max=40, memory_mb=52000):
    return PoolJob(
        job_id=job_id,
        cores_min=cores_min,
        cores_optimal=cores_optimal,
        cores_max=cores_max,
        memory_mb=memory_mb,
        priority=JobPriority.HIGH,
        execute_func=lambda *a, **k: None,
        args=(),
        kwargs={},
    )


def test_job_exceeding_pool_budget_still_gets_scheduled_into_empty_pool():
    """The exact shape of job 4646143: PAL=40 x maxcore=4000 vs 136000 MB budget.

    Before the guard, _calculate_optimal_allocation returned None forever.
    """
    pool = _make_pool(total_cores=40, total_memory_mb=136000)
    job = _make_job(cores_min=40, cores_optimal=40, cores_max=40, memory_mb=160000)

    cores = pool._calculate_optimal_allocation(job)

    assert cores is not None, "job must be schedulable into an idle pool"
    assert cores == 40
    assert job.memory_mb == 136000


def test_submit_clamps_memory_above_pool_budget_and_logs_error(caplog):
    pool = _make_pool(total_memory_mb=136000)
    job = _make_job(memory_mb=160000)

    with caplog.at_level(logging.ERROR):
        pool.submit_job(job)

    assert job.memory_mb == 136000
    assert "more than the pool can ever provide" in caplog.text
    assert "160000 MB > pool budget 136000 MB" in caplog.text


def test_submit_clamps_cores_min_above_pool_size_and_keeps_core_invariant(caplog):
    pool = _make_pool(total_cores=40)
    job = _make_job(cores_min=64, cores_optimal=64, cores_max=64, memory_mb=1000)

    with caplog.at_level(logging.ERROR):
        pool.submit_job(job)

    assert job.cores_min == 40
    assert job.cores_min <= job.cores_optimal <= job.cores_max <= pool.total_cores
    assert "cores_min 64 > pool size 40 cores" in caplog.text


def test_clamp_is_logged_once_per_job(caplog):
    pool = _make_pool(total_memory_mb=136000)
    job = _make_job(memory_mb=160000)

    with caplog.at_level(logging.ERROR):
        pool._clamp_to_pool_capacity(job, "submit")
        job.memory_mb = 160000  # simulate a later mutation of the same job
        pool._clamp_to_pool_capacity(job, "allocation")

    assert job.memory_mb == 136000
    assert caplog.text.count("more than the pool can ever provide") == 1


def test_fitting_job_is_left_untouched_and_silent(caplog):
    pool = _make_pool(total_cores=40, total_memory_mb=136000)
    job = _make_job(cores_min=4, cores_optimal=13, cores_max=40, memory_mb=52000)

    with caplog.at_level(logging.ERROR):
        cores = pool._calculate_optimal_allocation(job)

    assert cores == 13
    assert (job.cores_min, job.cores_optimal, job.cores_max) == (4, 13, 40)
    assert job.memory_mb == 52000
    assert caplog.text == ""


def test_pal_wide_job_matching_the_budget_exactly_is_untouched(caplog):
    """The pre-cap production shape (Jerome, RSS_paper): PAL x maxcore == budget.

    40 x 6000 MB against a 240000 MB pool fits exactly. Those runs must keep
    behaving as before: no clamp, no ERROR, full PAL allocation.
    """
    pool = _make_pool(total_cores=40, total_memory_mb=240000)
    job = _make_job(cores_min=40, cores_optimal=40, cores_max=40, memory_mb=240000)

    with caplog.at_level(logging.ERROR):
        cores = pool._calculate_optimal_allocation(job)

    assert cores == 40
    assert job.memory_mb == 240000
    assert caplog.text == ""


def test_wedged_pool_with_waiting_job_is_reported(caplog):
    """The line that was missing from the run log for 57 hours."""
    pool = _make_pool(total_cores=40, total_memory_mb=136000)
    job = _make_job(cores_min=40, memory_mb=160000)
    job.queue_time = 0.0  # queued "long ago"
    pool._job_queue.put((1, 0.0, 0, job))

    with caplog.at_level(logging.ERROR):
        pool._prevent_starvation()

    assert "still does not fit" in caplog.text
    assert "cannot make progress" in caplog.text


def test_wedge_is_reported_once_not_every_monitor_tick(caplog):
    pool = _make_pool()
    job = _make_job(cores_min=40, memory_mb=160000)
    job.queue_time = 0.0
    pool._job_queue.put((1, 0.0, 0, job))

    with caplog.at_level(logging.ERROR):
        pool._prevent_starvation()
        pool._prevent_starvation()

    assert caplog.text.count("still does not fit") == 1


def test_leaked_core_accounting_on_an_idle_pool_is_reported(caplog):
    """The other way a pool wedges: allocation counters never released."""
    pool = _make_pool(total_cores=40, total_memory_mb=136000)
    pool._allocated_cores = 40  # nothing running, yet everything booked out
    job = _make_job(cores_min=4, memory_mb=52000)
    job.queue_time = 0.0
    pool._job_queue.put((1, 0.0, 0, job))

    with caplog.at_level(logging.ERROR):
        pool._prevent_starvation()

    assert "still does not fit" in caplog.text


def test_no_wedge_report_while_other_jobs_are_running(caplog):
    """A busy pool is normal queueing, not a deadlock."""
    pool = _make_pool()
    running = _make_job("running_job")
    pool._running_jobs = {running.job_id: running}

    waiting = _make_job("waiting_job")
    waiting.queue_time = 0.0
    pool._job_queue.put((1, 0.0, 0, waiting))

    with caplog.at_level(logging.ERROR):
        pool._prevent_starvation()

    assert "still does not fit" not in caplog.text


def test_no_wedge_report_for_a_long_waiting_job_that_fits(caplog):
    """Idle pool, job about to start on the next pass — must stay quiet."""
    pool = _make_pool(total_cores=40, total_memory_mb=136000)
    job = _make_job(cores_min=4, cores_optimal=13, memory_mb=52000)
    job.queue_time = 0.0
    pool._job_queue.put((1, 0.0, 0, job))

    with caplog.at_level(logging.ERROR):
        pool._prevent_starvation()

    assert caplog.text == ""


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
