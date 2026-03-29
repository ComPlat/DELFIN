import threading

from delfin.parallel_classic_manually import JobPriority, WorkflowJob, _WorkflowManager


def _make_manager(total_cores=30):
    manager = _WorkflowManager.__new__(_WorkflowManager)
    manager.config = {}
    manager.total_cores = total_cores
    manager.max_jobs = 4
    manager.maxcore_mb = 3800
    manager.label = "test"
    manager._parallel_enabled = True
    manager._lock = threading.RLock()
    manager._inflight = set()
    manager._jobs = {}
    manager._completed = set()
    manager._failed = {}
    manager._skipped = {}
    manager._check_exclusive_bottleneck_boost = lambda job: None
    return manager


def _make_job(job_id, *, deps=None, cores_min=4, cores_optimal=8, cores_max=30):
    return WorkflowJob(
        job_id=job_id,
        work=lambda cores: None,
        description=f"{job_id} optimization",
        dependencies=set(deps or ()),
        cores_min=cores_min,
        cores_optimal=cores_optimal,
        cores_max=cores_max,
        priority=JobPriority.NORMAL,
        memory_mb=cores_optimal * 3800,
    )


def test_single_ready_job_uses_full_available_cores_even_with_future_work_remaining():
    manager = _make_manager(total_cores=30)
    ready = _make_job("initial", deps=())
    blocked = _make_job("ox_step_1", deps={"initial"})
    manager._jobs = {ready.job_id: ready, blocked.job_id: blocked}

    allocations, caps = manager._plan_core_allocations([ready], {"running_jobs": 0, "queued_jobs": 0, "allocated_cores": 0})

    assert allocations == {"initial": 30}
    assert caps == {"initial": False}


def test_single_ready_job_consumes_currently_available_pool_capacity():
    manager = _make_manager(total_cores=30)
    ready = _make_job("initial", deps=(), cores_max=30)
    manager._jobs = {ready.job_id: ready}
    manager._inflight = {"other_running_job"}

    allocations, caps = manager._plan_core_allocations(
        [ready],
        {"running_jobs": 1, "queued_jobs": 0, "allocated_cores": 12},
    )

    assert allocations == {"initial": 18}
    assert caps == {"initial": False}


def test_aggressive_default_does_not_wait_speculatively_for_bottlenecks():
    manager = _make_manager(total_cores=30)
    ready = _make_job("initial", deps=())
    manager._jobs = {ready.job_id: ready}
    manager._inflight = {"other_running_job"}

    should_wait = manager._should_wait_for_exclusive_bottleneck(
        [ready],
        {"running_jobs": 1, "queued_jobs": 0, "allocated_cores": 20},
    )

    assert should_wait is False


def test_stage_profile_file_provides_duration_and_core_fallbacks():
    manager = _make_manager(total_cores=40)
    job = _make_job("ox_step_1", deps=(), cores_max=40)

    duration_hint = manager._get_duration_hint(job)
    core_hint = manager._get_stage_core_hint(job)

    assert duration_hint is not None
    assert duration_hint > 600
    assert core_hint is None


def test_observed_stage_core_hints_can_be_enabled_explicitly():
    manager = _make_manager(total_cores=40)
    manager.config = {"scheduler_stage_core_learning": "observed"}
    job = _make_job("ox_step_1", deps=(), cores_max=40)

    core_hint = manager._get_stage_core_hint(job)

    assert core_hint is not None
    assert 6.0 <= core_hint <= 10.0
