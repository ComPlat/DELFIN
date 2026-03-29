from delfin.global_manager import GlobalJobManager


def _make_manager():
    manager = GlobalJobManager.__new__(GlobalJobManager)
    manager.total_cores = 1
    manager.maxcore_per_job = 1000
    return manager


def test_pal_jobs_defaults_prioritize_makespan_for_medium_pal():
    manager = _make_manager()

    sanitized = manager._sanitize_resource_config({"PAL": 30, "maxcore": 3800})

    assert sanitized["pal_jobs"] == 3


def test_pal_jobs_defaults_prioritize_makespan_for_large_pal():
    manager = _make_manager()

    sanitized = manager._sanitize_resource_config({"PAL": 64, "maxcore": 3800})

    assert sanitized["pal_jobs"] == 6


def test_parallel_disable_forces_single_slot_even_without_explicit_pal_jobs():
    manager = _make_manager()

    sanitized = manager._sanitize_resource_config(
        {"PAL": 64, "maxcore": 3800, "parallel_workflows": "disable"}
    )

    assert sanitized["pal_jobs"] == 1


def test_resolve_job_resources_caps_to_active_global_limits():
    manager = _make_manager()
    manager.total_cores = 24
    manager.maxcore_per_job = 3500

    cores, maxcore = manager.resolve_job_resources(requested_cores=40, requested_maxcore=6000)

    assert cores == 24
    assert maxcore == 3500


def test_resolve_job_resources_uses_global_defaults_when_request_missing():
    manager = _make_manager()
    manager.total_cores = 12
    manager.maxcore_per_job = 2800

    cores, maxcore = manager.resolve_job_resources()

    assert cores == 12
    assert maxcore == 2800
