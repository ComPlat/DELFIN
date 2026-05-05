from delfin.global_manager import GlobalJobManager
from delfin.workflows.scheduling.manager import (
    _budgeted_memory_mb,
    _slurm_node_memory_mb,
    _SLURM_MEM_HEADROOM,
)


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


def test_slurm_memory_cap_clamps_pool_budget(monkeypatch):
    """When SLURM allocates less mem than PAL × MaxCore, pool must be capped to prevent OOM.

    Scenario from JEW-R465 bug: PAL=60 × MaxCore=6000 MB = 360 GB requested,
    but SLURM allocation is 360 000 MB (~352 GB). Pool budget must respect SLURM,
    otherwise the pool admits two parallel ORCA jobs and OOM-kills hit.
    """
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "360000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)

    budget, slurm = _budgeted_memory_mb(pal=60, maxcore_mb=6000)

    assert slurm == 360000
    assert budget == int(360000 * _SLURM_MEM_HEADROOM)
    assert budget < 60 * 6000


def test_slurm_memory_unset_keeps_request(monkeypatch):
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.delenv("SLURM_CPUS_ON_NODE", raising=False)
    monkeypatch.delenv("SLURM_NTASKS", raising=False)

    budget, slurm = _budgeted_memory_mb(pal=12, maxcore_mb=4000)

    assert slurm is None
    assert budget == 12 * 4000


def test_slurm_memory_per_cpu_fallback(monkeypatch):
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "5000")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "40")

    mem = _slurm_node_memory_mb()

    assert mem == 5000 * 40


def test_slurm_memory_request_below_cap_keeps_request(monkeypatch):
    """Pool with low memory request must not be inflated to the SLURM cap."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "360000")

    budget, slurm = _budgeted_memory_mb(pal=4, maxcore_mb=2000)

    assert slurm == 360000
    assert budget == 8000  # request < cap, take request
