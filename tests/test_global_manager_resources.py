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


def test_sanitize_derates_maxcore_when_request_exceeds_slurm_cap(monkeypatch):
    """PAL × maxcore > SLURM × headroom must derate maxcore so a single PAL-wide
    job still fits the pool. Without this, the scheduler's memory gate blocks the
    job forever and the SLURM allocation walltimes out at exit 138 (job 4233067
    bug, 2026-05-08): pool 91800 MB, esd_S0 wants 12 × 9000 = 108000 MB, no
    progress for 47 h, CPU utilised 29 s.
    """
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "108000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    manager = _make_manager()

    sanitized = manager._sanitize_resource_config({"PAL": 12, "maxcore": 9000})

    cap = int(108000 * _SLURM_MEM_HEADROOM)
    assert sanitized["PAL"] == 12
    assert sanitized["maxcore"] == cap // 12
    assert sanitized["PAL"] * sanitized["maxcore"] <= cap


def test_sanitize_keeps_maxcore_when_request_fits_slurm_cap(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "108000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    manager = _make_manager()

    sanitized = manager._sanitize_resource_config({"PAL": 12, "maxcore": 4000})

    assert sanitized["maxcore"] == 4000


def test_sanitize_derate_idempotent_under_slurm_cap(monkeypatch):
    """Re-sanitizing an already-derated config must not derate again."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "108000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    manager = _make_manager()

    first = manager._sanitize_resource_config({"PAL": 12, "maxcore": 9000})
    second = manager._sanitize_resource_config(first)

    assert first["maxcore"] == second["maxcore"]


def test_initialize_propagates_derated_maxcore_to_caller_dict(monkeypatch):
    """Manager must mutate the caller's config dict with the sanitized maxcore.

    Without this, ClassicEngine.maxcore_mb (read from config.get('maxcore')
    in classic.py:199) keeps the original 9000 MB while the pool was capped at
    91800 MB. _resolve_memory(forced=12) then returns 12 * 9000 = 108000 MB,
    pool gate 91800 < 108000 blocks forever, esd_S0 deadlocks until walltime
    (job 4649394 bug, 2026-05-11). Covers all initialize() callers: cli.py
    (initial submit + recalc), cli_imag.py, occupier.py (nested).
    """
    from delfin.global_manager import GlobalJobManager

    monkeypatch.setenv("SLURM_MEM_PER_NODE", "108000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)

    # Reset singleton so initialize creates a fresh pool
    GlobalJobManager._instance = None
    mgr = GlobalJobManager()

    config = {"PAL": 12, "maxcore": 9000, "pal_jobs": 2}
    mgr.initialize(config)

    cap = int(108000 * _SLURM_MEM_HEADROOM)
    expected_derated = cap // 12
    assert config["maxcore"] == expected_derated, (
        f"Caller's config dict not propagated: got maxcore={config['maxcore']}, "
        f"expected derated {expected_derated}. Engine.maxcore_mb would stay "
        f"at the un-derated value and trigger esd_S0 deadlock."
    )
    assert mgr.maxcore_per_job == expected_derated
    assert config["PAL"] == 12
    # cleanup so other tests get a fresh singleton
    mgr.shutdown()
    GlobalJobManager._instance = None


def test_initialize_preserves_caller_dict_when_no_derate_needed(monkeypatch):
    """When PAL × maxcore fits the SLURM cap, the caller's maxcore is unchanged."""
    from delfin.global_manager import GlobalJobManager

    monkeypatch.setenv("SLURM_MEM_PER_NODE", "108000")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)

    GlobalJobManager._instance = None
    mgr = GlobalJobManager()

    config = {"PAL": 12, "maxcore": 4000, "pal_jobs": 2}
    mgr.initialize(config)

    assert config["maxcore"] == 4000  # 12 * 4000 = 48000 < 91800, no derate
    mgr.shutdown()
    GlobalJobManager._instance = None


def test_resolve_headroom_default_is_90_percent(monkeypatch):
    """Default headroom raised from 0.85 → 0.90 on 2026-05-12. Pin it.

    Anyone who needs the older conservative 0.85 should set the env var or
    Dashboard Settings; the in-code default reflects the validated Fritz-
    archive measurement that DFT-OPT real RSS is ~80 % of promised.
    """
    from delfin.workflows.scheduling.manager import (
        _resolve_slurm_mem_headroom,
        _SLURM_MEM_HEADROOM_DEFAULT,
    )
    monkeypatch.delenv("DELFIN_SLURM_MEM_HEADROOM", raising=False)
    assert _SLURM_MEM_HEADROOM_DEFAULT == 0.90
    # When no env override and no settings file with override, default holds
    assert abs(_resolve_slurm_mem_headroom() - 0.90) < 1e-6 or \
           0.50 <= _resolve_slurm_mem_headroom() <= 0.99  # (settings could override)


def test_resolve_headroom_env_override(monkeypatch):
    """ENV var beats settings + default — for power-users + ad-hoc tuning."""
    from delfin.workflows.scheduling.manager import _resolve_slurm_mem_headroom

    monkeypatch.setenv("DELFIN_SLURM_MEM_HEADROOM", "0.85")
    assert abs(_resolve_slurm_mem_headroom() - 0.85) < 1e-6


def test_resolve_headroom_clamps_above_max(monkeypatch):
    """Values >= 1.0 risk OOM — clamp to 0.99."""
    from delfin.workflows.scheduling.manager import _resolve_slurm_mem_headroom

    monkeypatch.setenv("DELFIN_SLURM_MEM_HEADROOM", "1.5")
    assert abs(_resolve_slurm_mem_headroom() - 0.99) < 1e-6


def test_resolve_headroom_clamps_below_min(monkeypatch):
    """Values < 0.50 waste resources — clamp to 0.50."""
    from delfin.workflows.scheduling.manager import _resolve_slurm_mem_headroom

    monkeypatch.setenv("DELFIN_SLURM_MEM_HEADROOM", "0.10")
    assert abs(_resolve_slurm_mem_headroom() - 0.50) < 1e-6


def test_resolve_headroom_invalid_falls_back(monkeypatch):
    """Garbage env value silently falls back to settings/default."""
    from delfin.workflows.scheduling.manager import _resolve_slurm_mem_headroom

    monkeypatch.setenv("DELFIN_SLURM_MEM_HEADROOM", "not-a-number")
    value = _resolve_slurm_mem_headroom()
    # Either the persisted settings value, or the default — must be within range
    assert 0.50 <= value <= 0.99


def test_user_settings_round_trip_headroom(monkeypatch, tmp_path):
    """Settings can persist scheduling.slurm_mem_headroom and load it back."""
    from delfin import user_settings as us

    settings_file = tmp_path / "delfin_settings.json"
    monkeypatch.setattr(us, "get_settings_path", lambda *a, **k: settings_file)

    payload = {"scheduling": {"slurm_mem_headroom": 0.85}}
    us.save_settings(payload)
    loaded = us.load_settings()
    assert loaded["scheduling"]["slurm_mem_headroom"] == 0.85


def test_user_settings_rejects_invalid_headroom(monkeypatch, tmp_path):
    """Out-of-range headroom in settings file is rejected with a clear error."""
    from delfin import user_settings as us

    settings_file = tmp_path / "delfin_settings.json"
    monkeypatch.setattr(us, "get_settings_path", lambda *a, **k: settings_file)

    import pytest as _pt
    with _pt.raises(ValueError, match="headroom"):
        us.save_settings({"scheduling": {"slurm_mem_headroom": 1.5}})
