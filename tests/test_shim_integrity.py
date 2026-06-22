"""Guard against shim drift: every backward-compat shim must re-export
the same objects as the canonical module.

If this test fails, a shim has diverged from its canonical source —
exactly the bug class that caused the GlobalJobManager dual-singleton
crash (Tilman, 2026-04).

The test is intentionally exhaustive: it checks identity (``is``), not
equality, so that accidental copy-paste of classes is caught immediately.
"""

import importlib
import pytest


# (shim_module, canonical_module)
SHIM_PAIRS = [
    ("delfin.global_manager", "delfin.workflows.scheduling.manager"),
    ("delfin.global_scheduler", "delfin.workflows.engine.scheduler"),
    ("delfin.job_priority", "delfin.workflows.scheduling.priority"),
    ("delfin.dynamic_pool", "delfin.workflows.scheduling.pool"),
    ("delfin.parallel_classic_manually", "delfin.workflows.engine.classic"),
    ("delfin.parallel_occupier", "delfin.workflows.engine.occupier"),
    ("delfin.pipeline", "delfin.workflows.pipeline"),
]


# ---------------------------------------------------------------------------
# Names that external code actually imports through the old (shim) path.
# If a name is added here, the corresponding shim MUST re-export it.
# ---------------------------------------------------------------------------
REQUIRED_SHIM_EXPORTS = {
    "delfin.global_manager": [
        "GlobalJobManager",
        "get_global_manager",
        "bootstrap_global_manager_from_env",
        "_TrackedProcess",
        "_normalize_parallel_token",
        "_safe_int",
    ],
    "delfin.global_scheduler": [
        "GlobalOrcaScheduler",
    ],
    "delfin.job_priority": [
        "adjust_job_priorities",
        "is_exclusive_bottleneck",
    ],
    "delfin.dynamic_pool": [
        "DynamicCorePool",
        "JobPriority",
        "PoolJob",
        "get_current_job_id",
        "get_current_job_cores",
        "_set_current_job_id",
        "_set_current_job_cores",
    ],
    "delfin.parallel_classic_manually": [
        "WorkflowJob",
        "WorkflowRunResult",
        "_WorkflowManager",
        "execute_classic_workflows",
        "execute_manually_workflows",
        "normalize_parallel_token",
        "_update_pal_block",
    ],
    "delfin.parallel_occupier": [
        "OccupierExecutionContext",
        "run_occupier_orca_jobs",
        "build_flat_occupier_fob_jobs",
    ],
    "delfin.pipeline": [
        "FileBundle",
        "PipelineContext",
        "run_occuper_phase",
        "run_classic_phase",
        "run_manual_phase",
        "run_esd_phase",
        "run_hyperpol_xtb_phase",
        "run_tadf_xtb_phase",
        "collect_gibbs_energies",
        "compute_summary",
        "interpret_method_alias",
        "normalize_input_file",
        "_run_guppy_for_smiles",
        "_resolve_smiles_converter",
        "_skip_xtb_goat_after_guppy",
    ],
}


# ---------------------------------------------------------------------------
# 1. Every required name must exist in the shim
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shim_path", list(REQUIRED_SHIM_EXPORTS.keys()))
def test_shim_exports_required_names(shim_path):
    """Shim must re-export every name that external code imports."""
    shim_mod = importlib.import_module(shim_path)
    missing = [
        name
        for name in REQUIRED_SHIM_EXPORTS[shim_path]
        if not hasattr(shim_mod, name)
    ]
    assert not missing, (
        f"{shim_path} is missing required exports: {missing}\n"
        f"Add explicit imports to the shim file."
    )


# ---------------------------------------------------------------------------
# 2. Every required name must be the SAME OBJECT as in the canonical module
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shim_path,canonical_path", SHIM_PAIRS)
def test_shim_identity(shim_path, canonical_path):
    """Shim names must be identical objects (not copies) to canonical names."""
    shim_mod = importlib.import_module(shim_path)
    canonical_mod = importlib.import_module(canonical_path)

    names_to_check = REQUIRED_SHIM_EXPORTS.get(shim_path, [])
    mismatches = []
    for name in names_to_check:
        shim_obj = getattr(shim_mod, name, None)
        canonical_obj = getattr(canonical_mod, name, None)
        if shim_obj is not canonical_obj:
            mismatches.append(name)

    assert not mismatches, (
        f"Shim {shim_path} returns DIFFERENT objects than {canonical_path} "
        f"for: {mismatches}\n"
        f"This means the shim has been replaced with a copy — "
        f"restore it to a thin re-export."
    )


# ---------------------------------------------------------------------------
# 3. GlobalJobManager singleton must be unique across all import paths
# ---------------------------------------------------------------------------
def test_global_manager_singleton_identity():
    """The singleton that caused the original crash: importing via different
    paths must return the exact same instance."""
    from delfin.global_manager import get_global_manager as get_via_shim
    from delfin.workflows.scheduling.manager import get_global_manager as get_via_canonical

    mgr_shim = get_via_shim()
    mgr_canonical = get_via_canonical()

    assert mgr_shim is mgr_canonical, (
        "GlobalJobManager singleton differs between import paths!\n"
        f"  via shim:      {id(mgr_shim):#x} ({type(mgr_shim).__module__})\n"
        f"  via canonical: {id(mgr_canonical):#x} ({type(mgr_canonical).__module__})\n"
        "This is the exact bug that caused the cluster crash."
    )


# ---------------------------------------------------------------------------
# 4. All CLI entry points must be importable
# ---------------------------------------------------------------------------
CLI_ENTRY_POINTS = [
    ("delfin.main", "main"),
    ("delfin.cli_step", "main"),
    ("delfin.cli_pipeline", "main"),
    ("delfin.cli_voila", "main"),
    ("delfin.build_up_complex", "main"),
    ("delfin.guppy_sampling", "main"),
    ("delfin.cli_esd_report", "main"),
    ("delfin.cli_ir_report", "main"),
    ("delfin.cli_nmr_report", "main"),
    ("delfin.cli_delfin_collect", "main"),
]


@pytest.mark.parametrize("module_path,func_name", CLI_ENTRY_POINTS)
def test_cli_entry_point_importable(module_path, func_name):
    """Every console_script entry point must resolve to a callable."""
    mod = importlib.import_module(module_path)
    func = getattr(mod, func_name, None)
    assert func is not None, f"{module_path} has no attribute '{func_name}'"
    assert callable(func), f"{module_path}.{func_name} is not callable"


# ---------------------------------------------------------------------------
# 5. Cross-module singleton sharing (the real-world scenario)
# ---------------------------------------------------------------------------
def test_workflow_manager_shares_pool_with_global_manager():
    """WorkflowManager (via any import path) must use the same pool
    as the GlobalJobManager — otherwise resource allocation breaks."""
    from delfin.global_manager import get_global_manager
    from delfin.parallel_classic_manually import _WorkflowManager

    mgr = get_global_manager()
    mgr.initialize({
        "PAL": "2",
        "maxcore": "500",
        "parallel_workflows": "auto",
    })

    try:
        wm = _WorkflowManager(
            {"PAL": "2", "maxcore": "500", "parallel_workflows": "auto"},
            label="shim_test",
        )
        assert wm.pool is mgr.pool, (
            "WorkflowManager created via shim uses a different pool "
            "than GlobalJobManager — resource allocation will be broken."
        )
        wm.shutdown()
    finally:
        mgr.shutdown()
