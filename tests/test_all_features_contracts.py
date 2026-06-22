"""Master contract test for all user-visible DELFIN features.

The ONE place that says "DELFIN exposes these features, they are still
importable, their signatures are still these, their dataclasses still
have these fields, their CLI entry points still work". Every refactor
has to either keep these intact or explicitly update the contract.

This file is intentionally shallow-and-wide: many small tests, each
pointing at a distinct feature. Depth tests (parsing a full
CONTROL.txt, exercising FoB job graphs) live in:

  - test_co2_coordination_workflow.py  (GUPPY / OCCUPIER / CO2)
  - test_orca_workflow_contracts.py     (classic / hyperpol / TADF /
                                         NMR / raw ORCA)

Here we sweep across the entire feature surface so no subsystem is
silently lost. Each @pytest.mark.parametrize block is a feature
category — adding a new subsystem means adding a row, not rewriting
test logic.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect

import pytest


# ===========================================================================
# Helper: feature-registry-like table of contracts
# ===========================================================================
#
# Each table row is one "promise DELFIN makes about this feature".
# Categorised so failures point at the right owner.


# ---------------------------------------------------------------------------
# 1. Workflows registered in delfin.workflows.registry
# ---------------------------------------------------------------------------

_EXPECTED_WORKFLOWS = {
    "tadf_xtb",
    "hyperpol",
    "occupier",
    "esd",
    "imag",
    "co2_coordinator",
}


def test_workflow_registry_lists_all_expected_entries():
    """The workflows registry is the canonical catalogue — a refactor
    that forgets to re-register a workflow silently removes it."""
    from delfin.workflows import registry

    # Registry API exposes list_workflows or similar
    if hasattr(registry, "list_workflows"):
        registered = set(registry.list_workflows().keys())
    elif hasattr(registry, "WORKFLOWS"):
        registered = set(registry.WORKFLOWS.keys())
    else:
        pytest.skip("no recognisable registry accessor; update this test")

    missing = _EXPECTED_WORKFLOWS - registered
    assert not missing, (
        f"Workflows dropped from registry: {sorted(missing)}. "
        f"Users invoking these via CONTROL.txt toggles will silently no-op."
    )


# ---------------------------------------------------------------------------
# 2. Tools adapter registry — 20+ tool bindings
# ---------------------------------------------------------------------------

_EXPECTED_TOOL_ADAPTER_MODULES = [
    "delfin.tools.adapters.afp",
    "delfin.tools.adapters.cclib",
    "delfin.tools.adapters.censo",
    "delfin.tools.adapters.crest",
    "delfin.tools.adapters.genarris",
    "delfin.tools.adapters.guppy",
    "delfin.tools.adapters.imag",
    "delfin.tools.adapters.mlp",
    "delfin.tools.adapters.morfeus",
    "delfin.tools.adapters.multiwfn",
    "delfin.tools.adapters.orca",
    "delfin.tools.adapters.packmol",
    "delfin.tools.adapters.reporting",
    "delfin.tools.adapters.smiles",
    "delfin.tools.adapters.uv_vis",
    "delfin.tools.adapters.xtb",
    "delfin.tools.adapters.xtb_solvator",
]


@pytest.mark.parametrize("module_name", _EXPECTED_TOOL_ADAPTER_MODULES)
def test_tool_adapter_module_importable(module_name):
    """Every tool adapter must stay importable. If one crashes on
    import (missing `logger`, bad module path in a `from X import Y`),
    every PipelineTemplate that uses this step silently fails."""
    importlib.import_module(module_name)


# ---------------------------------------------------------------------------
# 3. Tools runtime: run_step / StepResult / StepStatus
# ---------------------------------------------------------------------------


def test_tools_run_step_is_importable():
    """The run_step public API from delfin.tools._runner is the
    uniform entrypoint for invoking any tool. If this disappears,
    PipelineTemplate-based user code has nowhere to go."""
    from delfin.tools._runner import run_step
    assert callable(run_step)


def test_tools_step_result_has_essential_fields():
    from delfin.tools._types import StepResult
    fields = {f.name for f in dataclasses.fields(StepResult)}
    expected = {"step_name", "status"}
    missing = expected - fields
    assert not missing, (
        f"StepResult lost {sorted(missing)} — every tool adapter's return "
        f"shape depends on these."
    )


def test_tools_step_status_is_importable_and_has_states():
    from delfin.tools._types import StepStatus
    # Must at minimum distinguish success vs failure
    values = {v.name for v in StepStatus}
    expected_min = {"SUCCESS", "FAILED"}
    missing = expected_min - values
    assert not missing, (
        f"StepStatus lost {sorted(missing)} — adapters can no longer "
        f"report success/failure to the pipeline."
    )


# ---------------------------------------------------------------------------
# 4. CLI commands — every pyproject.toml-declared entry point works
# ---------------------------------------------------------------------------

_CLI_ENTRY_POINTS = [
    ("delfin", "delfin.main", "main"),
    ("delfin_ESD", "delfin.cli_esd_report", "main"),
    ("delfin_IR", "delfin.cli_ir_report", "main"),
    ("delfin_NMR", "delfin.cli_nmr_report", "main"),
    ("delfin-json", "delfin.cli_delfin_collect", "main"),
    ("delfin-build", "delfin.build_up_complex", "main"),
    ("delfin-guppy", "delfin.guppy_sampling", "main"),
    ("delfin-voila", "delfin.cli_voila", "main"),
    ("delfin-step", "delfin.cli_step", "main"),
    ("delfin-pipeline", "delfin.cli_pipeline", "main"),
    ("delfin-docs-index", "delfin.doc_server.indexer", "main"),
]


@pytest.mark.parametrize(
    "cli_name,module,func", _CLI_ENTRY_POINTS,
    ids=[c[0] for c in _CLI_ENTRY_POINTS],
)
def test_cli_command_target_resolves(cli_name, module, func):
    m = importlib.import_module(module)
    assert hasattr(m, func), f"CLI `{cli_name}` → {module}:{func} broken"
    assert callable(getattr(m, func))


# ---------------------------------------------------------------------------
# 5. SMILES → XYZ converters (GUPPY, NORMAL, hapto)
# ---------------------------------------------------------------------------


def test_smiles_to_xyz_public_api():
    """The quick SMILES → XYZ conversion is used by the NORMAL
    converter path and as GUPPY fallback."""
    from delfin.smiles_converter import smiles_to_xyz
    assert callable(smiles_to_xyz)

    sig = inspect.signature(smiles_to_xyz)
    params = list(sig.parameters.keys())
    assert params[0] == "smiles", (
        f"smiles_to_xyz first param must be 'smiles', got {params[0]!r}"
    )


def test_guppy_sampling_main_entry():
    """GUPPY is triggered by CONTROL.txt smiles_converter=GUPPY. Its
    main() is the script entry point `delfin-guppy` resolves to."""
    from delfin.guppy_sampling import main
    assert callable(main)


# ---------------------------------------------------------------------------
# 6. Stability Constant workflow (thermodynamics)
# ---------------------------------------------------------------------------


def test_stability_constant_dataclasses():
    """The stability-constant workflow publishes data classes that
    downstream reporting consumes. Pin the field contracts."""
    from delfin.stability_constant import (
        LigandInfo, StabilityResult, ReactionStabilityAnalysis,
    )

    lig_fields = {f.name for f in dataclasses.fields(LigandInfo)}
    assert "charge" in lig_fields, (
        "LigandInfo.charge is gone — stability-constant workflow can no "
        "longer reconstruct the neutral ligand charge balance."
    )

    res_fields = {f.name for f in dataclasses.fields(StabilityResult)}
    assert res_fields, "StabilityResult has no fields — definitely broken"


def test_stability_analyze_complex_signature():
    """analyze_complex(smiles, solvent, n_explicit_solvent) is the
    main analysis helper. Renaming any arg breaks callers."""
    from delfin.stability_constant import analyze_complex

    sig = inspect.signature(analyze_complex)
    params = list(sig.parameters.keys())
    assert params[:3] == ["smiles", "solvent", "n_explicit_solvent"], (
        f"analyze_complex signature changed: {params[:3]}"
    )


# ---------------------------------------------------------------------------
# 7. ANMR (analysis tools)
# ---------------------------------------------------------------------------


def test_anmr_availability_helpers_exist():
    """The ANMR pipeline checks `anmr_available()` before running.
    Losing it silently breaks ensemble NMR workflows."""
    from delfin.analysis_tools import anmr_available, get_anmr_path
    assert callable(anmr_available)
    assert callable(get_anmr_path)


# ---------------------------------------------------------------------------
# 8. Build-up complex (3D metal-complex assembly)
# ---------------------------------------------------------------------------


def test_build_up_complex_public_api():
    """`delfin-build` CLI entry + the run_build_up / run_goat helpers
    are what the CLI invokes."""
    from delfin.build_up_complex import main, run_build_up, run_goat_optimization
    assert callable(main)
    assert callable(run_build_up)
    assert callable(run_goat_optimization)


# ---------------------------------------------------------------------------
# 9. ORCA helpers: smart_recalc / orca_recovery / qm_runtime
# ---------------------------------------------------------------------------


def test_smart_recalc_api():
    """Smart recalc is what lets users re-submit failed DELFIN jobs and
    skip already-finished ORCA steps. Its public API is should_skip()."""
    from delfin.smart_recalc import should_skip, smart_mode_enabled
    assert callable(should_skip)
    assert callable(smart_mode_enabled)


def test_orca_recovery_module_importable():
    """Intelligent error-recovery wrapper around run_orca. If this
    module breaks, enable_auto_recovery=yes silently does nothing."""
    mod = importlib.import_module("delfin.orca_recovery")
    # Should export the main recovery coordinator
    public = [n for n in dir(mod) if not n.startswith("_")]
    assert public, "orca_recovery exposes no public API — broken"


def test_qm_runtime_api():
    """qm_runtime knows how to find ORCA / xTB binaries on disk. The
    pipeline depends on find_tool_executable() to locate them before
    launching jobs."""
    from delfin.qm_runtime import find_tool_executable
    assert callable(find_tool_executable)
    sig = inspect.signature(find_tool_executable)
    params = list(sig.parameters.keys())
    assert params, f"find_tool_executable has no parameters: {params}"


# ---------------------------------------------------------------------------
# 10. Reporting sub-package: the write-out side of every workflow
# ---------------------------------------------------------------------------


def test_reporting_occupier_summary_api():
    from delfin.reporting import generate_summary_report_OCCUPIER
    assert callable(generate_summary_report_OCCUPIER)

    sig = inspect.signature(generate_summary_report_OCCUPIER)
    params = list(sig.parameters.keys())
    expected_first_args = {"duration", "fspe_values", "is_even", "charge"}
    missing = expected_first_args - set(params)
    assert not missing, (
        f"generate_summary_report_OCCUPIER lost params {sorted(missing)}"
    )


def test_reporting_subpackage_importable():
    import delfin.reporting
    # These submodules are the ones reporting CLIs hit
    for sub in ("nmr_report", "occupier_reports"):
        importlib.import_module(f"delfin.reporting.{sub}")


# ---------------------------------------------------------------------------
# 11. Runtime & dashboard bootstrapping
# ---------------------------------------------------------------------------


def test_runtime_setup_public_api():
    """apply_runtime_environment / detect_local_runtime_limits are
    called by submit_delfin.sh indirectly via the dashboard. If the
    symbols move, SLURM jobs fail at env-bootstrap."""
    from delfin.runtime_setup import (
        apply_runtime_environment,
        detect_local_runtime_limits,
        resolve_backend_choice,
    )
    assert callable(apply_runtime_environment)
    assert callable(detect_local_runtime_limits)
    assert callable(resolve_backend_choice)


def test_user_settings_loader_api():
    """load_settings is called at every dashboard startup. Contract:
    must take an optional settings_path and return a dict."""
    from delfin.user_settings import load_settings, load_remote_archive_enabled
    assert callable(load_settings)
    assert callable(load_remote_archive_enabled)


# ---------------------------------------------------------------------------
# 12. Scheduling primitives — used by every pipeline
# ---------------------------------------------------------------------------


def test_scheduling_primitives_still_on_shims():
    """Old-path imports `from delfin.global_manager import ...` etc.
    must keep working — thousands of user scripts hardcode them."""
    from delfin.global_manager import get_global_manager, GlobalJobManager
    from delfin.dynamic_pool import DynamicCorePool, JobPriority, PoolJob
    from delfin.job_priority import (
        count_downstream_jobs,
        is_exclusive_bottleneck,
        adjust_job_priorities,
    )

    for fn in (get_global_manager, count_downstream_jobs,
               is_exclusive_bottleneck, adjust_job_priorities):
        assert callable(fn)


# ---------------------------------------------------------------------------
# 13. Config parser: both soft (OCCUPIER_parser) and strict paths
# ---------------------------------------------------------------------------


def test_config_parsers_importable():
    """OCCUPIER_parser / read_control_file / parse_control_text are
    the three CONTROL.txt entry points. All three are used at
    different places in the pipeline."""
    from delfin.config import (
        OCCUPIER_parser, read_control_file, parse_control_text,
    )
    assert callable(OCCUPIER_parser)
    assert callable(read_control_file)
    assert callable(parse_control_text)


# ---------------------------------------------------------------------------
# 14. Pipeline phases — master list of entry points
# ---------------------------------------------------------------------------

_PIPELINE_PHASES = [
    "run_occuper_phase",
    "run_classic_phase",
    "run_manual_phase",
    "run_esd_phase",
    "run_hyperpol_xtb_phase",
    "run_tadf_xtb_phase",
]


@pytest.mark.parametrize("phase_name", _PIPELINE_PHASES)
def test_pipeline_phase_entry_point(phase_name):
    """Master list of pipeline phase entry points — each is called
    from delfin.cli based on a CONTROL.txt toggle. Losing any breaks
    a subset of workflows silently."""
    from delfin.workflows import pipeline
    phase = getattr(pipeline, phase_name, None)
    assert callable(phase), (
        f"Pipeline phase {phase_name} is missing from "
        f"delfin.workflows.pipeline. CONTROL.txt toggles routing to it "
        f"will silently skip their phase."
    )
    # Signature: all phase functions take (ctx: PipelineContext)
    sig = inspect.signature(phase)
    params = list(sig.parameters.keys())
    assert params, f"{phase_name} takes no args — broken signature"


# ---------------------------------------------------------------------------
# 15. Submit-templates infrastructure (SLURM + local)
# ---------------------------------------------------------------------------


def test_submit_templates_directory_exists():
    """The SLURM submit script lives under delfin/submit_templates/
    and is copied by the dashboard on submit. Losing the directory
    breaks every SLURM submission."""
    from delfin.runtime_setup import get_packaged_submit_templates_dir
    d = get_packaged_submit_templates_dir()
    assert d.exists(), f"{d} is missing — SLURM submit will fail"
    # The main script must be there
    main_script = d / "submit_delfin.sh"
    assert main_script.exists(), (
        f"{main_script} is missing — all SLURM jobs will fail to submit"
    )


# ---------------------------------------------------------------------------
# 16. Logging infrastructure (every module uses get_logger)
# ---------------------------------------------------------------------------


def test_common_logging_configurable():
    """configure_logging and get_logger are the two log-infrastructure
    functions every module imports. Pin both."""
    from delfin.common.logging import get_logger, configure_logging
    lg = get_logger("contract_test")
    assert lg is not None
    assert hasattr(lg, "info")
    assert callable(configure_logging)
