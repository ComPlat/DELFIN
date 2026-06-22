"""Guard tests against re-introduction of the dual-implementation drift
risks identified in the post-refactor cleanup audit.

Each test pins one consolidation that the cleanup performed. The pattern
is always the same:

  - A previous refactor left two implementations of the same public name
    in two different modules with diverging signatures.
  - The cleanup picked one canonical implementation and made the other
    site delegate to it (or renamed one to remove the collision).
  - These tests ensure the consolidation stays consolidated.

If a new dual-implementation appears, add it here and consolidate it
the same way: pick a canonical home, make every other site delegate.
"""

from __future__ import annotations

import importlib
import inspect
import subprocess
from pathlib import Path

import pytest


# ===========================================================================
# 1. run_esd_phase: pipeline-phase vs. inner worker (renamed)
# ===========================================================================
#
# Before cleanup: `delfin.esd_module.run_esd_phase` and
# `delfin.workflows.pipeline.run_esd_phase` both existed with incompatible
# signatures. The pipeline phase was `(ctx) -> bool`; the inner worker was
# `(config, charge, solvent, metals, main_basisset, metal_basisset) ->
# WorkflowRunResult`. The inner worker is now `execute_esd_jobs`.


def test_esd_module_no_longer_exposes_run_esd_phase():
    """The inner worker was renamed to execute_esd_jobs to remove the
    name collision with the pipeline phase. If `run_esd_phase` reappears
    on esd_module, two functions with the same public name and different
    signatures coexist again — exactly the drift risk we removed."""
    mod = importlib.import_module("delfin.esd_module")
    assert not hasattr(mod, "run_esd_phase"), (
        "delfin.esd_module.run_esd_phase reappeared. The inner worker is "
        "now `execute_esd_jobs`. Either consolidate the new function into "
        "execute_esd_jobs or rename it to a non-colliding name."
    )


def test_esd_module_execute_esd_jobs_is_the_canonical_worker():
    from delfin.esd_module import execute_esd_jobs
    sig = inspect.signature(execute_esd_jobs)
    params = list(sig.parameters.keys())
    expected = ["config", "charge", "solvent", "metals",
                "main_basisset", "metal_basisset"]
    for i, want in enumerate(expected):
        assert params[i] == want, (
            f"execute_esd_jobs param #{i}: contract says '{want}', got "
            f"'{params[i]}'. Full sig: {sig}"
        )


def test_pipeline_run_esd_phase_keeps_ctx_signature():
    """The pipeline-phase entry point signature stays (ctx) -> bool."""
    from delfin.workflows.pipeline import run_esd_phase
    sig = inspect.signature(run_esd_phase)
    params = list(sig.parameters.keys())
    assert params == ["ctx"], (
        f"workflows.pipeline.run_esd_phase signature changed: {sig}"
    )


# ===========================================================================
# 2. read_control_file: full pipeline parser vs. CO2 standalone parser
# ===========================================================================


def test_co2_coordinator_does_not_export_read_control_file():
    """The CO2 Coordinator's lightweight CONTROL parser was renamed to
    `_minimal_read_control_file` so it no longer collides with
    `delfin.config.read_control_file`. If `read_control_file` reappears
    on the CO2 module, the two parsers risk drifting on new CONTROL keys."""
    mod = importlib.import_module("delfin.co2.CO2_Coordinator6")
    assert not hasattr(mod, "read_control_file"), (
        "delfin.co2.CO2_Coordinator6.read_control_file reappeared — this "
        "name collides with the canonical parser in delfin.config and "
        "risks silent drift. Use `_minimal_read_control_file` for the "
        "standalone-friendly variant."
    )


def test_co2_coordinator_has_minimal_parser_helper():
    """The renamed standalone parser must still be available."""
    from delfin.co2.CO2_Coordinator6 import _minimal_read_control_file
    assert callable(_minimal_read_control_file)


def test_canonical_read_control_file_is_in_config():
    from delfin.config import read_control_file
    assert callable(read_control_file)


# ===========================================================================
# 3. collect_gibbs_energies: pipeline wrapper vs. reporting wrapper, both
#     delegate to delfin.energies.collect_gibbs_energies_from_dir
# ===========================================================================


def test_collect_gibbs_energies_helper_lives_in_energies():
    """The single source of truth for the redox-state file map and the
    Gibbs-energy collection logic must live in delfin.energies. Both the
    pipeline wrapper and the reporting wrapper delegate to it."""
    from delfin.energies import (
        collect_gibbs_energies_from_dir,
        GIBBS_STATE_FILES,
    )
    assert callable(collect_gibbs_energies_from_dir)
    # Schema check: 7 redox states are pinned
    assert set(GIBBS_STATE_FILES.keys()) == {"0", "+1", "+2", "+3",
                                             "-1", "-2", "-3"}


def test_pipeline_collect_gibbs_energies_delegates_to_helper():
    """The pipeline-side wrapper must be a delegator: same result as the
    helper for the same working directory and esd-enabled flag."""
    src = Path(__file__).resolve().parents[1] / "delfin" / "workflows" / "pipeline.py"
    text = src.read_text()
    # Contract: the wrapper imports and calls the helper
    assert "collect_gibbs_energies_from_dir" in text, (
        "pipeline.collect_gibbs_energies must call the helper in "
        "delfin.energies — re-implementing the file map locally re-creates "
        "the drift risk we removed."
    )


def test_reporting_collect_gibbs_energies_delegates_to_helper():
    src = (
        Path(__file__).resolve().parents[1]
        / "delfin" / "reporting" / "delfin_collector.py"
    )
    text = src.read_text()
    assert "collect_gibbs_energies_from_dir" in text, (
        "reporting.delfin_collector.collect_gibbs_energies must delegate "
        "to delfin.energies — re-implementing the state-file table here "
        "would let the JSON exporter look at different files than the "
        "live pipeline."
    )


# ===========================================================================
# 4. Dead modules stay deleted
# ===========================================================================
#
# These were 1900+ lines of unmaintained / unused code, removed during
# the cleanup. If a future commit reintroduces them, this fails.

_DELETED_DEAD_MODULES = [
    "delfin/safe.py",  # 1838 LOC, 4 _safe variants of OCCUPIER summary
]


@pytest.mark.parametrize("relpath", _DELETED_DEAD_MODULES)
def test_dead_module_stays_deleted(relpath):
    p = Path(__file__).resolve().parents[1] / relpath
    assert not p.exists(), (
        f"{relpath} reappeared in the tree. It was deleted in the post-"
        f"refactor cleanup because it was unmaintained dead code. If you "
        f"need to bring it back, consolidate the live functionality into "
        f"its canonical home (delfin.reporting.occupier_reports) instead."
    )


def test_verify_global_manager_moved_to_scripts():
    """Standalone Dev-script for verifying global-manager singleton
    behaviour was moved from delfin/ (where it counted as a package
    member but no one imported it) to scripts/ (where it lives as a
    Dev-tool). Make sure it stays out of the package."""
    pkg_path = Path(__file__).resolve().parents[1] / "delfin" / "verify_global_manager.py"
    assert not pkg_path.exists(), (
        "delfin/verify_global_manager.py is back in the package — it was "
        "moved to scripts/ because it is a standalone diagnostic, not "
        "package code."
    )
    scripts_path = Path(__file__).resolve().parents[1] / "scripts" / "verify_global_manager.py"
    assert scripts_path.exists(), (
        "scripts/verify_global_manager.py is missing. If the script is no "
        "longer wanted, delete it and remove this assertion."
    )


# ===========================================================================
# 5. stability_constant uses canonical workflows.pipeline imports
# ===========================================================================


def test_stability_constant_uses_canonical_pipeline_imports():
    """stability_constant.py used to import internal pipeline helpers
    via the backward-compat shim `delfin.pipeline`. The cleanup migrated
    these to the canonical `delfin.workflows.pipeline` so the codebase
    is internally consistent (only user code should still rely on the
    shim)."""
    src = Path(__file__).resolve().parents[1] / "delfin" / "stability_constant.py"
    text = src.read_text()
    assert "from delfin.pipeline import" not in text, (
        "stability_constant.py imports via the old shim path. Replace "
        "`from delfin.pipeline import X` with "
        "`from delfin.workflows.pipeline import X`."
    )
