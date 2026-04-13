"""Regression test that catches NameErrors, missing imports, and scope bugs
introduced by refactoring.

This is the safety net Tilman and Jan's bug reports made necessary. The
`cf78545` / `3ed6566` / `39b1919` refactorings all introduced latent
NameErrors that existing tests did not catch because they didn't
exercise the exact code path. This test shifts left: instead of waiting
for users to hit the bug in production, pyflakes inspects every module
at build-time and fails the test if it finds undefined names.

What counts as a failure here:
- `undefined name 'X'` where X is a module-level reference (real bug)
- A module that cannot even be imported (ImportError / NameError at
  module load)

What is explicitly ignored:
- `undefined name '<typing-name>'` inside function bodies — these are
  variable annotations which Python does not evaluate at runtime (see
  PEP 526), so pyflakes flags a false positive.
- Forward-reference strings like Optional["SomeClass"] which pyflakes
  sometimes flags but Python never evaluates.
- Unused imports / unused variables — noise, not a bug.
- Shim files that do `from X import *` (backward-compat layers).
"""

from __future__ import annotations

import importlib
import pkgutil
import subprocess
import sys
from pathlib import Path

import pytest


_ROOT = Path(__file__).resolve().parents[1]
_DELFIN = _ROOT / "delfin"


# Modules we don't try to import because they have heavy side effects
# (start subprocesses, open dashboards, etc.) or they require runtime
# extensions that aren't installed in the test environment. Importing
# them is a nice-to-have, not required to prove the refactor is sound.
_SKIP_IMPORT = {
    # Heavy optional deps — missing on CI / dev machines
    "delfin.dashboard.tab_agent",
    "delfin.dashboard",  # pulls in ipywidgets, voila, big GUI stack
    # CLI entry points that do sys.argv parsing on import
    "delfin.cli",
    "delfin.cli_nmr_report",
    "delfin.cli_esd_report",
    "delfin.cli_ir_report",
    "delfin.cli_voila",
    "delfin.cli_imag",
    "delfin.cli_recalc",
    "delfin.cli_delfin_collect",
    # Agent subsystem needs model providers that may not be installed
    "delfin.agent",
    # Build tools that expect external binaries
    "delfin.build_up_complex",
    "delfin.build_up_complex2",
    "delfin.guppy_sampling",
}


# Typing names that pyflakes falsely flags as "undefined" inside function
# bodies because they're part of PEP 526 variable annotations (never
# evaluated at runtime). Also forward-reference strings in annotations
# that resolve correctly at runtime.
_PYFLAKES_FALSE_POSITIVES = {
    "'Tuple'",
    "'List'",
    "'Dict'",
    "'Set'",
    "'Any'",
    "'Callable'",
    "'Optional'",
    "'Iterator'",
    "'Iterable'",
    "'Union'",
    "'Sequence'",
    "'Mapping'",
    "'GlobalOrcaScheduler'",   # forward-ref string
    "'_WorkflowManager'",       # forward-ref string
    "'PipelineContext'",        # forward-ref string
    "'JobPriority'",            # forward-ref string
    "'WorkflowJob'",            # forward-ref string
}


def _iter_delfin_modules():
    """Yield all importable dotted module names under delfin/."""
    for py in _DELFIN.rglob("*.py"):
        # Skip build artifacts, tests, vendored code
        if any(part in py.parts for part in (
            "__pycache__", ".build", "_deps", "ComPlat_DELFIN.egg-info"
        )):
            continue
        if py.name == "__init__.py":
            dotted = ".".join(py.parent.relative_to(_ROOT).parts)
        else:
            dotted = ".".join(py.relative_to(_ROOT).with_suffix("").parts)
        yield dotted, py


def test_every_delfin_module_imports_cleanly():
    """Every non-skipped delfin module must import without raising.

    This is the fastest way to catch refactor-time bugs like:
      - missing `import copy` when `copy.deepcopy(...)` is called
      - missing `logger = get_logger(...)` when `logger.info(...)` is called
      - missing import of a helper after a module split
      - module-level reference to a name not yet defined

    Each of the five bugs from Tilman's and Jan's reports would fail
    this test (except the inline-branch cores bug which needed a call
    to trigger).
    """
    failures = []
    for dotted, path in _iter_delfin_modules():
        if dotted in _SKIP_IMPORT:
            continue
        # Skip anything under a skipped package
        if any(dotted.startswith(pkg + ".") for pkg in _SKIP_IMPORT):
            continue
        try:
            importlib.import_module(dotted)
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{dotted}: {type(exc).__name__}: {exc}")

    if failures:
        failed_list = "\n  ".join(failures)
        pytest.fail(
            f"{len(failures)} module(s) failed to import:\n  {failed_list}\n\n"
            "Each failure is a refactor that broke an import chain. "
            "Fix the import or adjust _SKIP_IMPORT if the module is "
            "intentionally not importable in the test env."
        )


def _run_pyflakes_on_all_delfin() -> list[str]:
    """Return the list of pyflakes warnings that are real bugs."""
    result = subprocess.run(
        [sys.executable, "-m", "pyflakes", str(_DELFIN)],
        capture_output=True,
        text=True,
    )
    real_bugs: list[str] = []
    for line in result.stdout.splitlines():
        # Only care about "undefined name 'X'" lines
        if "undefined name" not in line:
            continue
        # Skip false positives (typing names, forward refs)
        if any(fp in line for fp in _PYFLAKES_FALSE_POSITIVES):
            continue
        # Skip backward-compat shim stars
        if "from" in line and "import *" in line:
            continue
        # Skip the legacy safe.py module (dead code, no imports from it)
        if "/delfin/safe.py:" in line:
            continue
        real_bugs.append(line)
    return real_bugs


def test_no_undefined_names_in_delfin_tree():
    """pyflakes finds no module-level undefined names in any delfin module.

    This is how we would have caught:
      - classic.py: `_parse_step_set` used without import
      - occupier.py: `copy` used without import
      - occupier.py: `folder_path` referenced outside its closure
      - tools/adapters/orca.py: `logger` used without definition
      - runtime_setup.py: `logger` typo (should be _logger)
      - smiles_converter.py: `np`/`deque`/`Geometry` used without import
      - reporting/occupier_reports.py: `occ_method` undefined in safe variant
      - tab_calculations_browser.py: `_calc_update_png_frame` nonexistent call
      - tab_calculations_browser.py: `delete_current` undefined local
    """
    try:
        import pyflakes  # noqa: F401
    except ImportError:
        pytest.skip("pyflakes not installed in test env")

    bugs = _run_pyflakes_on_all_delfin()
    if bugs:
        formatted = "\n  ".join(bugs)
        pytest.fail(
            f"pyflakes found {len(bugs)} real undefined-name bug(s):\n  "
            f"{formatted}\n\n"
            "Each of these will raise NameError at runtime when the code "
            "path is exercised. Either import the missing name, fix the "
            "scope, or delete the dead reference."
        )


# ---------------------------------------------------------------------------
# Shim export consistency
# ---------------------------------------------------------------------------
#
# The workflows/ and tools/ refactor in cf78545 left behind backward-compat
# shim files at the old flat paths (delfin/global_manager.py,
# delfin/dynamic_pool.py, etc.). Every import in user code and tests still
# goes through these shims. If a shim stops re-exporting a symbol the
# target module defines, user code breaks silently at import.
#
# This test pins the shim contract: for every shim, every public name the
# target module exports must also be accessible on the shim.

_BACKWARD_COMPAT_SHIMS = {
    "delfin.dynamic_pool": "delfin.workflows.scheduling.pool",
    "delfin.global_manager": "delfin.workflows.scheduling.manager",
    "delfin.global_scheduler": "delfin.workflows.engine.scheduler",
    "delfin.job_priority": "delfin.workflows.scheduling.priority",
    "delfin.parallel_classic_manually": "delfin.workflows.engine.classic",
    "delfin.parallel_occupier": "delfin.workflows.engine.occupier",
    "delfin.pipeline": "delfin.workflows.pipeline",
}


@pytest.mark.parametrize("shim_name,target_name", sorted(_BACKWARD_COMPAT_SHIMS.items()))
def test_shim_reexports_all_target_public_names(shim_name, target_name):
    """Backward-compat shims must re-export every public name of their target.

    If the target module adds a new export, the shim's `from X import *`
    must also pick it up. Missing re-exports break any caller that
    imports via the old flat path.
    """
    import importlib

    shim = importlib.import_module(shim_name)
    target = importlib.import_module(target_name)

    target_exports = set(getattr(target, "__all__", None) or [
        name for name in dir(target) if not name.startswith("_")
    ])
    shim_names = set(dir(shim))

    missing = target_exports - shim_names
    assert not missing, (
        f"{shim_name} is missing re-exports from {target_name}: "
        f"{sorted(missing)}. The shim uses `from X import *` — if a name "
        "is in the target's __all__ but not on the shim, the star-import "
        "is filtering it out or the target added it after the shim was "
        "last regenerated."
    )


# ---------------------------------------------------------------------------
# Smoke import of critical public entry points
# ---------------------------------------------------------------------------

_CRITICAL_IMPORTS = [
    # Workflow engine — anything below here breaks every CONTROL.txt run
    ("delfin.workflows.pipeline", "run_occuper_phase"),
    ("delfin.workflows.pipeline", "run_classic_phase"),
    ("delfin.workflows.engine.classic", "_WorkflowManager"),
    ("delfin.workflows.engine.classic", "WorkflowJob"),
    ("delfin.workflows.engine.classic", "_parse_step_set"),
    ("delfin.workflows.engine.occupier", "build_flat_occupier_fob_jobs"),
    ("delfin.workflows.engine.occupier", "OccupierExecutionContext"),
    ("delfin.workflows.engine.occupier", "run_occupier_orca_jobs"),
    ("delfin.workflows.scheduling.manager", "get_global_manager"),
    ("delfin.workflows.scheduling.pool", "DynamicCorePool"),
    ("delfin.workflows.scheduling.pool", "JobPriority"),
    # Tools adapter layer
    ("delfin.tools._registry", "register"),
    ("delfin.tools._types", "StepResult"),
    ("delfin.tools._types", "StepStatus"),
    # Core delfin
    ("delfin.orca", "run_orca"),
    ("delfin.orca", "run_orca_with_intelligent_recovery"),
    ("delfin.user_settings", "load_settings"),
    ("delfin.runtime_setup", "apply_runtime_environment"),
    ("delfin.runtime_setup", "detect_local_runtime_limits"),
    ("delfin.smart_recalc", "should_skip"),
    ("delfin.common.logging", "get_logger"),
    # Reporting
    ("delfin.reporting", "generate_summary_report_OCCUPIER"),
]


@pytest.mark.parametrize("module_name,symbol_name", _CRITICAL_IMPORTS)
def test_critical_public_api_is_importable(module_name, symbol_name):
    """Public names that other modules rely on must stay importable.

    This is the contract-test companion to the shim test: if a refactor
    renames or removes a public function without updating callers, this
    test fails loudly before users see a traceback.
    """
    import importlib

    module = importlib.import_module(module_name)
    assert hasattr(module, symbol_name), (
        f"{module_name}.{symbol_name} is missing. This is a public API "
        f"other modules/shims depend on — a refactor either removed it "
        f"or forgot to export it. Either restore the symbol or update "
        f"_CRITICAL_IMPORTS after consciously dropping the API."
    )
