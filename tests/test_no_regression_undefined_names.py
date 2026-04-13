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
