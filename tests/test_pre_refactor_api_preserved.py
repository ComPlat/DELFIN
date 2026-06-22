"""Pre-refactor API preservation tests.

The workflows/tools refactor in commit cf78545 (22.03.2026) relocated
7 top-level modules into delfin/workflows/ and delfin/tools/, leaving
backward-compat shims at the old flat paths. This test file verifies
that NO user-visible public API was lost or silently changed in that
move — which is the single strongest guarantee we can give users that
"DELFIN still works the way it did before".

How it works:
    For each shim path, we git-show the version at cf78545^ (the commit
    RIGHT BEFORE the refactor), parse out every top-level public
    `def`/`class`, and require that:

      a) the name is still accessible on the shim today
      b) the function's positional-arg names still match position-by-
         position
      c) no keyword-only argument was removed

If a refactor legitimately drops an old API, that has to be explicit:
remove the symbol from _WAIVED_DROPS and document why.

This is the companion to test_no_regression_undefined_names.py and
test_functional_contracts.py. Those check that the current code is
well-formed and matches its pinned contracts. THIS file checks that
we didn't lose anything users were relying on from the old code.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import subprocess
from pathlib import Path

import pytest


# The commit RIGHT BEFORE the big refactor — cf78545 is the refactor,
# cf78545^ is the last pre-refactor commit.
_PRE_REFACTOR_REF = "cf78545^"


# Shim files: every module that cf78545 hollowed out. After the refactor,
# each of these is a 1–8 line stub that does `from <target> import *`.
# Every public name they used to define must still be reachable via the
# shim, so old `from delfin.X import Y` imports in user code keep working.
_SHIM_MODULES = [
    "delfin.parallel_occupier",
    "delfin.parallel_classic_manually",
    "delfin.dynamic_pool",
    "delfin.global_manager",
    "delfin.global_scheduler",
    "delfin.job_priority",
    "delfin.pipeline",
]


# ---------------------------------------------------------------------------
# Intentionally-dropped or -renamed APIs (waiver list)
# ---------------------------------------------------------------------------
#
# This set is the ONE place where we document "yes, we know this was a
# public API before cf78545 and we deliberately removed/renamed it".
# Any entry here should be accompanied by:
#
#   1. A comment explaining the rationale and the replacement (if any)
#   2. A reference to the commit that intentionally dropped the symbol
#   3. A note about migration path for users of the old API
#
# The test below will SKIP the symbols in this set instead of failing.
# If you don't add it here, the test will fail loudly with a message
# telling the next developer why they're breaking user code.
#
# In other words:
#   - New public API → just add it to the target module. No test change.
#     The new name shows up as a superset on the shim, old tests still
#     pass, users get more functionality for free.
#   - Intentional rename/drop of old API → ADD the entry here, document
#     why, the test then skips that symbol. The diff of this file tells
#     code review exactly what user-visible behaviour is being broken.
#   - Accidental rename/drop of old API → test fails, blocks the merge.
#
# This is the "direkt ersichtlich wenn alte Funktionalität kaputt geht"
# guarantee users and reviewers rely on.
_WAIVED_DROPS: set[tuple[str, str]] = set()


def _git_show(ref: str, path: str) -> str | None:
    """Return the contents of ``path`` at ``ref`` or None if missing."""
    result = subprocess.run(
        ["git", "show", f"{ref}:{path}"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    return result.stdout


def _extract_top_level_public_symbols(src: str) -> dict[str, dict]:
    """Parse top-level `def`/`class` names and their arg signatures.

    Returns dict mapping ``name`` → {"kind": "def"|"class",
    "pos_args": [...], "kwonly_args": [...], "has_vararg": bool,
    "has_kwarg": bool}. Only public names (no leading underscore).
    """
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return {}

    out: dict[str, dict] = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            if node.name.startswith("_"):
                continue
            out[node.name] = {"kind": "class"}
        elif isinstance(node, ast.FunctionDef):
            if node.name.startswith("_"):
                continue
            args = node.args
            out[node.name] = {
                "kind": "def",
                "pos_args": [a.arg for a in args.args],
                "kwonly_args": [a.arg for a in args.kwonlyargs],
                "has_vararg": args.vararg is not None,
                "has_kwarg": args.kwarg is not None,
            }
    return out


# ---------------------------------------------------------------------------
# Build the expected vs. actual comparison once (at test collection time)
# ---------------------------------------------------------------------------


def _build_comparison_cases():
    """Yield (shim_name, symbol_name, old_meta) tuples for every symbol
    the pre-refactor module publicly exported."""
    cases: list[tuple[str, str, dict]] = []
    for shim_name in _SHIM_MODULES:
        path = shim_name.replace(".", "/") + ".py"
        old_src = _git_show(_PRE_REFACTOR_REF, path)
        if old_src is None:
            continue
        for name, meta in _extract_top_level_public_symbols(old_src).items():
            cases.append((shim_name, name, meta))
    return cases


_COMPARISON_CASES = _build_comparison_cases()


# ===========================================================================
# 1. Every pre-refactor public symbol still exists on its shim today
# ===========================================================================


@pytest.mark.parametrize(
    "shim_name,symbol_name,old_meta",
    _COMPARISON_CASES,
    ids=[f"{s}.{n}" for s, n, _ in _COMPARISON_CASES],
)
def test_pre_refactor_public_symbol_still_importable(shim_name, symbol_name, old_meta):
    """Every name that was public BEFORE cf78545 must still be reachable.

    This catches silent API loss: a refactor deletes a function, the
    shim's `from X import *` doesn't complain, and users only notice
    months later when an old import suddenly raises ImportError.
    """
    if (shim_name, symbol_name) in _WAIVED_DROPS:
        pytest.skip(f"{shim_name}.{symbol_name}: drop explicitly waived")

    mod = importlib.import_module(shim_name)
    assert hasattr(mod, symbol_name), (
        f"{shim_name}.{symbol_name} existed before the refactor (at "
        f"{_PRE_REFACTOR_REF}) but is no longer accessible via the shim. "
        f"Either re-export it from the new module (so the shim's "
        f"`from X import *` picks it up) or, if the drop is intentional, "
        f"add ('{shim_name}', '{symbol_name}') to _WAIVED_DROPS with a "
        f"comment explaining why."
    )


# ===========================================================================
# 2. Function signatures are position-name-compatible
# ===========================================================================


_DEF_CASES = [
    (s, n, m) for s, n, m in _COMPARISON_CASES if m.get("kind") == "def"
]


@pytest.mark.parametrize(
    "shim_name,func_name,old_meta",
    _DEF_CASES,
    ids=[f"{s}.{n}" for s, n, _ in _DEF_CASES],
)
def test_pre_refactor_function_signature_unchanged(shim_name, func_name, old_meta):
    """Positional argument names and order must survive the refactor.

    Users have code like `run_orca(inp, out, timeout=60)`. If a refactor
    renames the first positional to `input_path`, that caller now fails
    with TypeError. Adding NEW params at the end or new keyword-only
    args is fine (backward-compatible); removing or renaming an
    existing positional is NOT.
    """
    if (shim_name, func_name) in _WAIVED_DROPS:
        pytest.skip(f"{shim_name}.{func_name}: change explicitly waived")

    mod = importlib.import_module(shim_name)
    func = getattr(mod, func_name, None)
    if func is None:
        pytest.fail(
            f"{shim_name}.{func_name} is gone — covered by the existence test"
        )
    if not callable(func):
        pytest.fail(f"{shim_name}.{func_name} exists but is not callable")

    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        pytest.skip(f"{shim_name}.{func_name} has no introspectable signature")

    # Current positional-compatible arg names, in order
    new_pos = [
        p.name
        for p in sig.parameters.values()
        if p.kind in (p.POSITIONAL_OR_KEYWORD, p.POSITIONAL_ONLY)
    ]

    old_pos = list(old_meta["pos_args"])
    # Drop `self`/`cls` if it slipped in (shouldn't — these are top-level)
    old_pos = [a for a in old_pos if a not in ("self", "cls")]

    # Rule 1: every old positional must appear at the same position
    for i, old_name in enumerate(old_pos):
        assert i < len(new_pos), (
            f"{shim_name}.{func_name}: used to have positional arg "
            f"#{i} ({old_name}) but the new signature has only "
            f"{len(new_pos)} positional args: {new_pos}. "
            f"This breaks callers that pass `{old_name}=...` or positional."
        )
        assert new_pos[i] == old_name, (
            f"{shim_name}.{func_name}: positional arg #{i} was "
            f"'{old_name}' before cf78545 but is '{new_pos[i]}' now. "
            f"Callers using kwarg syntax will fail with TypeError."
        )

    # Rule 2: no old keyword-only arg disappeared (new name would be missing)
    new_all = set(p.name for p in sig.parameters.values())
    lost_kwonly = set(old_meta["kwonly_args"]) - new_all
    if lost_kwonly and not old_meta.get("has_kwarg"):
        # If the function had **kwargs the removal may be absorbed, otherwise
        # every caller passing `foo=bar` will fail.
        pytest.fail(
            f"{shim_name}.{func_name}: keyword-only args "
            f"{sorted(lost_kwonly)} were dropped. Callers will fail "
            f"with TypeError."
        )


# ===========================================================================
# 3. Smoke test — the comparison harness itself works
# ===========================================================================


def test_pre_refactor_diff_harness_found_some_symbols():
    """If we suddenly have zero cases, something broke in the git-show
    logic and we're silently passing every test."""
    assert len(_COMPARISON_CASES) > 20, (
        f"Expected many pre-refactor symbols to compare, got "
        f"{len(_COMPARISON_CASES)}. The git-show call is probably "
        f"failing silently — tests are giving a false-green."
    )
