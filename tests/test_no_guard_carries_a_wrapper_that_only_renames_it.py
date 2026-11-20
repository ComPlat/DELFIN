"""A guard is named by what it refuses, not by a wrapper that renames it.

Two module-level callables in the security core were unreachable: nothing in
``delfin/`` and nothing in ``tests/`` named them.

* ``principles_guard.expected_digest()`` returned the module constant and
  nothing else. The comparison it looks like it performs runs in
  ``check()`` (every pinned digest against the file's digest), and
  ``enforce()`` raises unless ``check()`` is ok; ``enforce()`` is called
  from ``cli.guard_principles_or_exit`` and ``AgentEngine.__init__``, so
  every entry point passes through it.
* ``refusal_memory._redirect_targets_safe()`` was an alias kept for a
  caller its docstring named (``api_client._bash_named_write_paths``) and
  which never existed -- that hook folds in ``_redirect_targets`` directly.
  The directory-change guard both names advertised lives inside
  ``_redirect_targets`` itself.

Neither added a refusal. Both were removed, and this test is what keeps a
thin forwarding wrapper from reappearing in either module: a second name
for one check is a second place for the answer to drift, and a name that
looks like a check but only forwards invites a caller that believes it has
gained a guarantee.

What is refused here: a module-level function in either guard module whose
whole body is one ``return`` of something a sibling already owns -- either a
module-level constant, or a call to another module-level function passing
through its own parameters unchanged. Nothing is said about functions with a
body of their own, which is the normal case.

This is a source check, so it measures the tree the interpreter loaded:
the paths are derived from the imported modules' ``__file__``, not from a
working directory.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from delfin.agent import principles_guard, refusal_memory

# The two modules that each held one such wrapper. Adding a module here is
# a decision about where a guard may be renamed twice; the test is written
# to run over any module given to it.
_GUARD_MODULES = (principles_guard, refusal_memory)


def _module_level_names(tree: ast.Module) -> tuple[set[str], set[str]]:
    """(function names, constant names) assigned at module level."""
    funcs, consts = set(), set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs.add(node.name)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    consts.add(target.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            consts.add(node.target.id)
    return funcs, consts


def _forwarding_wrappers(source: str) -> list[str]:
    """Names whose whole body is one ``return`` of what a sibling owns.

    A wrapper is reported when its body (docstring aside) is exactly one
    ``Return`` whose value is either

    * a bare module-level constant -- the zero-argument re-statement, or
    * a call to another module-level function, forwarding this function's
      own parameters in order -- the alias.

    A body that does any work of its own is not reported.
    """
    tree = ast.parse(source)
    funcs, consts = _module_level_names(tree)
    found: list[str] = []
    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        body = list(node.body)
        if body and isinstance(body[0], ast.Expr) and isinstance(
                getattr(body[0], "value", None), ast.Constant) and isinstance(
                body[0].value.value, str):
            body = body[1:]          # drop the docstring
        if len(body) != 1 or not isinstance(body[0], ast.Return):
            continue
        value = body[0].value
        if isinstance(value, ast.Name) and value.id in consts:
            found.append(f"{node.name}() -> {value.id}")
            continue
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Name):
            callee = value.func.id
            own = [a.arg for a in node.args.args]
            passed = [a.id for a in value.args if isinstance(a, ast.Name)]
            if (callee in funcs and callee != node.name
                    and passed == own and not value.keywords):
                found.append(f"{node.name}() -> {callee}()")
    return found


@pytest.mark.parametrize("module", _GUARD_MODULES,
                         ids=[m.__name__.rsplit(".", 1)[-1]
                              for m in _GUARD_MODULES])
def test_a_guard_module_carries_no_forwarding_wrapper(module):
    source = Path(module.__file__).read_text(encoding="utf-8")
    found = _forwarding_wrappers(source)
    assert found == [], (
        f"{module.__name__} carries a wrapper that only renames what it "
        f"forwards: {found}. Call the name it forwards to, or give it a "
        f"body that refuses something of its own -- a second name for one "
        f"check is a second place for the answer to drift.")


def test_the_checker_reports_both_shapes():
    """The detector is measured on source it must catch, so a green result
    above is not the detector having found nothing to look at."""
    source = (
        "CONST = 'x'\n"
        "def _real(cmd):\n"
        "    return [cmd]\n"
        "def _restated():\n"
        "    return CONST\n"
        "def _alias(cmd):\n"
        "    return _real(cmd)\n"
        "def _own_work(cmd):\n"
        "    return [cmd.upper()]\n"
    )
    assert _forwarding_wrappers(source) == [
        "_restated() -> CONST", "_alias() -> _real()"]


def test_the_checker_leaves_a_body_of_its_own_alone():
    """A function that does anything besides forward is not a wrapper: a
    guard with a condition must not be reported, or the check would punish
    the shape it is trying to encourage."""
    source = (
        "def _guarded(cmd):\n"
        "    if not cmd:\n"
        "        return []\n"
        "    return [cmd]\n"
    )
    assert _forwarding_wrappers(source) == []
