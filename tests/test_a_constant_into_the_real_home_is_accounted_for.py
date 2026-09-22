"""Every module constant into ``~/.delfin`` is accounted for.

On 2026-09-21 a full suite run on a compute node left two ``$ ls``
confirmations "waiting at a terminal" for the operator: no session, the
host column naming the node. A test had asked through
``TerminalConfirmBroker`` and the question was published under the user's
REAL ``~/.delfin/terminal_confirmations`` -- the room's path is a module
constant resolved at import, and import-time constants are exactly what
the suite's redirect cannot reach (see the comment above the sink table
in ``delfin/agent/state_paths.py``: it was the same class of failure as
the 203 job records in the user's real locator index).

The tables are the product's, and they only work as a fence when nothing
stands outside them. This file walks ``delfin/`` at test time, finds every
module-level constant whose value is a path under ``Path.home() /
".delfin"``, and requires each one to be one of exactly three things:

  an entry in ``USER_STATE_SINKS``      a sink redirected by the suite
  an entry in ``USER_STATE_RESOLVERS``  a per-call resolver, ditto
  an entry in ``USER_STATE_REAL_HOME``  a deliberate reader of the real
                                        home, with a reason

A new module that grows such a constant therefore fails here first --
not as a question on somebody's terminal a day later.
"""

from __future__ import annotations

import ast
import importlib
from pathlib import Path

import pytest

from delfin.agent import state_paths

# The package root this file's assertions walk. Resolved from the imported
# product module so a relocated checkout needs no edit here.
_PKG_ROOT = Path(state_paths.__file__).resolve().parent.parent


def _calls_path_home(node: ast.AST) -> bool:
    """Whether the expression reads ``Path.home()`` at any depth."""
    for sub in ast.walk(node):
        if (isinstance(sub, ast.Call)
                and isinstance(sub.func, ast.Attribute)
                and sub.func.attr == "home"
                and isinstance(sub.func.value, ast.Name)
                and sub.func.value.id == "Path"):
            return True
    return False


def _names_delfin_dir(node: ast.AST) -> bool:
    """Whether the path chain names ``.delfin`` as a segment.

    Only ``~/.delfin`` is the user's state tree; ``Path.home() / "calc"``
    and its siblings are working directories with a meaning of their own
    and are outside what the tables govern.
    """
    for sub in ast.walk(node):
        if (isinstance(sub, ast.Constant)
                and isinstance(sub.value, str)
                and sub.value.strip() == ".delfin"):
            return True
    return False


def _module_constants_into_delfin() -> list[tuple[str, str]]:
    """Every ``(module, attribute)`` whose constant value sits in ~/.delfin.

    Read from the SOURCE, not from the imported attribute: the product
    also carries a function's worth of aliases (``_PID_PATH`` and friends
    built off ``_DELFIN_DIR``), and the constant that matters is the one
    at the head of the chain -- ``_DELFIN_DIR`` itself is an expression
    mentioning ``.delfin`` only through the constant it is built from,
    which is exactly why the walk finds the head, not the derived leaves
    (those are ``_DELFIN_DIR / "name"`` and name no ``.delfin`` of their
    own).
    """
    out: list[tuple[str, str]] = []
    for py in sorted(_PKG_ROOT.rglob("*.py")):
        rel = py.relative_to(_PKG_ROOT.parent)
        if not rel.as_posix().startswith("delfin/"):
            continue
        mod_name = rel.with_suffix("").as_posix().replace("/", ".")
        if mod_name.endswith(".__init__"):
            mod_name = mod_name[: -len(".__init__")]
        tree = ast.parse(py.read_text(encoding="utf-8"))
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            if not _names_delfin_dir(node.value) or not _calls_path_home(node.value):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    out.append((mod_name, target.id))
    return out


def test_every_constant_into_the_real_home_is_accounted_for():
    found = _module_constants_into_delfin()

    known_sinks = {(m, a) for m, a, _ in state_paths.USER_STATE_SINKS}
    known_resolvers = {
        (m, a) for m, a, _ in state_paths.USER_STATE_RESOLVERS}
    deliberate = {
        (m, a) for m, a, _ in state_paths.USER_STATE_REAL_HOME}

    unaccounted = [
        f"{m}.{a}" for m, a in found
        if (m, a) not in known_sinks
        and (m, a) not in known_resolvers
        and (m, a) not in deliberate
    ]
    assert not unaccounted, (
        "module constants into the user's real ~/.delfin that no table "
        "accounts for:\n  " + "\n  ".join(unaccounted) +
        "\neach needs an entry in USER_STATE_SINKS, USER_STATE_RESOLVERS "
        "or (with a reason) USER_STATE_REAL_HOME")


def test_every_real_home_entry_justifies_itself():
    """The deliberate list is not a dumping ground: an entry whose module
    no longer carries the constant is stale, and one whose reason is
    empty is not deliberate at all."""
    found = {(m, a) for m, a in _module_constants_into_delfin()}
    for mod_name, attr, reason in state_paths.USER_STATE_REAL_HOME:
        assert reason.strip(), (
            f"{mod_name}.{attr} reads the real home 'with intent' but "
            "carries no reason for that intent")
        assert (mod_name, attr) in found, (
            f"{mod_name}.{attr} is listed as reading the real home, but "
            "no such constant exists any more -- a stale entry")


def test_a_new_constant_outside_the_tables_is_caught(tmp_path):
    """The guard fires, not just passes: a file written into the package
    tree with a fresh constant must be found by the walk. Written to a
    tmp copy? No -- the walk reads the real tree, so this writes a small
    helper module next to the product and removes it again. That is a
    mutation of the package dir, so instead the check runs on an in-memory
    variant: a tree the walk would find, asserted by hand.
    """
    src = (
        "from pathlib import Path\n"
        "_PENDING_DIR = Path.home() / \".delfin\" / \"guard_probe\"\n"
    )
    tree = ast.parse(src)
    node = tree.body[1].value
    assert _names_delfin_dir(node) and _calls_path_home(node), (
        "the walk no longer recognizes a plain constant into ~/.delfin")


def test_the_measured_case_is_in_a_table():
    """The entry this file exists for, named explicitly, so removing it
    from its table fails here and not on somebody's terminal."""
    assert ("delfin.agent.terminal_confirm", "_PENDING_DIR") in {
        (m, a) for m, a, _ in state_paths.USER_STATE_SINKS
    } or ("delfin.agent.terminal_confirm", "_PENDING_DIR") in {
        (m, a) for m, a, _ in state_paths.USER_STATE_RESOLVERS
    }, "terminal_confirm._PENDING_DIR is not in any table -- the measured "
    "incident (2026-09-21) returns with the next suite run"
