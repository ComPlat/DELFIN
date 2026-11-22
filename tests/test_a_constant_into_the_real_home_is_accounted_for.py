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

import _package_ast
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
        tree = _package_ast.parse(py)
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


# --- the other half: a real-home constant must not be DESTROYED -----------
#
# The table above answers "is this path accounted for?". It does not ask
# what the module DOES with it, and the two questions came apart on
# 2026-11-21: model_capabilities._CACHE_PATH is on USER_STATE_REAL_HOME as
# "read-mostly ... a redirect would hide the real cache from a run that
# reads it", and the module called `_CACHE_PATH.unlink()`. A suite run
# removed the operator's accumulated cache through it -- the only caller
# anywhere was a test that had not redirected the path.
#
# USER_STATE_REAL_HOME states the rule itself ("a constant on this list
# that starts WRITING as a matter of course belongs in USER_STATE_SINKS
# instead"); nothing checked it. Removal is the half that cannot be
# undone, so it is the half checked here.

_DESTRUCTIVE_METHODS = ("unlink", "rmdir")


def _destructive_calls_on(attr: str, tree: ast.AST) -> list[str]:
    """``attr.unlink()`` / ``attr.rmdir()`` / ``rmtree(attr)`` in *tree*.

    Matched on the NAME the table names, so a path derived from it (a
    leaf, which has its own table entry) is not this check's business.
    """
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Attribute):
            if (func.attr in _DESTRUCTIVE_METHODS
                    and isinstance(func.value, ast.Name)
                    and func.value.id == attr):
                out.append(f"{attr}.{func.attr}() at line {node.lineno}")
            elif func.attr == "rmtree":
                for arg in node.args:
                    if isinstance(arg, ast.Name) and arg.id == attr:
                        out.append(f"rmtree({attr}) at line {node.lineno}")
    return out


def test_no_real_home_constant_is_removed_by_the_product():
    for mod_name, attr, _reason in state_paths.USER_STATE_REAL_HOME:
        path = _PKG_ROOT.parent / (mod_name.replace(".", "/") + ".py")
        assert path.is_file(), f"{mod_name} has no file at {path}"
        hits = _destructive_calls_on(attr, ast.parse(path.read_text()))
        assert not hits, (
            f"{mod_name}.{attr} points at the operator's real ~/.delfin "
            f"and the module removes it: {'; '.join(hits)}. A read-mostly "
            "store is not a cache to drop: either the removal becomes "
            "opt-in (a parameter the caller asks for), or the constant "
            "belongs in USER_STATE_SINKS where the suite redirects it.")


def test_the_removal_check_sees_a_planted_unlink():
    """The guard fires, not just passes -- the measured shape, verbatim."""
    tree = ast.parse(
        "def clear_cache():\n"
        "    try:\n"
        "        _CACHE_PATH.unlink()\n"
        "    except Exception:\n"
        "        pass\n")
    assert _destructive_calls_on("_CACHE_PATH", tree)
    # And the two ways it must stay quiet: another name, and a read.
    assert not _destructive_calls_on("_OTHER_PATH", tree)
    assert not _destructive_calls_on(
        "_CACHE_PATH", ast.parse("_CACHE_PATH.read_text()\n"))


def test_the_measured_case_keeps_the_file_it_used_to_delete():
    """The entry this half exists for, driven rather than read: clearing
    the model-capability cache leaves the stored file where it is.

    The path is redirected into tmp_path, so this test cannot repeat the
    incident it is about.
    """
    import json

    from delfin.agent import model_capabilities as mc

    with pytest.MonkeyPatch.context() as mp, \
            __import__("tempfile").TemporaryDirectory() as tmp:
        stored = Path(tmp) / "model_capabilities_cache.json"
        stored.write_text(json.dumps({"_v": mc._CACHE_VERSION,
                                      "entries": {}}))
        mp.setattr(mc, "_CACHE_PATH", stored)
        mc.clear_cache()
        assert stored.is_file(), (
            "clear_cache removed the stored cache; on an unredirected "
            "path that is the operator's accumulated file")
