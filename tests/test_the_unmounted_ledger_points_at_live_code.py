"""The unmounted-capability ledger must point at code that is there.

`tests/fixtures/unmounted_capabilities.txt` is the ledger that decides
whether a capability counts as unreachable, and its comment blocks are where
the reason for each row lives -- they are read as fact by the next person or
agent working the row. Two ways it can rot, both measured on 2026-11-25:

  * a comment names a file that was later DELETED. `shared_tmp.py` and its
    five test files went in a92ed6e36; the `map_to_host` and
    `bwrap_substitution` notes still describe `tests/test_shared_tmp_wiring.py`
    as a live strict-xfail contract waiting for wiring. A wave then plans
    work against a contract that does not exist.
  * a row's `file:line` anchor drifts as the file grows, so the pointer lands
    on an unrelated line (6 of 9 rows did: `role_tool_surface_report` said
    12107, the definition is at 12667).

Input: the ledger file and the source tree -- no model, no cage, nothing
installed. Output: which pointers do not resolve.

A path is checked only when it is written root-relative (`delfin/...`,
`tests/...`): a bare `build.py` or `engine.py` in prose is an example command,
not a pointer, and `/tmp/...` / `../x.txt` are shell operands inside a quoted
command. That keeps the check to the one claim it can actually make.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = Path(__file__).resolve().parent / "fixtures" / "unmounted_capabilities.txt"

#: A root-relative path to a tracked source file, as the notes cite one.
_PATH_RE = re.compile(
    r"\b((?:delfin|tests|docs|scripts|tools)/[A-Za-z0-9_./-]*\.(?:py|md|txt))\b")


def _commented_paths() -> list[str]:
    out: list[str] = []
    for line in LEDGER.read_text().splitlines():
        if not line.lstrip().startswith("#"):
            continue
        out += _PATH_RE.findall(line)
    return out


def _rows() -> list[tuple[str, str, str, int]]:
    """``(tier, name, file, line)`` for every data row."""
    rows = []
    for line in LEDGER.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tier, name, where = line.split("\t")
        path, _, lineno = where.rpartition(":")
        rows.append((tier, name, path, int(lineno)))
    return rows


def _definition_lines(path: Path, name: str) -> set[int]:
    """Lines where ``name`` is bound at module level or as a class method."""
    tree = ast.parse(path.read_text())
    found: set[int] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            if node.name == name:
                found.add(node.lineno)
            if isinstance(node, ast.ClassDef):
                for sub in node.body:
                    if getattr(sub, "name", None) == name and isinstance(
                            sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        found.add(sub.lineno)
        elif isinstance(node, ast.Assign):
            found.update(t.id for t in node.targets
                         if isinstance(t, ast.Name) and t.id == name)
        elif isinstance(node, ast.AnnAssign):
            if getattr(node.target, "id", None) == name:
                found.add(node.lineno)
    return found


def test_the_ledger_exists_and_holds_rows():
    """Control on both sides of the fix: the file is there and readable, so
    the two tests below assert something."""
    assert LEDGER.is_file(), LEDGER
    assert _rows(), "the ledger holds no data rows"


def test_a_note_cites_a_file_that_still_exists():
    """A reason that points at a deleted file cannot be checked, and reads
    as a live obligation. Failing input: a comment naming a test file that a
    later commit deleted."""
    dead = sorted({p for p in _commented_paths() if not (ROOT / p).exists()})
    assert not dead, (
        f"{len(dead)} cited path(s) no longer exist: "
        + ", ".join(dead)
        + " -- say what happened to them instead of leaving the pointer.")


@pytest.mark.parametrize("tier,name,path,lineno", _rows(),
                         ids=[r[1] for r in _rows()])
def test_a_row_anchor_lands_on_its_definition(tier, name, path, lineno):
    """The `file:line` column is a pointer, so it must resolve. Failing
    input: the definition moved when the file grew, and the anchor now names
    an unrelated line."""
    target = ROOT / path
    assert target.is_file(), f"{path} does not exist"
    lines = target.read_text().splitlines()
    assert 0 < lineno <= len(lines), (
        f"{path}:{lineno} is past the end of a {len(lines)}-line file; "
        f"definition line(s) now: {sorted(_definition_lines(target, name))}")
    got = _definition_lines(target, name)
    assert lineno in got, (
        f"{name} is not defined at {path}:{lineno} "
        f"(that line reads {lines[lineno - 1].strip()!r}); "
        f"it is defined at {sorted(got) or 'no module-level binding'}")
