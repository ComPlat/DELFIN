"""A tracked file may not point a reader at a path git cannot produce.

``delfin/agent/commit_attribution.py`` and its two contract tests name the
patch that wires the executor half of the attribution decision, and tell the
reader to apply it. ``scratch/`` is gitignored (``.gitignore:178``), so a path
under it exists in no clone and no merge of this branch: the one instruction
those files give cannot be followed from the repository, and the switch they
promise stays a promise.

Scope is the attribution surface this module owns -- the module and the tests
that pin its contract -- not the whole tree: other ``*.patch`` strings in the
repository are historical prose (a CHANGELOG entry for a patch removed in
June) or ``unittest.mock.patch``, and a tree-wide scan would report those
whether or not this change lands.

The shape is checked, not the name: every ``*.patch`` path one of these files
mentions must resolve from the repository root and must be a unified diff. A
path git cannot produce fails whether it is called
``h1-trailer-off-for-members.patch`` or anything else.

Input: the text of the attribution module and its contract tests, and the
working tree. Output: the references whose target is missing or not a diff.
"""

from __future__ import annotations

import pathlib
import re

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]

#: The files that carry the executor half's contract, relative to the root.
_SURFACE = (
    "delfin/agent/commit_attribution.py",
    "tests/test_the_executor_asks_the_same_attribution_question.py",
    "tests/test_the_swarm_member_decides_its_commit_attribution.py",
)

#: A repository-relative path ending in ``.patch`` as those files spell it.
_PATCH_REF = re.compile(r"(?<![\w/.-])([\w./-]+\.patch)\b")

#: Forms that describe a patch without naming a file in this repository: the
#: landing rule quotes ``scratch/<file>.patch`` as a placeholder.
_NOT_A_PATH = frozenset({"<file>.patch", "*.patch", "scratch/<file>.patch"})


def _patch_references() -> list[tuple[str, int, str]]:
    """``(file, line number, the path it names)`` for every reference."""
    found: list[tuple[str, int, str]] = []
    for rel in _SURFACE:
        path = _ROOT / rel
        assert path.is_file(), f"{rel} is part of the contract surface but is absent"
        for lineno, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            for match in _PATCH_REF.finditer(line):
                ref = match.group(1)
                if ref in _NOT_A_PATH or "<file>" in ref:
                    continue
                found.append((rel, lineno, ref))
    return found


def test_every_patch_the_attribution_contract_names_is_in_the_repository():
    """The pointer a reader is told to follow must resolve in a clone."""
    refs = _patch_references()
    assert refs, (
        "no attribution-surface file names a patch any more; either the "
        "executor half landed (drop this test) or the instruction was lost")
    missing = [f"{rel}:{lineno} -> {ref}"
               for rel, lineno, ref in refs if not (_ROOT / ref).is_file()]
    assert not missing, (
        "tracked text names a patch that is not in the repository, so a clone "
        "cannot follow the instruction to apply it:\n  " + "\n  ".join(missing))


def test_the_patch_those_files_name_is_a_unified_diff():
    """And the file that resolves must be the diff it claims to be."""
    for ref in sorted({r for _rel, _lineno, r in _patch_references()}):
        text = (_ROOT / ref).read_text(errors="replace")
        assert text.lstrip().startswith(("--- ", "diff ")), (
            f"{ref} is named as a patch but is not a unified diff")
        assert "+++ " in text, f"{ref} carries no diff hunk header"


if __name__ == "__main__":  # pragma: no cover
    pytest.main([__file__, "-q"])
