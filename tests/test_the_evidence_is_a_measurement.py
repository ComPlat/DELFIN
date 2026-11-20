"""A file under test-evidence/ carries a measurement, not prose.

``test-evidence/`` holds the measurements CI cannot reproduce — a run
against a real service, cluster or device, a browser, a pty, a cage, or a
before/after on a host CI does not have. The suite is its own evidence for
everything else, which is why a CI result does not belong there: a copy of
output that runs on every commit is a copy that goes stale, and a directory
of stale copies is a changelog.

The failure mode this guards is the one that kills such a directory: files
that read like a report and contain no number, so nobody can tell afterwards
whether anything was driven. Measured against the reason it exists: the
evidence for these changes previously lived in a pull-request body and in a
coordination file in another repository whose log had grown past two hundred
entries, and neither could be reached from a clone.

Input: the files under ``test-evidence/``. Output: whether each carries a
commit reference and at least one figure.

Semantics: the README is documentation and exempt, named by that role rather
than by its filename so a second guide is exempt for the same reason. No
check on wording, length or structure — a rule about prose would be a rule
about taste, and the point here is the number.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "test-evidence"

#: A commit, a range, or a branch tip the measurement is attributed to.
_COMMIT = re.compile(r"\b[0-9a-f]{7,40}\b")
#: A figure with or without separators: 12, 11,678, 2.4, 35%.
_FIGURE = re.compile(r"\b\d[\d,.]*\s*%?\b")


def _evidence_files() -> list[Path]:
    if not EVIDENCE.is_dir():
        return []
    return sorted(p for p in EVIDENCE.glob("*.md")
                  if not _is_guide(p))


def _is_guide(path: Path) -> bool:
    """A file that explains the directory rather than reporting a run.

    Decided by its first heading, not by its name: a second guide added
    later is exempt for the same reason, and a measurement cannot become
    exempt by being renamed.
    """
    head = path.read_text(errors="replace").lstrip().splitlines()[:1]
    first = (head[0] if head else "").strip().lower().lstrip("# ").strip()
    return first in {"test-evidence", "evidence", "readme"}


def test_the_directory_exists_and_holds_something():
    """A glob that matches nothing would make every assertion below pass
    while checking no file at all."""
    assert EVIDENCE.is_dir(), "test-evidence/ is missing"
    assert (EVIDENCE / "README.md").is_file(), "the directory has no guide"
    assert _evidence_files(), (
        "test-evidence/ holds no measurement, so this file asserts nothing")


@pytest.mark.parametrize("path", _evidence_files() or [None],
                         ids=lambda p: p.name if p is not None else "none")
def test_an_evidence_file_names_the_commit_it_measured(path):
    """A measurement without a commit cannot be attributed, and
    ``origin/main`` moves under it."""
    if path is None:
        pytest.skip("no evidence file to check")
    text = path.read_text(errors="replace")
    assert _COMMIT.search(text), (
        f"{path.name} names no commit; a number that cannot be attributed to "
        "a tree is not evidence")


@pytest.mark.parametrize("path", _evidence_files() or [None],
                         ids=lambda p: p.name if p is not None else "none")
def test_an_evidence_file_carries_a_figure(path):
    """Prose without a number is a report, and a directory of reports is a
    changelog."""
    if path is None:
        pytest.skip("no evidence file to check")
    text = path.read_text(errors="replace")
    figures = _FIGURE.findall(text)
    assert len(figures) >= 3, (
        f"{path.name} carries {len(figures)} figure(s); state what was "
        "measured, before and after")


@pytest.mark.parametrize("path", _evidence_files() or [None],
                         ids=lambda p: p.name if p is not None else "none")
def test_an_evidence_file_says_what_it_could_not_verify(path):
    """The part that makes it worth reading. A measurement that claims no
    limits is claiming the wrong thing, and every real run on this project
    has had a boundary -- a device that was not there, a service with no
    key, a platform this host is not."""
    if path is None:
        pytest.skip("no evidence file to check")
    text = path.read_text(errors="replace").lower()
    assert ("not verifiable" in text or "could not verify" in text
            or "not measured" in text), (
        f"{path.name} names no limit of its own measurement")
