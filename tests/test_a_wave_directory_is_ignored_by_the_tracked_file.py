"""Wave working directories are ignored by the TRACKED .gitignore.

A wave leaves its working directories under ``scratch/``: probe scripts, suite
logs and raw real-run transcripts. Those hold absolute paths inside somebody's
home directory and complete transcripts of tool calls. They were kept out of
git only by ``.git/info/exclude`` -- a file that is per-checkout and does not
travel to a clone -- while ``scratch/record_final_geometry.py`` is tracked on
main, so ``scratch/`` read as an ordinary committable directory. One
``git add -A`` in a fresh checkout, or in a landing checkout, published the
traces.

The rule now lives in ``.gitignore``. This pins two things so it cannot
silently regress:

* a path under ``scratch/`` is ignored, and the reason git names is
  ``.gitignore`` rather than ``.git/info/exclude``;
* a file under ``scratch/`` that is already tracked stays tracked -- an ignore
  rule must not hide a file the repository already follows, or
  ``record_final_geometry.py`` would quietly stop being versioned.

``git check-ignore`` is asked, not the file text: the text of a rule and the
effect of one are different things (a later ``!`` line, a global excludes file
or an ordering rule all change the effect).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(
    not (_REPO / ".git").exists(),
    reason="needs a git checkout to ask git itself",
)


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=_REPO, capture_output=True, text=True,
        check=False).stdout


def test_a_wave_working_directory_is_ignored_by_the_tracked_gitignore():
    """The rule that keeps transcripts out of the repository is versioned."""
    out = _git("check-ignore", "-v", "scratch/arms/zz-land/realrun.json")

    assert out.strip(), (
        "scratch/arms/... is not ignored at all; a `git add -A` would publish "
        "the raw run traces")
    source = out.split(":", 1)[0]
    assert source == ".gitignore", (
        f"the ignore rule comes from {source!r}, not the tracked .gitignore; "
        "a clone does not inherit that rule")


def test_an_already_tracked_scratch_file_stays_tracked():
    """An ignore rule never hides a file git already follows."""
    tracked = _git("ls-files", "scratch/").split()

    assert "scratch/record_final_geometry.py" in tracked, (
        "the tracked helper under scratch/ disappeared from the index; the "
        "ignore rule must not cost us a versioned file")
