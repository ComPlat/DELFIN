"""The dashboard header on a DELFIN that is installed, not checked out.

The header shows the git state of the installed package's directory. A pip
install (site-packages) and an unpacked deploy (scripts/pi) have no .git:
the header said "git status failed" on every page, and Pull and Rollback
pointed at a repository that was not there.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from delfin import dashboard as D


def _checkout(path: Path) -> Path:
    """A directory git itself calls a checkout.

    ``(path / ".git").mkdir()`` is not one: git answers "not a git
    repository" for an empty ``.git``, and a stray empty one on a
    workstation (measured: ``/tmp/.git``) made every path beneath it read as
    a checkout. The predicate under test reads the layout, so the stand-in
    has to have it.
    """
    path.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q"], cwd=str(path), check=True,
                   capture_output=True)
    return path


def test_an_installed_copy_is_said_to_be_one(tmp_path):
    (tmp_path / "delfin").mkdir()
    hint = D._installed_copy_hint(tmp_path)
    assert hint is not None
    assert hint["text"] == "installed" and hint["state"] == "ok"
    assert str(tmp_path) in hint["title"] and "nothing to pull" in hint["title"]


def test_a_checkout_is_left_to_its_git_status(tmp_path):
    _checkout(tmp_path)
    assert D._installed_copy_hint(tmp_path) is None


def test_a_worktree_checkout_counts_too(tmp_path):
    # In a `git worktree` checkout .git is a file pointing at the real one.
    (tmp_path / ".git").write_text("gitdir: /elsewhere/.git/worktrees/x\n")
    assert D._installed_copy_hint(tmp_path) is None
