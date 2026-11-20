"""A directory named .git is not the same question as a git checkout.

Eight places asked ``(path / ".git").exists()`` to decide whether a path is
a repository checkout. An empty directory named ``.git`` answers yes, and
git answers "not a git repository" for the same path.

Measured on a workstation where something had left an empty ``/tmp/.git``
behind: every path under ``/tmp`` was classified as a repository checkout.
The refusal for a write into the calculation archive then named the
repository rule instead of the archive rule -- two suites failed on that
host and passed on CI, which is what made it visible at all. The refusal
still refused, so nothing was let through; the reason was wrong.

Input: a directory. Output: ``safe_git.is_checkout``, and
``safe_git.is_repository`` beside it.

Semantics pinned here: the layout predicate accepts what git accepts and
refuses what git refuses, over the cases that occur -- a repository, a
linked worktree, a submodule-style ``gitdir:`` file, a plain directory, an
empty ``.git`` directory, a ``.git`` that is neither, and a path that does
not exist. The two predicates are checked against each other on the same
trees, because two answers to one question drift apart when nothing holds
them together: ``is_repository`` asks git and costs a subprocess,
``is_checkout`` reads the layout and can be called for every parent of a
path, and only the second is cheap enough for the loops that need it.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import safe_git


def _git(*args: str, cwd: Path) -> None:
    subprocess.run(["git", *args], cwd=str(cwd), check=True,
                   capture_output=True)


@pytest.fixture
def repo(tmp_path):
    """A real repository with one commit, so a worktree can be added."""
    root = tmp_path / "repo"
    root.mkdir()
    _git("init", "-q", cwd=root)
    _git("config", "user.email", "t@example.org", cwd=root)
    _git("config", "user.name", "t", cwd=root)
    (root / "f.txt").write_text("x\n")
    _git("add", "f.txt", cwd=root)
    _git("commit", "-qm", "first", cwd=root)
    return root


def test_a_repository_is_a_checkout(repo):
    assert safe_git.is_checkout(repo) is True
    assert safe_git.is_repository(repo) is True


def test_a_linked_worktree_is_a_checkout(repo, tmp_path):
    """Its ``.git`` is a FILE holding ``gitdir:``, not a directory. Every
    swarm agent works in one of these, so a predicate that only knows the
    directory form calls their workspace no checkout."""
    wt = tmp_path / "wt"
    _git("worktree", "add", "-q", str(wt), "-b", "side", cwd=repo)
    assert (wt / ".git").is_file()
    assert safe_git.is_checkout(wt) is True
    assert safe_git.is_repository(wt) is True


def test_an_empty_dot_git_directory_is_not_a_checkout(tmp_path):
    """The case that was measured. git refuses it; so does this."""
    (tmp_path / ".git").mkdir()
    assert safe_git.is_checkout(tmp_path) is False
    assert safe_git.is_repository(tmp_path) is False


def test_a_dot_git_file_that_is_not_a_gitdir_pointer_is_not_a_checkout(tmp_path):
    (tmp_path / ".git").write_text("notes to self\n")
    assert safe_git.is_checkout(tmp_path) is False


def test_a_plain_directory_is_not_a_checkout(tmp_path):
    assert safe_git.is_checkout(tmp_path) is False
    assert safe_git.is_repository(tmp_path) is False


def test_a_path_that_does_not_exist_is_not_a_checkout(tmp_path):
    assert safe_git.is_checkout(tmp_path / "nowhere") is False
    assert safe_git.is_repository(tmp_path / "nowhere") is False


def test_a_subdirectory_of_a_checkout_is_not_the_top_of_one(repo):
    """Where the two predicates differ, and it is deliberate. git works
    from anywhere inside a checkout; the callers here walk UP the parents
    looking for the top, and a yes at every level would stop the walk at
    the first directory it looked at."""
    sub = repo / "pkg" / "inner"
    sub.mkdir(parents=True)
    assert safe_git.is_checkout(sub) is False
    assert safe_git.is_repository(sub) is True


def test_no_caller_asks_the_neighbouring_question_again(tmp_path):
    """A tree-wide check, because this is a shape and not one line.

    An importer grep would not catch a new call site; scanning the tree
    does. The predicate's own module is where the spelling is allowed to
    appear, in its docstring.
    """
    root = Path(__file__).resolve().parents[1]
    offenders = []
    for path in (root / "delfin").rglob("*.py"):
        if path.name == "safe_git.py":
            continue                 # where the spelling is quoted, in prose
        text = path.read_text(errors="replace")
        for n, line in enumerate(text.splitlines(), 1):
            if '".git")' in line or "'.git')" in line:
                if ".exists()" in line:
                    offenders.append(f"{path.relative_to(root)}:{n}")
    assert not offenders, (
        "a .git entry is being tested with .exists(), which an empty "
        f"directory passes; use safe_git.is_checkout: {offenders}")
