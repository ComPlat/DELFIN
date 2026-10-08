"""A missing file that exists on a sibling branch says how it reaches you.

The third shape from the field (2026-10-08): a session built its
component, found `delfin/chemdarwin_interface.py` absent from its own
worktree, and cleanly reverted 475 lines -- the module was on another
session's branch. "File not found" was true and said nothing. The shared
handover directory (#119) carries documents; code reaches a worktree
through main, and the clause says so, naming the worktree and the way.

Only existence of the same relative path in a sibling worktree of the
SAME repository is looked at -- never content. That a file exists on a
sibling branch is what `git branch -a` would tell anyone.
"""

from __future__ import annotations

import json
import subprocess

import pytest

from delfin.agent import api_client as A


def _git(where, *args):
    return subprocess.run(["git", "-C", str(where), *args],
                          capture_output=True, text=True, timeout=30)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    from delfin.agent import session_presence as P
    monkeypatch.setattr(P, "_git_cache", {})
    root = tmp_path / "project"
    root.mkdir()
    _git(root, "init", "-q")
    _git(root, "config", "user.email", "t@example.invalid")
    _git(root, "config", "user.name", "t")
    (root / "README.md").write_text("x\n")
    _git(root, "add", "README.md")
    _git(root, "commit", "-qm", "first")
    trees = root / ".delfin" / "worktrees"
    trees.mkdir(parents=True)
    for name in ("wt-a", "wt-b"):
        _git(root, "worktree", "add", "-q", "-b", f"session/{name}", str(trees / name))
    (trees / "wt-b" / "delfin").mkdir()
    (trees / "wt-b" / "delfin" / "chemdarwin_interface.py").write_text("# theirs\n")
    return root


def test_the_clause_names_the_worktree_and_the_way(repo):
    mine = repo / ".delfin" / "worktrees" / "wt-a"
    said = A._why_that_path_is_not_there("delfin/chemdarwin_interface.py", mine)
    assert "another session's worktree" in said, said
    assert "wt-b" in said and "session/wt-b" in said, said
    assert "pull request" in said and "merge origin/main" in said, said
    assert "handover directory" in said, "documents have another way, and it is named"


def test_a_path_in_no_sibling_says_nothing(repo):
    mine = repo / ".delfin" / "worktrees" / "wt-a"
    assert A._why_that_path_is_not_there("delfin/nowhere.py", mine) == ""


def test_an_absolute_path_is_not_searched(repo):
    """An absolute path names a place; the clause is for a relative path
    that the model expects beside it."""
    mine = repo / ".delfin" / "worktrees" / "wt-a"
    theirs = repo / ".delfin" / "worktrees" / "wt-b" / "delfin" / "chemdarwin_interface.py"
    assert "another session" not in A._why_that_path_is_not_there(str(theirs), mine)


def test_outside_a_repository_nothing_is_said(tmp_path, monkeypatch):
    from delfin.agent import session_presence as P
    monkeypatch.setattr(P, "_git_cache", {})
    loose = tmp_path / "loose"
    loose.mkdir()
    assert A._why_that_path_is_not_there("delfin/x.py", loose) == ""


def test_it_reaches_read_file(repo):
    """Through the real tool, with the workspace the gate hands over."""
    mine = repo / ".delfin" / "worktrees" / "wt-a"
    perms = A.KitToolPermissions(workspace=mine)
    perms.mode = "bypassPermissions"
    eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
    eng._permissions = perms
    out = json.loads(eng._execute_read_file(
        {"path": "delfin/chemdarwin_interface.py"}, perms))
    err = str(out.get("error") or "")
    assert "File not found" in err
    assert "wt-b" in err and "merge origin/main" in err, err


def test_the_existing_shapes_are_untouched():
    assert "abbreviation" in A._why_that_path_is_not_there("a/.../b.md")
    assert A._why_that_path_is_not_there("notes/readme.md") == ""
