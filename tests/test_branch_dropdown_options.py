"""The dashboard branch-switch dropdown lists ``main`` + this checkout's LOCAL
branches (the ones you / a contributor like Jerome create with DELFIN) — never
the dozens of shared remote feature branches. Replaces the old hardcoded
``main / tools-and-workflows / GUPPY`` list.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.dashboard import _compute_branch_options


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    d = tmp_path / "checkout"
    d.mkdir()
    _git(d, "init", "-b", "main")
    _git(d, "config", "user.email", "t@t")
    _git(d, "config", "user.name", "t")
    (d / "f.txt").write_text("x")
    _git(d, "add", "-A")
    _git(d, "commit", "-m", "init")
    return d


def test_main_only_on_fresh_checkout(repo):
    # A clone with no feature branches yet → just main.
    assert _compute_branch_options(repo) == [("main", "main")]


def test_lists_local_branches_main_first(repo):
    _git(repo, "branch", "jerome/balance-fix")
    _git(repo, "branch", "max/experiment")
    opts = _compute_branch_options(repo)
    values = [v for _label, v in opts]
    assert values[0] == "main"                       # main always first
    assert "jerome/balance-fix" in values            # contributor branch shows
    assert "max/experiment" in values                # the user's own shows
    # deterministic: main, then the rest sorted
    assert values == ["main", "jerome/balance-fix", "max/experiment"]


def test_remote_branches_are_not_listed(repo, tmp_path):
    # A remote-tracking ref must NOT leak into the menu (we scope to refs/heads).
    origin = tmp_path / "origin.git"
    _git(repo, "init", "--bare", str(origin))  # noqa: not a real remote workflow,
    # create a remote-tracking ref by hand so we can assert it is filtered out
    _git(repo, "update-ref", "refs/remotes/origin/GUPPY", "HEAD")
    values = [v for _label, v in _compute_branch_options(repo)]
    assert "origin/GUPPY" not in values
    assert "GUPPY" not in values
    assert values == ["main"]


def test_missing_repo_dir_still_offers_main():
    # No repo / unreadable path → always still offer main so the user isn't stuck.
    assert _compute_branch_options(None) == [("main", "main")]
    assert _compute_branch_options("/nonexistent/path/xyz") == [("main", "main")]
