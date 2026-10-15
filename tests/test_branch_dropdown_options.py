"""The dashboard branch dropdown lists ``main`` + this checkout's local
branches that carry work not yet on main, never the shared remote feature
branches and never a branch already merged. (The full rule, with pull
requests and the user's own pins, is in
test_the_branch_list_shows_what_is_not_on_main.py.)
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
    assert [v for _l, v in _compute_branch_options(repo)] == ["main"]


def _branch_with_commit(repo, name):
    _git(repo, "checkout", "-q", "-b", name, "main")
    (repo / (name.replace("/", "_") + ".txt")).write_text(name)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-m", name)
    _git(repo, "checkout", "-q", "main")


def test_lists_local_branches_main_first(repo):
    _branch_with_commit(repo, "contributor/balance-fix")
    _branch_with_commit(repo, "max/experiment")
    _git(repo, "branch", "already-merged")           # no commit of its own
    opts = _compute_branch_options(repo)
    values = [v for _label, v in opts]
    assert values[0] == "main"                       # main always first
    assert "contributor/balance-fix" in values            # contributor branch shows
    assert "max/experiment" in values                # the user's own shows
    assert "already-merged" not in values            # nothing to test there
    assert sorted(values[1:]) == ["contributor/balance-fix", "max/experiment"]


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
    assert [v for _l, v in _compute_branch_options(None)] == ["main"]
    assert [v for _l, v in _compute_branch_options("/nonexistent/path/xyz")] == ["main"]
