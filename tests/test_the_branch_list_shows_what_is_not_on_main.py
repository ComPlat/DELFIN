"""The header's branch list shows what is not on main yet, and the user's own picks.

Input: a checkout whose ``origin/main`` contains some local branches and not
others, an optional ``gh`` that reports pull requests, and a per-user list
of pinned and hidden entries under ``~/.delfin``. Output: the dropdown
entries -- ``main`` first, then pinned entries, then branches with an open
pull request (number and CI state), then local branches with commits not in
``origin/main`` (ahead count). A branch fully contained in ``origin/main``
is not listed; neither is one whose pull request was merged. Without ``gh``
(missing, failing, no network) the list is built from git alone.

Universal: every repository is made in ``tmp_path`` with ``git init`` and a
bare ``origin`` beside it; ``gh`` is a script this test writes; no network.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import sys

import pytest

from delfin.dashboard import branch_list as BL

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args) -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True, text=True)
    return out.stdout.strip()


def _commit(cwd, name):
    (cwd / name).write_text(name + "\n")
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", name)
    return _git(cwd, "rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    _commit(co, "a.txt")
    _git(co, "push", "-q", "-u", "origin", "main")
    # merged: its commit is on origin/main
    _git(co, "branch", "merged-one")
    # ahead: two commits not on main
    _git(co, "checkout", "-q", "-b", "feat-ahead")
    _commit(co, "f1.txt")
    _commit(co, "f2.txt")
    _git(co, "checkout", "-q", "main")
    # a branch whose pull request was squash-merged: commits still ahead
    _git(co, "checkout", "-q", "-b", "squashed")
    _commit(co, "s.txt")
    _git(co, "checkout", "-q", "main")
    return co


def _fake_gh(tmp_path, prs, *, rc=0):
    script = tmp_path / "gh"
    payload = json.dumps(prs)
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        f"sys.stdout.write({payload!r})\n"
        f"sys.exit({rc})\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return str(script)


def _values(entries):
    return [e.ref for e in entries]


def test_merged_branches_are_hidden_and_ahead_ones_shown(repo):
    entries = BL.branch_entries(repo, gh=None)
    values = _values(entries)
    assert values[0] == "main"
    assert "merged-one" not in values
    assert "feat-ahead" in values
    ahead = {e.ref: e.ahead for e in entries}
    assert ahead["feat-ahead"] == 2
    label = next(e.label for e in entries if e.ref == "feat-ahead")
    assert "+2" in label


def test_without_gh_the_list_is_git_only(repo, tmp_path):
    missing = str(tmp_path / "no-such-gh")
    entries = BL.branch_entries(repo, gh=missing)
    assert "feat-ahead" in _values(entries)
    assert all(e.pr is None for e in entries)
    failing = _fake_gh(tmp_path, [], rc=1)
    assert "feat-ahead" in _values(BL.branch_entries(repo, gh=failing,
                                                     refresh=True))


def test_open_pr_shows_number_and_ci_and_merged_pr_hides(repo, tmp_path):
    gh = _fake_gh(tmp_path, [
        {"number": 7, "headRefName": "feat-ahead", "state": "OPEN",
         "url": "u7", "statusCheckRollup": [
             {"conclusion": "SUCCESS", "status": "COMPLETED"}]},
        {"number": 5, "headRefName": "squashed", "state": "MERGED",
         "url": "u5", "statusCheckRollup": []},
        {"number": 9, "headRefName": "someone-else", "state": "OPEN",
         "url": "u9", "statusCheckRollup": [
             {"conclusion": "FAILURE", "status": "COMPLETED"}]},
    ])
    entries = BL.branch_entries(repo, gh=gh, refresh=True)
    by = {e.ref: e for e in entries}
    assert by["feat-ahead"].pr == 7
    assert by["feat-ahead"].ci == "passing"
    assert "#7" in by["feat-ahead"].label
    # squash-merged: still "ahead" in git terms, hidden by its merged PR
    assert "squashed" not in by
    # an open PR without a local branch is offered through origin
    assert "origin/someone-else" in by
    assert by["origin/someone-else"].ci == "failing"


def test_the_shared_checkout_on_another_branch_stays_visible(repo):
    _git(repo, "checkout", "-q", "merged-one")
    values = _values(BL.branch_entries(repo, gh=None))
    assert "merged-one" in values


def test_filter_matches_substrings(repo):
    entries = BL.branch_entries(repo, gh=None)
    assert _values(BL.filter_entries(entries, "AHEAD")) == ["main", "feat-ahead"]
    assert _values(BL.filter_entries(entries, "")) == _values(entries)


# -- the user's own list -----------------------------------------------------

def test_pin_shows_a_hidden_rule_entry_and_hide_removes_one(repo):
    ok, msg, ref = BL.add_to_list(repo, "merged-one")
    assert ok, msg
    assert ref == "merged-one"
    values = _values(BL.branch_entries(repo, gh=None))
    assert "merged-one" in values
    BL.hide(repo, "feat-ahead")
    values = _values(BL.branch_entries(repo, gh=None))
    assert "feat-ahead" not in values
    # still in git, only not shown
    assert _git(repo, "rev-parse", "--verify", "feat-ahead")
    assert "feat-ahead" in BL.load_prefs(repo)["hidden"]
    BL.unhide(repo, "feat-ahead")
    assert "feat-ahead" in _values(BL.branch_entries(repo, gh=None))


def test_main_cannot_be_hidden(repo):
    ok, _msg = BL.hide(repo, "main")
    assert not ok
    assert _values(BL.branch_entries(repo, gh=None))[0] == "main"


def test_the_list_is_per_user_under_home_not_in_the_repo(repo, tmp_path):
    BL.add_to_list(repo, "merged-one")
    path = BL._prefs_path()
    assert str(path).startswith(str(tmp_path / "home"))
    assert path.is_file()
    assert oct(path.stat().st_mode & 0o777) == "0o600"
    assert _git(repo, "status", "--porcelain") == ""
    # another user (another HOME) has their own list
    os.environ["HOME"] = str(tmp_path / "other")
    (tmp_path / "other").mkdir()
    assert BL.load_prefs(repo) == {"pinned": [], "hidden": []}


def test_add_fetches_a_remote_only_branch(repo, tmp_path):
    other = tmp_path / "other-clone"
    _git(tmp_path, "clone", "-q", str(tmp_path / "origin.git"), str(other))
    _git(other, "checkout", "-q", "-b", "colleague/x")
    _commit(other, "c.txt")
    _git(other, "push", "-q", "origin", "colleague/x")
    ok, msg, ref = BL.add_to_list(repo, "colleague/x")
    assert ok, msg
    assert ref == "origin/colleague/x"
    assert _git(repo, "rev-parse", "origin/colleague/x") == _git(
        other, "rev-parse", "HEAD")
    assert "origin/colleague/x" in _values(BL.branch_entries(repo, gh=None))


def test_add_takes_a_tag_or_commit(repo):
    sha = _git(repo, "rev-parse", "feat-ahead~1")
    _git(repo, "tag", "v1", "feat-ahead")
    ok, _m, ref = BL.add_to_list(repo, "v1")
    assert ok and ref == "v1"
    ok, _m, ref = BL.add_to_list(repo, sha[:10])
    assert ok and ref == sha[:10]
    values = _values(BL.branch_entries(repo, gh=None))
    assert "v1" in values and sha[:10] in values


def test_add_refuses_what_does_not_exist_or_reads_as_an_option(repo):
    ok, msg, _ref = BL.add_to_list(repo, "no-such-branch")
    assert not ok and "not found" in msg
    ok, msg, _ref = BL.add_to_list(repo, "--upload-pack=touch /tmp/x")
    assert not ok


def test_new_branch_from_main_does_not_switch_the_checkout(repo, tmp_path):
    before = _git(repo, "branch", "--show-current")
    (repo / "dirty.txt").write_text("uncommitted\n")
    ok, msg, ref = BL.create_branch_from_main(repo, "max/new-idea")
    assert ok, msg
    assert ref == "max/new-idea"
    assert _git(repo, "branch", "--show-current") == before
    assert (repo / "dirty.txt").read_text() == "uncommitted\n"
    assert _git(repo, "rev-parse", "max/new-idea") == _git(
        repo, "rev-parse", "origin/main")
    assert "max/new-idea" in _values(BL.branch_entries(repo, gh=None))
    ok, msg, _ = BL.create_branch_from_main(repo, "max/new-idea")
    assert not ok and "exists" in msg
    ok, msg, _ = BL.create_branch_from_main(repo, "bad..name")
    assert not ok
