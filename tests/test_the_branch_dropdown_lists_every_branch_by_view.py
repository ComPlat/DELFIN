"""The branch dropdown's views: Mine, All, Hidden, Merged, Remote, and a typed ref.

Input: a checkout made here with a bare origin -- local branches with and
without commits of their own, a branch only origin has, a tag, a pinned
branch and a hidden one -- and the per-user list under a temporary HOME.
Output of :func:`branch_list.all_rows`: one row per branch that exists
(local, and ``origin/<x>`` when no local ``x`` exists), with its last commit
date, author, whether ``origin/main`` contains it, and whether it is pinned,
hidden or in the default view. :func:`branch_list.view_rows` selects one
view and narrows it by a case-insensitive substring of the name or author.
:func:`branch_list.lookup_row` turns a typed tag, commit or unknown branch
name into a candidate row. Nothing here changes git.
"""

from __future__ import annotations

import json
import os
import subprocess

import pytest

from delfin.dashboard import branch_list as BL

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args, check=True, author="t") -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", f"user.name={author}", *args],
        cwd=str(cwd), env=_ENV, check=check, capture_output=True, text=True)
    return out.stdout.strip()


def _commit(cwd, name, author="t"):
    (cwd / name).write_text(name + "\n")
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", name, author=author)


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
    _git(co, "branch", "old/merged")
    _git(co, "branch", "mine/pinned")
    _git(co, "checkout", "-q", "-b", "feat")
    _commit(co, "f.txt", author="Lena S")
    _git(co, "checkout", "-q", "-b", "hidden/one")
    _commit(co, "h.txt")
    _git(co, "checkout", "-q", "-b", "remote/only")
    _commit(co, "r.txt", author="Tim B")
    _git(co, "push", "-q", "origin", "remote/only")
    _git(co, "checkout", "-q", "main")
    _git(co, "branch", "-q", "-D", "remote/only")
    _git(co, "tag", "v1", "main")
    _git(co, "fetch", "-q", "origin")
    BL.pin(co, "mine/pinned")
    BL.hide(co, "hidden/one")
    before = (_git(co, "for-each-ref"), _git(co, "status", "--porcelain"))
    yield co
    assert (_git(co, "for-each-ref"), _git(co, "status", "--porcelain")) == before


def _refs(rows):
    return [r.ref for r in rows]


def test_all_rows_cover_local_and_remote_only_with_date_and_author(repo):
    rows = {r.ref: r for r in BL.all_rows(repo, gh=None)}
    assert set(rows) >= {"main", "feat", "old/merged", "mine/pinned",
                         "hidden/one", "origin/remote/only"}
    assert "origin/main" not in rows and "origin/HEAD" not in rows
    assert rows["origin/remote/only"].kind == "remote"
    assert rows["origin/remote/only"].author == "Tim B"
    assert rows["feat"].author == "Lena S"
    assert len(rows["feat"].date) == 10 and rows["feat"].date[4] == "-"
    assert rows["old/merged"].merged and not rows["feat"].merged
    assert rows["mine/pinned"].pinned and rows["mine/pinned"].mine
    assert rows["hidden/one"].hidden and not rows["hidden/one"].mine
    assert rows["feat"].mine and rows["feat"].ahead == 1
    assert not rows["old/merged"].mine


def test_each_view_selects_its_rows(repo):
    rows = BL.all_rows(repo, gh=None)
    assert _refs(BL.view_rows(rows, "mine")) == ["main", "mine/pinned", "feat"]
    assert set(_refs(BL.view_rows(rows, "all"))) == {r.ref for r in rows}
    assert _refs(BL.view_rows(rows, "hidden")) == ["hidden/one"]
    assert set(_refs(BL.view_rows(rows, "merged"))) == {"old/merged", "mine/pinned"}
    assert _refs(BL.view_rows(rows, "remote")) == ["origin/remote/only"]


def test_search_narrows_by_name_or_author_case_insensitive(repo):
    rows = BL.all_rows(repo, gh=None)
    assert _refs(BL.view_rows(rows, "all", "REMOTE")) == ["origin/remote/only"]
    assert _refs(BL.view_rows(rows, "all", "lena")) == ["feat"]
    # main stays reachable in the default view whatever is typed
    assert _refs(BL.view_rows(rows, "mine", "zzz")) == ["main"]
    assert _refs(BL.view_rows(rows, "all", "zzz")) == []


def test_a_typed_tag_commit_or_unknown_name_becomes_a_candidate(repo):
    rows = BL.all_rows(repo, gh=None)
    known = {r.ref for r in rows}
    tag = BL.lookup_row(repo, "v1", known)
    assert tag is not None and tag.kind == "tag" and tag.ref == "v1"
    sha = _git(repo, "rev-parse", "feat")[:10]
    commit = BL.lookup_row(repo, sha, known)
    assert commit is not None and commit.kind == "commit"
    unknown = BL.lookup_row(repo, "someone/new-work", known)
    assert unknown is not None and unknown.kind == "lookup"
    assert BL.lookup_row(repo, "feat", known) is None        # already listed
    assert BL.lookup_row(repo, "--upload-pack=x", known) is None
    assert BL.lookup_row(repo, "", known) is None


def test_rows_carry_the_pull_request_and_its_ci_state(repo, tmp_path):
    gh = tmp_path / "gh"
    gh.write_text("#!/bin/sh\necho '" + json.dumps(
        [{"number": 7, "headRefName": "feat", "state": "OPEN", "url": "",
          "statusCheckRollup": [{"conclusion": "SUCCESS"}]}]) + "'\n")
    gh.chmod(0o755)
    rows = {r.ref: r for r in BL.all_rows(repo, gh=str(gh), refresh=True)}
    assert rows["feat"].pr == 7 and rows["feat"].ci == "passing"
