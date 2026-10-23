"""Clean up lists merged branches and stale worktrees, deletes them on request,
and writes what it deleted so every deletion can be undone.

Input: a checkout with local branches (merged into ``origin/main``, with a
merged or closed pull request, unmerged) and worktrees (clean and merged,
dirty, carrying commits on no branch, old). Output of
:func:`branch_cleanup.candidates`: what may be deleted, with the reason and
the count of commits not in ``origin/main``; and what is kept, with why.
:func:`branch_cleanup.delete` removes only selected candidates that still
hold: a branch with unmerged commits needs ``confirm_unmerged``; ``main``,
the checked-out branch, a branch checked out in any worktree, and anything
on the remote are never deleted. Each deletion is a line in
``~/.delfin/branch_cleanup.jsonl`` with the ref, its sha and the command
that brings it back.

Universal: repositories in ``tmp_path``; ``gh`` is a script written here.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
import time

import pytest

from delfin.dashboard import branch_cleanup as BC

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args, check=True) -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=check, capture_output=True, text=True)
    return out.stdout.strip()


def _commit(cwd, name):
    (cwd / name).write_text(name + "\n")
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", name)
    return _git(cwd, "rev-parse", "HEAD")


def _fake_gh(tmp_path, prs):
    script = tmp_path / "gh"
    script.write_text(
        f"#!{sys.executable}\nimport sys\n"
        f"sys.stdout.write({json.dumps(prs)!r})\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return str(script)


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
    _git(co, "branch", "merged-one")
    _git(co, "checkout", "-q", "-b", "unmerged")
    _commit(co, "u.txt")
    _git(co, "checkout", "-q", "-b", "squashed", "main")
    _commit(co, "s.txt")
    _git(co, "checkout", "-q", "main")
    return co


def _ids(cands):
    return {c.name for c in cands}


def test_lists_merged_and_pr_merged_branches_only(repo, tmp_path):
    gh = _fake_gh(tmp_path, [{"number": 3, "headRefName": "squashed",
                              "state": "MERGED", "statusCheckRollup": []}])
    cands, kept = BC.candidates(repo, gh=gh)
    names = _ids(c for c in cands if c.kind == "branch")
    assert names == {"merged-one", "squashed"}
    sq = next(c for c in cands if c.name == "squashed")
    assert sq.unmerged == 1 and "PR #3" in sq.reason
    mo = next(c for c in cands if c.name == "merged-one")
    assert mo.unmerged == 0
    assert "main" not in names and "unmerged" not in names


def test_deleting_a_merged_branch_logs_its_sha_and_the_undo_restores_it(repo):
    sha = _git(repo, "rev-parse", "merged-one")
    cands, _ = BC.candidates(repo, gh=None)
    report = BC.delete(repo, [c for c in cands if c.name == "merged-one"])
    assert [d["name"] for d in report.deleted] == ["merged-one"]
    assert _git(repo, "rev-parse", "--verify", "--quiet", "merged-one",
                check=False) == ""
    log = BC._log_path()
    assert str(log).startswith(os.environ["HOME"])
    rows = [json.loads(x) for x in log.read_text().splitlines()]
    assert rows[-1]["name"] == "merged-one" and rows[-1]["sha"] == sha
    undo = report.undo[0]
    assert sha in undo
    subprocess.run(undo, shell=True, check=True, env=_ENV)
    assert _git(repo, "rev-parse", "merged-one") == sha


def test_unmerged_commits_need_the_second_confirmation(repo, tmp_path):
    gh = _fake_gh(tmp_path, [{"number": 3, "headRefName": "squashed",
                              "state": "MERGED", "statusCheckRollup": []}])
    cands, _ = BC.candidates(repo, gh=gh)
    sq = [c for c in cands if c.name == "squashed"]
    report = BC.delete(repo, sq)
    assert report.deleted == []
    assert "second confirmation" in report.skipped[0][1]
    assert _git(repo, "rev-parse", "--verify", "squashed")
    report = BC.delete(repo, sq, confirm_unmerged=True)
    assert [d["name"] for d in report.deleted] == ["squashed"]


def test_main_current_and_held_branches_are_never_deleted(repo, tmp_path):
    wt = tmp_path / "wt-held"
    _git(repo, "worktree", "add", "-q", str(wt), "merged-one")
    (wt / "work.txt").write_text("in progress\n")       # keeps the worktree
    _git(repo, "checkout", "-q", "-b", "current-one", "main")
    cands, kept = BC.candidates(repo, gh=None)
    names = _ids(cands)
    assert "merged-one" not in names            # held by a kept worktree
    assert "current-one" not in names
    # forged candidates are refused at deletion time too
    forged = [BC.Candidate(kind="branch", name=n, sha=_git(repo, "rev-parse", n),
                           reason="forged") for n in ("main", "current-one",
                                                      "merged-one")]
    report = BC.delete(repo, forged, confirm_unmerged=True)
    assert report.deleted == []
    for n in ("main", "current-one", "merged-one"):
        assert _git(repo, "rev-parse", "--verify", n)


def test_nothing_on_the_remote_changes(repo):
    _git(repo, "push", "-q", "origin", "merged-one", "unmerged")
    before = _git(repo, "ls-remote", "origin")
    cands, _ = BC.candidates(repo, gh=None)
    BC.delete(repo, cands, confirm_unmerged=True)
    assert _git(repo, "ls-remote", "origin") == before
    assert _git(repo, "rev-parse", "--verify", "origin/merged-one")


def test_a_branch_that_moved_since_the_listing_is_skipped(repo):
    cands, _ = BC.candidates(repo, gh=None)
    _git(repo, "branch", "-f", "merged-one", "unmerged")
    report = BC.delete(repo, [c for c in cands if c.name == "merged-one"],
                       confirm_unmerged=True)
    assert report.deleted == [] and "moved" in report.skipped[0][1]


def test_stale_worktrees_are_listed_and_live_ones_kept(repo, tmp_path):
    merged = tmp_path / "wt-merged"
    _git(repo, "worktree", "add", "-q", str(merged), "merged-one")
    dirty = tmp_path / "wt-dirty"
    _git(repo, "worktree", "add", "-q", "--detach", str(dirty), "main")
    (dirty / "notes.txt").write_text("not committed\n")
    orphan = tmp_path / "wt-orphan"
    _git(repo, "worktree", "add", "-q", "--detach", str(orphan), "main")
    _commit(orphan, "lonely.txt")                      # a commit on no branch
    active = tmp_path / "wt-active"
    _git(repo, "worktree", "add", "-q", str(active), "unmerged")
    running = tmp_path / "wt-running"
    _git(repo, "worktree", "add", "-q", "--detach", str(running), "main")
    cands, kept = BC.candidates(repo, gh=None, protect_paths=[str(running)])
    paths = {c.path for c in cands if c.kind == "worktree"}
    assert str(merged) in paths
    kept_why = {k.name: k.why for k in kept}
    assert "uncommitted" in kept_why[str(dirty)]
    assert "no branch" in kept_why[str(orphan)]
    assert str(active) not in paths and str(active) in kept_why
    assert str(running) not in paths and "preview" in kept_why[str(running)]
    assert str(repo) not in paths
    # the branch the stale worktree holds is deleted after it, in one go
    br = next(c for c in cands if c.name == "merged-one")
    assert br.after == str(merged)
    report = BC.delete(repo, [c for c in cands
                              if c.path == str(merged) or c.name == "merged-one"])
    done = [(d["kind"], d["name"]) for d in report.deleted]
    assert done == [("worktree", str(merged)), ("branch", "merged-one")]
    assert not merged.exists()
    assert dirty.exists() and orphan.exists() and active.exists()
    assert any("worktree add" in u for u in report.undo)


def test_an_old_clean_worktree_is_stale(repo, tmp_path):
    old = tmp_path / "wt-old"
    _git(repo, "worktree", "add", "-q", str(old), "unmerged")
    cands, _ = BC.candidates(repo, gh=None)
    assert str(old) not in {c.path for c in cands}
    later = time.time() + 30 * 86400
    cands, _ = BC.candidates(repo, gh=None, stale_days=14, now=later)
    assert str(old) in {c.path for c in cands}


def test_a_worktree_whose_directory_is_gone_is_pruned(repo, tmp_path):
    gone = tmp_path / "wt-gone"
    _git(repo, "worktree", "add", "-q", "--detach", str(gone), "main")
    import shutil
    shutil.rmtree(gone)
    cands, _ = BC.candidates(repo, gh=None)
    c = [c for c in cands if c.path == str(gone)]
    assert c and "gone" in c[0].reason
    BC.delete(repo, c)
    assert str(gone) not in _git(repo, "worktree", "list")


def test_delete_branch_from_the_list_uses_the_same_rules(repo):
    ok, msg, report = BC.delete_branch(repo, "unmerged")
    assert not ok and "second confirmation" in msg
    ok, msg, report = BC.delete_branch(repo, "unmerged", confirm_unmerged=True)
    assert ok, msg
    assert report.undo and "unmerged" in report.undo[0]
    ok, msg, _ = BC.delete_branch(repo, "main", confirm_unmerged=True)
    assert not ok


def test_nothing_is_deleted_when_the_undo_log_cannot_be_written(repo):
    from pathlib import Path
    (Path(os.environ["HOME"]) / ".delfin").write_text("a file, not a dir\n")
    cands, _ = BC.candidates(repo, gh=None)
    report = BC.delete(repo, [c for c in cands if c.name == "merged-one"])
    assert report.deleted == []
    assert "undo log" in report.skipped[0][1]
    assert _git(repo, "rev-parse", "--verify", "merged-one")


def test_unmerged_branches_are_offered_only_when_asked_for(repo, tmp_path):
    """The route to delete a branch with commits of its own: Clean up with
    ``include_unmerged``; the row says how many commits are not in
    origin/main, and deletion still needs the second confirmation."""
    wt = tmp_path / "wt-held"
    _git(repo, "branch", "held-unmerged", "unmerged")
    _git(repo, "worktree", "add", "-q", str(wt), "held-unmerged")
    (wt / "work.txt").write_text("in progress\n")
    cands, _ = BC.candidates(repo, gh=None)
    assert "unmerged" not in _ids(cands)
    cands, _ = BC.candidates(repo, gh=None, include_unmerged=True)
    un = next(c for c in cands if c.name == "unmerged")
    assert un.unmerged == 1 and "not in" in un.reason
    assert "held-unmerged" not in _ids(cands)       # held: never offered
    assert "main" not in _ids(cands)
    report = BC.delete(repo, [un])
    assert report.deleted == [] and "second confirmation" in report.skipped[0][1]


def test_delete_refusals_answers_for_many_names_at_once(repo, tmp_path):
    wt = tmp_path / "wt-held"
    _git(repo, "worktree", "add", "-q", str(wt), "squashed")
    why = BC.delete_refusals(repo, ["main", "merged-one", "squashed",
                                    "unmerged", "origin/main", "no-such"])
    assert why["merged-one"] == "" and why["unmerged"] == ""
    assert why["main"] and why["squashed"] and why["origin/main"]
    assert why["no-such"]
    _git(repo, "checkout", "-q", "merged-one")
    assert BC.delete_refusals(repo, ["merged-one"])["merged-one"]
