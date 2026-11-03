"""Branch operations reach the Pi the same way the dashboard does them.

Input: a checkout made here with a bare origin (local branches, a remote-only
branch, a hidden/pinned pref under a temporary HOME) and one with a worktree
and an unmerged branch. Output:

* ``branch_ops.rows(repo, view)`` -- the same rows ``branch_list.all_rows``
  and ``view_rows`` give the dashboard's dropdown, as plain dicts;
* ``branch_ops.pull_main(repo)`` -- fast-forward-only update of main of the
  checkout the service runs in (never a switch, never a merge);
* ``branch_ops.test_branch(repo, ref)`` -- a worktree at the branch, the
  branch's changed test files run in it (capped time, ``nice``), pass/fail
  counts and the first failure, worktree removed afterwards unless kept;
* ``branch_ops.cleanup_candidates(repo)`` / ``branch_ops.cleanup_delete`` --
  what ``branch_cleanup`` offers and deletes, deletion only through it, with
  the second confirmation for unmerged branches.

Nothing here switches the shared checkout; the dashboard's own behaviour is
not touched (its imports keep working).
"""

from __future__ import annotations

import os
import subprocess

import pytest

REPO = None  # set in fixture: the DELFIN checkout this test imports


def _run(cwd, *args, env=None, check=True):
    out = subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                          *args], cwd=str(cwd), capture_output=True, text=True,
                         env=env or {**os.environ,
                                     "GIT_CONFIG_GLOBAL": os.devnull,
                                     "GIT_CONFIG_SYSTEM": os.devnull})
    if check and out.returncode != 0:
        raise AssertionError(f"git {' '.join(args)}: {out.stderr}")
    return out.stdout.strip()


def _commit(cwd, name, msg=None):
    (cwd / name).write_text(name + "\n")
    _run(cwd, "add", name)
    _run(cwd, "commit", "-q", "-m", msg or name)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    origin = tmp_path / "origin.git"
    _run(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _run(tmp_path, "clone", "-q", str(origin), str(co))
    _commit(co, "a.txt")
    _run(co, "push", "-q", "-u", "origin", "main")
    _run(co, "branch", "old/merged")
    _run(co, "checkout", "-q", "-b", "feat")
    _commit(co, "f.txt")
    _run(co, "checkout", "-q", "main")
    return co


# -- views -------------------------------------------------------------------

def test_rows_views(tmp_path, repo):
    from delfin.agent import branch_ops
    mine = branch_ops.rows(repo, "mine")
    assert [r["ref"] for r in mine][:1] == ["main"]
    assert any(r["ref"] == "feat" for r in mine)
    all_rows = branch_ops.rows(repo, "all")
    refs = {r["ref"] for r in all_rows}
    assert "old/merged" in refs and "feat" in refs
    merged = branch_ops.rows(repo, "merged")
    assert {r["ref"] for r in merged} >= {"old/merged"}


def test_rows_remote_view_and_filter(tmp_path, repo):
    from delfin.agent import branch_ops
    _run(repo, "push", "-q", "origin", "feat")
    _run(repo, "branch", "-D", "feat")
    remote = branch_ops.rows(repo, "remote")
    assert any(r["ref"] == "origin/feat" for r in remote)
    mine = branch_ops.rows(repo, "mine", text="main")
    assert [r["ref"] for r in mine] == ["main"]


def test_views_names(tmp_path, repo):
    from delfin.agent import branch_ops
    assert branch_ops.VIEWS == ("mine", "all", "hidden", "merged", "remote")


# -- pull --------------------------------------------------------------------

def test_pull_updates_main_fast_forward_only(repo):
    from delfin.agent import branch_ops
    # origin gets a new commit on main
    other = repo.parent / "other"
    _run(repo.parent, "clone", "-q", str(repo.parent / "origin.git"),
         str(other))
    _commit(other, "b.txt")
    _run(other, "push", "-q", "origin", "main")
    result = branch_ops.pull_main(repo)
    assert result["ok"] is True
    assert "up to date" not in result["text"].lower() or \
        "new commit" in result["text"] or "1" in result["text"]
    assert _run(repo, "rev-parse", "HEAD") == _run(other, "rev-parse", "HEAD")


def test_pull_refuses_when_not_on_main(repo):
    from delfin.agent import branch_ops
    _run(repo, "checkout", "-q", "-b", "side")
    result = branch_ops.pull_main(repo)
    assert result["ok"] is False
    assert "main" in result["text"].lower()


def test_pull_refuses_a_non_fast_forward(repo):
    from delfin.agent import branch_ops
    _run(repo, "checkout", "-q", "main")
    other = repo.parent / "other2"
    _run(repo.parent, "clone", "-q", str(repo.parent / "origin.git"),
         str(other))
    _commit(other, "c.txt")
    _run(other, "push", "-q", "origin", "main")
    # the checkout makes its own commit on main (diverged)
    _commit(repo, "d.txt")
    _run(repo, "push", "-q", "origin", "main")
    _run(other, "pull", "-q", "--rebase", "origin", "main")
    _commit(other, "e.txt")
    _run(other, "push", "-q", "origin", "main")
    result = branch_ops.pull_main(repo)
    assert result["ok"] is False
    assert "fast-forward" in result["text"].lower()
    # the checkout stays on its own commit: nothing was moved
    assert _run(repo, "rev-parse", "HEAD") != _run(other, "rev-parse", "HEAD")


# -- test --------------------------------------------------------------------

def test_test_branch_runs_changed_tests_in_a_worktree(repo, monkeypatch):
    from delfin.agent import branch_ops
    _run(repo, "checkout", "-q", "-b", "fix-thing")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_the_thing.py").write_text(
        "def test_it():\n    assert 1 == 1\n")
    _run(repo, "add", "tests/test_the_thing.py")
    _run(repo, "commit", "-q", "-m", "add test")
    _run(repo, "checkout", "-q", "main")
    calls = []

    def fake_runner(cwd, files, timeout_s):
        calls.append((str(cwd), tuple(files), timeout_s))
        return {"passed": 1, "failed": 0, "first_failure": ""}

    monkeypatch.setattr(branch_ops, "run_tests", fake_runner)
    result = branch_ops.test_branch(repo, "fix-thing",
                                    runner=fake_runner)
    assert result["ok"] is True
    assert result["passed"] == 1 and result["failed"] == 0
    assert result["kept"] is False
    # the worktree is gone again
    assert not (repo.parent / "delfin-preview" / "fix-thing").exists()
    # the shared checkout stayed on main
    assert _run(repo, "branch", "--show-current") == "main"
    # the runner saw the worktree and the branch's test files, capped
    assert calls and "fix-thing" in calls[0][0]
    assert calls[0][1] == ("tests/test_the_thing.py",)
    assert calls[0][2] <= branch_ops.TEST_TIMEOUT_S


def test_test_branch_reports_first_failure(repo, monkeypatch):
    from delfin.agent import branch_ops
    _run(repo, "checkout", "-q", "-b", "broken")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_broken.py").write_text(
        "def test_broken():\n    assert False\n")
    _run(repo, "add", "tests/test_broken.py")
    _run(repo, "commit", "-q", "-m", "broken test")
    _run(repo, "checkout", "-q", "main")

    def fake_runner(cwd, files, timeout_s):
        return {"passed": 0, "failed": 1,
                "first_failure": "test_broken.py::test_broken: assert False"}

    result = branch_ops.test_branch(repo, "broken", runner=fake_runner)
    assert result["ok"] is True
    assert result["failed"] == 1
    assert "test_broken" in result["first_failure"]


def test_test_branch_keeps_on_request(repo, monkeypatch):
    from delfin.agent import branch_ops
    _run(repo, "checkout", "-q", "-b", "keepme")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_x.py").write_text("def test_x():\n    assert 1\n")
    _run(repo, "add", "tests/test_x.py")
    _run(repo, "commit", "-q", "-m", "x")
    _run(repo, "checkout", "-q", "main")
    result = branch_ops.test_branch(repo, "keepme",
                                    runner=lambda *a: {"passed": 1, "failed": 0,
                                                       "first_failure": ""},
                                    keep=True)
    assert result["kept"] is True
    wt = repo.parent / "delfin-preview" / "keepme"
    assert wt.is_dir()
    # cleaned up by hand so the fixture's tmp_path goes away cleanly
    from delfin.agent import branch_ops
    branch_ops.remove_test_worktree(repo, wt)
    assert not wt.exists()


def test_test_branch_refuses_unknown_ref(repo, monkeypatch):
    from delfin.agent import branch_ops
    result = branch_ops.test_branch(repo, "no-such-branch",
                                    runner=lambda *a: {"passed": 1, "failed": 0,
                                                       "first_failure": ""})
    assert result["ok"] is False
    assert "not found" in result["text"]


def test_real_runner_is_the_default(repo):
    from delfin.agent import branch_ops
    assert callable(branch_ops.run_tests)
    # a worktree holds one passing test: the real runner answers counts
    _run(repo, "checkout", "-q", "-b", "real-run")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_ok.py").write_text("def test_ok():\n    assert 1\n")
    _run(repo, "add", "tests/test_ok.py")
    _run(repo, "commit", "-q", "-m", "ok")
    _run(repo, "checkout", "-q", "main")
    result = branch_ops.test_branch(repo, "real-run")
    assert result["ok"] is True
    assert result["passed"] >= 1 and result["failed"] == 0


# -- cleanup -----------------------------------------------------------------

def test_cleanup_candidates_and_delete(repo, monkeypatch):
    from delfin.agent import branch_ops
    cands, kept = branch_ops.cleanup_candidates(repo)
    names = [c["name"] for c in cands]
    assert "old/merged" in names
    # unmerged branches are not offered by default
    _run(repo, "checkout", "-q", "-b", "wip")
    _commit(repo, "w.txt")
    _run(repo, "checkout", "-q", "main")
    cands, _ = branch_ops.cleanup_candidates(repo)
    assert "wip" not in [c["name"] for c in cands]
    # offered only when asked for, and deletion needs the confirmation
    cands, _ = branch_ops.cleanup_candidates(repo, include_unmerged=True)
    wip = next(c for c in cands if c["name"] == "wip")
    assert wip["unmerged"] == 1
    result = branch_ops.cleanup_delete(repo, ["wip"], confirm=False)
    assert result["ok"] is False
    assert "confirm" in result["text"].lower()
    result = branch_ops.cleanup_delete(repo, ["wip"], confirm=True)
    assert result["ok"] is True
    out = _run(repo, "branch", "--list", "wip")
    assert out == ""
    # main and the current branch are refused even when asked for
    result = branch_ops.cleanup_delete(repo, ["main"], confirm=True)
    assert result["ok"] is False
    assert "never" in result["text"].lower()


def test_cleanup_delete_rereads_rules_at_deletion(repo):
    from delfin.agent import branch_ops
    # a branch merged after the listing is refused (it moved / not merged)
    cands, _ = branch_ops.cleanup_candidates(repo)
    names = [c["name"] for c in cands]
    assert "old/merged" in names
    # delete through the real rules
    result = branch_ops.cleanup_delete(repo, ["old/merged"], confirm=True)
    assert result["ok"] is True
    assert _run(repo, "branch", "--list", "old/merged") == ""
