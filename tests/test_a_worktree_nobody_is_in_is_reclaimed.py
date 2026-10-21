"""A session that crashed left its worktree behind for good.

close_session releases a session's own worktree when it is spare, but a
session that never reaches it -- a killed kernel, a lost node, a browser
window that went away with the process -- leaves the checkout and its
orphan ``session/<hex>`` branch on disk. Nothing collected them, and
before the sidecar there was no record to collect them BY.

This is the pass that does, and because it REMOVES DIRECTORIES most of
what is asserted here is what it must refuse to touch. Five conditions,
all of which must hold:

  * the tree carries a sidecar this host wrote
  * its owning process is gone
  * no live session is working in it
  * no saved session would be reopened into it
  * it has no uncommitted changes and no background jobs inside it

The last two are ``exit_worktree``'s decision, asked by the caller rather
than restated, so "is this tree spare" has one answer and not two.

Universal: a repository made in ``tmp_path`` with ``git init`` and one
commit; liveness supplied by a process this test starts and waits for, so
no pid is assumed; the host read from ``socket.gethostname()`` on both
sides.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.dashboard import agent_sessions as AS


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=str(repo), check=True,
                   capture_output=True)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.invalid")
    _git(root, "config", "user.name", "t")
    (root / "a.txt").write_text("one\n", encoding="utf-8")
    _git(root, "add", "a.txt")
    _git(root, "commit", "-qm", "first")
    return root


@pytest.fixture(autouse=True)
def _no_live_sessions(monkeypatch):
    """No presence records unless a test adds them."""
    from delfin.agent import session_presence as P
    monkeypatch.setattr(P, "open_sessions", lambda **kw: [])
    monkeypatch.setattr(AS, "_saved_session_workspaces", lambda: [])


def _a_dead_pid() -> int:
    p = subprocess.Popen([sys.executable, "-c", "pass"])
    p.wait()
    return p.pid


def _make_tree(repo, *, pid=None, host=None):
    """A session worktree whose sidecar names ``pid`` on ``host``."""
    ws = Path(AS.session_worktree(str(repo)))
    side = json.loads((ws / AS._WORKTREE_SIDECAR).read_text())
    if pid is not None:
        side["pid"] = pid
    if host is not None:
        side["host"] = host
    (ws / AS._WORKTREE_SIDECAR).write_text(json.dumps(side))
    return ws


def _branches(repo) -> set[str]:
    out = subprocess.run(["git", "branch", "--format=%(refname:short)"],
                         cwd=str(repo), capture_output=True, text=True)
    return {b.strip() for b in out.stdout.splitlines() if b.strip()}


# ---------------------------------------------------------------------------
# What it reclaims
# ---------------------------------------------------------------------------

def test_a_tree_whose_owner_is_gone_is_reclaimed(repo):
    ws = _make_tree(repo, pid=_a_dead_pid())
    branch = AS.read_worktree_sidecar(ws)["branch"]

    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == [str(ws)], out
    assert not ws.exists()
    assert branch not in _branches(repo), "the orphan branch was left"


def test_the_report_says_why_a_tree_was_kept(repo):
    ws = _make_tree(repo, pid=os.getpid())
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert "still running" in out["kept"][str(ws)]


def test_nothing_to_do_is_an_empty_report(repo):
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out == {"released": [], "kept": {}}


def test_a_directory_with_no_sidecar_is_not_its_business(repo):
    """It touches only trees whose sidecar says DELFIN made them for a
    session on this host. A folder someone put there by hand is not a
    candidate at all."""
    stray = repo / ".delfin" / "worktrees" / "by-hand"
    stray.mkdir(parents=True)
    (stray / "keep.txt").write_text("mine\n", encoding="utf-8")
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out == {"released": [], "kept": {}}
    assert (stray / "keep.txt").is_file()


# ---------------------------------------------------------------------------
# What it must never touch
# ---------------------------------------------------------------------------

def test_a_tree_whose_owner_is_alive_is_kept(repo):
    ws = _make_tree(repo, pid=os.getpid())
    AS.reclaim_orphaned_worktrees(str(repo))
    assert ws.is_dir()


def test_a_tree_belonging_to_another_host_is_kept(repo):
    """A pid written on another login node names nothing here, so it
    cannot be checked. Each host reclaims its own."""
    ws = _make_tree(repo, pid=_a_dead_pid(), host="some-other-node")
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert ws.is_dir()


def test_a_tree_a_live_session_works_in_is_kept(repo, monkeypatch):
    from delfin.agent import session_presence as P

    ws = _make_tree(repo, pid=_a_dead_pid())
    monkeypatch.setattr(P, "open_sessions", lambda **kw: [
        {"key": "k", "title": "Session A", "workspace": str(ws)}])
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert "live session" in out["kept"][str(ws)]
    assert ws.is_dir()


def test_a_tree_a_saved_session_would_reopen_into_is_kept(repo, monkeypatch):
    """A reopened session reuses its tree; removing it would turn a resume
    into a session that starts somewhere else without saying so."""
    ws = _make_tree(repo, pid=_a_dead_pid())
    monkeypatch.setattr(AS, "_saved_session_workspaces", lambda: [str(ws)])
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert "saved session" in out["kept"][str(ws)]
    assert ws.is_dir()


def test_a_saved_workspace_inside_the_tree_also_keeps_it(repo, monkeypatch):
    """The session's workspace is the path at the same place INSIDE the
    tree, not the tree root."""
    ws = _make_tree(repo, pid=_a_dead_pid())
    inner = ws / "sub" / "dir"
    inner.mkdir(parents=True)
    monkeypatch.setattr(AS, "_saved_session_workspaces", lambda: [str(inner)])
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert ws.is_dir()


def test_uncommitted_work_is_never_reclaimed(repo):
    """exit_worktree's decision, and this pass is covered by it."""
    ws = _make_tree(repo, pid=_a_dead_pid())
    (ws / "unsaved.txt").write_text("do not lose me\n", encoding="utf-8")
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert "uncommitted" in out["kept"][str(ws)]
    assert (ws / "unsaved.txt").read_text() == "do not lose me\n"


def test_a_tree_a_job_is_running_in_is_kept(repo, monkeypatch):
    from delfin.agent import worktree as WT

    ws = _make_tree(repo, pid=_a_dead_pid())
    monkeypatch.setattr(WT, "jobs_holding_worktree",
                        lambda path: [{"id": "job-1"}])
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out["released"] == []
    assert ws.is_dir()


# ---------------------------------------------------------------------------
# Bounded, defensive, and visible
# ---------------------------------------------------------------------------

def test_one_pass_removes_at_most_the_limit(repo):
    trees = [_make_tree(repo, pid=_a_dead_pid()) for _ in range(4)]
    out = AS.reclaim_orphaned_worktrees(str(repo), limit=2)
    assert len(out["released"]) == 2
    assert sum(1 for t in trees if t.is_dir()) == 2
    # and the next pass picks up the rest
    out2 = AS.reclaim_orphaned_worktrees(str(repo), limit=2)
    assert len(out2["released"]) == 2
    assert all(not t.is_dir() for t in trees)


def test_a_path_that_is_not_a_repository_is_a_no_op(tmp_path):
    assert AS.reclaim_orphaned_worktrees(str(tmp_path)) == {"released": [],
                                                            "kept": {}}
    assert AS.reclaim_orphaned_worktrees("") == {"released": [], "kept": {}}


def test_an_unreadable_sidecar_is_skipped_not_raised(repo):
    ws = _make_tree(repo, pid=_a_dead_pid())
    (ws / AS._WORKTREE_SIDECAR).write_text("{not json")
    out = AS.reclaim_orphaned_worktrees(str(repo))
    assert out == {"released": [], "kept": {}}
    assert ws.is_dir()


def test_the_sweep_runs_after_the_sessions_are_restored():
    """Order matters: a tree a restored session works in must read as live
    when the question is asked."""
    import inspect

    src = inspect.getsource(AS)
    i_restore = src.index('view["restoring"] = False')
    i_sweep = src.index("_swept = reclaim_orphaned_worktrees")
    assert i_restore < i_sweep


def test_the_reclamation_is_announced():
    import inspect

    src = inspect.getsource(AS)
    arm = src[src.index("_swept = reclaim_orphaned_worktrees"):]
    assert "_say(" in arm[:600], "a removed checkout must not vanish quietly"


def test_the_host_is_read_the_same_way_on_both_sides(repo):
    """The sidecar writer and the sweep must agree on what this host is
    called, or every tree reads as foreign and nothing is ever reclaimed."""
    ws = _make_tree(repo)
    assert (AS.read_worktree_sidecar(ws)["host"]
            == (socket.gethostname() or "").strip())
