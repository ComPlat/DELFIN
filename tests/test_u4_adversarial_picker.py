"""Adversarial tests for the U4 session-worktree picker.

These defend the OBSERVABLE contract of the phase-2 fix (busy worktrees in
the new-session picker), NOT the shape of any classifier the builder invents.
The release-offer for a worktree whose session is gone must route through the
existing ``release_session_worktree`` / ``_worktree_is_orphaned`` gate, and a
tree whose owning process is alive, which another host wrote, or which a live
session works in is NEVER releasable here.

Universal: a fake project in tmp_path, no host paths, no real git operations.
"""
from __future__ import annotations

import json
import os
import socket
from types import SimpleNamespace

from delfin.dashboard import agent_sessions as asm


def _fake_ctx(project):
    """A minimal dashboard context rooted at ``project`` (a tmp_path)."""
    return SimpleNamespace(repo_dir=project, agent_dir=project,
                           calc_dir=project)


def _sidecar(wt, *, host="somehost", pid=-1, key="k"):
    """Write a session-worktree sidecar for ``wt`` and return ``wt``."""
    side = wt / ".delfin"
    side.mkdir(parents=True, exist_ok=True)
    (side / "session_worktree.json").write_text(json.dumps({
        "path": str(wt), "branch": "session/abc123", "repo_dir": str(wt.parent),
        "base_ref": "deadbeef", "created_at": 0.0, "host": host, "pid": pid,
        "key": key,
    }), encoding="utf-8")
    return wt


def _worktree(project, name="sb"):
    wt = project / name
    wt.mkdir(parents=True, exist_ok=True)
    return wt


# ---------------------------------------------------------------------------
# Picker exclusion: a held worktree must not be offered as a start folder
# ---------------------------------------------------------------------------

def test_held_worktree_is_not_offered_by_the_picker(tmp_path, monkeypatch):
    """RED on current code: a worktree a live session holds is still offered.

    The picker (_folders_like) lists it as a plain subdirectory; nothing
    consults liveness, so a user can type/select a directory that _on_start
    then refuses. After the fix the held tree must be absent from the choices.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt)
    holder = "My Session on kitzbuhel"
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        holder if str(ws) == str(wt) else ""))
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(wt) not in choices, (
        "a worktree a live session holds must not be offered as a start folder")


def test_plain_folder_is_still_offered_when_a_worktree_is_held(tmp_path,
                                                               monkeypatch):
    """GREEN (regression guard): exclusion is specific to held worktrees.

    A plain directory (no sidecar) must remain offered even while a sibling
    worktree is held — an over-broad fix that drops any folder when some
    other session is busy would wrongly delete normal start options.
    """
    project = tmp_path / "proj"
    project.mkdir()
    held = _worktree(project, "held_wt")
    _sidecar(held)
    plain = _worktree(project, "plain_dir")
    holder = "Session A on kitzbuhel"
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        holder if str(ws) == str(held) else ""))
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(plain) in choices, "a plain directory must remain offered"
    assert str(held) not in choices, "the held worktree alone must be dropped"


def test_an_unheld_worktree_is_still_offered(tmp_path, monkeypatch):
    """GREEN (regression guard): only LIVE-held worktrees are excluded.

    A worktree whose session is gone (no holder) must stay offered — its
    parent dir is still a valid place to look, and the fix must not exclude
    every session worktree, only the busy ones.
    """
    project = tmp_path / "proj"
    project.mkdir()
    for name in ("sb_a", "sb_b"):
        _sidecar(_worktree(project, name))
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        "Session X on kitzbuhel" if str(ws).endswith("sb_a") else ""))
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(project / "sb_b") in choices, "a session-gone worktree stays offered"
    assert str(project / "sb_a") not in choices, "only the held one is dropped"


def test_alive_pid_worktree_is_not_offered_even_with_stale_presence(
        tmp_path, monkeypatch):
    """RED: a worktree whose owning process is alive must not be offered.

    A fix that keys the picker exclusion on session-presence alone
    (``_live_session_in``) misses a tree whose session vanished from the
    presence file but whose owning pid is still alive: stale presence would
    let the busy tree be offered as a start folder. The exclusion must
    consult the worktree sidecar's owning pid, not presence only.
    Uses os.getpid() so the check is deterministic on this machine.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, pid=os.getpid(), host=socket.gethostname() or "somehost")
    # Session gone from the presence file: liveness-by-presence says "".
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(wt) not in choices, (
        "a worktree an alive owning process holds must not be offered")


# ---------------------------------------------------------------------------
# Release safety: never release an alive / other-host / live-session tree
# ---------------------------------------------------------------------------

def test_live_pid_worktree_is_not_orphaned(tmp_path):
    """A tree whose owning process is still running is never releasable.

    ``_worktree_is_orphaned`` must answer with a reason (not ""), so the
    phase-2 release offer never calls release_session_worktree on it.
    Uses os.getpid() so the check is deterministic on this machine.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, pid=os.getpid(), host=socket.gethostname() or "somehost")
    reason = asm._worktree_is_orphaned(asm.read_worktree_sidecar(wt))
    assert reason, "a live-pid worktree must not be reported orphaned/releasable"


def test_other_host_worktree_is_not_orphaned(tmp_path):
    """A tree another host wrote cannot be judged or released from here."""
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, host="other-node.example")
    reason = asm._worktree_is_orphaned(asm.read_worktree_sidecar(wt))
    assert reason, "an other-host worktree must not be releasable from here"


def test_held_worktree_is_not_orphaned(tmp_path, monkeypatch):
    """A tree a live session works in is not releasable, even pid==gone."""
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, pid=-1)
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        "Session Y on kitzbuhel" if str(ws) == str(wt) else ""))
    reason = asm._worktree_is_orphaned(asm.read_worktree_sidecar(wt))
    assert reason, "a live-held worktree must not be reported releasable"


def test_alive_pid_worktree_is_not_releasable_via_state(tmp_path, monkeypatch):
    """RED: the release-offer classifier must not mark a live-pid tree releasable.

    ``_worktree_is_orphaned`` blocks a live owning pid, but the phase-2
    release-offer classifier ``session_worktree_state`` decides ``releasable``
    from session-presence alone (``_live_session_in``). With stale presence
    (session gone from the presence file) but an alive owning pid, it reports
    ``releasable=True`` -- inviting the user to release a worktree a running
    DELFIN process still owns. ``releasable`` must also be False while the
    owning pid is alive (the same liveness the offer side and
    ``_worktree_is_orphaned`` use). Uses os.getpid() so it is deterministic.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, pid=os.getpid(), host=socket.gethostname() or "somehost")
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")
    state = asm.session_worktree_state(str(wt))
    assert state["releasable"] is False, (
        "a worktree an alive owning process owns must not be offered for release")


def test_saved_session_worktree_is_not_releasable_via_state(tmp_path, monkeypatch):
    """RED: the release-offer must not mark a saved-session tree releasable.

    ``_worktree_is_orphaned`` protects a worktree a resumable (saved) session
    would reopen ("a saved session would be reopened into it"), and the
    reclaim path (reclaim_orphaned_worktrees) passes the real saved
    workspaces. But ``session_worktree_state`` calls it with
    ``saved_workspaces=()``, so its ``releasable`` -- which gates the
    picker's release offer, routed to release_session_worktree ->
    exit_worktree, which does NOT re-check saved sessions -- ignores the
    guard. Offering AND releasing a saved-session tree breaks a future
    resume (the session then starts elsewhere without saying so).
    ``releasable`` must consult the saved sessions, not pass an empty list.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _worktree(project)
    _sidecar(wt, pid=-1, host=socket.gethostname() or "somehost")  # dead pid, same host
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")  # no live holder
    # A resumable session would reopen this very worktree:
    monkeypatch.setattr(asm, "_saved_session_workspaces", lambda: [str(wt)])
    state = asm.session_worktree_state(str(wt))
    assert state["releasable"] is False, (
        "a worktree a saved session would reopen must not be offered for release")
