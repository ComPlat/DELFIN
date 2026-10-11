"""How the new-session picker decides that a session worktree is held.

The owning process recorded in ``.delfin/session_worktree.json`` is judged
by ``proc_identity.alive``: pid, start time and host together. A pid alone
names a different process once the number is handed out again, and a pid
written on another login node (shared home) names nothing on this one.

Every test runs against a presence directory and a session store under
tmp_path, never the user's ``~/.delfin``.
"""
from __future__ import annotations

import json
import os
import socket
from types import SimpleNamespace

import pytest

from delfin.agent import proc_identity
from delfin.agent import session_presence
from delfin.agent import session_store
from delfin.dashboard import agent_sessions as asm


@pytest.fixture(autouse=True)
def _private_stores(tmp_path, monkeypatch):
    presence = tmp_path / "presence"
    presence.mkdir()
    monkeypatch.setattr(session_presence, "_DIR", presence)
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "saved")
    return presence


def _ctx(project):
    return SimpleNamespace(repo_dir=project, agent_dir=project,
                           calc_dir=project)


def _tree(project, name, **fields):
    wt = project / name
    (wt / ".delfin").mkdir(parents=True)
    record = {"path": str(wt), "branch": "session/x", "repo_dir": str(project),
              "base_ref": "x", "created_at": 0.0,
              "host": socket.gethostname(), "pid": -1}
    record.update(fields)
    (wt / ".delfin" / "session_worktree.json").write_text(
        json.dumps(record), encoding="utf-8")
    return wt


def test_a_pid_written_on_another_host_does_not_hide_the_tree(tmp_path):
    project = tmp_path / "proj"
    project.mkdir()
    wt = _tree(project, "elsewhere", host="other-node.invalid",
               pid=os.getpid())
    assert str(wt) in asm.workspace_choices(_ctx(project), str(project) + "/")


def test_a_reused_pid_does_not_hide_the_tree(tmp_path):
    """Same host, the number is alive, but it is a later process."""
    project = tmp_path / "proj"
    project.mkdir()
    wt = _tree(project, "reused", pid=os.getpid(), proc_start="1")
    assert str(wt) in asm.workspace_choices(_ctx(project), str(project) + "/")
    side = asm.read_worktree_sidecar(wt)
    assert asm._worktree_is_orphaned(side) == ""


def test_the_live_owner_still_hides_the_tree(tmp_path):
    project = tmp_path / "proj"
    project.mkdir()
    wt = _tree(project, "mine", pid=os.getpid(),
               proc_start=proc_identity.process_start(os.getpid()))
    assert str(wt) not in asm.workspace_choices(_ctx(project),
                                                str(project) + "/")


def test_the_sidecar_records_the_owner_start_time(tmp_path):
    info = SimpleNamespace(path=tmp_path / "wt", repo_dir=tmp_path,
                           branch="session/x", base_ref="x", created_at=0.0)
    asm._write_worktree_sidecar(info)
    side = asm.read_worktree_sidecar(tmp_path / "wt")
    assert side["proc_start"] == proc_identity.process_start(os.getpid())


def test_the_presence_records_are_read_once_per_listing(tmp_path, monkeypatch):
    """Typing inside a worktree lists its subfolders; each of them finds the
    same sidecar. The presence directory is read once per call, not once
    per entry (measured before: 199 reads for 200 entries, 200 ms)."""
    project = tmp_path / "proj"
    project.mkdir()
    wt = _tree(project, "big")
    for i in range(40):
        (wt / f"d{i:02}").mkdir()
    calls = []
    real = session_presence.open_sessions
    monkeypatch.setattr(session_presence, "open_sessions",
                        lambda **kw: calls.append(1) or real(**kw))
    choices = asm.workspace_choices(_ctx(project), str(wt) + "/")
    assert len([c for c in choices if c.startswith(str(wt) + "/")]) == 40
    assert len(calls) <= 1


def test_a_tree_hidden_by_its_live_owner_says_why(tmp_path):
    """The start form explains a missing entry; the owner-alive case had no
    text at all because only a presence holder produced a hint."""
    project = tmp_path / "proj"
    project.mkdir()
    wt = _tree(project, "mine", pid=os.getpid(),
               proc_start=proc_identity.process_start(os.getpid()))
    state = asm.session_worktree_state(str(wt))
    assert state["holder"] == ""
    assert state["why"]
