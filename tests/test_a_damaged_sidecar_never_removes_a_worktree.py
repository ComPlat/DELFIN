"""A damaged sidecar must never send a worktree to its removal.

``AS.read_worktree_sidecar`` answers None for a sidecar that cannot be
parsed, and ``reclaim_orphaned_worktrees`` only ever releases a tree it
could positively identify as DELFIN's from a sidecar this host wrote.
So the lenient read is correct by design -- but the safety property is
what must be pinned: a worktree whose sidecar is damaged survives the
sweep, with nothing lost and nothing deleted.
"""
import json
import socket
from pathlib import Path

import pytest

from delfin.dashboard import agent_sessions as AS


def _tree(tmp_path):
    tree = tmp_path / "wt1"
    (tree / ".delfin").mkdir(parents=True)
    return tree


def test_a_damaged_sidecar_keeps_the_worktree(tmp_path, monkeypatch):
    tree = _tree(tmp_path)
    (tree / ".delfin" / "session_worktree.json").write_text(
        "{not json", encoding="utf-8")
    assert AS.read_worktree_sidecar(str(tree)) is None

    # Even with no owning process and no live session, the sweep must
    # not recognise a tree it cannot read a sidecar from.
    monkeypatch.setattr(socket, "gethostname", lambda: "h1")
    out = AS.reclaim_orphaned_worktrees(str(tmp_path))
    assert tree.is_dir(), "the tree was removed although its sidecar was unreadable"
    assert out["released"] == []


def test_read_worktree_sidecar_does_not_quarantine_a_damaged_record(tmp_path):
    """The sidecar has no writer that reads first: the bytes are the only
    copy of the record, and no caller rewrites the file, so nothing may
    move them aside either."""
    tree = _tree(tmp_path)
    side = tree / ".delfin" / "session_worktree.json"
    side.write_text("{not json", encoding="utf-8")
    AS.read_worktree_sidecar(str(tree))
    assert side.read_text(encoding="utf-8") == "{not json"
    assert not list((tree / ".delfin").glob("*.unreadable-*"))
