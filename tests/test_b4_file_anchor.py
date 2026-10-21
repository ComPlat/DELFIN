"""B4 P2 finding (reviewer nacht-s19): the file-wait root must be the
session's own workspace, not the state dir.

Verified: _resolve_workspace_path anchored at
state_paths.workspace_state_dir("").parent.parent — the STATE root. A
file wait on "watched.txt" watched <home>/.delfin/watched.txt. Target 3
says "a file in the session's own workspace": the session's workspace is
the process cwd (each DELFIN session runs in its own worktree; the
operator rule "relative paths only" is cwd-relative too).
"""

from __future__ import annotations

import os
import time

import pytest

from delfin.agent import state_paths
from delfin.agent.wait_registry import due, register


@pytest.fixture
def cage(tmp_path, monkeypatch):
    with state_paths.RedirectedUserState(tmp_path / "home"):
        os.makedirs(tmp_path / "ws", exist_ok=True)
        os.makedirs(tmp_path / "home", exist_ok=True)
        monkeypatch.chdir(tmp_path / "ws")  # the session's workspace
        yield tmp_path


def test_file_wait_watches_the_workspace_cwd_not_the_state_dir(cage):
    ws = cage / "ws"
    (ws / "watched.txt").write_text("v1", encoding="utf-8")
    register("file", "watched.txt", key="s16")
    assert due(time.time(), [], key="s16") == []
    (ws / "watched.txt").write_text("v2", encoding="utf-8")
    hits = due(time.time(), [], key="s16")
    assert [w["kind"] for w in hits] == ["file"]
