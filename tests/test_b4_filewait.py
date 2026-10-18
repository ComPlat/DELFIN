"""B4 phase 2 red control: a wait on a file wakes on change, not on touch.

Verified gap: wait_registry.due() (delfin/agent/wait_registry.py) matches
only incoming (kind, value) events and time waits; a "file" wait registered
in phase 1 stores a content signature (file_sig, size+sha256, mtime
excluded) but due() never evaluates it, so a file change can never wake a
waiting session. The wake route is due() itself, called at the idle path's
existing cadence (repl.py:4224) — no new process, no extra polling.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

from delfin.agent import state_paths
from delfin.agent.wait_registry import due, register


@pytest.fixture
def home(tmp_path, monkeypatch):
    # The file-wait anchor is the process cwd: every test gets its own
    # throwaway workspace so nothing leaks into the checkout (the
    # conftest guard counts new checkout paths).
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.chdir(ws)
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


def _workspace_file(rel: str, content: str) -> str:
    # The file-wait anchor is the process cwd (the session's own
    # workspace; reviewer nacht-s19's P2 finding moved it off the state
    # dir). Write there; tmp_path fixture is cwd-independent for the
    # redirected HOME only.
    p = Path(rel)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8")
    return rel


def test_file_wait_fires_when_content_changes(home):
    rel = _workspace_file("reports/state.md", "v1")
    register("file", rel, key="s16")
    assert due(time.time(), [], key="s16") == []  # unchanged: no fire
    _workspace_file("reports/state.md", "v2")
    hits = due(time.time(), [], key="s16")
    assert [w["kind"] for w in hits] == ["file"]
    assert hits[0]["value"] == rel
    # consumed: the wake happens once
    assert due(time.time(), [], key="s16") == []


def test_file_wait_does_not_fire_on_touch_without_change(home):
    rel = _workspace_file("reports/state.md", "v1")
    register("file", rel, key="s16")
    assert due(time.time(), [], key="s16") == []
    # touch: same content, new mtime (cwd anchor)
    os.utime(rel, (time.time() + 10, time.time() + 10))
    assert due(time.time(), [], key="s16") == []


def test_file_wait_fires_when_the_file_appears(home):
    # registered against a not-yet-existing file: its appearance IS a
    # change (signature 'missing' -> content signature)
    rel = "reports/later.md"
    register("file", rel, key="s16")
    assert due(time.time(), [], key="s16") == []
    _workspace_file(rel, "now I exist")
    assert [w["kind"] for w in due(time.time(), [], key="s16")] == ["file"]


def test_file_wait_rejects_escaping_paths(home):
    for bad in ("/etc/passwd", "../outside.md", "a/../../outside.md"):
        with pytest.raises(ValueError):
            register("file", bad, key="s16")


def test_file_wait_rejects_symlink_escape(home, tmp_path):
    root = Path.cwd()
    outside = tmp_path / "outside" / "b4_outside_target.txt"
    outside.parent.mkdir(parents=True, exist_ok=True)
    outside.write_text("secret", encoding="utf-8")
    link = root / "reports"
    link.mkdir(parents=True, exist_ok=True)
    (link / "link.md").symlink_to(outside)
    with pytest.raises(ValueError):
        register("file", "reports/link.md", key="s16")
