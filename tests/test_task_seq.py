"""Session-relative task numbering (bug 20260619-172400, ka_xn0397).

The global task `id` keeps climbing across sessions in a long-lived workspace
("bin bei task 90"). `list(with_seq=True)` adds a small 1-based, session-local
`seq` for display, while the global `id` stays the reference for
task_update/task_get. Display surfaces (ticker, reminder, create hint) lead with
`seq`; the id is preserved everywhere it's needed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.agent_tasks import TaskStore
from delfin.agent import task_ticker


def _store(tmp_path) -> TaskStore:
    return TaskStore(tmp_path)


def test_seq_is_session_relative_one_based(tmp_path):
    s = _store(tmp_path)
    # Two sessions interleaved; global ids climb 1..4.
    a1 = s.create("A first", session_id="sessA")
    b1 = s.create("B first", session_id="sessB")
    a2 = s.create("A second", session_id="sessA")
    b2 = s.create("B second", session_id="sessB")
    assert [a1["id"], b1["id"], a2["id"], b2["id"]] == [1, 2, 3, 4]

    a = s.list(session_id="sessA", with_seq=True)
    b = s.list(session_id="sessB", with_seq=True)
    # Each session restarts at 1 for display, regardless of global id.
    assert {t["id"]: t["seq"] for t in a} == {1: 1, 3: 2}
    assert {t["id"]: t["seq"] for t in b} == {2: 1, 4: 2}


def test_seq_stable_when_a_task_completes(tmp_path):
    s = _store(tmp_path)
    s.create("one", session_id="x")
    s.create("two", session_id="x")
    s.create("three", session_id="x")
    s.update(1, status="completed")
    # Completed task keeps its slot; the others don't renumber.
    seqs = {t["id"]: t["seq"] for t in s.list(session_id="x", with_seq=True)}
    assert seqs == {1: 1, 2: 2, 3: 3}


def test_list_without_with_seq_is_unchanged(tmp_path):
    s = _store(tmp_path)
    s.create("one", session_id="x")
    tasks = s.list(session_id="x")
    assert "seq" not in tasks[0]            # opt-in only → no behaviour change


def test_ticker_shows_seq_not_global_id(tmp_path):
    s = _store(tmp_path)
    s.create("old session task", session_id="prev")   # id 1
    s.create("real one", session_id="cur")            # id 2 → seq 1
    html = task_ticker.render_html(tmp_path, session_id="cur")
    assert "#1 real one" in html                       # seq, not "#2"
    txt = task_ticker.render_text(tmp_path, session_id="cur")
    assert "#1 real one" in txt


def test_create_hint_uses_seq_but_keeps_id():
    from delfin.agent.api_client import _doc_executor, KitToolPermissions
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="bypassPermissions")
        # simulate a long-lived workspace: pre-seed prior-session tasks
        store = _doc_executor._task_store(perms)
        for i in range(5):
            store.create(f"old {i}", session_id="old")
        perms.task_session_id = "new"
        out = json.loads(_doc_executor._execute_task_create(
            {"subject": "fresh task"}, perms))
        assert out["status"] == "created"
        # Global id is 6 (6th task in the workspace) but seq is 1 this session.
        assert out["task"]["id"] == 6
        assert out["task"]["seq"] == 1
        assert "task 1" in out["hint"] and "id 6" in out["hint"]
