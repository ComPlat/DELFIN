"""Waiting on something is a status, and completed obeys the DAG.

THE INCIDENT. The valid statuses were pending / in_progress / completed
/ deleted. A task waiting on a user answer, a missing credential or a
failed dependency had nowhere to live: it stayed "pending", identical to
one nobody had started, in the list the user reads and in the reminder
the model reads.

At the same time the dependency guard applied only to the ``in_progress``
transition, so ``update(child, status="completed")`` walked straight
around the DAG — no in_progress step required — and the edge that exists
to order the two tasks ordered nothing.

Two invariants that were prompt-only and are one-line data checks are
enforced here too: never two tasks in_progress at once, and no silent
pending → completed (the step is what makes the work window exist).
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.agent_tasks import (
    OPEN_STATUSES, TaskStore, get_store, open_task_summary,
)
from delfin.agent.api_client import KitToolPermissions, _doc_executor
from delfin.agent import task_ticker as TT


@pytest.fixture
def store(tmp_path) -> TaskStore:
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return TaskStore(tmp_path)


# ---------------------------------------------------------------------------
# blocked
# ---------------------------------------------------------------------------

def test_blocked_needs_a_reason(store):
    t = store.create("Deploy", session_id="s")
    with pytest.raises(ValueError) as exc:
        store.update(t["id"], status="blocked")
    assert "blocked_reason" in str(exc.value)


def test_blocked_keeps_what_it_waits_on(store):
    t = store.create("Deploy", session_id="s")
    out = store.update(t["id"], status="blocked",
                       blocked_reason="waiting for the KIT API key")
    assert out["status"] == "blocked"
    assert "KIT API key" in out["blocked_reason"]


def test_blocked_is_open_work(store):
    t = store.create("Deploy", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="user answer")
    assert "blocked" in OPEN_STATUSES
    summary = open_task_summary(store.base_dir, "s")
    assert summary["state"] == "open"
    assert summary["counts"]["blocked"] == 1


def test_a_blocked_task_can_be_started_once_it_is_unblocked(store):
    t = store.create("Deploy", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="key")
    assert store.update(t["id"], status="in_progress")["status"] == "in_progress"


def test_blocked_has_its_own_bucket_in_the_ticker(store):
    t = store.create("Deploy to the cluster", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="missing key")
    html = TT.render_html(store.base_dir, session_id="s")
    assert "Deploy to the cluster" in html
    assert "missing key" in html
    text = TT.render_text(store.base_dir, session_id="s")
    assert "[!]" in text


def test_blocked_has_its_own_bucket_in_the_per_turn_block(store, monkeypatch):
    from delfin.agent.engine import AgentEngine

    class _Perms:
        workspace = store.base_dir
        task_session_id = "s"
        mode = "default"

    t = store.create("Deploy to the cluster", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="missing key")
    eng = AgentEngine.__new__(AgentEngine)
    monkeypatch.setattr(AgentEngine, "kit_permissions",
                        property(lambda self: _Perms()))
    block = eng._build_open_tasks_block()
    assert "[blocked]" in block
    assert "missing key" in block


def test_the_tool_reports_blocked_as_its_own_group(tmp_path):
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default",
                               task_session_id="s")
    tid = json.loads(_doc_executor.execute(
        "task_create", {"subject": "Deploy"}, perms))["task"]["id"]
    out = json.loads(_doc_executor.execute("task_update", {
        "task_id": tid, "status": "blocked",
        "blocked_reason": "the KIT API key is missing",
    }, perms))
    assert out["task"]["status"] == "blocked"
    listed = json.loads(_doc_executor.execute("task_list", {}, perms))
    assert listed["by_status"]["blocked"] == 1


def test_the_tool_refuses_blocked_without_a_reason(tmp_path):
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    tid = json.loads(_doc_executor.execute(
        "task_create", {"subject": "Deploy"}, perms))["task"]["id"]
    out = json.loads(_doc_executor.execute(
        "task_update", {"task_id": tid, "status": "blocked"}, perms))
    assert "blocked_reason" in out["error"]


# ---------------------------------------------------------------------------
# completed obeys the DAG
# ---------------------------------------------------------------------------

def test_completed_is_refused_while_a_predecessor_is_open(store):
    parent = store.create("setup", session_id="s")
    child = store.create("execute", session_id="s",
                         blocked_by=[parent["id"]])
    store.update(child["id"], status="pending")
    with pytest.raises(ValueError) as exc:
        store.update(child["id"], status="completed")
    assert "blocked by unfinished" in str(exc.value)


def test_completed_is_allowed_once_the_predecessor_is_done(store):
    parent = store.create("setup", session_id="s")
    child = store.create("execute", session_id="s", blocked_by=[parent["id"]])
    store.update(parent["id"], status="in_progress")
    store.update(parent["id"], status="completed")
    store.update(child["id"], status="in_progress")
    assert store.update(child["id"], status="completed")["status"] == "completed"


# ---------------------------------------------------------------------------
# The two invariants that were prompt-only
# ---------------------------------------------------------------------------

def test_pending_cannot_jump_straight_to_completed(store):
    t = store.create("Wire the parser", session_id="s")
    with pytest.raises(ValueError) as exc:
        store.update(t["id"], status="completed")
    assert "in_progress" in str(exc.value)


def test_re_completing_a_completed_task_is_not_an_error(store):
    """Idempotence matters: a model that repeats the call must not be
    handed an error it will then try to 'fix'."""
    t = store.create("Wire the parser", session_id="s")
    store.update(t["id"], status="in_progress")
    store.update(t["id"], status="completed")
    assert store.update(t["id"], status="completed")["status"] == "completed"


def test_only_one_task_is_in_progress_per_session(store):
    a = store.create("first", session_id="s")
    b = store.create("second", session_id="s")
    store.update(a["id"], status="in_progress")
    with pytest.raises(ValueError) as exc:
        store.update(b["id"], status="in_progress")
    assert "already in_progress" in str(exc.value)
    # And the error says both ways out.
    assert "completed" in str(exc.value) and "blocked" in str(exc.value)


def test_parking_the_first_task_frees_the_slot(store):
    a = store.create("first", session_id="s")
    b = store.create("second", session_id="s")
    store.update(a["id"], status="in_progress")
    store.update(a["id"], status="blocked", blocked_reason="user answer")
    assert store.update(b["id"], status="in_progress")["status"] == "in_progress"


def test_two_sessions_may_each_have_one_in_progress(store):
    a = store.create("mine", session_id="s1")
    b = store.create("theirs", session_id="s2")
    store.update(a["id"], status="in_progress")
    assert store.update(b["id"], status="in_progress")["status"] == "in_progress"


def test_starting_a_task_records_the_window_start(store):
    t = store.create("Wire the parser", session_id="s")
    assert store.get(t["id"])["started_at"] == ""
    started = store.update(t["id"], status="in_progress")
    assert started["started_at"]
