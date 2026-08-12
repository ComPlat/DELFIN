"""The "Open items" section was reading a list nobody filled.

THE INCIDENT. Both consumers of the open-work list — the handoff brief
for a fresh agent, and the saved episode — key on ``todo_payload``. That
field had exactly one producer: the dashboard branch that handles an
external CLI's todo tool. Nothing converted the real task store into it,
the engine's exported state carries no task field, and the headless save
never passed the argument at all.

Consequence on the open-weights path, which is the one in use: every
handoff brief and every episode reported an empty "Open items" section
while the store still held in_progress work — a fresh agent was briefed
that nothing was outstanding, which is the one thing a handoff exists to
prevent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import session_store as ss
from delfin.agent.agent_tasks import get_store
from delfin.agent.episodes import build_episode_from_state


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions")
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return tmp_path


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def test_a_brief_after_one_task_create_names_that_task(fake_home, workspace):
    get_store(workspace).create(
        "Wire the parser into ops_server", session_id="sid-1")
    ss.save_session("sid-1", workspace=str(workspace),
                    chat_messages=[{"role": "user", "content": "do it"}])
    data = ss.load_session("sid-1")
    brief = ss.build_handoff_brief(data)
    assert "Wire the parser into ops_server" in brief


def test_the_saved_session_carries_the_task_list(fake_home, workspace):
    store = get_store(workspace)
    t = store.create("Run the migration", session_id="sid-2")
    store.update(t["id"], status="in_progress")
    ss.save_session("sid-2", workspace=str(workspace))
    payload = ss.load_session("sid-2")["todo_payload"]
    assert [p["subject"] for p in payload] == ["Run the migration"]
    assert payload[0]["status"] == "in_progress"


def test_an_explicit_payload_still_wins(fake_home, workspace):
    """The dashboard's external-CLI branch has a real list of its own; the
    fallback must not overwrite it."""
    get_store(workspace).create("from the store", session_id="sid-3")
    ss.save_session("sid-3", workspace=str(workspace),
                    todo_payload=[{"id": 9, "subject": "from the caller",
                                   "status": "pending"}])
    payload = ss.load_session("sid-3")["todo_payload"]
    assert [p["subject"] for p in payload] == ["from the caller"]


def test_a_brief_built_straight_from_a_workspace_dict_finds_the_tasks(
        fake_home, workspace):
    get_store(workspace).create("Collect the results", session_id="sid-4")
    brief = ss.build_handoff_brief(
        {"session_id": "sid-4", "workspace": str(workspace)})
    assert "Collect the results" in brief


def test_a_blocked_task_reaches_the_brief_with_its_reason(
        fake_home, workspace):
    store = get_store(workspace)
    t = store.create("Deploy to the cluster", session_id="sid-5")
    store.update(t["id"], status="blocked",
                 blocked_reason="the KIT API key is missing")
    ss.save_session("sid-5", workspace=str(workspace))
    brief = ss.build_handoff_brief(ss.load_session("sid-5"))
    assert "Deploy to the cluster" in brief
    assert "KIT API key" in brief


def test_the_episode_open_items_come_from_the_store_too(fake_home, workspace):
    """The headless save path: ``build_episode_from_state`` gets the
    engine's exported state, which has no task field at all."""
    store = get_store(workspace)
    t = store.create("Finish the wrapper", session_id="sid-6")
    store.update(t["id"], status="in_progress")
    fields = build_episode_from_state(
        {"session_id": "sid-6", "project_dir": str(workspace),
         "cost_usd": 0.5},
        [{"role": "user", "content": "go"}],
    )
    assert any("Finish the wrapper" in item for item in fields["open_items"])


def test_a_finished_list_produces_no_open_items(fake_home, workspace):
    store = get_store(workspace)
    t = store.create("Done and dusted", session_id="sid-7")
    store.update(t["id"], status="in_progress")
    store.update(t["id"], status="completed")
    ss.save_session("sid-7", workspace=str(workspace))
    brief = ss.build_handoff_brief(ss.load_session("sid-7"))
    assert "(nothing outstanding" in brief or "Open items" in brief
    assert "▶ #" not in brief


def test_the_fill_never_breaks_a_save(fake_home, tmp_path):
    """Bookkeeping must not be able to fail a save: a workspace that does
    not exist simply yields no items."""
    assert ss.tasks_as_todo_payload(tmp_path / "nope", "sid") == []
    assert ss.tasks_as_todo_payload("", "sid") == []
    ss.save_session("sid-8", workspace=str(tmp_path / "nope"))
    assert ss.load_session("sid-8")["todo_payload"] == []
