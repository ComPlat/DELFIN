"""Cross-session task visibility: open work from previous sessions.

The task store is workspace-scoped and survives across sessions, but every
listing surface filters to the CURRENT session id — a fresh session mints a
new id, so pending/in_progress tasks from earlier sessions were invisible
and died silently. `open_foreign_tasks` summarises that leftover work and
the engine surfaces it exactly once, on a session's first prompt build,
with explicit opt-in adoption (session scoping itself stays authoritative —
it fixed a task-leak bug and lists are never auto-merged).
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone

from delfin.agent.agent_tasks import TaskStore, get_store, open_foreign_tasks
from delfin.agent.engine import AgentEngine


# ---------------------------------------------------------------------------
# Helper: open_foreign_tasks
# ---------------------------------------------------------------------------

def test_helper_filters_to_open_foreign_only(tmp_path):
    s = get_store(tmp_path)
    s.create("mine open", session_id="cur")               # own → excluded
    s.create("theirs pending", session_id="old1")         # foreign, open
    started = s.create("theirs running", session_id="old2")
    s.update(started["id"], status="in_progress")         # foreign, open
    done = s.create("theirs finished", session_id="old1")
    s.update(done["id"], status="completed")              # closed → excluded
    gone = s.create("theirs dropped", session_id="old1")
    s.update(gone["id"], status="deleted")                # deleted → excluded

    out = open_foreign_tasks(tmp_path, "cur")
    assert out["count"] == 2
    subjects = {t["subject"] for t in out["tasks"]}
    assert subjects == {"theirs pending", "theirs running"}
    # Records keep their owning session id so the agent can tell them apart.
    assert {t["session_id"] for t in out["tasks"]} == {"old1", "old2"}
    assert all(t["status"] in ("pending", "in_progress") for t in out["tasks"])


def test_helper_empty_for_empty_session_id(tmp_path):
    s = get_store(tmp_path)
    s.create("left behind", session_id="old")
    # Empty current id → task_list already shows every workspace task, so
    # nothing is invisible and the summary must stay empty.
    assert open_foreign_tasks(tmp_path, "") == {
        "count": 0, "oldest_age_days": 0, "tasks": [],
    }


def test_helper_caps_list_but_counts_all(tmp_path):
    s = get_store(tmp_path)
    for i in range(7):
        s.create(f"leftover {i}", session_id="old")
    out = open_foreign_tasks(tmp_path, "cur")
    assert out["count"] == 7
    assert len(out["tasks"]) == 5          # default cap
    out2 = open_foreign_tasks(tmp_path, "cur", cap=2)
    assert out2["count"] == 7
    assert len(out2["tasks"]) == 2


def test_helper_age_math_and_oldest_first(tmp_path):
    s = get_store(tmp_path)
    a = s.create("three days old", session_id="old")
    b = s.create("fresh", session_id="old")
    # Rewrite the on-disk timestamps to known ages (the store re-reads the
    # file on every list, so this is authoritative).
    now = datetime.now(timezone.utc)
    stamp = (now - timedelta(days=3, hours=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    data = json.loads(s.path.read_text(encoding="utf-8"))
    for t in data["tasks"]:
        if t["id"] == a["id"]:
            t["updated_at"] = stamp
    s.path.write_text(json.dumps(data), encoding="utf-8")

    out = open_foreign_tasks(tmp_path, "cur")
    assert out["count"] == 2
    assert out["oldest_age_days"] == 3
    # Oldest first; the fresh task ages 0 days.
    assert [t["id"] for t in out["tasks"]] == [a["id"], b["id"]]
    assert [t["age_days"] for t in out["tasks"]] == [3, 0]


def test_helper_never_raises(tmp_path):
    # No store file at all → empty summary, no exception.
    assert open_foreign_tasks(tmp_path / "nowhere", "cur")["count"] == 0
    # Corrupt store file → empty summary, no exception.
    s = TaskStore(tmp_path)
    s.path.parent.mkdir(parents=True, exist_ok=True)
    s.path.write_text("{not json", encoding="utf-8")
    assert open_foreign_tasks(tmp_path, "cur")["count"] == 0


# ---------------------------------------------------------------------------
# Engine block: one-shot notice on the session's first prompt build
# ---------------------------------------------------------------------------

def _bare_engine() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = [{"role": "user", "content": "hello"}]  # first-turn state
    return eng


def _perms_for(tmp_path, sid):
    class _Perms:
        workspace = tmp_path
        task_session_id = sid
    return _Perms()


def test_foreign_block_on_first_build_only(tmp_path, monkeypatch):
    s = get_store(tmp_path)
    t1 = s.create("Finish BoTorch wrapper", session_id="prev")
    s.create("Write comparison notebook", session_id="prev")

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "cur")),
    )
    block = eng._build_open_foreign_tasks_block()
    assert "# Open work from previous sessions" in block
    assert "2 open task(s)" in block
    assert "task_adopt" in block and "task_list(all_sessions=true)" in block
    assert f"id {t1['id']}" in block and "Finish BoTorch wrapper" in block
    # One-shot: the very next build stays silent.
    assert eng._build_open_foreign_tasks_block() == ""


def test_foreign_block_empty_without_foreign_tasks(tmp_path, monkeypatch):
    s = get_store(tmp_path)
    s.create("current work", session_id="cur")   # own session only

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "cur")),
    )
    assert eng._build_open_foreign_tasks_block() == ""


def test_foreign_block_suppressed_for_restored_history(tmp_path, monkeypatch):
    s = get_store(tmp_path)
    s.create("left behind", session_id="prev")

    eng = _bare_engine()
    eng.messages = [
        {"role": "user", "content": "old"},
        {"role": "assistant", "content": "reply"},
        {"role": "user", "content": "new"},
    ]
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "cur")),
    )
    # Mid-conversation (e.g. a restored session) → never shown, now or later.
    assert eng._build_open_foreign_tasks_block() == ""
    eng.messages = [{"role": "user", "content": "new"}]
    assert eng._build_open_foreign_tasks_block() == ""


def test_foreign_block_caps_listing_at_five(tmp_path, monkeypatch):
    s = get_store(tmp_path)
    for i in range(8):
        s.create(f"leftover {i}", session_id="prev")

    eng = _bare_engine()
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: _perms_for(tmp_path, "cur")),
    )
    block = eng._build_open_foreign_tasks_block()
    assert "8 open task(s)" in block
    assert block.count("- id ") == 5


# ---------------------------------------------------------------------------
# Adoption: session_id becomes updatable (the task_adopt path)
# ---------------------------------------------------------------------------

def test_store_update_session_id_adopts_task(tmp_path):
    s = get_store(tmp_path)
    t = s.create("orphaned work", session_id="prev")
    adopted = s.update(t["id"], session_id="cur")
    assert adopted["session_id"] == "cur"
    # Now visible in the current session's list …
    assert [x["id"] for x in s.list(session_id="cur")] == [t["id"]]
    # … and no longer foreign.
    assert open_foreign_tasks(tmp_path, "cur")["count"] == 0
    assert open_foreign_tasks(tmp_path, "prev")["count"] == 1


# ---------------------------------------------------------------------------
# blocked_by round-trip through the store (DAG ordering stays enforced)
# ---------------------------------------------------------------------------

def test_blocked_by_round_trip_via_store(tmp_path):
    import pytest

    s = get_store(tmp_path)
    a = s.create("parent", session_id="x")
    b = s.create("child", session_id="x", blocked_by=[a["id"]])
    assert b["blocked_by"] == [a["id"]]
    assert s.get(a["id"])["blocks"] == [b["id"]]
    # Blocked → cannot start.
    with pytest.raises(ValueError, match="blocked by unfinished"):
        s.update(b["id"], status="in_progress")
    # Edge edits round-trip and keep the reverse index in sync.
    s.update(b["id"], remove_blocked_by=[a["id"]])
    assert s.get(b["id"])["blocked_by"] == []
    assert s.get(a["id"])["blocks"] == []
    s.update(b["id"], add_blocked_by=[a["id"]])
    assert s.get(b["id"])["blocked_by"] == [a["id"]]
    assert s.get(a["id"])["blocks"] == [b["id"]]
    # Completing the parent unblocks the child.
    s.update(a["id"], status="completed")
    assert s.update(b["id"], status="in_progress")["status"] == "in_progress"
