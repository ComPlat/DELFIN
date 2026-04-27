"""Tests for delfin.agent.session_store."""

import json
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent.session_store import (
    delete_session,
    fork_session,
    list_sessions,
    load_session,
    save_session,
)


@pytest.fixture
def sessions_dir(tmp_path):
    """Redirect session storage to a temp directory."""
    d = tmp_path / "agent_sessions"
    d.mkdir()
    with patch("delfin.agent.session_store._SESSIONS_DIR", d):
        yield d


def test_save_and_load_roundtrip(sessions_dir):
    chat = [{"role": "user", "content": "Hello"}, {"role": "assistant", "content": "Hi!"}]
    path = save_session(
        "sess-001",
        mode="quick",
        role_index=1,
        route=["session_manager", "builder_agent"],
        chat_messages=chat,
        engine_messages=[{"role": "user", "content": "Hello"}, {"role": "assistant", "content": "Hi!"}],
        cost_usd=0.05,
    )
    assert path.exists()
    assert path.name == "sess-001.json"

    data = load_session("sess-001")
    assert data is not None
    assert data["session_id"] == "sess-001"
    assert data["mode"] == "quick"
    assert data["role_index"] == 1
    assert data["route"] == ["session_manager", "builder_agent"]
    assert len(data["chat_messages"]) == 2
    assert data["cost_usd"] == 0.05
    assert data["title"] == "Hello"  # auto-generated from first user message


def test_auto_title_from_long_message(sessions_dir):
    long_msg = "A" * 200
    save_session("sess-002", chat_messages=[{"role": "user", "content": long_msg}])
    data = load_session("sess-002")
    assert len(data["title"]) <= 84  # 80 chars + "..."


def test_auto_title_preserves_explicit(sessions_dir):
    save_session("sess-003", title="My custom title", chat_messages=[{"role": "user", "content": "ignored"}])
    data = load_session("sess-003")
    assert data["title"] == "My custom title"


def test_save_preserves_created_at(sessions_dir):
    save_session("sess-004")
    data1 = load_session("sess-004")
    created = data1["created_at"]

    time.sleep(0.05)
    save_session("sess-004", cost_usd=0.1)
    data2 = load_session("sess-004")
    assert data2["created_at"] == created
    assert data2["updated_at"] > created


def test_load_nonexistent_returns_none(sessions_dir):
    assert load_session("does-not-exist") is None


def test_list_sessions_newest_first(sessions_dir):
    save_session("sess-a", chat_messages=[{"role": "user", "content": "First"}])
    time.sleep(0.05)
    save_session("sess-b", chat_messages=[{"role": "user", "content": "Second"}])
    time.sleep(0.05)
    save_session("sess-c", chat_messages=[{"role": "user", "content": "Third"}])

    sessions = list_sessions()
    assert len(sessions) == 3
    assert sessions[0]["session_id"] == "sess-c"
    assert sessions[1]["session_id"] == "sess-b"
    assert sessions[2]["session_id"] == "sess-a"


def test_list_sessions_respects_limit(sessions_dir):
    for i in range(10):
        save_session(f"sess-{i:03d}")
        time.sleep(0.01)

    sessions = list_sessions(limit=3)
    assert len(sessions) == 3


def test_list_sessions_summary_fields(sessions_dir):
    save_session(
        "sess-x",
        mode="reviewed",
        chat_messages=[
            {"role": "user", "content": "Fix bug"},
            {"role": "assistant", "content": "Done"},
            {"role": "user", "content": "Test it"},
        ],
    )
    sessions = list_sessions()
    assert len(sessions) == 1
    s = sessions[0]
    assert s["session_id"] == "sess-x"
    assert s["mode"] == "reviewed"
    assert s["message_count"] == 3
    # Full chat_messages should NOT be in the summary
    assert "chat_messages" not in s


def test_delete_session(sessions_dir):
    save_session("sess-del")
    assert load_session("sess-del") is not None
    assert delete_session("sess-del") is True
    assert load_session("sess-del") is None


def test_delete_nonexistent_returns_false(sessions_dir):
    assert delete_session("nope") is False


def test_save_empty_session_id_is_noop(sessions_dir):
    path = save_session("")
    assert list_sessions() == []


def test_list_sessions_empty_dir(sessions_dir):
    assert list_sessions() == []


def test_corrupt_file_skipped(sessions_dir):
    (sessions_dir / "bad.json").write_text("not json at all")
    save_session("good", chat_messages=[{"role": "user", "content": "ok"}])
    sessions = list_sessions()
    assert len(sessions) == 1
    assert sessions[0]["session_id"] == "good"


def test_load_corrupt_returns_none(sessions_dir):
    (sessions_dir / "broken.json").write_text("{invalid")
    assert load_session("broken") is None


# ---------------------------------------------------------------------------
# fork_session
# ---------------------------------------------------------------------------

def test_fork_creates_new_id(sessions_dir):
    save_session("src-1",
                 chat_messages=[{"role": "user", "content": "hello"}],
                 cost_usd=0.42, title="Original")
    new_id = fork_session("src-1")
    assert new_id is not None
    assert new_id != "src-1"
    # Both files exist
    assert (sessions_dir / "src-1.json").exists()
    assert (sessions_dir / f"{new_id}.json").exists()


def test_fork_preserves_chat_and_cost(sessions_dir):
    save_session("src-2",
                 chat_messages=[{"role": "user", "content": "hi"},
                                {"role": "assistant", "content": "hey"}],
                 cost_usd=1.23, token_usage={"input": 100, "output": 50})
    new_id = fork_session("src-2")
    forked = load_session(new_id)
    assert forked["chat_messages"] == [
        {"role": "user", "content": "hi"},
        {"role": "assistant", "content": "hey"},
    ]
    assert forked["cost_usd"] == 1.23
    assert forked["token_usage"] == {"input": 100, "output": 50}


def test_fork_records_forked_from_link(sessions_dir):
    save_session("origin",
                 chat_messages=[{"role": "user", "content": "x"}])
    new_id = fork_session("origin")
    forked = load_session(new_id)
    assert forked["forked_from"] == "origin"
    # session_id is the new one, not the source
    assert forked["session_id"] == new_id


def test_fork_appends_title_suffix(sessions_dir):
    save_session("titled",
                 chat_messages=[{"role": "user", "content": "x"}],
                 title="My Original Work")
    new_id = fork_session("titled")
    forked = load_session(new_id)
    assert forked["title"] == "My Original Work (fork)"


def test_fork_does_not_double_suffix(sessions_dir):
    """Forking a fork must not append the suffix twice."""
    save_session("src",
                 chat_messages=[{"role": "user", "content": "x"}],
                 title="Work (fork)")
    new_id = fork_session("src")
    forked = load_session(new_id)
    # Title still contains exactly one "(fork)"
    assert forked["title"] == "Work (fork)"


def test_fork_unknown_source_returns_none(sessions_dir):
    assert fork_session("does-not-exist") is None


def test_fork_with_explicit_new_id(sessions_dir):
    save_session("origin",
                 chat_messages=[{"role": "user", "content": "x"}])
    new_id = fork_session("origin", new_id="my-fork")
    assert new_id == "my-fork"
    assert (sessions_dir / "my-fork.json").exists()


def test_fork_resets_timestamps(sessions_dir):
    save_session("origin",
                 chat_messages=[{"role": "user", "content": "x"}])
    # Manually back-date the source's updated_at
    src_data = json.loads((sessions_dir / "origin.json").read_text())
    src_data["updated_at"] = src_data["created_at"] = 1000.0
    (sessions_dir / "origin.json").write_text(json.dumps(src_data, indent=1))

    new_id = fork_session("origin")
    forked = load_session(new_id)
    # Fork's timestamps are recent, not 1000.0
    assert forked["created_at"] > 1000.0
    assert forked["updated_at"] > 1000.0


def test_fork_corrupt_source_returns_none(sessions_dir):
    (sessions_dir / "corrupt.json").write_text("{not valid")
    assert fork_session("corrupt") is None
