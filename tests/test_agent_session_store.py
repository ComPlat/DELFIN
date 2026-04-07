"""Tests for delfin.agent.session_store."""

import json
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent.session_store import (
    delete_session,
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
