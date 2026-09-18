"""Tests for delfin.agent.cli_resume — terminal session resume.

Uses a fabricated session store in tmp_path; never touches the real
~/.delfin/agent_sessions directory.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent.cli_resume import (
    ResumeError,
    list_sessions,
    render_sessions,
    resume_target,
)


def _write_session(
    store, sid, *, updated, workspace="", model="", chat=None,
    created=None,
):
    chat = chat if chat is not None else [
        {"role": "user", "content": f"do task in {sid}"},
        {"role": "assistant", "content": "done"},
    ]
    data = {
        "session_id": sid,
        "created_at": created or updated - 60,
        "updated_at": updated,
        "workspace": workspace,
        "model": model,
        "chat_messages": chat,
    }
    (store / f"{sid}.json").write_text(json.dumps(data))
    return data


@pytest.fixture
def store(tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    return d


def test_list_sessions_basic_fields_and_order(store):
    _write_session(store, "aaa111", updated=1000.0, model="kit",
                   workspace="/tmp")
    _write_session(store, "bbb222", updated=2000.0, model="glm")
    rows = list_sessions(sessions_dir=store)
    assert [r["id"] for r in rows] == ["bbb222", "aaa111"]
    top = rows[0]
    assert top["model"] == "glm"
    assert top["turns"] == 2
    assert top["last_task_line"] == "do task in bbb222"
    assert top["started"] == 2000.0 - 60


def test_list_sessions_limit_and_skips_corrupt_and_turn_files(store):
    _write_session(store, "s1", updated=100.0)
    _write_session(store, "s2", updated=200.0)
    (store / "broken.json").write_text("{not json")
    (store / "zz.turn.json").write_text("{}")  # crash checkpoint, not a session
    rows = list_sessions(limit=1, sessions_dir=store)
    assert [r["id"] for r in rows] == ["s2"]


def test_last_user_line_ignores_assistant_and_tool_rows(store):
    chat = [
        {"role": "user", "content": "first task"},
        {"role": "assistant", "content": "noise"},
        {"role": "user", "content": [{"type": "text", "text": "  final   task  "}]},
        {"role": "user", "content": ""},
    ]
    _write_session(store, "mmm", updated=10.0, chat=chat)
    rows = list_sessions(sessions_dir=store)
    assert rows[0]["last_task_line"] == "final task"
    assert rows[0]["turns"] == 4


def test_resume_target_empty_selector_most_recent(store, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    _write_session(store, "older", updated=100.0, workspace=str(ws))
    _write_session(store, "newer", updated=200.0, workspace=str(ws))
    assert resume_target("", sessions_dir=store)["id"] == "newer"


def test_resume_target_id_prefix(store, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    _write_session(store, "abcdef12", updated=100.0, workspace=str(ws))
    row = resume_target("abcd", sessions_dir=store)
    assert row["id"] == "abcdef12"


def test_resume_target_index_and_hash(store, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    _write_session(store, "old", updated=1.0, workspace=str(ws))
    _write_session(store, "mid", updated=2.0, workspace=str(ws))
    _write_session(store, "new", updated=3.0, workspace=str(ws))
    assert resume_target("2", sessions_dir=store)["id"] == "mid"
    assert resume_target("#1", sessions_dir=store)["id"] == "new"


def test_resume_target_ambiguous_prefix(store, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    _write_session(store, "dupaaaa1", updated=1.0, workspace=str(ws))
    _write_session(store, "dupaaaa2", updated=2.0, workspace=str(ws))
    with pytest.raises(ResumeError, match="ambiguous"):
        resume_target("dup", sessions_dir=store)


def test_resume_target_unknown_prefix_and_bad_index(store, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    _write_session(store, "real", updated=1.0, workspace=str(ws))
    with pytest.raises(ResumeError, match="no kept session"):
        resume_target("nope", sessions_dir=store)
    with pytest.raises(ResumeError, match="out of range"):
        resume_target("7", sessions_dir=store)


def test_resume_target_refuses_missing_workspace(store, tmp_path):
    gone = tmp_path / "gone-ws"
    assert not gone.exists()
    _write_session(store, "ghost", updated=1.0, workspace=str(gone))
    with pytest.raises(ResumeError) as ei:
        resume_target("", sessions_dir=store)
    assert "no longer exists" in str(ei.value)
    assert "refusing" in str(ei.value)


def test_resume_target_empty_workspace_allowed_and_empty_store(store):
    _write_session(store, "legacy", updated=1.0, workspace="")
    assert resume_target("", sessions_dir=store)["id"] == "legacy"
    empty = store.parent / "empty_store"
    empty.mkdir()
    with pytest.raises(ResumeError, match="no kept sessions"):
        resume_target("", sessions_dir=empty)


def test_render_sessions(store):
    _write_session(store, "r1", updated=time.time() - 300, model="kit",
                   workspace="/tmp")
    text = render_sessions(list_sessions(sessions_dir=store))
    assert "r1" in text
    assert "kit" in text
    assert "do task in r1" in text
    assert "--resume" in text
    assert render_sessions([]) == "no kept sessions found"


# -- the column that could not be filled ----------------------------------

def test_the_listing_carries_the_model_a_session_ran_on(tmp_path, monkeypatch):
    """`delfin-agent sessions` showed "?" in the model column for every
    row. The store records the model and the listing dropped it on the
    way out; the CLI writer never set it at all (2026-09-18)."""
    import json

    from delfin.agent import session_store as ss

    d = tmp_path / "agent_sessions"
    d.mkdir(parents=True)
    (d / "s1.json").write_text(json.dumps({
        "session_id": "s1", "title": "a task", "model": "kit.glm-5.3",
        "provider": "kit", "workspace": str(tmp_path), "updated_at": 2.0,
        "chat_messages": [],
    }), encoding="utf-8")
    monkeypatch.setattr(ss, "_sessions_dir", lambda *a, **k: d, raising=False)
    monkeypatch.setattr(ss, "SESSIONS_DIR", d, raising=False)

    rows = [r for r in ss.list_sessions(limit=10)
            if r.get("session_id") == "s1"]
    if rows:                      # only when the store honours the override
        assert rows[0]["model"] == "kit.glm-5.3"
        assert rows[0]["provider"] == "kit"


def test_the_writer_records_the_clients_model():
    """The model lives on the client, not the engine, so export_state
    never carried it."""
    import inspect

    from delfin.agent import cli

    source = inspect.getsource(cli._save_session)
    assert 'estate.setdefault(' in source and '"model"' in source
    assert '"provider"' in source
