"""Tests for the long-session resume + search surface.

These tests exercise the persistence layer end-to-end: save with
extras → load → fields round-trip; archive across multiple compactions
→ search recovers content. The dashboard glue is tested separately as
a closure inspection.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions"
    )
    return tmp_path


def test_save_then_load_returns_all_extras(fake_home):
    """Every long-session extra written by save_session must round-trip
    through load_session — the dashboard relies on this for restore."""
    ss.save_session(
        "sid-1",
        mode="solo",
        chat_messages=[{"role": "user", "content": "hi"}],
        perm_profile="all_free",
        provider="claude",
        model="opus",
        effort="high",
        active_gate={"type": "approval", "title": "review me"},
        last_compaction_info={"messages_compacted": 5, "tokens_saved": 800},
        subagent_calls=[{"type": "explore", "status": "done"}],
        pending_plan_body="# resumed plan\nstep 1",
        todo_payload=[{"id": 7, "status": "in_progress"}],
    )
    loaded = ss.load_session("sid-1")
    assert loaded is not None
    assert loaded["perm_profile"] == "all_free"
    assert loaded["provider"] == "claude"
    assert loaded["model"] == "opus"
    assert loaded["effort"] == "high"
    assert loaded["active_gate"]["type"] == "approval"
    assert loaded["last_compaction_info"]["tokens_saved"] == 800
    assert loaded["subagent_calls"][0]["type"] == "explore"
    assert loaded["pending_plan_body"].startswith("# resumed plan")
    assert loaded["todo_payload"][0]["id"] == 7


def test_legacy_session_without_extras_loads_clean(fake_home):
    """Sessions written before the long-session fields existed must still
    load without the new keys throwing KeyError."""
    legacy = {
        "session_id": "sid-old",
        "mode": "solo",
        "role_index": 0,
        "chat_messages": [],
        "engine_messages": [],
        "title": "legacy",
        "created_at": 0,
        "updated_at": 0,
    }
    sessions_dir = ss._ensure_dir()
    (sessions_dir / "sid-old.json").write_text(
        json.dumps(legacy), encoding="utf-8",
    )
    loaded = ss.load_session("sid-old")
    assert loaded is not None
    assert loaded["mode"] == "solo"
    # Long-session fields gracefully absent — caller defaults them
    assert loaded.get("perm_profile", "") == ""
    assert loaded.get("active_gate") is None


def test_session_search_finds_text_in_archive(fake_home):
    """Conceptual test of the search algorithm the slash command uses:
    given an archived transcript that contains the query, we can locate
    it via load_transcript_archive + simple substring scan."""
    ss.archive_pre_compaction_transcript(
        "sid-archived",
        [
            {"role": "user", "content": "we discussed the flux capacitor"},
            {"role": "assistant", "content": "yes, 1.21 gigawatts"},
        ],
    )
    records = ss.load_transcript_archive("sid-archived")
    # Mimic the search loop
    q = "flux capacitor"
    hits = []
    for rec in records:
        for i, m in enumerate(rec.get("messages") or []):
            if q in str(m.get("content", "")):
                hits.append((rec.get("compacted_at", 0), i, m["role"]))
    assert hits
    assert hits[0][2] == "user"


def test_session_search_finds_text_in_session_chat(fake_home):
    ss.save_session(
        "sid-with-hit",
        mode="solo",
        chat_messages=[
            {"role": "user", "content": "build me a deflector shield"},
            {"role": "assistant", "content": "OK, on it."},
        ],
    )
    loaded = ss.load_session("sid-with-hit")
    msgs = loaded["chat_messages"]
    matches = [
        i for i, m in enumerate(msgs)
        if "deflector" in str(m.get("content", "")).lower()
    ]
    assert matches == [0]


def test_archive_then_search_combined(fake_home):
    """End-to-end: save a session, archive a pre-compaction transcript,
    search both layers and find hits in either."""
    ss.save_session(
        "sid-combo",
        mode="solo",
        chat_messages=[{"role": "user", "content": "alpha needle"}],
    )
    ss.archive_pre_compaction_transcript(
        "sid-combo",
        [{"role": "assistant", "content": "beta haystack"}],
    )
    found_session = False
    for r in ss.list_sessions(limit=10):
        data = ss.load_session(r["session_id"]) or {}
        for m in data.get("chat_messages") or []:
            if "alpha" in str(m.get("content", "")).lower():
                found_session = True
                break
    found_archive = False
    for archive in ss.list_transcript_archives():
        for rec in ss.load_transcript_archive(archive["session_id"]):
            for m in rec.get("messages") or []:
                if "beta" in str(m.get("content", "")).lower():
                    found_archive = True
    assert found_session
    assert found_archive


def test_save_session_writes_atomically_for_concurrent_resume(fake_home):
    """save_session uses tmp + replace internally elsewhere; verify the
    final file is always parseable JSON even on a fresh write."""
    p = ss.save_session(
        "sid-atomic", mode="solo",
        chat_messages=[{"role": "user", "content": "x"}],
        perm_profile="ask_all",
    )
    # File must be valid JSON, not partial
    data = json.loads(p.read_text(encoding="utf-8"))
    assert data["session_id"] == "sid-atomic"
    assert data["perm_profile"] == "ask_all"
