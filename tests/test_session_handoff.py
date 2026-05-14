"""Tests for the session-handoff surface:

- build_handoff_brief — structured recap for a fresh agent
- save_handoff_brief — persistence under ~/.delfin/handoffs/
- export_bundle / import_bundle — portable .delfin-session round-trip
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


def _sample_session() -> dict:
    return {
        "session_id": "sid-handoff-1",
        "title": "Refactor the flux capacitor",
        "chat_messages": [
            {"role": "user", "content": "Refactor flux_capacitor.py to use 1.21 GW."},
            {"role": "assistant", "content": "Read the file; the power constant is hard-coded.\n\nNext I'll extract it to config."},
            {"role": "user", "content": "go"},
            {"role": "assistant", "content": "Extracted POWER_GW to config.py.\n\nNext: update the 3 call sites."},
        ],
        "engine_messages": [
            {"role": "assistant", "content": "I read src/flux_capacitor.py and src/config.py."},
            {"role": "tool", "content": "edited src/config.py:12"},
        ],
        "token_usage": {"input": 4200, "output": 1100},
        "cost_usd": 0.0312,
        "todo_payload": [
            {"id": 1, "subject": "extract POWER_GW", "status": "completed"},
            {"id": 2, "subject": "update call sites", "status": "in_progress"},
            {"id": 3, "subject": "add a test", "status": "pending"},
        ],
        "active_gate": None,
        "pending_plan_body": "",
    }


# --- build_handoff_brief ----------------------------------------------------

def test_build_handoff_brief_has_all_sections():
    brief = ss.build_handoff_brief(_sample_session())
    for header in (
        "# Session Handoff",
        "## Goal",
        "## Key decisions & current state",
        "## Files touched",
        "## Open items",
        "## Recommended next step",
    ):
        assert header in brief, header


def test_build_handoff_brief_extracts_goal_from_first_user_msg():
    brief = ss.build_handoff_brief(_sample_session())
    assert "Refactor flux_capacitor.py" in brief


def test_build_handoff_brief_lists_open_tasks_not_completed_ones():
    brief = ss.build_handoff_brief(_sample_session())
    assert "update call sites" in brief      # in_progress
    assert "add a test" in brief              # pending
    assert "extract POWER_GW" not in brief    # completed — not an open item


def test_build_handoff_brief_extracts_files_touched():
    brief = ss.build_handoff_brief(_sample_session())
    assert "src/config.py" in brief or "src/flux_capacitor.py" in brief


def test_build_handoff_brief_surfaces_active_gate_as_open_item():
    data = _sample_session()
    data["active_gate"] = {"type": "approval", "title": "review the edit"}
    brief = ss.build_handoff_brief(data)
    assert "paused on approval gate" in brief


def test_build_handoff_brief_surfaces_pending_plan():
    data = _sample_session()
    data["pending_plan_body"] = "# Plan\n1. do X"
    brief = ss.build_handoff_brief(data)
    assert "plan-mode plan is awaiting approval" in brief


def test_build_handoff_brief_handles_empty_session():
    """A session with no messages must still produce a valid brief
    rather than crashing."""
    brief = ss.build_handoff_brief({"session_id": "empty", "chat_messages": []})
    assert "# Session Handoff" in brief
    assert "no explicit goal" in brief


# --- save_handoff_brief -----------------------------------------------------

def test_save_handoff_brief_writes_file(fake_home):
    brief = ss.build_handoff_brief(_sample_session())
    p = ss.save_handoff_brief("sid-handoff-1", brief)
    assert p.exists()
    assert p.suffix == ".md"
    assert "Session Handoff" in p.read_text(encoding="utf-8")


def test_save_handoff_brief_sanitises_session_id(fake_home):
    """A session id with slashes / dots must not escape the handoffs
    directory."""
    p = ss.save_handoff_brief("../../etc/passwd", "body")
    assert ".delfin/handoffs" in str(p)
    assert "etc/passwd" not in str(p)


# --- export_bundle / import_bundle round-trip -------------------------------

def test_export_bundle_then_import_round_trip(fake_home):
    # Save the source session so export_bundle can find it
    data = _sample_session()
    ss.save_session(
        data["session_id"], mode="solo",
        chat_messages=data["chat_messages"],
        engine_messages=data["engine_messages"],
        token_usage=data["token_usage"],
        cost_usd=data["cost_usd"],
        title=data["title"],
        todo_payload=data["todo_payload"],
    )
    # Also drop a transcript archive so the bundle picks it up
    ss.archive_pre_compaction_transcript(
        data["session_id"],
        [{"role": "user", "content": "older compacted content"}],
    )
    bundle = ss.export_bundle(data["session_id"])
    assert bundle is not None
    assert bundle.suffix == ".delfin-session" or bundle.name.endswith(".delfin-session")
    assert bundle.exists()

    # Import it — should land under a FRESH id
    new_sid = ss.import_bundle(bundle)
    assert new_sid is not None
    assert new_sid != data["session_id"]
    restored = ss.load_session(new_sid)
    assert restored is not None
    assert restored["title"].endswith("(imported)")
    assert len(restored["chat_messages"]) == len(data["chat_messages"])
    # Transcript archive was restored alongside
    archives = {a["session_id"] for a in ss.list_transcript_archives()}
    assert new_sid in archives


def test_export_bundle_missing_session_returns_none(fake_home):
    assert ss.export_bundle("does-not-exist") is None


def test_import_bundle_bad_path_returns_none(fake_home, tmp_path):
    assert ss.import_bundle(tmp_path / "nope.delfin-session") is None


def test_import_bundle_rejects_non_zip(fake_home, tmp_path):
    junk = tmp_path / "junk.delfin-session"
    junk.write_text("this is not a zip", encoding="utf-8")
    assert ss.import_bundle(junk) is None


def test_import_bundle_does_not_clobber_existing(fake_home):
    """Importing must always mint a fresh id so a local session with the
    same id is never overwritten."""
    data = _sample_session()
    ss.save_session(
        data["session_id"], mode="solo",
        chat_messages=data["chat_messages"], title=data["title"],
    )
    bundle = ss.export_bundle(data["session_id"])
    new_sid = ss.import_bundle(bundle)
    # Original still intact
    assert ss.load_session(data["session_id"]) is not None
    assert new_sid != data["session_id"]


def test_bundle_contains_handoff_brief(fake_home):
    import zipfile
    data = _sample_session()
    ss.save_session(
        data["session_id"], mode="solo",
        chat_messages=data["chat_messages"], title=data["title"],
    )
    bundle = ss.export_bundle(data["session_id"])
    with zipfile.ZipFile(bundle, "r") as z:
        names = set(z.namelist())
        assert "session.json" in names
        assert "handoff.md" in names
        assert "manifest.json" in names
        manifest = json.loads(z.read("manifest.json"))
        assert manifest["kind"] == "delfin-session-bundle"
