"""Continuity guarantees for long conversations.

Two promises the user relies on:
  1. Auto-compaction never drops the user's GOALS — the original intent
     survives verbatim (up to 400 chars) in the in-place summary.
  2. The full pre-compaction transcript is archived, so nothing is ever
     truly lost — it can be reloaded after the live history was trimmed.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def tmp_archive(monkeypatch, tmp_path):
    """Redirect the transcript archive to a throwaway dir."""
    d = tmp_path / "transcript_archive"
    d.mkdir()
    monkeypatch.setattr(ss, "_transcript_archive_dir", lambda: d)
    return d


# ---------------------------------------------------------------------------
# Archive round-trip — nothing is lost
# ---------------------------------------------------------------------------

def test_archive_round_trips_full_transcript(tmp_archive):
    msgs = [
        {"role": "user", "content": "GOAL: build feature X"},
        {"role": "assistant", "content": "ok, here is the plan"},
        {"role": "tool", "content": '{"result": 42}'},
    ]
    ss.archive_pre_compaction_transcript(
        "sess-A", msgs, info={"messages_compacted": 3},
    )
    records = ss.load_transcript_archive("sess-A")
    assert len(records) == 1
    archived = records[0]["messages"]
    assert archived[0]["content"] == "GOAL: build feature X"
    assert records[0]["info"]["messages_compacted"] == 3


def test_archive_is_append_only_across_compactions(tmp_archive):
    ss.archive_pre_compaction_transcript("sess-B", [{"role": "user", "content": "first"}])
    ss.archive_pre_compaction_transcript("sess-B", [{"role": "user", "content": "second"}])
    records = ss.load_transcript_archive("sess-B")
    assert len(records) == 2                       # both compactions kept
    assert records[0]["messages"][0]["content"] == "first"
    assert records[1]["messages"][0]["content"] == "second"


def test_empty_session_id_is_noop(tmp_archive):
    ss.archive_pre_compaction_transcript("", [{"role": "user", "content": "x"}])
    assert ss.load_transcript_archive("") == []


# ---------------------------------------------------------------------------
# Engine compaction preserves user goals
# ---------------------------------------------------------------------------

def _make_engine(tmp_path):
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine(repo_dir=str(tmp_path), mode="solo")
    eng.session_id = "sess-compact"
    return eng


def test_full_compaction_keeps_user_goal_and_recent(monkeypatch, tmp_path, tmp_archive):
    eng = _make_engine(tmp_path)
    # Force the deterministic extractive path (no API summariser).
    monkeypatch.setattr(eng, "_llm_summarize_old_messages", lambda *a, **k: "")

    goal = "GOALSENTINEL implement the dark-mode toggle end to end"
    msgs = [{"role": "user", "content": goal}]
    for i in range(13):                            # well over _COMPACTION_THRESHOLD
        role = "assistant" if i % 2 == 0 else "user"
        msgs.append({"role": role, "content": f"turn {i} chatter " * 5})
    last = "RECENTSENTINEL the very latest exchange"
    msgs[-1] = {"role": "user", "content": last}
    eng.messages = msgs

    assert len(eng.messages) >= eng._COMPACTION_THRESHOLD
    eng._compact_history()

    blob = "\n".join(
        m.get("content", "") for m in eng.messages
        if isinstance(m.get("content"), str)
    )
    # The original GOAL survived the compaction...
    assert "GOALSENTINEL" in blob
    # ...and the most-recent message is still intact (KEEP_RECENT window).
    assert "RECENTSENTINEL" in blob
    # Compaction was recorded for the /context view.
    assert eng.last_compaction_info
    assert eng.last_compaction_info.get("messages_compacted", 0) > 0


def test_full_compaction_archives_old_messages(monkeypatch, tmp_path, tmp_archive):
    eng = _make_engine(tmp_path)
    monkeypatch.setattr(eng, "_llm_summarize_old_messages", lambda *a, **k: "")

    goal = "ARCHIVEGOAL refactor the api client into modules"
    msgs = [{"role": "user", "content": goal}]
    for i in range(13):
        role = "assistant" if i % 2 == 0 else "user"
        msgs.append({"role": role, "content": f"detail {i} " * 5})
    eng.messages = msgs
    eng._compact_history()

    # The pre-compaction transcript is recoverable from the archive.
    records = ss.load_transcript_archive("sess-compact")
    assert records, "no archive written"
    archived_blob = "\n".join(
        m["content"] for rec in records for m in rec["messages"]
    )
    assert "ARCHIVEGOAL" in archived_blob
