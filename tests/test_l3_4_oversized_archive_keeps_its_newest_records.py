"""A transcript archive above the read cap is indexed from its newest end.

Input: an append-only archive (one JSON record per compaction) whose size
exceeds ``session_index.MAX_SOURCE_BYTES``.
Output: the newest whole records within the cap are searchable; the
oldest records and the torn record at the cut are not.

Before this, an archive above the cap was skipped whole, so the longest
sessions (the ones with the most compactions) became unsearchable, where
the unbounded read before the cap had indexed them.
"""

from __future__ import annotations

import json
import time

import pytest

import delfin.agent.session_index as si


@pytest.fixture
def index_home(tmp_path, monkeypatch):
    from delfin.agent import session_store

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(si, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


def _append_record(path, text):
    rec = {"compacted_at": time.time(), "n_messages": 1, "info": {},
           "messages": [{"role": "user", "content": text}]}
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")


def test_an_archive_above_the_cap_keeps_its_newest_records(
        index_home, monkeypatch):
    monkeypatch.setattr(si, "MAX_SOURCE_BYTES", 4096)
    d = index_home / "transcript_archive"
    d.mkdir()
    p = d / "longsession.jsonl"
    _append_record(p, "oldestmarker " + "x" * 3000)
    for i in range(3):
        _append_record(p, f"filler{i} " + "y" * 900)
    _append_record(p, "newestmarker")
    assert p.stat().st_size > si.MAX_SOURCE_BYTES

    assert si.index_session("longsession") is True
    assert [h.session_id for h in si.search("newestmarker")] == [
        "longsession"]
    assert si.search("oldestmarker") == []


def test_an_archive_within_the_cap_is_indexed_whole(index_home, monkeypatch):
    monkeypatch.setattr(si, "MAX_SOURCE_BYTES", 4096)
    d = index_home / "transcript_archive"
    d.mkdir()
    p = d / "short.jsonl"
    _append_record(p, "firstmarker")
    _append_record(p, "secondmarker")

    assert si.index_session("short") is True
    assert si.search("firstmarker")
    assert si.search("secondmarker")


def test_an_oversized_session_json_is_still_skipped(index_home, monkeypatch):
    """The session json is one JSON document: a tail of it parses to
    nothing, so it stays skipped above the cap."""
    monkeypatch.setattr(si, "MAX_SOURCE_BYTES", 256)
    d = index_home / "agent_sessions"
    d.mkdir()
    (d / "big.json").write_text(json.dumps(
        {"title": "t", "chat_messages": [
            {"role": "user", "content": "bigmarker " + "z" * 500}]}))
    assert si.index_session("big") is False
    assert si.search("bigmarker") == []
