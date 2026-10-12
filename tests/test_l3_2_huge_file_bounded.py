"""A huge session file does not stall a turn.

The prior index (``delfin/agent/session_index.py``) reads every source with a
whole-file ``p.read_text()`` (``_session_json_documents`` :148,
``_archive_documents`` :174): a multi-GB session json or transcript archive
is read fully into memory and would stall a search. Every source read
must be bounded to a fixed byte cap and return promptly on oversized
files -- a search must stay cheap no matter how large the transcripts were.
"""

from __future__ import annotations

import json
import os
import time

import pytest

import delfin.agent.session_index as sr


@pytest.fixture
def index_home(tmp_path, monkeypatch):
    """A fake ``~/.delfin`` holding the index and the session sources."""
    from delfin.agent import session_store

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(sr, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


def _write_session(home, session_id, *, title="A session", text="needle 42",
                   messages=None):
    d = home / "agent_sessions"
    d.mkdir(parents=True, exist_ok=True)
    msgs = messages or [{"role": "user", "content": text}]
    data = {
        "schema_version": 2,
        "session_id": session_id,
        "title": title,
        "chat_messages": msgs,
        "created_at": time.time() - 86400,
        "updated_at": time.time(),
    }
    p = d / f"{session_id}.json"
    p.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    return p


def _write_archive(home, session_id, messages):
    d = home / "transcript_archive"
    d.mkdir(parents=True, exist_ok=True)
    rec = {
        "compacted_at": time.time(),
        "n_messages": len(messages),
        "info": {},
        "messages": [{"role": r, "content": c} for r, c in messages],
    }
    p = d / f"{session_id}.jsonl"
    with p.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return p


def test_oversized_session_json_does_not_stall_the_search(index_home):
    """A session json larger than the read cap is never read whole: the
    index skips/truncates it and search still returns promptly."""
    _write_session(index_home, "small", text="the small session content")
    d = index_home / "agent_sessions"

    big = d / "huge.json"
    # Far beyond any sane read cap -- writing it as a real file would be
    # slow, so we fake a huge size observed at the stat gate: the module
    # must consult the size BEFORE opening (never open a many-GB file).
    _write_session(index_home, "huge", text="bone of many kilobytes")
    # Make the file's on-disk size lie far beyond the cap by truncating it
    # to a sparse 512 MiB -- cheap, and the stat gate sees a huge size.
    os.truncate(big, 512 * 1024 * 1024)

    # A backfill/index pass must neither hang nor read the huge file whole.
    sr.index_session("huge")
    sr.backfill_missing(batch=100)

    hits = sr.search("small")
    assert any(h.session_id == "small" for h in hits), (
        "the small session is still indexed and searchable")
    # The huge session must not have wedged the index: this search returns.
    sr.search("bone")


def test_oversized_archive_jsonl_is_bounded(index_home):
    """A transcript archive larger than the read cap is read bounded (or
    skipped), never fully loaded into a search's hot path."""
    _write_session(index_home, "core", text="core session text 4711")
    d = index_home / "transcript_archive"
    d.mkdir(parents=True, exist_ok=True)

    big = d / "hugearch.jsonl"
    _write_archive(index_home, "hugearch", [("user", "archive needle 99")])
    os.truncate(big, 256 * 1024 * 1024)  # pretend 256 MiB of archive

    sr.index_session("hugearch")
    sr.backfill_missing(batch=100)
    hits = sr.search("core")
    assert any(h.session_id == "core" for h in hits)
    # And the archive read itself is bounded: this call must return.
    sr.search("archive")


def test_backfill_returns_quickly_with_many_oversized_files(index_home):
    """A wall of oversized source files must not make a backfill stall:
    the size gate is checked before any open."""
    for i in range(5):
        _write_session(index_home, f"s{i}", text=f"session {i} marker")
        p = index_home / "agent_sessions" / f"s{i}.json"
        os.truncate(p, 64 * 1024 * 1024)
    # Must return promptly, indexing at most the bounded ones (which the
    # size gate rejects here), and never attempt a whole-file read.
    assert sr.backfill_missing(batch=100) == 0
    hits = sr.search("marker")
    assert hits == [], (
        "oversized sessions are skipped by the size gate, so none is "
        "searchable -- and backfill reports 0 indexed rather than stalling")
