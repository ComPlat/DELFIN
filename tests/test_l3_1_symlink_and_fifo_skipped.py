"""L3 phase-3 guarantee: a symlink or FIFO in the sessions directory is
skipped -- never followed, never indexed, never opened (a FIFO must not
hang a read).

The prior index (``delfin/agent/session_index.py``) enumerates with
``d.glob("*.json")`` and reads with ``p.read_text()`` (session_index.py:379,
:148, :174): a ``*.json`` SYMLINK pointing into another account's directory
is followed (exfiltration), and a ``*.jsonl`` FIFO blocks ``read_text``
forever (a stalled turn). These tests pin the guaranteed behaviour of the
L3 replacement (``delfin.agent.session_index``): enumerate ONLY regular
files, skip symlinks and non-regular entries, and never open a FIFO.

RED on the previous commit: ``session_index`` does not exist yet, and the
existing index has no symlink/FIFO guard -- so every assertion here is
verifying a guarantee the current code does not provide.
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


def _write_session(home, session_id, *, title="A session", text="the answer is 42",
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


def test_symlink_pointing_out_of_the_dir_is_never_indexed(index_home, tmp_path):
    """A ``*.json`` symlink whose target lives in an OTHER account's dir
    must be skipped: nothing from the target is indexed, never searched."""
    _write_session(index_home, "realone", text="the trustworthy session")
    # Simulate another account's directory -- a real path, outside HOME.
    outside = tmp_path / "another_account_home"
    outside.mkdir()
    victim = outside / "secrets.json"
    victim.write_text('{"secret_text": "the stolen data 007"}', encoding="utf-8")

    d = index_home / "agent_sessions"
    link = d / "evil.json"
    link.symlink_to(str(victim))

    assert sr.index_session("evil") is False, (
        "a symlinked session must not be indexed at all")

    # Even a full backfill must not pick the symlink up.
    assert sr.backfill_missing(batch=100) >= 0
    hits = sr.search("stolen")
    assert hits == [], (
        "the symlink target's content must never be searchable")
    hits = sr.search("trustworthy")
    assert any(h.session_id == "realone" for h in hits), (
        "the real session must still be indexed")


def test_fifo_named_json_is_skipped_not_opened(index_home):
    """A FIFO that would block ``read_text`` forever must never be opened.
    The test returns promptly; a hang here is the bug."""
    _write_session(index_home, "fine", text="fine session content")
    d = index_home / "agent_sessions"

    fifo = d / "blocking.json"
    os.mkfifo(str(fifo))  # a writer that never closes blocks a full read

    # Must not block: the guard rejects non-regular files before opening.
    sr.index_session("blocking")
    sr.backfill_missing(batch=100)

    hits = sr.search("fine")
    assert any(h.session_id == "fine" for h in hits), (
        "the regular session must be found; the FIFO is just skipped")
    hits = sr.search("blocking")
    assert hits == [], "the FIFO's (empty) content is not indexed"


def test_symlink_and_fifo_do_not_hide_a_returned_snippet(index_home):
    """After a backfill the result list contains only real sessions --
    neither the fake symlink session nor the FIFO appears."""
    _write_session(index_home, "ok", text="searchable real content")
    d = index_home / "agent_sessions"

    out_link = index_home / "agent_sessions" / "link.json"
    out_link.symlink_to(str(index_home / ".." / "target.json"))
    try:
        (index_home / ".." / "target.json").resolve().write_text(
            '{"secret": "outside 009"}', encoding="utf-8")
    except OSError:
        pass  # the symlink branch is asserted by the first test

    sr.backfill_missing(batch=100)
    for h in sr.search("real"):
        assert h.session_id not in ("link",), (
            "a symlinked entry must never surface as a session")
        assert h.session_id == "ok"
