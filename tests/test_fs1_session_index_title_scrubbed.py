# FS1 finding (nacht-s6 probe20): session_index stores the session TITLE
# unscrubbed on disk. documents INSERT scrubs text2 but passes title2 raw
# (session_index.py:397), and the title's provenance is the raw session
# json (:263) -- a secret in the title lands in the FTS index on disk and
# comes back verbatim through _search_locked (Hit.title, :532).
import json

import pytest

from delfin.agent import output_guard

SECRET = "sk-" + "live-" + "a1B2c3D4e5F6g7H8i9J0kLmN"


@pytest.fixture
def seeded_session(tmp_path, monkeypatch):
    """A session store + transcript archive pointing at tmp_path."""
    from delfin.agent import session_index as si
    from delfin.agent import session_store as ss

    monkeypatch.setattr(ss, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(
        ss, "_transcript_archive_path",
        lambda: tmp_path / "archive")
    (tmp_path / "sessions").mkdir(parents=True)
    monkeypatch.setattr(si, "_index_path", lambda: tmp_path / "idx.sqlite")
    yield tmp_path


def _write_session_json(path, title):
    path.write_text(json.dumps({
        "title": title,
        "created_at": 1700000000.0,
        "chat_messages": [{"role": "user", "content": "hello world"}],
    }))


def test_session_title_with_secret_is_scrubbed_on_disk(seeded_session):
    from delfin.agent import session_index as si

    _write_session_json(
        seeded_session / "sessions" / "s-title.json",
        f"My api key {SECRET}")
    si.index_session("s-title")

    import sqlite3
    con = sqlite3.connect(seeded_session / "idx.sqlite")
    rows = con.execute(
        "SELECT title FROM documents").fetchall()
    con.close()
    assert rows, "the session was indexed"
    for (title,) in rows:
        assert SECRET not in title, "the title hit disk unscrubbed"


def test_session_title_with_secret_is_scrubbed_when_redactor_fails(
        seeded_session, monkeypatch):
    # FS1 contract: the redactor failing closed must withhold, never leak.
    from delfin.agent import session_index as si

    def raising(text, **kwargs):
        raise RuntimeError("redactor down")

    monkeypatch.setattr(output_guard, "scrub_secrets", raising)
    _write_session_json(
        seeded_session / "sessions" / "s-title2.json",
        f"My api key {SECRET}")
    si.index_session("s-title2")

    import sqlite3
    con = sqlite3.connect(seeded_session / "idx.sqlite")
    rows = con.execute(
        "SELECT title FROM documents").fetchall()
    con.close()
    assert rows
    for (title,) in rows:
        assert SECRET not in title
