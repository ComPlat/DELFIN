"""The session-end re-index step never breaks the session end.

Paket 5 Phase 4: whatever the session-end path of Paket 2 looks like,
the indexing step it calls must be one safe line. These tests pin that:
a finished session becomes searchable through
``reindex_finished_session``, and no failure of the index -- not a
missing session, not a corrupt index file, not a broken import -- raises
to the caller.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import session_store


@pytest.fixture
def index_home(tmp_path, monkeypatch):
    from delfin.agent import session_index

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


class TestReindexFinishedSession:
    def test_a_finished_session_becomes_searchable(self, index_home):
        from delfin.agent import session_end_index, session_index

        d = index_home / "agent_sessions"
        d.mkdir(parents=True)
        (d / "sess-end.json").write_text(json.dumps({
            "session_id": "sess-end", "title": "the hafnium alkyl story",
            "created_at": 1760000000,
            "chat_messages": [
                {"role": "user", "content": "hafnium alkyl co-catalyst tuning"}],
        }), encoding="utf-8")
        assert session_end_index.reindex_finished_session("sess-end") is True
        assert session_index.search("hafnium alkyl")[0].session_id == "sess-end"

    def test_an_unknown_session_returns_false_quietly(self, index_home):
        from delfin.agent import session_end_index
        assert session_end_index.reindex_finished_session("nope") is False

    def test_an_empty_session_id_returns_false_quietly(self, index_home):
        from delfin.agent import session_end_index
        assert session_end_index.reindex_finished_session("") is False

    def test_a_corrupt_index_file_does_not_raise(self, index_home):
        from delfin.agent import session_end_index, session_index

        d = index_home / "agent_sessions"
        d.mkdir(parents=True)
        (d / "sess-corrupt.json").write_text(json.dumps({
            "session_id": "sess-corrupt", "title": "t",
            "created_at": 1.0, "chat_messages": [
                {"role": "user", "content": "ytterbium reduction"}],
        }), encoding="utf-8")
        # Corrupt the index so the write must fail.
        index_home.joinpath("session_index.sqlite").write_bytes(b"not sqlite")
        assert session_end_index.reindex_finished_session(
            "sess-corrupt") is False
        # The session itself is untouched.
        assert (d / "sess-corrupt.json").exists()

    def test_a_broken_session_index_module_does_not_raise(
            self, index_home, monkeypatch):
        from delfin.agent import session_end_index

        def _boom(session_id):
            raise RuntimeError("index exploded")

        monkeypatch.setattr(
            "delfin.agent.session_index.index_session", _boom)
        assert session_end_index.reindex_finished_session("x") is False
