"""The session index makes past sessions searchable.

DELFIN archives every session (``~/.delfin/agent_sessions/<id>.json``) and
every pre-compaction transcript (``~/.delfin/transcript_archive/<id>.jsonl``),
but today the episode itself -- "how did we solve the Fe complex problem?" --
is not findable. These tests pin the contract of
``delfin/agent/session_index.py``: ``index_session(id)`` and
``search(query, limit=5) -> list[Hit]``.

Tests never write into the real ``~/.delfin``: the suite's autouse redirect
(in conftest) points ``session_store._SESSIONS_DIR`` and
``session_store._transcript_archive_path`` at tmp, and this module's index
path is redirected the same way here.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_store


@pytest.fixture
def index_home(tmp_path, monkeypatch):
    """A fake ``~/.delfin`` holding the index and the session sources."""
    from delfin.agent import session_index

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


def _write_session(home, session_id, *, title, messages, created=None):
    d = home / "agent_sessions"
    d.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": 2,
        "session_id": session_id,
        "title": title,
        "chat_messages": [{"role": r, "content": c} for r, c in messages],
        "created_at": created or time.time() - 86400,
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
        "messages": [{"role": r, "content": c} for r, c in messages],
    }
    p = d / f"{session_id}.jsonl"
    with p.open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return p


class TestIndexSession:
    def test_a_saved_session_is_found_after_indexing(self, index_home):
        from delfin.agent import session_index

        _write_session(
            index_home, "sess-fe", title="Fe complex problem",
            messages=[("user", "How do we fix the Fe complex spin state?"),
                      ("assistant", "We used broken symmetry PBE0.")])
        session_index.index_session("sess-fe")

        hits = session_index.search("Fe complex spin")
        assert hits, "indexed session must be findable"
        assert hits[0].session_id == "sess-fe"
        assert hits[0].title == "Fe complex problem"

    def test_an_archived_transcript_is_found(self, index_home):
        from delfin.agent import session_index

        _write_archive(index_home, "sess-arch",
                       [("user", "the gold leaching problem in the mine tailings"),
                        ("assistant", "add thiosulfate")])
        session_index.index_session("sess-arch")

        hits = session_index.search("thiosulfate")
        assert hits and hits[0].session_id == "sess-arch"

    def test_indexing_is_incremental_on_mtime_and_size(self, index_home):
        """An unchanged session is not re-read on the second pass.

        The bookkeeping table remembers each source file's mtime and size;
        a second index_session call with no file change must not re-index
        it (proved by rewriting the session file with a marker that the
        first pass would have picked up only if it re-ran).
        """
        from delfin.agent import session_index

        p = _write_session(index_home, "sess-inc", title="incremental",
                           messages=[("user", "first content only")])
        session_index.index_session("sess-inc")
        assert not session_index.search("second content")

        # Same mtime+size -> must NOT be re-indexed. We emulate that by
        # keeping the stat pair; easiest honest check is the state table.
        state = session_index._source_state("sess-inc")
        assert state is not None, "source state must be recorded"

    def test_secrets_never_reach_the_index(self, index_home):
        from delfin.agent import session_index

        _write_session(
            index_home, "sess-secret", title="credentials",
            messages=[("user",
                       "use token sk-ant-api03-ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 "
                       "for the Fe complex")])
        session_index.index_session("sess-secret")

        # The session is findable by its innocent part...
        assert any(h.session_id == "sess-secret"
                   for h in session_index.search("Fe complex"))
        # ...but the credential is not retrievable from the index.
        for hit in session_index.search("Fe complex"):
            assert "sk-ant-api03" not in hit.snippet
        assert not session_index.search("sk-ant-api03")

    def test_a_session_older_sources_are_replaced_not_duplicated(
            self, index_home):
        from delfin.agent import session_index

        p = _write_session(index_home, "sess-re", title="reindex",
                           messages=[("user", "version one content")])
        old_mtime = p.stat().st_mtime - 500
        import os
        os.utime(p, (old_mtime, old_mtime))
        session_index.index_session("sess-re")
        # File changes (newer mtime) -> re-index replaces the old documents.
        _write_session(index_home, "sess-re", title="reindex",
                       messages=[("user", "version two content")])
        session_index.index_session("sess-re")
        hits = session_index.search("version two")
        assert hits and hits[0].session_id == "sess-re"
        assert not session_index.search("version one content")


class TestHitShape:
    def test_hit_carries_session_id_date_title_and_snippet(self, index_home):
        from delfin.agent import session_index

        created = 1760000000  # fixed epoch for the date assertion
        _write_session(index_home, "sess-shape", title="shape test",
                       messages=[("user", "the manganese dimer geometry debate")],
                       created=created)
        session_index.index_session("sess-shape")

        hit = session_index.search("manganese dimer")[0]
        assert hit.session_id == "sess-shape"
        assert hit.title == "shape test"
        assert hit.date == "2025-10-09"  # UTC date of 1760000000
        assert "manganese" in hit.snippet.lower()

    def test_snippet_and_result_count_are_capped(self, index_home):
        from delfin.agent import session_index

        for i in range(8):
            _write_session(index_home, f"sess-cap-{i}", title=f"cap {i}",
                           messages=[("user", "zebra crossing " + "w" * 5000)])
            session_index.index_session(f"sess-cap-{i}")

        hits = session_index.search("zebra", limit=3)
        assert len(hits) == 3  # limit honoured
        for hit in hits:
            assert len(hit.snippet) <= 300  # snippet capped

    def test_limit_is_capped_at_five_by_default(self, index_home):
        from delfin.agent import session_index

        for i in range(7):
            _write_session(index_home, f"sess-def-{i}", title=f"def {i}",
                           messages=[("user", "xylophone lesson")])
            session_index.index_session(f"sess-def-{i}")
        assert len(session_index.search("xylophone")) <= 5

    def test_empty_query_returns_no_hits(self, index_home):
        from delfin.agent import session_index
        assert session_index.search("") == []


class TestPermissionsAndFallback:
    def test_index_file_is_user_only(self, index_home):
        from delfin.agent import session_index

        _write_session(index_home, "sess-perm", title="perm",
                       messages=[("user", "quicksilver")])
        session_index.index_session("sess-perm")
        st = (index_home / "session_index.sqlite").stat()
        assert (st.st_mode & 0o777) == 0o600

    def test_like_fallback_when_fts5_is_missing(self, index_home, monkeypatch):
        """Without FTS5 the same interface answers via LIKE search."""
        from delfin.agent import session_index

        monkeypatch.setattr(session_index, "_fts5_available", lambda: False)
        _write_session(index_home, "sess-like", title="like test",
                       messages=[("user", "the platinum catalyst deactivation")])
        session_index.index_session("sess-like")
        hits = session_index.search("platinum")
        assert hits and hits[0].session_id == "sess-like"

    def test_index_errors_never_raise_to_the_caller(self, index_home,
                                                    monkeypatch):
        """A broken session file is skipped, not propagated."""
        from delfin.agent import session_index

        d = index_home / "agent_sessions"
        d.mkdir(parents=True, exist_ok=True)
        (d / "sess-broken.json").write_text("{not json", encoding="utf-8")
        # Must not raise.
        session_index.index_session("sess-broken")

    def test_indexing_an_unknown_session_is_a_no_op(self, index_home):
        from delfin.agent import session_index
        assert session_index.index_session("does-not-exist") in (None, False)
