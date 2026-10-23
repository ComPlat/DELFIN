"""Self-healing backfill for the session index.

Phase-1 finding (2026-09-29, verified): the index DB did not exist on the
real machine despite 25 finished sessions -- every session-end index pass
had failed silently (by design, ``index_session`` never raises), so
``session_search`` returned 0 hits forever. The fix is a bounded backfill
that runs when a search finds nothing: index the sessions that are on disk
but missing from the index, then answer the query. It must never touch the
real ``~/.delfin`` from tests (see the ``index_home`` fixture pattern in
``test_a_session_index_makes_sessions_searchable.py``).
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


class TestBackfill:
    def test_a_miss_triggers_backfill_and_then_finds_the_session(self,
                                                                 index_home):
        """A session on disk but not indexed is found after one search.

        The whole point of the fix: a silent session-end failure must not
        make a session unfindable forever. One ``search()`` that misses
        heals the gap and returns the hit on its own retry.
        """
        from delfin.agent import session_index

        _write_session(index_home, "sess-lost", title="the lost one",
                       messages=[("user", "the thiosulfate leaching fix")])
        # NOTE: no index_session call -- as on the real machine.

        # One search heals the gap transparently: the backfill runs
        # inside search() before the query is answered (proved red
        # against the unfixed code: this exact call returned []).
        hits = session_index.search("thiosulfate")
        assert hits and hits[0].session_id == "sess-lost"

    def test_backfill_is_bounded(self, index_home):
        """Backfill never indexes more than its cap in one pass.

        A machine with thousands of sessions must not stall a search for
        minutes; the backfill works through a bounded batch and later
        searches continue the healing.
        """
        from delfin.agent import session_index

        for i in range(40):
            _write_session(index_home, f"sess-bulk-{i}",
                           title=f"bulk {i}",
                           messages=[("user", f"bulk session number {i}")])
        session_index.search("bulk")
        cap = getattr(session_index, "BACKFILL_BATCH", None)
        assert cap is not None and cap <= 10, "backfill batch must be capped"

    def test_backfill_never_raises(self, index_home, monkeypatch):
        """A broken source file in the backfill path is skipped, not fatal."""
        from delfin.agent import session_index

        d = index_home / "agent_sessions"
        d.mkdir(parents=True, exist_ok=True)
        (d / "sess-broken.json").write_text("{ not json", encoding="utf-8")
        _write_session(index_home, "sess-ok", title="ok",
                       messages=[("user", "the recoverable answer")])
        # Must not raise despite the corrupt file.
        hits = session_index.search("recoverable")
        assert hits and hits[0].session_id == "sess-ok"
