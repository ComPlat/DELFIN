"""The elided-content store's silent failures do not erase silently.

Two contracts pinned here:

1. ``_enforce_elided_cap``: dropping the oldest records is best-effort —
   but a store that cannot be READ (OSError) must not be REWRITTEN as
   if it were empty. The old code caught the OSError from the read loop
   and silently left the file as it was, which is correct; what it must
   NEVER do is catch a write failure after the read succeeded and leave
   a store whose cap is silently no longer enforced while every caller
   believes it is. A failed rewrite is now said once (module logger), and
   the store stays whole — dropping-oldest either happened or did not.

2. ``append_elided_record``: returning None on OSError is the documented
   best-effort contract (compaction must not break). But the LOSS is not
   invisible at every consumer: ``history_get`` already returns a
   readable error for a missing ref. What is pinned here is that a
   record whose APPEND failed resolves to the same readable "may have
   been dropped" error, not to a KeyError or a crash.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    return d


def test_enforce_cap_failure_keeps_store_whole_and_is_said(
        fake_home, monkeypatch, caplog):
    """A failing rewrite after a successful read keeps the store whole;
    the failure is visible in the log instead of silent."""
    ref = ss.append_elided_record(
        "sid-capfail", index=0, role="assistant", content="niobium record")
    assert ref
    # A second record so the store exceeds the cap and the rewrite
    # actually fires (with one line, total == cap and nothing happens).
    ss.append_elided_record(
        "sid-capfail", index=1, role="assistant", content="rhodium record")
    store = ss.elided_store_path("sid-capfail")
    before = store.read_text()

    def boom(path, text):
        raise OSError("simulated rewrite failure")

    monkeypatch.setattr(ss, "_atomic_write_text", boom)
    with caplog.at_level("WARNING", logger="delfin.agent.session_store"):
        ss._enforce_elided_cap(store, cap=1)
    # The store is whole: the read succeeded, the rewrite did not, and
    # nothing truncated it.
    assert store.read_text() == before
    assert json.loads(before.splitlines()[0])["ref"] == ref
    # And the failure is not silent.
    assert any("elided store" in r.getMessage() for r in caplog.records)


def test_enforce_cap_healthy_run_silent(fake_home, caplog):
    """The normal path gains no new log speech."""
    for i in range(3):
        ss.append_elided_record(
            "sid-capok", index=i, role="assistant", content=f"rec {i}")
    with caplog.at_level("WARNING", logger="delfin.agent.session_store"):
        ss._enforce_elided_cap(ss.elided_store_path("sid-capok"), cap=2)
    assert not [r for r in caplog.records if "elided store" in r.message]


def test_append_failure_resolves_to_readable_error(fake_home, monkeypatch):
    """A record whose append failed resolves through history_get's error
    path — never a crash — when the caller keeps the marker."""
    from delfin.agent import history_search as hs

    # Direct OSError path: append returns None by contract.
    monkeypatch.setattr(ss, "_ensure_dir", lambda: None)
    p = ss.elided_store_path("sid-appendfail")

    def open_boom(*a, **k):
        raise OSError("simulated append failure")

    monkeypatch.setattr(type(p), "open", open_boom)
    assert ss.append_elided_record(
        "sid-appendfail", index=0, role="assistant", content="x") is None
    # The consumer's error path stays readable.
    got = hs.history_get("sid-appendfail", "elided:deadbeef")
    assert "error" in got
    assert "may have been dropped" in got["error"]
