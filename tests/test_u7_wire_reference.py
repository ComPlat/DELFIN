"""U7 / phase 3 — the stale marker is wired to the session's last activity.

QS decision (s21, Option B): a stale blocked reason only reaches the user on
the real panel if a call site passes the reference. This tests the helper
that supplies it: ``task_ticker.session_last_activity_at`` reads the
session's last recorded activity (the proxy for "the last message") from the
session store, defensively.

The session store is global keyed by session id; every test here points
``session_store._SESSIONS_DIR`` at a tmp dir so nothing touches the real
store. The store-write floor and ceiling are planted directly (a dict with a
controlled ``updated_at``), so the tests are deterministic — no wall-clock
dependence.

The integration case is the user-visible effect of the wiring: a blocked
task whose block time predates the session's last activity reads as stale on
the panel.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import agent_tasks, session_store
from delfin.agent import task_ticker as TT


@pytest.fixture
def store(tmp_path):
    agent_tasks._STORES.clear()
    return agent_tasks.get_store(tmp_path)


@pytest.fixture
def sessions(tmp_path, monkeypatch):
    """A throwaway session store dir, and a helper to plant a record in it."""
    d = tmp_path / "sessions"
    d.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", d)

    def plant(sid, updated_at):
        (d / f"{sid}.json").write_text(json.dumps({
            "session_id": sid, "updated_at": updated_at,
        }))

    return plant


# ---------------------------------------------------------------------------
# the helper, defensively
# ---------------------------------------------------------------------------

def test_the_helper_reads_the_saved_epoch(store, sessions):
    sessions("s", 1_700_000_000.0)
    assert TT.session_last_activity_at("ws", session_id="s") == 1_700_000_000.0


def test_no_session_record_yields_none(store, sessions):
    assert TT.session_last_activity_at("ws", session_id="ghost") is None


def test_a_missing_or_blank_session_id_yields_none(store, sessions):
    sessions("s", 1_700_000_000.0)
    assert TT.session_last_activity_at("ws") is None
    assert TT.session_last_activity_at("ws", session_id="") is None


def test_a_non_numeric_updated_at_yields_none(store, sessions):
    sessions("s", "not-a-number")
    assert TT.session_last_activity_at("ws", session_id="s") is None


def test_a_zero_updated_at_yields_none(store, sessions):
    sessions("s", 0)
    assert TT.session_last_activity_at("ws", session_id="s") is None


def test_a_corrupt_session_file_never_raises(store, tmp_path, monkeypatch):
    d = tmp_path / "sessions2"
    d.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", d)
    (d / "s.json").write_text("{not json")
    assert TT.session_last_activity_at("ws", session_id="s") is None


# ---------------------------------------------------------------------------
# the integration: marker on the panel when the block predates the activity
# ---------------------------------------------------------------------------

def _blocked(store, subject, reason="waiting on approval"):
    t = store.create(subject, session_id="s")
    return store.update(t["id"], status="blocked", blocked_reason=reason)


def test_marker_appears_when_the_block_predates_the_last_activity(store, sessions, tmp_path):
    _blocked(store, "push the branch")
    block_epoch = TT._block_epoch(store.get(1))
    sessions("s", block_epoch + 3600.0)          # the user said more since
    last = TT.session_last_activity_at(tmp_path, session_id="s")
    row = TT.rows(tmp_path, session_id="s", last_message_at=last)
    assert " (stale)" in row[0]["label"]


def test_no_marker_when_the_block_postdates_the_last_activity(store, sessions, tmp_path):
    _blocked(store, "push the branch")
    block_epoch = TT._block_epoch(store.get(1))
    sessions("s", max(0.0, block_epoch - 3600.0))  # the block is the newest thing
    last = TT.session_last_activity_at(tmp_path, session_id="s")
    row = TT.rows(tmp_path, session_id="s", last_message_at=last)
    assert " (stale)" not in row[0]["label"]


def test_no_marker_when_there_is_no_session_at_all(store, sessions, tmp_path):
    _blocked(store, "push the branch")
    assert session_store.load_session("s") is None
    last = TT.session_last_activity_at(tmp_path, session_id="s")
    assert last is None
    row = TT.rows(tmp_path, session_id="s", last_message_at=last)
    assert " (stale)" not in row[0]["label"]
