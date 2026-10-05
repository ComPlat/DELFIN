"""T3 phase 2 — reliable sessions: reserved operator mailbox, no refusal for
closed-but-known sessions.

A message to ``operator`` was refused with "no other open session 'operator'"
unless a heartbeat process faked presence, because the recipient is resolved
only among open sessions (api_client.py:17654). This package adds, on the
storage side, ``deliverable(to_key)``: a session is reachable when it is the
reserved operator mailbox, has ever announced itself, or already has an
inbox — so a known-but-closed session queued instead of refused.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def test_the_operator_is_always_deliverable(tmp_path):
    """The operator's mailbox is reserved: a message is queued for it whether
    or not a presence record exists (no heartbeat may fake one)."""
    assert M.deliverable("operator") is True
    # Presence is empty: the operator has not announced itself.
    assert P.open_sessions() == []
    # The queue still accepts the message and it waits to be taken.
    sent = M.send("operator", "phase 2 done", from_key="nacht-s16")
    assert sent["to"] == "operator"
    (msg,) = M.take("operator")
    assert msg["text"] == "phase 2 done"


def test_a_known_but_closed_session_is_deliverable(tmp_path, monkeypatch):
    """A session whose presence record exists but is stale (it closed, or is
    mid-restart) is still reachable: the message queues for its next start
    instead of being refused."""
    P.announce("nacht-s12", title="builder", workspace=str(tmp_path))
    # Age the record past the staleness window so it is no longer "open",
    # exactly a session that has closed or just restarted.
    record_path = P._path("nacht-s12")
    old = json.loads(record_path.read_text(encoding="utf-8"))
    old["updated_at"] = time.time() - P._STALE_S - 60
    record_path.write_text(json.dumps(old), encoding="utf-8")
    # No open_sessions() call here: it would run `_reap` and delete the same-
    # host stale record, which is exactly the mid-restart window we must cover.
    # (presence._reap reaps same-host dead records on any open_sessions check.)
    assert M.deliverable("nacht-s12") is True   # stale but known: queue it
    M.send("nacht-s12", "wake up after your restart", from_key="nacht-s16")
    (msg,) = M.take("nacht-s12")
    assert msg["text"] == "wake up after your restart"


def test_a_session_with_an_inbox_is_deliverable(tmp_path):
    """A session that has received before (its inbox exists) remains
    deliverable even once its presence has gone stale."""
    M.send("nacht-s13", "hello", from_key="nacht-s16")
    # No presence record: never announced, but it has an inbox.
    assert M.deliverable("nacht-s13") is True


def test_gibberish_is_not_deliverable(tmp_path):
    """An address that is neither reserved, announced nor has an inbox is
    unknown — a typo — and is refused rather than queued forever."""
    assert M.deliverable("no-such-session-xyz") is False


def test_reserved_mailbox_case_variant_reads_the_same_inbox(tmp_path):
    """The reserved check is case-insensitive (deliverable('OPERATOR') is True),
    so a case-variant address must write and read the SAME mailbox file. Without
    this, send('OPERATOR') writes OPERATOR.jsonl while take('operator') reads
    operator.jsonl and the message is silently dropped (reviewer finding B)."""
    assert M.deliverable("OPERATOR") is True
    sent = M.send("OPERATOR", "case-variant queued", from_key="nacht-s16")
    assert sent["to"] == "OPERATOR"
    (msg,) = M.take("operator")
    assert msg["text"] == "case-variant queued"
