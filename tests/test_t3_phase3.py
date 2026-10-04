"""T3 phase 3 — delivery receipts and an address-latch bound.

A sender must be able to ask what happened to a message it sent:
`status(message_id)` reports queued / delivered / read (or unknown), and
`ls()` lists the messages it knows about — the storage-side core a
`delfin-agent messages ls` / `status` CLI and the session_message tool burn
into. Every message gets an id at send; a delivered/read receipt survives the
inbox being taken, so the sender can look it up later.

Also the address lattice: reviewer nacht-s17 showed `_known_key` latches any
address that has an inbox file FOREVER (no TTL), so one message to a typo
(e.g. ``operatro``) makes it permanently deliverable and every later message
to it is silently accepted and dropped — the wave-11/12 "sender never knows
if it was read" bug re-created on the typo path. An inbox therefore proves an
address only while it is fresh.
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


def test_send_returns_a_message_id(tmp_path):
    sent = M.send("nacht-s17", "hello", from_key="nacht-s16")
    assert sent["to"] == "nacht-s17"
    mid = sent.get("id")
    assert mid, "send() must hand the sender a message id to ask status of"


def test_status_is_queued_then_delivered_then_read(tmp_path):
    sent = M.send("nacht-s17", "do you have the hash?", from_key="nacht-s16")
    mid = sent["id"]
    assert M.status(mid) == "queued"
    (msg,) = M.take("nacht-s17")
    assert msg["id"] == mid
    assert M.status(mid) == "delivered"
    M.mark_read(mid)
    assert M.status(mid) == "read"


def test_status_survives_the_inbox_being_taken(tmp_path):
    """The receipt outlives take(): take() unlinks the inbox, yet the sender
    can still ask status once the message has been delivered/read."""
    sent = M.send("nacht-s17", "heads-up", from_key="nacht-s16")
    mid = sent["id"]
    M.take("nacht-s17")
    assert M.status(mid) == "delivered"


def test_status_unknown_for_a_missing_id(tmp_path):
    assert M.status("definitely-not-a-real-id") == "unknown"


def test_ls_lists_sent_and_delivered_messages(tmp_path):
    M.send("nacht-s17", "first", from_key="nacht-s16")
    M.send("nacht-s18", "second", from_key="nacht-s16")
    rows = M.ls()
    assert {r["to"] for r in rows} == {"nacht-s17", "nacht-s18"}
    assert {r["status"] for r in rows} == {"queued"}
    M.take("nacht-s17")
    by_to = {r["id"]: r["status"] for r in M.ls()}
    delivered_ids = [i for i, s in by_to.items()
                     if M.status(i) == "delivered"]
    assert delivered_ids, "delivered message is still listed, with a receipt"


def test_a_typo_inbox_ages_out_of_known(tmp_path):
    """One message to a typo must not make it permanently deliverable: after
    its inbox is idle past the freshness window it is unknown again, so a
    later message is refused instead of silently accepted forever (reviewer
    finding: no-TTL latch re-creates 'sender never knows if read'). The
    window is session_presence._STALE_S (QS decision)."""
    M.send("operatro", "typo", from_key="nacht-s16")
    assert M.deliverable("operatro") is True   # freshly written inbox
    inbox = M._inbox("operatro")
    now = time.time()
    import os
    os.utime(inbox, (now - P._STALE_S - 60, now - P._STALE_S - 60))
    assert M.deliverable("operatro") is False   # aged past the window: unknown


def test_a_fresh_inbox_is_still_known(tmp_path):
    """The freshness TTL does not break the phase-2 promise: a session that
    has received recently (its inbox is fresh) stays deliverable even before
    any presence record exists."""
    M.send("nacht-s13", "hello", from_key="nacht-s16")
    assert M.deliverable("nacht-s13") is True
