"""T3 phase 3 — adversarial review of the delivery-receipt storage core.

Attacks session_messages.status()/ls()/mark_read()/_queued_sender_of() (the
merged, green phase-3 core that the pending ``delfin-agent messages`` CLI
patch burns into) beyond the builder's own test_t3_phase3.py. Focus:

- the security contract the operator set (receipts carry NO text; status is
  scoped to the sender) must hold even in the corner cases;
- the rows ls() returns are exactly metadata, never a message body;
- a delivered/read receipt is stable: a second take, a failed-ish concurrent
  read, or a later status() query must not corrupt the answer;
- ls()/status() never merge the queued and delivered views of one message
  into a row that leaks state to a non-sender.

Everything here drives the committed storage functions directly (the CLI
surface is a protected .gate/t3_cli.patch handed to the operator), so these
are orthogonal to the builder's red CLI controls.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def test_ls_rows_are_metadata_only_never_a_message_body(tmp_path):
    """A listed row is {id,to,from,status,sent_at} and nothing else: the
    sender must never read its own text back through the summary the CLI
    prints, nor any title."""
    M.send("nacht-s17", "TOP-SECRET-BODY", from_key="nacht-s16",
           from_title="TITLE-LEAK")
    rows = M.ls("nacht-s16")
    assert rows, "at least one queued row"
    for row in rows:
        assert set(row) == {"id", "to", "from", "status", "sent_at"}
        assert "TOP-SECRET-BODY" not in json.dumps(row)
        assert "TITLE-LEAK" not in json.dumps(row)
    # after delivery the receipt row is equally body-free
    M.take("nacht-s17")
    for row in M.ls("nacht-s16"):
        assert set(row) == {"id", "to", "from", "status", "sent_at"}
        assert "TOP-SECRET-BODY" not in json.dumps(row)


def test_status_never_prints_the_body_for_any_state(tmp_path):
    """status() returns exactly one of queued / delivered / read / unknown and
    never the message text, whatever state the message is in."""
    sent = M.send("nacht-s17", "SECRET-BODY", from_key="nacht-s16")
    for _ in range(3):  # repeated queries, queued state
        assert M.status(sent["id"], "nacht-s16") in (
            "queued", "delivered", "read", "unknown")
    out = M.status(sent["id"], "nacht-s16")
    assert "SECRET-BODY" not in str(out)


def test_a_second_take_is_idempotent_for_the_receipt(tmp_path):
    """Taking an already-taken mailbox returns [] and must not corrupt the
    existing delivered/read receipt: status stays delivered."""
    sent = M.send("nacht-s17", "hello", from_key="nacht-s16")
    M.take("nacht-s17")          # delivered
    assert M.status(sent["id"], "nacht-s16") == "delivered"
    assert M.take("nacht-s17") == []   # nothing left, idempotent
    assert M.status(sent["id"], "nacht-s16") == "delivered"
    # and mark_read then survives the second take too
    M.mark_read(sent["id"])
    assert M.status(sent["id"], "nacht-s16") == "read"
    M.take("nacht-s17")
    assert M.status(sent["id"], "nacht-s16") == "read"


def test_status_does_not_merge_two_senders_for_one_recipient(tmp_path):
    """Two senders targeting the same recipient keep separate receipts: A
    never learns the id B sent, even though both sit in the same inbox."""
    a = M.send("nacht-s17", "A-letter", from_key="nacht-s16")
    b = M.send("nacht-s17", "B-letter", from_key="nacht-s18")
    # A cannot get B's id from ls() -- B's queued row is not listed for A
    assert a["id"] in [r["id"] for r in M.ls("nacht-s16")]
    assert b["id"] not in [r["id"] for r in M.ls("nacht-s16")]
    M.take("nacht-s17")
    assert M.status(a["id"], "nacht-s16") == "delivered"
    assert M.status(a["id"], "nacht-s18") == "unknown"   # not B's message
    assert M.status(b["id"], "nacht-s16") == "unknown"   # not A's message
    assert M.status(b["id"], "nacht-s18") == "delivered"


def test_queued_message_is_unknown_to_an_empty_and_foreign_sender(tmp_path):
    """A queued message's id is not a handle for anyone but its sender: status
    with no sender, or the wrong sender, is 'unknown' -- never 'queued'."""
    sent = M.send("nacht-s17", "x", from_key="nacht-s16")
    assert M.status(sent["id"]) == "unknown"
    assert M.status(sent["id"], "") == "unknown"
    assert M.status(sent["id"], "nacht-s99") == "unknown"
    assert M.status(sent["id"], "nacht-s16") == "queued"


def test_mark_read_only_reads_an_existing_delivered_receipt(tmp_path):
    """mark_read on an unknown id is a silent no-op (no crash, no fabricated
    receipt), and never reads a message that is still queued in the inbox."""
    ghost = "0" * 32
    M.mark_read(ghost)                 # never sent
    assert M.status(ghost, "nacht-s16") == "unknown"
    live = M.send("nacht-s17", "still-queued", from_key="nacht-s16")
    M.mark_read(live["id"])            # queued, not yet delivered
    assert M.status(live["id"], "nacht-s16") == "queued"
    assert M._receipts().get(live["id"]) is None
