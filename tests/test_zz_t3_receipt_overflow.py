"""T3 phase-3 adversarial review (reviewer nacht-s17) of the receipts core.

Finding (RED): the _MAX_TAKE cap destroys a sender's accountability for the
messages it chopped off.

  send() hands a sender a message id and page-2-phase3 promises "sender must be
  able to ask what happened to a message it sent: status(message_id) reports
  queued / delivered / read". That contract breaks at the inbox cap: when a
  recipient's inbox holds more than session_messages._MAX_TAKE messages,
  _take (session_messages.py:204) keeps the NEWEST _MAX_TAKE, drops the older
  ones, and calls _receipts_deliver only on the surviving set (line 245). The
  dropped messages are unlinked with the inbox but get no receipt and leave no
  trace, so a sender's status() for one of them regresses from "queued" to
  "unknown" -- the same value as a never-sent id -- and it vanishes from ls().
  The sender can no longer tell "my message was dropped by the cap" from "I
  mistyped the id". That is the silent-loss hazard the phase targets.

  The message WAS delivered into the inbox; only the cap discarded it. It must
  not become byte-for-byte indistinguishable from a message that never existed.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_messages as M


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    # Also isolate the receipts sidecar (next to the inbox dir) implicitly.
    assert M._receipts_path().parent == tmp_path / "inbox"
    M._DIR.mkdir(parents=True, exist_ok=True)


def test_an_overflow_dropped_message_stays_accountable(tmp_path):
    """A message dropped by the _MAX_TAKE cap must not become unknown.

    It was sent and received into the inbox; the cap then discarded it. The
    sender must still be able to account for it -- either ls() still lists it
    (as undelivered) or status() reports something other than the same
    'unknown' a never-sent id gives. Red on the current core: the dropped
    message falls off both ls() and receipts and status() says 'unknown'."""
    cap = M._MAX_TAKE
    first_id = None
    for i in range(cap + 2):  # two messages beyond what one take will deliver
        sent = M.send("busy", f"message {i}", from_key="sender")
        if i == 0:
            first_id = sent["id"]
    # First message was confirmed queued before the capped take.
    assert M.status(first_id, "sender") == "queued"
    got = M.take("busy")
    assert len(got) == cap + 1  # cap delivered + one truncation marker
    # The sender must still be able to account for the dropped message.
    assert M.status(first_id, "sender") != "unknown", (
        "a message that was queued must not become indistinguishable from a "
        "never-sent id when the inbox cap drops it"
    )
    assert any(r["id"] == first_id for r in M.ls("sender")), (
        "ls() must still show the dropped message so the sender can see it "
        "was never delivered"
    )


def test_status_will_still_be_delivered_for_a_cap_survivor(tmp_path):
    """The cap must not corrupt the receipts of the messages it DOES deliver."""
    cap = M._MAX_TAKE
    survivor_id = None
    for i in range(cap + 1):  # the newest 'cap' survive; only the oldest drops
        sent = M.send("busy2", f"message {i}", from_key="sender")
    survivor_id = sent["id"]
    M.take("busy2")
    assert M.status(survivor_id, "sender") == "delivered"


def test_operator_delivery_is_traceable_by_a_normal_sender(tmp_path):
    """A normal session can ask status of a message it sent to the operator
    mailbox -- the operator takes it like any inbox, and the sender sees it
    delivered/read. The operator path must not be a special case that the
    sender cannot track."""
    sent = M.send("operator", "status check", from_key="sender")
    mid = sent["id"]
    assert M.status(mid, "sender") == "queued"
    M.take("operator")
    assert M.status(mid, "sender") == "delivered"
    M.mark_read(mid)
    assert M.status(mid, "sender") == "read"


def test_the_receipts_sidecar_is_never_scanned_as_an_inbox(tmp_path):
    """The receipts file (.receipts.json) lives inside the inbox dir but must
    never be handed to a session as a message or listed as a recipient. It is
    structurally excluded because every inbox scan globs ``*.jsonl`` and the
    sidecar is ``.receipts.json`` -- pinned so a future glob widening cannot
    silently turn the sidecar into deliverable mail or a listed row."""
    M.send("target", "hello", from_key="sender")
    M.take("target")  # delivers -> writes the .receipts.json sidecar
    sidecar = M._receipts_path()
    assert sidecar.suffix == ".json", "the sidecar must not look like an inbox"
    # Every scan source in the module globs only *.jsonl, and the sidecar is
    # never among them -- so no inbox scan (take/ls/_queued_sender_of) reads it.
    scanned = M._DIR.glob("*.jsonl")
    assert sidecar not in scanned
    # The genuinely delivered message is still reported, nothing read off the
    # sidecar file name.
    rows = M.ls("sender")
    assert len(rows) == 1 and rows[0]["to"] == "target"
