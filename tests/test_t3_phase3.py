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
    assert M.status(mid, "nacht-s16") == "queued"
    (msg,) = M.take("nacht-s17")
    assert msg["id"] == mid
    assert M.status(mid, "nacht-s16") == "delivered"
    M.mark_read(mid)
    assert M.status(mid, "nacht-s16") == "read"


def test_status_survives_the_inbox_being_taken(tmp_path):
    """The receipt outlives take(): take() unlinks the inbox, yet the sender
    can still ask status once the message has been delivered/read."""
    sent = M.send("nacht-s17", "heads-up", from_key="nacht-s16")
    mid = sent["id"]
    M.take("nacht-s17")
    assert M.status(mid, "nacht-s16") == "delivered"


def test_status_unknown_for_a_missing_id(tmp_path):
    assert M.status("definitely-not-a-real-id", "nacht-s16") == "unknown"


def test_ls_lists_sent_and_delivered_messages(tmp_path):
    M.send("nacht-s17", "first", from_key="nacht-s16")
    M.send("nacht-s18", "second", from_key="nacht-s16")
    rows = M.ls("nacht-s16")
    assert {r["to"] for r in rows} == {"nacht-s17", "nacht-s18"}
    assert {r["status"] for r in rows} == {"queued"}
    M.take("nacht-s17")
    by_id = {r["id"]: r["status"] for r in M.ls("nacht-s16")}
    delivered_ids = [i for i, s in by_id.items()
                     if M.status(i, "nacht-s16") == "delivered"]
    assert delivered_ids, "delivered message is still listed, with a receipt"


def test_ls_without_a_sender_sees_nothing(tmp_path):
    """A caller that names no key sees nothing -- status and ls never leak a
    message's existence to a caller who is not its sender."""
    M.send("nacht-s17", "secret", from_key="nacht-s16")
    mid = M.ls("nacht-s16")[0]["id"]
    assert M.ls() == []
    assert M.status(mid) == "unknown"
    assert M.status(mid, "") == "unknown"


def test_a_session_cannot_see_another_sessions_or_the_operators_mail(tmp_path):
    """Security: session A must not learn about B's or the operator's
    messages -- status(scoped to sender) returns unknown and ls(scoped to
    sender) omits them, and no receipt record ever persists a message body."""
    b_msg = M.send("nacht-s17", "B's secret to 17", from_key="nacht-s18")
    op_msg = M.send("operator", "secret to the operator", from_key="nacht-s18")
    # A asks status for B's and the operator's message ids:
    assert M.status(b_msg["id"], "nacht-s16") == "unknown"
    assert M.status(op_msg["id"], "nacht-s16") == "unknown"
    assert M.status(op_msg["id"], "nacht-s18") in ("queued", "delivered")
    # A's ls sees none of them:
    assert [r["id"] for r in M.ls("nacht-s16")] == []
    # B's own ls sees B's sent messages only (both messages B sent):
    assert {r["id"] for r in M.ls("nacht-s18")} == {b_msg["id"], op_msg["id"]}
    # And once delivered, no receipt ever stores the message body:
    M.take("nacht-s17")
    receipts = M._receipts()
    if b_msg["id"] in receipts:
        assert "text" not in receipts[b_msg["id"]]
        assert "from_title" not in receipts[b_msg["id"]]


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


def test_an_overflow_dropped_message_stays_accountable(tmp_path):
    """Reviewer finding (s17, absorbed from the withdrawn probe): when a
    recipient's inbox holds more than _MAX_TAKE messages, _take keeps the
    newest _MAX_TAKE and discards the older ones -- and, on the old core,
    discards them WITHOUT a receipt. A sender's status() for a dropped message
    then regressed from 'queued' to 'unknown', the same value as a never-sent
    id, and ls() lost it: the sender could no longer tell 'my message was
    dropped by the cap' from 'I mistyped the id'. That is the silent-loss
    hazard phase 3 targets. The cap must never make a queued message
    indistinguishable from one that never existed."""
    cap = M._MAX_TAKE
    first_id = None
    for i in range(cap + 2):  # two messages beyond what one take will deliver
        sent = M.send("nacht-s17", f"message {i}", from_key="nacht-s16")
        if i == 0:
            first_id = sent["id"]
    # First message was confirmed queued before the capped take.
    assert M.status(first_id, "nacht-s16") == "queued"
    got = M.take("nacht-s17")
    assert len(got) == cap + 1  # cap delivered + one truncation marker
    # The sender must still be able to account for the dropped message.
    assert M.status(first_id, "nacht-s16") != "unknown", (
        "a message that was queued must not become indistinguishable from a "
        "never-sent id when the inbox cap drops it"
    )
    assert any(r["id"] == first_id for r in M.ls("nacht-s16")), (
        "ls() must still show the dropped message so the sender can see it "
        "was never delivered"
    )


def test_an_overflow_dropped_message_reports_dropped(tmp_path):
    """The status a dropped message reports must name the drop honestly --
    'dropped', not the 'delivered' that a receipt-with-no-read_at would give.
    It was accepted into the inbox but never handed to the recipient's prompt,
    so it is not delivered; calling it delivered would lie to the sender."""
    cap = M._MAX_TAKE
    first_id = None
    for i in range(cap + 2):
        sent = M.send("nacht-s18", f"message {i}", from_key="nacht-s16")
        if i == 0:
            first_id = sent["id"]
    M.take("nacht-s18")
    assert M.status(first_id, "nacht-s16") == "dropped"
    row = next(r for r in M.ls("nacht-s16") if r["id"] == first_id)
    assert row["status"] == "dropped"


def test_status_will_still_be_delivered_for_a_cap_survivor(tmp_path):
    """The cap must not corrupt the receipts of the messages it DOES deliver:
    the newest _MAX_TAKE that survive a capped take still report delivered."""
    cap = M._MAX_TAKE
    survivor_id = None
    for i in range(cap + 1):  # the newest 'cap' survive; only the oldest drops
        sent = M.send("nacht-s17", f"message {i}", from_key="nacht-s16")
        survivor_id = sent["id"]
    M.take("nacht-s17")
    assert M.status(survivor_id, "nacht-s16") == "delivered"


def test_operator_delivery_is_traceable_by_a_normal_sender(tmp_path):
    """A normal session can ask status of a message it sent to the operator
    mailbox -- the operator takes it like any inbox, and the sender sees it
    delivered/read. The operator path must not be a special case that the
    sender cannot track."""
    sent = M.send("operator", "status check", from_key="nacht-s16")
    mid = sent["id"]
    assert M.status(mid, "nacht-s16") == "queued"
    M.take("operator")
    assert M.status(mid, "nacht-s16") == "delivered"
    M.mark_read(mid)
    assert M.status(mid, "nacht-s16") == "read"


def test_the_receipts_sidecar_is_never_scanned_as_an_inbox(tmp_path):
    """The receipts file (.receipts.json) lives inside the inbox dir but must
    never be handed to a session as a message or listed as a recipient. It is
    structurally excluded because every inbox scan globs *.jsonl and the
    sidecar is .receipts.json -- pinned so a future glob widening cannot
    silently turn the sidecar into deliverable mail or a listed row."""
    M.send("nacht-s17", "hello", from_key="nacht-s16")
    M.take("nacht-s17")  # delivers -> writes the .receipts.json sidecar
    sidecar = M._receipts_path()
    assert sidecar.suffix == ".json", "the sidecar must not look like an inbox"
    # Every scan source in the module globs only *.jsonl, and the sidecar is
    # never among them -- so no inbox scan (take/ls/_queued_sender_of) reads it.
    scanned = M._DIR.glob("*.jsonl")
    assert sidecar not in scanned
    # The genuinely delivered message is still reported, nothing read off the
    # sidecar file name.
    rows = M.ls("nacht-s16")
    assert len(rows) == 1 and rows[0]["to"] == "nacht-s17"
