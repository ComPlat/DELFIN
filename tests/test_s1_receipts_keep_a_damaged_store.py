"""The session-message receipts store is kept, not silently re-wiped.

``_receipts()`` answered "absent" and "unreadable" with the same empty dict,
and the three writers that follow it (``_receipts_deliver``,
``_receipts_dropped``, ``mark_read``) save the whole store back with
``O_TRUNC``. A receipts file torn the way a full disk or a kill mid-write
tears it made ``_receipts()`` return ``{}``, the writer wrote ``{}`` plus the
new delivery back over it, and every earlier delivery/read record was gone
--- no message, no copy. Delivery receipts are not reconstructable (they
record events that passed), which is why they are persisted at all.

The read used to DISPLAY (``status`` / ``ls``) stays lenient on purpose; only
the writers go strict here.
"""

from __future__ import annotations

from delfin.agent import session_messages as SM

_MSG = {
    "to": "k", "from": "alice", "from_title": "alice", "text": "hi",
    "sent_at": 1.0, "id": "aaaa1111",
}


def _inbox(tmp_path, monkeypatch) -> "object":
    inbox = tmp_path / "inbox"
    inbox.mkdir()
    monkeypatch.setattr(SM, "_DIR", inbox)
    return inbox


def _aside(store):
    return sorted(q for q in store.parent.iterdir()
                  if q.name.startswith(store.name + ".unreadable-"))


def test_delivering_keeps_a_torn_receipts_store_and_reports(
        tmp_path, monkeypatch, caplog):
    import logging
    inbox = _inbox(tmp_path, monkeypatch)
    caplog.set_level(logging.WARNING, logger="delfin.agent.session_messages")

    # Real receipts, then the failure that truncates them.
    SM._receipts_deliver([{**_MSG, "id": "aaaa1111"}])
    SM._receipts_deliver([{**_MSG, "id": "bbbb2222"}])
    store = inbox / ".receipts.json"
    original = store.read_text()
    assert "aaaa1111" in original
    store.write_text(original[: len(original) // 2])

    # The writer runs again (a later delivery).
    SM._receipts_deliver([{**_MSG, "id": "cccc3333"}])

    # The damaged original is not silently gone: it sits aside.
    kept = _aside(store)
    assert kept, (
        "a receipts store that could not be read must be moved aside, "
        "not written over and lost")
    assert "aaaa1111" in kept[0].read_text(), (
        "the delivery/read records the damaged store still held must be "
        "recoverable from the aside")

    # The loss is reported, not swallowed.
    assert any("receipts" in r.message for r in caplog.records), (
        "the writer must say that a damaged receipts store was kept aside")

    # The new receipt is still recorded --- the write is not refused.
    assert SM._receipts().get("cccc3333"), "the fresh delivery was stored"


def test_an_intact_receipts_store_is_untouched(tmp_path, monkeypatch):
    """The control: healthy stores save exactly as before, no aside."""
    inbox = _inbox(tmp_path, monkeypatch)
    SM._receipts_deliver([{**_MSG, "id": "aaaa1111"}])
    store = inbox / ".receipts.json"
    assert not _aside(store)
    assert "aaaa1111" in SM._receipts()


def test_a_reader_still_answers_an_empty_store(tmp_path, monkeypatch):
    """Lenient on purpose: a damaged store reads as nothing recorded."""
    inbox = _inbox(tmp_path, monkeypatch)
    store = inbox / ".receipts.json"
    store.write_text('{"aaaa1111": {"to": "k"')   # torn
    assert SM._receipts() == {}
    assert _aside(store) == [], "a reader does not move the file"
