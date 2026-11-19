"""A damaged session record is kept, not silently rewritten, and the phone
must not re-enable a session whose record it cannot read.

Measured here on the activation record (one JSON file per session under the
phone-messenger per-user state): the writer ``activation.configure`` reads the
record with ``record()``, whose ``except (OSError, ValueError, AttributeError):
return {}`` conflates "nothing is recorded" with "what is recorded cannot be
read". A file truncated the way a kill mid-write truncates it then reads as {}
-- "the session is not switched off in the session" -- and:

* ``switched_off_in_session`` returns False, so the phone's pull re-enables a
  session the person deliberately switched off in the session itself (the
  invariant "switched off in the session cannot be switched on from the
  phone" is silently bypassed by corruption);
* ``configure`` builds the new record from that empty read and writes it in
  place, destroying the person's ``by`` / ``generation`` decision with no
  copy anywhere.

Both are the "no restriction is recorded" direction: the file existed to
withhold, and the empty read grants.
"""
from __future__ import annotations

import pytest

from delfin.messenger import activation


TRUNCATED = '{"session_id": "s1", "enabled": true, "by": "session", "generation"'


@pytest.fixture(autouse=True)
def _records(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "act")
    return tmp_path / "act"


def _record_file(sid="s1"):
    return activation._path(sid)


def test_the_phone_cannot_re_enable_a_session_whose_record_is_damaged(_records):
    """A session switched off in the session stays off from the phone even
    when its record file stops parsing."""
    activation.configure("s1", "off")
    # The person switched off IN the session: the phone may not undo it.
    assert activation.switched_off_in_session("s1") is True
    # Damage the record the way a crash mid-write does.
    _record_file().write_text(TRUNCATED)
    # The guard must not read damage as "no restriction recorded": the file
    # may exist exactly to withhold the phone, so it fails closed.
    assert activation.switched_off_in_session("s1") is True, (
        "an unreadable record must not grant the phone a re-enable")


def test_configure_keeps_a_damaged_record_aside_never_overwriting_it(_records):
    """Writing a new state after a damaged record must not destroy the old
    bytes: they are quarantined, and the write no longer replaces the file
    in place."""
    activation.configure("s1", "on")
    store = _record_file()
    store.write_text(TRUNCATED)
    activation.configure("s1", "off")
    # The damaged bytes must survive somewhere -- never silently overwritten
    # in place. They are moved aside, not lost.
    aside = [p for p in store.parent.glob("*.unreadable-*")]
    assert aside, "the damaged record must be kept aside, not overwritten"
    assert any(p.read_bytes() == TRUNCATED.encode() for p in aside), (
        "the kept bytes are what the damaged file held")


def test_configure_refuses_the_phone_re_enable_on_an_unreadable_record(
        _records):
    """A phone turning the session ON must not overwrite a record it could
    not read: that record may exist to withhold it."""
    activation.configure("s1", "on")
    store = _record_file()
    store.write_text(TRUNCATED)
    out = activation.configure("s1", "on", by="phone", phone="+49123")
    assert "not switched on from the phone" in out, (
        "the phone re-enable is refused, and the person is told why")
    # The damaged file stays exactly where it is. Moving it aside on this
    # path would be the hole the next test measures.
    assert store.read_text() == TRUNCATED, (
        "the refused path must not touch the record at all")
    assert not list(store.parent.glob("*.unreadable-*"))
    assert activation.enabled("s1") is False


def test_the_refusal_holds_for_a_second_attempt_from_the_phone(_records):
    """A refusal that lasts one attempt is not a refusal.

    Driven on 2026-11-19 against the first version of this guard, which
    quarantined the damaged record before refusing the phone. A quarantine
    is a rename, so after it the live path does not exist -- and an absent
    record reads as "nothing was recorded", which is the one answer this
    guard exists never to give. The phone's next pull, a second later, was
    granted and wrote an enabled record:

        guard before: True -> attempt 1 refused -> guard after: False
        -> attempt 2: "switched on" -> enabled: True

    A phone that is refused retries. So the damaged file is left in place:
    it is the only thing that still withholds, nothing the phone does can
    clear it, and the person at the computer resolves it (their own write
    quarantines it and stores the new state).
    """
    activation.configure("s1", "off")
    store = _record_file()
    store.write_text(TRUNCATED)

    for attempt in (1, 2, 3):
        out = activation.configure("s1", "on", by="phone", phone="+49123")
        assert "not switched on from the phone" in out, (attempt, out)
        assert activation.switched_off_in_session("s1") is True, attempt
        assert activation.enabled("s1") is False, attempt


def test_the_person_at_the_computer_is_still_the_way_back(_records):
    """Failing closed must not lock the session out for good: a write from
    the session itself quarantines the damaged record and stores the new
    state, so there is exactly one route back and it needs a person."""
    activation.configure("s1", "off")
    store = _record_file()
    store.write_text(TRUNCATED)

    out = activation.configure("s1", "on", endpoint="http://127.0.0.1:1")

    assert activation.enabled("s1") is True, out
    aside = list(store.parent.glob("*.unreadable-*"))
    assert any(p.read_bytes() == TRUNCATED.encode() for p in aside), (
        "the damaged bytes are kept when the person's write replaces them")
    assert "moved" in out and "aside" in out, (
        "and the person is told, in the answer they get: " + out)


def test_the_phone_may_still_switch_off_an_unreadable_record(_records):
    """Switching OFF is withheld by nothing, so it is not refused. The
    opposite would let a damaged file keep a session reachable from a phone
    the person wants to stop."""
    activation.configure("s1", "on")
    _record_file().write_text(TRUNCATED)

    out = activation.configure("s1", "off", by="phone")

    assert activation.enabled("s1") is False, out
