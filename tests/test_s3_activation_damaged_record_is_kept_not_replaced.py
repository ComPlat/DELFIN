"""A damaged session record is kept, not silently rewritten, and the phone
must not re-enable a session whose record it cannot read.

Measured here on the activation record (~/.delfin/messenger_sessions/<hash>.json):
the writer ``activation.configure`` reads the record with ``record()``, whose
``except (OSError, ValueError, AttributeError): return {}`` conflates "nothing
is recorded" with "what is recorded cannot be read". A file truncated the way
a kill mid-write truncates it then reads as {} -- "the session is not switched
off in the session" -- and:

* ``switched_off_in_session`` returns False, so the phone's ``hol`` re-enables
  a session the person deliberately switched off in the session itself
  (the invariant "switched off in the session cannot be switched on from the
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
    activation.configure("s1", "on")
    assert activation.switched_off_in_session("s1")
    # Damage the record the way a crash mid-write does.
    store = _record_file()
    store.write_text(TRUNCATED)
    # The guard must not read damage as "no restriction recorded".
    assert activation.switched_off_in_session("s1") is True, (
        "an unreadable record must not grant the phone a re-enable")


def test_configure_keeps_a_damaged_record_aside_never_overwriting_it(_records):
    """Writing a new state after a damaged record must not destroy the old
    bytes: they are quarantined, and the write no longer replaces the file
    in place."""
    activation.configure("s1", "on")
    store = _record_file()
    original = store.read_bytes()
    store.write_text(TRUNCATED)
    activation.configure("s1", "off")
    # The old bytes must survive somewhere -- never silently overwritten.
    assert any(p.read_bytes() == original for p in store.parent.glob("*.unreadable-*")), \
        "the damaged record must be kept, not overwritten in place"
