"""A damaged MCP-pins record is kept aside and the write is refused, never
silently retired by the next pin write.

Measured here on the trust-store record ``mcp_pins.json``, which holds the
hash each MCP tool's description had when it was first seen (a change from
that pin is flagged, never silently served). The pins are written by
``mcp_client._pins_update``, which builds the next record on
``trust_store.read_record`` -- ``except ValueError: return None`` conflates
"no tool was ever pinned" with "what was pinned cannot be read". A file
truncated the way a kill mid-write truncates it then reads as None -- "no
restriction is recorded" -- and the write replaces the record, retiring every
previously-pinned hash with no copy anywhere.

This is the dangerous direction: a tool whose description changed from its
pin would present as 'first seen' on the next write instead of flagged.
The fix keeps the writer from silently replacing a store it could not read:
the write is refused and the damaged bytes are moved aside for the person.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import trust_store as T

MARKER = "DELFIN_AGENT_SESSION"


@pytest.fixture(autouse=True)
def _store(tmp_path, monkeypatch):
    store = tmp_path / "store"
    monkeypatch.setattr(T, "store_dir", lambda: store)
    monkeypatch.delenv(MARKER, raising=False)
    return store


def _seed():
    return {"servers": {"srv#1": {
        "name": "srv", "first_seen": 1000.0,
        "tools": {"tool1": "deadbeef" * 8}, "pending": {}}}}


def _change(data):
    data["servers"].setdefault("srv#2", {"name": "srv2", "first_seen": 2000.0,
                                         "tools": {}, "pending": {}})


def test_a_damaged_pins_record_is_refused_and_kept_aside(_store):
    """The next pin write must not silently replace a record it could not
    read: the write is refused and the damaged bytes are kept aside."""
    M._PIN_CACHE.clear()
    M._PIN_MEMORY.clear()
    M._PIN_ANSWERS.clear()
    M._PIN_NOTED.clear()
    store = _store / M.PIN_RECORD
    store.parent.mkdir(parents=True, exist_ok=True)
    store.write_text(json.dumps(_seed()))
    store.write_text('{"servers": {"srv#1"'[:32])  # truncate mid-document
    before = store.read_bytes()
    # Drive the real writer on the damaged record.
    ok = M._pins_update(_change)
    # The write must not land on top of a store the writer could not read.
    assert ok is False, (
        "a pin write that cannot read the record must be refused, not "
        "repair it into an empty store")
    # The damaged bytes are kept aside -- never overwritten in place.
    aside = [p for p in _store.glob(M.PIN_RECORD + ".unreadable-*")]
    assert aside, "the damaged pins record must be kept aside, not overwritten"
    assert any(p.read_bytes() == before for p in aside), (
        "the kept bytes are what the damaged record held")
