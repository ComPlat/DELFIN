"""A damaged occupier auto-state store is quarantined, never emptied by a write.

``delfin.occupier_auto`` remembers which configuration won for a given
delta (``.delfin_occ_auto_state.json`` in the workspace root); the winner
drives the next redox step. ``record_auto_preference`` is the one writer:
it read the store, recorded one delta's preference, and wrote the whole
file back; the shared read answered "nothing is recorded" for both an
absent file and one that cannot be read, so a store truncated the way a
full disk or a kill mid-write truncates it became, on the next record, a
store holding only the one delta just written -- every other delta's
preferred index, multiplicity and basis set gone.

The fix splits the two reads. Readers keep the lenient answer -- a
damaged store must not crash the step that consumes it. The writer goes
through ``_load_state_for_write``, which on damage moves the file to the
dated ``.unreadable-`` quarantine of ``agent.state_file`` (bytes kept)
and raises ``StateFileUnreadable``: a store that cannot be read cannot be
safely replaced, so the write is refused and the reason carried to the
caller (the report steps log it). The store is gone, so the next record
writes a fresh one safely.

The store is an explicit ``root=`` argument into the test's tmp dir;
nothing here ever touches the real workspace.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.state_file import StateFileUnreadable
from delfin.occupier_auto import STATE_FILENAME, record_auto_preference


def _path(root):
    return root / STATE_FILENAME


def _aside(root):
    return sorted(p for p in root.iterdir() if "unreadable" in p.name)


def test_record_refuses_to_replace_a_damaged_store_and_keeps_its_bytes(tmp_path):
    # Two real preferences (two deltas), then the file truncated the way a
    # full disk or a kill mid-write truncates it.
    record_auto_preference("even", 3, 1, m_value=1, bs_value="def2-SVP", root=tmp_path)
    record_auto_preference("odd", 1, 2, m_value=3, root=tmp_path)
    damaged = _path(tmp_path).read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    _path(tmp_path).write_text(damaged, encoding="utf-8")

    with pytest.raises(StateFileUnreadable):
        record_auto_preference("even", 0, 3, m_value=1, root=tmp_path)

    # The write is refused: the damaged bytes are not rewritten, truncated
    # further, or silently replaced by a store that drops the other deltas.
    moved = _aside(tmp_path)
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")
    assert not _path(tmp_path).exists(), (
        "nothing may be written over a store that could not be read")

    # The store is gone, so the next record writes a fresh one safely.
    record_auto_preference("even", 0, 3, m_value=1, root=tmp_path)
    fresh = _path(tmp_path).read_text(encoding="utf-8")
    assert '"3"' in fresh, "the new store holds only what was just recorded"
    assert len(_aside(tmp_path)) == 1, "the second record quarantines nothing"


def test_an_intact_store_is_untouched_by_any_of_this(tmp_path):
    """The control: nothing above may change the normal path."""
    record_auto_preference("even", 3, 1, m_value=1, bs_value="def2-SVP", root=tmp_path)
    record_auto_preference("odd", 1, 2, m_value=3, root=tmp_path)
    # Overwrite the even preference for delta 1: the delta keeps one entry.
    record_auto_preference("even", 5, 1, m_value=1, root=tmp_path)
    data = json.loads(_path(tmp_path).read_text(encoding="utf-8"))
    assert set(data) == {"1", "2"}, "both deltas are recorded"
    assert data["1"]["even"]["index"] == 5, "delta 1 holds the latest even preference"
    assert data["2"]["odd"]["index"] == 1, "delta 2 keeps its own preference"
    assert not _aside(tmp_path), "nothing is quarantined when the store parses"
