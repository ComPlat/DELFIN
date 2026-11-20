"""A damaged observed-window store is quarantined, never emptied.

The observed-window store (``delfin.agent.observed_window``) persists the
largest input the backend accepted per model. ``note`` is the writer: it
read the store, raised one model's floor, and wrote the whole file back.
Its read answered "nothing was learned" for both "the file is absent" and
"the file cannot be read", so a store truncated the way a full disk or a
kill mid-write truncates it became, on the next ``note``, a file holding
only the one model just written. The floors every other run had measured
were gone, with no copy anywhere.

The fix routes the writer's read through ``agent.state_file``'s
``read_json_for_write``: a damaged store is moved to a dated
``.unreadable-`` quarantine (bytes and mode kept) and the write is
REFUSED -- a store that cannot be read cannot be safely replaced, so
nothing lands on top of it. The store is gone, so the next note writes a
fresh one safely, and whatever the damaged file still held is recoverable
from the quarantine; the loss is logged.
"""

from __future__ import annotations

import json

import delfin.agent.observed_window as ow


def _store_path():
    # The suite redirects this sink into the test's own directory.
    p = ow._PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _aside(path):
    return sorted(q for q in path.parent.iterdir()
                  if "unreadable" in q.name)


def test_note_refuses_to_replace_a_damaged_store_and_keeps_its_bytes():
    p = _store_path()
    p.write_text(json.dumps({"model-a": 100, "model-b": 200}), encoding="utf-8")
    # Truncate the way a full disk or a kill mid-write truncates it.
    damaged = '{"model-a": 100, "model-b": 2'
    p.write_text(damaged, encoding="utf-8")

    ow.note("model-a", 500)

    # The write is refused: the damaged bytes are not rewritten, truncated
    # further, or silently replaced by a store that drops the other floors.
    moved = _aside(p)
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")
    assert not p.exists(), (
        "nothing may be written over a store that could not be read")
    assert ow.floor("model-a") == 0, "the damaged floors are not read back"

    # The store is gone, so the next note writes a fresh one safely.
    ow.note("model-a", 500)
    assert json.loads(p.read_text(encoding="utf-8")) == {"model-a": 500}
    assert len(_aside(p)) == 1, "the second note quarantines nothing"


def test_note_does_not_quarantine_an_intact_store():
    """The control: nothing above may change the normal path."""
    p = _store_path()
    p.write_text(json.dumps({"model-a": 100, "model-b": 200}), encoding="utf-8")

    ow.note("model-a", 500)

    live = json.loads(p.read_text(encoding="utf-8"))
    assert live["model-a"] == 500
    assert live["model-b"] == 200, "an ordinary write keeps the other floors"
    assert not _aside(p), "nothing is quarantined when the store parses"


def test_a_reader_tolerates_a_damaged_store_without_raising():
    p = _store_path()
    p.write_text('{"model-a": 100, "model-b": ', encoding="utf-8")
    # A lookup must not stop a turn on a broken file.
    assert ow.floor("model-a") == 0
    assert ow.at_least("model-a", 4096) == 4096
