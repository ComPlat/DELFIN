"""A damaged user-model store is quarantined, never silently wiped.

Measured in the same class as the memory store and credential store
(2026-11-17): one read answers "absent" and "unreadable" with the same
value, and a writer that builds the file from that answer deletes
everything else in it with no copy and no message.

``user_model.UserModel`` is the accumulated, non-reproducible picture of
the person the agent works with -- role, expertise, preferences, working
style, corrections, learned across sessions. ``_read`` (a reader) may
treat a damaged file as "nothing learned"; the WRITERS (``update_from_user_message``,
``forget``, ``forget_all``) may not build the new file from that answer,
because the old facts are then gone.

The repair follows ``delfin.agent.state_file.py`` and the memory store's
precedent:

* ``update_from_user_message`` stores its new entry, quarantines the
  damaged file (bytes kept, dated name reserved with O_EXCL) and raises
  an exception naming the moved file, so the loss is never silent.
* ``forget`` / ``forget_all`` raise ``StateFileUnreadable`` instead of
  writing a store built from nothing -- deleting entries of a model you
  could not read removes the entries you could have kept.
"""

from __future__ import annotations

import pytest

from delfin.agent import state_file as SF
from delfin.agent.user_model import UserModel

TRUNCATED = ('{"schema": 1, "fields": {"preferences": ['
             '{"value": "the user prefers minimal diffs", '
             '"source_role": "user"}')


def _aside(path):
    return sorted(q for q in path.parent.iterdir()
                  if "unreadable" in q.name)


def _seeded_model(tmp_path):
    """A real store holding two typed facts."""
    p = tmp_path / "user_model.json"
    m = UserModel(path=p)
    m.update_from_user_message(
        "the user prefers minimal diffs", date="2026-10-11")
    m.update_from_user_message(
        "the user prefers short answers", date="2026-10-11")
    assert len(m.entries()["preferences"]) == 2
    return p


def test_a_damaged_model_is_quarantined_not_wiped(tmp_path):
    """update_from_user_message keeps the old bytes and reports the loss."""
    p = _seeded_model(tmp_path)
    # Truncate the way a kill mid-write or a full disk truncates it.
    p.write_text(TRUNCATED, encoding="utf-8")

    m = UserModel(path=p)
    with pytest.raises(Exception) as exc:
        m.update_from_user_message(
            "the user prefers concise report headers", date="2026-10-12")

    # The damaged file is preserved aside, named as unreadable.
    asides = _aside(p)
    assert asides, "the damaged store must be preserved, not overwritten"
    kept = asides[0].read_text(encoding="utf-8")
    assert "the user prefers minimal diffs" in kept, (
        "the original facts must survive"
    )
    # The message names the file the person must look at.
    assert "user_model.json" in str(exc.value)

    # The new entry still landed.
    fresh = UserModel(path=p)
    assert any(e["value"] == "the user prefers concise report headers"
               for e in fresh.entries()["preferences"])


def test_forget_raises_rather_than_wiping_a_damaged_model(tmp_path):
    """forget must not build a store from nothing: it raises instead."""
    p = _seeded_model(tmp_path)
    p.write_text(TRUNCATED, encoding="utf-8")

    m = UserModel(path=p)
    with pytest.raises(SF.StateFileUnreadable):
        m.forget("preferences", 0)

    # Nothing was written over the damaged store; bytes still there.
    asides = _aside(p)
    assert not asides, "forget writes nothing, so nothing to quarantine"
    raw = p.read_text(encoding="utf-8")
    assert "prefers minimal diffs" in raw


def test_forget_all_raises_rather_than_wiping_a_damaged_model(tmp_path):
    """forget_all on an unreadable store must not empty the file on disk."""
    p = _seeded_model(tmp_path)
    p.write_text(TRUNCATED, encoding="utf-8")

    m = UserModel(path=p)
    with pytest.raises(SF.StateFileUnreadable):
        m.forget_all()

    assert "prefers minimal diffs" in p.read_text(encoding="utf-8")


def test_the_reader_stays_lenient_for_display(tmp_path):
    """entries()/select() (display) must not raise on a damaged model."""
    p = _seeded_model(tmp_path)
    p.write_text(TRUNCATED, encoding="utf-8")
    m = UserModel(path=p)
    # A dashboard panel loads the model for display; it must show empty,
    # not raise.
    assert m.entries()["preferences"] == []
