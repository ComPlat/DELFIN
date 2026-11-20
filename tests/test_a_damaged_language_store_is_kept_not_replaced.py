"""A damaged contact-language store is quarantined, never emptied by a write.

``delfin.messenger.language`` keeps one entry per contact (their language,
mode, streak) in ``languages.json`` in the messenger state dir. The
writers ``set_contact`` and ``observe`` read the store, changed one
contact, and wrote the whole file back; the shared read answered "nothing
is known" for both an absent file and one that cannot be read, so a store
truncated the way a full disk or a kill mid-write truncates it became, on
the next message or the next ``language en`` command, a store holding
only the one contact just touched -- every other contact's language gone.

The fix splits the two reads. Readers (``of_contact``/``of_person``) keep
the lenient answer -- a damaged store must not crash message handling.
Writers go through ``_load_for_write``, which on damage moves the file to
the dated ``.unreadable-`` quarantine of ``agent.state_file`` (bytes
kept) and raises ``StateFileUnreadable``: a store that cannot be read
cannot be safely replaced, so the write is refused and the reason carried
to the caller (the daemon loop answers the person and retries next
round). The store is gone, so the next message writes a fresh one safely.

The store is redirected with ``DELFIN_MESSENGER_DIR`` to the test's temp
dir; nothing here ever touches the real state.
"""

from __future__ import annotations

import pytest

from delfin.agent.state_file import StateFileUnreadable
from delfin.messenger import language


@pytest.fixture()
def store_dir(monkeypatch, tmp_path):
    target = tmp_path / "messenger"
    target.mkdir()
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(target))
    return target


def _file(store_dir):
    return store_dir / "languages.json"


def _aside(store_dir):
    return sorted(p for p in store_dir.iterdir()
                  if "unreadable" in p.name)


def test_set_contact_refuses_to_replace_a_damaged_store_and_keeps_its_bytes(store_dir):
    language.set_contact("c1", "de")
    language.set_contact("c2", "en")
    damaged = _file(store_dir).read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    _file(store_dir).write_text(damaged, encoding="utf-8")

    # Readers stay lenient while the damaged file is still present: a
    # broken store does not crash message handling, it answers the
    # default rather than hiding the failure.
    assert language.of_contact("c1") in language.LANGS

    with pytest.raises(StateFileUnreadable):
        language.set_contact("c3", "en")

    moved = _aside(store_dir)
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")
    assert not _file(store_dir).exists(), (
        "nothing may be written over a store that could not be read")

    # The store is gone, so the next write starts a fresh one safely.
    language.set_contact("c3", "en")
    assert language.of_contact("c3") == "en"
    assert len(_aside(store_dir)) == 1, "the second write quarantines nothing"


def test_observe_refuses_to_write_over_a_damaged_store(store_dir):
    language.set_contact("c1", "en")
    damaged = _file(store_dir).read_text(encoding="utf-8")[:30]
    _file(store_dir).write_text(damaged, encoding="utf-8")

    with pytest.raises(StateFileUnreadable):
        language.observe("c1", "hallo")

    assert len(_aside(store_dir)) == 1
    assert _aside(store_dir)[0].read_text(encoding="utf-8") == damaged
    assert not _file(store_dir).exists()


def test_an_intact_store_is_untouched_by_any_of_this(store_dir):
    """The control: nothing above may change the normal path."""
    language.set_contact("c1", "de")
    language.set_contact("c2", "en")
    assert language.observe("c2", "hello") in language.LANGS
    assert language.of_contact("c1") == "de"
    assert language.of_contact("c2") == "en"
    assert not _aside(store_dir), "nothing is quarantined when the store parses"
