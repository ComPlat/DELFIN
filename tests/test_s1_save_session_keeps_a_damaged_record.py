"""A session record damaged on disk is kept, not silently overwritten.

The writer path of ``session_store.save_session`` reads the existing
``<id>.json`` only to preserve two fields (``created_at`` and the pending
plan's "already announced" mark), then writes the WHOLE file back with the
live turn's data. A single damaged byte in that record made the reader throw,
the two fields were rebuilt as if the session were new, and the full
conversation history was then silently overwritten --- no message, no copy.

Measured on 2026-11-19 by driving the real writer against a torn file: the
``.unreadable-*`` quarantine that ``state_file`` exists to produce was never
created, because the reader answered "absent" and "damaged" with the same
empty result.

The reader used by display (``load_session``) stays lenient on purpose, an
empty list is the honest answer for a reader. Only the writer is strict here.
"""

from __future__ import annotations

import delfin.agent.session_store as SS


def _torn(record: str) -> str:
    """Cut the file the way a kill mid-write does, keeping the first half."""
    return record[: len(record) // 2]


def test_save_session_keeps_a_torn_record_and_reports(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(SS, "_SESSIONS_DIR", tmp_path)
    monkeypatch.setattr(SS, "_startup_maintenance", lambda: None)

    import logging
    caplog.set_level(logging.WARNING, logger="delfin.agent.session_store")

    # A real session, then the failure that truncates it.
    SS.save_session(
        "s1",
        chat_messages=[{"role": "user", "content": "the earlier conversation"}],
        title="the work",
    )
    stored = tmp_path / "s1.json"
    original = stored.read_text()
    stored.write_text(_torn(original))

    # The writer runs again (the next turn's save).
    SS.save_session(
        "s1",
        chat_messages=[{"role": "user", "content": "the later turn"}],
        title="the work",
    )

    # The damaged original is not silently gone: it sits aside, byte for byte.
    kept = [q for q in tmp_path.iterdir() if q.name.startswith("s1.json.unreadable-")]
    assert kept, (
        "a session record that could not be read must be moved aside, "
        "not replaced out from under its history"
    )
    assert "the earlier conversation" in kept[0].read_text(), (
        "what the damaged record still held must be recoverable from the aside")

    # The loss is reported, not swallowed: the warning names the file.
    assert any("s1.json.unreadable-" in r.message for r in caplog.records), (
        "the save must say that a damaged record was kept aside, "
        "so the person can act on it")

    # The live session still carries the new turn --- a save is not refused.
    live = SS.load_session("s1")
    assert live is not None
    assert any(m.get("content") == "the later turn"
               for m in live.get("chat_messages", []))


def test_an_intact_record_is_untouched_by_all_of_this(tmp_path, monkeypatch):
    """The control: a healthy session saves exactly as before, no aside."""
    monkeypatch.setattr(SS, "_SESSIONS_DIR", tmp_path)
    monkeypatch.setattr(SS, "_startup_maintenance", lambda: None)

    SS.save_session(
        "c1",
        chat_messages=[{"role": "user", "content": "a turn"}],
        title="control",
    )
    assert not [q for q in tmp_path.iterdir()
                if q.name.startswith("c1.json.unreadable-")]
    data = SS.load_session("c1")
    assert data.get("title") == "control"
