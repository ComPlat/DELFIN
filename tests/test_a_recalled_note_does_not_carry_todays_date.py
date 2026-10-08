"""The header of a recalled memory is the same today as yesterday.

The prompt head is cached as a prefix: the first byte that differs from
the previous request makes everything after it cold. Between two builds
of the solo prompt one second apart, the ONLY difference was
"last recalled 2026-10-07" becoming "-08" inside a recalled note, at
95.6% of the way in -- so the tail of the system prompt and the whole
conversation after it were re-read at full price every day, for a date
the model does nothing with. The written date stays: a two-year-old note
should read as one.
"""

from __future__ import annotations

from delfin.agent.prompt_loader import _memory_entry_header


def _note(created: int, updated: int) -> str:
    return (f"---\nsource: user\ncreated_at: {created}\n"
            f"updated_at: {updated}\n---\nUser prefers short answers.\n")


def test_the_recalled_date_is_not_in_the_header():
    day = 86_400
    created = 1_700_000_000
    a = _memory_entry_header("note", "notes/a.md", _note(created, created + 30 * day))
    b = _memory_entry_header("note", "notes/a.md", _note(created, created + 31 * day))
    assert a == b, (a, b)
    assert "last recalled" not in a


def test_the_written_date_survives():
    import time
    created = 1_700_000_000
    head = _memory_entry_header("note", "notes/a.md", _note(created, created + 99))
    assert time.strftime("%Y-%m-%d", time.localtime(created)) in head
    assert "notes/a.md" in head and head.startswith("# note")
