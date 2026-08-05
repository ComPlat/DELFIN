"""The undo the framework records is one only the model can reach.

change_journal.revert is a careful piece of work: it refuses when the
file's current hash no longer matches what the agent left, refuses an
oversize pre-image, verifies the stored pre-image's own hash before
restoring, preserves the target's permission bits, and unwinds chained
edits newest-first. It has exactly one caller in the whole tree — the
model's own `undo_changes` tool.

What the user can reach: `/undo`, which drops messages from context and
touches no file; and a hidden "Undo Edit" button whose handler shells out
to `git checkout` and, for a workspace that is not a git repository,
prints "the original content was not saved before the edit" — while the
correct pre-image sits in ~/.delfin/undo/<sid>/.

So after the agent overwrites an uncommitted file in a calc folder, which
is precisely the case the journal was built for, the only remaining move
is to ask the model that just made the mistake to please call its own
undo tool. list_changes, which would let the user see what an undo would
do first, has zero production callers.
"""

from __future__ import annotations

import pathlib

import pytest


_TAB = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
        / "dashboard" / "tab_agent.py")


def _source() -> str:
    return _TAB.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# The commands exist and are discoverable
# ---------------------------------------------------------------------------

def test_the_command_is_registered_in_the_palette():
    src = _source()
    assert '"/undo-file"' in src, (
        "the journal-backed undo is still unreachable from the UI")


def test_the_palette_entry_distinguishes_it_from_the_context_undo():
    """Two commands called undo-something must not read the same."""
    src = _source()
    i = src.index('"/undo-file"')
    entry = src[i:i + 260]
    assert "file" in entry.lower()
    j = src.index('("Session", "/undo"')
    assert "context" in src[j:j + 200].lower()


def test_the_handler_calls_the_journal():
    src = _source()
    assert "change_journal" in src, "the dashboard still never imports it"
    i = src.index('if cmd.startswith("/undo-file")')
    body = src[i:i + 3000]
    assert ".revert(" in body
    assert ".list_changes(" in body


def test_the_preview_is_reachable_without_reverting():
    """A user must be able to see what an undo would do first."""
    src = _source()
    i = src.index('if cmd.startswith("/undo-file")')
    body = src[i:i + 3000]
    j = body.index(".list_changes(")
    k = body.index(".revert(")
    assert j < k, "the listing branch must come before the acting branch"


def test_the_scope_is_passed_through():
    src = _source()
    i = src.index('if cmd.startswith("/undo-file")')
    body = src[i:i + 3000]
    for scope in ("last", "turn", "session"):
        assert scope in body, scope


def test_the_result_is_rendered_not_swallowed():
    """revert returns reverted / conflicts / skipped; a conflict is the
    whole point of the safety contract and has to be shown."""
    src = _source()
    i = src.index('if cmd.startswith("/undo-file")')
    body = src[i:i + 3000]
    for key in ("reverted", "conflicts", "skipped"):
        assert key in body, key


def test_it_refuses_without_a_session():
    src = _source()
    i = src.index('if cmd.startswith("/undo-file")')
    body = src[i:i + 3000]
    assert "session" in body.lower()


# ---------------------------------------------------------------------------
# The journal side still behaves
# ---------------------------------------------------------------------------

def test_the_journal_functions_are_what_the_handler_expects(tmp_path):
    from delfin.agent import change_journal as cj

    sid = "s-undo-test"
    target = tmp_path / "a.txt"
    target.write_text("original\n", encoding="utf-8")
    cj.record_change(sid, tool="write_file", path=str(target),
                     old_text="original\n", new_text="changed\n")
    target.write_text("changed\n", encoding="utf-8")

    listed = cj.list_changes(sid)
    assert listed and listed[-1].get("path", "").endswith("a.txt")

    out = cj.revert(sid, scope="last", workspace=tmp_path)
    assert set(out) >= {"reverted", "conflicts", "skipped"}
    assert target.read_text(encoding="utf-8") == "original\n"


def test_a_file_changed_since_the_edit_is_a_conflict_not_an_overwrite(tmp_path):
    from delfin.agent import change_journal as cj

    sid = "s-undo-conflict"
    target = tmp_path / "b.txt"
    target.write_text("original\n", encoding="utf-8")
    cj.record_change(sid, tool="write_file", path=str(target),
                     old_text="original\n", new_text="agent wrote this\n")
    target.write_text("the user then edited it\n", encoding="utf-8")

    out = cj.revert(sid, scope="last", workspace=tmp_path)
    assert out["conflicts"], out
    assert target.read_text(encoding="utf-8") == "the user then edited it\n"
