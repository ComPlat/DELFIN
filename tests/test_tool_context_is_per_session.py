"""ToolContext keeps ONE session's per-turn tool-executor state apart from
every other session's.

A process runs several sessions (the dashboard's several sessions, a parent
agent with a background sub-agent) through one shared ``_DocToolExecutor``.
Before this change the undo-window seqs, the pending ``view_image`` images and
the ``live_messages`` for ``history_search`` all lived on the executor itself,
so two sessions shared them: one session starting a turn emptied the other's
undo window and, on the next write, wiped the other session's per-session
journal entry.

``ToolContext`` is the object the session's client owns: one instance per
session. These tests pin the per-session semantics of the context from the
outside — the state an undo window and a history-search binding must keep
once per session.
"""
from __future__ import annotations

from delfin.agent.tool_context import ToolContext


class _Perms:
    """Minimal stand-in for KitToolPermissions: the executor reads only
    ``task_session_id`` out of perms in this context."""
    task_session_id = ""

    def __init__(self, sid: str) -> None:
        self.task_session_id = sid


def test_two_contexts_keep_separate_undo_windows_even_across_a_turn():
    """Session A's window must hold A's seqs after session B runs a whole turn
    (begin_turn + two writes): the turn boundary of B touches nothing of A's
    window. This is the exact fix the step-3 executor test exercises."""
    a = ToolContext()
    b = ToolContext()

    # A's turn: two writes.
    a.note_change_seq({"seq": 11}, _Perms("A"))
    a.note_change_seq({"seq": 12}, _Perms("A"))

    # B starts a turn while A's turn is still open, then writes.
    b.begin_turn()
    b.note_change_seq({"seq": 21}, _Perms("B"))
    b.note_change_seq({"seq": 22}, _Perms("B"))

    assert a.turn_seqs_for("A") == [11, 12], (
        "session A's undo window leaked B's turn into it")
    assert b.turn_seqs_for("B") == [21, 22], "session B's own window is wrong"


def test_begin_turn_starts_a_fresh_window_for_that_session_only():
    """begin_turn resets THIS context's window; it never touches a sibling."""
    a = ToolContext()
    b = ToolContext()
    a.note_change_seq({"seq": 11}, _Perms("A"))
    b.note_change_seq({"seq": 21}, _Perms("B"))
    a.begin_turn()
    a.note_change_seq({"seq": 12}, _Perms("A"))
    # A's new turn has only seq 12.
    assert a.turn_seqs_for("A") == [12]
    # B's earlier seq survived A's turn boundary.
    assert b.turn_seqs_for("B") == [21]


def test_pending_images_are_taken_per_context():
    """The image session A opened is drained from A's context, not handed to
    whichever session's round finishes first."""
    a = ToolContext()
    b = ToolContext()
    a._pending_view_images.append("img-a")
    assert b.take_pending_images() == []
    assert a.take_pending_images() == ["img-a"]


def test_undo_window_falls_back_to_flat_list_when_no_session_entry():
    """A turn that only produced document writes (whose capture path feeds
    the flat list) still has a window to undo, as before."""
    ctx = ToolContext()
    ctx._turn_change_seqs = [5, 6]
    # No per-session entry: falls back to the flat list.
    assert ctx.turn_seqs_for("solo") == [5, 6]


def test_note_change_seq_ignores_non_int_or_missing_seq():
    """Robustness carried over from the executor: a rec without an int seq is
    a no-op and must not advance the window."""
    ctx = ToolContext()
    ctx.note_change_seq(None, _Perms("A"))
    ctx.note_change_seq({"seq": "bogus"}, _Perms("A"))
    assert ctx.turn_seqs_for("A") == []
