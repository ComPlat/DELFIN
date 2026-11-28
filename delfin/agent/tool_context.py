"""Per-session tool-executor state, lifted off the process-wide ``_doc_executor``.

``_DocToolExecutor`` in :mod:`delfin.agent.api_client` is a module-level
singleton shared by every session in a process: the dashboard's several
sessions, and a parent agent with a background sub-agent. Three pieces of
state that describe a single turn used to live directly on it and were
therefore shared by every session:

* ``_turn_change_seqs`` / ``_turn_seqs_by_session`` — the undo journal seqs
  recorded this turn, read by undo_changes(scope="turn");
* ``_pending_view_images`` — images a ``view_image`` call opened, drained by
  the next model round that finishes;
* ``live_messages`` — the conversation ``history_search`` resolves ``live:``
  refs against.

Sharing those between sessions was a bug: one session starting a turn emptied
the other's undo window and, on the next write, wiped the other session's
per-session journal entry.

This module moves that state into a :class:`ToolContext`: one context per
session, owned (created, reset at turn start, drained) by the session's
client. The executor routes every per-turn / per-session field access through
a context, so two sessions in one process no longer share any of it.

Compatibility: the executor keeps every old attribute name and method working
(``_turn_change_seqs``, ``_turn_seqs_by_session``, ``_turn_seqs_seen``,
``live_messages``, ``_pending_view_images``, ``_note_change_seq``,
``_turn_seqs_for``) by delegating them to a context. Callers that pass no
context to ``execute`` fall back to a single default context per process, as
before.
"""
from __future__ import annotations

import contextvars as _contextvars
from typing import Optional

# The ToolContext set active for the duration of ONE execute() call. When a
# client passes its own context it is authoritative for that call; methods that
# run outside any execute() call (a direct `_note_change_seq` or
# `_turn_seqs_for`, the round-loop drain) read through the executor's per-session
# registry instead (_ctx_for). A ContextVar keeps a nested subagent stream on the
# same thread from leaking one session's context into another -- the same
# mechanism _MODE_PIN in api_client.py uses.
ACTIVE_TOOL_CONTEXT: "_contextvars.ContextVar[Optional[object]]" = (
    _contextvars.ContextVar("delfin_active_tool_context", default=None))

_ACTIVE_TOOL_CONTEXT = ACTIVE_TOOL_CONTEXT


class ToolContext:
    """One session's per-turn tool-executor state.

    An instance owns exactly one session's undo-window seqs, its pending
    ``view_image`` results and its live conversation for ``history_search``.
    The owner (the session's client) resets the turn and drains the pending
    images; the executor routes its per-session field access here.

    The undo logic is the one that used to live on ``_DocToolExecutor``
    (``_note_change_seq`` / ``_turn_seqs_for``), kept verbatim so behaviour
    does not change: the flat list still tracks the turn boundary, because the
    turn reset only ever touches this one context's window now.
    """

    __slots__ = (
        "_turn_change_seqs",
        "_turn_seqs_by_session",
        "_turn_seqs_seen",
        "_pending_view_images",
        "live_messages",
    )

    def __init__(self) -> None:
        # Undo journal: seqs of file changes captured THIS turn (reset by the
        # client at turn start) so undo_changes(scope="turn") knows the turn
        # boundary. Per-session: this context belongs to one session.
        self._turn_change_seqs: list[int] = []
        # The session-keyed view of the same window. Kept keyed by session id
        # for compatibility with the old flat attribute, though a context
        # already belongs to one session.
        self._turn_seqs_by_session: dict[str, list[int]] = {}
        self._turn_seqs_seen: int = 0
        # Images a view_image call opened this turn, drained by the next round.
        self._pending_view_images: list = []
        # The live conversation the session's history_search resolves live:
        # refs against, bound just-in-time before each tool call.
        self.live_messages: Optional[list] = None

    # -- undo-window journal -------------------------------------------------

    def note_change_seq(
        self,
        rec: "Optional[dict]",
        perms: "Optional[object]" = None,
    ) -> None:
        """Remember a journal seq as belonging to the CURRENT turn.

        Keyed by session: this context belongs to ONE session, so a turn-scoped
        undo never reverts a sibling's writes. The flat list is still fed
        because the turn reset and the document-write capture path use it.
        """
        if not rec:
            return
        seq = rec.get("seq")
        if not isinstance(seq, int):
            return
        if not self._turn_change_seqs:
            self._turn_change_seqs = []
        if not self._turn_seqs_by_session:
            self._turn_seqs_by_session = {}
        # The turn boundary is marked by the client emptying the flat list.
        # Detecting that here keeps the per-session cursor honest without a
        # second reset point that could be forgotten.
        if len(self._turn_change_seqs) < self._turn_seqs_seen:
            self._turn_seqs_by_session = {}
        self._turn_change_seqs.append(seq)
        self._turn_seqs_seen = len(self._turn_change_seqs)
        sid = getattr(perms, "task_session_id", "") or ""
        self._turn_seqs_by_session.setdefault(sid, []).append(seq)

    def turn_seqs_for(self, session_id: str) -> list[int]:
        """The journal seqs written this turn BY THIS SESSION.

        Falls back to the shared flat list only when the per-session cursor
        holds nothing for it, so a turn that only produced document writes
        (whose capture path feeds the flat list) still has a turn to undo.
        """
        per = self._turn_seqs_by_session or {}
        seqs = list(per.get(session_id or "", []) or [])
        if seqs:
            return seqs
        return list(self._turn_change_seqs or [])

    def begin_turn(self) -> None:
        """Reset the turn boundary: the next write starts a clean undo window.

        BOTH views of the window, not just the flat list. ``turn_seqs_for``
        reads the session-keyed dict FIRST and only falls back to the flat
        list, so clearing the list alone left the previous turn's seqs
        standing: a turn that wrote no file still reported the seqs of the
        turn before it, and ``undo_changes(scope="turn")`` reverted a file
        this turn never touched -- then reported that as this turn's undo.
        The terminal's ``/undo-file turn`` reads the same call and did the
        same thing.

        The lazy reset in :meth:`note_change_seq` stays: it is the guard for
        a caller that empties ``_turn_change_seqs`` through the executor's
        property without starting a turn here (a document write appends to
        that list directly), which this method cannot see.
        """
        self._turn_change_seqs = []
        self._turn_seqs_by_session = {}

    # -- pending view images -------------------------------------------------

    def take_pending_images(self) -> list:
        """Drain and return the images opened this turn (the round loop's view)."""
        imgs = self._pending_view_images
        self._pending_view_images = []
        return list(imgs)


def for_client(client) -> "ToolContext":
    """The client's own context, created on first use.

    A client built without ``__init__`` raises ``AttributeError`` in the
    round loop otherwise. That is not a corner case here: 31 test files
    build one with ``OpenAIClient.__new__`` and set only the attributes
    they need, and any subclass that skips ``__init__`` behaves the same.
    Creating it on demand keeps the per-session state exactly as it is for
    a client that was built normally.
    """
    ctx = getattr(client, "_tool_context", None)
    if ctx is None:
        ctx = ToolContext()
        try:
            client._tool_context = ctx
        except Exception:      # a client that refuses attributes
            pass
    return ctx
