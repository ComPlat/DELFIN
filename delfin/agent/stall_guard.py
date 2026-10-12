"""Stall guard: end a turn that repeats the same harness notice.

A model can get stuck writing a tool call as literal text instead of
calling the tool: the harness emits the same corrective notice ("was
written as text ... it did not run"), the model writes the same text
again, and without a guard the turn burns every tool round on the same
loop before the round cap ends it.

This module is the pure decision helper behind that guard. It is fed each
harness notice of the current turn as ``(kind, fingerprint)``, where
``kind`` is the class of notice (for example ``"dsml_leak"``) and
``fingerprint`` is a short stable digest of the model text that caused it,
so identical text maps to the same value and different text maps to
another. It decides, from that sequence alone, whether the turn should
keep going or end with a status line. There is no timing and no wall clock
anywhere in it: only the identity of consecutive notices matters.

Rule: the same ``(kind, fingerprint)`` N times in a row (default 3,
setting ``agent.stall_guard.repeats``) ends the turn. A different
``(kind, fingerprint)`` in between resets the streak, and explicit
progress (a real tool call, :meth:`StallGuard.note_progress`) resets it
too, so a turn that is actually making progress is never stopped on a
stale streak. A notice with an empty or missing fingerprint is never
counted toward a stall: without a fingerprint the caller could not prove
what text it was, and a guard whose job is to detect *identical* text
must not assume identity it cannot prove.
"""

from __future__ import annotations

#: Decision returned by :meth:`StallGuard.advance`: the turn should keep going.
CONTINUE = "continue"
#: Decision returned by :meth:`StallGuard.advance`: end the turn with the status.
STOP = "stop"

#: Default number of identical consecutive notices that ends a turn. The
#: caller may read ``agent.stall_guard.repeats`` and pass it in.
DEFAULT_REPEATS = 3


class StallGuard:
    """Decide, from the sequence of harness notices, whether the turn stalls.

    An instance is per-turn: create one at the start of a turn, feed every
    harness notice through :meth:`advance`, and call :meth:`note_progress`
    whenever a real tool call executes (genuine work resets the streak).
    When :meth:`advance` returns :data:`STOP` the turn should end;
    :meth:`status_line` carries the reason and what the model should do
    differently next turn.
    """

    def __init__(self, repeats: int = DEFAULT_REPEATS) -> None:
        """``repeats`` is the number of identical consecutive notices that
        end the turn. Read defensively so a missing or corrupt setting never
        throws inside the stream loop (mirrors how ``agent.max_tool_rounds``
        is resolved)."""
        try:
            self.repeats = int(repeats)
        except (TypeError, ValueError):
            self.repeats = DEFAULT_REPEATS
        if self.repeats < 1:
            self.repeats = DEFAULT_REPEATS
        # Identity of the previous notice and how many times the current one
        # has repeated consecutively, including the call that set it.
        self._last: tuple[str, str] | None = None
        self._streak = 0
        # The (kind, fingerprint) that ended the turn, if any.
        self._stalled: tuple[str, str] | None = None

    def advance(self, kind: str, fingerprint: str) -> str:
        """Record one harness notice; return :data:`CONTINUE` or :data:`STOP`.

        ``kind`` is the class of notice (for example ``"dsml_leak"``).
        ``fingerprint`` is a short stable digest of the model text that
        caused it. A notice that differs from the previous one (in kind OR
        fingerprint) resets the streak to 1. The same ``(kind, fingerprint)``
        as the previous one increments it; at ``self.repeats`` the turn
        stalls and :data:`STOP` is returned.
        """
        # A fingerprint that is empty or missing means the caller could not
        # prove what text this was. The guard exists to detect IDENTICAL text
        # repeated, so when identity cannot be proven it must never assume
        # identity: reset the streak and keep going. A non-identifiable
        # notice also separates two otherwise-identical notices, breaking
        # their consecutiveness.
        if fingerprint is None or str(fingerprint).strip() == "":
            self._last = None
            self._streak = 0
            self._stalled = None
            return CONTINUE
        key = (str(kind), str(fingerprint))
        if key == self._last:
            self._streak += 1
        else:
            self._last = key
            self._streak = 1
            # A different notice ends any earlier stall: the last advance no
            # longer decided the turn should end.
            self._stalled = None
        if self._streak >= self.repeats:
            self._stalled = key
            return STOP
        return CONTINUE

    def note_progress(self) -> None:
        """Record that the turn made real progress (a tool actually ran).

        Resets the streak so a turn that works between two identical notices
        is never stopped on a stale streak.
        """
        self._streak = 0
        self._stalled = None

    @property
    def streak(self) -> int:
        """How many consecutive identical notices have been seen so far."""
        return self._streak

    def has_stalled(self) -> bool:
        """Whether the last :meth:`advance` decided the turn should end."""
        return self._stalled is not None

    def status_line(self) -> str:
        """The status to end the turn with, or ``""`` when not stalled.

        Names the notice class, that it repeated, and what the model should
        do differently instead of emitting the same text a fourth time.
        """
        if not self._stalled:
            return ""
        kind, fingerprint = self._stalled
        return (
            f"The harness notice (kind={kind!r}, text={fingerprint!r}) was "
            f"repeated {self.repeats} times in a row this turn without real "
            "progress, so the turn was ended to avoid burning rounds on a "
            "stall. Do not repeat the same text again: call the tool "
            "properly (as a real tool call, not as text), change approach, "
            "or answer the user plainly."
        )
