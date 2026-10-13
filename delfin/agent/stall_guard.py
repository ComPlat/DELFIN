"""Stall guard: count identical consecutive turn endings.

Input: one ``(kind, fingerprint)`` per turn that ended dead, through
:meth:`StallGuard.advance`, and :meth:`StallGuard.note_progress` whenever
something real happened in between (a tool ran, a turn ended on an
answer, a person sent a message, a conversation began). The engine feeds
``("dsml_leak", <tool name>)`` for a turn whose answer was a tool call
written as text, re-requested and written as text again.

Output: :data:`STOP` once the same pair arrived ``repeats`` times in a row
(default 3, setting ``agent.stall_guard.repeats``), :data:`CONTINUE`
otherwise. A different pair restarts the streak at 1; an empty or missing
fingerprint restarts it at 0, because identity that cannot be shown is
never assumed. No clock: only the order of calls matters.

One instance per engine, so each conversation counts for itself.
"""

from __future__ import annotations

#: :meth:`StallGuard.advance`: the session goes on.
CONTINUE = "continue"
#: :meth:`StallGuard.advance`: the session stops asking.
STOP = "stop"

#: Identical consecutive endings that stop the session.
DEFAULT_REPEATS = 3


class StallGuard:
    """Streak of identical consecutive ``(kind, fingerprint)`` endings."""

    def __init__(self, repeats: int = DEFAULT_REPEATS) -> None:
        """``repeats`` is parsed here, so a corrupt setting falls back to
        :data:`DEFAULT_REPEATS` instead of raising in the caller."""
        try:
            self.repeats = int(repeats)
        except (TypeError, ValueError):
            self.repeats = DEFAULT_REPEATS
        if self.repeats < 1:
            self.repeats = DEFAULT_REPEATS
        self._last: tuple[str, str] | None = None
        self._streak = 0

    def advance(self, kind: str, fingerprint: str) -> str:
        """Record one dead ending; :data:`STOP` at ``repeats`` in a row."""
        if fingerprint is None or str(fingerprint).strip() == "":
            self._last = None
            self._streak = 0
            return CONTINUE
        key = (str(kind), str(fingerprint))
        if key == self._last:
            self._streak += 1
        else:
            self._last = key
            self._streak = 1
        return STOP if self._streak >= self.repeats else CONTINUE

    def note_progress(self) -> None:
        """Something real happened: the streak starts over."""
        self._streak = 0

    @property
    def streak(self) -> int:
        """Identical consecutive endings seen since the last reset."""
        return self._streak
