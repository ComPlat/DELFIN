"""Per-user-turn policy for how many DSML re-requests a turn may make.

A leaked DSML tool call (written as literal ``<invoke name=...>...</invoke>``
text instead of calling the tool) means nothing ran and the turn is re-asked
to make the real call.  The original U5 one-shot latch allowed exactly one
re-request per user turn.  That was too strict: a long autonomous turn with
dozens of tool calls can legitimately leak a second call after real progress,
and the one-shot latch then ended the turn — the 06.11 defect (a session sat
idle ~4 hours with uncommitted work because its second ``git commit`` leak was
not re-requested).

This module holds the re-request DECISION as a small, pure, tested state
object.  The engine/api_client hooks in this package only *call* it; every
bound and every rule lives here so it can be tested in isolation.

Policy::

    first re-request  -> allowed iff the leaked call dominates the answer
                         (today's behaviour, unchanged)
    later re-requests  -> allowed iff the leaked call dominates AND at least
                         one REAL tool call has completed since the last
                         re-request (progress) AND fewer than
                         ``max_per_turn`` re-requests happened this turn

Two leaks back-to-back without a real tool call between them stay annotated
and stopped, exactly as today — the re-request budget stops a loop, never
recovery after real progress.
"""
from __future__ import annotations


class DSMLReplayPolicy:
    """The per-turn DSML re-request budget.

    The engine owns one instance per conversation and re-arms it (or builds a
    fresh one) on every real user turn; ``max_per_turn`` bounds how many
    re-requests a single user turn may make in total.
    """

    #: Module-level default so the hook wiring can read one number.
    DEFAULT_MAX_PER_TURN = 3

    def __init__(self, max_per_turn: int = DEFAULT_MAX_PER_TURN) -> None:
        self._max = max_per_turn
        self._used = 0
        #: Real tool calls that completed since the last re-request.
        self._real_since = 0

    @property
    def used(self) -> int:
        """How many re-requests this turn has consumed so far."""
        return self._used

    @property
    def exhausted(self) -> bool:
        """True when the per-turn budget is spent and no more re-requests are
        allowed this turn."""
        return self._used >= self._max

    def note_real_tool_call(self) -> None:
        """Record that a REAL tool call completed (progress).  Call this each
        time a tool result arrives; it is what lets a later leak be
        re-requested again."""
        self._real_since += 1

    def allowed_for(self, dominates: bool) -> bool:
        """Whether a re-request may be issued for the current leak.

        ``dominates`` is what ``leaked_tool_call_dominates`` already decided:
        the whole answer IS the leaked call.  Refused re-requests consume
        nothing, so a non-dominant answer never spends budget.
        """
        if not dominates:
            return False
        if self._used >= self._max:
            # Budget spent -- covers max_per_turn <= 0 too, where the budget
            # is exhausted before any re-request.  Must agree with exhausted.
            return False
        if self._used == 0:
            # The first re-request of the turn: allowed exactly as before,
            # no progress required.
            return True
        # A later re-request needs real progress since the last one AND room
        # in the per-turn budget.
        return self._real_since > 0

    def spend(self) -> None:
        """Consume one re-request (the caller is about to re-request)."""
        self._used += 1
        # Progress is counted "since the last re-request": after a re-request
        # the next one needs NEW progress.
        self._real_since = 0


def dsml_leak_notice(name: str, *, limit_reached: bool = False) -> str:
    """The user-visible notice for a leaked DSML tool call written as text.

    One source of truth for both consumers: the api_client notice (emitted
    per leaking request) and the engine's exhaustion caveat (LK1 target 5:
    the notice says when the per-turn re-request bound was reached, not a
    bare "did not run" that implies a re-request is still coming).

    ``name`` is the tool the model wrote as ``<invoke name=...>`` text.
    ``limit_reached``: the per-turn re-request budget is spent, so no
    re-request will happen this turn — the notice names the bound.
    """
    base = (f"\n\n⚠️ a call to {name} was written as text instead of "
            f"calling the tool; it did not run")
    if limit_reached:
        return base + ("; the re-request limit for this turn is reached; "
                       "the call was not run again.\n")
    return base + ".\n"
