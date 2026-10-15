"""One step of a permission-mode cycle (Shift+Tab, the mode chip).

Input: the ring of modes in cycle order, the current mode, and whether a
turn is running. Output: the next mode. While a turn runs the step passes
over ``plan``: every step takes effect at once, so walking the ring from
bypass back round to acceptEdits put a running turn into plan for a
moment, and the tool call that landed then was refused (bug archive
2026-09-16..10-14: 15 plan-mode refusals in 9 sessions; one push refused
mid-cycle, after which the model stopped). An explicit choice of plan
(/mode plan, the dropdown) is not a cycle step and still applies at once.

Rejected: applying the cycled mode only at the next turn boundary (a
user who cycles to a stricter mode mid-turn wants it now); a debounce
timer (the posture between two key presses would still be live).
"""

from __future__ import annotations

from typing import Sequence

__all__ = ["next_in_cycle"]

PLAN = "plan"


def next_in_cycle(order: Sequence[str], current: str, *,
                  turn_running: bool, plan: str = PLAN) -> str:
    """The mode after ``current`` in ``order``; ``plan`` is passed over
    while ``turn_running``. A ``current`` not in the ring starts it at its
    first entry that may be used. A ring with nothing but ``plan`` is
    returned unchanged."""
    ring = [str(m) for m in order]
    if not ring:
        return current
    usable = [m for m in ring if not (turn_running and m == plan)] or ring
    if current in ring:
        i = ring.index(current)
        for step in range(1, len(ring) + 1):
            cand = ring[(i + step) % len(ring)]
            if cand in usable:
                return cand
    return usable[0]
