"""The in-turn context gap (control for the design, not yet the fix).

RED today: the elision that governs context WITHIN one long turn derives
its budget from the model's RAW context window and never sees
DELFIN_CONTEXT_WINDOW_CAP / agent.context_window_cap. The engine caps its
own window (engine.py _capped_context_window) but nothing tells the
client, so with a 40 000-token cap on a 131 000-token model the elision
lets tool output grow to ~45% of 131k — three times the whole capped
window — before it cuts anything, and the request dies with
context_length_exceeded instead.
"""

from __future__ import annotations

import os

from delfin.agent import api_client
from delfin.agent.engine import _capped_context_window


class _Caps:
    """Minimal caps stand-in: what _tool_context_char_budget reads."""
    def __init__(self, window: int) -> None:
        self.context_window = window


def test_engine_cap_lowers_the_engine_window(monkeypatch):
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    assert _capped_context_window(131_000) == 40_000


def test_elision_budget_respects_the_same_cap(monkeypatch):
    """The in-turn elision must budget against the CAPPED window.

    45% of a 40k-token window is 18k tokens ~ 72k chars. Budgeting
    against the raw 131k window gives ~235k chars — three times the
    whole capped window; the turn overflows before the elision fires.
    """
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    budget = api_client._tool_context_char_budget(_Caps(131_000))
    assert budget <= 40_000 * 0.45 * 4, (
        f"elision budget {budget} chars ignores the 40k cap; "
        f"a capped turn overflows before eliding"
    )


def test_no_cap_keeps_the_scaled_budget():
    """Without a cap the bug-172455 scaling stays as it is."""
    budget = api_client._tool_context_char_budget(_Caps(131_000))
    assert budget >= 60_000
