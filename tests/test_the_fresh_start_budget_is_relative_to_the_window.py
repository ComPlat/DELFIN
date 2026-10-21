"""The fresh-start budget cannot exceed the model's context window.

`_DEFAULT_FRESH_BUDGET` is 900,000 tokens -- the operator's number, taken
from a session that re-sent more than that per TURN (the sum over every
round). The budget is compared against the CURRENT context, which the
model's window bounds. GLM 5.3's window is 131,072, so a 900k budget was
6.9x the whole window and could never fire: the fresh start was dead for
every KIT model while one session paid 20.5M input tokens over 22 turns.

The window itself, not a fraction: the sliding trim and compaction run
first; a fresh start is the harder reset for a context still over the
window after both.
"""

from __future__ import annotations

from delfin.agent import repl as R


class _Engine:
    def __init__(self, window=0, explicit=None):
        self.context_window_tokens = window
        if explicit is not None:
            self.context_budget = explicit


def test_a_small_window_bounds_the_budget():
    assert R._resolve_fresh_budget(None, _Engine(window=131_072)) == 131_072


def test_a_huge_window_keeps_the_operators_ceiling():
    assert R._resolve_fresh_budget(None, _Engine(window=2_000_000)) == R._DEFAULT_FRESH_BUDGET


def test_no_window_known_falls_back_to_the_ceiling():
    assert R._resolve_fresh_budget(None, _Engine(window=0)) == R._DEFAULT_FRESH_BUDGET


def test_an_explicit_setting_still_wins():
    assert R._resolve_fresh_budget(50_000, _Engine(window=131_072)) == 50_000
    assert R._resolve_fresh_budget(None, _Engine(window=131_072, explicit=70_000)) == 70_000


def test_zero_still_disables():
    assert R._resolve_fresh_budget(0, _Engine(window=131_072)) == 0


def test_it_now_fires_for_a_kit_sized_context():
    """The point: a 166k context on a 131k window -- observed -- starts
    fresh, where before it was 734k short of the budget."""
    eng = _Engine(window=131_072)
    eng.messages = [{"role": "user", "content": "x" * 4 * 166_000}]
    budget = R._resolve_fresh_budget(None, eng)
    assert R._context_tokens(eng) > budget
