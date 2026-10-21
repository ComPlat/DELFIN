"""A trim that shaves to the trigger shaves on every turn.

The sliding-window trim fires above `_SLIDING_WINDOW_PCT` (0.70 of the
model's context window) and shortens the oldest long messages IN PLACE.
It used to stop as soon as the estimate was back under that same 0.70 --
so the next turn's growth put it over again, and it trimmed again. Every
turn, for the rest of the session.

Why that is expensive: a mutated earlier message is a changed prefix, and
a changed prefix is a cold prefix cache for everything after it. Measured
in the field across three sessions (2026-10-08): rounds INSIDE one turn
reported 84-92% of their input cached, while the first round of each new
turn reported 0-9% -- about the size of the system prompt and nothing
more. One session paid 20.5 M input tokens over 22 turns.

So the trim now cuts to a FLOOR (0.50) and the session grows from there.
The number this file exists to produce is how many turns out of N mutate
the history under each rule.

The cost side, measured too, because a cut that takes more history is not
free: this drops more per cut. The dropped middles stay retrievable
(`history_get('elided:…')`), so the price is a lookup rather than a loss.
"""

from __future__ import annotations

import pytest

from delfin.agent import engine as engine_mod


_WINDOW = 100_000          # round numbers: 0.70 -> 70k, 0.50 -> 50k


class _Engine:
    """A real AgentEngine, with only what the trim reads filled in.

    Built with `__new__` rather than assembled from borrowed methods: a
    stand-in that copies five methods measures those five, and the trim
    reaches further than that (the candidate order, the irreducible
    system-prompt term, the machine-turn test). The real class is the
    thing under measurement.
    """

    def __new__(cls, floor_pct):
        eng = engine_mod.AgentEngine.__new__(engine_mod.AgentEngine)
        eng.context_window_tokens = _WINDOW
        eng._SLIDING_WINDOW_FLOOR_PCT = floor_pct
        eng.messages = [{"role": "user", "content": "build the oracle"}]
        eng._trimmed_chars_since_floor = 0
        eng._last_input_tokens = 0
        eng.last_system_prompt = "x" * 40_000      # ~10k tokens, fixed
        eng._system_prompt_chars = 40_000
        eng.auto_compact_pct = 0.95
        eng.session_id = "measure"
        eng._elide_original = lambda *a, **k: ""
        return eng


def _append_turn(eng, size: int) -> None:
    eng.messages.append({"role": "assistant", "content": "a" * size})
    eng.messages.append({"role": "user", "content": "next step"})


def _run(floor_pct, turns=30, growth=24_000):
    """N turns of growth. Returns (trimming turns, final estimate)."""
    eng = _Engine(floor_pct)
    trimmed_on = 0
    for _ in range(turns):
        _append_turn(eng, growth)
        if eng._should_slide():
            before = [str(m.get("content")) for m in eng.messages]
            eng._shorten_oldest_non_goal_messages()
            after = [str(m.get("content")) for m in eng.messages]
            if before != after:
                trimmed_on += 1
    return trimmed_on, eng._estimate_context_tokens()


class TestHowOftenTheHistoryIsRewritten:
    def test_stopping_at_the_trigger_rewrites_almost_every_turn(self):
        """The previous rule, run as a control."""
        at_trigger, _ = _run(engine_mod.AgentEngine._SLIDING_WINDOW_PCT)
        assert at_trigger >= 20, (
            f"the control should rewrite on most turns; got {at_trigger}")

    def test_cutting_to_the_floor_rewrites_far_fewer(self):
        at_floor, _ = _run(engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT)
        at_trigger, _ = _run(engine_mod.AgentEngine._SLIDING_WINDOW_PCT)
        assert at_floor < at_trigger / 2, (
            f"cutting to the floor rewrote on {at_floor} turns of 30, "
            f"stopping at the trigger on {at_trigger}")

    def test_the_shipped_floor_is_below_the_trigger(self):
        """The two numbers are what make this hysteresis rather than a
        threshold; equal numbers would restore the old behaviour in a
        form nobody would notice."""
        assert (engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT
                < engine_mod.AgentEngine._SLIDING_WINDOW_PCT)

    def test_it_still_gets_under_the_trigger(self):
        """A cut that does not reach its goal is worse than none: the
        next turn trims again and nothing was gained."""
        _, final = _run(engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT)
        assert final <= int(_WINDOW
                            * engine_mod.AgentEngine._SLIDING_WINDOW_PCT), final


class TestWhatItCosts:
    def test_one_cut_drops_more_history(self):
        """Stated rather than hidden: the floor takes more per cut. The
        dropped middles are elided, not deleted."""
        kept = {}
        for name, pct in (("trigger",
                           engine_mod.AgentEngine._SLIDING_WINDOW_PCT),
                          ("floor",
                           engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT)):
            eng = _Engine(pct)
            # Enough growth to be over the TRIGGER before the single cut:
            # four turns leaves the estimate at 34k, under both marks, so
            # neither rule trimmed and the first run of this measured
            # nothing.
            for _ in range(12):
                _append_turn(eng, 24_000)
            eng._shorten_oldest_non_goal_messages()
            kept[name] = eng._estimate_context_tokens()
        assert kept["floor"] < kept["trigger"], kept

    def test_a_user_goal_is_never_trimmed(self, ):
        """The limit that must survive any threshold change."""
        eng = _Engine(engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT)
        goal = "the goal: " + "g" * 5_000
        eng.messages[0]["content"] = goal
        for _ in range(8):
            _append_turn(eng, 24_000)
        eng._shorten_oldest_non_goal_messages()
        assert eng.messages[0]["content"] == goal

    def test_the_recent_messages_are_untouched(self):
        eng = _Engine(engine_mod.AgentEngine._SLIDING_WINDOW_FLOOR_PCT)
        for _ in range(8):
            _append_turn(eng, 24_000)
        tail = [str(m.get("content")) for m in eng.messages[-eng._KEEP_RECENT:]]
        eng._shorten_oldest_non_goal_messages()
        assert [str(m.get("content"))
                for m in eng.messages[-eng._KEEP_RECENT:]] == tail
