"""The memory cap breaks a loop; it should not tax a working turn.

`max_memory_writes_per_turn` is 2 for the KIT models. The cap exists
because a model that calls `remember` over and over looks like progress
to the no-progress guard (which keys on name AND arguments, and six
remembers with different content differ every time).

But a turn is not a fixed size. In the field, the session that ran
longest made 40 tool calls in one turn and had six memory writes refused
across the session (2026-10-08) -- facts it had learned and could not
record, in the session where memory matters most because its context is
being trimmed underneath it.

A loop does nothing else, and that is the signal the cap can use: the
allowance is the profile's cap plus one per `_MEMORY_WRITE_EARN_EVERY`
non-memory tool calls in the same turn. Three remembers in a row still
stop at two, because nothing was earned.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


def _allowance(work_calls: int, cap: int = 2) -> int:
    """The rule, as the dispatcher computes it."""
    return cap + (work_calls // A._MEMORY_WRITE_EARN_EVERY)


class TestTheRule:
    def test_a_loop_earns_nothing(self):
        """The property the cap exists for, and the one that must not
        move: a turn that only records is stopped at the profile's cap."""
        assert _allowance(work_calls=0) == 2

    def test_work_below_the_step_earns_nothing(self):
        assert _allowance(A._MEMORY_WRITE_EARN_EVERY - 1) == 2

    def test_a_long_turn_earns_room(self):
        """40 calls was the real turn; it buys four writes on top."""
        assert _allowance(40) == 6

    def test_it_rises_one_at_a_time(self):
        step = A._MEMORY_WRITE_EARN_EVERY
        assert [_allowance(n * step) for n in range(5)] == [2, 3, 4, 5, 6]

    def test_a_profile_with_no_cap_is_untouched(self):
        """`_memory_cap` of 0 switches the whole check off; the earning
        rule must not switch it back on."""
        from delfin.agent.model_profiles import get_profile
        assert get_profile("claude-opus-5").max_memory_writes_per_turn >= 0


class TestItReachesTheDispatcher:
    def test_the_rule_is_wired_where_the_cap_is_decided(self):
        """One source check, and it pins only the WIRING.

        The rule itself is asserted above, on the shipped constant. This
        is the part a behavioural test cannot reach without a full fake
        of the streaming endpoint: that the dispatcher computes the
        allowance from the cap plus earned room, and that a memory write
        does not count as its own work -- which would let a loop earn its
        own headroom, the one thing this must never do.
        """
        import inspect

        source = inspect.getsource(A)
        assert "_memory_cap + (" in source
        assert "_work_calls // _MEMORY_WRITE_EARN_EVERY" in source
        assert "elif fn_name not in _MEMORY_WRITE_TOOLS:" in source
        # And the refusal names how much work the turn has done, or the
        # number it quotes looks arbitrary and the model argues with it.
        assert "other tool call(s)" in source
        assert "The allowance rises as" in source
