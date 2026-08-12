"""The trim loops chased a target that the trimming could not reach.

``_estimate_context_tokens`` is message characters plus SYSTEM PROMPT
characters plus a provider floor. The trim loops then run

    for msg in self.messages[:protected_from]:
        if self._estimate_context_tokens() <= budget:
            break

against a budget of 70 % of the window. Two of those three terms are
constants that trimming messages cannot reduce, so on any model whose
prompt is a large share of the window the condition is unsatisfiable
before a single message is touched.

Reproduced on a 32k window with a 90 kB prompt:

    budget            22400
    system prompt     22500   ← already over, on its own
    estimate before   31535
    trimmed           8 messages, to 1057 chars each
    estimate after    27627   ← still over

The break is never reached, so the loop runs to the end of the eligible
prefix and guillotines everything in it. Next turn the prompt is the
same size, the verdict is the same, and it fires again on whatever grew
since. The agent's own answers -- the only place its conclusions,
decisions and file lists survive -- are permanently reduced to stubs.

The rule this restores is not "trim less". It is that a budget must be
compared against the part that trimming can actually change. What is
left after the prompt is the history's share; when that share is zero or
negative the window is too small for the prompt, and shredding the
history is pure loss that does not fix it. That is worth saying rather
than doing silently, because the remedy is a smaller prompt, and only
the user can decide which part of it to give up.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine


def _engine(*, window: int, prompt_chars: int, n_msgs: int = 12,
            msg_chars: int = 3000) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.context_window_tokens = window
    eng.last_system_prompt = "x" * prompt_chars
    eng._system_prompt_chars = prompt_chars
    eng._last_input_tokens = 0
    eng._trimmed_chars_since_floor = 0
    eng.auto_compact_pct = 0.95
    eng.messages = [{"role": "user", "content": "the goal"}] + [
        {"role": "assistant", "content": f"finding {i}: " + "y" * msg_chars}
        for i in range(n_msgs)
    ]
    return eng


def _assistant_lengths(eng) -> list[int]:
    return [len(m["content"]) for m in eng.messages if m["role"] == "assistant"]


# ---------------------------------------------------------------------------
# A target trimming cannot reach is not chased
# ---------------------------------------------------------------------------

def test_a_prompt_larger_than_the_budget_does_not_shred_the_history():
    """The reproduction above: 22.5k of prompt against a 22.4k budget."""
    eng = _engine(window=32_000, prompt_chars=90_000)
    before = _assistant_lengths(eng)
    eng._slide_window_trim()
    assert _assistant_lengths(eng) == before


def test_it_says_why_instead_of_trimming_silently():
    eng = _engine(window=32_000, prompt_chars=90_000)
    eng._slide_window_trim()
    info = getattr(eng, "last_compaction_info", None) or {}
    text = str(info.get("note") or info.get("reason") or "")
    assert "prompt" in text.lower(), (
        "the window is too small for the prompt and nothing says so")


def test_it_reports_that_it_trimmed_nothing():
    eng = _engine(window=32_000, prompt_chars=90_000)
    assert eng._slide_window_trim() == 0


# ---------------------------------------------------------------------------
# ...while a reachable target is still reached
# ---------------------------------------------------------------------------

def test_a_small_prompt_still_trims_the_history():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=250,
                  msg_chars=6_000)
    before = sum(_assistant_lengths(eng))
    eng._slide_window_trim()
    assert sum(_assistant_lengths(eng)) < before


def test_trimming_stops_once_the_history_fits():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=250,
                  msg_chars=6_000)
    eng._slide_window_trim()
    budget = int(eng.context_window_tokens * AgentEngine._SLIDING_WINDOW_PCT)
    assert eng._estimate_context_tokens() <= budget


def test_the_most_recent_messages_are_still_protected():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=250,
                  msg_chars=6_000)
    eng._slide_window_trim()
    kept = _assistant_lengths(eng)[-AgentEngine._KEEP_RECENT:]
    assert all(n > 3_000 for n in kept)


def test_a_history_that_already_fits_is_untouched():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=3,
                  msg_chars=500)
    before = _assistant_lengths(eng)
    assert eng._slide_window_trim() == 0
    assert _assistant_lengths(eng) == before


def test_the_user_goal_is_never_trimmed():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=250,
                  msg_chars=6_000)
    eng._slide_window_trim()
    assert eng.messages[0]["content"] == "the goal"


# ---------------------------------------------------------------------------
# The same rule for the harder stage
# ---------------------------------------------------------------------------

def test_the_hard_clear_still_runs_when_only_the_tighter_budget_is_out_of_reach():
    """Its budget is 95% of the window, the sliding window's is 70%. A
    prompt over the tighter one can still be under this one, and then
    stubbing DOES help -- refusing here would leave the session at the
    cliff for nothing."""
    eng = _engine(window=32_000, prompt_chars=90_000, n_msgs=60,
                  msg_chars=6_000)
    before = sum(_assistant_lengths(eng))
    eng._hard_clear_old_tool_results(eng.messages[:-AgentEngine._KEEP_RECENT])
    assert sum(_assistant_lengths(eng)) < before, (
        "it refused a target it could actually have reached")


def test_the_hard_clear_refuses_a_target_out_of_reach_for_it_too():
    """Prompt over 95% of the window: now nothing it stubs can help."""
    eng = _engine(window=32_000, prompt_chars=130_000)
    before = _assistant_lengths(eng)
    eng._hard_clear_old_tool_results(eng.messages[:-AgentEngine._KEEP_RECENT])
    assert _assistant_lengths(eng) == before


def test_the_hard_clear_still_works_when_it_can_help():
    eng = _engine(window=200_000, prompt_chars=4_000, n_msgs=250,
                  msg_chars=6_000)
    before = sum(_assistant_lengths(eng))
    eng._hard_clear_old_tool_results(eng.messages[:-AgentEngine._KEEP_RECENT])
    assert sum(_assistant_lengths(eng)) <= before
