"""Wave 14 U11 — the long-turn measurement: 300 rounds, request size per round.

Phase-1/2 mandate: a long autonomous session is ONE turn of hundreds of tool
rounds, and nothing records whether each round's built request stays within
the window the engine plans with. This drives a simulated 300-round turn
through the REAL in-turn machinery that the round loop uses every round
(_tool_context_char_budget, _elide_old_tool_results, _compact_between_rounds)
and records the whole-request size after every round.

RED on current code: without the phase-2 floor-aware budget and the
whole-request guard, assistant reasoning + tool rows accumulate unbounded —
the request grows past the planning window the engine plans with and nothing
stops it. GREEN once the operator builds .gate/u11_one_window.patch (the
in-turn budget plans against planning_window = floor-raised-then-capped, and
_compact_between_rounds compacts a whole request that exceeds the window
even when its tool output is small).
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.agent import api_client
from delfin.agent import observed_window as OW


@pytest.fixture
def floors(tmp_path, monkeypatch):
    monkeypatch.setattr(OW, "_PATH", tmp_path / "observed_windows.json")
    return tmp_path / "observed_windows.json"


def _caps(context_window: int) -> "SimpleNamespace":
    return SimpleNamespace(context_window=context_window)


def _request_tokens(msgs: list, system_chars: int = 0) -> int:
    """Whole-request token estimate exactly as the engine reads it
    (engine.py:5100-5133: message chars / 4 + system-prompt chars / 4)."""
    total = sum(len(str(m.get("content") or "")) for m in msgs)
    return (total + int(system_chars or 0)) // 4


def _simulate_loop(n_rounds, *, window, model, system_chars=1000,
                   reasoning_chars=1800, tool_chars=200):
    """Drive one long turn the way api_client.stream_message does: per round
    re-derive the tool budget, elide old tool results, compact between rounds,
    then append this round's assistant + tool rows. Records the whole-request
    token size after each round. Returns (per_round_sizes, final_messages)."""
    msgs = [{"role": "user", "content": "the user's goal for this session"}]
    caps = _caps(window)
    sizes = []
    for i in range(n_rounds):
        budget = api_client._tool_context_char_budget(caps, model=model)
        api_client._elide_old_tool_results(
            msgs, char_budget=budget, protect_before=0)
        api_client._compact_between_rounds(
            msgs, budget, 0, "s",
            planning_window=window, system_chars=system_chars)
        msgs.append({"role": "assistant",
                     "content": "x" * reasoning_chars,
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": f"{{}} // {i}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}",
                     "content": "x" * tool_chars})
        sizes.append(_request_tokens(msgs, system_chars))
    return sizes, msgs


def test_a_long_turn_records_request_size_and_stays_within_the_window(floors):
    """300 rounds, whole-request size recorded per round — every round must
    stay within the window the engine plans with (131072 tokens, no floor).

    RED on current code: the tool-only budget does not account for assistant
    reasoning + system prompt, so the request grows past the window and no
    round ever records a value that would trip a guard."""
    sizes, _ = _simulate_loop(
        300, window=131_072, model="kit.glm-5.3-flash", system_chars=36577)
    assert len(sizes) == 300
    assert max(sizes) <= 131_072                 # RED: exceeds without a guard
    # The growth is monotone until a compaction fires; record the peak too.
    assert sizes[-1] <= 131_072


def test_a_floor_raises_the_plan_and_the_budget_must_follow(floors):
    """With an observed floor the engine plans against a LARGER window; the
    in-turn budget and the per-round guard must plan against that SAME raised
    window (finding A). RED on current code: the in-turn budget stays anchored
    to the raw window, so a turn that the engine plans to hold (166k) gets
    compacted far earlier than it needs to."""
    OW.note("kit.deepseek-v4-flash", 166_379)
    caps = _caps(131_072)
    budget = api_client._tool_context_char_budget(caps, model="kit.deepseek-v4-flash")
    assert budget == int(166_379 * 0.45 * 4)     # RED: today 235929 (raw)


def test_the_whole_request_guard_covers_assistant_heavy_turns(floors):
    """Finding B: a turn whose assistant reasoning (not its tool output) is the
    dominant term must still be compacted against the planning window. The
    whole-request guard is separate from the tool-only elider.

    Uses a 4-round turn (goal + 4 complete tool rounds), because the in-turn
    compactor only fires when there are rounds older than the kept tail to
    summarise; a single-round toy turn could never compact and would prove
    nothing.

    RED on current code: _compact_between_rounds has no planning_window
    parameter and no whole-request guard, so an assistant-heavy turn reports no
    compaction and the request stays over the window."""
    msgs = [{"role": "user", "content": "goal", "tool_calls": None}]
    for i in range(4):
        msgs.append({"role": "assistant",
                     "content": "x" * 170_000,
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": "{}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}",
                     "content": "x" * 100})
    note = api_client._compact_between_rounds(
        msgs, 100_000, 0, "s",
        planning_window=131_072, system_chars=0)
    assert note != ""                            # RED: no whole-request guard today
    assert "compacted" in note
