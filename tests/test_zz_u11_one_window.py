"""One planning window for the engine and the in-turn budgets.

  A -- the in-turn elision budget plans with the floor-raised-then-capped
       window the engine uses (context_window.planning_window), not the raw
       caps window.
  B -- a whole-request guard compacts a turn whose assistant text exceeds
       the planning window while its tool output is small; the tool-only
       elider is unchanged.
  C -- the engine learns a model's floor from every accepted request;
       observed_window.note stores only a new maximum, so a smaller request
       never lowers it.
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


# ---------------------------------------------------------------------------
# FINDING A — the in-turn budget rises with the observed floor
# ---------------------------------------------------------------------------

def test_the_in_turn_budget_rises_with_the_observed_floor(floors):
    """With the flash-model floor (166379) the in-turn budget must plan
    against the floor-raised window, exactly like the engine does."""
    OW.note("kit.deepseek-v4-flash", 166_379)
    caps = _caps(131_072)
    budget = api_client._tool_context_char_budget(caps, model="kit.deepseek-v4-flash")
    assert budget == int(166_379 * 0.45 * 4)   # 299482, not 235929


def test_without_a_floor_the_in_turn_budget_tracks_the_engine(floors):
    caps = _caps(131_072)
    assert api_client._tool_context_char_budget(caps, model="kit.glm-5.3") \
        == int(131_072 * 0.45 * 4)


def test_a_floor_does_not_poison_a_capped_turn(floors, monkeypatch):
    """The cap (operator's hard ceiling) stays winning over a stale floor."""
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    OW.note("m", 166_379)
    caps = _caps(131_072)
    assert api_client._tool_context_char_budget(caps, model="m") \
        == int(40_000 * 0.45 * 4)


# ---------------------------------------------------------------------------
# FINDING B — a whole-request guard accounts for assistant + system
# ---------------------------------------------------------------------------

def _msgs(n_rounds, result_chars=100, reasoning_chars=500_000):
    msgs = [{"role": "user", "content": "the user's goal for this turn"}]
    msgs.append({"role": "assistant", "content": "x" * reasoning_chars})
    for i in range(n_rounds):
        msgs.append({"role": "assistant", "content": f"round {i}",
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": "{}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}",
                     "content": "x" * result_chars})
    return msgs


def test_an_assistant_heavy_turn_is_compacted_against_the_planning_window(floors):
    """The tool-only budget sees little tool output, but the WHOLE request
    (assistant reasoning + system prompt) exceeds the planning window: the
    whole-request guard must fire and report it."""
    msgs = _msgs(3, result_chars=100, reasoning_chars=500_000)
    note = api_client._compact_between_rounds(
        msgs, 100_000, 0, "s", planning_window=131_072, system_chars=0)
    assert note != ""
    assert "compacted" in note


def test_the_tool_only_elider_is_not_silently_changed(floors):
    """Finding B is a SEPARATE guard: the tool-only elider's return for a
    small tool payload must stay 0 (a turn with fine tool rows is not
    rewritten underneath)."""
    msgs = _msgs(3, result_chars=100, reasoning_chars=0)
    n = api_client._elide_old_tool_results(
        msgs, char_budget=100_000)
    assert n == 0


# ---------------------------------------------------------------------------
# FINDING C — a floor is learned for any accepted request, bounded to the max
# ---------------------------------------------------------------------------

def test_the_engine_bounds_the_floor_raise_to_what_was_accepted(floors):
    """An accepted request is proof of the window, whether or not it stayed
    within the static table window: 166379 (over the 100k plan) IS learned,
    because the backend genuinely held it. What is BOUNDED is the raise —
    observed_window.note stores only a new maximum actually accepted and never
    extrapolates, so a single real request cannot push the floor past a number
    the backend really held, and a lower request never lowers it."""
    from delfin.agent import engine as E

    def _engine():
        eng = E.AgentEngine.__new__(E.AgentEngine)
        eng.client = type("_C", (), {"model": "kit.m"})()
        eng.context_window_tokens = 100_000
        eng._floor_captured_this_turn = False
        eng._last_input_tokens = None
        eng._trimmed_chars_since_floor = 0
        return eng

    over = _engine()
    over._capture_floor(166_379)   # 166k over the 100k plan window
    assert OW.floor("kit.m") == 166_379

    later = _engine()                         # a fresh turn, same model
    later._capture_floor(90_000)   # a smaller accepted request
    assert OW.floor("kit.m") == 166_379       # the raise is never lowered
