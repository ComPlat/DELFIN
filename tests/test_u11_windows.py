"""Wave 14, package U11 — phase 1: the window that decides each trim/compaction.

The whole point of the observed floor (delfin/agent/observed_window.py) is
that the static window table understates: a request the backend actually
accepted is proof the window is at least that large, and every protection
of a long session fires at a share of the window. So the floor must feed
EVERY decision, not just the engine's.

Measured on three live dashboard sessions (2026-10-08): session C sent a
190,618-token request through, session A 88,851. Both became floors. Yet
the in-turn elision/compaction budget (api_client._tool_context_char_budget)
reads `capped_window(caps.context_window)` and never consults the observed
floor — so once a floor raises the engine's planning window, the two
protections of the SAME turn are deciding against DIFFERENT windows.

These tests are the measurement-first gate for phase 1: they pin the
window the in-turn budget OUGHT to plan with (the same window the engine
plans with), and record the request size the turn actually builds. They are
RED on the current code because the in-turn budget is floor-blind and
measures tool output only.
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
# Finding A — two windows decide the same turn
# ---------------------------------------------------------------------------

def test_the_engine_plans_with_the_floor_raised_window(floors):
    """The engine's planning window is at_least(table, floor).

    engine.py:1280-1281 already plans with
    _capped_context_window(_ow.at_least(model, caps.context_window)). This
    is what every slide/autocompact budget locks to.
    """
    OW.note("kit.deepseek-v4-flash", 166_379)
    assert OW.at_least("kit.deepseek-v4-flash", 131_072) == 166_379


def test_the_in_turn_budget_is_floor_blind(floors):
    """RED (current): the in-turn elision budget must be floor-aware.

    The engine plans this same turn against the floor-raised window
    (166_379 for the flash model), so the in-turn budget must be derived
    from that same window. Today _tool_context_char_budget returns
    int(capped_window(caps.context_window) * 0.45 * 4) — it reads the RAW
    table window and never the observed floor, so a floor raised the engine
    window while the in-turn elision stays anchored ~35k tokens lower.
    """
    OW.note("kit.deepseek-v4-flash", 166_379)
    caps = _caps(131_072)
    planning = OW.at_least("kit.deepseek-v4-flash", caps.context_window)
    # The floor is per model, so the budget is asked for THIS model.
    assert api_client._tool_context_char_budget(
        caps, model="kit.deepseek-v4-flash") == int(planning * 0.45 * 4)
    # Without a model there is no floor to apply: the table window stands.
    assert api_client._tool_context_char_budget(caps) == int(
        caps.context_window * 0.45 * 4)


def test_the_floor_does_not_poison_the_in_turn_budget_above_the_cap(floors, monkeypatch):
    """The cap still bounds the in-turn budget even when the floor is higher.

    A floor may raise the window only up to what the backend accepted, and
    the operator's cap must keep WINNING over it: _capped_context_window
    applies the cap AFTER at_least. The in-turn budget must do the same —
    otherwise a stale enormous floor would make the in-turn budget balloon
    past the cap that the engine still enforces.
    """
    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    OW.note("m", 166_379)
    caps = _caps(131_072)
    # The engine's window here: cap wins over the floor (40k < 166k).
    planning = min(OW.at_least("m", caps.context_window), 40_000)
    assert api_client._tool_context_char_budget(caps) == int(planning * 0.45 * 4)


# ---------------------------------------------------------------------------
# Finding B — the budget counts tool output only, not the whole request
# ---------------------------------------------------------------------------

def _msg(role, content, tool_call_id=None):
    m = {"role": role, "content": content}
    if tool_call_id:
        m["tool_call_id"] = tool_call_id
    return m


def test_the_tool_only_budget_lets_the_whole_request_overflow(floors):
    """RED (current): the guard must account for the whole request.

    The in-turn budget is sized as ~45% of the window for TOOL OUTPUT
    (api_client.py:10953-10955). The elision measures only role=="tool"
    content (_tool_chars, api_client.py:11036-11038), so assistant
    reasoning + the system prompt never count against it. Build a turn whose
    REASONING is the dominant term: the tool-only budget reports it "under
    budget", while the whole request (system + assistant + tool) already
    exceeds the window the engine plans with. RED because the whole-request
    guard is what must fire, and today only the tool-only one exists.
    """
    # Four rounds: a compaction summarises OLDER rounds and keeps the last
    # ones, so a turn of a single round has nothing it could give up.
    msgs = [_msg("user", "the user's goal for this long turn")]
    for i in range(4):
        msgs.append({
            "role": "assistant", "content": "x" * 150_000,  # reasoning
            "tool_calls": [{"id": f"c{i}", "type": "function",
                            "function": {"name": "bash", "arguments": "{}"}}]})
        msgs.append(_msg("tool", "x" * 1000, tool_call_id=f"c{i}"))
    before = sum(len(str(m.get("content") or "")) for m in msgs)

    # The turn, as the round loop would treat it: elide tool output only.
    n_elided = api_client._elide_old_tool_results(
        msgs, char_budget=api_client._tool_context_char_budget(_caps(131_072)))
    assert n_elided == 0  # tool output (1000 chars) is well under budget

    # The WHOLE-REQUEST guard is the fix (api_client._compact_between_rounds
    # with a planning window), not the tool-only elider above. An
    # assistant-heavy turn (500k chars of reasoning, tiny tool output) exceeds
    # the planning window, so the guard must compact it. RED on current code
    # (_compact_between_rounds has no planning_window parameter), GREEN once
    # the phase-2 patch lands — the guard reports the compaction.
    note = api_client._compact_between_rounds(
        msgs, 100_000, 0, "s",
        planning_window=131_072, system_chars=0)
    assert note != ""               # RED: no whole-request guard today
    assert "summarised" in note
    # ... and the request really is smaller, not only announced as such.
    assert sum(len(str(m.get("content") or "")) for m in msgs) < before


def test_request_size_per_round_is_recorded_for_the_long_turn(floors, monkeypatch):
    """RED (current): a long turn records nothing about per-round request size.

    Phase-1 mandate: "records the request size per round" and asserts no
    request exceeds the window the engine plans with. There is no hook today
    that records each round's built request size (system + history + tool)
    against the planning window; this test documents the missing surface.
    """
    import inspect
    # The in-turn loop calls _elide_old_tool_results and _compact_between_rounds
    # with a char budget, but neither returns the built request size nor checks
    # it against the engine's planning-window share. Pin the shape: there must
    # be a per-round check that yields a max_context notice when the WHOLE
    # request would exceed the planning window.
    sig = inspect.signature(api_client._compact_between_rounds)
    assert "planning_window" in sig.parameters  # RED: not present today


# ---------------------------------------------------------------------------
# Phase 2 — ONE window for every decision: context_window.planning_window
# ---------------------------------------------------------------------------

def test_planning_window_uses_the_floor_raised_then_capped_window(floors):
    """The single planning window = floor applied, then cap.

    The engine plans compaction with _capped_context_window(_ow.at_least(
    model, window)); the in-turn budget must read the SAME number. This is
    the owned, single source both sides will consult.
    """
    from delfin.agent import context_window as CW

    OW.note("kit.deepseek-v4-flash", 166_379)
    assert CW.planning_window(131_072, "kit.deepseek-v4-flash") == 166_379


def test_planning_window_applies_the_cap_after_the_floor(floors, monkeypatch):
    """A cap is a hard ceiling and must clamp a floor above it.

    The engine does _capped_context_window(at_least(...)) with the cap
    applied AFTER the floor. _tool_context_char_budget must not balloon past
    the cap just because a stale floor is large.
    """
    from delfin.agent import context_window as CW

    monkeypatch.setenv("DELFIN_CONTEXT_WINDOW_CAP", "40000")
    OW.note("m", 166_379)
    assert CW.planning_window(131_072, "m") == 40_000


def test_planning_window_without_a_floor_keeps_the_raw_window(floors):
    from delfin.agent import context_window as CW

    assert CW.planning_window(131_072, "kit.glm-5.3") == 131_072
    assert CW.planning_window(131_072, "") == 131_072  # empty model -> no floor
