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
    assert api_client._tool_context_char_budget(caps) == int(planning * 0.45 * 4)


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
    user = _msg("user", "the user's goal for this long turn")
    reasoning = _msg("assistant", "x" * (500_000))  # ~125k tokens of reasoning
    tool = _msg("tool", "x" * 1000, tool_call_id="c0")
    msgs = [user, reasoning, tool]

    # The turn, as the round loop would treat it: elide tool output only.
    n_elided = api_client._elide_old_tool_results(
        msgs, char_budget=api_client._tool_context_char_budget(_caps(131_072)))
    assert n_elided == 0  # tool output (1000 chars) is well under budget

    # The whole request (chars / 4) is what the backend receives, and — with
    # the system prompt and this turn's reasoning included — already exceeds
    # the whole window. RED: the tool-only guard reports "under budget"
    # (n_elided == 0) while the request the backend would get is over-window.
    whole_tokens = sum(len(str(m.get("content", ""))) for m in msgs) // 4
    assert whole_tokens <= api_client._tool_context_char_budget(_caps(131_072)) // 4


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
