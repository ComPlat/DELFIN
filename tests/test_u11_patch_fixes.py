"""Wave 14 U11 — the three operator findings on .gate/u11_one_window.patch.

The operator rejected the phase-2 patch with three findings; each is fixed in
the re-delivered patch and pinned here by a test that is RED on the buggy code
and GREEN once the fixed patch is built:

  1. The whole-request guard was a no-op: it compared a TOKEN estimate
     (``len(content) // 4``) to the window, then handed that same number to
     ``in_turn_compaction.compact`` as a *char* budget. But ``compact``
     measures tool-CHAR content, so a token-sized budget was always larger
     than the (small) tool output and nothing was ever summarised — the guard
     never fired on the assistant-heavy turns it exists for. The fix measures
     the whole request in chars against the window's char capacity and gives
     ``compact`` a real tool-char budget that lies under the current tool
     chars, so it provably fires.
  2. ``system_chars`` is counted in the whole-request total when the rows
     carry no system/developer row. The client loop carries the prompt as
     row 0 and passes no ``system_chars`` (see
     test_u11_whole_request_is_counted_once.py).
  3. The floor throttle forbade any floor above the static table window from
     being learned, so ``131072 -> 166379`` (the exact case observed_window.py
     exists to capture) was never recorded. The fix LEARNS the floor
     unconditionally — ``observed_window.note`` already bounds the raise to
     the max actually accepted, so it is bounded, not forbidden.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import api_client
from delfin.agent import engine
from delfin.agent import observed_window as OW


@pytest.fixture
def home(tmp_path, monkeypatch):
    """Isolate the persisted observed-window store (lives under Path.home())."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(Path, "home", classmethod(lambda _cls: tmp_path))
    return tmp_path


_WINDOW = 131_072          # table window for kit.deepseek-v4-flash
_CAPACITY = _WINDOW * 4    # char capacity (4 chars/token, as the tool budget)
_BUDGET = int(_CAPACITY * 0.95)  # ≈ 498073 chars


def _messages(*, n_rounds: int = 4, tool_chars: int = 0,
              assistant_chars: int = 0) -> list[dict]:
    """A realistic turn for _compact_between_rounds: one goal + ``n_rounds``
    complete tool rounds (each an assistant call message + its matching tool
    result). ``in_turn_compaction.compact`` only fires when there are STRICTLY
    MORE than ``keep_last_rounds`` (2) rounds to drop, so a single-round toy
    turn would always no-op and prove nothing -- the guard's job is to summarise
    OLD rounds of a long turn, and these complete rounds model that."""
    msgs = [{"role": "user", "content": "goal", "tool_calls": None}]
    for i in range(n_rounds):
        msgs.append({"role": "assistant",
                     "content": "x" * assistant_chars,
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": "{}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}",
                     "content": "x" * tool_chars})
    return msgs


def test_whole_request_guard_fires_when_assistant_dominates_even_with_tiny_tool():
    """Finding 1: assistant-heavy turn with tiny tool output must compact.

    On the buggy patch the guard handed the tool-only compactor a token-sized
    budget (window*0.95 / 4 ≈ 124518) as its tool-CHAR budget; the tool output
    is far under that, so compact() returned None -- the guard was a silent
    no-op and the turn kept its overflowing assistant text.
    On the fix the whole request is measured in chars (4*170000 + 4*100 > 498073)
    and the tool budget is set to `budget - non_tool` (1), which is under the
    400 tool chars present, and there are 4 complete rounds (> 2 kept), so
    compact must summarise the older rounds.
    """
    msgs = _messages(n_rounds=4, assistant_chars=170_000, tool_chars=100)
    note = api_client._compact_between_rounds(
        msgs, 200_000, 0, "s", planning_window=_WINDOW, system_chars=0)
    assert note, "whole-request guard must fire on an assistant-heavy turn"
    # The whole-request branch reports "summarised" (its note names the window
    # crossing, not the tool compactor's wording) -- assert on that, which is
    # what the fix actually returns.
    assert "summarised" in note.lower()


def test_whole_request_guard_counts_the_system_prompt_in_chars():
    """Finding 2: a large system prompt alone must push the guard to fire.

    The buggy guard read ``system_chars`` but then token-estimated it
    (``// 4``) and handed a token number to the char-measuring compactor, so
    the system prompt was never actually booked against the window. The fix
    counts the raw ``system_chars`` toward the whole-request char total.
    """
    # Whole = 4 (goal) + (1000 assistant + 2000 tool) * 4 rounds = 12004;
    #        + system 510000 = 522004 > 498073 -> guard fires, and the tool
    #        budget drops below the 8000 tool chars present -> compacts.
    msgs = _messages(n_rounds=4, assistant_chars=1_000, tool_chars=2_000)
    note = api_client._compact_between_rounds(
        msgs, 200_000, 0, "s",
        planning_window=_WINDOW, system_chars=510_000)
    assert note, "guard must book the system prompt against the window"


def test_floor_learned_even_when_request_exceeds_the_table_window(home):
    """Finding 3: any accepted request is proof of the window.

    The floor must be learned even when the request is ABOVE the static table
    window (166379 > 131072) -- that is exactly the case observed_window.py
    documents. The buggy patch gated learning on ``within_plan`` and dropped
    the 166379 floor; the fix learns it unconditionally, and
    observed_window.note bounds the raise (stores only the max accepted).
    """
    model = "kit.deepseek-v4-flash"
    eng = engine.AgentEngine.__new__(engine.AgentEngine)
    eng.context_window_tokens = _WINDOW
    eng.client = SimpleNamespace(model=model)
    eng._floor_captured_this_turn = False
    eng._last_input_tokens = 0
    eng._trimmed_chars_since_floor = 0

    eng._capture_floor_if_in_plan(166_379)

    floor = OW.at_least(model, _WINDOW)
    assert floor >= 166_379, (
        f"floor must be learned even above the table window (got {floor})")
