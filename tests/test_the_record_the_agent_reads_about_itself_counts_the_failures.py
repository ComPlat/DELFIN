"""The agent was told a pass rate computed only over turns that worked.

``_record_turn_outcome`` has one call site, and it is guarded by

    if chunks and engine.mode in ("solo", "dashboard", "quick", ...):

``chunks`` is the streamed text. So an outcome row is written only when
the turn produced text. A turn killed by the first-token watchdog ("No
tokens were produced"), a turn that raised out of the stream loop, a
silent exit — none of them are recorded. The same turn IS written to
``agent_metrics.jsonl`` with ``silent_exit=True``, so the two logs
disagree about the same event, and only one of them is read back.

That record does not stay in a file. ``briefing.generate_briefing``
writes "Task class: X | History: N outcomes, 96% pass rate" INTO THE
SYSTEM PROMPT, and the provider profile uses the same history for
adaptive routing. So the model is handed, as evidence about itself, a
rate from which every zero-output failure has been removed — and the
worse the backend behaves, the cleaner the record looks.

The mode list has a second hole the audit did not name: the three modes
that exist are dashboard, solo and office, and **office is not in the
tuple**. Every office turn — the mode this user actually works in — has
been recorded nowhere at all.

Both are the same mistake in different clothes: a measurement whose
sampling is decided by the thing being measured.
"""

from __future__ import annotations

import pathlib
import re

import pytest


_TAB = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
        / "dashboard" / "tab_agent.py")


def _call_site() -> str:
    src = _TAB.read_text(encoding="utf-8")
    # The CALL, not the definition: the definition sits earlier in the
    # file and carries none of the gating this test is about.
    i = src.index("                        _record_turn_outcome(")
    return src[max(0, i - 1200):i + 400]


def _recorded_modes() -> set[str]:
    """The modes whose turns reach the outcome log."""
    block = _call_site()
    return set(re.findall(r'"(dashboard|solo|office)"', block))


# ---------------------------------------------------------------------------
# Every mode that exists is recorded
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["dashboard", "solo", "office"])
def test_every_selectable_mode_reaches_the_record(mode):
    assert mode in _recorded_modes(), (
        f"{mode} turns are recorded nowhere, so they cannot appear in the "
        f"pass rate the agent is shown")


# ---------------------------------------------------------------------------
# A turn that produced nothing is a failure, not an absence
# ---------------------------------------------------------------------------

def test_a_turn_without_output_is_still_recorded():
    """`if chunks` made silence unobservable. The row has to be written;
    what it says is the verdict's business, not the gate's."""
    block = _call_site()
    assert not re.search(r"if chunks and", block), (
        "the outcome is still gated on the turn having produced text")


def test_the_verdict_can_see_that_nothing_was_produced():
    from delfin.dashboard.tab_agent import _compute_turn_verdict
    import inspect
    src = inspect.getsource(_compute_turn_verdict)
    assert "error_type" in src


def test_an_empty_answer_is_not_a_pass():
    from delfin.dashboard.tab_agent import _compute_turn_verdict
    verdict = _compute_turn_verdict(
        engine=None, response_text="", denied_commands=[], state={},
        error_type="no_output")
    assert verdict == "FAIL"


def test_a_real_answer_is_still_a_pass():
    from delfin.dashboard.tab_agent import _compute_turn_verdict
    verdict = _compute_turn_verdict(
        engine=None,
        response_text="Done — the file now holds the corrected total.",
        denied_commands=[], state={}, error_type="")
    assert verdict == "PASS"


def test_the_two_logs_no_longer_disagree():
    """agent_metrics records silent_exit for exactly the turns the outcome
    log dropped; both should describe the same event."""
    src = _TAB.read_text(encoding="utf-8")
    assert "silent_exit" in src
    i = src.index("                        _record_turn_outcome(")
    assert "if chunks and" not in src[max(0, i - 1200):i]
