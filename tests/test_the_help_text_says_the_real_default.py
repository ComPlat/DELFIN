"""The Tool-rounds help text says what actually happens.

Field report 2026-10-06: a fresh installation stopped after 20 tool
calls while the help text under the control read "Default 500". The
reader had no way to find the real number from the dashboard, and the
one printed was off by a factor of 25.

The text was true until 2026-08-28. Before that commit the control
showed 500 for an unset setting and wrote it back on save, so 500 WAS
what everybody had -- including the installations that had never chosen
it. The fix made the field default to -1 ("let the profile decide") and
left the sentence behind.

What is true now:

  * unset (the field shows -1) -> the per-model profile decides;
  * KIT models carry the lowest budgets, Sonnet/Opus/GPT higher ones;
  * 500 is only the last-resort fallback when the profile registry
    cannot be read at all;
  * 0 disables the cap, leaving the cost circuit-breaker and the
    repeated-error abort as the only stops.

This file pins the text against the mechanism rather than against a
number: a profile may be re-measured, and a help text that hardcodes
today's figure would be wrong again the next time one is.
"""

from __future__ import annotations

import pathlib
import re

_TAB = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_settings.py")


def _help_text() -> str:
    """The paragraph rendered under the Tool-rounds control."""
    src = _TAB.read_text(encoding="utf-8")
    start = src.index("max_tool_rounds_input, max_tool_rounds_hint")
    end = src.index("Memory store key", start)
    block = src[start:end]
    return " ".join(re.findall(r"'([^']*)'", block))


def test_the_text_does_not_claim_a_numeric_default():
    text = _help_text()
    assert "Default 500" not in text, (
        "the control defaults to -1 (per-model); printing 500 is what sent "
        "a user looking for a limit they could not find")
    assert not re.search(r"\bDefault\s+\d+", text), (
        "any fixed number here is a figure that drifts when a profile is "
        "re-measured")


def test_the_text_names_the_three_states_the_control_has():
    """-1, a number, and 0 mean three different things, and the person
    setting it has to be able to tell which one they are choosing."""
    text = _help_text().lower()
    # Any wording that hands the decision to the profile. Checked as
    # meaning rather than as one phrase: the text reads "leaves it to the
    # model profile", which says it better than the "per-model" this
    # assertion first demanded, and a test that pins my wording instead
    # of the meaning would have blocked the clearer sentence.
    assert any(k in text for k in ("per-model", "per model", "model profile")), (
        "an unset control hands the decision to the model profile; the "
        "text has to say so")
    assert "uncapped" in text or "no cap" in text
    assert "continue" in text, "it must say what happens when the cap is hit"


def test_the_text_says_a_small_model_gets_less():
    """The whole surprise was that the number depends on the model. A
    reader who knows that looks in the right place."""
    text = _help_text().lower()
    assert "kit" in text or "smaller model" in text or "small model" in text, (
        "the text must point at what decides the number when it is unset")


def test_the_hint_beside_the_field_still_states_the_two_sentinels():
    src = _TAB.read_text(encoding="utf-8")
    assert "per-model default" in src
    assert "uncapped" in src
