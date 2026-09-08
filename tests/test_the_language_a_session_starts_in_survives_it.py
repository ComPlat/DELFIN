"""The session language, at both ends of a save.

Two faults, one in the detector and one in the record.

**Set.** `_note_session_language` asks `detect_language`, which refuses
below twelve words. That floor is right where it was written — judging
an ANSWER, where a wrong verdict forces a correction turn on something
that was fine. For the pin it is a second guard on top of one that
already works: three function-word hits for one language with a 2×
margin over the other. "Merk dir: die Kennzahl ist 5309. Antworte nur
mit OK." is nine words with four German hits and zero English, and the
length gate refused it anyway.

Measured over the 85 benchmark prompts: 54 pinned at twelve, 58 at
eight, and not one English verdict on a German prompt at any floor.

**Kept.** The field dumped with `str`, so an unset pin was written as the
string "None" — which is TRUTHY. A resumed session then had a pin that
named no language: `_note_session_language` saw it as already set and
never looked again, and the block that renders it looked up "None".
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as VG


@pytest.mark.parametrize("text", [
    "Merk dir: die Kennzahl ist 5309. Antworte nur mit OK.",
    "stell den orca functional auf b3lyp ein",
    "Wie hoch ist der Saldo in Gutschriften.xlsx?",
])
def test_a_short_but_unambiguous_opener_pins_the_session(text):
    assert VG.detect_language(text, min_words=8) == "de", text


def test_judging_an_answer_keeps_the_stricter_floor():
    """The default is unchanged: a wrong verdict there costs a correction
    turn on an answer that was fine."""
    assert VG.detect_language("Merk dir: die Kennzahl ist 5309.") == ""


@pytest.mark.parametrize("text", [
    "wechsel zu Submit",
    "/tab calc",
    "b3lyp def2-TZVP",
    "42",
])
def test_a_fragment_still_says_nothing(text):
    assert VG.detect_language(text, min_words=6) == ""


def test_the_margin_rule_still_applies_below_the_floor():
    """Technical German quotes English identifiers and vice versa, and
    the 2x margin is what stops one borrowed word deciding. Lowering the
    length floor does not touch it.

    A genuine 3:1 majority IS a verdict and always was — what the margin
    refuses is a near-tie, and a near-tie is refused at eight words for
    the same reason it is at twelve.
    """
    # near-tie: nothing is claimed
    assert VG.detect_language(
        "the input für die calculation is der test", min_words=6) == ""
    # clear majority: claimed, at either floor
    mostly_german = "the builder shows den input für die calculation"
    assert VG.detect_language(mostly_german, min_words=6) == "de"


# ---------------------------------------------------------------------------
# ...and it has to survive the save
# ---------------------------------------------------------------------------

def _field(name: str):
    from delfin.agent.engine import AgentEngine
    for f in AgentEngine._SESSION_FIELDS:
        if f.attr == name:
            return f
    raise AssertionError(f"no session field {name}")


def test_an_unset_pin_is_not_written_as_the_word_none():
    """`str(None)` is "None", and "None" is truthy. A resumed session had
    a pin that named no language and would never be set again."""
    f = _field("_session_language")
    assert f.dump(None) == ""
    assert f.load("") == ""


def test_a_real_pin_round_trips():
    f = _field("_session_language")
    assert f.load(f.dump("de")) == "de"
    assert f.load(f.dump("en")) == "en"


def test_a_resumed_session_with_no_pin_can_still_be_pinned(tmp_path):
    from unittest.mock import MagicMock, patch

    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                          model="kit.glm-5.3", mode="solo")
    eng._session_language = _field("_session_language").load(
        _field("_session_language").dump(None))
    eng._note_session_language(
        "Bau mir bitte ein kleines Werkzeug, das die Tags aus der Datei "
        "zaehlt und mir das Ergebnis nennt.")
    assert eng._session_language == "de"
