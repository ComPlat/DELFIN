"""A session opened with "wie funktioniert aktuell der co2 coordinator"
ran with no language: the detector knew "der" and none of the other
five words, and the session pin needed eight words anyway. The role
answered in English (field report 2026-09-11; "Warum immer noch
englisch und deutsch?" on 2026-09-04 was the same gap). The function
words of a user's opener -- question words, pronouns, "bitte" -- count
now, and the pin takes five words; the three-hit, two-to-one rule is
what keeps a fragment from pinning anything. Measured over the 99
benchmark prompts: 77 pinned German at both floors, one more at five,
none English.
"""

from delfin.agent import verify_guard as vg
from delfin.agent.engine import AgentEngine


def test_short_german_openers_pin_german_and_english_ones_pin_english():
    assert vg.detect_language("wie funktioniert aktuell der co2 coordinator", min_words=5) == "de"
    assert vg.detect_language("kannst du die tabelle prüfen?", min_words=5) == "de"
    assert vg.detect_language("hier ist die datei 07a", min_words=5) == "de"
    assert vg.detect_language("what does the co2 coordinator do?", min_words=5) == "en"


def test_a_fragment_still_pins_nothing():
    for text in ("hallo", "ok start job xyz", "/calc search co2", "PBE0 def2-TZVP CPCM(water)",
                 "run the tests please"):
        assert vg.detect_language(text, min_words=5) == "", text


def test_the_answer_judge_keeps_its_own_higher_floor():
    """Judging an ANSWER's language forces a correction turn when wrong;
    it keeps the twelve-word caution."""
    assert vg.MIN_WORDS_FOR_LANGUAGE == 12
    assert vg.detect_language("kannst du die tabelle prüfen?") == ""


def test_the_session_pin_reads_five_words(tmp_path):
    import inspect
    src = inspect.getsource(AgentEngine._note_session_language)
    assert "min_words=5" in src and "min_words=8" not in src
