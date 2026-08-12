"""A guard that any phrasing can switch off is not a guard.

Every disarm in the output guardrails used to be searched over the WHOLE
answer, so a word that had nothing to do with the claim silenced it. All
four pairs below were produced by executing the shipped functions:

    "Gesamtsumme der Spalte Anschaffungswert: 45.231,50 EUR."   flagged
    ... + " Soll ich noch nach Kostenstelle gruppieren?"        NOT flagged
    "The S1 energy is 2.31 eV."                                 flagged
    "The S1 energy is 2.31 eV. The geometry is around ..."      NOT flagged
    "Das Skript app.py läuft jetzt fehlerfrei."                 flagged
    ... + ", wenn Python 3.11 installiert ist."                 NOT flagged
    "31 PDF-Dateien" over a 29-item list                        flagged
    "31 Rechnungen" over the same list                          NOT flagged

Nothing about the claims had changed. The rule now is that a disarm
excuses the claim it stands next to: same sentence for a hedge or a
question, same-or-earlier clause for a condition, and a general noun
shape instead of a closed list nobody agreed to write in.

The second half of each section is the other half of the deal: a guard
that over-fires teaches the model to write around it, so every widening
here is paired with a correct answer that must stay untouched.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg


_AMBIGUOUS = ["Anschaffungswert"]
_TOTAL = "Gesamtsumme der Spalte Anschaffungswert: 45.231,50 EUR."


# ---------------------------------------------------------------------------
# A question somewhere else in the answer is not a question about the figure
# ---------------------------------------------------------------------------

def test_a_closing_question_does_not_excuse_a_total():
    with_offer = _TOTAL + " Soll ich noch nach Kostenstelle gruppieren?"
    assert vg.scan_for_totals_over_ambiguous_columns(
        _TOTAL, _AMBIGUOUS) == ["Anschaffungswert"]
    assert vg.scan_for_totals_over_ambiguous_columns(
        with_offer, _AMBIGUOUS) == ["Anschaffungswert"]


def test_naming_the_ambiguity_beside_the_figure_still_passes():
    for text in (
        "Gesamtsumme der Spalte Anschaffungswert: 45.231,50 EUR — welche "
        "Lesart gilt?",
        "Die Spalte Anschaffungswert ist mehrdeutig. Die Summe wäre "
        "45.231,50 EUR.",
        "Summe der Spalte Anschaffungswert: 45.231,50 EUR oder 45231,50 EUR.",
    ):
        assert vg.scan_for_totals_over_ambiguous_columns(
            text, _AMBIGUOUS) == [], text


def test_both_readings_of_one_figure_are_recognised():
    assert vg._states_both_readings("8.986 = 8986") is True
    assert vg._states_both_readings("25.136 EUR bzw. 25,136 EUR") is True
    # Two unrelated numbers are not two readings of one.
    assert vg._states_both_readings("45.231,50 EUR aus 31 Zeilen") is False


# ---------------------------------------------------------------------------
# A hedge belongs to its own sentence
# ---------------------------------------------------------------------------

def test_a_hedge_one_sentence_away_does_not_excuse_a_quantity():
    bare = "The S1 energy is 2.31 eV."
    elsewhere = "The S1 energy is 2.31 eV. The geometry is around the minimum."
    assert [f.quantity for f in vg.scan_for_unsourced_quantities(bare)] == [
        "2.31 eV"]
    assert [f.quantity for f in
            vg.scan_for_unsourced_quantities(elsewhere)] == ["2.31 eV"]


def test_a_hedge_on_the_claim_itself_still_excuses_it():
    for text in (
        "The S1 energy is around 2.31 eV.",
        "Die S1-Energie liegt bei ca. 2.31 eV.",
        "Vermutlich sind es 2.31 eV, geprüft habe ich es nicht.",
        "I have not checked this: the S1 energy is 2.31 eV.",
    ):
        assert vg.scan_for_unsourced_quantities(text) == [], text


def test_a_hedge_one_sentence_away_does_not_excuse_a_location():
    bare = "Zeile 26: class AgentEngine"
    elsewhere = "Zeile 26: class AgentEngine. Der Rest ist vermutlich anders."
    assert len(vg.scan_for_ungrounded_location_claims(
        bare, observed_files=set())) == 1
    assert len(vg.scan_for_ungrounded_location_claims(
        elsewhere, observed_files=set())) == 1


def _sentence_around(text: str, needle: str) -> str:
    start = text.index(needle)
    return vg._claim_sentence(text, start, start + len(needle))


def test_an_abbreviation_never_splits_a_sentence():
    """"ca." ends in a period, and splitting there would put the hedge in a
    different sentence from the number it hedges — the exact false positive
    this guard must not produce."""
    assert _sentence_around("Der Wert liegt bei ca. 2.31 eV.", "2.31 eV") == (
        "Der Wert liegt bei ca. 2.31 eV.")
    # A real sentence break does split.
    assert _sentence_around(
        "Das Ergebnis ist 2.31 eV. Der Rest ist etwa gleich.",
        "2.31 eV") == "Das Ergebnis ist 2.31 eV"


def test_a_file_name_never_splits_a_sentence():
    assert _sentence_around(
        "Das Skript app.py läuft fehlerfrei.", "läuft") == (
        "Das Skript app.py läuft fehlerfrei.")


# ---------------------------------------------------------------------------
# A trailing condition names a precondition; it does not retract the claim
# ---------------------------------------------------------------------------

def test_a_trailing_condition_does_not_retract_a_functional_claim():
    bare = "Das Skript app.py läuft jetzt fehlerfrei."
    trailing = ("Das Skript app.py läuft jetzt fehlerfrei, wenn Python 3.11 "
                "installiert ist.")
    assert [f.kind for f in vg.scan_for_unexercised_functional_claims(
        bare, exec_commands=[])] == ["unexercised"]
    assert [f.kind for f in vg.scan_for_unexercised_functional_claims(
        trailing, exec_commands=[])] == ["unexercised"]


def test_a_leading_condition_still_means_it_is_no_assertion():
    for text in (
        "Wenn Python 3.11 installiert ist, läuft das Skript app.py "
        "fehlerfrei.",
        "Damit app.py fehlerfrei läuft, brauchst du Python 3.11.",
        "Sobald du es startest, läuft app.py fehlerfrei.",
        "Der Test prüft, ob app.py fehlerfrei läuft.",
    ):
        assert vg.scan_for_unexercised_functional_claims(
            text, exec_commands=[]) == [], text


def test_an_honest_disclosure_holds_wherever_it_stands():
    """A disclosure is the phrasing the framework asks for. Demanding it
    come FIRST would punish the answer it is trying to elicit."""
    for text in (
        "app.py läuft fehlerfrei, getestet habe ich es allerdings nicht.",
        "Ich konnte nicht verifizieren, dass app.py fehlerfrei läuft.",
        "app.py works now, though I could not verify it.",
        "Wie du bestätigt hast, läuft app.py jetzt fehlerfrei.",
    ):
        assert vg.scan_for_unexercised_functional_claims(
            text, exec_commands=[]) == [], text


def test_a_negation_in_a_later_clause_does_not_clear_the_claim():
    flags = vg.scan_for_unexercised_functional_claims(
        "Das Skript app.py läuft jetzt fehlerfrei, Fehler gibt es keine mehr.",
        exec_commands=[])
    assert [f.kind for f in flags] == ["unexercised"]
    # while a negation ON the predicate still does
    assert vg.scan_for_unexercised_functional_claims(
        "Das Skript app.py läuft nicht fehlerfrei.", exec_commands=[]) == []


# ---------------------------------------------------------------------------
# Counting things does not depend on which noun the user happens to use
# ---------------------------------------------------------------------------

_LIST_OF_29 = "\n".join(f"{i}. Beleg_{i}" for i in range(1, 30))


@pytest.mark.parametrize("noun", [
    "PDF-Dateien", "Rechnungen", "Belege", "Formulare", "Datensätze",
    "Positionen", "Dateien", "Vorgänge", "invoices", "receipts",
])
def test_every_countable_noun_is_treated_the_same(noun):
    text = f"Ich habe 31 {noun} verifiziert.\n{_LIST_OF_29}"
    assert vg.scan_for_counts_over_truncated_output(text, ["list_files"])
    assert vg.scan_for_count_vs_enumeration(text) == [(31, 29)]


@pytest.mark.parametrize("phrase", [
    "Der Lauf dauerte 45 Minuten.",
    "Die Datei ist 12000 Zeichen lang.",
    "Der Anteil liegt bei 31 Prozent.",
    "Die Messung ergab 25 Grad.",
    "Das kostet 45 Euro.",
    "Die Datei ist 20 MB groß.",
    "Der Wert war 31 mal höher.",
    "Version 3.11 ist installiert.",
])
def test_measures_are_not_counts_of_things(phrase):
    """Durations, sizes, money and percentages are not enumerations, and
    caveating them would train the model to avoid stating them."""
    assert vg.scan_for_counts_over_truncated_output(
        phrase, ["list_files"]) == [], phrase


def test_a_singular_noun_is_not_a_count():
    assert vg.scan_for_counts_over_truncated_output(
        "Seite 31 Absatz 2 beschreibt es.", ["list_files"]) == []


def test_the_reported_claim_names_the_whole_phrase():
    """"31 neue Dateien" reads back as itself, not as "31 neue"."""
    assert vg.scan_for_counts_over_truncated_output(
        "Ich habe 31 neue Dateien angelegt.", ["list_files"]) == [
        "31 neue Dateien"]


def test_nothing_truncated_still_means_no_count_caveat():
    text = f"Ich habe 31 Rechnungen verifiziert.\n{_LIST_OF_29}"
    assert vg.scan_for_counts_over_truncated_output(text, []) == []


def test_a_count_matching_its_list_is_not_a_contradiction():
    text = f"Ich habe 29 Rechnungen verifiziert.\n{_LIST_OF_29}"
    assert vg.scan_for_count_vs_enumeration(text) == []


# ---------------------------------------------------------------------------
# The marker that says a correction verified something
# ---------------------------------------------------------------------------

def test_the_marker_names_what_was_read():
    marker = vg.verification_marker({"calc/tddft.out", "delfin/engine.py"})
    assert marker.startswith("\n\n[verify] Self-check")
    assert "tddft.out" in marker
    assert "engine.py" in marker


def test_no_new_file_means_no_marker():
    assert vg.verification_marker(set()) == ""
    assert vg.verification_marker(None) == ""
