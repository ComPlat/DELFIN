"""Two defects in how a German figure is judged, and they interact.

**The total cue missed German inflection.** The figure guard reads a
sentence for a claim about a total and checks whether any tool produced
the number. Measured on the shipped cue list, with an empty ledger and a
fabricated amount:

    flags=0 | Die Kosten betragen 99.999,99 EUR.
    flags=0 | Die Summen der Konten betragen 99.999,99 EUR.
    flags=0 | Die Personalkosten liegen bei 99.999,99 EUR.
    flags=1 | Die Gesamtsumme beträgt 99.999,99 EUR.

The list had the third-person singular ``beträgt`` and the compound noun
``Gesamtsumme``, and neither the plural ``betragen`` — the most ordinary
way German states a total — nor the plural ``Summen``. The comment above
it claimed compounding was the hard part; inflection was.

**The hedges disagreed with each other.** Two guards read the same
sentence and only one of them saw the hedge:

    hedged=False | Die Summe liegt bei rund 45.000 EUR.
    hedged=False | Ich schätze die Summe auf 45.000 EUR.
    hedged=False | Ich glaube, es sind 45.000 EUR.
    hedged=False | Soweit ich weiß sind es 45.000 EUR.

The answer guard had a full English first-person vocabulary ("i think",
"from memory", "iirc") and covered German with adverbs only. And ``rund``
WAS in the office guard's non-assertion set — so one sentence was an
honest hedge to one guard and a confident claim to the other, inside one
answer.
"""

from __future__ import annotations

import pytest

from delfin.agent import german as G
from delfin.agent import office as O
from delfin.agent import verify_guard as VG


def _flags(sentence):
    return O.scan_answer_for_unledgered_figures(sentence, ledger=[])


# ---------------------------------------------------------------------------
# A fabricated total, said the way German says it
# ---------------------------------------------------------------------------

_FABRICATED = (
    "Die Kosten betragen 99.999,99 EUR.",
    "Die Summen der Konten betragen 99.999,99 EUR.",
    "Die Personalkosten liegen bei 99.999,99 EUR.",
    "Die Gesamtsumme beträgt 99.999,99 EUR.",
    "Die Ausgaben belaufen sich auf 99.999,99 EUR.",
    "Die Reisekosten lagen bei 99.999,99 EUR.",
    "Der Aufwand betrug 99.999,99 EUR.",
    "Die Beträge der Belege betragen 99.999,99 EUR.",
    "Die Einnahmen liegen bei 99.999,99 EUR.",
    "Der Endbetrag ist 99.999,99 EUR.",
)


@pytest.mark.parametrize("sentence", _FABRICATED)
def test_a_total_nothing_produced_is_flagged_however_german_states_it(
        sentence):
    assert _flags(sentence), sentence


@pytest.mark.parametrize("sentence", _FABRICATED)
def test_the_same_total_is_silent_once_a_tool_produced_it(sentence):
    """The widening must not turn the guard into noise: a figure the
    ledger holds is an ordinary answer."""
    ledger = [O.ToolFigure(value=99999.99, kind="total", label="Summe",
                           tool="sum_column")]
    assert not O.scan_answer_for_unledgered_figures(sentence, ledger=ledger), \
        sentence


@pytest.mark.parametrize("sentence", [
    "Die Rechnung steht in Zeile 4711.",
    "Der Beleg trägt die Nummer 99999.",
    "Die Datei hat 3 Spalten.",
    "Version 2.31 ist installiert.",
])
def test_a_sentence_that_states_no_total_is_left_alone(sentence):
    assert not _flags(sentence), sentence


# ---------------------------------------------------------------------------
# One hedge vocabulary, two guards
# ---------------------------------------------------------------------------

_HEDGED = (
    "Die Summe liegt bei rund 45.000 EUR.",
    "Ich schätze die Summe auf 45.000 EUR.",
    "Ich glaube, es sind 45.000 EUR.",
    "Soweit ich weiß sind es 45.000 EUR.",
    "Ich nehme an, es sind 45.000 EUR.",
    "Ich vermute, die Summe beträgt 45.000 EUR.",
    "Ich gehe davon aus, dass es 45.000 EUR sind.",
    "Meines Wissens beträgt die Summe 45.000 EUR.",
    "Wenn ich mich recht erinnere, sind es 45.000 EUR.",
    "Aus dem Gedächtnis: die Summe beträgt 45.000 EUR.",
    "Die Summe beträgt ca. 45.000 EUR.",
    "Die Summe beträgt ungefähr 45.000 EUR.",
    "Die Summe wäre dann 45.000 EUR.",
    "Die Summe ist vermutlich 45.000 EUR.",
    "Die Summe ist nicht geprüft: 45.000 EUR.",
    "I think the total is 45,000 EUR.",
    "From memory, the total is 45,000 EUR.",
    "The total is roughly 45,000 EUR.",
)

_ASSERTED = (
    "Die Summe beträgt 45.000 EUR.",
    "Die Kosten betragen 45.000 EUR.",
    "Die Gesamtsumme der Belege ist 45.000 EUR.",
    "Der Saldo liegt bei 45.000 EUR.",
    "The total is 45,000 EUR.",
)


@pytest.mark.parametrize("sentence", _HEDGED)
def test_both_guards_call_the_same_sentence_hedged(sentence):
    """The contract that was broken: one answer, one verdict. It is the
    same regex object in both places, and this asserts the wiring."""
    assert VG._is_hedged(sentence, 0, len(sentence)), sentence
    assert O._NOT_ASSERTED_RE.search(sentence), sentence


@pytest.mark.parametrize("sentence", _ASSERTED)
def test_both_guards_call_the_same_sentence_asserted(sentence):
    assert not VG._is_hedged(sentence, 0, len(sentence)), sentence
    assert not O._NOT_ASSERTED_RE.search(sentence), sentence


def test_the_two_guards_share_one_vocabulary():
    assert O._NOT_ASSERTED_RE is G.HEDGE_RE
    assert VG._HEDGE_MARKERS is G.HEDGE_RE


@pytest.mark.parametrize("sentence", _HEDGED)
def test_a_hedged_figure_is_not_flagged_as_fabricated(sentence):
    """The two halves together: an answer that says it is guessing has
    already disclosed its grounding, and caveating it would punish the
    honest wording."""
    assert not _flags(sentence), sentence
