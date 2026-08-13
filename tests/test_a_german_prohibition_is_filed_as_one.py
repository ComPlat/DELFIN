"""The strongest German prohibition was filed as background.

``_classify_by_heuristic`` picks a memory type when the user gives no
prefix. ``feedback`` is the corrective class — the one that is protected
from pruning and shown back to the model. Measured on the shipped hints:

    type='user'     | Committe niemals ohne Test
    type='user'     | Nie ohne Test committen
    type='feedback' | Bitte nicht ohne Test committen

Two causes, both mechanical. ``nie`` and ``niemals`` were absent from the
list outright. And every German hint carried a TRAILING SPACE — "nicht ",
"kein " — so a hint could only match with a word after it, which is
exactly where German does NOT put its negation: at the end of the
sentence.
"""

from __future__ import annotations

import pytest

from delfin.agent.memory_store import _classify_by_heuristic as classify


@pytest.mark.parametrize("text", [
    "Committe niemals ohne Test",
    "Nie ohne Test committen",
    "Bitte nicht ohne Test committen",
    "Committe nicht ohne Test",
    "Das machen wir grundsätzlich nicht",
    "Auf keinen Fall force-pushen",
    "Immer erst die Tests laufen lassen",
    "Bitte stets die Belege prüfen",
    "Vermeide bulk git add",
    "Der Branch darf nie direkt auf main",
    "Never commit without a test",
    "Don't use bulk add",
    "Always run the tests first",
])
def test_a_correction_is_filed_as_feedback(text):
    assert classify(text) == "feedback", text


@pytest.mark.parametrize("text,expected", [
    # The other classes must still win where they should: widening the
    # feedback hints must not swallow every note.
    ("Die Abgabe ist nächste Woche", "project"),
    ("Der Release ist am Freitag", "project"),
    ("Siehe https://example.org/doku", "reference"),
    ("Das Dashboard läuft in Grafana", "reference"),
    # ... and an ordinary note stays an ordinary note.
    ("Der Kunde heißt Müller GmbH", "user"),
    ("Die Kostenstelle für das Institut ist 4711", "user"),
    ("Die Buchungen liegen in Buchungen.xlsx", "user"),
])
def test_what_is_not_a_correction_is_not_filed_as_one(text, expected):
    assert classify(text) == expected, text


@pytest.mark.parametrize("text", [
    "Die Priorität ist niedrig",           # "nie" inside "niedrig"
    "Die Immersion war gut",               # "immer" inside "Immersion"
    "Der Nichtraucherbereich ist hinten",  # "nicht" inside a compound
])
def test_a_hint_inside_a_longer_word_does_not_count(text):
    """Word-bounded, not substring. Dropping the trailing space is what let
    a hint match at the end of a sentence; it must not also let it match
    in the middle of a word."""
    assert classify(text) == "user", text
