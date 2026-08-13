"""A German write task was signed off with nothing written.

``check_completion_claim`` decides what a finished task can show for
itself. Measured on the shipped matchers, the same task passed in German
and failed in English:

    verified  path_read      | Trage die Werte in Buchungen.csv ein
    verified  path_read      | Übertrage die Beträge nach Journal.xlsx
    unmet     path_unwritten | Write the values into Buchungen.csv

with nothing written in any of the three. Two causes:

* German puts the prefix of a separable verb at the end of the clause, so
  a pattern written for ``anpass\\w*`` cannot see "Passe die Tabelle an";
* the German half of the write vocabulary was a CODING vocabulary
  (``erstell|implementier|refaktor``) with no office verb in it at all.

Ten of ten realistic German office subjects failed.

The second half of this file is the reverse defect in the same function:
the artefact branch matched its nouns as SUBSTRINGS and ran before any
read/write distinction, so "Summarize the findings briefly" promised a
letter ("brief" inside "briefly") and "Prüfe den Bericht auf Fehler" — a
read-only task — was unmet for having written no PDF.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import (
    _READ_VERB_RE, _WRITE_VERB_RE, _artifact_word, check_completion_claim,
)


# ---------------------------------------------------------------------------
# The write vocabulary: German office work, however the verb is split
# ---------------------------------------------------------------------------

_WRITE_SUBJECTS = (
    # verb-second: the prefix is at the end of the clause
    "Trage die Werte in Buchungen.csv ein",
    "Trag das bitte in die Liste ein",
    "Passe die Tabelle an",
    "Lege ein neues Tabellenblatt an",
    "Fülle das Formular aus",
    "Füge die fehlende Zeile hinzu",
    "Gib die Beträge in die Maske ein",
    "Lege die Belege im Archiv ab",
    "Buche die Rechnung auf 4711 um",
    "Ordne die Belege den Kostenstellen zu",
    "Trage die Summe nach und speichere die Datei",
    # prefix attached, and the perfect
    "Die Werte eintragen",
    "Das Blatt anlegen",
    "Die Beträge übertragen",
    "Die Rechnung verbuchen",
    "Das Formular ausfüllen",
    "Die Liste sortieren",
    "Die Zeile ist bereits eingetragen",
    "Das Blatt wurde angelegt",
    # inseparable office verbs
    "Übertrage die Beträge nach Journal.xlsx",
    "Verbuche die Rechnung auf der Kostenstelle",
    "Sortiere die Liste nach Datum",
    "Storniere die doppelte Buchung",
    "Beschrifte die Spalten neu",
    # Found by running the SHIPPED matcher over a corpus written from the
    # intent rather than from the table. Each of these was read as a task
    # that wants no write, so a session that had merely opened the file
    # reported it verified.
    "Setze das Datum auf den 31.07.2026",
    "Hinterlege die Bankverbindung im Stammblatt",
    "Streiche die Position 7 aus dem Angebot",
    # A rename is a write, and only the joined form was ever matched.
    "Benenne die Datei um",
    "Benenne das Tabellenblatt in Q3 um",
)


@pytest.mark.parametrize("subject", _WRITE_SUBJECTS)
def test_a_german_office_task_promises_a_change(subject):
    assert _WRITE_VERB_RE.search(subject), subject


_NOT_WRITE_SUBJECTS = (
    # reading and computing: the answer is a number in the chat
    "Prüfe den Bericht auf Fehler",
    "Analysiere die Tabelle",
    "Rechne die Summe der Spalte aus",
    "Werte die Belege aus",
    "Fasse den Bericht zusammen",
    "Sieh dir die Tabelle an",
    "Gleiche die beiden Listen ab",
    "Vergleiche die Beträge",
    "Lies die Datei",
    "Zähle die offenen Posten",
    # sentences whose words merely LOOK like a split verb
    "Ich trage ein Kleid",
    "Die Angelegenheit ist erledigt",
    "Das passiert an manchen Tagen",
    "Ein Beispiel für eine Tabelle",
    "Das Ergebnis liegt an der Konfiguration",
    # "setzen" is the finite half of several verbs and only one of them
    # writes. The prefix at the end of the clause is what tells them apart.
    "Setze die Recherche fort",
    "Das ist eine Voraussetzung",
    "Setze die Auswertung morgen fort",
)


@pytest.mark.parametrize("subject", _NOT_WRITE_SUBJECTS)
def test_a_task_that_changes_nothing_is_not_read_as_a_write(subject):
    """A matcher that fires on everything is worse than one that fires on
    nothing, because the user stops reading it."""
    assert not _WRITE_VERB_RE.search(subject), subject


@pytest.mark.parametrize("subject", [
    "Rechne die Summe der Spalte aus",
    "Werte die Belege aus",
    "Fasse den Bericht zusammen",
    "Sieh dir die Tabelle an",
    "Gleiche die beiden Listen ab",
    "Sieh in Buchungen.csv nach",
    # Read verbs the tables did not have. A read task that matches no read
    # verb is not merely unrecognised: the artefact branch then runs on it
    # and reports it unmet for having produced no document.
    "Zeig mir die Umsatzliste",
    "Öffne die Kostenstellenübersicht",
    "Nenne mir die zehn größten Posten",
    # A question is a read intent even without a read verb, and it is the
    # most ordinary way somebody asks for a figure.
    "Wie hoch ist die Gesamtsumme?",
    "Wie viele Belege fehlen?",
    "Was steht in der Spalte D?",
    "Berichte mir über die Kosten",
])
def test_a_german_compute_task_is_read_as_reading(subject):
    """"Rechne die Summe aus" is answered by reading and arithmetic. Calling
    it a write would tell a user their honest, finished task wrote
    nothing."""
    assert _READ_VERB_RE.search(subject), subject


# ---------------------------------------------------------------------------
# ... and what the verdict becomes
# ---------------------------------------------------------------------------

def _verdict(subject, *, written=(), read=()):
    changes = [{"path": p, "ts": 1.0, "created": False} for p in written]
    return check_completion_claim(subject, changes=changes, observed=list(read))


@pytest.mark.parametrize("subject", [
    "Trage die Werte in Buchungen.csv ein",
    "Übertrage die Beträge nach Journal.xlsx",
    "Passe Buchungen.csv an",
    # The verbs added after measuring: the point is not that the regex
    # matches, it is that opening the file no longer finishes the task.
    "Setze das Datum in Buchungen.csv auf den 31.07.2026",
    "Hinterlege die Bankverbindung in Journal.xlsx",
    "Streiche die Position 7 in Buchungen.csv",
    "Benenne Buchungen.csv um",
])
def test_a_read_does_not_finish_a_german_write_task(subject):
    assert _verdict(subject, read=["/w/Buchungen.csv", "/w/Journal.xlsx"]
                    )["verdict"] == "unmet"


@pytest.mark.parametrize("subject,path", [
    ("Trage die Werte in Buchungen.csv ein", "/w/Buchungen.csv"),
    ("Übertrage die Beträge nach Journal.xlsx", "/w/Journal.xlsx"),
])
def test_the_write_that_did_happen_finishes_it(subject, path):
    out = _verdict(subject, written=[path])
    assert out["verdict"] == "verified" and out["kind"] == "path_write"


def test_a_german_read_task_is_finished_by_the_read():
    out = _verdict("Prüfe Buchungen.csv auf Fehler",
                   read=["/w/Buchungen.csv"])
    assert out["verdict"] == "verified" and out["kind"] == "path_read"


def test_a_german_compute_task_over_a_file_is_finished_by_the_read():
    out = _verdict("Rechne die Summe in Buchungen.csv aus",
                   read=["/w/Buchungen.csv"])
    assert out["verdict"] == "verified", out


# ---------------------------------------------------------------------------
# The artefact noun: a word, not a substring
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("subject", [
    "Summarize the findings briefly",
    "Debrief the team",
    "Explain it briefly",
    "The consumption is high",
])
def test_an_english_word_that_merely_contains_a_german_noun_promises_nothing(
        subject):
    assert _artifact_word(subject) == "", subject


@pytest.mark.parametrize("subject,word", [
    ("Schreibe einen Brief an den Kunden", "brief"),
    ("Schreibe die Briefe an alle Kunden", "brief"),
    ("Erstelle den Serienbrief", "serienbrief"),
    ("Erstelle das Anschreiben", "anschreiben"),
    ("Erstelle den Bericht", "bericht"),
    ("Erstelle die Berichte", "bericht"),
    ("Erstelle eine Tabelle", "tabelle"),
    ("Erstelle die Tabellen", "tabelle"),
])
def test_the_german_noun_and_its_plural_still_promise_the_artefact(
        subject, word):
    assert _artifact_word(subject) == word, subject


def test_a_read_only_task_is_not_asked_to_produce_an_artefact():
    """"Prüfe den Bericht auf Fehler" names a report and promises to make
    none. The artefact branch used to run before the read/write
    distinction and reported it unmet for having written no PDF."""
    out = check_completion_claim("Prüfe den Bericht auf Fehler",
                                 changes=[], observed=[])
    assert out["verdict"] != "unmet", out
    assert out["kind"] != "artifact", out


def test_a_write_task_that_names_an_artefact_still_needs_one():
    out = check_completion_claim("Erstelle den Bericht als PDF",
                                 changes=[], observed=[])
    assert out["verdict"] == "unmet" and out["kind"] == "artifact"
