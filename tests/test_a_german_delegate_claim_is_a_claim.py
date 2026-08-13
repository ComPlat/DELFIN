"""A delegate reporting in German claimed nothing the guard could see.

``scan_report_completion_claims`` reads a delegate's report and tags what
it claims. The ``mutation`` family — "I changed something" — is the one
that has to be checked against the write ledger. Measured on the shipped
pattern:

    []           | Ich habe die Datei geschrieben.
    []           | Ich habe die Tabelle erstellt.
    []           | Ich habe die Werte eingetragen.
    []           | Ich habe die Beträge übertragen.
    []           | Ich habe ein Blatt angelegt.
    []           | Ich habe die Spalte ergänzt.
    []           | Ich habe die Datei gespeichert.
    ['mutation'] | Ich habe die Datei aktualisiert.

One participle out of eight. The same gap ran the other way in English:
the list carried ``rewrote`` and not ``wrote``, and no ``saved`` at all,
so "I wrote the file" and "I saved the file" were no claim either.
"""

from __future__ import annotations

import pytest

from delfin.agent.subagents import scan_report_completion_claims


def _families(text) -> list[str]:
    return [c.family for c in scan_report_completion_claims(text)]


_GERMAN_MUTATIONS = (
    "Ich habe die Datei geschrieben.",
    "Ich habe die Tabelle erstellt.",
    "Ich habe die Werte eingetragen.",
    "Ich habe die Beträge übertragen.",
    "Ich habe ein Blatt angelegt.",
    "Ich habe die Spalte ergänzt.",
    "Ich habe die Datei gespeichert.",
    "Ich habe die Rechnung verbucht.",
    "Ich habe das Formular ausgefüllt.",
    "Ich habe die Liste sortiert.",
    "Ich habe die Zeile nachgetragen.",
    "Ich habe die Belege abgelegt.",
    "Ich habe die Buchung storniert.",
    "Wir haben die Spalten neu beschriftet.",
    "Ich habe die Datei aktualisiert.",
)

_ENGLISH_MUTATIONS = (
    "I wrote the file.",
    "I saved the file.",
    "I rewrote the file.",
    "I appended the row.",
    "I filled in the form.",
    "We exported the table.",
    "I sorted the list.",
)


@pytest.mark.parametrize("report", _GERMAN_MUTATIONS + _ENGLISH_MUTATIONS)
def test_a_delegate_saying_it_changed_a_file_makes_a_mutation_claim(report):
    assert "mutation" in _families(report), report


@pytest.mark.parametrize("report", [
    # A FINDING about existing code is not a claim of authorship — the
    # first-person frame is what separates the two, and widening the
    # vocabulary must not weaken that.
    "Die Spalte wurde in einer früheren Sitzung eingetragen.",
    "Der Wert steht bereits in der Datei.",
    # Reading, not writing.
    "Ich habe die Datei gelesen.",
    "Ich habe die Tabelle geprüft.",
    "Ich habe die Summe nachgerechnet.",
    # Explicitly hedged or negated: already disclosed.
    "Ich habe vermutlich die falsche Datei gespeichert.",
    "Ich habe die Datei nicht geschrieben.",
])
def test_what_is_not_a_claim_of_authorship_stays_silent(report):
    assert "mutation" not in _families(report), report


def test_a_bare_report_bullet_still_counts():
    """A delegate that drops the pronoun is still reporting its own work."""
    assert "mutation" in _families("- Eingetragen: die Werte aus der Liste")


@pytest.mark.xfail(
    reason="the agentive test asks whether the verb OPENS the line, which "
           "is English word order. German puts the participle last, so "
           "'- Werte eingetragen' still makes no claim. The helper that "
           "decides it is outside this change.",
    strict=True,
)
def test_a_german_bullet_puts_its_verb_last():
    assert "mutation" in _families("- Werte in Buchungen.csv eingetragen")


def test_a_verification_claim_is_still_its_own_family():
    """Widening the mutation family must not swallow the neighbouring one:
    they are judged against different evidence."""
    assert _families("Ich habe die Summe verifiziert.") == ["verification"]
    assert _families("Die Datei ist unverändert.") == ["verification"]
