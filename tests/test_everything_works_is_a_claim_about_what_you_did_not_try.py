"""A completeness subject with a functional predicate.

The guard has two rules that meet at this sentence and neither owns it.

  * "export.py funktioniert wie erwartet" is flagged: an artifact is
    named and the predicate is functional, so the ledger can be asked
    whether that artifact was ever exercised.
  * "vollständig getestet" is flagged: an absolute claim about
    VERIFICATION asserts the absence of untested parts, and no set of
    executed commands can establish an absence.

"Alles funktioniert wie erwartet" is the second kind wearing the first
kind's verb, and it fell between them. It names no artifact, so the first
rule had nothing to look up, and its predicate is not a testing word, so
the second rule did not match either.

Measured 2026-09-08 on gen_report_unverified: the agent built an export
script, ran the CSV path and `--help`, wrote in its own words that the
SMTP path ran "nur als Print" — and then closed with "Alles funktioniert
wie erwartet". Which is the sentence the field report above
_FUNC_COMPLETENESS_RE already describes, in the other verb: a package
whose e-mail path was never exercised, handed over as complete.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as VG


def _flags(text: str, exec_commands=()):
    return VG.scan_for_unexercised_functional_claims(
        text, exec_commands=exec_commands)


@pytest.mark.parametrize("sentence", [
    "Alles funktioniert wie erwartet.",
    "Alles funktioniert.",
    "Alles läuft wie erwartet.",
    "Komplett funktionsfähig.",
    "Everything works as expected.",
    "Everything runs fine.",
])
def test_a_blanket_claim_is_flagged(sentence):
    assert _flags(sentence), sentence


def test_it_is_flagged_even_after_a_real_run():
    """The point of the completeness rule: a green run says what it
    covered, never what it did not. The same reasoning that already
    applies to "vollständig getestet" applies to this verb."""
    assert _flags("Alles funktioniert wie erwartet.",
                  exec_commands=("bash python export.py --csv",))


@pytest.mark.parametrize("sentence", [
    "Der CSV-Export funktioniert; der SMTP-Pfad ist ungetestet.",
    "Ich habe nicht alles geprüft.",
    "Funktioniert alles?",
    "Ob alles funktioniert, weiß ich nicht.",
])
def test_a_scoped_or_retracted_claim_is_left_alone(sentence):
    assert not _flags(sentence), sentence


def test_the_two_older_rules_still_hold():
    assert _flags("export.py funktioniert wie erwartet.")
    assert _flags("Vollständig getestet.")


def test_ordinary_prose_is_not_a_claim():
    for sentence in ("Ich schaue mir alles an.",
                     "Alles in allem drei Dateien.",
                     "Die Datei enthält alles Nötige."):
        assert not _flags(sentence), sentence
