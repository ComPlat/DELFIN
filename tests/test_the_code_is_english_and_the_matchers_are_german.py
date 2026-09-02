"""Which language a string is in, decided by who reads it.

One German answer about one spreadsheet used to come back with three
caveats — German, English, German — two of them from the SAME file, so
the split was accidental rather than decided. The English one even quoted
the German noun back at the reader ("This answer states 31 Belege").

**The decision taken here** is not "user-facing text is English". It is
who READS the string and whether anything translates it afterwards:

* A tool result, an error, a note the MODEL gets back — ENGLISH. The
  model reads it and answers the user in the user's language, so it is
  translated on the way out. Every message in ``office.py`` is that kind.
* A caveat — GERMAN. It is appended to the FINISHED answer, after the
  model's last turn. Nothing follows it to translate it, and it asks the
  reader to do something ("bitte nachrechnen, bevor die Zahl
  weitergegeben wird"). An instruction the reader skips is not a warning,
  and this framework's users write German.

The ``[verify] Caveat:`` / ``[verify] Self-check:`` tags stay as they
are: they are machine markers the engine and the tests key on.

The second half of the file is a census. The rule this project runs on is
that code and comments are English and only the matchers that recognise
what the USER writes are German. A census that only ever goes up is what
stops the German half being quietly refactored away.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent import german as G
from delfin.agent import office as O
from delfin.agent import verify_guard as VG


# ---------------------------------------------------------------------------
# Every caveat a finished answer can carry, in one language
# ---------------------------------------------------------------------------

def _every_caveat() -> dict[str, str]:
    """Each caveat producer, called with input that makes it speak."""
    figure_flags = O.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 99.999,99 EUR.", ledger=[])
    location_flags = VG.scan_for_ungrounded_location_claims(
        "Zeile 26: class AgentEngine", observed_files=set())
    functional_flags = VG.scan_for_unexercised_functional_claims(
        "Das Paket ist funktionsfähig und vollständig getestet.",
        exec_commands={"bash pytest -q"}, exec_ledger_available=True)
    return {
        "figure": O.figure_caveat(figure_flags),
        "grounding": VG.grounding_caveat(location_flags, []),
        "functional": VG.functional_claim_caveat(functional_flags),
        "count_vs_enumeration": VG.count_vs_enumeration_caveat([(31, 29)]),
        "truncated_output": VG.truncated_output_caveat(
            ["31 Belege"], ["list_files"]),
        "ambiguous_column": VG.ambiguous_column_caveat(["Anschaffungswert"]),
    }


# Words that only appear in a German sentence, and English function words
# that only appear in an English one. Cheap, and enough to tell a
# translated caveat from an untranslated one.
_GERMAN_MARK = re.compile(
    r"(?i)\b(?:die|der|das|und|nicht|bitte|wurde|werden|ist|sind|dieser|"
    r"diese|dieses|eine|einen|kein|keine)\b")
# Only words that cannot also be German: "was" and "die" are both, and a
# census that counts them measures nothing.
_ENGLISH_MARK = re.compile(
    r"(?i)\b(?:the|and|this|these|those|with|from|which|there|treat|"
    r"before|passing|following|claims?|answer|states?|entries|"
    r"unverified|unconfirmed|estimated|counted|truncated)\b")


@pytest.mark.parametrize("name", sorted(_every_caveat()))
def test_every_caveat_reaches_the_reader_in_german(name):
    text = _every_caveat()[name]
    assert text.strip(), f"{name} produced nothing to check"
    assert _GERMAN_MARK.search(text), f"{name}: {text!r}"


@pytest.mark.parametrize("name", sorted(_every_caveat()))
def test_no_caveat_is_half_translated(name):
    """The failure this replaces: one answer, three caveats, two
    languages. The `[verify]` tag is a marker and does not count."""
    text = _every_caveat()[name].replace("[verify] Caveat:", "")
    text = text.replace("[verify] Self-check:", "")
    assert not _ENGLISH_MARK.search(text), f"{name}: {text!r}"


def test_the_machine_markers_are_untouched():
    """The engine, the CLI and the tests key on these."""
    caveats = _every_caveat()
    assert caveats["grounding"].startswith("\n\n[verify] Caveat")
    assert caveats["functional"].startswith("\n\n[verify] Caveat")


@pytest.mark.parametrize("call", [
    lambda: O.fill_docx_template("nope.docx", {}, output="out.docx"),
    lambda: O.read_document("nope.xlsx"),
    lambda: O.sum_column("nope.xlsx", "Betrag"),
])
def test_a_tool_message_stays_english(call):
    """The other half of the rule: what the MODEL reads is English,
    because the model translates it into the user's language on the way
    out. Only what nothing translates is German."""
    with pytest.raises(O.OfficeError) as excinfo:
        call()
    assert not _GERMAN_MARK.search(str(excinfo.value)), str(excinfo.value)


# ---------------------------------------------------------------------------
# The census: how much German matching this framework carries
# ---------------------------------------------------------------------------

# The modules that are allowed to hold German, and what they hold it for.
_MATCHER_MODULES = ("german", "office", "verify_guard", "api_client",
                    "subagents", "memory_store", "prompt_loader")

_UMLAUT = re.compile(r"[äöüßÄÖÜ]")

# German words that only ever appear in a matcher, never in prose a
# comment would carry. Counting these rather than "any umlaut" keeps the
# census honest when a comment quotes a German sentence.
_GERMAN_MATCHER_WORDS = (
    "summe", "gesamt", "betrag", "beträgt", "betragen", "rechnung",
    "beleg", "buchung", "kostenstelle", "eintrag", "übertrag", "verbuch",
    "erstell", "schreib", "speicher", "prüf", "vermutlich", "ungefähr",
    "geschätzt", "niemals", "tabelle", "spalte", "vorlage", "hintergrund",
    "abhängigkeit", "erlauben", "notizbuch", "lösungsmittel", "zusammen",
)


def _matcher_hits(name: str) -> int:
    import importlib
    import inspect
    source = inspect.getsource(importlib.import_module(
        f"delfin.agent.{name}")).lower()
    return sum(source.count(word) for word in _GERMAN_MATCHER_WORDS)


# The floor, not the target. Raise it when a change adds German matching;
# never lower it to make a refactor pass — the German half is what makes
# this framework usable by the people it was built for.
#
# Measured: 124 before the German-matcher audit was fixed, 236 after, and
# every one of the seven modules went up. (238 at the peak; two of those
# were dead triggers that could never decide anything — see the
# no-superstring rule in the prompt-module test. Deleting dead weight is
# not the same as dropping coverage, which is why the number is written
# down with its reason rather than just tracked.)
_CENSUS_FLOOR = 236


def test_the_german_matchers_are_not_quietly_dropped():
    total = sum(_matcher_hits(name) for name in _MATCHER_MODULES)
    assert total >= _CENSUS_FLOOR, (
        f"German matching fell to {total} occurrences, below the "
        f"{_CENSUS_FLOOR} this framework shipped with")


@pytest.mark.parametrize("name", _MATCHER_MODULES)
def test_every_matcher_module_still_carries_german(name):
    import importlib
    import inspect
    source = inspect.getsource(importlib.import_module(f"delfin.agent.{name}"))
    assert _UMLAUT.search(source), name


def test_the_shared_module_is_where_the_vocabulary_lives():
    """Three defects came from the same vocabulary existing twice. These
    are the ones that must stay single."""
    assert O._NOT_ASSERTED_RE is G.HEDGE_RE
    assert VG._HEDGE_MARKERS is G.HEDGE_RE
    from delfin.agent.deliverables import _slugify as report_slug
    from delfin.agent.memory_store import _slugify as memory_slug
    assert memory_slug("Müller GmbH") == report_slug("Müller GmbH")


# ---------------------------------------------------------------------------
# The premise this whole file rests on: the model answers in the user's
# language, so an English string aimed at the model is translated on the
# way out. Nothing ever checked that the prompt pack actually says so.
# ---------------------------------------------------------------------------
#
# Field report 20260902-080635: the user wrote the task in English —
# "Extract the hyperpolarizability tensor from all folders … and calculate
# the Hyper-Rayleigh scattering hyperpolarizability βHRS from it" — and the
# answer came back in German, headed "## Zusammenfassung".
#
# The rule existed. It was the last clause of a bullet about CODE comments
# ("You still talk to the user in their language"), which is the weakest
# place to put a rule about answers, and two skills contradicted it
# outright by prescribing "German chat summary" and "German prose" as
# their output format regardless of who was asking.

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"


def test_the_pack_states_the_answer_language_as_its_own_rule():
    text = (_PACK / "shared" / "work_cycle_rules.md").read_text()
    assert "Answer in the language the user wrote in" in text
    # Decided per message — not from the project, not from habit.
    assert "latest\n  message" in text or "latest message" in text


def test_no_skill_prescribes_the_language_of_its_own_output():
    """A skill that pins German overrides the user's choice silently,
    and the user never sees the skill file that did it."""
    offenders = []
    for path in sorted((_PACK / "skills").glob("*.md")):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"(?i)\b(?:german|deutsch)\b", line) and \
                    re.search(r"(?i)output|format|summary|prose|antwort", line):
                offenders.append(f"{path.name}:{i}: {line.strip()}")
    assert not offenders, "skills pinning an output language:\n" + \
        "\n".join(offenders)


def test_code_is_still_english_regardless_of_the_conversation():
    """The load-bearing half of the same bullet, which the rewrite must
    not have loosened: what goes INTO code stays English."""
    text = (_PACK / "shared" / "work_cycle_rules.md").read_text()
    assert "Everything you write INTO code is English" in text
    assert "regardless of the conversation" in text
