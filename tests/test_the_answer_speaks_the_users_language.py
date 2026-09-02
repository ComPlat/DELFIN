"""The prompt rule was measured to be insufficient, so a mechanism asks.

Field report 20260902-080635: the task was written in English and the
answer came back headed "## Zusammenfassung". Three things were wrong in
sequence, each found only by running it:

  1. The rule lived in work_cycle_rules.md, which reaches only the four
     full-context roles — so it was never in the solo prompt at all.
  2. Against a real HOME a stored memory said "User prefers German
     language for communication" (2026-07-28) and outranked it.
  3. With the rule delivered, no memory in the sandbox, and five of 1109
     prompt lines carrying any German — all of them examples of what a
     USER might type — Qwen still answered in German four times in four
     runs.

(1) and (2) are fixed in the pack. (3) is why this exists: a rule the
model does not follow is not a mechanism, and the reader can see the
language even when the model cannot.

FUNCTION WORDS, not vocabulary. They are the part of a text that does not
move with the subject, so an answer about `beta_HRS_au` and `xTB` carries
the same scaffolding as one about anything else, and they survive a code
block. A noun list would be beaten by every technical term.

Silent whenever it cannot tell. That is most of the time — a path, a
number, a short acknowledgement, a code block — and it has to be, because
a wrong verdict here forces a correction turn on an answer that was fine.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg


# The two texts this was built from, quoted from the reports.
FIELD_QUESTION = (
    "Extract the hyperpolarizability tensor from all folders in "
    "/pfs/data6/home/ka/ka_ibcs/ka_ew7404/archive and calculate the "
    "Hyper-Rayleigh scattering hyperpolarizability betaHRS from it.")
FIELD_ANSWER = (
    "## Zusammenfassung\nIch habe den Hyperpolarisierbarkeitstensor aus "
    "allen 404 Ordnern im Archive extrahiert und die Hyper-Rayleigh-"
    "Streuung-Hyperpolarisierbarkeit berechnet. Die Werte sind in der "
    "Tabelle unten aufgeführt und wurden nach der Standardformel ermittelt.")
PROBE_ANSWER = (
    "Die Funktion `add` nimmt zwei Parameter `a` und `b` entgegen, führt "
    "jedoch eine Subtraktion aus und gibt damit nicht die Summe zurück, "
    "sondern die Differenz der beiden Werte.")


# ---------------------------------------------------------------------------
# The reported case
# ---------------------------------------------------------------------------

def test_the_reported_turn_is_caught():
    assert vg.scan_for_language_mismatch(FIELD_ANSWER, FIELD_QUESTION) == "en"


def test_the_live_probe_turn_is_caught():
    question = ("Read calc.py and write one short paragraph to answer.txt "
                "explaining what the function add actually does.")
    assert vg.scan_for_language_mismatch(PROBE_ANSWER, question) == "en"


def test_the_feedback_asks_only_for_the_language():
    msg = vg.language_mismatch_feedback("en")
    assert "English" in msg
    assert "the content is not in question" in msg
    # The precedence that the field case turned on.
    assert "remembered preference does not override" in msg
    # And the half that must not be swept along.
    assert "code stays English" in msg


def test_the_caveat_reaches_the_reader_in_their_own_language():
    assert "Frage war auf Deutsch" in vg.language_mismatch_caveat("de")
    assert "question was in English" in vg.language_mismatch_caveat("en")
    assert vg.language_mismatch_caveat("") == ""


# ---------------------------------------------------------------------------
# What must stay silent — most of it
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("question, answer, why", [
    ("Kannst du mir bitte sagen, was die Funktion add in dieser Datei "
     "eigentlich macht?", PROBE_ANSWER, "German asked, German given"),
    (FIELD_QUESTION,
     "I extracted the tensor from all of the folders and computed the value "
     "with the formula that is stored in the file, which is the one that "
     "should be used here.", "English asked, English given"),
    (FIELD_QUESTION, "171232.01 au", "the answer is a number"),
    (FIELD_QUESTION, "/pfs/data6/home/ka/calc/TADFs-Jeneesh_3/ESD/S1.out",
     "the answer is a path"),
    (FIELD_QUESTION, "```python\nfor x in items:\n    print(x)\n```",
     "the answer is code"),
    ("ok", PROBE_ANSWER, "the question is too short to have a language"),
    ("", PROBE_ANSWER, "there is no question"),
    (FIELD_QUESTION, "", "there is no answer"),
])
def test_it_says_nothing_when_it_cannot_say(question, answer, why):
    assert vg.scan_for_language_mismatch(answer, question) == "", why


def test_a_technical_answer_is_not_mistaken_for_english():
    """German prose quoting English identifiers stays German."""
    answer = ("Die Werte für `beta_HRS_au` und `beta_zzz_au` sind in der "
              "Datei DELFIN_Data.json gespeichert und werden von der "
              "Funktion parse_orca_output gelesen, nicht neu berechnet.")
    assert vg.detect_language(answer) == "de"


def test_english_prose_quoting_german_names_stays_english():
    answer = ("The values are stored in the Ordner that the user named, and "
              "the file DELFIN_Data.json is the one that has been read for "
              "this, so the number does not have to be computed again.")
    assert vg.detect_language(answer) == "en"


def test_a_short_text_has_no_language():
    for text in ("OK.", "Fertig.", "Done.", "42", "beta_HRS_au = 171232.01"):
        assert vg.detect_language(text) == "", text


# ---------------------------------------------------------------------------
# The guard's own input must be the USER's message
# ---------------------------------------------------------------------------

def test_the_guards_own_feedback_is_not_read_as_the_user():
    """`[Verify] …` is fed back in AS a user message by the correction
    turns, and those are always English. If the guard read one back it
    would demand English of every German session — the very defect it
    exists to catch, built into the catcher.

    The engine records the user's message at the top of the turn in
    `_last_user_message`, and a nested correction overwrites it. So the
    check ignores a body that starts with the marker, and this pins that
    the marker survives to be recognised.
    """
    from delfin.agent import verify_guard as vg
    guard_feedback = ("[Verify] The following cited paths do not exist in "
                      "the workspace: 'a.out'. Read or grep the actual "
                      "files and correct the answer.")
    assert guard_feedback.lstrip().startswith("[Verify]")
    # Read as a question it would look English and force a switch.
    assert vg.detect_language(guard_feedback) == "en"
    assert vg.scan_for_language_mismatch(PROBE_ANSWER, guard_feedback) == "en"
    # Ignored, nothing is claimed.
    assert vg.scan_for_language_mismatch(PROBE_ANSWER, "") == ""
