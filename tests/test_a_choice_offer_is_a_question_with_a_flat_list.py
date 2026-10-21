"""``detect_question`` offers numbered choices only for a choice question.

Input: the text of a finished answer. Output of interest: whether
``detect_question`` returns ``{"type": "numbered", ...}``, which draws one
button per option under the answer ("Pick one or select multiple").

Semantics pinned here. A numbered offer is

* the LAST list of the answer, numbered 1..n with 2 <= n <= 8;
* flat: no indented sub-item, no lettered or bulleted child under an item;
* final: after the list stands nothing, or one short paragraph ending
  in "?";
* asked for: the line directly above it is a choice question ("Welches
  Format möchtest du?") or a request to pick ("Möchtest du:", "tell me
  which to use:"), or the closing paragraph is a question naming a
  choice ("Welche Variante soll ich umsetzen?", "Which one?").

Rejected alternative (the previous rule): any two lines anywhere in the
last 2000 characters that start with a number, plus one "?" anywhere in
the same tail. It turned plans, reports, rankings and questionnaires into
buttons, and a nested plan contributed its sub-items as options too.

Corpus: ``tests/fixtures/choice_offer_corpus.json``, 70 benchmark
answers (``~/.delfin/benchmark_runs/*.jsonl``, field ``text_excerpt``;
3349 distinct texts in all), labelled by reading each one: ``offer``
(must offer), ``none`` (must not), ``either`` (a real choice the rule
declines because more text follows, or an excerpt cut mid-list; not
asserted). Rows: the 56 texts the previous rule offered as numbered
choices (46 none, 4 offer, 6 either); 11 whose line above the list only
introduces it ("Sichere Alternativen:", "Typische Beispiele:"), which a
looser draft of this rule offered; 3 "**Option 1:** ... Welchen Weg
bevorzugst du?" answers the previous rule drew as an open question.

Counts over all 3349 texts, previous rule -> this rule: numbered 56 -> 8
(7 offers + 1 either), open 59 -> 71, yes/no 38 -> 45, none 3196 -> 3225.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.dashboard.tab_agent import detect_question as _detect


_CORPUS = json.loads(
    (Path(__file__).parent / "fixtures" / "choice_offer_corpus.json")
    .read_text(encoding="utf-8"))


def _numbered(text):
    got = _detect(text)
    return got if got and got["type"] == "numbered" else None


# -- the corpus ------------------------------------------------------------

@pytest.mark.parametrize(
    "row", [r for r in _CORPUS if r["label"] == "none"],
    ids=lambda r: r["id"])
def test_no_plan_report_or_questionnaire_is_offered_as_choices(row):
    assert _numbered(row["text"]) is None, row["id"]


@pytest.mark.parametrize(
    "row", [r for r in _CORPUS if r["label"] == "offer"],
    ids=lambda r: r["id"])
def test_a_real_offer_in_the_corpus_is_still_offered(row):
    got = _numbered(row["text"])
    assert got, row["id"]
    nums = [n for n, _ in got["options"]]
    assert nums == [str(i) for i in range(1, len(nums) + 1)], nums
    assert 2 <= len(nums) <= 8


def test_the_corpus_holds_what_the_docstring_says():
    labels = [r["label"] for r in _CORPUS]
    assert len(labels) == 70
    assert labels.count("offer") == 7
    assert labels.count("none") == 57


# -- the shapes, one at a time --------------------------------------------

_NESTED_PLAN = """Mein Vorschlag:

1. **Erst lesen**: die beiden Module ansehen.
2. **Auswertung definieren**:
   1. Häufigkeit jedes Tags
   2. Bookmarks pro Tag
3. **Ausgabe**: Tabelle im Chat oder CSV?
"""


def test_an_answer_ending_in_a_nested_numbered_list_is_not_an_offer():
    assert _numbered(_NESTED_PLAN) is None


def test_bulleted_children_also_make_a_list_nested():
    text = ("Welche Variante willst du?\n\n"
            "1. Tags zählen\n   - pro Bookmark\n2. Tags gruppieren\n")
    assert _numbered(text) is None


def test_lettered_children_also_make_a_list_nested():
    text = ("Welche Variante willst du?\n\n"
            "1. Tags zählen\n a) pro Bookmark\n b) gesamt\n2. Tags gruppieren\n")
    assert _numbered(text) is None


def test_numbered_steps_of_a_report_are_not_an_offer():
    text = ("Alle Schritte wurden ausgeführt:\n\n"
            "1. Functional auf BP86 gesetzt\n2. Basis auf def2-TZVP gesetzt\n"
            "3. Zum Submit-Tab gewechselt\n\n"
            "Möchtest du die Berechnung jetzt einreichen?")
    assert _numbered(text) is None


def test_a_list_followed_by_more_text_is_not_an_offer():
    text = ("Welche Option soll ich nehmen?\n\n"
            "1. CSV\n2. JSON\n\n"
            "Die Datei liegt danach im Workspace. Ich lege außerdem einen "
            "Test an.\n\nDanach melde ich mich.")
    assert _numbered(text) is None


def test_a_list_with_no_question_is_not_an_offer():
    text = "Reihenfolge nach Gap:\n\n1. run_c — 4.401 eV\n2. run_a — 4.073 eV\n"
    assert _numbered(text) is None


def test_a_request_to_run_a_command_sequence_is_not_a_choice():
    text = ("**Frage:** Darf ich folgende Befehle ausführen?\n\n"
            "1. `ss -tlnp | grep 8899` — PID finden\n"
            "2. `kill <PID>` — Prozess beenden\n")
    assert _numbered(text) is None


def test_numbering_that_does_not_start_at_one_is_not_an_offer():
    text = "Welche nimmst du?\n\n2. foo\n3. bar\n"
    assert _numbered(text) is None


def test_a_list_inside_a_code_fence_is_not_an_offer():
    text = "Welche Option?\n\n```\n1. foo\n2. bar\n```\n"
    assert _numbered(text) is None


@pytest.mark.parametrize("text, want", [
    ("Welches Format möchtest du?\n\n1. CSV\n2. JSON\n3. XLSX\n",
     ["CSV", "JSON", "XLSX"]),
    ("Which one should I use?\n1) the cage\n2) the gate\n",
     ["the cage", "the gate"]),
    ("Options — tell me which to use:\n\n1. **Current time**\n"
     "2. **File mtime**\n\nWhich one?", ["Current time", "File mtime"]),
    ("Mögliche Varianten:\n\n1. Nur anzeigen\n2. Zählen pro Tag\n\n"
     "Welche Variante soll ich umsetzen?", ["Nur anzeigen", "Zählen pro Tag"]),
    ("Welche Option?\n\n- **Option 1**: schnell\n- **Option 2**: gründlich\n",
     ["schnell", "gründlich"]),
    ("Welche Option?\n\n(1) schnell\n(2) gründlich\n", ["schnell", "gründlich"]),
    ("Möchtest du:\n\n1. die ganze Suite laufen lassen?\n"
     "2. einen Test benennen?\n", ["die ganze Suite laufen lassen?",
                                    "einen Test benennen?"]),
])
def test_a_flat_list_under_a_choice_question_is_offered(text, want):
    got = _numbered(text)
    assert got, text
    assert [label for _, label in got["options"]] == want


def test_a_wrapped_item_continues_without_making_the_list_nested():
    text = ("Welche Variante willst du?\n\n"
            "1. Tags zählen und die Häufigkeit\n   pro Tag ausgeben\n"
            "2. Tags gruppieren\n")
    got = _numbered(text)
    assert got and [n for n, _ in got["options"]] == ["1", "2"]


def test_a_rejected_list_still_reaches_the_yes_no_and_open_checks():
    """Declining the buttons must not swallow the question under them."""
    got = _detect(_NESTED_PLAN)
    assert got and got["type"] == "open", got
