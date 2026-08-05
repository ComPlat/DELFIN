"""Two claims a finished turn makes about the world, now checked.

Both from field case 20260803-143354, an office run costing $2.83.

A TASK MARKED DONE. The agent completed a task titled "PDF-Bericht" while
create_pdf had failed for a missing library and only a .docx existed. The
task said done, the summary said done, and the report's own file table
lists the PDF as missing. Nothing in the framework compared the two.

A COUNT FROM CUT-SHORT OUTPUT. fill_series reported 29 complete of 31;
the listing that would have settled it was truncated; the answer asserted
"31 PDF-Dateien verifiziert" and then named 29. The truncation marker was
present in the result and simply not honoured — and this was the answer
produced BY the framework's own correction turn, so a second retry was
never going to help.

Both take the caveat consequence rather than a forced correction. A retry
cannot count what is still not in context, and it cannot produce a file
whose library is missing.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _unmet_artifact
from delfin.agent.engine import AgentEngine
from delfin.agent import verify_guard as vg


# ---------------------------------------------------------------------------
# A completed task must have produced what it names
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("subject,produced,missing", [
    ("PDF-Bericht für Juni erstellen", ["out/bericht.docx"], "pdf"),
    ("Create the PDF report", [], "pdf"),
    ("Excel-Tabelle abgleichen", ["notiz.txt"], "excel"),
    ("Brief an die Verwaltung", [], "brief"),
])
def test_a_promised_artifact_that_was_never_produced_is_named(
        subject, produced, missing):
    assert _unmet_artifact(subject, produced) == missing


@pytest.mark.parametrize("subject,produced", [
    ("PDF-Bericht für Juni erstellen", ["out/bericht.pdf"]),
    ("Excel-Tabelle abgleichen", ["abgleich.xlsx"]),
    ("Brief an die Verwaltung", ["brief.docx"]),
    ("Word-Vorlage füllen", ["gefuellt.docx"]),
])
def test_a_task_whose_artifact_exists_says_nothing(subject, produced):
    assert _unmet_artifact(subject, produced) == ""


@pytest.mark.parametrize("subject", [
    "Daten prüfen", "Rückfrage klären", "Die Zahlen abgleichen", "",
])
def test_a_subject_promising_nothing_checkable_is_left_alone(subject):
    """Advisory and narrow on purpose: only an unambiguous extension is
    checked, so the note cannot become noise the model learns to skip."""
    assert _unmet_artifact(subject, ["irgendwas.txt"]) == ""


def test_a_broken_input_does_not_raise():
    assert _unmet_artifact(None, None) == ""
    assert _unmet_artifact("PDF", object()) in ("", "pdf")


# ---------------------------------------------------------------------------
# A count needs output that was not cut short
# ---------------------------------------------------------------------------

def _engine_with_truncation(*tools):
    eng = AgentEngine.__new__(AgentEngine)
    eng._truncated_tools_turn = list(tools)
    return eng


@pytest.mark.parametrize("answer", [
    "Ich habe 31 PDF-Dateien verifiziert.",
    "Die Liste enthält 29 Einträge.",
    "All 31 files were created.",
    "Es wurden 250 Zeilen geprüft.",
])
def test_a_count_after_a_truncated_result_is_flagged(answer):
    assert _engine_with_truncation("list_files")._scan_truncated_counts(answer)


@pytest.mark.parametrize("answer", [
    "Alles erledigt.",
    "Es sind 7 Dateien.",          # single digit: below the floor
    "Der Bericht liegt bereit.",
])
def test_an_answer_without_a_real_count_is_not_flagged(answer):
    assert not _engine_with_truncation("list_files")._scan_truncated_counts(
        answer)


def test_nothing_truncated_means_no_caveat():
    """The whole point is the SOURCE, not the number: a count from output
    the model saw in full is an ordinary answer."""
    assert not _engine_with_truncation()._scan_truncated_counts(
        "Ich habe 31 PDF-Dateien verifiziert.")


def test_the_caveat_names_the_number_and_the_tool():
    text = vg.truncated_output_caveat(["31 PDF-Dateien"], ["list_files"])
    assert "31 PDF-Dateien" in text
    assert "list_files" in text
    assert "gekürzt" in text


def test_the_caveat_says_what_the_number_is_worth():
    """"Estimated, not counted" is the sentence the user needs; "warning"
    on its own tells them nothing they can act on."""
    text = vg.truncated_output_caveat(["31 Dateien"], ["list_files"])
    assert "geschätzt" in text


def test_the_ledger_is_per_turn_and_bounded():
    eng = _engine_with_truncation()
    for i in range(40):
        eng._note_truncated_tool(f"tool_{i}")
    assert len(eng._truncated_tools_turn) <= 20


def test_the_same_tool_is_recorded_once():
    eng = _engine_with_truncation()
    for _ in range(5):
        eng._note_truncated_tool("list_files")
    assert eng._truncated_tools_turn == ["list_files"]
