"""The evidence module's own edge cases, beyond the characterization.

The characterization file pins the shapes the old check judged and the
ones it judged right; this one pins the boundaries the NEW heuristics
draw -- where a word or path is at the edge of being a claim.
"""

from __future__ import annotations

from delfin.agent.task_evidence import (
    check_completion_claim,
    _format_word,
    _is_mentioned,
    _unmet_format,
)


def _changes(*paths):
    return [{"path": p, "ts": 10.0, "created": True} for p in paths]


# ---------------------------------------------------------------------------
# The format word
# ---------------------------------------------------------------------------

def test_bare_english_word_is_not_a_format_promise():
    assert _format_word("Explain the word Tabelle to the user") == ""


def test_a_compound_word_dokument_is_a_format_promise():
    assert _format_word("Erstelle das Word-Dokument") == "word-"


def test_pdf_compound_and_suffix_form_both_promise():
    assert _format_word("PDF-Bericht erstellen") == "pdf"
    assert _format_word("Erstelle den Bericht als PDF") == "pdf"


def test_excel_compound_promises_a_spreadsheet():
    assert _unmet_format(
        "Lege eine Excel-Tabelle an", ["/w/t.md"]) == "excel"
    assert _unmet_format(
        "Lege eine Excel-Tabelle an", ["/w/tabelle.xlsx"]) == ""


def test_a_written_docx_satisfies_a_word_task():
    res = check_completion_claim(
        "Erstelle das Word-Dokument", changes=_changes("/w/brief.docx"))
    assert res["verdict"] == "verified"
    assert res["kind"] == "artifact"


# ---------------------------------------------------------------------------
# The mention guard
# ---------------------------------------------------------------------------

def test_a_mention_construction_guards_only_its_own_clause():
    """A path that is BOTH the object and mentioned: object wins."""
    assert not _is_mentioned(
        "foo.py", "Fix the bug in foo.py, as described in foo.py")


def test_a_pure_mention_is_recognised_in_german_and_english():
    assert _is_mentioned("foo.py", "Bau es wie in foo.py beschrieben")
    assert _is_mentioned("foo.py", "Build it as described in foo.py")


def test_an_unrelated_path_is_not_a_mention():
    assert not _is_mentioned("bar.py", "Fix the bug in bar.py")


def test_a_mentioned_path_does_not_accuse_but_a_written_one_verifies():
    """The write of the WORK file finishes the task; the reference is
    context either way."""
    res = check_completion_claim(
        "Bau die Prüflogik wie in foo.py beschrieben",
        changes=_changes("/w/prueflogik.py"),
        observed=["/w/foo.py"])
    assert res["verdict"] == "verified"
    # Without the write: the write-verb branch still demands a mutation,
    # and a session that wrote nothing at all is honestly unmet.
    res = check_completion_claim(
        "Bau die Prüflogik wie in foo.py beschrieben",
        changes=[], observed=["/w/foo.py"])
    assert res["verdict"] == "unmet"
    assert res["kind"] == "no_change"


def test_a_read_task_over_a_mentioned_file_still_verifies_by_reading_it():
    """For a read task the mention guard must not bite: the file IS the
    object being looked at, however the sentence carries it."""
    res = check_completion_claim(
        "Lies den Stand, wie in foo.py beschrieben",
        changes=[], observed=["/w/foo.py"])
    assert res["verdict"] == "verified"


# ---------------------------------------------------------------------------
# Description paths are context, not claims
# ---------------------------------------------------------------------------

def test_a_path_that_only_the_description_names_never_accuses():
    """Subject: fix bar.py; description mentions where the bug came
    from. The description cannot make the task unmet."""
    res = check_completion_claim(
        "Behebe den Fehler in bar.py",
        description="Der Fehler kam ursprünglich aus foo.py.",
        changes=_changes("/w/bar.py"),
        observed=["/w/foo.py"])
    assert res["verdict"] == "verified"
    assert res["kind"] == "path_write"


def test_a_description_path_still_helps_when_the_subject_names_none():
    """A write task with no path in the subject is judged by the
    journal; the description mentioning a file changes nothing about
    that."""
    res = check_completion_claim(
        "Bereinige die Prüflogik",
        description="Analog zu check_base.py, aber ohne foo.py.",
        changes=_changes("/w/validate.py"))
    assert res["verdict"] == "verified"


# ---------------------------------------------------------------------------
# Bare artefact nouns describe content, not containers
# ---------------------------------------------------------------------------

def test_a_bare_tabelle_word_is_decided_by_the_journal():
    """The headline s4 case: the session appended a markdown section.
    The journal shows the mutation; no spreadsheet is demanded."""
    res = check_completion_claim(
        "Ergänze die Übersicht um eine Tabelle der Testfälle",
        changes=_changes("/w/README.md"), observed=[])
    assert res["verdict"] in ("verified", "unchecked"), res
    assert res["kind"] != "artifact", res


def test_a_bare_tabelle_word_with_no_change_is_unmet_via_the_journal():
    """Not demanding a spreadsheet does not mean demanding nothing: a
    write task that wrote nothing at all is still unmet."""
    res = check_completion_claim(
        "Ergänze die Übersicht um eine Tabelle der Testfälle",
        changes=[], observed=[])
    assert res["verdict"] == "unmet"
    assert res["kind"] == "no_change"
