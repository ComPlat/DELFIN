"""What ``check_completion_claim`` actually does today, per shape.

Characterization, not specification: every assertion here is the CURRENT
behaviour of ``delfin.agent.api_client.check_completion_claim``, green on
the unchanged code. The supervised run of 2026-09-22 (s4) reported
``unmet``/``unchecked`` notes that keyed on WORDING of the task text
rather than on the work done; the cases marked MISJUDGMENT reproduce
those shapes and pin today's wrong answer, so a rewrite
(``delfin.agent.task_evidence``) has something to be measured against.
The CONTROL cases are the ones that judge correctly today and must keep
doing so.

The distinction the old check misses: a word or a filename can be the
OBJECT of the task (then it is a claim) or it can merely occur in the
text (then it promises nothing).
"""

from __future__ import annotations

from delfin.agent.api_client import check_completion_claim


def _changes(*paths):
    return [{"path": p, "ts": 10.0, "created": True} for p in paths]


# ---------------------------------------------------------------------------
# CONTROL — correct today, must stay correct
# ---------------------------------------------------------------------------

def test_control_a_named_file_that_was_written_is_verified():
    res = check_completion_claim(
        "Ergänze eine Tabelle der Testfälle in REIBUNG.md",
        changes=_changes("/w/REIBUNG.md"))
    assert res["verdict"] == "verified"
    assert res["kind"] == "path_write"


def test_control_a_named_file_that_was_not_written_is_unmet():
    res = check_completion_claim(
        "Add mylib/optimizers/wrapper.py",
        changes=_changes("/proj/other.py"),
        observed=["/proj/mylib/optimizers/wrapper.py"])
    assert res["verdict"] == "unmet"
    assert res["kind"] == "path_unwritten"


def test_control_a_read_task_is_finished_by_the_read():
    res = check_completion_claim(
        "analysiere mylib/core.py", changes=[], observed=["/p/mylib/core.py"])
    assert res["verdict"] == "verified"
    assert res["kind"] == "path_read"


def test_control_an_artifact_word_with_a_format_qualifier_still_demands_it():
    """The original incident: a PDF task with only a docx on disk."""
    res = check_completion_claim(
        "Erstelle den Bericht als PDF",
        changes=_changes("/w/bericht.docx"))
    assert res["verdict"] == "unmet"
    assert res["kind"] == "artifact"


def test_control_a_written_pdf_satisfies_a_pdf_task():
    res = check_completion_claim(
        "PDF-Bericht für Juni erstellen",
        changes=_changes("/w/bericht.pdf"))
    assert res["verdict"] == "verified"


def test_control_an_edit_task_with_no_mutation_is_unmet():
    res = check_completion_claim("Refactor the optimizer wrapper", changes=[])
    assert res["verdict"] == "unmet"
    assert res["kind"] == "no_change"


def test_control_a_read_only_task_is_not_asked_to_produce_an_artifact():
    res = check_completion_claim(
        "Prüfe den Bericht auf Fehler", changes=[], observed=[])
    assert res["verdict"] != "unmet"


def test_control_a_subject_with_nothing_to_key_on_is_unchecked():
    res = check_completion_claim(
        "Rückfrage mit Jerome klären", changes=[], observed=[], tests=[])
    assert res["verdict"] == "unchecked"


# ---------------------------------------------------------------------------
# MISJUDGMENT — the shapes from the supervised run of 2026-09-22 (s4).
# Each pins today's answer; the target module must judge them right.
# The rule to judge by: a word or filename is a CLAIM only when it is the
# object of the task. Mentioned in passing, it promises nothing — and
# when the evidence cannot decide, the answer is "unchecked", never a
# false "unmet" that makes the agent redo finished work.
# ---------------------------------------------------------------------------

# 1. A table in a markdown file: "Tabelle" in a subordinate role does not
#    promise a spreadsheet. Done right: the .md write is the claim.
def test_misjudgment_a_markdown_table_is_not_a_spreadsheet_promise():
    res = check_completion_claim(
        "Ergänze eine Tabelle der Testfälle in REIBUNG.md",
        changes=_changes("/w/REIBUNG.md"),
        observed=[],
    )
    # Wrong today: artifact/Tabelle demands .xlsx/.csv even though the
    # named file was written. (The path branch answers first only when a
    # path is matched; when the write went to a differently-named file the
    # word alone decides — and accuses.)
    assert res["verdict"] in ("verified", "unchecked"), res
    assert res["kind"] not in ("artifact",), res


# 2. The variant that really fired: the .md write is recorded, so the path
#    branch saves it — but only by luck of ordering. The same subject
#    where the table is the object and the .md is only mentioned:
def test_misjudgment_a_table_word_in_a_subordinate_clause():
    """Subject: append a table of test cases to the overview (no file
    named). The session appended a markdown section; the word "Tabelle"
    fired and demanded a .xlsx/.csv."""
    res = check_completion_claim(
        "Ergänze die Übersicht um eine Tabelle der Testfälle",
        changes=_changes("/w/README.md"),
        observed=[],
    )
    assert res["verdict"] != "unmet", res


# 3. A file only MENTIONED is not a file to change.
def test_misjudgment_a_mentioned_file_is_not_a_promised_write():
    """Build the check logic as described in foo.py — the work went into
    a new module, foo.py was only the reference. Today foo.py in the
    text makes the task unmet until foo.py itself is written."""
    res = check_completion_claim(
        "Bau die Prüflogik wie in foo.py beschrieben",
        changes=_changes("/w/prueflogik.py"),
        observed=["/w/foo.py"],
    )
    assert res["verdict"] != "unmet", res


# 4. The same, comparative: "wie in foo.py" vs. a task that genuinely
#    targets foo.py. The second must stay unmet when foo.py is untouched.
def test_misjudgment_a_compared_to_file_is_not_a_promised_write():
    res = check_completion_claim(
        "Passe die Validierung an, analog zu check_base.py",
        changes=_changes("/w/validate_new.py"),
        observed=["/w/check_base.py"],
    )
    assert res["verdict"] != "unmet", res


# 5. A read task whose subject names a file the session has a ledger for
#    must not be judged by the write branch just because the description
#    carries a write verb ("notiere deine Erkenntnisse in Stichpunkten").
def test_misjudgment_a_read_task_with_a_note_instruction_in_the_body():
    res = check_completion_claim(
        "Lies den Stand in engine.py und notiere Stichpunkte",
        description="Ergebnis als Notiz in der Antwort, keine Datei.",
        changes=[],
        observed=["/w/engine.py"],
    )
    assert res["verdict"] != "unmet", res


# 6. A write task about one file that mentions a second in the
#    description must not be judged on the mentioned one.
def test_misjudgment_the_second_named_file_in_the_description_wins():
    """Subject: fix the bug in bar.py; the description says the bug
    originally came from foo.py. Today the first unmatched candidate
    decides — the ORDER of the paths in the text, not the object of the
    task."""
    res = check_completion_claim(
        "Behebe den Fehler in bar.py",
        description="Der Fehler kam ursprünglich aus foo.py.",
        changes=_changes("/w/bar.py"),
        observed=["/w/foo.py"],
    )
    assert res["verdict"] == "verified", res
    assert res["kind"] == "path_write", res

