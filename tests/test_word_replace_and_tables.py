"""Find and replace, paragraphs and tables in the Word view.

The mechanics are Word's, not invented ones: replacing a word inside a bold
sentence leaves the sentence bold, Enter puts the new paragraph *after* the one
the cursor is in, and a table cannot lose its last row or column.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import docx_view

docx = pytest.importorskip("docx")


@pytest.fixture
def document(tmp_path):
    doc = docx.Document()
    doc.add_paragraph("The ligand binds the metal.")
    formatted = doc.add_paragraph()
    formatted.add_run("The ").bold = True
    formatted.add_run("ligand").bold = True
    formatted.add_run(" is bulky.").bold = True
    doc.add_paragraph("LIGAND in capitals, and ligands in plural.")
    table = doc.add_table(rows=2, cols=3)
    table.cell(0, 0).text = "ligand"
    path = tmp_path / "report.docx"
    doc.save(path)
    return path


def _texts(path):
    doc = docx.Document(str(path))
    return [p.text for p in doc.paragraphs]


# ---------------------------------------------------------------------------
# find and replace
# ---------------------------------------------------------------------------
def test_replace_all_replaces_every_occurrence(document):
    result = docx_view.replace_all(document, "ligand", "chelate")
    docx_view.save(result["document"], document)

    assert result["replaced"] >= 3
    text = "\n".join(_texts(document))
    assert "chelate" in text
    assert "ligand" not in text.lower().replace("chelate", "")


def test_replace_reaches_inside_tables(document):
    result = docx_view.replace_all(document, "ligand", "chelate")
    docx_view.save(result["document"], document)

    doc = docx.Document(str(document))
    assert doc.tables[0].cell(0, 0).text == "chelate"


def test_matching_case_replaces_only_what_matches(document):
    result = docx_view.replace_all(document, "LIGAND", "CHELATE", match_case=True)
    docx_view.save(result["document"], document)

    text = "\n".join(_texts(document))
    assert "CHELATE in capitals" in text
    assert "The ligand binds" in text, "the lower-case one must be left alone"


def test_whole_word_does_not_replace_inside_a_longer_word(document):
    result = docx_view.replace_all(document, "ligand", "chelate", whole_word=True)
    docx_view.save(result["document"], document)

    text = "\n".join(_texts(document))
    assert "ligands in plural" in text, "the plural is a different word"


def test_replacing_a_word_leaves_the_formatting_around_it(document):
    result = docx_view.replace_all(document, "ligand", "chelate")
    docx_view.save(result["document"], document)

    doc = docx.Document(str(document))
    bold_paragraph = doc.paragraphs[1]
    assert "chelate" in bold_paragraph.text
    assert all(run.bold for run in bold_paragraph.runs if run.text.strip()), (
        "replacing a word flattened the formatting of the sentence"
    )


def test_nothing_is_written_until_the_caller_saves(document):
    before = document.read_bytes()

    docx_view.replace_all(document, "ligand", "chelate")

    assert document.read_bytes() == before


def test_searching_for_nothing_is_refused(document):
    with pytest.raises(docx_view.DocxError):
        docx_view.replace_all(document, "", "x")


def test_the_two_switches_change_which_matches_exist(document):
    read = docx_view.read_document(document)

    loose = docx_view.search_matches(read, "ligand")
    cased = docx_view.search_matches(read, "ligand", match_case=True)
    whole = docx_view.search_matches(read, "ligand", whole_word=True)

    assert len(loose) > len(cased), "case-insensitive must find the capitals too"
    assert len(loose) > len(whole), "whole word must not find the plural"


# ---------------------------------------------------------------------------
# paragraphs
# ---------------------------------------------------------------------------
def test_a_new_paragraph_lands_after_the_one_it_was_asked_for(document):
    doc = docx.Document(str(document))

    docx_view.insert_paragraph(doc, "p:0", "Inserted.")
    docx_view.save(doc, document)

    assert _texts(document)[1] == "Inserted."


def test_a_new_paragraph_can_be_asked_to_land_before(document):
    doc = docx.Document(str(document))

    docx_view.insert_paragraph(doc, "p:0", "First.", before=True)
    docx_view.save(doc, document)

    assert _texts(document)[0] == "First."


def test_a_paragraph_can_be_taken_out(document):
    doc = docx.Document(str(document))
    before = len(doc.paragraphs)

    docx_view.delete_paragraph(doc, "p:0")
    docx_view.save(doc, document)

    assert len(_texts(document)) == before - 1
    assert "binds the metal" not in "\n".join(_texts(document))


def test_an_address_that_is_not_there_is_reported(document):
    doc = docx.Document(str(document))
    with pytest.raises(docx_view.DocxError):
        docx_view.delete_paragraph(doc, "p:999")


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def test_a_table_can_be_added(document):
    doc = docx.Document(str(document))

    docx_view.add_table(doc, 3, 4)
    docx_view.save(doc, document)

    tables = docx.Document(str(document)).tables
    assert len(tables) == 2
    assert len(tables[1].rows) == 3 and len(tables[1].columns) == 4


def test_a_table_gains_and_loses_rows(document):
    doc = docx.Document(str(document))

    docx_view.add_table_row(doc, 0)
    docx_view.save(doc, document)
    assert len(docx.Document(str(document)).tables[0].rows) == 3

    doc = docx.Document(str(document))
    docx_view.drop_table_row(doc, 0, 0)
    docx_view.save(doc, document)
    assert len(docx.Document(str(document)).tables[0].rows) == 2


def test_a_table_gains_and_loses_columns(document):
    doc = docx.Document(str(document))

    docx_view.add_table_column(doc, 0)
    docx_view.save(doc, document)
    assert len(docx.Document(str(document)).tables[0].columns) == 4

    doc = docx.Document(str(document))
    docx_view.drop_table_column(doc, 0, 3)
    docx_view.save(doc, document)
    assert len(docx.Document(str(document)).tables[0].columns) == 3


def test_a_table_keeps_its_last_row_and_column(document):
    doc = docx.Document(str(document))
    docx_view.add_table(doc, 1, 1)

    with pytest.raises(docx_view.DocxError):
        docx_view.drop_table_row(doc, 1, 0)
    with pytest.raises(docx_view.DocxError):
        docx_view.drop_table_column(doc, 1, 0)


def test_a_table_that_would_be_enormous_is_refused(document):
    doc = docx.Document(str(document))
    with pytest.raises(docx_view.DocxError):
        docx_view.add_table(doc, 500, 500)
