"""The rest of the document-writer audit: five places that state a fact
they did not check.

Measured against the shipped module, 2026-08-13:

* ``read_docx`` walks ``doc.paragraphs`` and ``doc.tables`` and nothing
  else. A .docx whose Rechnungsnummer sits in the Kopfzeile — which is
  where letterhead puts it — comes back with ``header found: False``,
  ``paragraphs: 1`` stated as a fact and ``notes: []``. Text boxes and
  content controls are invisible the same way, and the paragraph count
  now feeds the figure ledger, so a wrong count is not merely cosmetic.
* ``create_pdf`` verifies by searching the written file for the longest
  word of ≥4 letters. ``_PROBE_RE`` excludes digits, so a Kostenauf-
  stellung of dates and amounts yields NO probe, ``found`` stays True by
  default and ``verified`` reports only that the file opens as a PDF.
* ``fill_pdf_form`` and ``fill_docx_template`` raise on a write error and
  leave the part-written file behind. ``create_pdf`` beside them unlinks.
* ``draft_email`` collects ``attachments_skipped`` and then says "Draft
  written." with nothing about them — the data is there and the sentence
  does not carry it.
* ``create_docx`` is the only writer with no read-back at all: it returns
  path, blocks and counts, and no ``verified``.

One shape, five times: the report is written from what the code MEANT to
do rather than from what the file turned out to contain.
"""

from __future__ import annotations

import pytest

from delfin.agent import office

docx = pytest.importorskip("docx")
pytest.importorskip("reportlab")
pypdf = pytest.importorskip("pypdf")


# ---------------------------------------------------------------------------
# read_docx: everything the object model does not reach
# ---------------------------------------------------------------------------

def _letterhead(path):
    doc = docx.Document()
    doc.sections[0].header.paragraphs[0].text = "Rechnungsnummer RE-2026-0042"
    doc.sections[0].footer.paragraphs[0].text = "Seite 1 von 1"
    doc.add_paragraph("Sehr geehrte Damen und Herren,")
    doc.save(str(path))
    return path


def test_the_number_in_the_letterhead_is_read(tmp_path):
    """A .docx Rechnung read as not containing its own Rechnungsnummer."""
    out = office.read_docx(_letterhead(tmp_path / "rechnung.docx"))
    assert "RE-2026-0042" in out["text"]


def test_the_footer_is_read_too(tmp_path):
    out = office.read_docx(_letterhead(tmp_path / "rechnung.docx"))
    assert "Seite 1 von 1" in out["text"]


def test_the_body_is_still_read(tmp_path):
    out = office.read_docx(_letterhead(tmp_path / "rechnung.docx"))
    assert "Sehr geehrte Damen und Herren," in out["text"]


def test_the_header_is_labelled_so_it_is_not_taken_for_body_text(tmp_path):
    """A number lifted out of a running header is not a body figure, and
    an answer that cites it should be able to say where it came from."""
    out = office.read_docx(_letterhead(tmp_path / "rechnung.docx"))
    assert "header" in out["text"].lower() or "kopfzeile" in out["text"].lower()


def test_a_document_with_no_header_says_nothing_extra(tmp_path):
    doc = docx.Document()
    doc.add_paragraph("nur der Fließtext")
    doc.save(str(tmp_path / "plain.docx"))
    out = office.read_docx(tmp_path / "plain.docx")
    assert out["text"].strip() == "nur der Fließtext"


def test_a_passage_the_object_model_cannot_reach_is_reported(tmp_path):
    """Text boxes and content controls live outside ``doc.paragraphs``.
    python-docx cannot create them, so the part is written directly —
    which is exactly how Word stores a letterhead address block."""
    import zipfile

    src = tmp_path / "with_box.docx"
    doc = docx.Document()
    doc.add_paragraph("Fließtext")
    doc.save(str(src))

    # Splice a text box into the body part.
    with zipfile.ZipFile(src) as bundle:
        parts = {n: bundle.read(n) for n in bundle.namelist()}
    body = parts["word/document.xml"].decode("utf-8")
    box = (
        '<w:p><w:r><w:pict><v:shape><v:textbox><w:txbxContent>'
        '<w:p><w:r><w:t>Empfänger: Musterfirma GmbH</w:t></w:r></w:p>'
        '</w:txbxContent></v:textbox></v:shape></w:pict></w:r></w:p>'
    )
    body = body.replace("</w:body>", box + "</w:body>")
    parts["word/document.xml"] = body.encode("utf-8")
    with zipfile.ZipFile(src, "w") as bundle:
        for name, data in parts.items():
            bundle.writestr(name, data)

    out = office.read_docx(src)
    assert "Musterfirma GmbH" in out["text"]


def test_the_paragraph_count_is_not_stated_as_a_fact_about_the_file(
        tmp_path):
    """It counts what was extracted. When parts of the file are carried
    separately the record has to say which count it is."""
    out = office.read_docx(_letterhead(tmp_path / "rechnung.docx"))
    assert out["paragraphs"] >= 1
    assert "counted" in " ".join(out["notes"]).lower() or out.get(
        "paragraphs_source")


def test_a_broken_zip_does_not_take_the_read_down(tmp_path):
    """The extra reader is an aid; it must not become a reason not to
    read the document at all."""
    path = _letterhead(tmp_path / "rechnung.docx")
    import delfin.agent.office as mod

    def _boom(*a, **k):
        raise OSError("zip unreadable")

    original = mod._docx_unreached_text
    mod._docx_unreached_text = _boom
    try:
        out = office.read_docx(path)
    finally:
        mod._docx_unreached_text = original
    assert "Sehr geehrte Damen und Herren," in out["text"]


# ---------------------------------------------------------------------------
# create_pdf: a document of nothing but figures
# ---------------------------------------------------------------------------

def _kostenaufstellung(path):
    return office.create_pdf(path, [
        {"table": [["12.03.2026", "1.240,50"],
                   ["19.03.2026", "880,00"]]},
    ])


def test_a_document_of_dates_and_amounts_is_actually_checked(tmp_path):
    out = _kostenaufstellung(tmp_path / "kosten.pdf")
    assert out["verified"] is True
    # ... and it checked something, rather than defaulting to True.
    text = " ".join((p.extract_text() or "")
                    for p in pypdf.PdfReader(str(tmp_path / "kosten.pdf")).pages)
    assert "1.240,50" in text.replace(" ", "") or "1.240,50" in text


def test_a_figure_that_did_not_survive_the_write_is_reported(tmp_path,
                                                             monkeypatch):
    """The point of the probe: if the file does not read back with what
    was put in it, verified must be False and say so."""
    import delfin.agent.office as mod
    real = mod._squash
    monkeypatch.setattr(mod, "_squash",
                        lambda t: "MISMATCH" if "1.240" in str(t) else real(t))
    out = _kostenaufstellung(tmp_path / "kosten.pdf")
    assert out["verified"] is False
    assert any("read back" in n for n in out["notes"])


def test_a_document_with_no_checkable_content_does_not_claim_it_checked(
        tmp_path):
    """Nothing to probe for is not the same as probed and found.

    Both halves asserted, not either: a mutation that dropped the
    ``verified`` half survived an ``or`` here, because the note fired on
    its own and the report still read verified.
    """
    out = office.create_pdf(tmp_path / "empty.pdf", [{"paragraph": ""}])
    assert out["verified"] is False
    assert any("nothing about its CONTENT" in n for n in out["notes"])


def test_an_ordinary_document_still_verifies(tmp_path):
    out = office.create_pdf(tmp_path / "brief.pdf", [
        {"heading": "Zahlungserinnerung", "level": 1},
        {"paragraph": "Wir bitten um Überweisung des offenen Betrags."},
    ])
    assert out["verified"] is True
    assert out["notes"] == [] or all("read back" not in n for n in out["notes"])


# ---------------------------------------------------------------------------
# create_docx: the only writer that never looked
# ---------------------------------------------------------------------------

def test_a_written_document_is_read_back(tmp_path):
    out = office.create_docx(tmp_path / "brief.docx", [
        {"heading": "Zahlungserinnerung", "level": 1},
        {"paragraph": "Wir bitten um Überweisung des offenen Betrags."},
    ])
    assert out["verified"] is True


def test_a_document_that_did_not_survive_is_reported(tmp_path, monkeypatch):
    import delfin.agent.office as mod
    real = mod.read_docx
    monkeypatch.setattr(
        mod, "read_docx",
        lambda p, **k: {"text": "", "paragraphs": 0, "tables": 0,
                        "notes": [], "path": str(p), "kind": "word"})
    out = office.create_docx(tmp_path / "brief.docx",
                             [{"paragraph": "Überweisung"}])
    assert out["verified"] is False
    assert out["notes"]
    assert real is not None


def test_a_table_only_document_verifies_on_its_figures(tmp_path):
    out = office.create_docx(tmp_path / "kosten.docx", [
        {"table": [["12.03.2026", "1.240,50"]]},
    ])
    assert out["verified"] is True


# ---------------------------------------------------------------------------
# The write error paths leave nothing behind
# ---------------------------------------------------------------------------

def _save_that_dies_mid_write(monkeypatch):
    """What a real failure looks like: the file exists, then the write
    dies. A stub that raises BEFORE creating anything leaves nothing
    behind whatever the code does — it cannot tell an unlink from its
    absence, and two of these tests proved nothing until a mutation check
    said so."""
    def _half_written(self, target):
        from pathlib import Path
        Path(str(target)).write_bytes(b"PK\x03\x04partial")
        raise OSError("disk full")

    monkeypatch.setattr(docx.document.Document, "save", _half_written)


def test_a_failed_docx_fill_leaves_no_file(tmp_path, monkeypatch):
    template = tmp_path / "vorlage.docx"
    doc = docx.Document()
    doc.add_paragraph("Sehr geehrte(r) {{name}},")
    doc.save(str(template))

    out = tmp_path / "brief.docx"
    _save_that_dies_mid_write(monkeypatch)
    with pytest.raises(office.OfficeError):
        office.fill_docx_template(template, {"name": "Frau Meier"},
                                  output=out)
    assert not out.exists(), "a part-written document was left behind"


def test_a_failed_create_docx_leaves_no_file(tmp_path, monkeypatch):
    out = tmp_path / "brief.docx"
    _save_that_dies_mid_write(monkeypatch)
    with pytest.raises(office.OfficeError):
        office.create_docx(out, [{"paragraph": "Text"}])
    assert not out.exists()


# ---------------------------------------------------------------------------
# draft_email says what it dropped
# ---------------------------------------------------------------------------

def test_a_skipped_attachment_is_named_in_the_note(tmp_path):
    out = office.draft_email(
        tmp_path / "mail.eml",
        to=["empfaenger@example.org"], subject="Rechnung",
        body="Anbei die Rechnung.",
        attachments=[str(tmp_path / "gibt-es-nicht.pdf")],
    )
    assert out["attachments_skipped"]
    note = out["note"]
    assert "1" in note
    assert "attach" in note.lower()


def test_a_draft_with_every_attachment_says_nothing_extra(tmp_path):
    doc = tmp_path / "rechnung.pdf"
    office.create_pdf(doc, [{"paragraph": "Rechnungsbetrag"}])
    out = office.draft_email(
        tmp_path / "mail.eml",
        to=["empfaenger@example.org"], subject="Rechnung",
        body="Anbei.", attachments=[str(doc)],
    )
    assert out["attachments_skipped"] == []
    assert "skipped" not in out["note"].lower()


def test_the_draft_still_says_it_was_not_sent(tmp_path):
    """The sentence that matters most must survive the addition."""
    out = office.draft_email(
        tmp_path / "mail.eml", to=["a@example.org"], subject="s", body="b",
        attachments=[str(tmp_path / "fehlt.pdf")])
    assert "NOT sent" in out["note"]
