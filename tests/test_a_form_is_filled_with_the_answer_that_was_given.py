"""Three ways the document writers said yes and meant something else.

1. A Ja/Nein dropdown recorded the opposite answer. pypdf fills
   ``/_States_`` for a CHOICE field out of its ``/Opt``, and the filler
   branched on that list being present rather than on the field's type --
   so a dropdown went through the check-box mapper, which resolves any
   truthy word to the FIRST non-/Off state. On a field declared
   ``["Nein", "Ja"]``, "ja" was written as "Nein", reported verified
   (the read-back compares against what was just written) under a note
   calling it a check box. Einverständnis, Widerspruch, Datenweitergabe:
   every German consent field has that shape.

2. The template filler verified with the traversal that filled it. Both
   the writer and the read-back walk paragraphs and tables, so a
   placeholder in a text box -- which is how Word letterhead builds an
   address window -- was invisible to both, and the check could never
   catch the writer's blind spot. A Serienbrief printed {{name}} on every
   Bescheid and fill_series stamped each row ok.

3. Filling over an existing document destroyed it without a word. Every
   sibling writer refuses that; these two did not, and the PDF handler
   even reads the pre-image for the undo journal, so it knew the file was
   there.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

pypdf = pytest.importorskip("pypdf")
docx = pytest.importorskip("docx")

from delfin.agent import office                       # noqa: E402
from delfin.agent.office import OfficeError           # noqa: E402


# ---------------------------------------------------------------------------
# A PDF with a Ja/Nein choice, built here so the shape is explicit
# ---------------------------------------------------------------------------

def _consent_pdf(path: Path, options=("Nein", "Ja")) -> Path:
    from pypdf.generic import (ArrayObject, DictionaryObject, NameObject,
                               NumberObject, TextStringObject)
    writer = pypdf.PdfWriter()
    writer.add_blank_page(width=595, height=842)
    field = DictionaryObject({
        NameObject("/FT"): NameObject("/Ch"),
        NameObject("/T"): TextStringObject("Einverstaendnis"),
        NameObject("/V"): TextStringObject(str(options[0])),
        NameObject("/Opt"): ArrayObject(
            [TextStringObject(str(o)) for o in options]),
        NameObject("/Ff"): NumberObject(131072),
        NameObject("/Type"): NameObject("/Annot"),
        NameObject("/Subtype"): NameObject("/Widget"),
        NameObject("/Rect"): ArrayObject(
            [NumberObject(v) for v in (50, 700, 250, 720)]),
    })
    ref = writer._add_object(field)
    writer.pages[0][NameObject("/Annots")] = ArrayObject([ref])
    writer._root_object[NameObject("/AcroForm")] = writer._add_object(
        DictionaryObject({NameObject("/Fields"): ArrayObject([ref])}))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        writer.write(fh)
    return path


def _value_of(path: Path, field: str = "Einverstaendnis"):
    return (pypdf.PdfReader(str(path)).get_fields() or {})[field].get("/V")


def test_ja_is_written_as_ja(tmp_path):
    """The defect, at its plainest. This wrote 'Nein'."""
    src = _consent_pdf(tmp_path / "einwilligung.pdf")
    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(src, {"Einverstaendnis": "ja"}, output=str(out))
    assert str(_value_of(out)) == "Ja"


def test_nein_is_written_as_nein(tmp_path):
    src = _consent_pdf(tmp_path / "einwilligung.pdf")
    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(src, {"Einverstaendnis": "NEIN"}, output=str(out))
    assert str(_value_of(out)) == "Nein"


def test_the_order_of_the_options_does_not_decide_the_answer(tmp_path):
    """Same question, options the other way round. The old mapper took
    whichever came first, so this pair is what tells a real answer from a
    coincidence."""
    src = _consent_pdf(tmp_path / "andersrum.pdf", options=("Ja", "Nein"))
    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(src, {"Einverstaendnis": "nein"}, output=str(out))
    assert str(_value_of(out)) == "Nein"


def test_an_answer_the_form_does_not_offer_is_refused_with_the_options(
        tmp_path):
    """It used to be written verbatim with no note at all."""
    src = _consent_pdf(tmp_path / "einwilligung.pdf")
    with pytest.raises(OfficeError) as caught:
        office.fill_pdf_form(src, {"Einverstaendnis": "vielleicht"},
                             output=str(tmp_path / "x.pdf"))
    message = str(caught.value)
    assert "Nein" in message and "Ja" in message


def test_an_option_pair_is_read_as_one_option(tmp_path):
    """``/Opt`` may hold ``[export, display]`` pairs -- a real form says
    ``["J", "Ja"]``. Flattening the pair with str() produced the text
    "['J', 'Ja']", which no comparison and no message could use."""
    from pypdf.generic import (ArrayObject, DictionaryObject, NameObject,
                               NumberObject, TextStringObject)
    path = tmp_path / "paare.pdf"
    writer = pypdf.PdfWriter()
    writer.add_blank_page(width=595, height=842)
    field = DictionaryObject({
        NameObject("/FT"): NameObject("/Ch"),
        NameObject("/T"): TextStringObject("Einverstaendnis"),
        NameObject("/V"): TextStringObject("N"),
        NameObject("/Opt"): ArrayObject([
            ArrayObject([TextStringObject("N"), TextStringObject("Nein")]),
            ArrayObject([TextStringObject("J"), TextStringObject("Ja")]),
        ]),
        NameObject("/Ff"): NumberObject(131072),
        NameObject("/Type"): NameObject("/Annot"),
        NameObject("/Subtype"): NameObject("/Widget"),
        NameObject("/Rect"): ArrayObject(
            [NumberObject(v) for v in (50, 700, 250, 720)]),
    })
    ref = writer._add_object(field)
    writer.pages[0][NameObject("/Annots")] = ArrayObject([ref])
    writer._root_object[NameObject("/AcroForm")] = writer._add_object(
        DictionaryObject({NameObject("/Fields"): ArrayObject([ref])}))
    with path.open("wb") as fh:
        writer.write(fh)

    # The label is what a person answering the form reads, so it is
    # accepted -- and the EXPORT value is what the document stores.
    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(path, {"Einverstaendnis": "Ja"}, output=str(out))
    assert str(_value_of(out)) == "J"

    listed = office.pdf_form_fields(path)["fields"][0]["states"]
    assert listed == ["Nein (N)", "Ja (J)"], listed


def test_a_check_box_still_takes_a_loose_yes(tmp_path):
    """The other half: /Btn keeps the mapper that made it convenient."""
    from pypdf.generic import (ArrayObject, DictionaryObject, NameObject,
                               NumberObject, TextStringObject)
    path = tmp_path / "kasten.pdf"
    writer = pypdf.PdfWriter()
    writer.add_blank_page(width=595, height=842)
    field = DictionaryObject({
        NameObject("/FT"): NameObject("/Btn"),
        NameObject("/T"): TextStringObject("gelesen"),
        NameObject("/V"): NameObject("/Off"),
        NameObject("/Type"): NameObject("/Annot"),
        NameObject("/Subtype"): NameObject("/Widget"),
        NameObject("/AP"): DictionaryObject({
            NameObject("/N"): DictionaryObject({
                NameObject("/Ja"): DictionaryObject({}),
                NameObject("/Off"): DictionaryObject({}),
            })}),
        NameObject("/Rect"): ArrayObject(
            [NumberObject(v) for v in (50, 650, 70, 670)]),
    })
    ref = writer._add_object(field)
    writer.pages[0][NameObject("/Annots")] = ArrayObject([ref])
    writer._root_object[NameObject("/AcroForm")] = writer._add_object(
        DictionaryObject({NameObject("/Fields"): ArrayObject([ref])}))
    with path.open("wb") as fh:
        writer.write(fh)

    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(path, {"gelesen": True}, output=str(out))
    assert str(_value_of(out, "gelesen")) != "/Off"


# ---------------------------------------------------------------------------
# The template filler checks its own blind spot
# ---------------------------------------------------------------------------

_TEXTBOX = (
    '<w:p><w:r><w:pict>'
    '<v:shape xmlns:v="urn:schemas-microsoft-com:vml" '
    'style="width:200pt;height:60pt">'
    '<v:textbox><w:txbxContent><w:p><w:r><w:t>{{name}}</w:t></w:r></w:p>'
    '</w:txbxContent></v:textbox></v:shape></w:pict></w:r></w:p>')


def _letterhead_template(tmp_path: Path) -> Path:
    """A template whose address window is a text box, as Word builds it."""
    plain = tmp_path / "plain.docx"
    document = docx.Document()
    document.add_paragraph("Sehr geehrte(r) {{anrede}},")
    document.save(str(plain))

    tpl = tmp_path / "briefvorlage.docx"
    with zipfile.ZipFile(plain) as src, \
            zipfile.ZipFile(tpl, "w", zipfile.ZIP_DEFLATED) as dst:
        for item in src.infolist():
            data = src.read(item.filename)
            if item.filename == "word/document.xml":
                data = data.decode("utf-8").replace(
                    "</w:body>", _TEXTBOX + "</w:body>").encode("utf-8")
            dst.writestr(item, data)
    return tpl


def test_a_placeholder_the_writer_cannot_reach_is_reported(tmp_path):
    """filled: 1, unfilled: [], complete: True -- and {{name}} printed on
    every Bescheid."""
    tpl = _letterhead_template(tmp_path)
    out = tmp_path / "bescheid.docx"
    result = office.fill_docx_template(
        tpl, {"anrede": "Frau Müller"}, output=str(out))
    assert result["complete"] is False
    assert "name" in result["unfilled"]


def test_a_template_it_can_reach_still_completes(tmp_path):
    """The half that keeps this from being a way to fail every fill."""
    tpl = tmp_path / "einfach.docx"
    document = docx.Document()
    document.add_paragraph("Sehr geehrte(r) {{anrede}},")
    document.save(str(tpl))
    out = tmp_path / "fertig.docx"
    result = office.fill_docx_template(
        tpl, {"anrede": "Frau Müller"}, output=str(out))
    assert result["complete"] is True and result["unfilled"] == []


# ---------------------------------------------------------------------------
# Nothing is overwritten
# ---------------------------------------------------------------------------

def test_filling_a_pdf_over_an_existing_document_is_refused(tmp_path):
    src = _consent_pdf(tmp_path / "einwilligung.pdf")
    out = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(src, {"Einverstaendnis": "Ja"}, output=str(out))
    with pytest.raises(OfficeError, match="already exists"):
        office.fill_pdf_form(src, {"Einverstaendnis": "Nein"}, output=str(out))
    # The first answer is still the one in the file.
    assert str(_value_of(out)) == "Ja"


def test_filling_a_template_over_an_existing_document_is_refused(tmp_path):
    tpl = tmp_path / "einfach.docx"
    document = docx.Document()
    document.add_paragraph("{{anrede}}")
    document.save(str(tpl))
    out = tmp_path / "brief.docx"
    office.fill_docx_template(tpl, {"anrede": "Frau Müller"}, output=str(out))
    with pytest.raises(OfficeError, match="already exists"):
        office.fill_docx_template(tpl, {"anrede": "Herr Meier"},
                                  output=str(out))
    assert "Frau Müller" in office.read_docx(out)["text"]
