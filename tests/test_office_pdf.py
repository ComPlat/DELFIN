"""Assembling PDFs and producing one: merge, split, create.

The failure modes pinned here are the ones that look like success:

* a merge that stops at an unreadable input and leaves a file which
  reads as finished while documents are missing from it;
* an output name that already belonged to a document somebody has
  already sent, silently replaced;
* a batch reported as a count, with the one part that was never written
  hidden inside it;
* a page count taken from what was intended rather than from the file
  that now sits on disk;
* German names losing their umlauts on the way into a generated PDF,
  which nobody notices until the letter is printed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import change_journal as cj
from delfin.agent import office
from delfin.agent.api_client import (
    _DOC_TOOLS_OPENAI,
    _OFFICE_TOOL_NAMES,
    KitToolPermissions,
    ToolSurfaceContext,
    _DocToolExecutor,
    tool_unavailable_reason,
)

pypdf = pytest.importorskip("pypdf")
pytest.importorskip("reportlab")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "OFFICE"
    d.mkdir()
    return d


@pytest.fixture
def home(monkeypatch, tmp_path):
    """Keep the undo journal out of the real home directory."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", lambda: h)
    return h


def _perms(ws: Path, mode: str = "acceptEdits") -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = mode
    perms.task_session_id = "office-pdf"
    return perms


def _pdf(path: Path, pages: int, label: str) -> Path:
    """A PDF whose every page says which document and page it is."""
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    c = canvas.Canvas(str(path), pagesize=A4)
    for number in range(1, pages + 1):
        c.drawString(60, 780, f"{label} Seite {number}")
        c.showPage()
    c.save()
    return path


def _page_text(path: Path) -> list[str]:
    return [(page.extract_text() or "").strip()
            for page in pypdf.PdfReader(str(path)).pages]


@pytest.fixture
def anlage_a(ws):
    return _pdf(ws / "anlage_a.pdf", 2, "A")


@pytest.fixture
def anlage_b(ws):
    return _pdf(ws / "anlage_b.pdf", 3, "B")


@pytest.fixture
def verschluesselt(ws, anlage_a):
    """A PDF nobody can read without the password."""
    writer = pypdf.PdfWriter(clone_from=str(anlage_a))
    writer.encrypt("geheim")
    path = ws / "gesperrt.pdf"
    with path.open("wb") as fh:
        writer.write(fh)
    return path


@pytest.fixture
def form(ws):
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    path = ws / "antrag.pdf"
    c = canvas.Canvas(str(path), pagesize=A4)
    c.drawString(60, 780, "Antrag")
    c.acroForm.textfield(name="Name", x=180, y=740, width=280, height=18)
    c.save()
    return path


# ---------------------------------------------------------------------------
# merge_pdfs
# ---------------------------------------------------------------------------

def test_merge_reports_the_pages_of_every_input_and_of_the_result(
        ws, anlage_a, anlage_b):
    result = office.merge_pdfs([anlage_a, anlage_b], output=ws / "gesamt.pdf")
    assert [entry["pages"] for entry in result["inputs"]] == [2, 3]
    assert result["expected_pages"] == 5
    assert result["pages"] == 5
    assert result["verified"] is True


def test_the_reported_page_count_comes_from_the_written_file(
        ws, anlage_a, anlage_b):
    out = ws / "gesamt.pdf"
    result = office.merge_pdfs([anlage_a, anlage_b], output=out)
    assert result["pages"] == len(pypdf.PdfReader(str(out)).pages)


def test_merge_keeps_the_order_it_was_given(ws, anlage_a, anlage_b):
    out = ws / "gesamt.pdf"
    office.merge_pdfs([anlage_b, anlage_a], output=out)
    assert _page_text(out) == ["B Seite 1", "B Seite 2", "B Seite 3",
                               "A Seite 1", "A Seite 2"]


def test_merge_refuses_an_output_that_already_exists(ws, anlage_a, anlage_b):
    out = ws / "gesamt.pdf"
    out.write_bytes(b"%PDF-1.4 fremd")
    with pytest.raises(office.OfficeError, match="already exists"):
        office.merge_pdfs([anlage_a, anlage_b], output=out)
    # The file that was already there is the one somebody may have sent.
    assert out.read_bytes() == b"%PDF-1.4 fremd"


def test_an_encrypted_input_names_the_file_and_writes_nothing(
        ws, anlage_a, verschluesselt):
    out = ws / "gesamt.pdf"
    with pytest.raises(office.OfficeError) as excinfo:
        office.merge_pdfs([anlage_a, verschluesselt], output=out)
    assert "gesperrt.pdf" in str(excinfo.value)
    assert "encrypted" in str(excinfo.value)
    # Not a partial merge: nothing exists under the requested name.
    assert not out.exists()


def test_an_unreadable_input_names_the_file_and_writes_nothing(ws, anlage_a):
    kaputt = ws / "kaputt.pdf"
    kaputt.write_bytes(b"das ist kein PDF")
    out = ws / "gesamt.pdf"
    with pytest.raises(office.OfficeError) as excinfo:
        office.merge_pdfs([anlage_a, kaputt], output=out)
    assert "kaputt.pdf" in str(excinfo.value)
    assert not out.exists()


def test_a_missing_input_is_named(ws, anlage_a):
    with pytest.raises(office.OfficeError, match="fehlt.pdf"):
        office.merge_pdfs([anlage_a, ws / "fehlt.pdf"],
                          output=ws / "gesamt.pdf")


def test_a_non_pdf_input_is_refused(ws, anlage_a):
    notiz = ws / "notiz.txt"
    notiz.write_text("kein PDF", encoding="utf-8")
    with pytest.raises(office.OfficeError, match="not a PDF"):
        office.merge_pdfs([anlage_a, notiz], output=ws / "gesamt.pdf")


def test_merging_needs_at_least_two_documents(ws, anlage_a):
    with pytest.raises(office.OfficeError, match="at least two"):
        office.merge_pdfs([anlage_a], output=ws / "gesamt.pdf")


def test_the_output_suffix_has_to_be_pdf(ws, anlage_a, anlage_b):
    with pytest.raises(office.OfficeError, match="use .pdf"):
        office.merge_pdfs([anlage_a, anlage_b], output=ws / "gesamt.docx")


def test_merging_forms_warns_that_field_values_collapse(ws, form, anlage_a):
    result = office.merge_pdfs([form, anlage_a], output=ws / "gesamt.pdf")
    assert result["verified"] is True
    assert any("antrag.pdf" in note and "collapse" in note
               for note in result["notes"])


# ---------------------------------------------------------------------------
# split_pdf
# ---------------------------------------------------------------------------

def test_without_a_selection_every_page_becomes_its_own_file(ws, anlage_b):
    result = office.split_pdf(anlage_b, output_dir=ws / "teile")
    assert result["counts"] == {"ok": 3, "incomplete": 0, "failed": 0}
    assert sorted(p.name for p in (ws / "teile").iterdir()) == [
        "anlage_b_p1.pdf", "anlage_b_p2.pdf", "anlage_b_p3.pdf"]
    assert _page_text(ws / "teile" / "anlage_b_p2.pdf") == ["B Seite 2"]


def test_a_range_stays_one_document(ws, anlage_b):
    result = office.split_pdf(anlage_b, output_dir=ws / "teile", pages="2-3,1")
    names = [entry["output"] for entry in result["files"]]
    assert names == ["anlage_b_p2-3.pdf", "anlage_b_p1.pdf"]
    assert len(pypdf.PdfReader(str(ws / "teile" / names[0])).pages) == 2
    assert len(pypdf.PdfReader(str(ws / "teile" / names[1])).pages) == 1
    assert result["verified"] is True


def test_the_page_syntax_is_the_one_the_reader_uses(ws, anlage_b):
    # Same parser as read_pdf: a page beyond the end is named as such
    # rather than quietly dropped.
    with pytest.raises(office.OfficeError, match="has 3 page"):
        office.split_pdf(anlage_b, output_dir=ws / "teile", pages="7")
    with pytest.raises(office.OfficeError, match="bad page number"):
        office.split_pdf(anlage_b, output_dir=ws / "teile", pages="zwei")


def test_a_taken_name_fails_that_part_and_leaves_the_file_alone(ws, anlage_b):
    out_dir = ws / "teile"
    out_dir.mkdir()
    (out_dir / "anlage_b_p2.pdf").write_bytes(b"%PDF-1.4 alt")

    result = office.split_pdf(anlage_b, output_dir=out_dir)
    assert result["counts"] == {"ok": 2, "incomplete": 0, "failed": 1}
    failed = [e for e in result["files"] if e["status"] == "failed"]
    assert failed[0]["output"] == "anlage_b_p2.pdf"
    assert "already exists" in failed[0]["detail"]
    assert (out_dir / "anlage_b_p2.pdf").read_bytes() == b"%PDF-1.4 alt"
    # The batch does not report itself as done while a part is missing.
    assert result["verified"] is False


def test_every_part_is_reported_with_its_own_outcome(ws, anlage_b):
    out_dir = ws / "teile"
    out_dir.mkdir()
    (out_dir / "anlage_b_p1.pdf").write_bytes(b"%PDF-1.4 alt")
    result = office.split_pdf(anlage_b, output_dir=out_dir)
    assert len(result["files"]) == 3
    assert [e["status"] for e in result["files"]] == ["failed", "ok", "ok"]


def test_splitting_a_form_warns_that_the_fields_lose_their_document(ws, form):
    result = office.split_pdf(form, output_dir=ws / "teile")
    assert any("form" in note for note in result["notes"])


def test_split_refuses_a_non_pdf(ws):
    notiz = ws / "notiz.txt"
    notiz.write_text("kein PDF", encoding="utf-8")
    with pytest.raises(office.OfficeError, match="not a PDF"):
        office.split_pdf(notiz, output_dir=ws / "teile")


def test_split_refuses_an_encrypted_source(ws, verschluesselt):
    with pytest.raises(office.OfficeError, match="encrypted"):
        office.split_pdf(verschluesselt, output_dir=ws / "teile")


# ---------------------------------------------------------------------------
# create_pdf
# ---------------------------------------------------------------------------

_BRIEF = [
    {"heading": "Rechnungsübersicht", "level": 1},
    {"paragraph": "Sehr geehrte Frau Müller, anbei die Aufstellung."},
    {"table": [["Beleg", "Kostenstelle", "Betrag"],
               ["R-001", "4711", "1.234,50"],
               ["R-002", "4712", "89,90"]],
     "header_row": True},
    {"page_break": True},
    {"paragraph": "Mit freundlichen Grüßen, Özdemir"},
]


def test_umlauts_survive_into_the_pdf_text(ws):
    """German names are the everyday case; a mangled one is found in print."""
    target = ws / "brief.pdf"
    office.create_pdf(target, [
        {"heading": "Prüfbericht für Größe und Straße", "level": 1},
        {"paragraph": "Frau Müller und Herr Özdemir, Bearbeitungsgebühr 5 €."},
        {"table": [["Name", "Ort"], ["Weiß", "Köln"], ["Schäfer", "Zürich"]]},
    ])
    text = office.read_pdf(target)["text"]
    for word in ("Prüfbericht", "Größe", "Straße", "Müller", "Özdemir",
                 "Weiß", "Köln", "Zürich", "€"):
        assert word in text, word


def test_markup_characters_in_the_text_are_written_as_text(ws):
    target = ws / "brief.pdf"
    result = office.create_pdf(target, [
        {"paragraph": "Meier & Söhne <AG> — Anteil > 50%"},
    ])
    assert result["verified"] is True
    assert "Meier & Söhne <AG>" in office.read_pdf(target)["text"]


def test_blocks_become_headings_paragraphs_tables_and_pages(ws):
    result = office.create_pdf(ws / "brief.pdf", _BRIEF)
    assert result["headings"] == 1
    assert result["paragraphs"] == 2
    assert result["tables"] == 1
    assert result["pages"] == 2
    assert result["verified"] is True


def test_the_table_content_is_really_in_the_document(ws):
    target = ws / "brief.pdf"
    office.create_pdf(target, _BRIEF)
    text = office.read_pdf(target)["text"]
    assert "Kostenstelle" in text and "1.234,50" in text


def test_create_pdf_verifies_against_the_written_file(ws):
    target = ws / "brief.pdf"
    result = office.create_pdf(target, _BRIEF)
    assert result["pages"] == len(pypdf.PdfReader(str(target)).pages)


def test_create_pdf_refuses_an_existing_file(ws):
    target = ws / "brief.pdf"
    target.write_bytes(b"%PDF-1.4 alt")
    with pytest.raises(office.OfficeError, match="already exists"):
        office.create_pdf(target, _BRIEF)
    assert target.read_bytes() == b"%PDF-1.4 alt"


def test_create_pdf_refuses_a_non_pdf_suffix(ws):
    with pytest.raises(office.OfficeError, match="use .pdf"):
        office.create_pdf(ws / "brief.docx", _BRIEF)


def test_a_block_that_is_not_content_is_named_by_its_index(ws):
    with pytest.raises(office.OfficeError, match="block 2"):
        office.create_pdf(ws / "brief.pdf",
                          [{"paragraph": "ok"}, {"absatz": "falsch"}])
    assert not (ws / "brief.pdf").exists()


def test_an_empty_block_list_is_refused(ws):
    with pytest.raises(office.OfficeError, match="non-empty"):
        office.create_pdf(ws / "brief.pdf", [])


def test_the_same_blocks_work_for_word_and_for_pdf(ws):
    """One block shape for both writers, so a caller can switch format."""
    word = office.create_docx(ws / "brief.docx", _BRIEF)
    pdf = office.create_pdf(ws / "brief.pdf", _BRIEF)
    assert (word["headings"], word["paragraphs"], word["tables"]) == (
        pdf["headings"], pdf["paragraphs"], pdf["tables"])


# ---------------------------------------------------------------------------
# The tool layer
# ---------------------------------------------------------------------------

def test_the_merge_tool_reports_every_input_and_is_undoable(
        ws, home, anlage_a, anlage_b):
    ex = _DocToolExecutor()
    target = ws / "gesamt.pdf"
    out = json.loads(ex._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(anlage_b)], "output": str(target),
    }, _perms(ws)))
    assert out["status"] == "ok"
    assert [entry["pages"] for entry in out["inputs"]] == [2, 3]
    assert out["pages"] == 5
    assert out["verified"] is True

    assert cj.revert("office-pdf", scope="last")["reverted"] == [str(target)]
    assert not target.exists()


def test_the_merge_tool_refuses_an_existing_output(ws, anlage_a, anlage_b):
    target = ws / "gesamt.pdf"
    target.write_bytes(b"%PDF-1.4 alt")
    out = json.loads(_DocToolExecutor()._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(anlage_b)], "output": str(target),
    }, _perms(ws)))
    assert "already exists" in out["error"]
    assert target.read_bytes() == b"%PDF-1.4 alt"


def test_the_merge_tool_names_the_input_it_could_not_read(
        ws, anlage_a, verschluesselt):
    out = json.loads(_DocToolExecutor()._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(verschluesselt)],
        "output": str(ws / "gesamt.pdf"),
    }, _perms(ws)))
    assert "gesperrt.pdf" in out["error"]
    assert not (ws / "gesamt.pdf").exists()


def test_the_split_tool_names_the_parts_that_failed(ws, home, anlage_b):
    out_dir = ws / "teile"
    out_dir.mkdir()
    (out_dir / "anlage_b_p1.pdf").write_bytes(b"%PDF-1.4 alt")
    out = json.loads(_DocToolExecutor()._dispatch("split_pdf", {
        "path": str(anlage_b), "output_dir": str(out_dir),
    }, _perms(ws)))
    assert out["counts"] == {"ok": 2, "incomplete": 0, "failed": 1}
    assert out["problems"][0]["output"] == "anlage_b_p1.pdf"
    assert out["verified"] is False
    assert out["written"] == ["anlage_b_p2.pdf", "anlage_b_p3.pdf"]


def test_the_split_tool_takes_a_page_selection(ws, home, anlage_b):
    out = json.loads(_DocToolExecutor()._dispatch("split_pdf", {
        "path": str(anlage_b), "output_dir": str(ws / "teile"),
        "pages": "1-2",
    }, _perms(ws)))
    assert out["written"] == ["anlage_b_p1-2.pdf"]
    assert out["verified"] is True


def test_the_create_pdf_tool_writes_and_is_undoable(ws, home):
    ex = _DocToolExecutor()
    target = ws / "brief.pdf"
    out = json.loads(ex._dispatch("create_pdf", {
        "path": str(target), "blocks": _BRIEF,
    }, _perms(ws)))
    assert out["status"] == "ok" and out["verified"] is True
    assert out["pages"] == 2
    assert "Müller" in office.read_pdf(target)["text"]

    assert cj.revert("office-pdf", scope="last")["reverted"] == [str(target)]
    assert not target.exists()


def test_the_create_pdf_tool_refuses_an_existing_file(ws):
    target = ws / "brief.pdf"
    target.write_bytes(b"%PDF-1.4 alt")
    out = json.loads(_DocToolExecutor()._dispatch("create_pdf", {
        "path": str(target), "blocks": _BRIEF,
    }, _perms(ws)))
    assert "already exists" in out["error"]
    assert target.read_bytes() == b"%PDF-1.4 alt"


@pytest.mark.parametrize("name,arguments", [
    ("merge_pdfs", {"inputs": ["a.pdf", "b.pdf"], "output": "c.pdf"}),
    ("split_pdf", {"path": "a.pdf", "output_dir": "teile"}),
    ("create_pdf", {"path": "brief.pdf", "blocks": [{"paragraph": "x"}]}),
])
def test_the_assembly_tools_are_refused_in_plan_mode(ws, name, arguments):
    out = json.loads(_DocToolExecutor()._dispatch(
        name, arguments, _perms(ws, mode="plan")))
    assert "plan mode" in out["error"]


def test_a_diff_approval_decline_leaves_the_document_unwritten(
        ws, anlage_a, anlage_b):
    perms = _perms(ws, mode="diff_approval")
    perms.confirm_callback = lambda name, args, preview: False
    target = ws / "gesamt.pdf"
    out = json.loads(_DocToolExecutor()._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(anlage_b)], "output": str(target),
    }, perms))
    assert "declined" in out["error"]
    assert not target.exists()


def test_a_diff_approval_preview_names_what_would_be_written(
        ws, anlage_a, anlage_b):
    seen: dict = {}

    def confirm(name, args, preview):
        seen["preview"] = preview
        return True

    perms = _perms(ws, mode="diff_approval")
    perms.confirm_callback = confirm
    _DocToolExecutor()._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(anlage_b)],
        "output": str(ws / "gesamt.pdf"),
    }, perms)
    assert "anlage_a.pdf" in seen["preview"]
    assert "gesamt.pdf" in seen["preview"]


def test_the_assembly_tools_need_permissions(ws, anlage_a, anlage_b):
    out = json.loads(_DocToolExecutor()._dispatch("merge_pdfs", {
        "inputs": [str(anlage_a), str(anlage_b)],
        "output": str(ws / "gesamt.pdf"),
    }, None))
    assert "permissions" in out["error"]


# ---------------------------------------------------------------------------
# Surface and audit registration
# ---------------------------------------------------------------------------

def test_the_new_tools_are_in_the_catalogue():
    names = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
    assert {"merge_pdfs", "split_pdf", "create_pdf"} <= names


def test_the_new_tools_are_dropped_without_a_document_backend():
    ctx = ToolSurfaceContext(has_office_libs=False)
    for name in ("merge_pdfs", "split_pdf", "create_pdf"):
        assert name in _OFFICE_TOOL_NAMES
        assert tool_unavailable_reason(name, ctx) is not None


def test_the_new_tools_land_in_the_audit_trail():
    from delfin.agent import audit_log

    for name in ("merge_pdfs", "split_pdf", "create_pdf"):
        assert name in audit_log._WRITE_TOOLS
        assert name in _DocToolExecutor._AUDITED_TOOLS


def test_a_form_without_default_resources_can_still_be_filled(tmp_path):
    """An AcroForm may legally omit /DR, and some producers do. pypdf reads
    it while writing a value and calls get_object() on the plain dictionary
    it substitutes -- so a form that is merely missing an optional entry
    failed with an AttributeError naming neither the form nor the field."""
    fitz = pytest.importorskip("fitz")

    doc = fitz.open()
    page = doc.new_page()
    widget = fitz.Widget()
    widget.field_name = "begruendung"
    widget.field_type = fitz.PDF_WIDGET_TYPE_TEXT
    widget.field_flags = 1 << 12          # several lines
    widget.rect = fitz.Rect(60, 100, 400, 220)
    page.add_widget(widget)
    source = tmp_path / "antrag.pdf"
    doc.save(str(source))
    doc.close()

    # The premise: this document really has no /DR.
    writer = pypdf.PdfWriter(clone_from=str(source))
    assert "/DR" not in writer._root_object["/AcroForm"]

    target = tmp_path / "gefuellt.pdf"
    office.fill_pdf_form(source, {"begruendung": "Zeile eins\nZeile zwei"},
                         output=target)
    written = {f["name"]: f["value"]
               for f in office.pdf_form_fields(target)["fields"]}
    assert written["begruendung"] == "Zeile eins\nZeile zwei"


# ---------------------------------------------------------------------------
# Bold and colour
# ---------------------------------------------------------------------------

def test_a_paragraph_can_carry_emphasis_without_becoming_all_bold(ws):
    """What a report needs: the label plain, the figure red."""
    target = ws / "bericht.pdf"
    office.create_pdf(target, [
        {"paragraph": [
            {"text": "Abweichung bei "},
            {"text": "R-2026-002", "bold": True},
            {"text": " 89,90 gegen 98,90", "bold": True, "color": "red"},
        ]},
    ])
    text = office.read_pdf(target)["text"]
    assert "Abweichung bei R-2026-002" in text.replace("\n", " ")
    # The markup must not survive as visible characters.
    assert "<b>" not in text and "<font" not in text


def test_markup_characters_in_formatted_text_stay_text(ws):
    target = ws / "escape.pdf"
    office.create_pdf(target, [
        {"paragraph": "Hinweis & <Sonderzeichen>", "color": "grey"},
    ])
    assert "<Sonderzeichen>" in office.read_pdf(target)["text"]


def test_an_unknown_colour_is_refused_with_the_ones_that_work(ws):
    with pytest.raises(office.OfficeError) as exc:
        office.create_pdf(ws / "x.pdf", [
            {"paragraph": "text", "color": "himbeere"}])
    assert "himbeere" in str(exc.value)
    assert "#RRGGBB" in str(exc.value)


def test_colours_parse_by_name_and_by_hex():
    assert office.parse_colour("red") == (192, 0, 0)
    assert office.parse_colour("#00A0FF") == (0, 160, 255)
    assert office.parse_colour("00a0ff") == (0, 160, 255)
    assert office.parse_colour(None) is None
    assert office.parse_colour("") is None


def test_plain_text_still_works_unchanged(ws):
    """Every existing caller passes a string; none of them may break."""
    target = ws / "einfach.pdf"
    result = office.create_pdf(target, [
        {"heading": "Titel", "level": 1},
        {"paragraph": "Ein ganz normaler Absatz."},
        {"table": [["A", "B"], ["1", "2"]], "header_row": True},
    ])
    assert result["verified"] is True
    assert "normaler Absatz" in office.read_pdf(target)["text"]
