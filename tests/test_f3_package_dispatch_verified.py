"""Package f3: tools verified through the product's own dispatch.

Each tool is called via execute(name, arguments, permissions) on the real
executor, in a temporary HOME and workspace. Valid calls re-open the produced
file (python-docx, pypdf, openpyxl) to prove the content is there; each tool
also gets one invalid call whose refusal names the cause.

Nothing touches the network, mail transport or any device.
"""

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

docx = pytest.importorskip("docx")
openpyxl = pytest.importorskip("openpyxl")
pypdf = pytest.importorskip("pypdf")
reportlab = pytest.importorskip("reportlab")


@pytest.fixture
def home(monkeypatch, tmp_path):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", lambda: h)
    return h


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "werkraum"
    d.mkdir()
    return d


@pytest.fixture
def perms(ws):
    p = KitToolPermissions(workspace=str(ws))
    p.mode = "acceptEdits"
    p.task_session_id = "f3-package"
    return p


@pytest.fixture
def ex(home):
    return _DocToolExecutor()


def call(ex, perms, name, args):
    return json.loads(ex.execute(name, args, perms))


@pytest.fixture
def book(ws, ex, perms):
    """A workbook with three data rows and one formula cell."""
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Kosten"
    for row in (["Monat", "Betrag"], ["Jan", 100], ["Feb", 140], ["Mrz", 90]):
        sheet.append(row)
    sheet["D1"] = "=SUM(B2:B4)"
    path = ws / "kosten.xlsx"
    wb.save(path)
    wb.close()
    # edit_sheet demands a prior read in this session: the fixture hands
    # back a file the session has looked at.
    ex.execute("read_document", {"path": str(path)}, perms)
    return path


@pytest.fixture
def form(ws):
    """A PDF with two text fields."""
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    path = ws / "antrag.pdf"
    c = canvas.Canvas(str(path), pagesize=A4)
    c.drawString(60, 780, "Antrag")
    acro = c.acroForm
    acro.textfield(name="name", x=180, y=735, width=300, height=18)
    acro.textfield(name="betrag", x=180, y=705, width=300, height=18)
    c.save()
    return path


@pytest.fixture
def two_pdfs(ws):
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    paths = []
    for name, text in (("a.pdf", "erster"), ("b.pdf", "zweiter")):
        p = ws / name
        c = canvas.Canvas(str(p), pagesize=A4)
        c.drawString(60, 700, text)
        c.save()
        paths.append(p)
    return paths


# create_docx ---------------------------------------------------------------

def test_create_docx_content_reopened(ex, perms, ws):
    """The case's name promised a reopen, and it made none.

    Measured before this rewrite: the body asserted ``status == "ok"`` and
    opened nothing, so an executor answering ok for ``create_docx`` while
    writing no file at all passed it (scratch/integrate/cannot_fail.py,
    which is the input this assertion exists to fail on). Both blocks are
    now looked for in the reopened document, and the heading's style is
    checked, because a paragraph the writer silently dropped is exactly
    the defect this package found.
    """
    out = call(ex, perms, "create_docx", {
        "path": str(ws / "bericht.docx"),
        "blocks": [
            {"heading": "Versuch 1", "level": 1},
            {"paragraph": "Der Katalysator wurde erhitzt."},
        ],
    })
    assert out["status"] == "ok", out
    written = docx.Document(str(ws / "bericht.docx"))
    texts = [p.text for p in written.paragraphs]
    assert "Versuch 1" in texts, texts
    assert any("Katalysator" in t for t in texts), texts
    styles = [p.style.name for p in written.paragraphs if p.text]
    assert styles[0].startswith("Heading"), styles


def test_create_docx_refuses_empty_blocks(ex, perms, ws):
    out = call(ex, perms, "create_docx",
               {"path": str(ws / "leer.docx"), "blocks": []})
    assert "blocks" in out["error"]


# fill_docx_template --------------------------------------------------------

def _make_template(ws, name, paragraph):
    doc = docx.Document()
    doc.add_paragraph(paragraph)
    doc.save(str(ws / name))
    return ws / name


def test_fill_docx_template_value_holds(ex, perms, ws):
    src = _make_template(ws, "vorlage.docx", "Lieferant: {{name}}")
    out = call(ex, perms, "fill_docx_template", {
        "path": str(src),
        "output": str(ws / "ausgefuellt.docx"),
        "values": {"name": "Chemikalie"}})
    assert out.get("status") == "ok", out
    texts = [p.text for p in docx.Document(str(ws / "ausgefuellt.docx")).paragraphs]
    assert any("Chemikalie" in t for t in texts)
    assert not any("{{name}}" in t for t in texts)


def test_fill_docx_template_refuses_unknown_placeholder(ex, perms, ws):
    src = _make_template(ws, "vorlage2.docx", "Wert {{unbekannt}}")
    out = call(ex, perms, "fill_docx_template", {
        "path": str(src), "output": str(ws / "out.docx"),
        "values": {"etwasanderes": "x"}})
    assert "error" in out
    assert "unbekannt" in out["error"]


# create_pdf ----------------------------------------------------------------

def test_create_pdf_text_reopened(ex, perms, ws):
    """Both blocks, not just the heading.

    The heading alone was checked here, so the PDF half of the two-kind
    defect (a paragraph dropped, ``status: ok``) passed: measured at
    origin/main, create_pdf writes ``{'headings': 1, 'paragraphs': 0}`` for a
    two-kind block, and "Messreihe" is in the text either way.
    """
    out = call(ex, perms, "create_pdf", {
        "path": str(ws / "bericht.pdf"),
        "blocks": [
            {"heading": "Messreihe"},
            {"paragraph": "Zwoelf Grad."},
        ],
    })
    assert out.get("status") == "ok", out
    text = pypdf.PdfReader(str(ws / "bericht.pdf")).pages[0].extract_text()
    assert "Messreihe" in text, text
    assert "Zwoelf" in text, text


def test_create_pdf_refuses_empty_blocks(ex, perms, ws):
    out = call(ex, perms, "create_pdf",
               {"path": str(ws / "leer.pdf"), "blocks": []})
    assert "blocks" in out["error"]


# fill_pdf_form -------------------------------------------------------

def test_fill_pdf_form_field_holds(ex, perms, ws, form):
    out = call(ex, perms, "fill_pdf_form", {
        "path": str(form),
        "output": str(ws / "ausgefuellt.pdf"),
        "values": {"name": "Salzsaeure", "betrag": "1700"}})
    assert out.get("status") == "ok", out
    fields = pypdf.PdfReader(str(ws / "ausgefuellt.pdf")).get_fields()
    assert fields["name"].get("/V") == "Salzsaeure"
    assert fields["betrag"].get("/V") == "1700"


def test_fill_pdf_form_refuses_empty_values(ex, perms, ws, form):
    out = call(ex, perms, "fill_pdf_form", {
        "path": str(form), "output": str(ws / "x.pdf"), "values": {}})
    assert "values" in out["error"]


# merge_pdfs / split_pdf ----------------------------------------------

def test_merge_pdfs_result_reopened(ex, perms, ws, two_pdfs):
    out = call(ex, perms, "merge_pdfs", {
        "inputs": [str(two_pdfs[0]), str(two_pdfs[1])],
        "output": str(ws / "gesamt.pdf")})
    assert out.get("status") == "ok", out
    reader = pypdf.PdfReader(str(ws / "gesamt.pdf"))
    assert len(reader.pages) == 2
    both = (reader.pages[0].extract_text() + reader.pages[1].extract_text())
    assert "erster" in both and "zweiter" in both


def test_merge_pdfs_refuses_one_input(ex, perms, ws, two_pdfs):
    out = call(ex, perms, "merge_pdfs", {
        "inputs": [str(two_pdfs[0])], "output": str(ws / "x.pdf")})
    assert "at least two" in out["error"]


def test_split_pdf_output_reopened(ex, perms, ws):
    src = ws / "drei.pdf"
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas
    c = canvas.Canvas(str(src), pagesize=A4)
    for i in range(3):
        c.drawString(60, 700, f"Seite {i}")
        c.showPage()
    c.save()
    out = call(ex, perms, "split_pdf", {
        "path": str(src), "output_dir": str(ws / "teile"), "pages": "1-2"})
    assert out.get("status") == "ok", out
    written = out["written"][0]
    split = pypdf.PdfReader(ws / "teile" / written)
    # The page COUNT alone passes on a two-page file made of the wrong
    # pages, or of two blank ones. The range asked for is 1-2, so both of
    # those pages' own text has to be inside the result and the third's
    # must not.
    assert len(split.pages) == 2, len(split.pages)
    both = "".join(page.extract_text() for page in split.pages)
    assert "Seite 0" in both, both
    assert "Seite 1" in both, both
    assert "Seite 2" not in both, both


def test_split_pdf_refuses_unknown_page_syntax(ex, perms, ws):
    src = ws / "zwei.pdf"
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas
    c = canvas.Canvas(str(src), pagesize=A4)
    c.drawString(60, 700, "x")
    c.save()
    out = call(ex, perms, "split_pdf", {
        "path": str(src), "output_dir": str(ws / "teile2"),
                                "pages": "nichts"})
    assert "error" in out


# edit_sheet: the formula survives an edit to another cell -------------

def test_edit_sheet_formula_survives(ex, perms, ws, book):
    out = call(ex, perms, "edit_sheet", {
        "path": str(book), "edits": [{"cell": "B2", "value": 7}]})
    assert out.get("status") == "ok", out
    wb = openpyxl.load_workbook(str(book))
    sheet = wb["Kosten"]
    assert sheet["B2"].value == 7
    # the formula in D1 is still a formula, not a spilled number
    assert str(sheet["D1"].value).startswith("=")


def test_edit_sheet_refuses_without_prior_read(ex, perms, ws):
    wb = openpyxl.Workbook()
    wb.active["A1"] = "x"
    path = ws / "niegelesen.xlsx"
    wb.save(path)
    wb.close()
    out = call(ex, perms, "edit_sheet", {
        "path": str(path), "edits": [{"cell": "A1", "value": 2}]})
    assert "read" in out["error"].lower()


# sum_column / compare_tables -----------------------------------------

def test_sum_column_total_matches(ex, perms, ws, book):
    raw = ex.execute("sum_column",
                                    {"path": str(book), "column": "Betrag"}, perms)
    assert "total of 'Betrag': 330" in raw, raw
    assert "counted 3 of 3 data row(s)" in raw, raw


def test_sum_column_refuses_unknown_column(ex, perms, ws, book):
    out = call(ex, perms, "sum_column", {"path": str(book), "column": "NichtDa"})
    assert "error" in out


def test_compare_tables_reconciles(ex, perms, ws, book):
    second = ws / "kosten2.xlsx"
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Tabelle1"
    sheet.append(["Monat", "Betrag"])
    sheet.append(["Jan", 100])
    wb.save(second)
    wb.close()
    raw = ex.execute("compare_tables",
                     {"left": str(book), "right": str(second),
                      "key": "Monat"}, perms)
    assert "matched on 'Monat'" in raw, raw
    assert "equal:" in raw and "differing:" in raw, raw


def test_compare_tables_refuses_unknown_key(ex, perms, ws, book):
    out = call(ex, perms, "compare_tables",
                {"left": str(book), "right": str(book), "key": "Unfug"})
    assert "error" in out


# fill_series -----------------------------------------------------------

def test_fill_series_one_document_per_row(ex, perms, ws):
    table = ws / "adressen.xlsx"
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.append(["Name"])
    sheet.append(["Ada"])
    sheet.append(["Grace"])
    wb.save(table)
    wb.close()

    template = ws / "brief.docx"
    doc = docx.Document()
    doc.add_paragraph("Hallo {{Name}}")
    doc.save(str(template))

    raw = ex.execute("fill_series", {
        "table": str(table), "template": str(template),
        "output_dir": str(ws / "briefe")}, perms)
    assert "2 document(s) complete" in raw, raw
    assert "mapping: Name <- Name" in raw, raw
    # Output names follow the template stem + row number: brief_2.docx,
    # brief_3.docx (documented naming, not the row value).
    written = sorted((ws / "briefe").glob("*.docx"))
    assert len(written) == 2, raw
    letters = [docx.Document(str(out)).paragraphs[0].text for out in written]
    # Both rows, each in its own document. The weaker form of this check
    # ("Ada is in one of the two") passes when the writer drops a row and
    # duplicates the other, so each name is required on its own.
    assert any("Hallo Ada" in t for t in letters), letters
    assert any("Hallo Grace" in t for t in letters), letters
    for text in letters:
        assert "{{" not in text, text


def test_fill_series_refuses_empty_table(ex, perms, ws):
    empty = ws / "leer.xlsx"
    wb = openpyxl.Workbook()
    wb.save(empty)
    wb.close()
    template = ws / "t.docx"
    docx.Document().save(str(template))
    out = call(ex, perms, "fill_series", {
        "table": str(empty), "template": str(template),
        "output_dir": str(ws / "out")})
    assert "error" in out


# make_plot / publish_report / draft_email ------------------------------

def test_make_plot_writes_figure(ex, perms, ws):
    """The figure is a PNG, not merely a file with the right name.

    A name check passes on an empty file or on a traceback written into
    the target path, which is how a failing plot call used to look.
    """
    raw = ex.execute("make_plot", {
        "spec": {"kind": "line", "x": [1, 2, 3], "y": [2, 4, 9],
                 "filename": "messreihe.png"}}, perms)
    figure = ws / "messreihe.png"
    assert figure.exists(), raw
    assert "messreihe" in raw, raw
    data = figure.read_bytes()
    assert data.startswith(b"\x89PNG\r\n\x1a\n"), data[:16]
    assert len(data) > 1000, f"{len(data)} bytes is not a drawn figure"


def test_make_plot_refuses_unknown_kind(ex, perms, ws):
    out = call(ex, perms, "make_plot",
                                {"spec": {"kind": "kuchendiagramm"}})
    assert "error" in out
    assert "kind" in out["error"]


# publish_report -------------------------------------------------

def test_publish_report_writes_md_and_html(ex, perms, ws):
    """Both formats, because the default is ``both`` and both are claimed.

    The name of this case promised two files and the assertion checked one:
    ``_execute_publish_report`` defaults to ``fmt="both"`` (api_client.py)
    and ``deliverables.publish_report`` writes ``.md`` and ``.html`` on that
    default. A case that asserts the ``.md`` alone passes just as well when
    the HTML write is broken, so it measured nothing about half of its own
    name. Both files are now opened and each has to carry the text.
    """
    out = call(ex, perms, "publish_report", {
        "title": "F3 Paket", "markdown": "# F3\n\nGemessen."})
    assert out.get("status") == "written", out
    md = sorted((ws / "reports").glob("*.md"))
    html = sorted((ws / "reports").glob("*.html"))
    assert md, "no md report written"
    assert html, "no html report written"
    assert "Gemessen" in md[0].read_text(encoding="utf-8")
    page = html[0].read_text(encoding="utf-8")
    assert "Gemessen" in page, page[:400]
    assert "F3 Paket" in page, "the title did not reach the html page"


def test_publish_report_honours_an_explicit_md_only_format(ex, perms, ws):
    """The control for the case above: ``format`` still chooses.

    Without this, "both files exist" would also pass if the tool ignored
    ``fmt`` entirely and always wrote both. Asking for ``md`` alone must
    produce no ``.html`` -- that is the input which makes the assertion in
    the case above able to fail.
    """
    out = call(ex, perms, "publish_report", {
        "title": "Nur Markdown", "markdown": "Text.", "format": "md"})
    assert out.get("status") == "written", out
    assert sorted((ws / "reports").glob("*.md")), "no md report written"
    assert not sorted((ws / "reports").glob("*.html")), \
        "format=md wrote an html file anyway"


def test_publish_report_refuses_empty_title(ex, perms, ws):
    out = call(ex, perms, "publish_report",
               {"title": "", "markdown": "x"})
    assert "title is required" in out["error"]


# draft_email -------------------------------------------------

def test_draft_email_writes_eml_no_transport(ex, perms, ws):
    out = call(ex, perms, "draft_email", {
        "path": str(ws / "nachricht.eml"),
                                  "to": "person@example.org",
        "subject": "Messwerte",
        "body": "Die Messung ist fertig."})
    assert out.get("status") == "ok", out
    raw = (ws / "nachricht.eml").read_text(encoding="utf-8")
    assert "Messwerte" in raw
    assert "Die Messung ist fertig." in raw
    assert not out.get("sent")


def test_draft_email_refuses_missing_subject(ex, perms, ws):
    out = call(ex, perms, "draft_email", {
                                "path": str(ws / "n.eml"), "to": "person@example.org",
                                "body": "ohne Betreff"})
    assert "error" in out
