"""Series work and record-oriented edits — the two everyday office jobs.

Both replace something the model would otherwise improvise, and both
improvisations fail the same way: quietly.

A series built by calling the single filler in a loop ends in "done"
whether or not six rows failed. An edit addressed by cell coordinate
lands on the wrong row as soon as the sheet was paged, and the result
reads perfectly.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import office
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

openpyxl = pytest.importorskip("openpyxl")
pytest.importorskip("pypdf")


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "OFFICE"
    d.mkdir()
    return d


def _perms(ws) -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "office-series"
    return perms


@pytest.fixture
def form(ws):
    reportlab = pytest.importorskip("reportlab")
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    path = ws / "antrag.pdf"
    c = canvas.Canvas(str(path), pagesize=A4)
    c.drawString(60, 780, "Antrag")
    for index, name in enumerate(("Name", "Kostenstelle", "Betrag")):
        c.acroForm.textfield(name=name, x=180, y=740 - 30 * index,
                             width=280, height=18)
    c.save()
    return path


@pytest.fixture
def people(ws):
    """A German export: cp1252, semicolons, decimal commas."""
    path = ws / "personen.csv"
    path.write_bytes(
        ("Name;Kostenstelle;Betrag\n"
         "Müller;4711;1.234,50\n"
         "Özdemir;4712;89,90\n").encode("cp1252"))
    return path


@pytest.fixture
def book(ws):
    wb = openpyxl.Workbook()
    sheet = wb.active
    for row in (["Beleg", "Name", "Kostenstelle", "Betrag"],
                ["R-001", "Müller", "4711", "1.234,50"],
                ["R-002", "Özdemir", "4712", "89,90"],
                ["R-003", "Schmidt", "4711", "450,00"]):
        sheet.append(row)
    path = ws / "vorgaenge.xlsx"
    wb.save(path)
    wb.close()
    return path


# ---------------------------------------------------------------------------
# Series
# ---------------------------------------------------------------------------

def test_one_document_per_row(ws, form, people):
    result = office.fill_series(
        people, form, output_dir=ws / "out",
        name_pattern="Antrag_{Name}.pdf")
    assert result["counts"] == {"ok": 2, "incomplete": 0, "failed": 0}
    names = sorted(p.name for p in (ws / "out").iterdir())
    assert names == ["Antrag_Müller.pdf", "Antrag_Özdemir.pdf"]


def test_the_values_really_land_in_the_fields(ws, form, people):
    office.fill_series(people, form, output_dir=ws / "out",
                       name_pattern="Antrag_{Name}.pdf")
    values = {f["name"]: f["value"] for f in office.pdf_form_fields(
        ws / "out" / "Antrag_Müller.pdf")["fields"]}
    assert values["Name"] == "Müller"
    assert values["Kostenstelle"] == "4711"
    assert values["Betrag"] == "1.234,50"


def test_fields_are_mapped_to_columns_of_the_same_name(ws, form, people):
    result = office.fill_series(people, form, output_dir=ws / "out")
    assert result["mapping"] == {
        "Name": "Name", "Kostenstelle": "Kostenstelle", "Betrag": "Betrag"}


def test_an_explicit_mapping_is_honoured(ws, form):
    table = ws / "kurz.csv"
    table.write_text("Nachname,KST,Summe\nMeier,4711,10\n", encoding="utf-8")
    result = office.fill_series(
        table, form, output_dir=ws / "out",
        mapping={"Name": "Nachname", "Kostenstelle": "KST", "Betrag": "Summe"})
    assert result["counts"]["ok"] == 1


def test_a_mapping_that_names_a_missing_field_is_refused(ws, form, people):
    with pytest.raises(office.OfficeError) as exc:
        office.fill_series(people, form, output_dir=ws / "out",
                           mapping={"Nmae": "Name"})
    assert "Nmae" in str(exc.value)
    assert not (ws / "out").exists() or not list((ws / "out").iterdir())


def test_a_mapping_that_names_a_missing_column_is_refused(ws, form, people):
    with pytest.raises(office.OfficeError) as exc:
        office.fill_series(people, form, output_dir=ws / "out",
                           mapping={"Name": "Vorname"})
    assert "Vorname" in str(exc.value)


def test_nothing_is_written_when_validation_fails(ws, form, people):
    """A batch that fails halfway leaves a directory nobody can tell apart
    from a complete one."""
    out = ws / "out"
    with pytest.raises(office.OfficeError):
        office.fill_series(people, form, output_dir=out,
                           mapping={"Name": "gibtsnicht"})
    assert not out.exists() or list(out.iterdir()) == []


def test_a_row_missing_a_value_is_reported_as_incomplete(ws, form):
    table = ws / "luecke.csv"
    table.write_text("Name,Kostenstelle,Betrag\nMeier,,10\n", encoding="utf-8")
    result = office.fill_series(table, form, output_dir=ws / "out")
    assert result["counts"]["incomplete"] == 1
    entry = result["results"][0]
    assert entry["status"] == "incomplete"
    assert "Kostenstelle" in entry["detail"]
    assert any("not ready to hand over" in n for n in result["notes"])


def test_two_rows_wanting_one_file_name_do_not_overwrite(ws, form):
    """Two different records writing to one file is data loss that looks
    like success."""
    table = ws / "doppelt.csv"
    table.write_text(
        "Name,Kostenstelle,Betrag\nMeier,4711,10\nMeier,4712,20\n",
        encoding="utf-8")
    result = office.fill_series(table, form, output_dir=ws / "out",
                                name_pattern="Antrag_{Name}.pdf")
    assert result["counts"]["ok"] == 1
    failed = [r for r in result["results"] if r["status"] == "failed"]
    assert "already used by row 2" in failed[0]["detail"]
    assert len(list((ws / "out").iterdir())) == 1


def test_an_existing_file_is_never_overwritten(ws, form, people):
    out = ws / "out"
    out.mkdir()
    (out / "Antrag_Müller.pdf").write_bytes(b"vorhanden")
    result = office.fill_series(people, form, output_dir=out,
                                name_pattern="Antrag_{Name}.pdf")
    assert result["counts"]["failed"] == 1
    assert (out / "Antrag_Müller.pdf").read_bytes() == b"vorhanden"


def test_a_truncated_run_says_so(ws, form, people):
    """No silent caps: a limit that dropped rows has to be visible."""
    result = office.fill_series(people, form, output_dir=ws / "out", limit=1)
    assert result["processed"] == 1
    assert result["rows"] == 2
    assert any("only the first 1" in n for n in result["notes"])


def test_a_word_template_series(ws):
    docx = pytest.importorskip("docx")
    doc = docx.Document()
    para = doc.add_paragraph()
    for piece in ("Sehr geehrte(r) {{", "Na", "me}},"):
        para.add_run(piece)
    doc.add_paragraph("Betrag: {{Betrag}}")
    template = ws / "brief.docx"
    doc.save(template)

    table = ws / "leute.csv"
    table.write_text("Name,Betrag\nMüller,10\nÖzdemir,20\n", encoding="utf-8")
    result = office.fill_series(table, template, output_dir=ws / "briefe",
                                name_pattern="Brief_{Name}.docx")
    assert result["counts"]["ok"] == 2
    text = docx.Document(str(ws / "briefe" / "Brief_Müller.docx")).paragraphs
    assert "Sehr geehrte(r) Müller," in text[0].text


def test_a_template_without_fields_is_refused(ws, people):
    reportlab = pytest.importorskip("reportlab")
    from reportlab.pdfgen import canvas
    plain = ws / "plain.pdf"
    c = canvas.Canvas(str(plain))
    c.drawString(60, 780, "Kein Formular")
    c.save()
    with pytest.raises(office.OfficeError) as exc:
        office.fill_series(people, plain, output_dir=ws / "out")
    assert "no form fields" in str(exc.value)


def test_file_names_cannot_escape_the_output_directory(ws, form):
    table = ws / "boese.csv"
    table.write_text(
        "Name,Kostenstelle,Betrag\n../../etc/passwd,4711,10\n",
        encoding="utf-8")
    result = office.fill_series(table, form, output_dir=ws / "out",
                                name_pattern="{Name}.pdf")
    produced = list((ws / "out").iterdir())
    assert len(produced) == 1
    assert produced[0].parent == ws / "out"
    assert ".." not in produced[0].name


def test_the_series_tool_reports_every_group(ws, form, people):
    out = _DocToolExecutor()._dispatch("fill_series", {
        "table": str(people), "template": str(form),
        "output_dir": str(ws / "out"), "name_pattern": "A_{Name}.pdf",
    }, _perms(ws))
    assert "2 document(s) complete" in out
    assert "mapping:" in out


def test_the_series_tool_is_refused_in_plan_mode(ws, form, people):
    perms = _perms(ws)
    perms.mode = "plan"
    out = json.loads(_DocToolExecutor()._dispatch("fill_series", {
        "table": str(people), "template": str(form),
        "output_dir": str(ws / "out"),
    }, perms))
    assert "plan mode" in out["error"]


# ---------------------------------------------------------------------------
# Record-oriented edits
# ---------------------------------------------------------------------------

def test_a_row_is_changed_by_its_key(book):
    result = office.edit_sheet(book, key_column="Beleg", updates=[
        {"key": "R-002", "set": {"Kostenstelle": "4799", "Betrag": "98,90"}},
    ])
    assert result["verified"] is True
    cells = {c["cell"]: (c["old"], c["new"]) for c in result["applied"]}
    assert cells["C3"] == ("4712", "4799")
    assert cells["D3"] == ("89,90", "98,90")
    assert "Özdemir" in office.read_sheet(book)["grid"]


def test_the_change_names_the_record_not_just_the_cell(book):
    applied = office.edit_sheet(book, key_column="Beleg", updates=[
        {"key": "R-003", "set": {"Name": "Schmidt-Weber"}}])["applied"][0]
    assert applied["key"] == "R-003"
    assert applied["column"] == "Name"


def test_an_unknown_key_refuses_the_whole_call(book):
    before = book.read_bytes()
    with pytest.raises(office.OfficeError) as exc:
        office.edit_sheet(book, key_column="Beleg", updates=[
            {"key": "R-001", "set": {"Name": "Neu"}},
            {"key": "R-999", "set": {"Name": "Auch neu"}},
        ])
    assert "R-999" in str(exc.value)
    assert book.read_bytes() == before      # not even the valid one applied


def test_a_duplicate_key_is_refused_rather_than_guessed(ws):
    wb = openpyxl.Workbook()
    for row in (["Beleg", "Betrag"], ["R-001", "10"], ["R-001", "20"]):
        wb.active.append(row)
    path = ws / "doppelt.xlsx"
    wb.save(path)
    wb.close()
    with pytest.raises(office.OfficeError) as exc:
        office.edit_sheet(path, key_column="Beleg",
                          updates=[{"key": "R-001", "set": {"Betrag": "30"}}])
    assert "more than once" in str(exc.value)


def test_an_unknown_column_lists_the_real_ones(book):
    with pytest.raises(office.OfficeError) as exc:
        office.edit_sheet(book, key_column="Beleg", updates=[
            {"key": "R-001", "set": {"Gibtsnicht": "x"}}])
    assert "Gibtsnicht" in str(exc.value)
    assert "Kostenstelle" in str(exc.value)


def test_an_unknown_key_column_is_refused(book):
    with pytest.raises(office.OfficeError) as exc:
        office.edit_sheet(book, key_column="Belegnummer",
                          updates=[{"key": "R-001", "set": {"Name": "x"}}])
    assert "Belegnummer" in str(exc.value)


def test_updates_need_a_key_column(book):
    with pytest.raises(office.OfficeError) as exc:
        office.edit_sheet(book, updates=[{"key": "R-001", "set": {"Name": "x"}}])
    assert "key_column" in str(exc.value)


def test_a_record_is_appended_by_column_name(book):
    """A positional list lands in the wrong columns as soon as the sheet's
    order differs from the caller's assumption."""
    result = office.edit_sheet(book, append_records=[
        {"Betrag": "77,00", "Beleg": "R-004", "Name": "Neu"}])
    assert result["appended_rows"] == 1
    grid = office.read_sheet(book)["grid"]
    row = [line for line in grid.splitlines() if "R-004" in line][0]
    cells = [c.strip() for c in row.split("|")]
    assert cells[1] == "R-004" and cells[2] == "Neu" and cells[4] == "77,00"


def test_an_appended_record_with_an_unknown_column_is_refused(book):
    before = book.read_bytes()
    with pytest.raises(office.OfficeError):
        office.edit_sheet(book, append_records=[{"Gibtsnicht": "x"}])
    assert book.read_bytes() == before


def test_an_edit_is_verified_against_the_saved_file(book):
    """The form filler reads its result back; a sheet edit reported from
    memory would be claiming an outcome from the intent."""
    result = office.edit_sheet(book, edits=[{"cell": "B2", "value": "Geändert"}])
    assert result["verified"] is True


def test_the_tool_reports_the_record_and_the_verification(ws, book):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    ex.execute("read_document", {"path": str(book)}, perms)
    out = json.loads(ex._dispatch("edit_sheet", {
        "path": str(book), "key_column": "Beleg",
        "updates": [{"key": "R-002", "set": {"Betrag": "98,90"}}],
    }, perms))
    assert out["status"] == "ok"
    assert out["verified"] is True
    assert "R-002/Betrag" in out["summary"]


def test_the_tool_still_needs_a_prior_read_for_record_updates(ws, book):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    out = json.loads(ex._dispatch("edit_sheet", {
        "path": str(book), "key_column": "Beleg",
        "updates": [{"key": "R-001", "set": {"Name": "x"}}],
    }, perms))
    assert "without reading it first" in out["error"]
