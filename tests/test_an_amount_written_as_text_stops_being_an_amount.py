"""A correct edit turned a money column into a mixed one, silently.

`edit_sheet` put whatever it was handed straight into the cell. A person
correcting an amount writes it the way the rest of their file is written
— "1.265,85" — and a model passing that through did the reasonable thing.
The cell then held TEXT where every other row held a number, the tool
answered `"status": "ok", "verified": true`, the spreadsheet still showed
the right figure on screen, and the next total over that column refused
the whole thing:

    column 'Betrag' cannot be totalled as it stands: the column mixes
    conventions: 13 value(s) look like decimal_point, 1 like decimal_comma

Which is the shape the office prompt names in its first paragraph —
nothing here fails loudly, a wrong number raises no exception — arriving
through the tool that exists to keep it from happening.

Found by driving edit_sheet on the workbook fixture, because no benchmark
task had ever asked either model to change a record: all twelve office
tasks read. The tool had never been called in any recorded run.

What is pinned here: the amount goes in as a NUMBER whichever convention
it was written in; an amount that reads two ways is refused rather than
guessed, for the same reason the totals refuse it; deliberate text like
"n/a" is still allowed and now says what it will cost; and a batch
carrying one bad value changes nothing at all.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

_FIXTURE = (Path(__file__).resolve().parents[1] / "tests" / "fixtures"
            / "office_workspace" / "Buchungen_2026.xlsx")

# R-014 in the generated workbook, and the whole column's total.
_R014_BEFORE = 265.85
_TOTAL_BEFORE = 23716.25


@pytest.fixture
def sheet(tmp_path):
    if not _FIXTURE.exists():
        pytest.skip("the workbook fixture is not built in this checkout")
    shutil.copy(_FIXTURE, tmp_path / _FIXTURE.name)
    return tmp_path / _FIXTURE.name


def _set(sheet, value, *, key="R-014", column="Betrag"):
    from delfin.agent import office

    return office.edit_sheet(
        sheet, key_column="Beleg",
        updates=[{"key": key, "set": {column: value}}])


def _total(sheet):
    from delfin.agent import office

    return office.sum_column(sheet, "Betrag")


def test_the_premise_the_column_holds_numbers(sheet):
    res = _total(sheet)
    assert res["total"] == pytest.approx(_TOTAL_BEFORE)
    assert res["skipped"] == []


@pytest.mark.parametrize("written", ["1.265,85", "1,265.85", "1265.85"])
def test_an_amount_goes_in_as_a_number_however_it_was_written(sheet, written):
    _set(sheet, written)
    res = _total(sheet)
    assert res["total"] == pytest.approx(_TOTAL_BEFORE - _R014_BEFORE + 1265.85)
    assert res["skipped"] == [], (
        "the edited row is no longer readable as a number")


def test_a_real_number_is_still_passed_straight_through(sheet):
    _set(sheet, 1265.85)
    assert _total(sheet)["total"] == pytest.approx(
        _TOTAL_BEFORE - _R014_BEFORE + 1265.85)


def test_an_amount_that_reads_two_ways_is_refused_not_guessed(sheet):
    from delfin.agent.office import OfficeError

    with pytest.raises(OfficeError) as exc:
        _set(sheet, "1.265")
    message = str(exc.value)
    assert "two different numbers" in message
    # Both readings named, so the caller can settle it.
    assert "1.265,85" in message or "1265.85" in message
    assert _total(sheet)["total"] == pytest.approx(_TOTAL_BEFORE), (
        "the file was changed by a call that refused")


def test_deliberate_text_is_allowed_and_says_what_it_costs(sheet):
    res = _set(sheet, "n/a")
    joined = " ".join(res.get("notes") or [])
    assert "not a number" in joined
    assert "unreadable" in joined, joined
    after = _total(sheet)
    assert after["total"] == pytest.approx(_TOTAL_BEFORE - _R014_BEFORE)
    assert after["skipped"], "the row should now be named as unreadable"


def test_one_bad_value_leaves_the_whole_batch_unwritten(sheet):
    """The module's standing promise: a half-applied batch on someone's
    records is worse than a refused one, because nothing about the file
    afterwards says which half went through."""
    from delfin.agent import office
    from delfin.agent.office import OfficeError

    with pytest.raises(OfficeError):
        office.edit_sheet(
            sheet, key_column="Beleg",
            updates=[{"key": "R-013", "set": {"Betrag": "500,00"}},
                     {"key": "R-014", "set": {"Betrag": "1.265"}}])
    assert _total(sheet)["total"] == pytest.approx(_TOTAL_BEFORE), (
        "the first update of a refused batch was written")


def test_a_text_column_is_left_alone(sheet):
    """The coercion asks what the column IS. A name is not a number even
    when it could be parsed as one."""
    from delfin.agent import office

    office.edit_sheet(sheet, key_column="Beleg",
                      updates=[{"key": "R-014", "set": {"Kreditor": "1.265"}}])
    rows = office.read_document(sheet)
    assert "1.265" in json.dumps(rows, default=str)


def test_a_formula_is_still_written_as_a_formula(sheet):
    from delfin.agent import office

    res = office.edit_sheet(
        sheet, key_column="Beleg",
        updates=[{"key": "R-014", "set": {"Betrag": "=E14*2"}}])
    assert any("formula" in n.lower() for n in (res.get("notes") or [])), (
        res.get("notes"))


def test_the_same_holds_for_an_edit_addressed_by_cell(sheet):
    """The column does not care how the caller addressed the row."""
    from delfin.agent import office

    office.edit_sheet(sheet, edits=[{"cell": "E15", "value": "1.265,85"}])
    res = _total(sheet)
    assert res["total"] == pytest.approx(
        _TOTAL_BEFORE - _R014_BEFORE + 1265.85)
    assert res["skipped"] == []


# ---------------------------------------------------------------------------
# The other two tools nobody had ever called
# ---------------------------------------------------------------------------
#
# fill_series and draft_email were driven by hand at the same time and
# both did their job. draft_email had one thing worth changing: the
# Message-ID it generated ended in this machine's fully-qualified name,
# and the draft is a file the user opens and SENDS. The internal name of a
# shared compute host is not the sender's mail domain and travels to
# whoever receives the mail. The id only has to be unique.

def test_a_draft_does_not_carry_this_machines_name_to_the_recipient(tmp_path):
    import socket

    from delfin.agent import office

    office.draft_email(tmp_path / "anfrage.eml", to="einkauf@example.org",
                       subject="Rückfrage", body="Guten Tag")
    text = (tmp_path / "anfrage.eml").read_text(encoding="utf-8")
    host = socket.getfqdn()
    assert host not in text, f"the draft names this host: {host}"
    if "." in host:
        assert host.split(".", 1)[1] not in text, "the domain still travels"


def test_the_draft_is_still_a_valid_message_with_a_unique_id(tmp_path):
    """The id is reserved-TLD, not absent — a message without one is not
    a message, and two drafts must not collide."""
    from email import message_from_string

    from delfin.agent import office

    ids = []
    for name in ("a.eml", "b.eml"):
        office.draft_email(tmp_path / name, to="x@example.org",
                           subject="Test", body="Hallo")
        msg = message_from_string(
            (tmp_path / name).read_text(encoding="utf-8"))
        assert msg["To"] == "x@example.org"
        assert msg["Subject"] == "Test"
        ids.append(msg["Message-ID"])
    assert all(i and i.startswith("<") and i.endswith(">") for i in ids), ids
    assert ids[0] != ids[1], "two drafts share a Message-ID"


def test_a_series_fills_every_row_and_keeps_the_numbers_as_written(tmp_path):
    """fill_series, driven for the first time. Recorded because it is the
    tool the prompt's whole series section is about and nothing had ever
    run it."""
    from delfin.agent import office

    table = tmp_path / "kunden.csv"
    table.write_text("Nr;Name;Betrag\n1;Meier GmbH;1.234,50\n"
                     "2;Schulze AG;289,90\n", encoding="utf-8")
    template = tmp_path / "brief.docx"
    office.create_docx(template, [
        {"heading": "Zahlungserinnerung"},
        {"paragraph": "Sehr geehrte {{Name}},"},
        {"paragraph": "offen sind {{Betrag}} EUR."}])

    out = tmp_path / "briefe"
    res = office.fill_series(table, template, output_dir=out,
                             name_pattern="Brief_{Nr}.docx")
    assert res["counts"] == {"ok": 2, "incomplete": 0, "failed": 0}
    assert sorted(p.name for p in out.iterdir()) == [
        "Brief_1.docx", "Brief_2.docx"]

    body = json.dumps(office.read_document(out / "Brief_1.docx"),
                      ensure_ascii=False, default=str)
    assert "Meier GmbH" in body
    assert "{{" not in body, "a placeholder survived into the letter"
    assert "1.234,50" in body, "the amount was reformatted on the way in"
