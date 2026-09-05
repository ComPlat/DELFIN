"""A figure an answer states is checked against what the tools returned.

The defect, verified by execution before this existed: every guarantee in
the office path was a request to the model, not a mechanism. The tools
read, wrote and profiled correctly, and nothing compared a number an
ANSWER stated against what a tool had actually produced.

    scan_for_unsourced_quantities(
        "Die Gesamtsumme beträgt 45231.50 EUR bei 31 Belegen.",
        observed_files=set(), evidence_tools_used=set())  ->  []

Zero evidence acts, a total and a count asserted, nothing flagged. The
quantity guard covers eV, kcal/mol, Hartree and nm — the chemistry units.
No currency, no counts, no percentages. So a total that silently dropped a
row read exactly like a total that did not, which in administrative work
is the whole failure mode.

What is pinned here: the office tools record the figures they return, and
a finished answer is read against that ledger. A figure the tools DID
produce is passed in silence — a caveat on every answer is a caveat nobody
reads — and so is one the user typed, one from earlier in the
conversation, and one derivable from two tool figures. What is left is
named, with the reason, and never blocked.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import office
from delfin.agent import verify_guard as vg

openpyxl = pytest.importorskip("openpyxl")

BAD_ANSWER = "Die Gesamtsumme beträgt 45231.50 EUR bei 31 Belegen."


@pytest.fixture(autouse=True)
def clean_ledger():
    office.reset_figure_ledger()
    yield
    office.reset_figure_ledger()


@pytest.fixture
def workspace(tmp_path):
    """Two small German tables, written the way the users' files are."""
    left = tmp_path / "buchungen.xlsx"
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Buchungen"
    rows = [
        ("Beleg", "Kostenstelle", "Betrag"),
        ("R-001", "4711", "1.234,50"),
        ("R-002", "4712", "289,90"),
        ("R-003", "4711", "450,00"),
        ("R-004", "4713", "n/a"),
        ("R-005", "4711", "780,00"),
    ]
    for row in rows:
        sheet.append(row)
    wb.save(left)
    wb.close()
    return tmp_path, left


def _sum(path, column, **kw):
    result = office.sum_column(path, column, **kw)
    office.record_tool_figures("sum_column", result)
    return result


# ---------------------------------------------------------------------------
# The gap this closes
# ---------------------------------------------------------------------------

def test_the_existing_quantity_guard_says_nothing_about_a_money_total():
    """The measurement the mechanism was built from, kept as the record."""
    assert vg.scan_for_unsourced_quantities(
        BAD_ANSWER, observed_files=set(), evidence_tools_used=set()) == []


def test_a_total_and_a_count_after_zero_tool_calls_are_both_named():
    flags = office.scan_answer_for_unledgered_figures(BAD_ANSWER)
    assert [(f.figure, f.kind) for f in flags] == [
        ("45231.50", "total"), ("31", "count")]
    # In the order they are written, so the caveat reads like the answer.
    caveat = office.figure_caveat(flags)
    assert "45231.50" in caveat and "31" in caveat


def test_the_whole_check_is_one_call_for_the_caller_that_appends_it():
    """What the engine wires in: one call, "" for an answer whose figures
    the tools produced — which is nearly every answer."""
    assert office.figure_coverage_caveat(BAD_ANSWER).startswith("\n\n>")
    assert office.figure_coverage_caveat("Die Datei ist gelesen.") == ""
    assert office.figure_coverage_caveat(
        BAD_ANSWER, user_text="45231.50 EUR bei 31 Belegen, bitte prüfen"
    ) == ""


def test_the_caveat_annotates_and_does_not_block():
    """An answer naming a figure the tools did not produce is often
    legitimate. The mechanism marks it; it never withholds the answer."""
    assert office.figure_caveat([]) == ""
    caveat = office.figure_caveat(
        office.scan_answer_for_unledgered_figures(BAD_ANSWER))
    assert caveat.startswith("\n\n>")
    # It says which figure and why — the reader is the one who can tell an
    # invented total from one they typed themselves.
    assert "sum_column" in caveat
    assert "Werkzeug-Ergebnis" in caveat


# ---------------------------------------------------------------------------
# A figure a tool produced is passed in silence
# ---------------------------------------------------------------------------

def test_a_total_the_tool_computed_is_not_flagged(workspace):
    _root, book = workspace
    result = _sum(book, "Betrag")
    assert result["total"] == 2754.4
    answer = ("Die Summe der Spalte Betrag beträgt 2.754,40 EUR über 4 "
              "Werte.")
    assert office.scan_answer_for_unledgered_figures(answer) == []


def test_the_same_total_written_the_english_way_is_the_same_figure(workspace):
    _root, book = workspace
    _sum(book, "Betrag")
    for written in ("2754.40", "2.754,40", "2 754,40", "2754,4"):
        answer = f"Die Gesamtsumme beträgt {written} EUR."
        assert office.scan_answer_for_unledgered_figures(answer) == [], written


def test_a_thousands_group_reads_both_ways_rather_than_being_flagged():
    """'1.234' is 1234 to a German writer and 1.234 to an English one, and
    nothing in an answer settles it. Both readings count as a match, or a
    correctly written German total is reported as unsupported."""
    ledger = [office.ToolFigure(1234.0, "sum", "Summe", "sum_column")]
    for written in ("1.234", "1234", "1,234"):
        assert office.scan_answer_for_unledgered_figures(
            f"Die Summe beträgt {written} EUR.", ledger=ledger) == [], written


def test_a_wrong_total_over_the_same_column_is_still_named(workspace):
    _root, book = workspace
    _sum(book, "Betrag")
    flags = office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 2.854,40 EUR.")
    assert [f.figure for f in flags] == ["2.854,40"]


def test_a_total_is_matched_to_the_precision_it_was_written_to():
    ledger = [office.ToolFigure(12345.674, "sum", "Summe", "sum_column")]
    assert office.scan_answer_for_unledgered_figures(
        "Die Summe beträgt 12.345,67 EUR.", ledger=ledger) == []
    assert office.scan_answer_for_unledgered_figures(
        "Die Summe beträgt 12.345,68 EUR.", ledger=ledger)


# ---------------------------------------------------------------------------
# Counts come from the tool that counted
# ---------------------------------------------------------------------------

def test_the_counts_a_reconciliation_produced_ground_the_answer(tmp_path):
    left = tmp_path / "a.csv"
    right = tmp_path / "b.csv"
    left.write_text("Beleg;Betrag\nR-001;10,00\nR-002;20,00\nR-003;30,00\n",
                    encoding="utf-8")
    right.write_text("Beleg;Betrag\nR-001;10,00\nR-002;25,00\nR-004;40,00\n",
                     encoding="utf-8")
    result = office.compare_tables(left, right, key="Beleg")
    # The comparison states its figures outright, so the ledger does not
    # have to read them back out of a rendered report.
    assert {"kind", "value", "label"} <= set(result["figures"][0])
    office.record_tool_figures("compare_tables", result)

    answer = ("Der Abgleich ergibt 1 gleiche Zeile, 1 Abweichung, "
              "1 Beleg nur links und 1 Beleg nur rechts.")
    assert office.scan_answer_for_unledgered_figures(answer) == []
    flags = office.scan_answer_for_unledgered_figures(
        "Der Abgleich ergibt 7 Abweichungen.")
    assert [f.figure for f in flags] == ["7"]


def test_a_count_of_documents_read_is_grounded_by_the_reads(workspace):
    _root, book = workspace
    office.record_tool_figures("read_document", office.read_document(book))
    answer = "Ich habe 1 Datei gelesen und daraus berichtet."
    assert office.scan_answer_for_unledgered_figures(answer) == []


def test_a_turn_without_any_office_tool_only_reads_money():
    """The documented scope of the empty-ledger case. A coding session
    says "Insgesamt 3 Tests sind fehlgeschlagen" — a total cue and a
    count, and none of this mechanism's business. What it does read is a
    sentence stating an amount of money, because that is the claim the
    office tools exist to produce."""
    for answer in (
        "Insgesamt 3 Tests sind fehlgeschlagen.",
        "Ich habe 12 Dateien im Ordner gefunden.",
        "Der Anteil der fehlenden Docstrings liegt bei 12 %.",
        "Die Suite läuft in 355 Sekunden durch; insgesamt 9.029 Tests.",
    ):
        assert office.scan_answer_for_unledgered_figures(answer) == [], answer
    assert office.scan_answer_for_unledgered_figures(
        "Der Gesamtbetrag beträgt 12.345,67 EUR.")


def test_a_figure_is_named_after_the_more_specific_cue():
    ledger = [office.ToolFigure(1.0, "count", "1 Datei", "read_document")]
    flags = office.scan_answer_for_unledgered_figures(
        "Die Differenz zum Vorjahr beträgt 8.412,00 EUR.", ledger=ledger)
    assert [(f.figure, f.kind) for f in flags] == [("8.412,00", "derived")]
    assert "abgeleiteter Wert" in office.figure_caveat(flags)


# ---------------------------------------------------------------------------
# A derived figure is derived, not invented
# ---------------------------------------------------------------------------

def test_a_difference_between_two_tool_totals_counts_as_grounded():
    ledger = [
        office.ToolFigure(23716.25, "sum", "Summe 2026", "sum_column"),
        office.ToolFigure(19500.00, "sum", "Summe 2025", "sum_column"),
    ]
    answer = "Die Differenz zwischen beiden Jahren beträgt 4.216,25 EUR."
    assert office.scan_answer_for_unledgered_figures(
        answer, ledger=ledger) == []


def test_a_share_of_two_tool_counts_counts_as_grounded():
    ledger = [
        office.ToolFigure(31.0, "count", "31 Belege", "compare_tables"),
        office.ToolFigure(148.0, "count", "148 Zeilen", "read_document"),
    ]
    answer = "Das entspricht einem Anteil von 20,9 % der Belege."
    assert office.scan_answer_for_unledgered_figures(
        answer, ledger=ledger) == []


def test_a_derived_figure_that_derives_from_nothing_is_named():
    ledger = [
        office.ToolFigure(23716.25, "sum", "Summe 2026", "sum_column"),
        office.ToolFigure(25.0, "count", "25 Zeilen", "read_document"),
    ]
    flags = office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 8.412,00 EUR.", ledger=ledger)
    assert [f.figure for f in flags] == ["8.412,00"]


# ---------------------------------------------------------------------------
# Sources other than a tool
# ---------------------------------------------------------------------------

def test_a_figure_the_user_typed_is_grounded():
    assert office.scan_answer_for_unledgered_figures(
        "Der Gesamtbetrag von 45.231,50 EUR ist im Bescheid genannt.",
        user_text="Bitte prüfe die 45.231,50 EUR aus dem Bescheid.") == []


def test_a_figure_restated_from_earlier_in_the_conversation_is_grounded():
    assert office.scan_answer_for_unledgered_figures(
        "Wie oben: die Gesamtsumme beträgt 45.231,50 EUR.",
        prior_text="Die Summe war 45.231,50 EUR.") == []


def test_a_figure_the_answer_does_not_assert_is_left_alone():
    """Konjunktiv is how German says 'not measured'. Flagging it would
    caveat exactly the honest answer."""
    for answer in (
        "Die Summe wäre 25.136 EUR, wenn die Spalte deutsch gelesen wird.",
        "Die Gesamtsumme liegt bei ca. 45.000 EUR.",
        "Beispiel: eine Gesamtsumme von 1.000,00 EUR.",
        "Der Gesamtbetrag ist noch nicht berechnet — 0,00 EUR steht nur da.",
    ):
        assert office.scan_answer_for_unledgered_figures(answer) == [], answer


def test_a_figure_inside_code_is_not_the_answer_asserting_it():
    for answer in (
        "Die Summe kommt aus `=SUMME(D2:D26)` und ergibt 45.231,50 EUR.",
        "```\nGesamtsumme: 45231.50 EUR\n```",
        "> Gesamtsumme: 45231.50 EUR",
    ):
        flags = office.scan_answer_for_unledgered_figures(answer)
        assert all(f.figure != "45231.50" for f in flags), answer


# ---------------------------------------------------------------------------
# The ledger is per turn, and bounded
# ---------------------------------------------------------------------------

def test_the_ledger_is_forgotten_between_turns(workspace):
    _root, book = workspace
    _sum(book, "Betrag")
    assert office.figure_ledger()
    office.reset_figure_ledger()
    assert office.figure_ledger() == []
    # And a total grounded a moment ago is ungrounded again, which is the
    # point: it belonged to the turn that computed it.
    assert office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 2.754,40 EUR.")


def test_the_ledger_cannot_grow_without_limit(workspace):
    _root, book = workspace
    for _ in range(400):
        office.record_tool_figures(
            "read_document", {"path": str(book), "rows": 3, "columns": 2,
                              "numbers": list(range(50))})
    assert len(office.figure_ledger()) <= office.MAX_LEDGER_FIGURES + 1


def test_a_result_that_is_not_a_dict_is_ignored_rather_than_raising():
    assert office.record_tool_figures("read_document", "some text") == []
    assert office.record_tool_figures("read_document", None) == []


def test_the_scan_never_raises_on_junk():
    for text in ("", None, "…", "1.2.3.4.5", "€€€", "12,345,678,901"):
        assert isinstance(
            office.scan_answer_for_unledgered_figures(text), list)


# ---------------------------------------------------------------------------
# The false-positive rate on answers a careful person would write
# ---------------------------------------------------------------------------

FIXTURE_WS = pathlib.Path("tests/fixtures/office_workspace")


@pytest.fixture
def fixtures():
    from delfin.agent.benchmark_fixtures import ensure_office_fixtures

    _written, reason = ensure_office_fixtures()
    if reason:
        pytest.skip(reason)
    return FIXTURE_WS


def _fixture_ledger(root: pathlib.Path) -> None:
    """The tool calls a careful answer about these files rests on."""
    for name in ("Buchungen_2026.xlsx", "Gutschriften.xlsx",
                 "Verbrauch_2026.xlsx"):
        office.record_tool_figures(
            "read_document", office.read_document(root / name))
    office.record_tool_figures(
        "read_document", office.read_document(root / "Kostenstellen.xlsx"))
    _sum(root / "Buchungen_2026.xlsx", "Betrag")
    _sum(root / "Buchungen_2026.xlsx", "Betrag", group_by="Kostenstelle")
    _sum(root / "Gutschriften.xlsx", "Betrag")
    _sum(root / "Kostenstellen.xlsx", "Budget", header_row=2)
    office.record_tool_figures(
        "compare_tables",
        office.compare_tables(root / "buchungen.csv",
                              root / "rechnungen.csv", key="Beleg",
                              right_key="Belegnummer",
                              columns={"Betrag": "Rechnungsbetrag"}))


# Ten answers a careful person would give about these workbooks, written
# the way the users write: German figures, per-group totals, counts of
# things, and the caveats the tools themselves reported. None of them may
# be flagged — that is the measurement, not an aspiration.
REALISTIC_ANSWERS = [
    "Die Summe der Spalte Betrag in Buchungen_2026.xlsx beträgt 23.716,25 "
    "EUR. Gezählt wurden 24 Werte; eine Zeile hat keinen Betrag.",

    "Achtung: 7 Zeilen der Datei sind ausgeblendet. Sie sind im "
    "Gesamtbetrag von 23.716,25 EUR enthalten, den Sie am Bildschirm so "
    "nicht sehen.",

    "Summe je Kostenstelle: 4711 = 11.735,00 EUR, 4712 = 10.843,95 EUR, "
    "4713 = 1.137,30 EUR.",

    "Die offenen Posten summieren sich auf 17.076,70 EUR. 4 von 10 Werten "
    "sind Gutschriften in einer Schreibweise, die nicht als Zahl gelesen "
    "werden kann, und fehlen in dieser Summe.",

    "Buchungen_2026.xlsx hat 26 Zeilen und 5 Spalten; Verbrauch_2026.xlsx "
    "hat 261 Zeilen.",

    "Der Abgleich von buchungen.csv gegen rechnungen.csv ergibt 3 gleiche "
    "Belege, 1 Abweichung und 2 nicht vergleichbare Zeilen.",

    "Das Gesamtbudget über alle Kostenstellen beträgt 89.050,00 EUR, davon "
    "18.500,00 EUR auf Verbrauchsmaterial der Kostenstelle 4711.",

    "Der Anteil der Kostenstelle 4711 am Gesamtbudget liegt bei 20,8 %.",

    "Ich habe 4 Dateien gelesen: Buchungen_2026.xlsx, Gutschriften.xlsx, "
    "Verbrauch_2026.xlsx und Kostenstellen.xlsx.",

    "Zusammenfassung: 23.716,25 EUR auf 24 Buchungen, 17.076,70 EUR offene "
    "Posten, Differenz 6.639,55 EUR.",

    "Die Spalte Anschaffungswert in inventar.csv ist nicht eindeutig "
    "lesbar (8.986 bedeutet 8986 oder 8,986), deshalb wurde keine Summe "
    "gebildet.",

    "Kostenstelle 4711 trägt mit 11.735,00 EUR den größten Teil der "
    "Gesamtsumme von 23.716,25 EUR.",
]


@pytest.mark.parametrize("answer", REALISTIC_ANSWERS,
                         ids=range(len(REALISTIC_ANSWERS)))
def test_a_careful_answer_about_the_fixtures_is_not_caveated(
        fixtures, answer):
    _fixture_ledger(fixtures)
    assert office.scan_answer_for_unledgered_figures(answer) == [], answer


def test_the_false_positive_rate_over_all_of_them_is_zero(fixtures):
    """Measured, not asserted per case: the whole set at once, because a
    caveat on every answer is a caveat nobody reads."""
    _fixture_ledger(fixtures)
    flagged = [a for a in REALISTIC_ANSWERS
               if office.scan_answer_for_unledgered_figures(a)]
    assert flagged == []


def test_the_same_ledger_still_catches_a_figure_nothing_produced(fixtures):
    """The other half of the measurement: silence on correct answers is
    only worth something if the mechanism still speaks up."""
    _fixture_ledger(fixtures)
    flags = office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme aller Buchungen beträgt 31.402,80 EUR.")
    assert [f.figure for f in flags] == ["31.402,80"]


def test_no_money_total_slips_through_by_coincidence(fixtures):
    """The sharpness of the total check, measured rather than assumed:
    over the fixture ledger (295 figures, 1370 derivable values) not one
    of 300 random amounts is grounded by accident. A guard that says yes
    to anything is not a guard, and the tolerance — half a unit in the
    last place written — is what keeps it that way."""
    import random

    _fixture_ledger(fixtures)
    rng = random.Random(20260812)
    passed = [
        value for value in
        (round(rng.uniform(100, 100_000), 2) for _ in range(300))
        if not office.scan_answer_for_unledgered_figures(
            f"Die Gesamtsumme beträgt {value:.2f} EUR.".replace(".", ","))
    ]
    assert passed == []
