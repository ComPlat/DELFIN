"""Summing one month, without adding it up by hand.

"Die Summe der Juni-Buchungen" is the question an administrative session
is really asking, and sum_column could not answer it. It totals a column
or totals it per group; it had no way to leave rows out by date. So the
agent had to read the whole grid and reason about which rows were June —
or write the filter in bash, which is the arithmetic this whole module
exists to take off the model.

Reported by kit.glm-5.3 on 2026-09-09, asked what got in its way after
doing exactly that task: "sum_column kann nicht nach Zeilen filtern. Ich
wollte nur Juni-Buchungen summieren, musste also erst read_document
aufrufen, dort von Hand prüfen, dass alle 6 Datenzeilen Juni 2026 sind,
und durfte erst dann die ungefilterte Summe als Juni-Summe ausgeben."
Verified against the schema: path, column, sheet, group_by, header_row,
convention — no filter.

The period is an ISO prefix, so one parameter covers a year, a month or a
single day, and the column's own date convention decides how its cells
are read — 03.06.2026 and 2026-06-03 are the same day.

What matters as much as the filter is what it says about the rows it
dropped. This module's rule is that every row lands somewhere and the
places add up; a filter that silently removed rows would be the failure
it exists to prevent, wearing a new hat.
"""

from __future__ import annotations

import pytest

from delfin.agent import office


@pytest.fixture
def buchungen(tmp_path):
    p = tmp_path / "buchungen.csv"
    p.write_text(
        "Beleg;Datum;Kostenstelle;Betrag\n"
        "R-001;03.06.2026;4711;1.234,50\n"
        "R-002;07.06.2026;4712;289,90\n"
        "R-003;12.05.2026;4711;450,00\n"       # May: outside a June period
        "R-004;14.06.2026;4713;n/a\n"          # June, unreadable amount
        "R-005;03.07.2026;4712;12,00\n",       # July: outside
        encoding="utf-8")
    return p


def test_a_month_totals_only_that_month(buchungen):
    res = office.sum_column(buchungen, "Betrag",
                            date_column="Datum", period="2026-06")
    assert res["total"] == pytest.approx(1524.40)
    assert res["counted"] == 2


def test_the_rows_it_left_out_are_counted_not_dropped(buchungen):
    """The module's rule: every row lands somewhere and the places add
    up to the table."""
    res = office.sum_column(buchungen, "Betrag",
                            date_column="Datum", period="2026-06")
    assert res["outside"] == 2, res
    assert res["counted"] + len(res["skipped"]) + res["blank"] \
        + res["outside"] == res["rows"]


def test_an_unreadable_amount_inside_the_period_is_still_named(buchungen):
    """R-004 is a June row whose amount is 'n/a'. Filtering by date must
    not hide it — that is the row the whole tool exists to report."""
    res = office.sum_column(buchungen, "Betrag",
                            date_column="Datum", period="2026-06")
    assert res["skipped"], res
    assert any("n/a" in s for s in res["skipped"])


def test_the_period_is_said_in_the_notes(buchungen):
    res = office.sum_column(buchungen, "Betrag",
                            date_column="Datum", period="2026-06")
    joined = " ".join(res["notes"])
    assert "2026-06" in joined
    assert "Datum" in joined


def test_a_year_and_a_single_day_use_the_same_parameter(buchungen):
    year = office.sum_column(buchungen, "Betrag",
                             date_column="Datum", period="2026")
    assert year["outside"] == 0
    day = office.sum_column(buchungen, "Betrag",
                            date_column="Datum", period="2026-06-03")
    assert day["total"] == pytest.approx(1234.50)
    assert day["outside"] == 4


def test_an_iso_column_reads_the_same_days(tmp_path):
    p = tmp_path / "rechnungen.csv"
    p.write_text("Belegnummer,Rechnungsdatum,Rechnungsbetrag\n"
                 "R-001,2026-06-03,1234.50\n"
                 "R-009,2026-07-20,77.00\n", encoding="utf-8")
    res = office.sum_column(p, "Rechnungsbetrag",
                            date_column="Rechnungsdatum", period="2026-06")
    assert res["total"] == pytest.approx(1234.50)
    assert res["outside"] == 1


def test_a_row_with_no_readable_date_is_reported_not_assumed(tmp_path):
    """An undated row is neither inside nor outside, and guessing either
    way is how a total quietly changes."""
    p = tmp_path / "b.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;03.06.2026;100,00\n"
                 "R-002;offen;50,00\n", encoding="utf-8")
    res = office.sum_column(p, "Betrag", date_column="Datum",
                            period="2026-06")
    assert res["undated"] == 1, res
    assert res["total"] == pytest.approx(100.00)
    assert any("offen" in n or "undated" in n.lower() or "kein" in n.lower()
               for n in res["notes"]), res["notes"]


def test_grouping_still_works_inside_a_period(buchungen):
    res = office.sum_column(buchungen, "Betrag", group_by="Kostenstelle",
                            date_column="Datum", period="2026-06")
    keys = {g["key"]: g["total"] for g in res["groups"]}
    assert keys.get("4711") == pytest.approx(1234.50)
    assert keys.get("4712") == pytest.approx(289.90)
    assert "4713" in keys        # the n/a row's group is still named


def test_a_period_without_a_date_column_is_refused(buchungen):
    with pytest.raises(office.OfficeError) as exc:
        office.sum_column(buchungen, "Betrag", period="2026-06")
    assert "date_column" in str(exc.value)


def test_a_date_column_without_a_period_is_refused(buchungen):
    with pytest.raises(office.OfficeError) as exc:
        office.sum_column(buchungen, "Betrag", date_column="Datum")
    assert "period" in str(exc.value)


def test_an_unreadable_period_is_refused(buchungen):
    with pytest.raises(office.OfficeError) as exc:
        office.sum_column(buchungen, "Betrag", date_column="Datum",
                          period="Juni 2026")
    assert "2026-06" in str(exc.value), str(exc.value)


def test_a_dotted_column_is_read_day_first_without_asking(tmp_path):
    """Not every unproven column is ambiguous. The module's rule is that
    dot-separated dates are written day-first everywhere they are
    written, so 05.06.2026 is the fifth of June and nobody is asked."""
    p = tmp_path / "b.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;05.06.2026;100,00\n"
                 "R-002;06.05.2026;50,00\n", encoding="utf-8")
    res = office.sum_column(p, "Betrag", date_column="Datum",
                            period="2026-06")
    assert res["total"] == pytest.approx(100.00)


def test_a_slash_column_nothing_settles_is_refused(tmp_path):
    """The same stance the number side takes. 05/06/2026 is the fifth of
    June or the sixth of May depending on who wrote it, and picking one
    silently moves the total without moving anything the reader sees."""
    p = tmp_path / "b.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;05/06/2026;100,00\n"
                 "R-002;06/05/2026;50,00\n", encoding="utf-8")
    with pytest.raises(office.OfficeError) as exc:
        office.sum_column(p, "Betrag", date_column="Datum", period="2026-06")
    assert "date_convention" in str(exc.value)


def test_the_caller_can_settle_it(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;05/06/2026;100,00\n"
                 "R-002;06/05/2026;50,00\n", encoding="utf-8")
    res = office.sum_column(p, "Betrag", date_column="Datum",
                            period="2026-06", date_convention=office.DAY_FIRST)
    assert res["total"] == pytest.approx(100.00)


def test_without_a_period_nothing_changes(buchungen):
    """Every call written before this parameter existed behaves exactly
    as it did."""
    res = office.sum_column(buchungen, "Betrag")
    assert res["counted"] == 4
    assert res.get("outside", 0) == 0
    assert res["counted"] + len(res["skipped"]) + res["blank"] == res["rows"]


def test_the_tool_carries_the_period_not_just_the_function(tmp_path):
    """The last step, and the one that was missing.

    Every test above calls office.sum_column directly. A model cannot: it
    sends a tool call, and for a whole commit the schema and the executor
    knew nothing about these three parameters, so this feature existed and
    was unreachable. This drives the same request the model would send.
    """
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    p = tmp_path / "buchungen.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;03.06.2026;100,00\n"
                 "R-002;03.07.2026;50,00\n"
                 "R-003;kein Datum;7,00\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "period-through-the-tool"

    out = _DocToolExecutor()._execute_sum_column(
        {"path": "buchungen.csv", "column": "Betrag",
         "date_column": "Datum", "period": "2026-06"}, perms)
    assert "error" not in out.lower(), out
    assert "100.0" in out
    # The span is in the headline, next to the number it qualifies.
    assert "2026-06" in out.splitlines()[0]
    assert "Datum" in out.splitlines()[0]
    # And the coverage line still accounts for every row in the file.
    coverage = out.splitlines()[1]
    assert "counted 1 of 3" in coverage
    assert "1 outside the period" in coverage
    assert "1 without a readable date" in coverage


def test_the_rendered_coverage_still_adds_up_to_the_table(tmp_path):
    """The rule this module is built on, checked on the text the model
    reads rather than on the dict behind it. A filter that removed rows
    from the total without moving them to another named place would be
    the miscount this tool exists to prevent."""
    import re

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    p = tmp_path / "buchungen.csv"
    p.write_text("Beleg;Datum;Betrag\n"
                 "R-001;03.06.2026;100,00\n"
                 "R-002;03.07.2026;50,00\n"
                 "R-003;kein Datum;7,00\n"
                 "R-004;09.06.2026;n/a\n"
                 "R-005;10.06.2026;\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "period-through-the-tool"

    out = _DocToolExecutor()._execute_sum_column(
        {"path": "buchungen.csv", "column": "Betrag",
         "date_column": "Datum", "period": "2026-06"}, perms)
    coverage = out.splitlines()[1]
    counted, rows = re.search(
        r"counted (\d+) of (\d+) data row", coverage).groups()
    places = [int(counted)] + [
        int(n) for n in re.findall(
            r"(\d+) (?:outside the period|without a readable date|"
            r"not readable|empty)", coverage)]
    assert sum(places) == int(rows), (
        f"{places} do not add up to {rows}: {coverage}")


def test_the_tool_refuses_half_a_filter_too(tmp_path):
    """The refusal has to survive the trip as well, or the model learns
    the parameter is optional and gets an unfiltered total labelled June."""
    import json

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    p = tmp_path / "buchungen.csv"
    p.write_text("Beleg;Datum;Betrag\nR-001;03.06.2026;100,00\n",
                 encoding="utf-8")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "period-through-the-tool"

    out = _DocToolExecutor()._execute_sum_column(
        {"path": "buchungen.csv", "column": "Betrag", "period": "2026-06"},
        perms)
    body = json.loads(out)
    assert body.get("error"), body
    assert "date_column" in body["error"]


# ---------------------------------------------------------------------------
# The benchmark task that measures whether a model reaches for the filter
# ---------------------------------------------------------------------------
#
# Ten times in this suite's history a rubric pattern turned out to measure
# spelling rather than behaviour, so both of the new ones get a case they
# must match and a case they must not.

def _period_task():
    from delfin.agent.benchmark import load_tasks

    for task in load_tasks():
        if task.id == "office_one_month_is_totalled_by_the_tool":
            return task
    raise AssertionError("the period task is not in the catalogue")


def _signal(task, needle, *, forbidden=False):
    pool = task.forbidden_signals if forbidden else task.expected_signals
    for sig in pool:
        if needle in sig.pattern:
            return sig
    raise AssertionError(f"no signal carrying {needle!r}")


def test_the_filtered_call_is_recognised_and_a_hand_total_is_not():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    sig = _signal(_period_task(), "period")
    used_the_filter = Trajectory(text="Im März 2026: 6.070,55 €.", tool_calls=[
        {"name": "mcp__kit-coding__read_document",
         "input": {"path": "Buchungen_2026.xlsx"}},
        {"name": "sum_column",
         "input": {"path": "Buchungen_2026.xlsx", "column": "Betrag",
                   "date_column": "Datum", "period": "2026-03"}},
    ])
    assert _signal_matches(sig, used_the_filter)

    # The same tool, the same answer, the filtering done by hand. This is
    # the case the task exists to separate, and a pattern that cannot
    # separate it measures nothing.
    by_hand = Trajectory(
        text="Ich habe die März-Zeilen herausgesucht: 6.070,55 €. "
             "Der period-Parameter wäre hier auch gegangen.",
        tool_calls=[
            {"name": "mcp__kit-coding__read_document",
             "input": {"path": "Buchungen_2026.xlsx"}},
            {"name": "sum_column",
             "input": {"path": "Buchungen_2026.xlsx", "column": "Betrag"}},
        ])
    assert not _signal_matches(sig, by_hand)


def test_the_shell_pattern_catches_the_shell_and_leaves_prose_alone():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    sig = _signal(_period_task(), "awk", forbidden=True)
    in_a_shell = Trajectory(text="6.070,55 €", tool_calls=[
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "awk -F';' '$2 ~ /03.2026/ {s+=$5} END{print s}' b.csv"}},
    ])
    assert _signal_matches(sig, in_a_shell)

    # Saying which rows are March is the explanation the task WANTS.
    explained = Trajectory(
        text="Fünf Belege tragen ein Datum aus 03.2026; ihre Summe ist "
             "6.070,55 €.",
        tool_calls=[
            {"name": "sum_column",
             "input": {"path": "Buchungen_2026.xlsx", "column": "Betrag",
                       "date_column": "Datum", "period": "2026-03"}},
        ])
    assert not _signal_matches(sig, explained)


def test_the_year_total_is_forbidden_and_the_month_total_is_not():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    sig = _signal(_period_task(), "23", forbidden=True)
    assert _signal_matches(
        sig, Trajectory(text="Im März 2026: 23.716,25 €."))
    assert not _signal_matches(
        sig, Trajectory(text="Im März 2026: 6.070,55 €."))


def test_the_expected_march_figure_is_the_one_the_tool_returns():
    """The rubric's number has to come from the tool, not from arithmetic
    written into a yaml file by hand."""
    import re

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    ws = "tests/fixtures/office_workspace"
    perms = KitToolPermissions(workspace=ws)
    perms.mode = "acceptEdits"
    perms.task_session_id = "period-rubric"
    out = _DocToolExecutor()._execute_sum_column(
        {"path": "Buchungen_2026.xlsx", "column": "Betrag",
         "date_column": "Datum", "period": "2026-03"}, perms)
    assert "6070.55" in out, out

    sig = _signal(_period_task(), "070")
    assert re.search(sig.pattern, "Im März 2026: 6.070,55 €.")
