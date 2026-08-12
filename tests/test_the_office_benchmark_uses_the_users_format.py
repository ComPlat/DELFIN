"""The office suite scored a model against four small CSV files.

That is not the format the work arrives in. Real administrative files are
workbooks: a merged block for a heading, a filter somebody left switched
on, a total row at the bottom, a money column typed by hand, and more
rows than any reader shows at once. Every trap the office module learned
to detect this year lives in that format — and not one of them could fire
against a CSV, so the suite was measuring the parts that were never in
question while the new checks sat unexercised.

Committing the workbooks was the obvious fix and the wrong one: a binary
fixture is a file nobody can review in a diff, and the workspace README
says so. So the content lives in a spec that reads as a table and the
runner materialises it before the snapshot is taken.

Three properties that have to hold, and one that has to be loud:

* the workbooks are byte-stable, or a baseline compares two different
  files and calls the difference a model regression;
* a changed spec rebuilds them, or a run scores against a file the tasks
  no longer describe;
* they are built BEFORE the fixture snapshot, because the guard restores
  the workspace to whatever it held when the snapshot ran — a workbook
  written afterwards is deleted again after the first attempt;
* and when they cannot be built at all, the run says so. A task that
  fails for a missing library otherwise scores as a failing model.
"""

from __future__ import annotations

import pytest

from delfin.agent import benchmark_fixtures as BF

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / BF.OFFICE_WS_REL
    ws.mkdir(parents=True)
    return tmp_path, ws


# ---------------------------------------------------------------------------
# They get built, and they carry the traps
# ---------------------------------------------------------------------------

def test_the_workbooks_are_written(workspace):
    root, ws = workspace
    written, reason = BF.ensure_office_fixtures(root)
    assert reason == ""
    assert set(written) == set(BF.WORKBOOK_NAMES)
    for name in BF.WORKBOOK_NAMES:
        assert (ws / name).is_file()


def test_the_grouping_lives_in_merged_cells(workspace):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    wb = openpyxl.load_workbook(ws / "Kostenstellen.xlsx")
    sheet = wb.active
    assert len(sheet.merged_cells.ranges) >= 4
    # The row below a merged cost centre is blank in its own right --
    # that is the whole trap.
    assert sheet.cell(row=4, column=1).value is None
    assert sheet.cell(row=3, column=1).value == 4711


def test_a_filter_hides_rows_that_still_count(workspace):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    sheet = openpyxl.load_workbook(ws / "Buchungen_2026.xlsx").active
    hidden = [i for i in range(1, sheet.max_row + 1)
              if sheet.row_dimensions[i].hidden]
    assert len(hidden) >= 5


def test_the_total_carries_no_computed_value(workspace):
    """openpyxl never evaluates, so the number the row promises is not in
    the file. A model that reports it as computed is wrong."""
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    cached = openpyxl.load_workbook(
        ws / "Buchungen_2026.xlsx", data_only=True).active
    formulas = openpyxl.load_workbook(ws / "Buchungen_2026.xlsx").active
    last = formulas.max_row
    assert str(formulas.cell(row=last, column=5).value).startswith("=")
    assert cached.cell(row=last, column=5).value is None


def test_the_money_column_mixes_two_minus_conventions(workspace):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    sheet = openpyxl.load_workbook(ws / "Gutschriften.xlsx").active
    values = [sheet.cell(row=r, column=4).value
              for r in range(2, sheet.max_row + 1)]
    assert any(str(v).endswith("-") for v in values)
    assert any(str(v).startswith("(") for v in values)
    assert all(isinstance(v, str) for v in values), (
        "the column has to be text, or the trap is not a trap")


def test_the_maximum_sits_outside_the_default_window(workspace):
    from delfin.agent import office as O
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    sheet = openpyxl.load_workbook(ws / "Verbrauch_2026.xlsx").active
    column = [(r, sheet.cell(row=r, column=3).value)
              for r in range(2, sheet.max_row + 1)]
    peak_row = max(column, key=lambda rv: rv[1])[0]
    assert sheet.max_row > O.DEFAULT_MAX_ROWS
    assert peak_row > O.DEFAULT_MAX_ROWS, (
        "the point of the file is that a default read cannot see the answer")


# ---------------------------------------------------------------------------
# Stable, and rebuilt when the spec changes
# ---------------------------------------------------------------------------

def test_building_twice_gives_the_same_bytes(workspace):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    first = {n: (ws / n).read_bytes() for n in BF.WORKBOOK_NAMES}
    (ws / BF.STAMP_NAME).unlink()
    BF.ensure_office_fixtures(root)
    for name, blob in first.items():
        assert (ws / name).read_bytes() == blob, (
            f"{name} is not reproducible — a baseline would compare two "
            f"different files and call the difference a regression")


def test_a_second_call_does_no_work(workspace):
    root, _ws = workspace
    BF.ensure_office_fixtures(root)
    written, reason = BF.ensure_office_fixtures(root)
    assert written == [] and reason == ""


def test_a_changed_spec_rebuilds_them(workspace, monkeypatch):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    monkeypatch.setattr(BF, "spec_digest", lambda: "0123456789abcdef")
    assert not BF.fixtures_are_current(ws)
    written, _ = BF.ensure_office_fixtures(root)
    assert set(written) == set(BF.WORKBOOK_NAMES)


def test_a_missing_workbook_is_rebuilt_even_with_a_matching_stamp(workspace):
    root, ws = workspace
    BF.ensure_office_fixtures(root)
    (ws / "Gutschriften.xlsx").unlink()
    assert not BF.fixtures_are_current(ws)


# ---------------------------------------------------------------------------
# The run says when it cannot build them
# ---------------------------------------------------------------------------

def test_a_missing_dependency_gives_a_reason_not_a_silence(
    workspace, monkeypatch,
):
    root, _ws = workspace
    monkeypatch.setattr(
        BF, "missing_dependency_reason",
        lambda: "openpyxl is not importable — install the spreadsheet extra")
    written, reason = BF.ensure_office_fixtures(root)
    assert written == []
    assert "openpyxl" in reason


def test_no_workspace_is_a_reason_too(tmp_path):
    written, reason = BF.ensure_office_fixtures(tmp_path)
    assert written == [] and "office workspace" in reason


def test_the_runner_builds_them_before_it_snapshots():
    """After the snapshot is too late: the guard restores the workspace to
    what it held when it ran, so a workbook written later is deleted
    again after the first attempt, and the task then fails for a reason
    that has nothing to do with the model."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "benchmark_runner.py").read_text(encoding="utf-8")
    i = src.index("def run_task(")
    j = src.index("def run_suite(")
    assert "_prepare_office_fixtures()" in src[i:j]
    # ...and the guard is entered from _run_task_once, which run_task calls
    # afterwards.
    assert src.index("_prepare_office_fixtures()", i) < j


def test_the_runner_says_it_out_loud(capsys, monkeypatch):
    from delfin.agent import benchmark_runner as BR
    monkeypatch.setattr(BR, "_fixture_notice_shown", False)
    monkeypatch.setattr(
        BR, "ensure_office_fixtures",
        lambda *_a, **_k: ([], "openpyxl is not importable"))
    BR._prepare_office_fixtures()
    out = capsys.readouterr().out
    assert "openpyxl" in out
    assert "not because of the model" in out


def test_it_says_it_once(capsys, monkeypatch):
    from delfin.agent import benchmark_runner as BR
    monkeypatch.setattr(BR, "_fixture_notice_shown", False)
    monkeypatch.setattr(
        BR, "ensure_office_fixtures", lambda *_a, **_k: ([], "no dice"))
    BR._prepare_office_fixtures()
    BR._prepare_office_fixtures()
    assert capsys.readouterr().out.count("no dice") == 1


def test_a_broken_builder_does_not_take_the_run_down(monkeypatch):
    from delfin.agent import benchmark_runner as BR
    monkeypatch.setattr(BR, "_fixture_notice_shown", True)
    monkeypatch.setattr(
        BR, "ensure_office_fixtures",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("boom")))
    assert "boom" in BR._prepare_office_fixtures()


# ---------------------------------------------------------------------------
# The suite asks about them
# ---------------------------------------------------------------------------

def test_every_workbook_is_named_by_a_task():
    """A fixture no task reads measures nothing."""
    import pathlib
    import yaml
    tasks = yaml.safe_load(
        (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
         / "pack" / "benchmark" / "tasks_office.yaml").read_text(
            encoding="utf-8"))
    prompts = " ".join(t.get("prompt", "") for t in tasks["tasks"])
    for name in BF.WORKBOOK_NAMES:
        assert name in prompts, f"no office task reads {name}"


def test_the_workbooks_are_not_committed():
    """The whole reason they are generated: a binary cannot be reviewed."""
    import subprocess
    out = subprocess.run(
        ["git", "ls-files", "tests/fixtures/office_workspace"],
        capture_output=True, text=True, check=False).stdout
    assert ".xlsx" not in out, "a workbook was committed after all"
