"""The figure ledger belongs to the turn that filled it.

The defect, measured before this existed: the ledger was one module-level
list for the whole process, and the module said that was safe because "a
figure another session's tool produced can only GROUND a claim, never
create a flag". Both halves were false.

  * ``reset_figure_ledger()`` cleared that one list and was driven by ONE
    engine's turn boundary. Anything else in the process starting a turn
    wiped the evidence a turn in flight was still collecting.
  * The 2000-figure cap was shared. A background sub-agent runs in the
    same process through the same recording path, and one document read
    contributes up to 600 figures — four foreign reads filled the cap and
    the parent's own ``sum_column`` result was dropped in silence.

Either way the answer came out the same: a correct, tool-computed German
total carrying ``⚠️ Diese Zahl stammt nicht aus einem Werkzeug-Ergebnis
dieses Zuges``. That is the false positive the module itself names as the
thing that makes the guard unreadable — a caveat on a right answer. And a
figure a sub-agent recorded after the parent's reset survived into the
parent's next turn, which is the stale grounding the reset exists to
prevent.

What is pinned here: a figure recorded under one token is invisible under
another, a foreign fill cannot evict a turn's own figures, and the
correct case — a total the tools really computed, through the real tool
path — is still passed in silence.
"""

from __future__ import annotations

import threading

import pytest

from delfin.agent import office
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

openpyxl = pytest.importorskip("openpyxl")

TOTAL_ANSWER = "Die Summe der Spalte Betrag beträgt 2.754,40 EUR."


@pytest.fixture(autouse=True)
def clean_ledger():
    """Each test starts outside any turn, and leaves the scope as it
    found it — a token set here would otherwise follow this thread into
    the next test file."""
    outer = office._LEDGER_SCOPE.set("")
    office.reset_figure_ledger()
    yield
    office.reset_figure_ledger()
    office._LEDGER_SCOPE.reset(outer)


@pytest.fixture
def book(tmp_path):
    """One small German table, written the way the users' files are."""
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Buchungen"
    for row in (
        ("Beleg", "Betrag"),
        ("R-001", "1.234,50"),
        ("R-002", "289,90"),
        ("R-003", "450,00"),
        ("R-004", "780,00"),
    ):
        sheet.append(row)
    path = tmp_path / "buchungen.xlsx"
    wb.save(path)
    wb.close()
    return path


def _foreign_reads(token: str, reads: int = 4) -> None:
    """What a background sub-agent's document reads put in the ledger."""
    for _ in range(reads):
        office.record_tool_figures("read_document", {
            "path": "/foreign/doc.xlsx", "rows": 600, "columns": 3,
            "numbers": list(range(1000, 1600))}, token)


# ---------------------------------------------------------------------------
# The two directions, asserted separately
# ---------------------------------------------------------------------------

def test_a_figure_recorded_under_one_token_is_invisible_under_another(book):
    office.record_tool_figures("sum_column", office.sum_column(book, "Betrag"),
                               "token-A")
    assert [f.value for f in office.figure_ledger("token-A")] == [
        pytest.approx(2754.4), 4.0, 4.0, 1.0]
    assert office.figure_ledger("token-B") == []
    # Which is the whole point: B's answer is not grounded by A's tools.
    assert office.scan_answer_for_unledgered_figures(
        TOTAL_ANSWER, token="token-B")


def test_a_foreign_fill_cannot_evict_a_turns_own_figures(book):
    """The measured false positive: four foreign reads, and the parent's
    own total is dropped and then flagged."""
    _foreign_reads("token-B")
    assert len(office.figure_ledger("token-B")) >= office.MAX_LEDGER_FIGURES

    office.record_tool_figures("sum_column", office.sum_column(book, "Betrag"),
                               "token-A")
    assert any(f.value == pytest.approx(2754.4)
               for f in office.figure_ledger("token-A"))
    assert office.scan_answer_for_unledgered_figures(
        TOTAL_ANSWER, token="token-A") == []
    assert office.figure_coverage_caveat(TOTAL_ANSWER, token="token-A") == ""


def test_a_foreign_turn_boundary_does_not_clear_this_turns_ledger(book):
    office.record_tool_figures("sum_column", office.sum_column(book, "Betrag"),
                               "token-A")

    def foreign_turn_start() -> None:
        office.begin_figure_turn("another-session")
        office.reset_figure_ledger()

    thread = threading.Thread(target=foreign_turn_start)
    thread.start()
    thread.join()

    assert office.scan_answer_for_unledgered_figures(
        TOTAL_ANSWER, token="token-A") == []


def test_a_sub_agents_figures_do_not_ground_the_parents_answer():
    """A figure recorded after the parent's turn opened must not stand in
    for a total the parent's own tools never produced."""
    parent = office.begin_figure_turn("parent-session")

    def sub_agent() -> None:                    # its own thread, own scope
        office.record_tool_figures("sum_column", {"figures": [
            {"value": 99999.99, "kind": "sum", "label": "foreign"}]})

    thread = threading.Thread(target=sub_agent)
    thread.start()
    thread.join()

    flags = office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 99.999,99 EUR.", token=parent)
    assert [f.figure for f in flags] == ["99.999,99"]


def test_a_sub_agent_recording_without_a_token_lands_in_its_own_scope():
    """No wiring needed for the isolation to hold: a thread that never
    opened a turn records against itself, not against whoever did."""
    parent = office.begin_figure_turn("parent-session")
    scopes: list[str] = []

    def sub_agent() -> None:
        scopes.append(office.figure_ledger_scope())
        _foreign_reads("")

    thread = threading.Thread(target=sub_agent)
    thread.start()
    thread.join()

    assert scopes and scopes[0] != parent
    assert office.figure_ledger(parent) == []


# ---------------------------------------------------------------------------
# The correct case is still silent — through the real tool path
# ---------------------------------------------------------------------------

def test_the_turn_that_computed_the_total_is_still_passed_in_silence(
        tmp_path, book):
    """The sequence the engine actually runs: open the turn, call the
    tool through the executor, check the finished answer against the
    token the turn was opened with."""
    token = office.begin_figure_turn("office-session")
    perms = KitToolPermissions(workspace=str(book.parent))
    perms.mode = "acceptEdits"
    perms.task_session_id = "office-session"

    out = _DocToolExecutor()._dispatch(
        "sum_column", {"path": str(book), "column": "Betrag"}, perms)
    assert "2754.4" in out

    assert office.figure_coverage_caveat(TOTAL_ANSWER, token=token) == ""
    # And with a sub-agent flooding the process meanwhile, still silent.
    _foreign_reads("some-other-turn")
    assert office.figure_coverage_caveat(TOTAL_ANSWER, token=token) == ""


def test_the_guard_still_names_a_figure_its_own_turn_never_produced(book):
    token = office.begin_figure_turn("office-session")
    office.record_tool_figures("sum_column", office.sum_column(book, "Betrag"))
    flags = office.scan_answer_for_unledgered_figures(
        "Die Gesamtsumme beträgt 31.402,80 EUR.", token=token)
    assert [f.figure for f in flags] == ["31.402,80"]


def test_a_turn_opened_twice_does_not_reuse_the_previous_ledger(book):
    first = office.begin_figure_turn("office-session")
    office.record_tool_figures("sum_column", office.sum_column(book, "Betrag"))
    assert office.figure_coverage_caveat(TOTAL_ANSWER, token=first) == ""

    second = office.begin_figure_turn("office-session")
    assert second != first
    assert office.figure_ledger(second) == []
    assert office.figure_coverage_caveat(TOTAL_ANSWER, token=second)


# ---------------------------------------------------------------------------
# Bounded, in both dimensions
# ---------------------------------------------------------------------------

def test_the_figures_of_one_scope_are_still_capped():
    _foreign_reads("token-A", reads=40)
    assert len(office.figure_ledger("token-A")) <= (
        office.MAX_LEDGER_FIGURES + 1)


def test_the_number_of_scopes_is_bounded_and_the_live_one_survives():
    """A long session must not accumulate a ledger per turn for ever, and
    the eviction must never reach the turn doing the evicting."""
    live = office.begin_figure_turn("office-session")
    office.record_tool_figures("sum_column", {"figures": [
        {"value": 4711.0, "kind": "sum", "label": "live"}]}, live)
    for i in range(office.MAX_LEDGER_SCOPES * 3):
        office.record_tool_figures("sum_column", {"figures": [
            {"value": float(i), "kind": "sum", "label": "old"}]}, f"old-{i}")
        office.figure_ledger(live)          # the live turn keeps working
    assert len(office._LEDGERS) <= office.MAX_LEDGER_SCOPES
    assert [f.value for f in office.figure_ledger(live)] == [4711.0]


# ---------------------------------------------------------------------------
# The engine hands its own token to the check
# ---------------------------------------------------------------------------

def test_the_engine_opens_a_turn_ledger_and_checks_against_that_token():
    import inspect

    from delfin.agent import engine as engine_mod

    src = inspect.getsource(engine_mod.AgentEngine)
    assert "begin_figure_turn(" in src
    assert "_figure_ledger_token" in src
    caveat = inspect.getsource(
        engine_mod.AgentEngine._append_figure_coverage_caveat)
    assert "token=" in caveat
