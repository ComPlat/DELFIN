"""A count the agent's own script printed is not flagged as invented.

Measured on a Raspberry Pi install: asked to split a 737-row workbook into
one sheet per collaborator, the agent wrote and ran an openpyxl script that
printed the sheet and row counts, then reported them. Every answer carried
"This figure did not come from a tool result of this turn" for 192 and 49,
the numbers the script had counted. A command that succeeded is a call of
this turn, and its printed summary now grounds the counts in it; a failed
command and a figure nothing printed are still flagged.
"""

from __future__ import annotations

import pytest

from delfin.agent import office
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def perms(tmp_path, monkeypatch):
    """A turn that has read a workbook, as the measured one had."""
    monkeypatch.setattr(office, "_LEDGERS", {}, raising=False)
    book = tmp_path / "platten.xlsx"
    wb = openpyxl.Workbook()
    wb.active.append(["Platte", "Collaborator"])
    wb.active.append(["P1", "Anna"])
    wb.save(book)
    p = KitToolPermissions(workspace=str(tmp_path))
    p.mode = "default"
    p.task_session_id = "script-figures"
    _DocToolExecutor().execute("read_document", {"path": str(book)}, p)
    return p


ANSWER = ("Fertig: Ich habe 192 Tabellenblätter angelegt. "
          "49 Zeilen ohne Collaborator stehen auf einem eigenen Blatt.")


def _flags():
    return office.scan_answer_for_unledgered_figures(ANSWER)


def test_counts_a_successful_script_printed_are_grounded(perms):
    _DocToolExecutor().execute(
        "bash", {"command": "printf 'sheets 192\\nempty 49\\n'",
                 "description": "summary"}, perms)
    assert [f.figure for f in _flags()] == []


def test_a_count_nothing_printed_is_still_flagged(perms):
    _DocToolExecutor().execute(
        "bash", {"command": "printf 'sheets 192\\n'",
                 "description": "summary"}, perms)
    assert [f.figure for f in _flags()] == ["49"]


def test_a_failed_command_grounds_nothing(perms):
    _DocToolExecutor().execute(
        "bash", {"command": "printf 'sheets 192\\nempty 49\\n'; exit 3",
                 "description": "summary"}, perms)
    assert sorted(f.figure for f in _flags()) == ["192", "49"]
