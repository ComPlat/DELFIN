"""The task line shows the plan on every backend.

The Anthropic CLI backend keeps its plan in its own TodoWrite list, not in
DELFIN's task store, and the dashboard's task line read only the store: on
that backend the line never appeared. It now reads the TodoWrite list when
the engine has no KIT permissions. The plan also no longer outlives its
session: a new session clears it, and a resume takes the saved one even
when that is empty.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import task_ticker as TT

SRC = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
       / "tab_agent.py").read_text(encoding="utf-8")


def test_todo_rows_carry_order_label_and_subject():
    rows = TT.rows_from_todos([
        {"content": "Write the test", "status": "completed",
         "activeForm": "Writing the test"},
        {"content": "Fix the parser", "status": "in_progress",
         "activeForm": "Fixing the parser"},
        {"content": "Ship it", "status": "pending"},
        {"content": "", "status": "pending"},          # dropped
        "not a dict",                                  # dropped
    ])
    assert [r["status"] for r in rows] == [
        "in_progress", "pending", "completed"]
    assert rows[0]["label"] == "Fixing the parser"
    assert rows[0]["subject"] == "Fix the parser"
    assert TT.title_for(rows) == (
        "Tasks  ▶ 1  ☐ 1  ☑ 1 · Fixing the parser")


def test_a_finished_or_empty_plan_has_no_line():
    done = TT.rows_from_todos([{"content": "a", "status": "completed"}])
    assert TT.title_for(done) == ""
    assert TT.title_for(TT.rows_from_todos([])) == ""


def test_an_unknown_status_reads_as_pending():
    rows = TT.rows_from_todos([{"content": "a", "status": "weird"}])
    assert rows[0]["status"] == "pending"


def test_the_dashboard_reads_the_todo_list_without_kit_permissions():
    i = SRC.index("def _refresh_task_ticker():")
    body = SRC[i:SRC.index("def _", i + 30)]
    assert "if eng is not None and kp is None:" in body
    assert 'rows_from_todos' in body and 'state.get("current_todos")' in body


def test_a_todowrite_refreshes_the_line():
    i = SRC.index('elif tool_name == "TodoWrite":')
    assert "_refresh_task_ticker()" in SRC[i:i + 600]


def test_the_plan_does_not_outlive_its_session():
    # Resume: the saved plan, empty or not.
    assert ('state["current_todos"] = list(data.get("todo_payload") or [])'
            in SRC)
    assert 'if todo_payload:\n            state["current_todos"]' not in SRC
    # New session: cleared next to the new session id.
    i = SRC.index('state["active_session_id"] = _new_sid')
    assert 'state["current_todos"] = []' in SRC[i:i + 200]
