"""Package U3 — repair.plan / repair.apply over a doctor report.

The doctor report (delfin.agent.doctor.run_doctor) already marks the few
prerequisites that DELFIN can fix itself with a machine-actionable remedy:
a ``command`` row (test runner: pip install the test extra) and a
``setting`` row (MCP containment: ``agent.mcp_isolation = "builtin"``).
repair turns those rows into ordered steps, each carrying what it changes
and how to undo it, and applies one step at a time -- never without an
explicit approval, and never by overwriting a settings file (it is moved
to a dated backup first).

Tests use an explicit ``user_settings_path`` into tmp so the real
``~/.delfin/settings.json`` is never read or written.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import repair


def _report_row(**kw) -> dict:
    row = {"check": "test runner", "status": "WARN",
           "detail": "pytest is not installed",
           "fix": "install the test extra"}
    row.update(kw)
    return row


def _command_row(**kw) -> dict:
    row = _report_row(check="test runner",
                      command=f"{__import__('sys').executable} -m pip install 'x'")
    row.update(kw)
    return row


def _setting_row(**kw) -> dict:
    row = _report_row(check="mcp servers", status="PASS",
                      detail="2 without declared roots",
                      setting=("agent.mcp_isolation", "builtin"))
    row.update(kw)
    return row


# ---------------------------------------------------------------------------
# plan() -- steps ordered, each with what it changes and how to undo it
# ---------------------------------------------------------------------------


def test_plan_ignores_rows_without_a_machine_actionable_remedy():
    """Prose-only rows (system package, login a helper) have no declared
    command/setting -- the doctor says so on purpose, and repair must not
    improvise a remedy out of prose."""
    plain = _report_row(check="credentials", status="FAIL",
                        detail="no keys", fix="set a provider key")
    steps = repair.plan([plain, _command_row(), _setting_row()])
    checks = [s["check"] for s in steps]
    assert checks == ["test runner", "mcp servers"], steps


def test_plan_a_command_row_names_what_and_undo():
    steps = repair.plan([_command_row()])
    assert len(steps) == 1
    step = steps[0]
    assert step["kind"] == "command"
    assert step["command"]  # the pip install line is carried
    assert step["what"]     # a human sentence, not empty
    assert step["undo"]     # a stated way back, not empty


def test_plan_a_setting_row_names_key_value_what_and_undo():
    steps = repair.plan([_setting_row()])
    assert len(steps) == 1
    step = steps[0]
    assert step["kind"] == "setting"
    assert step["setting_key"] == "agent.mcp_isolation"
    assert step["setting_value"] == "builtin"
    assert step["what"]
    assert step["undo"]


def test_plan_steps_are_ordered_and_uniquely_identified():
    steps = repair.plan([_command_row(), _setting_row(),
                         _command_row(check="test runner 2")])
    ids = [s["id"] for s in steps]
    assert ids == sorted(ids), "steps must be stable and uniquely identified"
    assert len(set(ids)) == len(ids)


# ---------------------------------------------------------------------------
# apply() -- nothing runs without an approval
# ---------------------------------------------------------------------------


def test_apply_refuses_without_approval():
    step = repair.plan([_setting_row()])[0]
    with pytest.raises(ValueError):
        repair.apply(step, approved=False)


def test_apply_setting_moves_the_old_file_to_a_dated_backup_then_writes(tmp_path):
    """The core invariant: a settings file is never overwritten or deleted.
    The old file is moved (renamed) to a dated backup first, and the new
    file is written fresh from the backed-up content plus the setting."""
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    step = repair.plan([_setting_row()])[0]

    result = repair.apply(step, approved=True,
                          user_settings_path=settings)

    assert result["ok"] is True, result
    # old file still exists, under a dated name
    backups = [p for p in tmp_path.iterdir()
               if p.name.startswith("settings.json.") and p.name != settings.name]
    assert len(backups) == 1, [p.name for p in tmp_path.iterdir()]
    backed = json.loads(backups[0].read_text())
    assert backed == {"agent": {"other": 1}}, "backup keeps the old content"
    # new file is old content plus the setting (nothing else lost or added)
    assert json.loads(settings.read_text()) == \
        {"agent": {"other": 1, "mcp_isolation": "builtin"}}


def test_apply_setting_records_the_backup_so_undo_restores_it(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    step = repair.plan([_setting_row()])[0]
    result = repair.apply(step, approved=True, user_settings_path=settings)
    assert result["backup_path"]
    assert result["backup_path"] != str(settings)
    assert result["backup_path"].endswith(".bak")
    assert "restore" in result["undo"].lower(), result["undo"]


def test_apply_setting_with_no_prior_file_writes_fresh(tmp_path):
    settings = tmp_path / "settings.json"
    step = repair.plan([_setting_row()])[0]
    result = repair.apply(step, approved=True, user_settings_path=settings)
    assert result["ok"] is True, result
    assert json.loads(settings.read_text()) == {"agent": {"mcp_isolation": "builtin"}}
    # with no prior file there is nothing to back up
    assert not any(p.name.startswith("settings.json.") and p != settings
                   for p in tmp_path.iterdir())
