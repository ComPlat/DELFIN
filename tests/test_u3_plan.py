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


def test_two_setting_applies_keep_a_unique_backup_per_step(tmp_path, monkeypatch):
    """The phase-2 invariant "settings are never overwritten or deleted"
    must hold across two applies in the same second: each apply must move
    its OWN dated backup, so the FIRST (the true original) survives with the
    exact prior bytes. A second-resolved backup name lets the second apply
    overwrite the first's backup and silently delete the original."""
    import datetime as _dt

    frozen = _dt.datetime(2026, 10, 6, 12, 0, 0)  # fixed clock, same second
    class _FrozenClock(_dt.datetime):
        @classmethod
        def now(cls, tz=None):
            return frozen

    monkeypatch.setattr(repair, "datetime", _FrozenClock)

    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))

    step = repair.plan([_setting_row()])[0]
    repair.apply(step, approved=True, user_settings_path=settings)
    repair.apply(step, approved=True, user_settings_path=settings)

    backups = sorted(p.name for p in tmp_path.iterdir()
                     if p.name.startswith("settings.json.") and p != settings)
    assert len(backups) == 2, (
        f"two applies in one second must produce two distinct backups, "
        f"got {backups}")
    # the true original must be preserved verbatim in one of the backups
    preserved = any(
        json.loads((tmp_path / name).read_text()) == {"agent": {"other": 1}}
        for name in backups)
    assert preserved, (
        "the true original settings must survive in a backup; "
        f"backups: {[(n, (tmp_path / n).read_text()) for n in backups]}")


# ---------------------------------------------------------------------------
# The injectable apply seam the dashboard /doctor repair (and its tests) use
# ---------------------------------------------------------------------------
# doctor_repair_plan / apply_repair_step / doctor_repair_apply are the pure
# core behind the dashboard path: plan the steps, apply one step under an
# explicit approval, and drive an id -> step lookup + approval gate in one
# call. `authorize` is the injectable broker answer (real broker callback on
# the dashboard; a stub here), so these tests pin the approval semantics
# without any UI thread.


def test_doctor_repair_plan_returns_the_same_ordered_fixable_steps():
    steps = repair.doctor_repair_plan([_command_row(), _setting_row()])
    assert steps == repair.plan([_command_row(), _setting_row()])
    assert len(steps) == 2
    assert [s["check"] for s in steps] == ["test runner", "mcp servers"]


def test_apply_repair_step_refuses_without_approval(tmp_path):
    step = repair.plan([_setting_row()])[0]
    settings = tmp_path / "settings.json"
    with pytest.raises(ValueError):
        repair.apply_repair_step(step, approved=False,
                                 user_settings_path=settings)
    assert not settings.exists(), (
        "a refused step must not write the settings file")


def test_apply_repair_step_applies_a_setting_after_a_backup(tmp_path):
    step = repair.plan([_setting_row()])[0]
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    res = repair.apply_repair_step(step, approved=True,
                                   user_settings_path=settings)
    assert res["ok"] is True
    assert res["backup_path"], "a setting step must move a dated backup first"
    backup = tmp_path / res["backup_path"]
    assert backup.exists()
    assert json.loads(backup.read_text()) == {"agent": {"other": 1}}, (
        "the dated backup must hold the true prior settings")
    assert json.loads(settings.read_text()) == {
        "agent": {"mcp_isolation": "builtin", "other": 1}}, (
        "the rewritten settings file must carry the new value and keep the rest")


def test_doctor_repair_apply_unknown_id_never_asks():
    called = []
    def authorize(tool_name, args, preview):
        called.append((tool_name, args))
        return True
    res = repair.doctor_repair_apply("nope", results=[_setting_row()],
                                     authorize=authorize)
    assert res["ok"] is False
    assert "unknown step" in res["message"]
    assert called == [], "an unknown id must never reach the authorizer"


def test_doctor_repair_apply_declined_writes_nothing(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    seen = {}
    def authorize(tool_name, args, preview):
        seen["tool"] = tool_name
        seen["args"] = args
        seen["preview"] = preview
        return False
    res = repair.doctor_repair_apply(0, results=[_setting_row()],
                                     authorize=authorize,
                                     user_settings_path=settings)
    assert res["applied"] is False
    assert res["ok"] is True
    assert seen["tool"] == "repair"
    assert seen["args"] == {"step": 0}
    assert "undo" in seen["preview"], "the preview must carry the undo"
    # nothing written, no backup -- decline is a no-op for the disk
    assert json.loads(settings.read_text()) == {"agent": {"other": 1}}
    backups = [p for p in tmp_path.iterdir()
               if p.name.startswith("settings.json.")]
    assert backups == [], "a declined step must not create a backup"


def test_doctor_repair_apply_approved_backs_up_then_writes(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    def authorize(tool_name, args, preview):
        return True
    res = repair.doctor_repair_apply(0, results=[_setting_row()],
                                     authorize=authorize,
                                     user_settings_path=settings)
    assert res["applied"] is True
    assert res["result"]["ok"] is True
    backup = tmp_path / res["backup_path"]
    assert backup.exists()
    assert json.loads(backup.read_text()) == {"agent": {"other": 1}}, (
        "the dated backup must hold the true prior before the rewrite")
    assert json.loads(settings.read_text()) == {
        "agent": {"mcp_isolation": "builtin", "other": 1}}, (
        "the new value must be written and the rest preserved")
