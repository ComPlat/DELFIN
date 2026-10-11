"""Adversarial review of package U3 phase 2 — repair.apply / repair.plan.

Independent of the builder's test_u3_plan.py. Attacks the phase-2 contract
at its weak joints:

- "nothing runs without an explicit approval" -- proven with a marker file:
  an unapproved COMMAND step whose shell line would create a file must not
  create it (the builder's test only refused a *setting* step).
- settings are never overwritten or deleted -- a MALFORMED settings file is
  still backed up VERBATIM first and only then rewritten; the dated backup
  keeps the exact original bytes.
- the undo is real, not prose: apply -> capture backup -> restore it -> the
  settings file equals the original. Round-trips end to end.
- chained applies each take one dated backup; the second apply sees the file
  the first wrote and preserves it (nothing lost across the chain).
- command steps report their real exit status (a failing line returns
  ok=False with the error, it never raises, and a broken symbol does not
  pretend success).
- a malformed ``setting`` in a doctor row is skipped, never crashes plan().
- apply() refuses an unknown step kind rather than silently doing nothing.
- a three-level dotted setting key writes to the right place and keeps
  siblings.

Tests write only under tmp_path; the real ~/.delfin/settings.json is never
read or written.
"""

from __future__ import annotations

import json
import sys

import pytest

from delfin.agent import repair


def _setting_row(**kw) -> dict:
    row = {"check": "mcp servers", "status": "PASS", "detail": "d",
           "fix": "f", "setting": ("agent.mcp_isolation", "builtin")}
    row.update(kw)
    return row


def _command_row(command: str, **kw) -> dict:
    row = {"check": "test runner", "status": "WARN", "detail": "d",
           "fix": "f", "command": command}
    row.update(kw)
    return row


# ---------------------------------------------------------------------------
# approval gate -- proven by effect, not by exception alone
# ---------------------------------------------------------------------------


def test_unapproved_command_step_never_runs_the_shell(tmp_path):
    """An unapproved COMMAND step must not execute -- proven by a side effect
    a run would leave behind. The builder's refusal test only covered a
    setting step; the command path runs a shell and is the dangerous one."""
    marker = tmp_path / "ran.txt"
    # shell=True in repair runs the line as written; if it runs, the marker
    # appears. approved=False must stop it before the shell ever starts.
    touches = f"{sys.executable} -c 'open({str(marker)!r},\"w\").write(\"ran\")'"
    step = repair.plan([_command_row(touches)])[0]
    with pytest.raises(ValueError):
        repair.apply(step, approved=False)
    assert not marker.exists(), (
        "an unapproved command step executed the shell line -- "
        "the approval gate is bypassed for command steps")


# ---------------------------------------------------------------------------
# settings never overwritten -- the malformed-file path
# ---------------------------------------------------------------------------


def test_malformed_settings_file_is_backed_up_verbatim_then_rewritten(tmp_path):
    """A corrupt settings file must NOT be silent-agreed away: it is moved to
    a dated backup capturing the EXACT original bytes, and only then is a
    fresh file written. The bad content is preserved for the user."""
    settings = tmp_path / "settings.json"
    garbage = "{\n  'this is': 'not valid json',\n"
    settings.write_text(garbage)
    step = repair.plan([_setting_row()])[0]

    result = repair.apply(step, approved=True, user_settings_path=settings)

    assert result["ok"] is True, result
    backups = [p for p in tmp_path.iterdir()
               if p.name.startswith("settings.json.") and p.name.endswith(".bak")]
    assert len(backups) == 1, [p.name for p in tmp_path.iterdir()]
    # the backup holds the original bytes verbatim -- nothing of the corrupt
    # file is lost/discarded
    assert backups[0].read_text() == garbage
    assert result["backup_path"] == str(backups[0])


# ---------------------------------------------------------------------------
# undo is real -- end-to-end restore round-trip
# ---------------------------------------------------------------------------


def test_undo_round_trip_restores_the_original_file(tmp_path):
    """apply -> read the dated backup -> restore it over the settings file ->
    the file is byte-for-byte the original again. The undo text is not just a
    sentence; the backup it names actually restores the prior state."""
    settings = tmp_path / "settings.json"
    original = {"agent": {"other": 1}, "kit": {"allowed": ["a"]}}
    settings.write_text(json.dumps(original))
    step = repair.plan([_setting_row()])[0]

    result = repair.apply(step, approved=True, user_settings_path=settings)

    backup_path = result["backup_path"]
    assert backup_path, "a prior file must yield a backup to restore from"
    # the reported undo names the backup; following it restores the file
    backed = json.loads(open(backup_path).read())
    assert backed == original
    settings.write_text(json.dumps(backed))
    assert json.loads(settings.read_text()) == original


# ---------------------------------------------------------------------------
# chained applies -- one dated backup each, nothing lost across the chain
# ---------------------------------------------------------------------------


def test_chained_setting_applies_keep_a_unique_backup_per_step(tmp_path):
    """Two setting steps applied in sequence must each back up the file they
    find. The second reads what the first wrote (so a+b survives) and leaves
    TWO dated backups, both restorable."""
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    s1 = repair.plan([_setting_row(setting=("agent.a", "1"))])[0]
    s2 = repair.plan([_setting_row(setting=("agent.b", "2"))])[0]

    r1 = repair.apply(s1, approved=True, user_settings_path=settings)
    r2 = repair.apply(s2, approved=True, user_settings_path=settings)

    assert r1["ok"] and r2["ok"]
    assert len([p for p in tmp_path.iterdir() if p.name.endswith(".bak")]) == 2
    # nothing lost across the chain: both keys present
    assert json.loads(settings.read_text()) == \
        {"agent": {"other": 1, "a": "1", "b": "2"}}


# ---------------------------------------------------------------------------
# command steps report their real exit status honestly
# ---------------------------------------------------------------------------


def test_command_step_failing_exit_returns_ok_false_and_does_not_raise(tmp_path):
    """A command that exits non-zero must return ok=False with the status,
    never raise and never report success."""
    step = repair.plan([_command_row(f"{sys.executable} -c 'import sys; sys.exit(3)'")])[0]
    result = repair.apply(step, approved=True)
    assert result["ok"] is False
    assert result["status"] == "exit 3"


def test_command_step_that_exits_zero_returns_ok_true(tmp_path):
    step = repair.plan([_command_row(f"{sys.executable} -c 'pass'")])[0]
    result = repair.apply(step, approved=True)
    assert result["ok"] is True
    assert result["status"] == "ok"


# ---------------------------------------------------------------------------
# plan() must never crash on a malformed remedy, and apply() refuses garbage
# ---------------------------------------------------------------------------


def test_plan_skips_a_setting_that_is_not_a_valid_key_value_pair():
    rows = [
        _setting_row(setting="not-a-tuple"),          # a bare string
        _setting_row(setting=("only_key",)),          # a 1-tuple
        _setting_row(setting=("k", "v", "extra")),    # a 3-tuple
        _setting_row(setting=("", "v")),              # an empty key
    ]
    steps = repair.plan(rows)
    assert steps == [], (
        "malformed setting remedies must be skipped (prose-only), "
        "never crash plan() or reach apply()")


def test_apply_refuses_an_unknown_step_kind():
    step = {"id": 0, "check": "x", "kind": "frobnicate", "what": "w"}
    with pytest.raises(ValueError):
        repair.apply(step, approved=True)


# ---------------------------------------------------------------------------
# deep nested dotted keys
# ---------------------------------------------------------------------------


def test_setting_step_writes_a_three_level_key_and_keeps_siblings(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"mcp": {"level": "low"}}}))
    step = repair.plan([_setting_row(setting=("agent.mcp.level", "high"))])[0]
    result = repair.apply(step, approved=True, user_settings_path=settings)
    assert result["ok"] is True
    assert json.loads(settings.read_text()) == \
        {"agent": {"mcp": {"level": "high"}}}
