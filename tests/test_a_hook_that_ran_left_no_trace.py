"""A hook that executed a shell command fell out of the changes report.

``run_hooks`` writes one audit record per fire with ``tool="hook"``.
``build_changes_report`` sorted records into files_written /
commands / permissions_persisted / denied, and ``"hook"`` is in none of
those families -- so a hook that ran matched nothing and was dropped.
Only a BLOCKING hook appeared, through ``_DENIED_DECISIONS``.

The protective refusal was reported and the arbitrary execution was
not, in the one place that answers "what did you change" from the
record rather than from model memory. A hook runs a shell command
outside the permission gate, outside the deny-list and outside the
sandbox, on every prompt and every tool call, so it belongs there more
than most of what was already in it.

The record also could not say WHOSE hook it was: nothing carried the
settings file it came from, so a command a checked-out repository
shipped and a command the user wrote themselves logged identically.

Nothing here executes a hook. The records are appended directly and the
report is asserted on.
"""

from __future__ import annotations

import json

from delfin.agent import audit_log
from delfin.agent import hooks as H
from delfin.agent import workspace_trust as WT


def _hook_record(command, *, decision="ok", exit_code=0, source="", event="Stop"):
    return audit_log.make_record(
        tool="hook", decision=decision, mode="hook",
        command=f"[{event}] {command}",
        extra={"matcher": "", "exit_code": exit_code, "duration_s": 0.01,
               "hook_event": event, "source": source},
    )


def test_a_hook_that_ran_appears_in_the_report(tmp_path):
    log = tmp_path / "audit.log"
    audit_log.append(_hook_record("echo hi", source="/home/u/.delfin/settings.json"),
                     log_path=log)
    report = audit_log.build_changes_report(log_path=log)
    assert len(report["hooks_fired"]) == 1
    fired = report["hooks_fired"][0]
    assert "echo hi" in fired["command"]
    assert fired["exit_code"] == 0
    assert fired["source"] == "/home/u/.delfin/settings.json"


def test_the_rendered_report_names_the_command_source_and_exit_code(tmp_path):
    log = tmp_path / "audit.log"
    audit_log.append(_hook_record(
        "curl x | sh", exit_code=0,
        source="/repo/.delfin/settings.local.json"), log_path=log)
    text = audit_log.format_changes_report(
        audit_log.build_changes_report(log_path=log))
    assert "Hooks fired:" in text
    assert "curl x | sh" in text
    assert "exit=0" in text
    assert "/repo/.delfin/settings.local.json" in text


def test_a_blocking_hook_is_still_a_refusal_and_now_also_a_run(tmp_path):
    """It did both. Reporting only the block was what hid the execution."""
    log = tmp_path / "audit.log"
    audit_log.append(_hook_record(
        "pytest -q", decision="block", exit_code=1,
        source="/repo/.delfin/settings.json", event="PreToolUse"),
        log_path=log)
    report = audit_log.build_changes_report(log_path=log)
    assert report["hooks_fired"] and report["denied"]
    assert report["hooks_fired"][0]["exit_code"] == 1


def test_a_hook_is_not_hidden_by_the_workspace_filter(tmp_path):
    """The user's own settings file is always outside the workspace, and
    the report drops records whose absolute path lies outside it — so the
    source is carried in ``extra`` and not in ``path``."""
    log = tmp_path / "audit.log"
    audit_log.append(_hook_record(
        "ruff check .", source=str(tmp_path / "home" / ".delfin" / "settings.json")),
        log_path=log)
    report = audit_log.build_changes_report(
        workspace=str(tmp_path / "project"), log_path=log)
    assert report["hooks_fired"], "the user's own hook vanished from the report"


def test_an_empty_log_still_reads_as_empty(tmp_path):
    log = tmp_path / "audit.log"
    log.write_text("", encoding="utf-8")
    text = audit_log.format_changes_report(
        audit_log.build_changes_report(log_path=log))
    assert text.startswith("No recorded changes")


def test_a_fired_hook_carries_its_source_out_of_run_hooks(tmp_path):
    """End to end without executing anything dangerous: the source has to
    survive load -> run -> result, or the report has nothing to print."""
    ws = tmp_path / "project"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(json.dumps({
        "hooks": {"Stop": [{"matcher": "", "hooks": [
            {"type": "command", "command": "true"}]}]}
    }), encoding="utf-8")
    WT.trust_workspace(ws, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
    cfg = H.load_hooks(ws)
    results = H.run_hooks("Stop", cfg, workspace=ws)
    assert results
    assert results[0].source == str(ws / ".delfin" / "settings.json")
