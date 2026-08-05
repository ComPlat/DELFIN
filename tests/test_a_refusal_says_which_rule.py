"""A refusal was recorded without the reason it was refused for.

The gate writes a sentence naming the rule and what to do instead --
"command rejected by deny-pattern", the secret-path text, the egress
text, the not-auto-allowed text. That sentence went to the model and was
dropped on the way to the durable record, which has no reason field at
all. So `/changes` said:

    - bash: python analyse.py [denied]

and nothing else. After an unattended overnight run that stopped having
produced nothing, the user cannot tell whether it was the deny-list, the
secret scan, the egress scan, the workspace boundary, an auto-allow miss
or a hook -- and therefore cannot know what to allow. The one place the
rule name survived was a process-local deque that dies with the process.

Also here: a PreToolUse hook writes decision "block", which was in
neither the denied set nor any of the write/command/persist families, so
a hook that stopped an edit dropped out of the report entirely.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import audit_log as A


def _record(**kw):
    base = {"tool": "bash", "decision": "denied", "command": "rm -rf build"}
    base.update(kw)
    return A.make_record(**base)


# ---------------------------------------------------------------------------
# The record carries it
# ---------------------------------------------------------------------------

def test_the_reason_is_recorded():
    rec = _record(reason="command rejected by deny-pattern 'rm -rf'")
    assert "deny-pattern" in rec["reason"]


def test_a_record_without_one_is_unchanged():
    assert "reason" not in _record()


def test_the_reason_is_bounded():
    rec = _record(reason="x" * 5000)
    assert len(rec["reason"]) <= 300


# ---------------------------------------------------------------------------
# ...and the report shows it
# ---------------------------------------------------------------------------

def _report(tmp_path, records):
    log = tmp_path / "audit.log"
    for rec in records:
        A.append(rec, log_path=log)
    return A.build_changes_report(session_id="s1", log_path=log)


def test_the_changes_report_carries_the_reason(tmp_path):
    rep = _report(tmp_path, [_record(
        session_id="s1",
        reason="command rejected by deny-pattern 'rm -rf'")])
    assert rep["denied"], rep
    assert "deny-pattern" in rep["denied"][0]["reason"]


def test_the_rendered_line_names_the_rule(tmp_path):
    rep = _report(tmp_path, [_record(
        session_id="s1",
        reason="command rejected by deny-pattern 'rm -rf'")])
    text = A.format_changes_report(rep)
    assert "deny-pattern" in text, text


def test_a_denial_without_a_reason_still_renders(tmp_path):
    rep = _report(tmp_path, [_record(session_id="s1")])
    text = A.format_changes_report(rep)
    assert "bash" in text and "denied" in text


def test_a_hook_block_appears_in_the_report(tmp_path):
    """It matched no family before and vanished from the report."""
    rep = _report(tmp_path, [A.make_record(
        tool="hook", decision="block", session_id="s1",
        command="edit_file", reason="blocked by PreToolUse hook: protected path")])
    assert rep["denied"], rep
    assert "PreToolUse" in rep["denied"][0]["reason"]


# ---------------------------------------------------------------------------
# The gate's own message is what gets recorded
# ---------------------------------------------------------------------------

def test_the_executor_passes_the_gate_message_through(tmp_path):
    from delfin.agent import api_client as A2

    seen = {}

    def spy(record, **kw):
        seen.update(record)

    ex = A2._DocToolExecutor.__new__(A2._DocToolExecutor)
    perms = A2.KitToolPermissions(workspace=tmp_path)
    import delfin.agent.audit_log as _al
    orig = _al.append
    _al.append = spy
    try:
        ex._audit_call(
            "bash", {"command": "rm -rf /"}, perms,
            '{"error": "command rejected by deny-pattern \'rm -rf\'"}')
    finally:
        _al.append = orig

    assert seen.get("decision") == "denied"
    assert "deny-pattern" in seen.get("reason", "")


def test_a_successful_call_records_no_reason(tmp_path):
    from delfin.agent import api_client as A2

    seen = {}
    ex = A2._DocToolExecutor.__new__(A2._DocToolExecutor)
    perms = A2.KitToolPermissions(workspace=tmp_path)
    import delfin.agent.audit_log as _al
    orig = _al.append
    _al.append = lambda record, **kw: seen.update(record)
    try:
        ex._audit_call("bash", {"command": "ls"}, perms, "a.py\nb.py\n")
    finally:
        _al.append = orig

    assert seen.get("decision") == "ok"
    assert "reason" not in seen
