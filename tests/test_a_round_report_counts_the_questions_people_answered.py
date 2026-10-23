"""The round report counts the questions a person had to answer.

First real run of `delfin-agent report --since` (night of 2026-09-25):
every session read "dialogues: 0" although the operator had answered
dozens of permission dialogs, and `--name nacht-` matched nothing -- the
session name was stored nowhere the report could read. The broker now
writes one audit record per answered dialog (with the -n name), and the
report labels a session by that name or, failing it, its workspace.
"""

from __future__ import annotations

import json
import threading

from delfin.agent import round_report
from delfin.agent import terminal_confirm as tc


def test_dialog_records_count_and_label_the_session(tmp_path):
    audit = tmp_path / "audit.log"
    audit.write_text("\n".join(json.dumps(r) for r in (
        {"ts": "x", "session_id": "abc123", "tool": "bash", "decision": "ok",
         "workspace": "/w/trees/runde3-s2"},
        {"ts": "x", "event": "dialog", "session_id": "abc123",
         "session_key": "nacht-s2", "tool": "bash", "answer": "approved",
         "by": "outside"},
        {"ts": "x", "event": "dialog", "session_id": "abc123",
         "session_key": "nacht-s2", "tool": "bash", "answer": "denied",
         "by": "outside"},
    )) + "\n")
    labels, dialogs = round_report._audit_facts(audit)
    assert labels["abc123"] == "nacht-s2"
    assert dialogs["abc123"] == 2


def test_a_session_without_dialogs_is_labelled_by_its_workspace(tmp_path):
    audit = tmp_path / "audit.log"
    audit.write_text(json.dumps(
        {"ts": "x", "session_id": "def456", "tool": "read_file",
         "decision": "ok", "workspace": "/w/trees/runde3-s5"}) + "\n")
    labels, dialogs = round_report._audit_facts(audit)
    assert labels["def456"] == "runde3-s5"
    assert dialogs == {}


def test_the_broker_audits_an_answered_dialog(monkeypatch):
    records = []
    monkeypatch.setattr("delfin.agent.audit_log.append",
                        lambda rec, **kw: records.append(rec))
    broker = tc.TerminalConfirmBroker(session_id="s1", session_key="nacht-s1")
    req = tc.ConfirmRequest(kind="confirm", tool="bash", args={"command": "ls"},
                            preview="$ ls")
    broker._enqueue(req)
    threading.Timer(0.05, lambda: broker.resolve(req, True)).start()
    broker._wait(req)
    dialog = [r for r in records if r.get("event") == "dialog"]
    assert dialog and dialog[0]["session_key"] == "nacht-s1"
    assert dialog[0]["answer"] == "approved"
