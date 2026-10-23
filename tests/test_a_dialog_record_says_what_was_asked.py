"""A dialog's audit record says what was asked, and why it was refused.

The first analysis of a round's dialogs (2026-09-26) found the records
carried tool and answer only: 4 of 12 could be matched to their command
by timestamps alone, 4 not at all.
"""

from __future__ import annotations

import threading

from delfin.agent import terminal_confirm as tc


def _answered(monkeypatch, args, decision, reason=""):
    records = []
    monkeypatch.setattr("delfin.agent.audit_log.append",
                        lambda rec, **kw: records.append(rec))
    broker = tc.TerminalConfirmBroker(session_id="s1", session_key="w-s1")
    req = tc.ConfirmRequest(kind="confirm", tool="bash", args=args,
                            preview="$ x")
    broker._enqueue(req)

    def answer():
        if reason:
            broker.last_refusal_reason = reason
        broker.resolve(req, decision)

    threading.Timer(0.05, answer).start()
    broker._wait(req)
    return [r for r in records if r.get("event") == "dialog"][0]


def test_the_record_names_the_command(monkeypatch):
    rec = _answered(monkeypatch, {"command": "awk  '{print $1}'\nfile"}, True)
    assert rec["subject"] == "awk '{print $1}' file"
    assert rec["reason"] == ""


def test_a_long_subject_is_clipped(monkeypatch):
    rec = _answered(monkeypatch, {"path": "/x/" + "a" * 500}, True)
    assert len(rec["subject"]) == 200
