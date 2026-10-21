"""A terminal ask record names its tool and keeps its question and options.

``approvals ls/show/answer`` read the published record; without the
payload they showed a blank body for a question. Only a question's
payload is published -- other kinds may carry file contents.
"""

from __future__ import annotations

import json

from delfin.agent import terminal_confirm as tc

_QUESTION = {"question": "Restore how?",
             "options": [{"label": "Allow restore now"},
                         {"label": "I'll restore manually"}]}


def test_ask_user_names_the_tool(monkeypatch, tmp_path):
    room = tmp_path / "terminal_confirmations"
    room.mkdir()
    monkeypatch.setattr(tc, "_PENDING_DIR", room)
    payload = dict(_QUESTION)
    # Driving ask_user end to end would block on _wait; publish one side.
    req = tc.ConfirmRequest(kind=tc.ASK, tool="ask_user_question",
                            payload=payload)
    path = tc._publish_pending(req, "sess-x", "pane-x")
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record.get("tool") == "ask_user_question"
    assert record.get("payload") == payload


def test_a_non_question_publishes_no_payload(monkeypatch, tmp_path):
    room = tmp_path / "terminal_confirmations"
    room.mkdir()
    monkeypatch.setattr(tc, "_PENDING_DIR", room)
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file",
                            payload={"content": "secret file body"})
    path = tc._publish_pending(req, "sess-x", "pane-x")
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record.get("payload") == {}
