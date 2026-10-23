"""R2 finding 4 controls (display half): `approvals ls/show` render an
ask-user question.

An ask_user question is stored on the record under ``payload`` (question
+ options). ``approvals ls`` renders it through
``approval_answers.preview_line``; but ``approvals show`` printed only
the raw ``preview`` string, which an ask record leaves empty -- so even
a record that carries a full ``payload`` showed a blank question. These
controls pin that the question and its options are displayed.

A second, inseparable half -- that ``terminal_confirm._publish_pending``
persists ``payload`` so a *terminal* ask question reaches the listing at
all -- is a read-only security file, so it goes to the operator as a
patch (see .gate/r2_ask_persist.patch); its acceptance is
``test_published_ask_record_keeps_the_payload`` there.

The ``show`` members are RED on the unchanged code: the question never
appeared.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import file_confirm as _fc

_QUESTION = {"question": "Restore how?",
             "options": [{"label": "Allow restore now"},
                         {"label": "I'll restore manually"}],
             "multiSelect": False}


@pytest.fixture()
def headless_room(tmp_path, monkeypatch):
    room = tmp_path / "approvals"
    room.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: room)
    return room


def _write_ask(room, rid: str, **over):
    record = {"id": rid, "session_id": "sess-x", "kind": "ask",
              "tool": "ask_user_question", "command": "", "preview": "",
              "asked_at": time.time(),
              "payload": dict(_QUESTION)}
    record.update(over)
    p = room / f"{rid}.request.json"
    p.write_text(json.dumps(record), encoding="utf-8")
    os.chmod(p, 0o600)
    return p


def _run_approvals(capsys, action="ls", **kw):
    import argparse
    from delfin.agent import cli
    ns = argparse.Namespace(approvals_action=action, **kw)
    rc = cli.cmd_approvals(ns)
    out = capsys.readouterr()
    return rc, out.out, out.err


class TestAskQuestionVisibleInListing:
    def test_show_renders_the_question(self, headless_room, capsys):
        _write_ask(headless_room, "4001-aa")
        rc, out, _err = _run_approvals(capsys, "show", request_id="4001-aa")
        assert rc == 0
        assert "Restore how?" in out

    def test_show_renders_the_options(self, headless_room, capsys):
        _write_ask(headless_room, "4002-bb")
        rc, out, _err = _run_approvals(capsys, "show", request_id="4002-bb")
        assert rc == 0
        assert "Allow restore now" in out
        assert "I'll restore manually" in out

    def test_ls_renders_the_question(self, headless_room, capsys):
        _write_ask(headless_room, "4003-cc")
        rc, out, _err = _run_approvals(capsys, "ls")
        assert rc == 0
        assert "Restore how?" in out

    def test_answer_still_resolves_by_number(self, headless_room, capsys):
        _write_ask(headless_room, "4004-dd")
        rc, out, _err = _run_approvals(capsys, "answer", request_id="4004-dd",
                                       choice="1")
        assert rc == 0
        assert "Allow restore now" in out

    def test_answer_still_resolves_by_label(self, headless_room, capsys):
        _write_ask(headless_room, "4005-ee")
        rc, out, _err = _run_approvals(
            capsys, "answer", request_id="4005-ee",
            choice="I'll restore manually")
        assert rc == 0
        assert "I'll restore manually" in out
