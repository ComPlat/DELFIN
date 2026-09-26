"""Controls for the CLI wiring of `approvals answer` (order LF, phase 2).

These go through cmd_approvals the way the real parser calls it, and
were red on the previous commit: `approvals answer` did not exist, so
argparse rejected the choice before anything ran. Kept as the
integration half of the phase.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import approval_answers as aa
from delfin.agent import file_confirm as _fc


def _ask_record(room, rid: str):
    record = {
        "id": rid, "session_id": "sess-x", "session_key": "pane-x",
        "kind": "ask", "tool": "ask_user_question", "command": "",
        "preview": "", "asked_at": time.time(), "pid": os.getpid(),
        "host": "here",
        "payload": {"question": "Restore how?",
                    "options": [{"label": "Allow restore now"},
                                {"label": "I'll restore manually"}],
                    "multiSelect": False},
    }
    p = room / f"{rid}.request.json"
    p.write_text(json.dumps(record), encoding="utf-8")
    os.chmod(p, 0o600)
    return p


@pytest.fixture()
def headless_room(tmp_path, monkeypatch):
    room = tmp_path / "approvals"
    room.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: room)
    return room


def _run_cli(capsys, request_id, choice, reason=""):
    """cmd_approvals with the namespace the real parser builds."""
    import argparse
    from delfin.agent import cli
    args = argparse.Namespace(approvals_action="answer",
                              request_id=request_id, choice=choice,
                              reason=reason)
    rc = cli.cmd_approvals(args)
    out = capsys.readouterr()
    return rc, out.out, out.err


class TestApprovalsAnswerCli:
    def test_answer_by_number_writes_the_pick(self, headless_room, capsys):
        _ask_record(headless_room, "2001-dd")
        rc, out, err = _run_cli(capsys, "2001-dd", "1")
        assert rc == 0
        assert "Allow restore now" in out
        raw = json.loads(
            (headless_room / "2001-dd.answer.json").read_text(encoding="utf-8"))
        assert raw["answers"] == ["Allow restore now"]

    def test_answer_by_label(self, headless_room, capsys):
        _ask_record(headless_room, "2002-ee")
        rc, out, err = _run_cli(capsys, "2002-ee", "I'll restore manually")
        assert rc == 0
        raw = json.loads(
            (headless_room / "2002-ee.answer.json").read_text(encoding="utf-8"))
        assert raw["answers"] == ["I'll restore manually"]

    def test_answer_bad_choice_exits_two_and_keeps_the_question_open(
            self, headless_room, capsys):
        _ask_record(headless_room, "2003-ff")
        rc, out, err = _run_cli(capsys, "2003-ff", "nope")
        assert rc == 2
        assert "not one of the options" in err
        assert not (headless_room / "2003-ff.answer.json").exists()

    def test_answer_non_choice_refused(self, headless_room, capsys):
        rec = {"id": "2004-gg", "session_id": "s", "kind": "confirm",
               "tool": "bash", "command": "ls", "preview": "",
               "asked_at": time.time()}
        (headless_room / "2004-gg.request.json").write_text(
            json.dumps(rec), encoding="utf-8")
        rc, out, err = _run_cli(capsys, "2004-gg", "1")
        assert rc == 2
        assert "not a choice question" in err

    def test_parser_accepts_answer_subcommand(self):
        """The public surface: `delfin-agent approvals answer id choice`
        parses. Red before the wiring existed."""
        from delfin.agent.cli import build_parser
        p = build_parser()
        ns = p.parse_args(["approvals", "answer", "1234-ab", "2"])
        assert ns.approvals_action == "answer"
        assert ns.request_id == "1234-ab"
        assert ns.choice == "2"
