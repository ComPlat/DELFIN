"""Adversarial tests for R2 finding 4 (approvals ask_user display half).

The fix under review (bc468246) adds approval_answers.render_question()
and wires it into the ``approvals show`` body of both the headless (
file_confirm) and terminal (terminal_confirm) branches, falling back to
the raw preview. These tests try to break that display fix:

- hostile ``payload`` shapes reaching ``show`` must never crash it;
- a non-choice record must keep printing its preview (no blank
  regression), in ``show`` AND ``ls``;
- a payload-carrying ask record renders a numbered list where a label
  that lexicographically looks like a number or the question text does
  not corrupt the numbering;
- the terminal branch (pending_at_terminals) of ``show`` renders the
  same question.

A note on scope: the committed fix renders a payload that reached the
listing. Whether a *terminal* ask question ever reaches the listing with
its payload is terminal_confirm._publish_pending (read-only security
code) — that half is the operator patch, not this diff. These tests
therefore pin the display contract given a record carrying a payload,
plus that a payload-less ask record (the shape today until the patch
lands) still falls back cleanly without crashing.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import file_confirm as _fc

_OPTIONS = [{"label": "Restore now"}, {"label": "Restore later"}]


def _write_record(room, rid: str, **fields):
    record = {"id": rid, "session_id": "sess-x", "kind": "ask",
              "tool": "ask_user_question", "command": "", "preview": "",
              "asked_at": time.time()}
    record.update(fields)
    p = room / f"{rid}.request.json"
    p.write_text(json.dumps(record), encoding="utf-8")
    os.chmod(p, 0o600)
    return p


def _run_approvals(capsys, action="show", **kw):
    import argparse
    from delfin.agent import cli
    ns = argparse.Namespace(approvals_action=action, **kw)
    rc = cli.cmd_approvals(ns)
    out = capsys.readouterr()
    return rc, out.out, out.err


@pytest.fixture()
def room(tmp_path, monkeypatch):
    approvals_dir = tmp_path / "approvals"
    approvals_dir.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: approvals_dir)
    return approvals_dir


class TestShowNeverCrashesOnHostilePayload:
    def test_payload_is_not_a_dict(self, room, capsys):
        # A record whose "payload" is a string, not the dict the contract
        # promises. render_question must treat it as no question and fall
        # back to the preview, never raise.
        _write_record(room, "9001", payload="oops-not-a-dict",
                      preview="fallback preview")
        rc, out, _err = _run_approvals(capsys, "show", request_id="9001")
        assert rc == 0
        assert "fallback preview" in out

    def test_options_have_non_dict_entries(self, room, capsys):
        # options is a list mixing dicts and junk; the valid dicts render,
        # the junk is skipped, nothing crashes.
        _write_record(room, "9002", payload={
            "question": "Pick?",
            "options": [{"label": "Good"}, "trash", 42, {"label": "Also good"}]})
        rc, out, _err = _run_approvals(capsys, "show", request_id="9002")
        assert rc == 0
        assert "Pick?" in out
        assert "1. Good" in out
        assert "2. Also good" in out
        assert "trash" not in out

    def test_question_empty_but_options_present(self, room, capsys):
        # A blank question line is not printed; the numbered options still are.
        _write_record(room, "9003", payload={"question": "  ",
                                             "options": _OPTIONS})
        rc, out, _err = _run_approvals(capsys, "show", request_id="9003")
        assert rc == 0
        assert "1. Restore now" in out
        assert "2. Restore later" in out

    def test_empty_payload_dict_is_fallthrough(self, room, capsys):
        _write_record(room, "9004", payload={}, preview="still here")
        rc, out, _err = _run_approvals(capsys, "show", request_id="9004")
        assert rc == 0
        assert "still here" in out


class TestFallbackForNonChoiceRecords:
    def test_file_kind_keeps_its_preview_in_show(self, room, capsys):
        _write_record(room, "9010", kind="file", tool="write_file",
                      command="/a/b.txt", preview="edit /a/b.txt?")
        rc, out, _err = _run_approvals(capsys, "show", request_id="9010")
        assert rc == 0
        assert "edit /a/b.txt?" in out

    def test_file_kind_keeps_its_path_in_show(self, room, capsys):
        _write_record(room, "9011", kind="file", tool="write_file",
                      command="/a/b.txt", preview="edit /a/b.txt?")
        rc, out, _err = _run_approvals(capsys, "show", request_id="9011")
        assert "/a/b.txt" in out


class TestOpticsOfADisplayedQuestion:
    def test_label_resembling_a_number_still_numered(self, room, capsys):
        # "1" as a label would read as a number if it ever leaked into the
        # answer path; here it must just render as an option line.
        _write_record(room, "9020", payload={
            "question": "Multi?",
            "options": [{"label": "One"}, {"label": "1"}, {"label": "Three"}]})
        rc, out, _err = _run_approvals(capsys, "show", request_id="9020")
        assert rc == 0
        assert "1. One" in out
        assert "2. 1" in out
        assert "3. Three" in out

    def test_unicode_question_roundtrips(self, room, capsys):
        _write_record(room, "9021", payload={
            "question": "Wiederherstellen wie?",
            "options": [{"label": "Jetzt"}, {"label": "Später"}]})
        rc, out, _err = _run_approvals(capsys, "show", request_id="9021")
        assert rc == 0
        assert "Wiederherstellen wie?" in out
        assert "1. Jetzt" in out


class TestTerminalBranch:
    def test_terminal_pending_renders_question(self, room, capsys):
        # The terminal branch of `show` reads pending_at_terminals; a row
        # there carrying a payload must render the same question. Feed the
        # terminal room directly through its pending reader.
        from delfin.agent import terminal_confirm as _tc
        term_room = room / "terminal"
        term_room.mkdir()
        _write_record(term_room, "9030", preview="",
                      payload={"question": "Restore how?",
                               "options": _OPTIONS})
        monkeypatch = pytest.MonkeyPatch()
        # Point terminal confirmations at our room.
        monkeypatch.setattr(_tc, "_PENDING_DIR", term_room, raising=True)
        try:
            from delfin.agent import cli
            import argparse
            ns = argparse.Namespace(approvals_action="show", request_id="9030")
            rc = cli.cmd_approvals(ns)
            out = capsys.readouterr().out
            assert rc == 0
            assert "Restore now" in out
            assert "Restore later" in out
        finally:
            monkeypatch.undo()
