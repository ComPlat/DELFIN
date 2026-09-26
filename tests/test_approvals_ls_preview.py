"""Controls for the `approvals ls` preview line (order LF, phase 3).

`ls` showed id/tool/path and wait seconds -- never what the request
was ABOUT, so every question needed a second `show`. Wave 5 gap 2.

One line per entry now carries the subject (command, path, or the
question of a choice request), clipped to the terminal width.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import file_confirm as _fc


def _record(room, rid: str, **over):
    rec = {"id": rid, "session_id": "sess-x", "session_key": "pane-x",
           "kind": "confirm", "tool": "bash", "command": "", "path": "",
           "preview": "", "asked_at": time.time(), "pid": os.getpid(),
           "host": "here"}
    rec.update(over)
    p = room / f"{rid}.request.json"
    p.write_text(json.dumps(rec), encoding="utf-8")
    os.chmod(p, 0o600)
    return p


@pytest.fixture()
def headless_room(tmp_path, monkeypatch):
    room = tmp_path / "approvals"
    room.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: room)
    return room


def _ls(capsys):
    import argparse
    from delfin.agent import cli
    rc = cli.cmd_approvals(argparse.Namespace(approvals_action="ls"))
    out = capsys.readouterr()
    return rc, out.out


class TestLsPreview:
    def test_ls_shows_the_command(self, headless_room, capsys):
        _record(headless_room, "3001-hh", command="git push origin main")
        rc, out = _ls(capsys)
        assert "git push origin main" in out

    def test_ls_shows_the_path_when_no_command(self, headless_room, capsys):
        _record(headless_room, "3002-ii", tool="write_file",
                path="delfin/agent/cli.py")
        rc, out = _ls(capsys)
        assert "delfin/agent/cli.py" in out

    def test_ls_shows_the_question_of_a_choice_request(self, headless_room,
                                                       capsys):
        _record(headless_room, "3003-jj", kind="ask",
                payload={"question": "Restore now or manually?",
                         "options": [{"label": "now"}, {"label": "manually"}]})
        rc, out = _ls(capsys)
        assert "Restore now or manually?" in out

    def test_ls_clips_a_long_command_to_the_width(self, headless_room,
                                                  capsys, monkeypatch):
        _record(headless_room, "3004-kk",
                command="echo " + "x" * 300)
        # A fixed narrow terminal: the line must not exceed it.
        monkeypatch.setattr("shutil.get_terminal_size",
                            lambda: os.terminal_size((60, 24)))
        rc, out = _ls(capsys)
        assert "3004-kk" in out
        line = [l for l in out.splitlines() if "3004-kk" in l][0]
        assert len(line) <= 60, f"line is {len(line)} cols"
        # The ellipsis marks the clipping; the wait seconds stay visible
        # at the end -- they are the one number that keeps ticking.
        assert "…" in line
        assert line.rstrip().endswith("0s")

    def test_ls_width_falls_back_without_a_terminal(self, headless_room,
                                                    capsys, monkeypatch):
        _record(headless_room, "3005-ll",
                command="echo " + "y" * 300)
        monkeypatch.setattr("shutil.get_terminal_size",
                            lambda: (_f := None) or (_ for _ in ()).throw(
                                OSError("no terminal")))
        rc, out = _ls(capsys)
        assert rc == 0
        line = [l for l in out.splitlines() if "3005-ll" in l][0]
        assert len(line) <= 120

    def test_ls_still_shows_waiting_at_terminals_whole(
            self, tmp_path, monkeypatch, capsys):
        from delfin.agent import terminal_confirm as _tc
        room = tmp_path / "terminal_confirmations"
        room.mkdir()
        monkeypatch.setattr(_tc, "_PENDING_DIR", room)
        monkeypatch.setattr(_fc, "requests_dir",
                            lambda: tmp_path / "approvals")
        (tmp_path / "approvals").mkdir()
        _record(room, "3006-mm", command="ls -la")
        rc, out = _ls(capsys)
        assert "ls -la" in out
        assert "Waiting at a terminal" in out
