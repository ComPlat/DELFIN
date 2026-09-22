"""The report of a session is the session's own, and says what ran.

Two follow-ups to the session report, both found reading it rather than
running it:

  the file is owner-only     it collects every command of a session in
                             one place, while every store it reads from
                             is owner-only already
  the directory too          a shared home makes the difference visible
  a command, not its JSON    the trace stores a call's input as JSON, so
                             the first line of it was the blob, not the
                             command
  raw text still works       a recorder that stored the command itself
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import session_report


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(session_report.Path, "home", staticmethod(lambda: tmp_path))
    return tmp_path


def test_the_report_and_its_directory_are_owner_only(home, monkeypatch):
    monkeypatch.setattr(session_report, "collect_session_report",
                        lambda sid: session_report.SessionReport(session_id=sid))
    path = session_report.write_session_report("s1")
    assert path is not None
    assert os.stat(path).st_mode & 0o777 == 0o600
    assert os.stat(path.parent).st_mode & 0o777 == 0o700


def test_no_temporary_file_is_left_behind(home, monkeypatch):
    monkeypatch.setattr(session_report, "collect_session_report",
                        lambda sid: session_report.SessionReport(session_id=sid))
    path = session_report.write_session_report("s1")
    assert list(path.parent.glob("*.tmp")) == []


def test_the_commands_are_commands_and_not_json():
    entries = [
        {"tool": "mcp__kit-coding__bash",
         "input": json.dumps({"command": "git add a.py",
                              "description": "Stage the file"})},
        {"tool": "bash",
         "input": json.dumps({"command": "pytest -q\nsecond line",
                              "description": "run it"})},
    ]
    assert session_report._commands_run(entries) == ["git add a.py", "pytest -q"]


def test_a_recorder_that_stored_the_command_itself_still_works():
    entries = [{"tool": "bash", "input": "ls -la"},
               {"tool": "bash", "input": {"command": "whoami"}}]
    assert session_report._commands_run(entries) == ["ls -la", "whoami"]


def test_a_call_without_a_command_is_skipped_quietly():
    entries = [{"tool": "bash", "input": None},
               {"tool": "read_file", "input": json.dumps({"path": "x.py"})},
               {"tool": "bash", "input": json.dumps({"description": "no command"})}]
    assert session_report._commands_run(entries) == []


def test_the_rendered_report_shows_the_command(home, monkeypatch):
    report = session_report.SessionReport(
        session_id="s1",
        commands_run=session_report._commands_run(
            [{"tool": "bash", "input": json.dumps({"command": "git status"})}]),
    )
    text = session_report.render_markdown(report)
    assert "git status" in text
    assert "description" not in text
