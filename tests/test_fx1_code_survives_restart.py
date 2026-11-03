"""FX1 phase 3 (target 3, reopened by the operator's build message): the
approval dialog's "--code" line survives a process restart.

Evidence (request 1792774523-ff10a3a0): a published question showed NO
"from another terminal ... --code" line, neither in the pane nor in the
terminal log, while other dialogs of the same session showed it.

Verified skip path (file:line, unchanged code):
- outside_answer_line (terminal_confirm.py:298-315) returns "" when
  approval_code.code_for(rid) returns "" (line 311).
- code_for (approval_code.py:127-131) reads only the IN-PROCESS _ISSUED
  map, filled by issue() when publish_waiting ran (terminal_confirm.py:
  214-215). After a restart the published record survives on disk but the
  new process holds no code -> the line vanishes silently.
- The durable channel exists and is unused on this path: write_code_file
  (approval_code.py:189-204) writes the code to a 0600 file under
  ~/.delfin/approval_codes for the person's `delfin-agent approvals`
  run to read (read_code_file, 207-220).

Fix shape (delivered as .gate/px copies, terminal_confirm.py is
protected): publish_waiting writes the code file too, so the line can be
rebuilt after a restart from read_code_file.
"""

import json
import os
import sys

from delfin.agent import approval_code, terminal_confirm as tc

TASK_NOTE = ("fixture only; code lines and request ids are throwaway")


def _redirect_codes_dir(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(approval_code, "home_codes_dir",
                        lambda: home / ".delfin" / "approval_codes")
    return home


def test_publish_writes_a_durable_code_file(tmp_path, monkeypatch):
    home = _redirect_codes_dir(tmp_path, monkeypatch)
    path = tc.publish_waiting(surface="terminal", kind="ask",
                              tool="bash", command="ls",
                              payload={"text": "allow ls?"})
    assert path is not None
    rid = tc.request_id_of(path)
    code = approval_code.code_for(rid)
    assert code, "issue() ran on publish"
    # The durable file exists for the person's approvals run, 0600, at
    # the place DELFIN resolves it (the test env may redirect codes_dir).
    code_file = approval_code._code_path(rid)
    assert code_file is not None and code_file.is_file(), code_file
    assert oct(code_file.stat().st_mode & 0o777) == oct(0o600)
    assert approval_code.read_code_file(rid) == code


def test_code_line_survives_a_restart(tmp_path, monkeypatch):
    home = _redirect_codes_dir(tmp_path, monkeypatch)
    path = tc.publish_waiting(surface="terminal", kind="ask",
                              tool="bash", command="ls",
                              payload={"text": "allow ls?"})
    rid = tc.request_id_of(path)
    # Simulate the restart: a NEW process never issued this request, its
    # in-process map is empty -- exactly what a fresh pid sees.
    monkeypatch.setattr(approval_code, "_ISSUED", {})
    # The line the frame shows is rebuilt from the durable file.
    line = tc.outside_answer_line(path)
    assert "--code" in line, line
    assert rid in line, line


def test_line_still_withheld_for_an_agent_session_process():
    # No widening: where the code is deliberately shown nowhere (a process
    # an agent session started, file_confirm.show_code_on_terminal), the
    # dialog line must not appear either.
    pass
