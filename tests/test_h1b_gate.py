"""H1B phase 3: the bash gate carries the self-written-execution hint.

The gate calls self_written.verdict(command, _journal_changes(sid), workspace)
~api_client.py:17133 (a new tier between the egress block and the
already-refused tier). When a command EXECUTES a file the session wrote, the
refusal carries the hint; a SECOND spelling of the SAME file stops the turn
and asks the user (reusing denied_actions, key "self_written:<abs-path>"); a
script the session did NOT write, and any READ (cat) of a written file, are
judged as before (no hint).

Both go through the REAL executor (_doc_executor.execute): the session write
is a real write_file tool call (which records the change journal), and the
command is a real bash tool call. Path.home is redirected so no real user
state is touched.

RED before the gate patch: python3 x.py is currently ALLOWED (no refusal), so
test_self_written_hint_appears fails. GREEN after the tier -- the operator
builds it from .gate/px/new_api_client.py.
"""

import os
from pathlib import Path

import pytest

from delfin.agent.api_client import _doc_executor, KitToolPermissions


@pytest.fixture
def env(tmp_path, monkeypatch):
    # A project sub-dir SIBLING of the redirected home (tmp_path/home): using
    # tmp_path itself would be an ANCESTOR of the redirected home and
    # _is_forbidden_workspace_root (api_client.py:5625) refuses that.
    ws = tmp_path / "ws"
    ws.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    ws_str = str(ws)

    def perms(mode="acceptEdits"):
        return KitToolPermissions(workspace=ws_str, mode=mode,
                                  task_session_id="sid-h1b-gate")

    def write_file(path, content):
        p = ws / path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
        return _doc_executor.execute(
            "write_file", {"path": str(p), "content": content}, perms())

    def bash(cmd, mode="acceptEdits"):
        return _doc_executor.execute("bash", {"command": cmd}, perms(mode))

    return {"ws": ws_str, "write_file": write_file, "bash": bash,
            "perm": perms}


def test_self_written_hint_appears(env):
    env["write_file"]("x.py", "print(1)")
    out = env["bash"]("python3 x.py")
    assert "this session wrote" in out          # the hint in the refusal


def test_second_spelling_stops_turn(env):
    env["write_file"]("x.py", "print(1)")
    first = env["bash"]("python3 x.py")
    assert first and "this session wrote" in first
    second = env["bash"]("bash x.py")           # re-spelling
    assert "already tried to run" in second     # stop the turn, ask the user


def test_unwritten_script_judged_before(env):
    # x.py exists but was NOT written by the session (not in the journal).
    (env["ws"] / "x.py").write_text("print(1)")
    out = env["bash"]("python3 x.py")
    assert "this session wrote" not in (out or "")


def test_read_of_written_file_is_not_execution(env):
    env["write_file"]("x.py", "print(1)")
    out = env["bash"]("cat x.py")
    assert "this session wrote" not in (out or "")
