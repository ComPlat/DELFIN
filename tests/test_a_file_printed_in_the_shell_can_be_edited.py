"""A file printed through the shell counts as read for edit_file.

Report 20260915-112305: the agent read a test fixture with
``sed -n '122,175p'``, and edit_file refused it -- "call read_file before
editing" -- which cost two rounds and four minutes on GLM to print the
same lines again. The baseline guards against an edit landing on a file
that changed since it was looked at; a file whose text was printed has been
looked at. What proves nothing about the text still sets nothing.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "cfg.py").write_text("A = 1\nB = 2\n", encoding="utf-8")
    perms = A.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
    return perms, tmp_path


def _shown(perms, command, stdout="A = 1\n", exit_code=0):
    A._baseline_shell_reads(
        perms, {"command": command},
        json.dumps({"exit_code": exit_code, "stdout": stdout, "stderr": ""}))


def _has_baseline(perms, root):
    return str((root / "cfg.py").resolve()) in perms.read_tracker


@pytest.mark.parametrize("command", [
    "sed -n '1,2p' cfg.py",
    "cat cfg.py",
    "head -n 5 cfg.py",
    "grep -n A cfg.py; sed -n '1,2p' cfg.py",
])
def test_a_printed_file_has_a_baseline(ws, command):
    perms, root = ws
    _shown(perms, command)
    assert _has_baseline(perms, root)


@pytest.mark.parametrize("command", [
    "wc -l cfg.py",
    "stat cfg.py",
    "grep -n A cfg.py",
    "sed -i 's/A/B/' cfg.py",
    "sed -n -i 's/A/B/p' cfg.py",
    "cat /dev/null > cfg.py",
])
def test_what_does_not_show_the_text_sets_none(ws, command):
    perms, root = ws
    _shown(perms, command)
    assert not _has_baseline(perms, root)


def test_a_failed_or_silent_command_sets_none(ws):
    perms, root = ws
    _shown(perms, "cat cfg.py", exit_code=1)
    _shown(perms, "cat cfg.py", stdout="")
    assert not _has_baseline(perms, root)


def test_the_edit_goes_through_after_a_shell_read(ws):
    perms, root = ws
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    refused = executor._execute_edit_file(
        {"path": "cfg.py", "old_string": "A = 1", "new_string": "A = 3"}, perms)
    assert "before editing" in refused

    _shown(perms, "sed -n '1,2p' cfg.py")
    out = executor._execute_edit_file(
        {"path": "cfg.py", "old_string": "A = 1", "new_string": "A = 3"}, perms)
    assert "error" not in out[:20], out
    assert (root / "cfg.py").read_text(encoding="utf-8").startswith("A = 3")
