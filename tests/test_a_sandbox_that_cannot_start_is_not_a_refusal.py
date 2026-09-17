"""When the cage fails to start, the command did not run.

A session staged two files with `git diff --stat && git add ...` and got
back exit 1 and "bwrap: Can't bind mount /oldroot/...". It read that as a
refusal and rewrote the command; a plain retry was what worked
(2026-09-17). The line comes from bubblewrap, not from git: nothing ran,
nothing was refused, nothing was changed.

  the sandbox said it           bwrap / sandbox-exec on a failed exit
  the command's own failure      is left alone
  a successful run               is never annotated
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import _sandbox_start_note


def _result(exit_code: int, stderr: str = "", stdout: str = "") -> str:
    return json.dumps({"exit_code": exit_code, "stderr": stderr,
                       "stdout": stdout})


@pytest.mark.parametrize("stderr", [
    "bwrap: Can't bind mount /oldroot/pfs/data6 on /newroot/pfs/data6",
    "bwrap: Can't mkdir parents for /newroot/home/x/.ssh",
    "sandbox-exec: sandbox_apply: Operation not permitted",
])
def test_a_cage_that_could_not_start_says_so(stderr):
    note = _sandbox_start_note(_result(1, stderr))
    assert "never ran" in note
    assert "not a refusal" in note
    # And it says what to do: repeat it, do not invent a way around it.
    assert "again as it was" in note


@pytest.mark.parametrize("stderr", [
    "fatal: not a git repository",
    "error: pathspec 'x' did not match any file",
    "Traceback (most recent call last):",
    "",
])
def test_the_commands_own_failure_is_left_alone(stderr):
    assert _sandbox_start_note(_result(1, stderr)) == ""


def test_a_command_that_worked_is_not_annotated():
    assert _sandbox_start_note(_result(0, "bwrap: warning about something")) == ""


def test_junk_is_not_a_sandbox_failure():
    assert _sandbox_start_note("not json at all") == ""
    assert _sandbox_start_note(None) == ""
    assert _sandbox_start_note(json.dumps({"exit_code": None})) == ""
