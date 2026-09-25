"""Text between quotes is read the way bash reads it.

Four read-only commands from the supervised run of 2026-09-25 each cost
the operator a dialog, although nothing in them runs anything:

* a backtick between SINGLE quotes (`grep -n '^```' README.md`) was
  taken for command substitution -- bash expands nothing there;
* `\\|` inside a double-quoted grep pattern followed by `python` was
  taken for a pipe into an interpreter;
* an escaped `\\"` inside double quotes ended the quote early, so the
  rest of the pattern split into pseudo-commands.

The other half is the point of the gate and must not move: a
substitution between DOUBLE quotes runs, a real pipe into python is an
interpreter, and a backtick outside quotes is a command.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as ac
from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _asks(tmp_path, command: str) -> bool:
    asked: list[str] = []
    perms = KitToolPermissions(
        workspace=tmp_path, mode="acceptEdits",
        confirm_callback=lambda tool, args, reason: (asked.append(reason), False)[1])
    _doc_executor._run_permission_gate("bash", {"command": command}, perms)
    return bool(asked)


@pytest.mark.parametrize("command", [
    "grep -n '^```' README.md",
    "grep -c '```' README.md docs/USER_MANUAL.md",
    'grep -rln "sys.executable\\|python" tests/ | head -40',
    'grep -n "anyOf\\|not\\":\\|format\\":" delfin/agent/api_client.py | head -20',
])
def test_a_quoted_pattern_is_read_without_asking(tmp_path, command):
    assert not _asks(tmp_path, command)


@pytest.mark.parametrize("command", [
    'grep x "$(id)"',
    'grep x "`id`"',
    "grep x `id`",
    "ls $(touch x)",
    "grep -n x README.md | python3",
    "cat README.md | sh",
    "grep 'a' README.md; python3 -c 'import os'",
])
def test_what_runs_something_still_asks(tmp_path, command):
    assert _asks(tmp_path, command)


def test_an_escaped_quote_keeps_the_string_open():
    # bash: one argument `a" ; rm x` -- no second command.
    assert ac._split_shell_segments('echo "a\\" ; rm x"') == ['echo "a\\" ; rm x"']
    # an unquoted escaped operator is a character too
    assert ac._split_shell_segments("grep a\\|b f") == ["grep a\\|b f"]
    # the real operators still split
    assert ac._split_shell_segments('grep "a" f | head -3') == ['grep "a" f', "head -3"]


def test_blanking_keeps_offsets_and_double_quoted_substitution():
    cmd = "grep '`x`' \"$(y)\" z"
    blanked = ac._blank_quoted(cmd, "'")
    assert len(blanked) == len(cmd)
    assert "`" not in blanked
    assert "$(y)" in blanked
