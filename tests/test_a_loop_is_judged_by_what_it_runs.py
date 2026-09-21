"""A shell loop is judged by the commands in it, not by the word "for".

The auto-allow list reads a segment's first word, and a loop's first
word is `for`. So `for f in *.out; do grep ERROR "$f"; done` stopped for
approval although every command in it is allowed — and the split on `;`
turned the loop into "for f in *.out", "do grep ERROR" and "done", none
of which is a command at all. Measured over the whole trace history:
262 of 759 refusals of this kind were loops, more than a third.

  a loop of allowed commands       runs, like the commands would
  a loop with anything else        asks, exactly as before
  a header that runs something     asks: $(...) is a command nobody read
  a while that tests a command     the same
  the deny-list still decides      it runs before this, and stays first
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import (KitToolPermissions, _doc_executor,
                                     _loop_body)


@pytest.fixture
def perms(tmp_path):
    (tmp_path / "a.py").write_text("x = 1\n")
    return KitToolPermissions(workspace=tmp_path, mode="default",
                              confirm_callback=None)


@pytest.mark.parametrize("cmd", [
    'for f in *.py; do wc -l $f; done',
    'for f in *.out; do grep ERROR "$f"; done',
    "for i in 1 2 3; do echo $i; done",
    "for d in */; do ls $d; done",
    "while read line; do echo $line; done",
])
def test_a_loop_of_allowed_commands_runs(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is True


@pytest.mark.parametrize("cmd", [
    "for f in *; do rm -rf $f; done",                 # not an allowed command
    "for f in $(ls); do echo $f; done",               # header runs something
    "for f in `ls`; do echo $f; done",
    "while ps aux; do echo x; done",                  # condition is a command
    "for f in *; do echo $f; pip install evil; done",  # one bad command is enough
])
def test_anything_else_still_asks(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is False


def test_the_deny_list_decides_before_any_of_this(perms):
    """The list only says whether to ASK. What must never run is refused
    a layer earlier, inside a loop as anywhere else."""
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "for f in *; do cat $f > /etc/passwd; done"}, perms)
    assert out is not None and "deny-pattern" in out

    out = _doc_executor._run_permission_gate(
        "bash", {"command": "for f in *; do cat $f > ~/.ssh/authorized_keys; done"},
        perms)
    assert out is not None and "secret-deny path" in out


def test_an_allowed_loop_reaches_the_shell(perms):
    assert _doc_executor._run_permission_gate(
        "bash", {"command": 'for f in *.py; do wc -l "$f"; done'}, perms) is None


# -- the reading itself -----------------------------------------------------

def test_the_body_is_what_comes_back():
    assert _loop_body("for f in *.py; do wc -l $f; echo done; done") == \
        ["wc -l $f", "echo done"]


@pytest.mark.parametrize("cmd", [
    "ls -la",
    "grep -rn foo .",
    "",
    "for f in $(ls); do echo $f; done",
    "while ps aux; do echo x; done",
])
def test_what_is_not_a_readable_loop_says_so(cmd):
    assert _loop_body(cmd) is None
