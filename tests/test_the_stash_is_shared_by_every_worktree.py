"""The stash belongs to the repository. The worktrees do not know that.

``git worktree`` gives each session its own files and its own HEAD. It
does not give it its own stash: the stack lives in the repository, so
four sessions in four worktrees push onto ONE pile, and a pop takes
whatever is on top — another session's work, applied into the wrong
tree, with nothing on screen to say so.

On 2026-09-18 a session ran

    git stash push -m tmpbase -- delfin/agent/repl.py ...

while three others were working in worktrees of the same repository. It
was allowed: the mode was all_free, so nothing asked, and the deny list —
which applies in every mode — did not name it.

It does now, anchored on the subcommand so that a commit message
mentioning the word is not refused, and the refusal names the ways that
work: commit on your own branch, or ``enter_worktree`` for a control run.
Reading the stack stays allowed.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import (
    _DEFAULT_BASH_DENY_PATTERNS, _denied_command_hint,
)
import re


def _denied(cmd: str) -> str:
    """The pattern that refuses *cmd*, or ""."""
    for pat in _DEFAULT_BASH_DENY_PATTERNS:
        if re.search(pat, cmd):
            return pat
    return ""


@pytest.mark.parametrize("cmd", [
    "git stash",
    "git stash push -m tmpbase -- delfin/agent/repl.py",
    "git stash pop",
    "git stash apply",
    "git stash drop",
    "git stash clear",
    "git stash branch fix",
    "git -C /elsewhere stash push",
    "git --no-pager stash pop",
    "cd /tmp && git stash push",
])
def test_every_stash_that_writes_is_refused(cmd):
    assert _denied(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "git stash list",
    "git stash show",
    "git stash list | head -3",
    "git -C /elsewhere stash list",
])
def test_reading_the_stack_is_still_allowed(cmd):
    assert not _denied(cmd), cmd


@pytest.mark.parametrize("cmd", [
    'git commit -m "drop the stash usage"',
    'git commit -m "stash"',
    "grep -rn 'git stash' delfin/",
    "echo 'the stash is shared' > notes.txt",
])
def test_the_word_in_a_message_is_not_the_command(cmd):
    """The shutdown pattern once refused a commit because the word stood
    in its message, and cost the session the commit. Not twice."""
    hit = _denied(cmd)
    assert not hit, f"{cmd!r} refused by {hit!r}"


def test_the_refusal_names_the_ways_that_work():
    hint = _denied_command_hint("git\\s+stash\\b")
    assert "commit it on your own branch" in hint
    assert "enter_worktree" in hint
    assert "list" in hint, "and that reading the stack is still available"


def test_the_hint_says_why_the_stash_is_not_private():
    hint = _denied_command_hint("git\\s+stash\\b")
    assert "REPOSITORY" in hint and "worktree" in hint
