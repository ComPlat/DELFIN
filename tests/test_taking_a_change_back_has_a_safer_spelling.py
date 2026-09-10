"""`git checkout -- <file>` reverts by path. undo_changes reverts by author.

Four refusals in one suite tonight were the agent taking back a change
it had made to a fixture: `git checkout -- bookmarks.json`, `git
checkout -- bookmark_store.py bookmarks.json bookmarks_cli.py`. None of
them was told what to use instead.

It is the most dangerous spelling of that wish in a shared checkout,
which is why it is off the auto-allow list at all: it reverts by PATH,
not by authorship, so it discards whatever anybody else had uncommitted
under those paths — and this checkout has another agent working in it.

undo_changes reverts only what THIS session recorded writing, and
refuses a file whose content changed since. The refusal now says so.

`git stash` is deliberately not matched here. It is answered one branch
earlier by the control-run hint, which names enter_worktree(base_ref=…)
and states why that cannot destroy work in the user's tree.
"""

from __future__ import annotations

import json
import tempfile

import pytest

import delfin.agent.api_client as A

_n = [0]


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory() as tmp:
        yield tmp


def _err(ws, cmd, mode="default"):
    _n[0] += 1
    perms = A.KitToolPermissions(mode=mode, workspace=ws)
    perms.task_session_id = f"revert-{_n[0]}"
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms))
    return out.get("error", "")


@pytest.mark.parametrize("cmd", [
    "git checkout -- bookmarks.json",
    "git checkout -- bookmark_store.py bookmarks.json bookmarks_cli.py",
    "git checkout -- a.py b.py && ls",
    "git restore bookmarks.json",
    "git restore --staged --worktree a.py",
])
def test_taking_a_change_back_names_undo_changes(ws, cmd):
    err = _err(ws, cmd)
    assert err, cmd
    assert "undo_changes" in err, err


def test_the_hint_says_why_it_is_not_merely_the_allowed_spelling(ws):
    err = _err(ws, "git checkout -- bookmarks.json")
    assert "reverts by path" in err
    assert "anyone else had uncommitted" in err


def test_the_original_of_a_file_you_did_not_write_is_a_read(ws):
    """The other half of the wish: wanting the ORIGINAL is a read, not a
    restore over somebody's working tree."""
    err = _err(ws, "git checkout -- bookmarks.json")
    assert "read_file" in err


@pytest.mark.parametrize("cmd", [
    "git checkout main",
    "git checkout -b feature",
    "git status",
    "git diff",
    "git log --oneline -5",
])
def test_the_git_commands_that_were_always_fine_are_untouched(ws, cmd):
    assert not _err(ws, cmd), cmd


def test_git_stash_keeps_the_control_run_hint(ws):
    """Answered one branch up, and better: enter_worktree(base_ref=…) is
    the sanctioned control, and the stash is shared between agents."""
    err = _err(ws, "git stash")
    assert "enter_worktree" in err, err
    assert "undo_changes" not in err


def test_a_checkout_of_a_ref_is_not_a_revert(ws):
    """`git checkout <ref>` moves HEAD; it is not the shape this hint is
    about, and matching it would misdescribe what the command does."""
    err = _err(ws, "git checkout 74a890ad")
    assert "undo_changes" not in err
