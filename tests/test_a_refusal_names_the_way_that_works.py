"""A denied command is told what to do instead, where there is one.

A session made a temporary worktree to check a baseline, tried to remove
it again, and was told only "command rejected by deny-pattern ...:
refusing to run". The worktree verbs do exactly that job with the checks
a shell cannot make, and the refusal never mentioned them; the directory
stayed behind (2026-09-17).

A deny-list that names nothing leaves the model to invent a way around,
which is the one thing it must not provoke.

  the worktree verbs are named   for a denied `git worktree remove`
  so are the others with a route reset, clean, branch, push
  the refusal still refuses      the hint is added, never substituted
  no route, no invention         ending the machine gets no detour
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import (KitToolPermissions, _denied_command_hint,
                                     _doc_executor)


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path, mode="bypassPermissions",
                              confirm_callback=None)


def _gate(cmd, perms):
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


def test_a_denied_worktree_removal_names_the_verb(perms):
    out = _gate("git worktree remove /tmp/check-base", perms)
    assert "refusing to run" in out
    assert "worktree_remove" in out


@pytest.mark.parametrize("cmd, word", [
    ("git reset --hard HEAD~1", "revert"),
    ("git clean -xdf", "one by one"),
    ("git branch -D feature", "user's to delete"),
    ("git push --force origin main", "pull request"),
])
def test_the_others_with_a_route_name_it(perms, cmd, word):
    out = _gate(cmd, perms)
    assert "refusing to run" in out and word in out


def test_ending_the_machine_gets_no_detour(perms):
    out = _gate("shutdown -h now", perms)
    assert "refusing to run" in out
    assert out.rstrip().endswith("refusing to run."), "no route to offer"


def test_the_hint_is_added_and_never_replaces_the_refusal(perms):
    out = _gate("git worktree remove /tmp/x", perms)
    assert out.startswith("command rejected by deny-pattern")


def test_a_pattern_nobody_has_advice_for_is_silent():
    assert _denied_command_hint(r"\bsudo\b") == ""
    assert _denied_command_hint("") == ""
    assert _denied_command_hint(None) == ""
