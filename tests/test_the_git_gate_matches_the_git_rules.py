"""The git rules were prose; the gate permitted their opposite.

The prompt pack states the discipline plainly — branch before you build,
push the branch before any merge, stage what you changed and not the
whole tree, look before destructive git. A forensic pass compared each
sentence against the permission tables and found the tables silent or
inverted:

* ``git add -A`` and ``git add .`` were auto-allowed, unattended. Worse,
  the solo prompt *instructed* the bulk stage as a "checkpoint" and a
  test asserted it must stay auto-allowed. In a tree that also holds
  work this session did not do, that commit takes somebody else's
  changes with it and afterwards nothing can say which hunks were whose.
* ``git checkout -- <path>`` was auto-allowed while ``git reset --hard``
  was denied — two commands named as destructive in the same sentence of
  the project's own rules, landing on opposite lists.
* ``git clean -fd`` was denied and ``git clean -f -d`` was not. The
  deny pattern required the letters in one token and in that order, so
  the spellings people actually type went to the confirm gate — which
  under an unattended profile means they ran.
* ``git commit --author=`` was auto-allowed: the entry anchors on the
  message flag and ignores the rest of the line, so a commit under
  someone else's name needed no approval.
* Approving one ``git push -u origin feature`` persisted the pattern
  ``^\\s*git\\s+push\\b``. From that click on, ``git push origin main``
  auto-ran in every future session — the one invariant this area has a
  test for, revoked by a dialog that asked about a feature branch.

What is asserted here is the boundary, not the mood: the routine half of
each command still runs without asking, and only the half that can
destroy work this session did not create goes to the gate.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent import api_client as A
from delfin.agent import kit_settings as KS

_DENY = A._DEFAULT_BASH_DENY_PATTERNS
_ALLOW = A._DEFAULT_BASH_AUTO_ALLOW


def _denied(cmd: str) -> bool:
    return any(re.search(p, cmd, re.IGNORECASE) for p in _DENY)


def _auto(cmd: str) -> bool:
    return any(re.search(p, cmd, re.IGNORECASE) for p in _ALLOW)


def _verdict(cmd: str) -> str:
    if _denied(cmd):
        return "deny"
    return "auto" if _auto(cmd) else "gate"


# ---------------------------------------------------------------------------
# Staging: named paths run, the whole tree asks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git add -A",
    "git add .",
    "git add -u",
    "git add --all",
    "git add --update",
    'git add -A && git commit -m "checkpoint"',
])
def test_a_bulk_stage_is_not_unattended(cmd):
    assert _verdict(cmd) != "auto", (
        f"{cmd!r} stages work this session may not have done")


@pytest.mark.parametrize("cmd", [
    "git add delfin/agent/office.py",
    "git add tests/ delfin/",
    "git add README.md",
])
def test_staging_named_paths_still_runs(cmd):
    assert _verdict(cmd) == "auto"


def test_the_prompt_no_longer_teaches_the_bulk_stage():
    """The instruction and the permission were contradicting each other;
    fixing one without the other leaves the model asking for something
    it will not get."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "pack" / "agents" / "solo_agent.md").read_text(encoding="utf-8")
    assert "git add -A &&" not in src, "the checkpoint recipe is still there"
    # Lowercase both sides or neither: lowering only the source turns its
    # `-A` into `-a` and the needle never matches.
    assert "never `git add -a`" in src.lower(), (
        "the prompt has to say what to do instead, not merely stop saying "
        "the wrong thing")


# ---------------------------------------------------------------------------
# Moving between branches runs; discarding the working tree asks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git checkout -- .",
    "git checkout .",
    "git checkout -- delfin/agent/office.py",
    "git restore .",
    "git restore delfin/agent/office.py",
    "git switch --discard-changes main",
])
def test_discarding_uncommitted_work_is_not_unattended(cmd):
    assert _verdict(cmd) != "auto", (
        f"{cmd!r} overwrites uncommitted work in place")


@pytest.mark.parametrize("cmd", [
    "git checkout main",
    "git checkout -b topic",
    "git checkout feature/office-fixtures",
    "git switch main",
    "git switch -c topic",
])
def test_moving_between_branches_still_runs(cmd):
    assert _verdict(cmd) == "auto"


# ---------------------------------------------------------------------------
# The stash is not a scratch pad the agent may empty
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", ["git stash", "git stash drop", "git stash clear",
                                 "git stash pop"])
def test_touching_the_stash_is_not_unattended(cmd):
    assert _verdict(cmd) != "auto"


@pytest.mark.parametrize("cmd", ["git stash list", "git stash show"])
def test_looking_at_the_stash_still_runs(cmd):
    assert _verdict(cmd) == "auto"


# ---------------------------------------------------------------------------
# git clean: the flags, not one spelling of them
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git clean -fd",
    "git clean -f -d",
    "git clean -xdf",
    "git clean -d -f",
    "git clean --force -d",
    "git clean -f -x",
    "git clean -ffdx",
])
def test_every_spelling_of_a_destructive_clean_is_denied(cmd):
    assert _verdict(cmd) == "deny", f"{cmd!r} deletes untracked files"


@pytest.mark.parametrize("cmd", ["git clean -n", "git clean -nd", "git clean -d"])
def test_a_clean_that_deletes_nothing_is_not_denied(cmd):
    """`-d` without `-f` refuses on git's own side; a dry run is a read."""
    assert _verdict(cmd) != "deny"


# ---------------------------------------------------------------------------
# A commit carries this session's identity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    'git commit -m "x" --author="Someone <a@b.c>"',
    'git commit --author "Someone <a@b.c>" -m "x"',
    'git commit --date="2020-01-01" -m "x"',
])
def test_a_forged_identity_is_denied(cmd):
    assert _verdict(cmd) == "deny"


def test_an_ordinary_commit_still_runs():
    assert _verdict('git commit -m "fix the reader"') == "auto"


# ---------------------------------------------------------------------------
# One approval does not become a standing grant for a heavier command
# ---------------------------------------------------------------------------

def test_approving_a_feature_push_does_not_grant_pushing_to_main():
    pattern = KS.suggest_pattern_for_command("git push -u origin feature")
    assert pattern, "a suggestion is still expected"
    for heavier in ("git push origin main",
                    "git push --force-with-lease origin main"):
        assert not re.search(pattern, heavier), (
            f"remembering the feature push would also auto-allow {heavier!r}")


def test_approving_one_commit_does_not_grant_a_forged_author():
    pattern = KS.suggest_pattern_for_command('git commit -m "wip"')
    assert not re.search(
        pattern, 'git commit -m "x" --author="Someone <a@b.c>"')


def test_approving_one_install_does_not_grant_every_install():
    pattern = KS.suggest_pattern_for_command("pip install ./local-wheel.whl")
    assert not re.search(pattern, "pip install requests-from-the-internet")


@pytest.mark.parametrize("cmd,covers", [
    ("git status", "git status --porcelain"),
    ("git diff HEAD~1", "git diff --stat"),
    ("git log --oneline", "git log -20"),
])
def test_harmless_subcommands_still_generalise(cmd, covers):
    """Narrowing everything would bury the user in dialogs for `git log`."""
    pattern = KS.suggest_pattern_for_command(cmd)
    assert re.search(pattern, covers)


# ---------------------------------------------------------------------------
# ...and the invariants that already held still hold
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git push origin main", "git push", "git push -u origin feature",
    "git merge --no-ff topic", "git rebase -i HEAD~3",
])
def test_publishing_and_rewriting_still_go_through_the_gate(cmd):
    assert _verdict(cmd) != "auto"


@pytest.mark.parametrize("cmd", [
    "git push --force origin main", "git push -f origin main",
    "git reset --hard", "git branch -D topic", "git tag -d v1",
])
def test_the_destructive_set_is_still_denied(cmd):
    assert _verdict(cmd) == "deny"


@pytest.mark.parametrize("cmd", [
    "git status", "git diff", "git log --oneline -5", "git fetch",
    "git show HEAD", "git blame delfin/agent/office.py",
])
def test_reading_the_repository_still_runs(cmd):
    assert _verdict(cmd) == "auto"
