"""The prompt forbids `cd`, and the gate only understood `cd`'s absence.

`solo_agent.md` says it twice: never prepend `cd /pfad && …` to a bash
command; use the tool's `cwd`, or git's own `-C`. The auto-allow list
matched `git\\s+<subcommand>` with nothing in between, so `git -C . status`
— the spelling the prompt mandates — went to the confirm gate, and in a
headless run that is a refusal. So did `git --no-pager log`, which is how
you keep git from opening a pager in a shell that has no terminal.

Seen 2026-09-08 in a suite run: a compound whose every segment was
allowed on its own was refused for the `-C` in the middle of one of them.

Read-only subcommands only. `-C` names ANOTHER repository, so
`git -C /elsewhere commit` would be a write over there; the writing
subcommands keep requiring the bare form and the gate that goes with it.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


@pytest.mark.parametrize("cmd", [
    "git -C . status",
    "git -C . status --short",
    "git -C /tmp/other rev-parse --show-toplevel",
    "git -C .. diff",
    "git --no-pager log --oneline -5",
    "git --no-pager diff",
    "git --paginate show HEAD",
    "git -C . --no-pager log -3",
    # the bare forms keep working
    "git status",
    "git rev-parse --abbrev-ref HEAD",
    "git log --oneline -5",
])
def test_reading_a_repository_runs_however_it_is_spelled(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is True, cmd


@pytest.mark.parametrize("cmd", [
    # -C names another repository; writing there is not a read.
    "git -C /tmp/other commit -m x",
    "git -C /tmp/other push",
    "git -C /tmp/other add file.py",
    "git -C /tmp/other checkout -- .",
    "git -C /tmp/other init",
    # and push is never auto-allowed, in any spelling
    "git push",
    "git push -u origin main",
    "git --no-pager push",
])
def test_writing_still_goes_through_the_gate(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is False, cmd


def test_the_flag_prefix_is_not_a_place_to_hide_a_command(perms):
    """`-C` takes a path. A substitution there would be command text."""
    for cmd in ("git -C $(whoami) status",
                "git -C `id` status",
                "git -C ;whoami status"):
        assert perms.matches_bash_auto_allow(cmd) is False, cmd
