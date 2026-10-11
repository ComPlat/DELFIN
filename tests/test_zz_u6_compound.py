"""U6 phase-4 RED TEST (moved to .gate/ after the red control is confirmed).

A compound command is auto-allowed when EVERY segment (split outside quotes
at ; && |) is auto-allowed on its own and none of them writes outside the
workspace or reaches a remote. The gate already splits compounds
(matches_bash_auto_allow); the one real over-refusal is that read-only git
plumbing (`git merge-base`, `git rev-list`, `git cat-file`, ...) is not on
the auto-allow list even though _git_reads_only (api_client.py:1086
classifies the whole _GIT_READ_ONLY_SUBCOMMANDS set as read-only, line 1034).
That alone blocks the operator's headline compound. This test is RED on the
unpatched gate (merge-base -> auto=False); GREEN after the phase-4 patch adds
the missing read-only plumbing to the git allow alternation. The still-asked
shapes (git merge --no-edit writes the tree; git --output writes a file; $( )
runs a payload) must STAY not-auto-allowed.
"""
import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


@pytest.mark.parametrize("cmd", [
    "git merge-base HEAD X",
    "git log --oneline main -1; echo ---; git merge-base HEAD X",
    "git log -3; git merge-base main HEAD",
])
def test_read_only_git_plumbing_auto_allows(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "git merge X --no-edit 2>&1 | tail -5",   # git merge writes the tree
    "git log --output=/x",                     # --output writes a file
    "echo $(ls)",                              # command substitution runs a payload
])
def test_writes_and_payloads_stay_refused(perms, cmd):
    assert not perms.matches_bash_auto_allow(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "git grep -Ovim foo",
    "git grep -O'sh -c id' foo",
    "git grep --open-files-in-pager=sh foo",
    # git bundles short options and accepts a unique long-option prefix.
    "git grep -nOvim foo",
    "git grep -iOless foo",
    "git grep --open=less foo",
    "git grep --open-files=less foo",
])
def test_a_git_grep_that_runs_a_pager_program_is_not_auto_allowed(perms, cmd):
    assert not perms.matches_bash_auto_allow(cmd), cmd


def test_a_plain_git_grep_still_auto_allows(perms):
    assert perms.matches_bash_auto_allow("git grep -n foo -- delfin/")
