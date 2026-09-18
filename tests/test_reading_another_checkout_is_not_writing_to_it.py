"""`git -C <repo> log` reads. It was refused as a write.

``-C`` is on the list of options that carry a write destination, which is
right for ``tar -C`` and ``install -C`` and wrong for git: there it says
WHICH repository to work in, and whether that is a write is decided by the
subcommand.

Measured on 2026-09-18. A session asked

    git log --oneline -3 origin/main; ...; git -C <repo> log --oneline -3 main

and was told the repository "is in a READ-ONLY location" — true, and
beside the point, since nothing in that line writes anywhere. The call was
lost and the session was told the opposite of what was happening.

Fail closed: the mute applies only to a known read-only subcommand, and
only to ``-C``. Everything else about the command is judged as before.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _bash_write_targets, _git_subcommand


ELSEWHERE = "/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin"


# -- the reported command ---------------------------------------------------

def test_the_command_from_the_report_writes_nowhere():
    cmd = (f"git log --oneline -3 origin/main; git log --oneline -3 main; "
           f"git for-each-ref --count=5 refs/heads; "
           f"git -C {ELSEWHERE} log --oneline -3 main")
    assert _bash_write_targets(cmd) == []


@pytest.mark.parametrize("sub", [
    "log --oneline", "show HEAD", "status --porcelain", "diff HEAD~1",
    "rev-parse HEAD", "for-each-ref refs/heads", "ls-files", "blame x.py",
    "cat-file -p HEAD", "merge-base a b", "describe --tags",
])
def test_a_reading_subcommand_is_not_a_write(sub):
    assert _bash_write_targets(f"git -C {ELSEWHERE} {sub}") == []


# -- what must still be caught ----------------------------------------------

@pytest.mark.parametrize("sub", [
    "commit -m x", "add .", "checkout -- .", "reset --hard", "merge main",
    "rebase main", "clean -fd", "fetch origin", "pull", "init",
    # Ambiguous by design: each of these has a spelling that writes, so the
    # gate keeps treating the repository as a target.
    "branch -D old", "tag v1", "config user.name x", "stash pop",
    "worktree add /tmp/x", "remote add o url",
])
def test_a_writing_subcommand_still_names_the_repository(sub):
    assert ELSEWHERE in _bash_write_targets(f"git -C {ELSEWHERE} {sub}")


def test_an_output_option_on_a_reading_subcommand_still_counts():
    """The mute is for -C alone. A reading subcommand that is handed a
    write option writes."""
    out = _bash_write_targets(
        f"git -C {ELSEWHERE} diff --output=/etc/somewhere HEAD")
    assert "/etc/somewhere" in out


def test_tar_keeps_its_C():
    """-C is a real destination for the commands it was added for."""
    assert "/etc" in _bash_write_targets("tar -xf a.tar -C /etc")


# -- finding the subcommand past git's own options --------------------------

@pytest.mark.parametrize("args,expected", [
    (["log"], "log"),
    (["-C", "/x", "log"], "log"),
    (["--no-pager", "-C", "/x", "log"], "log"),
    (["-c", "core.pager=cat", "-C", "/x", "status"], "status"),
    (["--git-dir=/x/.git", "log"], "log"),
    (["--git-dir", "/x/.git", "commit"], "commit"),
    (["-C", "/x"], ""),
    ([], ""),
])
def test_the_subcommand_is_found_past_gits_own_options(args, expected):
    assert _git_subcommand(args) == expected


def test_a_value_that_looks_like_a_subcommand_is_not_one():
    """`-c log=1` carries a value; the value is not the subcommand."""
    assert _git_subcommand(["-c", "log", "commit"]) == "commit"
