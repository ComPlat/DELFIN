"""R2 finding 2 control: readonly_pipeline.is_safe recognises read-only
pipelines/sequences and refuses everything else.

Finding 2: harmless read-only forms cost approvals -- `grep x f 2>&1 | tail -5`,
`git log ...; git show ...`. The gate required every segment to match
individually, so a whole pipeline fell to the confirm gate. is_safe
answers "is this WHOLE command read-only" for the wiring (api_client
patch, operator-built).

Conservativeness is the point: it must approve the harmless read-only
forms it exists for, and refuse -- for the right reasons -- every write,
runner, command-substitution and redirection form in the hostile corpus.
"""

from __future__ import annotations

import pytest

from delfin.agent.readonly_pipeline import is_safe


# -- the forms the finding names, and kindred read-only pipelines ------
SAFE = [
    "ls -l",
    "cat file1 file2",
    "head -5 foo.out | tail -3",
    "grep ERROR *.out | wc -l",
    "grep ERROR file 2>&1 | tail -5",
    "git log --oneline -5; git show HEAD",
    "git status",
    "git diff",
    "git rev-parse --is-inside-work-tree",
    "git log --oneline | head -20",
    "git show HEAD:path/to/file.py",
    "git -C /someworktree log --oneline -3",
    "sort file | uniq -c",
    "ls | grep foo",
    "cat 'file with spaces'",
    "grep 'a|b' file",
    "LANG=C grep ERROR *.out",
    "echo hello",
    "grep -c ERROR file",
    "wc -l *.out | tail -1",
    "set -o pipefail; grep ERROR file | tail -5",
    "cat < file",
    "ls -la | sort -k5nr | head",
]

# -- every construct that writes, runs, or substitutes must lose ------
UNSAFE = [
    # write redirection to a file
    "cat file > out",
    "echo x >> log",
    "cat file 2> err",
    "cmd 2>&1 2> err2",
    "ls > /dev/null",
    # interpreters and command-runners
    "python -c 'x=1'",
    "python3 script.py",
    "ls | python3 script.py",
    "cat file | bash",
    "eval ls",
    "xargs echo hi",
    "env grep x file",
    "nohup ls",
    "sudo cat file",
    # sed / awk / find: write or run capability, never approved
    "sed -i 's/x/y/g' file",
    "awk '{print > \"out\"}' file",
    "find . -delete",
    "find . -exec rm {} ;",
    # command substitution, backticks, process substitution
    "grep $(cat file) x",
    "cat `ls`",
    "cat <(python script)",
    # writing git
    "git push",
    "git add file",
    "git commit -m x",
    "git checkout -- file",
    "git reset --hard",
    "git stash",
    # not whitelisted at all
    "rm file",
    "mv a b",
    "cp a b",
    "touch file",
    "mkdir -p dir",
    # malformed / empty
    "",
    "   ",
]
#: order-independent: (spec, expect) for a few that are not obvious
EXPLICIT = [
    ("cat file | tee out", False),      # tee writes
    ("ls | grep python", True),         # runner as argument is inert text
    ("git help log", False),            # not a read-only subcommand
    ("git log -p | grep -i fix", True),
]


@pytest.mark.parametrize("cmd", SAFE)
def test_readonly_pipelines_are_safe(cmd):
    assert is_safe(cmd) is True, cmd


@pytest.mark.parametrize("cmd", UNSAFE)
def test_writing_running_substituting_commands_are_refused(cmd):
    assert is_safe(cmd) is False, cmd


@pytest.mark.parametrize("cmd,expected", EXPLICIT)
def test_explicit_spec(cmd, expected):
    assert is_safe(cmd) is expected, cmd


def test_non_string_is_never_safe():
    assert is_safe(None) is False          # type: ignore[arg-type]
    assert is_safe(42) is False            # type: ignore[arg-type]


def test_quoted_separator_does_not_split():
    # a quoted pipe is data, not a pipeline split
    assert is_safe("grep 'a|b' file") is True
    # a quoted semicolon in an arg must not create a pseudo-command
    assert is_safe("git log -p ';'") is True
    assert is_safe("cat 'a;rm x'") is True


def test_sequence_safety_requires_every_member():
    assert is_safe("ls; rm x") is False
    assert is_safe("ls && python -c x") is False
    assert is_safe("grep x f || echo fallback") is True    # both read-only
    assert is_safe("grep x f || rm x") is False
