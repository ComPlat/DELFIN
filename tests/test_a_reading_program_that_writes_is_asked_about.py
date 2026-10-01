"""A reading program's name is not a reading command.

The gate auto-allows sort, uniq, sed -n, awk, find, tree, rg, git's
reading subcommands and env as READING tools. Each has a form that writes
a file the write-target gate never sees, or runs a program: found
2026-09-26 by probing the full permission gate with the spellings a
session's own read-only classifier had missed. Those forms go to the
confirm gate; the everyday reading forms stay free.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


WRITES_OR_RUNS = [
    "sed -n '1w /tmp/x' f.txt", "sed -n '/a/w out' f", "sed -n 's/a/b/w out' f",
    "sed -n '1e id' f", "sort -o out.txt f.txt", "sort -uo out f",
    "uniq f.txt out.txt", "env rm -rf f.txt", "env -i sh -c id",
    "rg --pre=sh foo .", "git -c core.fsmonitor=id status",
    "git diff --output=x.txt", "git log --ext-diff -p", "find . -fprint x",
    "find . -name '*.py' -exec cat {} +", "awk 'BEGIN{system(\"id\")}'",
    "awk '{print > \"f\"}' x", "awk '{print | \"sh\"}' x", "awk -f prog.awk x",
    "tree -o out .",
]

READS = [
    "sed -n '10,40p' f.txt", "sed -n '/def wait/,/^def /p' f.py", "sed -n '$p' f",
    "sed -n '1,5p;/west/p' f", "awk '{print $1}' f", "awk -F: '$2>5 {print $1}' f",
    "sort f", "uniq -c f", "find . -name '*.py'", "rg -n foo", "tree -L 2", "env",
    "git -C . log --oneline -5", "git diff", "git status", "grep -n wait f.py",
]


@pytest.mark.parametrize("cmd", WRITES_OR_RUNS)
def test_a_writing_or_running_form_is_not_auto_allowed(perms, cmd):
    assert not perms._segment_auto_allowed(cmd), cmd


@pytest.mark.parametrize("cmd", READS)
def test_the_reading_form_stays_free(perms, cmd):
    assert perms._segment_auto_allowed(cmd), cmd
