"""Operator review of readonly_pipeline: each line here was a read-only True.

A True skips the confirm dialog once the module is wired into the gate,
so every one of these would have run without asking.
"""
from __future__ import annotations

import pytest

from delfin.agent import readonly_pipeline as rp

WRITES_OR_RUNS = [
    "ls 3>f", "ls >&1.txt", "ls &>f", "cat <>f",
    "sort -of x", "sort -uof x", "sort --out=f x",
    "sort --compress-program=sh x", "sort -T /tmp x",
    "git diff --output=f", "git log --output f",
    "git reflog expire --all", "git branch -cnew",
    "git branch --set-upstream-to=o/x", "git -c core.pager=sh log",
    "git --exec-path=/x log", "git diff --ext-diff",
    "GIT_EXTERNAL_DIFF=evil git diff", "LD_PRELOAD=/x.so ls",
    "PAGER=sh git log", "file -C -m x", "date -s now",
]

TREE_WALKS = ["grep -r x ~ | head", "grep -rn x .", "grep --recursive x .",
              "ls -laR /"]

STILL_READ_ONLY = [
    "grep x f 2>&1 | tail -5", "git log --oneline -3; git show HEAD",
    "ls -la | sort -k5 -n | head", "sort -nr f | uniq -c", "git branch -a",
    "git tag -l", "ps aux | grep foo", "cut -d: -f1 f | sort -u",
    "echo a >&2", "git diff --stat HEAD~1",
]


@pytest.mark.parametrize("cmd", WRITES_OR_RUNS + TREE_WALKS)
def test_refused(cmd):
    assert rp.is_safe(cmd) is False


@pytest.mark.parametrize("cmd", STILL_READ_ONLY)
def test_still_read_only(cmd):
    assert rp.is_safe(cmd) is True
