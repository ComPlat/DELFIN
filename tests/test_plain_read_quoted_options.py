"""A quoted option is the same option to the shell.

``_bash_is_plain_read`` matched git's writing options (--output, -o,
--ext-diff, ...) and destructive flags against ``str.split()`` words.
The shell removes quotes before git sees the argument, so
``git log '--output=/tmp/x'`` writes /tmp/x while the split word
``'--output=/tmp/x'`` starts with a quote and matched no prefix: the
command passed as a plain read and never reached the classifier.
"""

import pytest

from delfin.agent import action_classifier as ac


@pytest.mark.parametrize("command", [
    "git log '--output=/tmp/x'",
    'git log "--output=/tmp/x"',
    "git diff '-o' out.txt",
    'git show "-oout.txt" HEAD',
    "git log '--ext-diff'",
    "git diff \"--textconv\"",
    "git branch '--set-upstream-to=origin/main'",
    "git branch '-D' old",
    "git tag \"-d\" v1",
    "git log --'output'=/tmp/x",
])
def test_quoted_writing_option_is_not_a_plain_read(command):
    assert ac._bash_is_plain_read(command) is False
    assert ac.is_plain_read({"tool": "bash", "command": command}) is False


@pytest.mark.parametrize("command", [
    "git log '--oneline' -3",
    "git log --grep 'fix this'",
    "git status",
    "grep -n 'two words' delfin/agent/action_classifier.py",
])
def test_quoted_read_stays_a_plain_read(command):
    assert ac._bash_is_plain_read(command) is True


def test_unbalanced_quote_is_not_a_plain_read():
    assert ac._bash_is_plain_read("git log '--output=/tmp/x") is False
