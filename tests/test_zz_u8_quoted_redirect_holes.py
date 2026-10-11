"""U8 (operator build): a quoted ``>`` is text only where nothing runs it.

The quoted-redirect fix must not hide a write that a nested shell, a command
substitution or a script-carrying argument performs.
"""

import pytest

from delfin.agent.api_client import _bash_write_targets


@pytest.mark.parametrize("cmd,target", [
    ('echo "$(echo hi > out.txt)"', "out.txt"),
    ("bash -c 'echo hi > out.txt'", "out.txt"),
    ('sh -c "echo hi > out.txt"', "out.txt"),
    ('eval "echo hi > out.txt"', "out.txt"),
    ('git commit -m "$(echo x > out.txt)"', "out.txt"),
    ("git -c alias.z='!echo > out.txt' z", "out.txt"),
    ("git config alias.z '!echo > out.txt'", "out.txt"),
    ("env git commit -m 'a > out.txt'", "out.txt"),
    ("git commit -m 'x => y' > out.txt", "out.txt"),
])
def test_a_write_behind_quotes_is_still_named(cmd, target):
    assert any(target in t for t in _bash_write_targets(cmd)), cmd


@pytest.mark.parametrize("cmd", [
    'git commit -m "a -> b"',
    "git commit -m 'x => new_value'",
    'echo "a -> b"',
    'grep -E "a=>b" notes.txt',
    'git log --grep "x -> y"',
])
def test_an_arrow_in_a_message_is_not_a_write(cmd):
    assert _bash_write_targets(cmd) == [], cmd
