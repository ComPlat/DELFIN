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


@pytest.mark.parametrize("cmd", [
    "echo x >&out.txt",
    "echo x >& out.txt",
    "echo x >|out.txt",
    "echo x >| out.txt",
    "echo x 2>&1 >&out.txt",
])
def test_a_redirect_spelled_with_an_ampersand_or_a_bar_names_its_file(cmd):
    assert _bash_write_targets(cmd) == ["out.txt"], cmd


@pytest.mark.parametrize("cmd", [
    "echo x >&2", "echo x 2>&1", "echo x >&-", "echo x 3>&1 1>&2 2>&3",
    "make 2>&1 | tail -5",
])
def test_a_descriptor_duplicate_is_not_a_file(cmd):
    assert _bash_write_targets(cmd) == [], cmd
