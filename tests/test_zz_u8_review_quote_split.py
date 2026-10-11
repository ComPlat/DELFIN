"""U8 review: a quoted ``>`` is text only when the quote is read as the shell
reads it.

The exemption looked for quotes inside each segment of a split that ignores
quotes. A ``;``, ``|``, ``&`` or newline inside a quoted argument cut that
argument in two, the half after the cut began with an unmatched quote, and
a real redirection after it was taken for quoted text. ``$'a\\'b'`` (an
escaped quote inside ANSI-C quoting) and a backtick, which bash runs even
across single quotes, did the same. Each command below writes the file in
bash; each was reported with no write target before this change.
"""

import pytest

from delfin.agent.api_client import _bash_write_targets


@pytest.mark.parametrize("cmd", [
    'echo "a;echo " > /etc/zz',
    'echo "a|echo " > /etc/zz',
    'echo "a&echo " > /etc/zz',
    "echo 'a;echo ' > /etc/zz",
    'echo "a\necho " > /etc/zz',
    'git commit -m "x; git commit -m " > /etc/zz',
    'echo "a;echo " > /etc/zz ";b"',
    "echo $'a\\'b' > /etc/zz",
    "echo `echo '` > /etc/zz `'`",
    "cat <<A\nit's\nA\necho hi > /etc/zz",
    "cat <<A <<B\nit's\nA\nb\nB\necho hi > /etc/zz",
])
def test_a_redirect_after_a_quote_with_a_separator_is_named(cmd):
    assert "/etc/zz" in _bash_write_targets(cmd), cmd


@pytest.mark.parametrize("cmd", [
    'git commit -m "a -> b; c"',
    "git commit -m 'a -> b | c'",
    'echo "a > b & c"',
    "echo $'a > b'",
    "cat <<'A'\nit's a -> b\nA\ngit commit -m 'x -> y'",
])
def test_an_arrow_in_a_quote_with_a_separator_stays_text(cmd):
    assert _bash_write_targets(cmd) == [], cmd
