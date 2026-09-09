"""Two ways past the interpreter guard, both a single word long.

The guard that keeps `python -c` off the bash auto-allow list anchors on
the head of a command or on what follows a shell operator. Nothing was
looking at what a wrapper word runs, and `_referenced_script_paths` --
which needs the same answer for a different question -- had had that list
for months:

    time python3 -c "open('/etc/passwd','a').write(1)"     auto-allowed

`time` is on the auto-allow list, `python3 -c …` after it was invisible
to the guard, and the write reached the file with no dialog. `nohup`,
`timeout 30`, `nice`, `setsid`, `stdbuf`, `command`, `watch` and a bare
`VAR=value` prefix all do the same thing.

The second was a word boundary in the wrong place. One alternative in
the guard's pattern ends in `=`:

    env\\s+[A-Za-z_]+=

and the `\\b` sat after the whole alternation, so it had to be satisfied
by whatever followed the `=`. `env X=1 python3 -c …` matched, because
`1` is a word character. `env PYTHONPATH=/elsewhere python3 -c …` did
not, because `/` is not -- and pointing PYTHONPATH somewhere else is the
entire reason the `env` alternative is in that list. The spelling that
matters was the one spelling it could not see. Same shape as the two
missing `\\b`s recorded on 2026-09-09, in the other direction: there a
boundary was absent, here one was present and load-bearing in a place
nobody meant it to be.

Both predate the payload analysis in
test_a_payload_in_the_command_is_not_hidden.py and are independent of
it: the analysis decides WHAT a readable payload does, this decides
whether the guard sees a payload at all.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


WRITE = "open('/etc/passwd','a').write(1)"


@pytest.mark.parametrize("prefix", [
    "time",
    "nohup",
    "timeout 30",
    "nice -n 5",
    "ionice -c 3",
    "setsid",
    "stdbuf -oL",
    "command",
    "watch",
    "env PYTHONPATH=/elsewhere",
    "env X=1",
    "PYTHONPATH=/elsewhere",
    "FOO=bar",
    "nohup time",                       # more than one
    "timeout 5 env PYTHONPATH=/x",
])
def test_a_wrapper_does_not_hide_the_interpreter(tmp_path, prefix):
    cmd = f"{prefix} python3 -c \"{WRITE}\""
    assert A._is_interpreter_invocation(cmd) is True, cmd
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    assert perms.matches_bash_auto_allow(cmd) is False, cmd


@pytest.mark.parametrize("cmd", [
    "env PYTHONPATH=/elsewhere python3 -c \"print(1)\"",
    "env PATH=/x:/y perl -e 'unlink $f'",
    "env HOME=~/other make build",
])
def test_the_env_alternative_matches_a_path_value(cmd):
    """The boundary bug, isolated: a value starting with `/` or `~`."""
    assert A._is_interpreter_invocation(cmd) is True, cmd


@pytest.mark.parametrize("cmd", [
    "ls -la",
    "time pytest -q",
    "timeout 30 pytest tests/",
    "nohup ./run.sh",
    "time python3 report.py",
])
def test_a_wrapper_in_front_of_something_ordinary_is_still_ordinary(cmd):
    """Peeling wrappers must not turn every wrapped command into an
    interpreter -- `time pytest` is `pytest`, and the whole point of the
    auto-allow list is that it runs.

    `env FOO=bar ls` is deliberately absent: `env VAR=` is its own entry
    in the guard's pattern and has always been opaque, whatever follows
    it. What changed is only that the entry now matches when the value
    starts with a `/`.
    """
    assert A._is_interpreter_invocation(cmd) is False, cmd


@pytest.mark.parametrize("cmd", [
    "time pytest -q",
    "timeout 30 pytest tests/",
    "time python3 report.py",
])
def test_the_wrapped_routine_case_still_runs_unattended(tmp_path, cmd):
    """`time` and `timeout` are the two wrappers the auto-allow list
    itself carries; nice/nohup/setsid are not on it and ask, before this
    change and after it."""
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    assert perms.matches_bash_auto_allow(cmd) is True, cmd


def test_a_wrapped_write_still_reaches_the_write_gate():
    """Seeing the interpreter is one thing; naming what it writes is the
    other, and the write gate reads the command wherever the payload sits
    in it."""
    assert A._bash_write_targets(
        f'time python3 -c "{WRITE}"') == ["/etc/passwd"]
    assert A._bash_write_targets(
        f'env PYTHONPATH=/x python3 -c "{WRITE}"') == ["/etc/passwd"]


def test_stripping_wrappers_terminates_on_junk():
    for junk in ("", "env", "time time time", "A=1 " * 500):
        A._strip_exec_wrappers(junk)
