"""IO2 adversarial tests (reviewer nacht-s2), round 2: the fix's own parsers.

Builder nacht-s1's 763698478 fixed F1-F5 of my first round and added three
new parsers: ``_strip_wrappers`` (wrapper prefixes), ``_collect_vars``
(script assignments) and a per-line ``cd`` rule. New parsers get the same
ALLOW-list attack as the old ones: the shell spellings that actually appear
in submit scripts, not the ones the regex was written against.

Each test states the value measured on 763698478. F1-F3 and F5 are NOT
re-tested here -- they are fixed and green in
``tests/test_io2_adversarial_detection.py``.

Reviewed: builder nacht-s1, commit 763698478 (IO2 phase 2).
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/pfs/home/user"}


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


def findings(script, env=None):
    return [(f.line, f.severity, f.rule) for f in hpc_io_check.check(
        script, env=ENV if env is None else env)]


def refusals(script, env=None):
    return [f for f in findings(script, env=env) if f[1] == "refuse"]


# --- 6. Wrappers that carry their own arguments --------------------------

@pytest.mark.parametrize("script", [
    # nice's normal spelling: an adjustment before the command.
    "nice -n 10 find /pfs/data6 -type f\n",
    "nice -10 find /pfs/data6 -type f\n",
    # env with an assignment is the common submit-script prefix; the
    # builder's `_strip_wrappers` strips `env` but only re-strips VAR=val
    # prefixes once, before its loop.
    "env X=1 find /pfs/data6 -type f\n",
    "env -i find /pfs/data6 -type f\n",
])
def test_a_wrapper_that_takes_arguments_does_not_hide_the_walk(lustre, script):
    """``_strip_wrappers`` (hpc_io_check.py:123-147) strips a bare
    ``nice``/``env`` but not the spelling that carries its own argument, so
    the command word no longer anchors the rule. Measured ``[]`` for all
    four; the bare forms refuse. ``nice -n 10`` is how nice is normally
    written, so this is the spelling that will actually appear."""
    assert refusals(script)


@pytest.mark.parametrize("script", [
    # A POSIX shell's `-c` form takes a name for $0 after the command.
    "bash -c 'find /pfs/data6 -type f' zero\n",
    "sh -c 'find /pfs/data6 -type f' \"$0\"\n",
])
def test_bash_c_with_a_zeroth_argument_still_refuses(lustre, script):
    """``bash -c 'cmd' name`` is a documented shell form (the name becomes
    $0). The wrapper regex (:132-133) requires the quoted command to end the
    line, so it does not match and the walk is invisible. Measured ``[]``."""
    assert refusals(script)


# --- 7. The per-line cd rule and chained cd ------------------------------

def test_a_second_cd_on_the_same_line_is_the_one_that_counts(lustre):
    """F2's fix looks at the last ``cd`` on or before the pytest line
    (763698478, commit message), but a chained cd is two cds on ONE line.
    ``cd /tmp && cd /pfs/data6/w`` leaves the shell on Lustre, so the pytest
    walks a Lustre tree. Measured ``[]`` -- the node-local /tmp is taken as
    the cwd."""
    assert refusals("cd /tmp && cd /pfs/data6/w\npytest tests/\n")


def test_the_last_cd_in_a_chain_to_lustre_first(lustre):
    """The mirror image, to show which direction the bug runs: Lustre first,
    then node-local, must stay allowed (the shell really is on /tmp)."""
    assert not refusals("cd /pfs/data6/w && cd /tmp\npytest tests/\n")


# --- 8. What the fix got right, pinned so it is not lost -----------------

@pytest.mark.parametrize("script", [
    'D=/pfs/data6\nfind "$D/x" -type f\n',
    'D="/pfs/data6"\nfind "$D/x" -type f\n',
    'export D=/pfs/data6\nfind "$D/x" -type f\n',
])
def test_a_script_assignment_no_longer_hides_the_path(lustre, script):
    """F3's fix: an assignment in the script itself now expands. Green on
    763698478; kept so a later refactor cannot silently lose it."""
    assert refusals(script)


def test_the_genuine_node_local_fix_is_still_allowed(lustre):
    assert not refusals("cd $TMPDIR\ntar -xf env.tar\npytest tests/\n")


def test_line_numbers_are_right_after_a_continuation(lustre):
    """F5's fix: the finding names the source line, past a backslash
    continuation. Green on 763698478; pinned because H2 depends on it."""
    script = ("#!/usr/bin/env bash\n"
              "A=1\n"
              "B=2 \\\n"
              "  C=3\n"
              "find /pfs/home/user -type f\n")
    assert [n for n, _, _ in refusals(script)] == [5]
