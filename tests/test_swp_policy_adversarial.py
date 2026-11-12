"""Package SWP: adversarial policy gaps (reviewer s15 probe, hardened).

These cases come from the reviewer's adversarial probe. The red control
is the group here that the committed policy.py (eeb76a05a) does NOT yet
refuse -- the gaps in the "no patch except --check" allowance, the patch
redirection boundaries, home traversal, and unicode command separators.
Each fix in policy.py turns the corresponding red case green.
"""

import re
from pathlib import Path

import pytest

import delfin.agent.swarm.policy as policy

HOME = Path("/home/alice")


def _refuses(patterns, cmd):
    return any(re.search(p, cmd) for p in patterns)


def _allows(patterns, cmd):
    return not _refuses(patterns, cmd)


def test_git_apply_dash_check_is_allowed():
    # The one patch-application form a member may run is the dry run.
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _allows(pats, "git apply --check change.patch")


def test_git_apply_attached_check_flag_is_refused():
    # "--checkfoo" is not the dry run; it must not ride the allowance.
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, "git apply --checkfoo change.patch")
    assert _refuses(pats, "git apply -3 series.patch")


def test_git_apply_excluded_file_is_refused():
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, "git apply --exclude=utils old.patch")


ENV_PREFIX_REFUSE = [
    "X=y git push origin main",
    "PYTHONPATH=. python -m pytest tests/test_x.py",
    "X=1 pytest tests/test_y.py",
]


@pytest.mark.parametrize("cmd", ENV_PREFIX_REFUSE)
def test_env_prefixed_command_still_refused(cmd):
    # A builder may legitimately `git add`, so this list has no staging;
    # an env assignment must not hide a verb the builder may not run.
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, cmd), f"env-prefixed must refuse: {cmd}"


def test_env_prefixed_product_staging_refused_for_reviewer():
    # `git add` is builder-legal; a reviewer may only stage its
    # adversarial tests, never product code -- env-prefixed too.
    pats = policy.deny_patterns_for("reviewer", home=HOME)
    assert _refuses(pats, "X=y git add -A")
    assert _refuses(pats, "REV=1 git add delfin/agent/swarm/launch.py")


REDIRECTION_REFUSE = [
    "(patch -p1 < change.patch)",
    "(patch -p1 2>err) ",
    "{ patch -p1 < change.patch; } >& out",
    "exec 3>f; patch -p1 < change.patch",
]


@pytest.mark.parametrize("cmd", REDIRECTION_REFUSE)
def test_patch_after_redirection_is_refused(cmd):
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, cmd), f"redirected patch must refuse: {cmd}"


UNICODE_SEPARATOR_REFUSE = [
    "git push origin main\uFF1B ls",          # FULLWIDTH SEMICOLON
    "git push origin main\uFD3E ls\uFD3F",    # ORNATE brackets
    "cd /tmp\u00A0&&\u00A0pytest tests/x.py",  # NBSP joins
    "git pull origin main\u200b",              # zero-width space
]


@pytest.mark.parametrize("cmd", UNICODE_SEPARATOR_REFUSE)
def test_unicode_lookalike_separators_are_refused(cmd):
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, cmd), f"unicode separator must refuse: {cmd}"


HOME_TRAVERSAL_REFUSE = [
    "cat $HOME/../bob/secret.txt",
    "cat ${HOME}/../bob/secret.txt",
    "cat ~bob/file.txt",
]


@pytest.mark.parametrize("cmd", HOME_TRAVERSAL_REFUSE)
def test_home_traversal_to_another_user_is_refused(cmd):
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _refuses(pats, cmd), f"home traversal must refuse: {cmd}"


HOME_OWN_OK = [
    "cat $HOME/project/README.md",
    "cd $HOME",
    "cat ~/project/README.md",
    "cd ~",
]


@pytest.mark.parametrize("cmd", HOME_OWN_OK)
def test_own_home_reads_stay_allowed(cmd):
    pats = policy.deny_patterns_for("builder", home=HOME)
    assert _allows(pats, cmd), f"own home must allow: {cmd}"
