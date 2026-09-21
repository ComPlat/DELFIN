"""A cost of zero is three different statements. The header said one.

The engine already sorts every turn into measured, non-billing (a quota
or local hardware: a real zero) and unpriced (no rate, nothing
observed). The dashboard header printed ``$0.00`` for all three, so a
whole session on a model with no published rate read as a free one --
on a surface whose rule is that it must not say anything untrue.

And a push now names the branch that actually went up: `git push -u
origin HEAD` pushes the branch the worktree is on, which is not always
the one the model believes it made. A session reported `diag-denials`
while the remote had `session/8c29dfaa` (2026-09-17).
"""

from __future__ import annotations

import json

import pytest


class _Engine:
    def __init__(self, cost=0.0, unpriced=0, non_billing=0):
        self.cost_usd = cost
        self._unpriced_turns = unpriced
        self._non_billing_turns = non_billing


def _label(engine) -> str:
    from delfin.dashboard.tab_agent import _engine_cost_label
    return _engine_cost_label(engine)


def test_a_measured_cost_is_the_number():
    assert _label(_Engine(cost=1.2345)) == "$1.23"


def test_a_zero_nobody_could_price_says_so():
    assert _label(_Engine(cost=0.0, unpriced=4)) == "cost not measured"


def test_a_quota_run_says_what_it_is():
    assert _label(_Engine(cost=0.0, non_billing=4)) == "no charge"


def test_an_unpriced_turn_outranks_a_measured_zero():
    """One turn nobody could price makes the total unmeasured, whatever
    the others were."""
    assert _label(_Engine(cost=0.0, unpriced=1, non_billing=9)) \
        == "cost not measured"


def test_a_session_that_ran_nothing_is_still_zero():
    assert _label(_Engine(cost=0.0)) == "$0.00"


def test_an_engine_with_nothing_to_say_does_not_raise():
    class _Bare:
        pass

    assert _label(_Bare()) == "$0.00"


# -- the branch a push actually created -------------------------------------

def test_the_push_note_names_the_branch_that_went_up():
    from delfin.agent.api_client import _after_push

    result = json.dumps({
        "exit_code": 0,
        "stdout": "branch 'session/8c29dfaa' set up to track it.\n",
        "stderr": ("remote: \n"
                   "To github.com:ComPlat/DELFIN.git\n"
                   " * [new branch]        HEAD -> session/8c29dfaa\n"),
    })

    class _Perms:
        workspace = None
        push_grants = {"push": 1}

    note = _after_push({"command": "git push -u origin HEAD"}, result, _Perms())
    assert "session/8c29dfaa" in note
    assert "It went up as" in note


def test_a_push_that_did_not_land_says_nothing():
    from delfin.agent.api_client import _after_push

    result = json.dumps({"exit_code": 1, "stdout": "",
                         "stderr": "error: failed to push some refs"})

    class _Perms:
        workspace = None
        push_grants = {"push": 1}

    assert _after_push({"command": "git push"}, result, _Perms()) == ""
