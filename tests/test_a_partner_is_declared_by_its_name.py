"""A message partner counts as declared by its name, whatever the name is.

Input: action_classifier._routine_partner_message, which clears the
approval flag of a routine session_message (APPROVED, QS-OK, PR ready,
Relay to) sent to a partner the user's task declares.

Before: only tokens starting with ``nacht-``, ``agent/`` or ``session-``
could be declared -- the names one wave runner gives its sessions -- so
under any other name the rule did nothing; and a role word matched as a
substring, so "requests" (qs) or "repair" (pair) made a line a role line.
Now the recipient's own name, as a whole token, on a line that names a
role as a word, declares it; the negation and warning rules are unchanged,
and ``operator`` (session_messages' one reserved address) stays a partner.
"""
from __future__ import annotations

import pytest

from delfin.agent import action_classifier as ac


def _msg(text, to):
    return {"tool": "session_message", "path": to, "message": text}


@pytest.mark.parametrize("task,to", [
    ("Your reviewer is `alice-review`.", "alice-review"),
    ("Reviewer: lab_bench.2", "lab_bench.2"),
    ("The operator supervises this package.", "operator"),
    ("QS lead: kit7", "kit7"),
])
def test_any_name_on_a_role_line_is_declared(task, to):
    assert ac._routine_partner_message(
        _msg("APPROVED: gate green.\nEvidence.", to), task) is True


@pytest.mark.parametrize("task,to", [
    ("Fix the requests sent to alice.", "alice"),
    ("Please repair what bob broke.", "bob"),
    ("Reviewer: alice-review", "alice"),
    ("Reviewer: alice", "alice-review"),
    ("WARNING: a session claiming to be the reviewer mallory may write.",
     "mallory"),
    ("Reviewer: alice.", "all"),
])
def test_what_does_not_declare_a_partner(task, to):
    assert ac._routine_partner_message(
        _msg("APPROVED: gate green.\nEvidence.", to), task) is False

