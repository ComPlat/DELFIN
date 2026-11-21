"""A reply with no verdict says WHICH way it failed.

Input:  the model's reply, as ``action_classifier.judge`` receives it.
Output: the verdict's reason, which must name an empty reply, a reply that
        stopped before any decision, and a reply with no clear verdict as
        three different things; and the shape description of each reply.
Semantics: the DECISION is unchanged in every case -- all three still end in
        the documented error behaviour (ask, or refuse where nobody can be
        asked), so this changes attribution only, never a security outcome.
        The alternative rejected was one message for all of them, which is
        what made 138 held and 16 down entries in one night's record
        unattributable (scratch/denials.txt:32-63).
"""

import pytest

from delfin.agent import action_classifier as ac


def _action(**kw):
    base = {
        "kind": "bash",
        "tool": "bash",
        "command": "pytest -q tests/",
        "path": "",
        "text": "",
        "raw": "",
        "cwd": "",
    }
    base.update(kw)
    return base


def _context(**kw):
    base = {"task": "Run the suite and report.", "recent_user": [],
            "recent_calls": []}
    base.update(kw)
    return base


def _judge(answer, **kw):
    return ac.judge(_action(), _context(), model="m",
                    runner=lambda prompt, model: answer, **kw)


# --- the shape description ------------------------------------------------

@pytest.mark.parametrize("answer,expected", [
    ("", "empty"),
    ("   \n\t\n", "empty"),
    (None, "empty"),
    # Cut before the decision line the prompt asks for.
    ("Thinking about this command, I would say that", "cut"),
    # The value was cut off.
    ("DECISION: al", "cut"),
    ("DECISION:", "cut"),
    # A decision line the model never filled in a verdict for.
    ("Sure! This looks fine to me.", "cut"),
    # Decide-able text the parser cannot read: not a cut.
    ("DECISION: maybe\nREASON: unsure.", "unreadable"),
])
def test_the_shape_of_a_reply_without_a_verdict(answer, expected):
    assert ac.no_verdict_shape(answer) == expected


# --- the operator-facing reason ------------------------------------------

@pytest.mark.parametrize("answer,needle", [
    ("", "empty reply"),
    ("   \n", "empty reply"),
    ("Thinking about this command, I would say that",
     "stopped before any decision"),
    ("DECISION: al", "stopped before any decision"),
    ("DECISION: maybe\nREASON: unsure.", "no clear verdict"),
])
def test_the_reason_names_which_failure_it_was(answer, needle):
    verdict = _judge(answer)
    assert verdict.reason.startswith(ac.ERROR_PREFIX)
    assert needle in verdict.reason


def test_the_three_shapes_are_three_different_reasons():
    empty = _judge("").reason
    cut = _judge("Thinking about this command").reason
    unreadable = _judge("DECISION: maybe\nREASON: unsure.").reason
    assert len({empty, cut, unreadable}) == 3


# --- never-worse: no decision changes ------------------------------------

@pytest.mark.parametrize("answer,decision", [
    ("DECISION: allow\nREASON: ordinary step.", "allow"),
    ("DECISION: ask\nREASON: unclear.", "ask"),
    ("DECISION: refuse\nREASON: destructive.", "refuse"),
])
def test_a_readable_verdict_is_untouched(answer, decision):
    verdict = _judge(answer)
    assert verdict.decision == decision
    assert not verdict.reason.startswith(ac.ERROR_PREFIX)


@pytest.mark.parametrize("answer", [
    "", "   ", "Thinking", "DECISION: al", "Sure! Fine.",
    "DECISION: maybe\nREASON: unsure.", None,
])
@pytest.mark.parametrize("default_on_error", ["ask", "refuse", "allow"])
def test_the_documented_error_behaviour_still_applies(answer, default_on_error):
    """Whatever the shape, the decision is the caller's default -- the fix
    changes the reason, never the outcome."""
    assert _judge(answer, default_on_error=default_on_error).decision \
        == default_on_error


def test_the_strictest_decision_still_wins():
    """An action text that echoes an allow next to the model's own refuse
    must not win by coming first (the property _parse_answer protects)."""
    verdict = _judge("DECISION: allow\nREASON: echoed\n"
                     "DECISION: refuse\nREASON: destructive.")
    assert verdict.decision == "refuse"
