"""A warning is justified by a decision, not by a word.

The grader tolerates a science-lint ``warning`` when the answer
"justifies" it. That test was: any word of more than three characters
from the finding's own message appearing anywhere in the answer.

Measured 2026-11-23 against the excited-state rungs, and each of these
passed the gate while deciding nothing:

    "This is the result."              -> justified tda-not-stated
                                          ("this" is in the message)
    "The results are listed here."     -> justified the CT warning
                                          ("here")
    the system prompt's ladder, quoted -> justified BOTH

The last one is the structural failure and the reason this file exists.
The ladder and the finding are rendered from one declaration -- which is
right, and is what keeps them from drifting -- so the grader was
checking for vocabulary it had handed the model itself. That is the
shape of "a rubric may score a vocabulary, not a behaviour".

Two routes replace it, and both are needed:

  * a rung whose decision lives in a SPEC FIELD is closed by the field.
    Prose cannot close it, so an agent cannot write "TDA was switched
    off" over a spec that is silent.
  * a rung where a counter-decision is legitimate in prose keeps prose,
    but only for what the PROMPT DID NOT SAY: the message's words minus
    the ladder's, which leaves exactly the case-specific interpolations
    -- the functional the spec chose, the state it declared.

A rule with no declared requirement keeps the word match with a floor
(two distinct words of five characters). That is an improvement and not
a solution; the asserted part is that it no longer accepts one short
connective.
"""
from __future__ import annotations

import pytest

from delfin.agent.workflow_grader import _lint_verdict
from delfin.tools._science_rules.excited import RUNGS, ladder_text
from delfin.tools.platform import scientific_lint

CT = "ct-needs-long-range-exchange"
TDA = "tda-not-stated"


def _spec(**extra):
    kwargs = {"charge": 0, "mult": 1, "method": "PBE0",
              "basis": "def2-TZVP", "excited_state": "charge_transfer"}
    kwargs.update(extra)
    return {"name": "m", "steps": [
        {"step": "orca_tddft", "label": "s", "kwargs": kwargs}]}


def _verdicts(spec, answer) -> dict:
    entries, _ = _lint_verdict(scientific_lint(spec), answer, spec)
    return {e["rule"]: e["justified"] for e in entries}


# --- the attacks that used to work ------------------------------------

@pytest.mark.parametrize("answer", [
    "This is the result.",
    "The results are listed here.",
    "Both functionals were used, from the same geometry.",
    "Done.",
])
def test_an_ordinary_report_sentence_justifies_nothing(answer):
    v = _verdicts(_spec(), answer)
    assert v.get(CT) is False, answer
    assert v.get(TDA) is False, answer


def test_the_ladder_quoted_back_justifies_nothing():
    """The first test the reviewer of this design asked for, and it
    failed before the change. The prompt hands the model this text; if
    repeating it satisfied the gate, the gate would measure reading."""
    v = _verdicts(_spec(), ladder_text())
    assert v.get(CT) is False
    assert v.get(TDA) is False


def test_prose_cannot_close_a_decision_that_lives_in_the_spec():
    """An agent must not be able to talk its way past a field. This is
    the one that matters most: the answer asserts the decision and the
    spec does not carry it."""
    v = _verdicts(_spec(), "TDA was switched off so the full response "
                           "is solved, as the ladder requires.")
    assert v.get(TDA) is False


# --- what must still pass ---------------------------------------------

def test_a_reasoned_counter_decision_is_accepted():
    """A global hybrid on a charge-transfer state is a legitimate choice
    when the answer says why. Refusing prose here would make the rung an
    error in disguise."""
    v = _verdicts(_spec(), "PBE0 is kept although the state has "
                           "charge-transfer character: the Lambda "
                           "diagnostic is 0.71, so the overlap is large.")
    assert v.get(CT) is True


@pytest.mark.parametrize("value", [True, False])
def test_the_field_closes_the_rung_whatever_the_answer_says(value):
    """Either value is a decision; the rung is about it being recorded.
    With the field set the lint does not fire at all, which is the
    cleanest form of closed."""
    v = _verdicts(_spec(tda=value), "Done.")
    assert TDA not in v


# --- the floor for rules with no declared requirement -----------------

def test_a_rule_without_a_requirement_needs_more_than_one_word():
    """Every other lint rule in the tree keeps the word match. What
    changed is the floor: one four-letter connective is no longer a
    justification."""
    finding = [{"step": "s", "level": "warning", "rule": "some-other-rule",
                "message": "this step does not state the solvation model"}]
    one_word, _ = _lint_verdict(finding, "This is it.", _spec())
    assert one_word[0]["justified"] is False
    two_words, _ = _lint_verdict(
        finding, "The solvation model is CPCM and the state is local.",
        _spec())
    assert two_words[0]["justified"] is True


def test_an_error_level_finding_is_never_justified():
    """Unchanged, and asserted because the routes above must not have
    opened a path around it."""
    finding = [{"step": "s", "level": "error", "rule": "x",
                "message": "electron count and multiplicity disagree"}]
    entries, reasons = _lint_verdict(
        finding, "electron count and multiplicity disagree", _spec())
    assert entries[0]["justified"] is False
    assert reasons


# --- the requirement is declared where the rung is --------------------

def test_every_rung_declares_how_it_may_be_justified():
    from delfin.tools._science import justification_requirement
    for rung in RUNGS:
        req = justification_requirement(rung.code)
        assert req, rung.code
        assert ("spec_field" in req) or ("prompt_vocabulary" in req), rung
        if rung.spec_field:
            assert req.get("spec_field") == rung.spec_field


def test_a_rule_nobody_declared_returns_no_requirement():
    from delfin.tools._science import justification_requirement
    assert justification_requirement("not-a-rule") is None


def test_the_reason_names_the_missing_field():
    """A reader of the verdict has to be able to act on it: "not
    justified" without the field name sends them to the message."""
    entries, reasons = _lint_verdict(
        scientific_lint(_spec()), "Done.", _spec())
    tda = next(e for e in entries if e["rule"] == TDA)
    assert "tda" in str(tda.get("justified_by", ""))
    assert any("tda-not-stated" in r for r in reasons)
