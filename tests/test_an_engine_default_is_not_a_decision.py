"""An excited-state result picked by a default was not decided.

The first rungs of the TD-DFT decision ladder. Both say the same thing
in two places: a spec that leaves a choice to the program records no
choice, and a level of theory nobody chose cannot be defended by the
person who publishes the number.

They are science-lint rules rather than a new mechanism, because the
enforcement already exists: ``agent/workflow_grader`` tolerates a
``warning`` only when the answer text justifies it. So a step that
decided deliberately and says why passes, and one that drifted into a
default does not — which is the pressure the ladder is for.

Why these two rungs and not a longer list: an outside measurement on
this task class (558 QUESTDB vertical transition energies) found an
unguided agent at 24.9 % coverage and the same agent with a structured
decision ladder at 44.1 %. The lever is structure, so the rungs have to
be ones a run can actually act on, and each is asserted here against
its own counter-case.

Input: a spec, through the single entry point ``scientific_lint``.
Output: the findings. Every assertion goes through that entry point and
not through the rule module, because a rule that fires only when called
directly is not in the product.
"""
from __future__ import annotations

import pytest

from delfin.tools.platform import scientific_lint

CT = "ct-needs-long-range-exchange"
TDA = "tda-not-stated"


def _rules(spec) -> list[str]:
    out = []
    for f in scientific_lint(spec):
        rule = f.get("rule") if isinstance(f, dict) else getattr(f, "rule", "")
        if rule in (CT, TDA):
            out.append(rule)
    return out


def _spec(step="orca_tddft", **kwargs):
    base = {"charge": 0, "mult": 1, "method": "PBE0", "basis": "def2-TZVP"}
    base.update(kwargs)
    return {"steps": [{"step": step, "label": "s", "kwargs": base}]}


# --- the positive control comes first ----------------------------------

def test_both_rungs_fire_on_the_case_they_were_written_for():
    """A charge-transfer state on a global hybrid, silent about TDA.

    First, and deliberately: four quiet cases below prove nothing unless
    one case is known to fire. A probe written for these rules reported
    "quiet" everywhere once because it read the wrong field name, and
    that looked exactly like a clean result.
    """
    assert sorted(_rules(_spec(excited_state="charge_transfer"))) == [CT, TDA]


# --- rung 1: charge transfer needs non-local exchange -------------------

@pytest.mark.parametrize("method", ["wB97X-D4", "CAM-B3LYP", "LC-BLYP",
                                    "wb97x-d3"])
def test_a_range_separated_functional_closes_the_first_rung(method):
    """The decision was made; the rule has nothing to say."""
    assert CT not in _rules(
        _spec(method=method, excited_state="charge_transfer", tda=True))


def test_a_state_with_no_declared_character_is_not_second_guessed():
    """The rule reads a DECLARED state, as the P1 Rydberg rule does. It
    does not infer charge-transfer character from a molecule, because
    nothing here can do that reliably and a wrong warning on every
    excited-state run would be worse than none."""
    assert CT not in _rules(_spec(tda=True))


@pytest.mark.parametrize("state", ["octahedral_field", "spectrum",
                                   "contracted"])
def test_a_word_that_merely_contains_ct_does_not_fire(state):
    """`ct` is matched as a whole token. Without that, "octahedral" and
    "spectrum" both claim charge-transfer character."""
    assert CT not in _rules(_spec(excited_state=state, tda=True))


@pytest.mark.parametrize("state", ["charge_transfer", "charge-transfer",
                                   "ct", "push_pull", "donor-acceptor",
                                   "MLCT", "intramolecular_ct"])
def test_the_spellings_a_spec_actually_uses_are_recognised(state):
    assert CT in _rules(_spec(excited_state=state, tda=True))


def test_a_placeholder_method_is_not_judged():
    """An unresolved spec is not a wrong spec."""
    from delfin.tools._science import _is_placeholder
    holder = next((v for v in ("<method>", "TBD", "?", "")
                   if _is_placeholder(v)), None)
    if holder is None:
        pytest.skip("no placeholder spelling to test with")
    assert CT not in _rules(
        _spec(method=holder, excited_state="charge_transfer", tda=True))


# --- rung 2: TDA is stated, not inherited ------------------------------

@pytest.mark.parametrize("value", [True, False])
def test_stating_tda_either_way_closes_the_second_rung(value):
    """The rung is about the decision being recorded, not about which
    way it went: TDA is a legitimate choice, and so is the full
    response. What is not legitimate is the record not saying."""
    assert TDA not in _rules(_spec(tda=value))


def test_a_silent_step_is_warned_because_orca_applies_tda_by_default():
    assert TDA in _rules(_spec())


# --- neither rung reaches past the steps it is for ---------------------

@pytest.mark.parametrize("step", ["orca_sp", "orca_opt", "orca_freq",
                                  "orca_opt_freq", "crest_conformers"])
def test_a_ground_state_step_is_left_alone(step):
    """These rules are about a TD-DFT response equation. A ground-state
    step has no TDA and no excited state, and warning about it would
    teach a reader to ignore the lint."""
    assert _rules(_spec(step=step, excited_state="charge_transfer")) == []


@pytest.mark.parametrize("step", ["orca_tddft", "esd", "esd_opi"])
def test_every_excited_state_step_is_covered(step):
    """The three steps that solve the response equation, named rather
    than pattern-matched -- so a new one is added deliberately."""
    assert TDA in _rules(_spec(step=step))


# --- the grader is what turns a warning into pressure -----------------

def test_the_rungs_are_warnings_so_the_answer_must_justify_them():
    """Not errors. A warning the answer justifies is tolerated by
    workflow_grader, which is exactly the behaviour wanted: a deliberate
    global hybrid with a stated reason passes, an unexamined one does
    not. An error would refuse the deliberate case too.
    """
    levels = {
        (f.get("rule") if isinstance(f, dict) else f.rule):
        (f.get("level") if isinstance(f, dict) else f.level)
        for f in scientific_lint(_spec(excited_state="charge_transfer"))}
    assert levels.get(CT) == "warning"
    assert levels.get(TDA) == "warning"


def test_every_finding_names_the_step_it_is_about():
    """A finding without a step cannot be acted on in a spec with four
    of them."""
    for f in scientific_lint(_spec(excited_state="ct")):
        step = f.get("step") if isinstance(f, dict) else f.step
        assert step == "s"
