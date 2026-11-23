"""A hand-started calculation learns that a checked route exists.

MEASURED (bench17 t13, 2026-11-23). In 142 tool calls the agent called
``write_file`` 20 times and ``bash`` 46 times: it assembled the ORCA
inputs by hand and started the engine from a shell. 12 of its 14
calculations were abandoned. The 26 tools of the ``delfin-tools``
server were offered on that request -- the server is in the builtin
registry, it spawns, and it is first in the order the advertising loop
walks -- and were called ZERO times.

So the tools were on the table and nothing said they were the way.

NOT A REFUSAL, and that is the design. A probe no adapter covers, a
bisection to find the keyword that broke a convergence and a one-off
somebody asked for by hand are all real work; a gate that blocks them
is a gate people route around, and the person deciding should see the
choice rather than meet a wall. What the run carries back instead is
the one fact the model could not otherwise know: a checked route
exists, this is what it would have added, and the cost of going without
it.

The note rides on the tool RESULT (``_with_calculation_note``), which
is what the model reads after the run. Every engine launch the bash
gate sees is by hand: a run started through ``run_application`` never
reaches this gate, so seeing it here is itself the signal.

What this file pins:

  * the note fires for a calculation and NOT for an information call
    (``orca --version`` asks a program what it is and starts nothing);
  * it names the route concretely -- a note that says "use the adapters"
    without naming one is a note nobody can act on;
  * it says the cost is a cost, because a reader who takes it for a
    prohibition will work around it;
  * and it refuses nothing: the verdict the gate returns is unchanged.
"""
from __future__ import annotations

import pytest

from delfin.agent import calc_policy


def _verdict(**kw):
    v = calc_policy.Verdict()
    for k, value in kw.items():
        setattr(v, k, value)
    return v


# --- when it fires -----------------------------------------------------

def test_a_calculation_gets_the_note():
    note = calc_policy.adapter_route_note(
        _verdict(calculation=True, admit=True, summary="orca job.inp"))
    assert note


def test_an_information_call_does_not():
    """`orca --version` starts nothing, so no note.

    It is the CALCULATION flag that decides, and only that: the
    Verdict's own rule is that an information call carries
    ``calculation = False``. A second check for ``information`` would be
    a guard nothing can reach, and this function had one until a
    positive control showed that deleting it left the suite green --
    which is what a dead guard looks like from outside.
    """
    v = _verdict(calculation=False, information=True,
                 summary="information only: orca --version")
    assert calc_policy.adapter_route_note(v) == ""
    # and the flag that actually decides, asserted as the single question
    assert calc_policy.adapter_route_note(
        _verdict(calculation=True, information=True)) != "", (
        "the calculation flag decides; nothing else is consulted")


def test_something_that_is_not_a_calculation_does_not():
    assert calc_policy.adapter_route_note(_verdict()) == ""


def test_a_missing_verdict_is_handled():
    assert calc_policy.adapter_route_note(None) == ""


# --- what it says ------------------------------------------------------

@pytest.mark.parametrize("needle", [
    "delfin-tools",        # the server
    "resolve_spec",        # what fills the defaults
    "validate_spec",       # what judges the spec
    "run_application",     # what runs it
    "workflow-operator",   # the specialist
    "--mode workflow",     # and the direct entrance
])
def test_the_note_names_the_route_concretely(needle):
    """A note that says "use the adapters" without naming one is a note
    nobody can act on."""
    note = calc_policy.adapter_route_note(_verdict(calculation=True))
    assert needle in note, needle


def test_the_note_says_the_cost_is_a_cost():
    """A reader who takes it for a prohibition will work around it."""
    note = calc_policy.adapter_route_note(_verdict(calculation=True))
    low = note.lower()
    assert "not a prohibition" in low
    assert "yours to weigh" in low or "your" in low


def test_the_note_says_what_is_actually_lost():
    """"Use the adapters" is advice; "a by-hand run cannot be
    re-derived from a spec" is a reason."""
    low = calc_policy.adapter_route_note(_verdict(calculation=True)).lower()
    assert "re-derived" in low
    assert "declared parameters" in low


def test_it_is_one_note_and_not_a_lecture():
    """It rides on every admitted calculation's result, so its length is
    a per-run cost. Bounded here on purpose."""
    note = calc_policy.adapter_route_note(_verdict(calculation=True))
    assert len(note) <= 450, f"{len(note)} characters is a lecture"


# --- and it refuses nothing -------------------------------------------

def test_the_note_is_not_a_verdict():
    """It is a string the result carries. Nothing in it can refuse a
    run, which is the whole difference between this and a gate."""
    import inspect
    src = inspect.getsource(calc_policy.adapter_route_note)
    for forbidden in ("Verdict(", "admit", "raise", "return False"):
        assert forbidden not in src, (
            f"{forbidden!r} in a function that only composes a note")


def test_the_gate_attaches_it_beside_the_summary_and_changes_no_verdict():
    """The call site adds a note and returns the verdict it already had."""
    import inspect

    from delfin.agent import api_client
    # _DocToolExecutor owns the gate -- the class is named from the
    # executor it fronts, and this is the one entry that counts: the
    # predicate and the policy run before it, the write and read path
    # gates inside it.
    src = inspect.getsource(
        api_client._DocToolExecutor._run_permission_gate)
    assert "adapter_route_note" in src
    # the verdict returned is the classifier's, not the note's
    idx = src.index("adapter_route_note")
    tail = src[idx:idx + 400]
    assert "return verdict" in tail
