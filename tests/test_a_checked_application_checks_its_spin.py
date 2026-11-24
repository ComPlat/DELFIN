"""The route called "checked" did not check what the rule says a check is.

The rule (2026-10-22): a calculation runs without a dialog when DELFIN
has checked that it is set up correctly -- "charge/spin parity, a valid
input, real resource figures" -- and its resources are moderate.
Otherwise the dialog states the reason, and with nobody to ask it is
refused. One policy for every route.

MEASURED at the real MCP entry, default mode, before this change:

    conformer_dG  charge=0 mult=2   admitted, no dialog
    conformer_dG  charge=1 mult=1   admitted, no dialog
    conformer_dG  no inputs at all  admitted, no dialog

Both of the first two are parity violations for the molecule a caller
would pass. The application judge read charge and mult from the call
but never the MOLECULE, so it could not count electrons, and it admitted
an unread setup while its own docstring said "An unread check is never
a pass; that rule is the whole reason this module is trusted".

Every other route already did this: the CONTROL.txt path counts the
geometry's electrons and refuses a multiplicity they cannot have, and
refuses when the geometry cannot be read at all. So the fix is the
existing check applied here -- ``electrons_of`` and ``spin_allowed``,
the same functions -- not a second one.

AND ONE THING IT MUST NOT DO, guarded below because it would have broken
the core of DELFIN: refuse a declared charge that differs from the
SMILES' own. A redox workflow passes the NEUTRAL molecule's SMILES with
charge=+1, mult=2 for its radical cation. Electrons are the nuclear
charge minus the DECLARED charge, as on the bash route, and that case is
admitted.
"""
from __future__ import annotations

import ast
import inspect
import tempfile

import pytest

from delfin.agent import api_client as ac
from delfin.agent import calc_policy as cp

APP = "conformer_dG"          # takes a {mult}
NO_MULT_APP = "redox_potential"   # the application chooses the spin
ETHANOL = "CCO"               # 26 electrons when neutral


def _v(**inputs):
    return cp.judge_application_run(APP, 1, inputs=inputs or None)


# --- the parity the rule names -----------------------------------------

@pytest.mark.parametrize("charge,mult", [(0, 2), (1, 1), (0, 4), (-1, 1)])
def test_a_multiplicity_the_electrons_cannot_have_is_not_admitted(charge,
                                                                  mult):
    v = _v(smiles=ETHANOL, charge=charge, mult=mult)
    assert v is not None and v.calculation
    assert not v.admit
    assert "multiplicity" in v.reason and "electrons" in v.reason


@pytest.mark.parametrize("charge,mult", [(0, 1), (0, 3), (1, 2), (-1, 2)])
def test_a_multiplicity_the_electrons_can_have_is_admitted(charge, mult):
    v = _v(smiles=ETHANOL, charge=charge, mult=mult)
    assert v is not None and v.admit, v.reason
    assert "electrons" in v.summary


def test_the_radical_cation_of_a_neutral_smiles_is_admitted():
    """The redox case. A check on 'SMILES charge == declared charge'
    would refuse this and break every redox workflow."""
    v = _v(smiles=ETHANOL, charge=1, mult=2)
    assert v is not None and v.admit, v.reason


def test_a_charged_smiles_counts_its_own_hydrogens():
    """[NH4+]: 11 nuclear charge, declared +1 -> 10 electrons, singlet."""
    ok = cp.judge_application_run(APP, 1, inputs={
        "smiles": "[NH4+]", "charge": 1, "mult": 1})
    bad = cp.judge_application_run(APP, 1, inputs={
        "smiles": "[NH4+]", "charge": 1, "mult": 2})
    assert ok.admit and not bad.admit


# --- an unread check is never a pass -----------------------------------

@pytest.mark.parametrize("inputs,why", [
    ({}, "molecule"),
    ({"charge": 0, "mult": 1}, "molecule"),
    ({"smiles": ETHANOL, "mult": 1}, "charge"),
    ({"smiles": ETHANOL, "charge": 0}, "multiplicity"),
    ({"smiles": "C(C(", "charge": 0, "mult": 1}, "molecule"),
])
def test_what_cannot_be_read_is_not_admitted(inputs, why):
    v = cp.judge_application_run(APP, 1, inputs=inputs or None)
    assert v is not None and not v.admit
    assert why in v.reason, v.reason


def test_an_application_that_chooses_its_own_spin_needs_no_mult():
    """redox_potential has no {mult}: the multiplicity is not an input,
    so its absence is not an unread check -- and the summary must not
    claim a parity it did not test."""
    v = cp.judge_application_run(NO_MULT_APP, 1, inputs={
        "smiles": ETHANOL, "charge": 0})
    assert v is not None and v.admit, v.reason
    assert "chosen by the application" in v.summary


def test_a_charge_beyond_the_nuclear_charge_is_not_admitted():
    v = cp.judge_application_run(NO_MULT_APP, 1, inputs={
        "smiles": "[H][H]", "charge": 3})
    assert v is not None and not v.admit


# --- one policy: the same functions as the other routes ----------------

def test_the_judge_uses_the_house_parity_functions():
    tree = ast.parse(inspect.getsource(cp.judge_application_run))
    called = {n.func.id for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert {"electrons_of", "spin_allowed"} <= called, called


# --- and at the entry that counts ---------------------------------------

def _entry(args):
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    perms.confirm_callback = None
    return ac._mcp_call_refusal("mcp__delfin-tools__run_application",
                                args, perms)


def test_unattended_a_parity_violation_is_refused_with_its_reason():
    err = _entry({"name": APP, "cores": 1, "inputs": {
        "smiles": ETHANOL, "charge": 0, "mult": 2}})
    assert err is not None and "multiplicity" in err


def test_unattended_a_correct_setup_still_runs():
    """The gradient this route was opened for must survive the check."""
    err = _entry({"name": APP, "cores": 1, "inputs": {
        "smiles": ETHANOL, "charge": 0, "mult": 1}})
    assert err is None, err
