"""WB0 phase 1 red control: the workflow grader grades a built SPEC, not prose.

The wave-21 goal: agents build their OWN physically realistic workflows and
the platform grades the spec they produced BEFORE anything runs. This file
pins the grading contract (delfin/agent/workflow_grader.py, new in this
package) against the physically wrong specs a careless builder produces,
each of which must fail grading, and against a correct spec, which must
pass.

Every graded spec uses the real adapter step names
(delfin/tools/adapters/*.py) and goes through the real call sites:
platform.resolve_spec / platform.validate_spec (delfin/tools/platform.py
:256 / :368) and the science lint (_science.scientific_lint via platform.py
:239). No shadow copy of anything is imported.
"""

from __future__ import annotations

import pytest

from delfin.agent import workflow_grader

# A correct gas-phase opt+freq of a small closed-shell molecule, one level
# of theory throughout — the answer justifies nothing because nothing needs
# justification.
GOOD_OPT_FREQ = {
    "name": "methanol_opt_freq",
    "steps": [
        {"step": "smiles_to_xyz", "kwargs": {"smiles": "CO"}},
        {"step": "orca_opt_freq",
         "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                    "basis": "def2-TZVP"}},
    ],
}

# A request a chemist would write; the grader checks that the spec covers
# what the request asks for (here: thermochemistry of the minimum).
GOOD_REQUEST = (
    "Optimise methanol and confirm the geometry is a minimum; I need the "
    "thermochemical correction at 298 K."
)


def grade(spec, request, answer=None):
    return workflow_grader.grade(spec, request, answer=answer)


# ---------------------------------------------------------------------------
# The correct spec passes
# ---------------------------------------------------------------------------

def test_a_correct_spec_passes():
    verdict = grade(GOOD_OPT_FREQ, GOOD_REQUEST)
    assert verdict["status"] == "pass", verdict["reasons"]


def test_the_verdict_is_structured():
    verdict = grade(GOOD_OPT_FREQ, GOOD_REQUEST)
    for key in ("status", "reasons", "findings"):
        assert key in verdict
    for key in ("resolve", "validate", "lint", "coverage", "dry_run", "cost"):
        assert key in verdict["findings"]


# ---------------------------------------------------------------------------
# Physically wrong specs fail — five designed defects
# ---------------------------------------------------------------------------

def test_an_energy_difference_across_two_levels_of_theory_is_refused():
    """DeltaE between two members of a couple must come from ONE method;
    mixing levels makes the difference meaningless."""
    spec = {
        "name": "redox_mixed",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": 1, "mult": 2,
                                           "method": "PBE0",
                                           "basis": "def2-SVP"}},
            {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                           "method": "B3LYP",
                                           "basis": "def2-TZVP"}},
        ],
    }
    verdict = grade(spec, "Compute the ionisation potential of the couple.")
    assert verdict["status"] == "fail"


def test_a_ts_without_frequency_check_or_irc_is_refused():
    spec = {
        "name": "ts_nocheck",
        "steps": [
            {"step": "orca_ts_opt", "kwargs": {"charge": 0, "mult": 1}},
        ],
    }
    verdict = grade(spec, "Find the TS and verify it.")
    assert verdict["status"] == "fail"


def test_an_anion_without_diffuse_functions_is_refused():
    spec = {
        "name": "anion_plain",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": -1, "mult": 1,
                                           "method": "PBE0",
                                           "basis": "def2-SVP"}},
        ],
    }
    verdict = grade(spec, "Energise the anion.")
    assert verdict["status"] == "fail"


def test_a_5d_metal_without_a_relativistic_hamiltonian_is_refused():
    spec = {
        "name": "gold_plain",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                           "method": "PBE0",
                                           "basis": "def2-SVP"}},
        ],
    }
    verdict = grade(spec, "A gold-containing complex: give me its energy.")
    assert verdict["status"] == "fail"


def test_a_flexible_molecule_straight_to_dft_without_a_conformer_step():
    spec = {
        "name": "flexible_noconf",
        "steps": [
            {"step": "orca_opt_freq",
             "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                        "basis": "def2-TZVP"}},
        ],
    }
    request = (
        "Ethanol is flexible. Find its lowest-energy conformer and then the "
        "thermochemistry of that conformer at one level of theory."
    )
    verdict = grade(spec, request)
    assert verdict["status"] == "fail"


# ---------------------------------------------------------------------------
# Justified warnings: the answer text carries the reason
# ---------------------------------------------------------------------------

def test_an_unjustified_lint_warning_fails_when_answer_is_graded():
    """A lint warning is a defensible choice ONLY with a stated reason.
    Without an answer (or with a non-answer) the warning is unjustified."""
    spec = {
        "name": "composite_redundant",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                           "method": "r2scan-3c",
                                           "basis": "def2-TZVP"}},
        ],
    }
    verdict = grade(spec, "Energy of methanol.", answer="done")
    assert verdict["status"] == "fail"


def test_a_justified_lint_warning_passes():
    spec = {
        "name": "composite_redundant",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                           "method": "r2scan-3c",
                                           "basis": "def2-TZVP"}},
        ],
    }
    verdict = grade(spec, "Energy of methanol.",
                    answer="I set basis def2-TZVP because the composite "
                           "default is too small for the diffuse tail I "
                           "need; I know it is ignored by r2scan-3c.")
    assert verdict["status"] == "pass"


def test_the_justification_branch_is_graded_per_finding():
    """s25 approved-package note: the 'justified' field of each lint
    finding must flip exactly with the answer text — with no answer the
    warning is unjustified and fails; with the substance named it passes;
    the same answer for a DIFFERENT unjustified deviation does not clear
    it."""
    spec = {
        "name": "composite_redundant",
        "steps": [
            {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                           "method": "r2scan-3c",
                                           "basis": "def2-TZVP"}},
        ],
    }
    no_answer = grade(spec, "Energy of methanol.", answer="")
    assert no_answer["status"] == "fail"
    assert no_answer["findings"]["lint"][0]["justified"] is False

    justified = grade(spec, "Energy of methanol.",
                      answer="The explicit basis def2-TZVP is ignored by "
                             "the composite method; I know this and I keep "
                             "it only as documentation for the metal_basis "
                             "override I plan later.")
    assert justified["status"] == "pass"
    assert justified["findings"]["lint"][0]["justified"] is True

    # the justification names the WRONG deviation: still unjustified
    wrong = grade(spec, "Energy of methanol.",
                  answer="I chose PBE0 because it is semi-local and fast.")
    assert wrong["status"] == "fail"
    assert wrong["findings"]["lint"][0]["justified"] is False


# ---------------------------------------------------------------------------
# Broken input is refused, never half-graded
# ---------------------------------------------------------------------------

def test_a_spec_that_does_not_resolve_fails():
    verdict = grade({"name": "bad", "steps": [{"step": "orca_opt"}]},
                    "Optimise something.")
    assert verdict["status"] == "fail"


def test_a_spec_without_a_name_fails():
    with pytest.raises(ValueError):
        grade({"steps": []}, "anything")


def test_an_unknown_step_name_fails():
    verdict = grade({"name": "m", "steps": [{"step": "classic_opt_freq"}]},
                    "Optimise methanol.")
    assert verdict["status"] == "fail"


# ---------------------------------------------------------------------------
# Reviewer s25 findings (A1-A4, B5) — false-positive physics demands
# ---------------------------------------------------------------------------

def test_reddox_request_does_not_demand_rhenium_relativity():
    """A1: the word 'redox' contains 're' — a word prefix must never
    trigger the 5d rule."""
    spec = {"name": "copper_redox", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 1, "mult": 2,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "redox potential of a copper complex")
    assert verdict["status"] == "pass", verdict["reasons"]


def test_ir_spectrum_does_not_demand_iridium_relativity():
    """A2: 'IR spectrum' contains 'ir'. The spec carries the IR parse step
    the request needs — only the false relativity demand must go."""
    spec = {"name": "ir_co", "steps": [
        {"step": "orca_opt_freq", "kwargs": {"charge": 0, "mult": 1,
                                             "method": "PBE0",
                                             "basis": "def2-TZVP"}},
        {"step": "ir_parse", "kwargs": {}}]}
    verdict = grade(spec, "the IR spectrum of carbon monoxide")
    assert verdict["status"] == "pass", verdict["reasons"]


def test_water_does_not_demand_tungsten_relativity():
    """A3: 'water' contains 'w'."""
    spec = {"name": "water_sp", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "the energy of water")
    assert verdict["status"] == "pass", verdict["reasons"]


def test_hf_symbol_case_insensitive_still_caught():
    """A4 positive control: a real 5d mention still fails without
    relativity (the symbol, upper case, as a chemist writes it)."""
    spec = {"name": "au_sp", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "the energy of the Au complex")
    assert verdict["status"] == "fail"


# ---------------------------------------------------------------------------
# QS finding F1 (s29, seconded s25): the relativity DECISION must come
# from DELFIN's own policy — delfin.utils._should_use_rel (:203-217,
# '4d5d'/'4f5f'/'mixed') — never from a grader-local element list. The
# grader's table is only the TEXT->SYMBOL resolver.
# ---------------------------------------------------------------------------

def test_a_4d_metal_is_caught_without_relativity():
    spec = {"name": "ru_sp", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "the energy of the Ru complex")
    assert verdict["status"] == "fail"


def test_a_lanthanide_is_caught_without_relativity():
    spec = {"name": "gd_sp", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "the energy of the gadolinium complex")
    assert verdict["status"] == "fail"


def test_pure_3d_iron_is_not_caught():
    """_should_use_rel returns False for a pure 3d system — the grader
    must agree (the spin task grades spin handling, not relativity)."""
    spec = {"name": "fe_sp", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 5,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    verdict = grade(spec, "the energy of the iron complex")
    assert verdict["status"] == "pass", verdict["reasons"]


def test_grader_and_policy_agree_on_the_pin_set():
    """The grader's verdict and DELFIN's _should_use_rel must agree for
    every element of the pin set {Fe, Au, Os, Gd, Rn}."""
    from delfin.utils import _should_use_rel

    spec = {"name": "pin", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0, "mult": 1,
                                       "method": "PBE0",
                                       "basis": "def2-TZVP"}}]}
    for request, metals in (
            ("the energy of the iron complex", ["Fe"]),
            ("the energy of the Au complex", ["Au"]),
            ("the energy of the osmium complex", ["Os"]),
            ("the energy of the gadolinium complex", ["Gd"]),
            ("the energy of the radon compound", ["Rn"])):
        verdict = grade(spec, request)
        needs = _should_use_rel(metals)
        rel_reasons = [r for r in verdict["reasons"] if "relativistic" in r]
        if needs:
            assert verdict["status"] == "fail", (request, verdict["reasons"])
            assert len(rel_reasons) == 1, (request, verdict["reasons"])
        else:
            # For not-needs, relativity is only one of several checks:
            # assert the specific reason is absent, not a blanket pass.
            assert not rel_reasons, (request, verdict["reasons"])


def test_rigid_molecule_frequency_request_does_not_demand_conformers():
    """B5: 'frequencies of the fixed chair conformation' mentions the
    word conformation — a rigid-molecule request must pass."""
    spec = {"name": "chair_freq", "steps": [
        {"step": "orca_opt_freq", "kwargs": {"charge": 0, "mult": 1,
                                             "method": "PBE0",
                                             "basis": "def2-TZVP"}}]}
    verdict = grade(
        spec, "frequencies of the fixed chair conformation of cyclohexane")
    assert verdict["status"] == "pass", verdict["reasons"]
