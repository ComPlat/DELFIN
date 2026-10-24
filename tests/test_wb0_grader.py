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
