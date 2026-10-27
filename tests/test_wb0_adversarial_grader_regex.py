"""WB0 adversarial tests for the grader (reviewer s25, package WB0).

Red tests here are FINDINGS against delfin/agent/workflow_grader.py
@860efa48; green tests pin behaviour that must not regress. Every spec
uses adapter step names and no template names.
"""
import pytest

from delfin.agent.workflow_grader import grade

GOOD = {"method": "pbe0", "basis": "def2-tzvp"}


def _spec(name, steps, **extra):
    return {"name": name, "steps": steps, **extra}


def _sp(label, charge=0, **kw):
    kwargs = dict(GOOD)
    kwargs.update(kw)
    if charge:
        kwargs["charge"] = charge
    return {"step": "orca_sp", "label": label, "kwargs": kwargs}


# --- A-group: false positives of the 5d/relativistic rule -----------------

def test_a1_redox_request_not_flagged_relativistic():
    # "redox" starts with "re"; the regex \b{sym}|{name}\b in
    # _relativistic_reason matches \bre in "redox" and demands a
    # relativistic Hamiltonian for a copper redox couple.
    spec = _spec("redox", [_sp("red", charge=-1, basis="ma-def2-tzvp"),
                           _sp("ox", charge=0)])
    verdict = grade(spec, "redox potential of a copper complex in water",
                    budget_core_hours=4.0)
    assert not any("relativistic" in r for r in verdict["reasons"]), verdict


def test_a2_ir_request_not_flagged_relativistic():
    # "IR" is the word "ir": \bir matches it, iridium is then demanded.
    spec = _spec("ir", [_sp("ir", charge=0)])
    verdict = grade(spec, "IR spectrum of benzene",
                    budget_core_hours=4.0)
    assert not any("relativistic" in r for r in verdict["reasons"]), verdict


def test_a3_water_request_not_flagged_relativistic():
    # "water" starts with "w" (tungsten's symbol).
    spec = _spec("aq", [_sp("aq", charge=0)])
    verdict = grade(spec, "solvation of glycine in water",
                    budget_core_hours=4.0)
    assert not any("relativistic" in r for r in verdict["reasons"]), verdict


def test_a4_hf_method_not_flagged_relativistic():
    # "HF" (Hartree-Fock) starts with "hf" (hafnium's symbol).
    spec = _spec("hf", [{"step": "orca_sp", "label": "sp",
                         "kwargs": {"method": "hf", "basis": "def2-tzvp"}}])
    verdict = grade(spec, "HF single point of water",
                    budget_core_hours=4.0)
    assert not any("relativistic" in r for r in verdict["reasons"]), verdict


def test_a5_gold_request_without_relativity_still_caught():
    # Positive control for the A-group: the rule must still FIRE when a 5d
    # metal really is named (fires-when-it-should side).
    spec = _spec("au", [_sp("au", charge=0)])
    verdict = grade(spec, "geometry of a gold chloride complex",
                    budget_core_hours=4.0)
    assert any("relativistic" in r for r in verdict["reasons"]), verdict


def test_a6_gold_formula_in_request_still_caught():
    spec = _spec("au2", [_sp("au2", charge=0)])
    verdict = grade(spec, "Au(I) chloride geometry, 20 atoms",
                    budget_core_hours=4.0)
    assert any("relativistic" in r for r in verdict["reasons"]), verdict


# --- B-group: mixed-level rule false positives -----------------------------

def test_b1_xtb_plus_orca_is_a_funnel_not_a_mixed_level():
    # A cheap semi-empirical step followed by a DFT step is standard
    # practice (crest_conformers -> orca_opt): xtb steps carry no method,
    # so they cannot create a second "level". _level_of only reads
    # orca-prefixed steps — pinning that xtb steps are correctly ignored.
    spec = _spec("funnel", [{"step": "crest_conformers", "label": "conf"},
                            _sp("sp")])
    verdict = grade(spec, "optimize structures of propane",
                    budget_core_hours=4.0)
    assert not any("level of theory" in r for r in verdict["reasons"]), verdict


def test_b2_same_level_different_capitalisation_is_not_mixed():
    spec = _spec("caps", [_sp("a", method="PBE0", basis="def2-TZVP"),
                          _sp("b", method="pbe0", basis="DEF2-TZVP")])
    verdict = grade(spec, "optimize propane and compare conformers",
                    budget_core_hours=4.0)
    assert not any("level of theory" in r for r in verdict["reasons"]), verdict


def test_b3_two_methods_in_two_orca_steps_is_mixed_red():
    # The true positive side: two orca steps at different levels compared
    # in one spec MUST be flagged.
    spec = _spec("mix", [_sp("a", method="pbe0"),
                         _sp("b", method="b3lyp")])
    verdict = grade(spec, "compare conformers of propane",
                    budget_core_hours=4.0)
    assert any("level of theory" in r for r in verdict["reasons"]), verdict


def test_b4_conformer_rule_demands_search_for_any_flexible_word():
    # R5 (the red spec the package must catch): flexible molecule, no
    # conformer step. This is the fires-when-it-should side.
    spec = _spec("flex", [_sp("opt")])
    verdict = grade(spec, "conformer landscape of flexible octanol",
                    budget_core_hours=4.0)
    assert any("conformer" in r for r in verdict["reasons"]), verdict


def test_b5_conformer_word_in_a_non_landscape_request_not_flagged():
    # A request that merely mentions a conformation in passing (asking for
    # a frequency on a rigid ring) must not be forced into a conformer
    # search — the rule keys on the word alone.
    spec = _spec("rigid", [_sp("opt")])
    verdict = grade(spec, "frequencies of a fixed chair conformation",
                    budget_core_hours=4.0)
    assert not any("conformer search" in r for r in verdict["reasons"]), verdict
