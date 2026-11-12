"""Adversarial P1 tests (reviewer nacht-s15) for WB1 rules 3/4/10.

Probes the real entry point scientific_lint. Each case states what the
physically right answer is; a failure means the rule mis-fires or misses.
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _step(label, **kw):
    return {"name": "s", "steps": [{"step": "orca_sp", "label": label, "kwargs": kw}]}


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


# ---------------------------------------------------------- rule 3 positives
def test_rule3_fires_for_anion_with_plain_basis():
    rules = _rules(_step("sp", charge=-1, basis="def2-TZVP"))
    assert ("diffuse-basis-required", "sp", "warning") in rules


def test_rule3_fires_for_rydberg_marker():
    rules = _rules(_step("sp", charge=0, basis="cc-pVDZ",
                         excited_state="Rydberg"))
    assert ("diffuse-basis-required", "sp", "warning") in rules


def test_rule3_silent_on_aug_marker():
    assert scientific_lint(_step("sp", charge=-1, basis="aug-cc-pVDZ")) == []


# ---------------------------------------------------------- rule 3 negatives
def test_rule3_uppercase_aug_marker_is_still_diffuse():
    # DELFIN users write AUG-cc-pVDZ; the marker check is case-sensitive.
    rules = _rules(_step("sp", charge=-1, basis="AUG-cc-pVDZ"))
    assert ("diffuse-basis-required", "sp", "warning") not in rules


def test_rule3_uppercase_ma_marker_is_still_minimally_augmented():
    rules = _rules(_step("sp", charge=-1, basis="MA-def2-TZVP"))
    assert ("diffuse-basis-required", "sp", "warning") not in rules


# ---------------------------------------------------------- rule 4
def test_rule4_fires_for_pt_without_rel():
    rules = _rules(_step("sp", metals=["Pt"], basis="def2-TZVP"))
    assert ("relativistic-required", "sp", "warning") in rules


def test_rule4_silent_when_relativity_set_with_rel_basis():
    # relativity=ZORA + a REL basis is the complete, correct setup.
    # (Note: relativity=ZORA with a non-rel basis def2-TZVP is correctly
    # flagged by the pre-existing relativity-without-rel-basis rule.)
    assert scientific_lint(_step("sp", metals=["Pt"], basis="zora-def2-TZVP",
                                 relativity="ZORA")) == []


def test_rule4_uppercase_zora_basis_is_a_rel_basis():
    # ZORA-def2-SVP (uppercase) is DELFIN's own rel default spelling
    # (utils.select_rel_and_aux).  The marker check must not be fooled.
    rules = _rules(_step("sp", metals=["Pt"], basis="ZORA-def2-SVP"))
    assert ("relativistic-required", "sp", "warning") not in rules


def test_rule4_uppercase_dkh_basis_is_a_rel_basis():
    rules = _rules(_step("sp", metals=["Ag"], basis="DKH-def2-TZVP"))
    assert ("relativistic-required", "sp", "warning") not in rules


def test_rule4_agrees_with_should_use_rel_on_non_rel_metals():
    # Zn, Fe are not 4d/5d/4f/5f: the rule must stay silent.
    assert scientific_lint(_step("sp", metals=["Zn", "Fe"],
                                 basis="def2-TZVP")) == []


# ---------------------------------------------------------- rule 10
def test_rule10_fires_on_same_parity():
    rules = _rules(_step("sp", mult=4, electrons=4))
    assert ("electron-multiplicity-parity", "sp", "error") in rules


def test_rule10_silent_on_legal_pair():
    assert scientific_lint(_step("sp", mult=2, electrons=3)) == []


def test_rule10_mult_1_with_odd_electrons_is_impossible():
    # mult=1 means S=0 (all paired) — impossible with an odd electron
    # count, yet the parities differ so parity alone passes it.
    # Physically right: this must be flagged.
    rules = _rules(_step("sp", mult=1, electrons=3))
    assert ("electron-multiplicity-parity", "sp", "error") in rules


def test_rule10_mult_exceeding_electron_capacity_is_impossible():
    # 3 electrons cannot give mult=7 (needs 4 unpaired > 3 electrons).
    rules = _rules(_step("sp", mult=7, electrons=3))
    assert ("electron-multiplicity-parity", "sp", "error") in rules
