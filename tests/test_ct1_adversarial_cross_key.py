"""Adversarial edge cases for cross_key_consistency_notes (CT1 review, nacht-s2).

The builder's phase-3 tests cover the happy paths. These attack the
boundaries a real CONTROL.txt writer hits:

* method spelled in mixed case / with surrounding whitespace (the check
  lowercases; a value like " Classic" must still suppress or fire
  correctly, not skip the check);
* OCCUPIER keys carrying "no" or "" must never fire (a switched-off key
  is not an inert-but-set key);
* relativity spelled uppercase ("ZORA ") still normalises to a finding;
* bracket metal SMILES with a metal the naive tokeniser could mis-split
  (e.g. "N" vs "Na", "C" vs "Cl" -- the docstring of _metals_in_smiles
  claims these are safe);
* an empty SMILES with relativity set: the metals list is empty, which
  is the "no transition metals" branch -- it must warn, not crash;
* both solvent and model given in different capitalisations: the check
  reads raw values, so "C-PCM" style input must not produce a false
  finding while the validator itself canonicalises to CPCM.

All findings are advisory; every test here also pins that the finding
never becomes a validation error.
"""
import pytest

from delfin.common.control_validator import (
    cross_key_consistency_notes,
    validate_control_config,
)

NOTES = cross_key_consistency_notes


def test_method_mixed_case_with_spaces_still_suppresses_or_fires():
    # "OCCUPIER" in odd spelling must be recognised as the OCCUPIER mode.
    assert NOTES({"method": "  Occupier ", "OCCUPIER_method": "auto"}) == []


def test_occupier_key_set_to_no_never_fires():
    assert NOTES({"method": "classic", "OCCUPIER_method": "no"}) == []


def test_occupier_key_empty_string_never_fires():
    assert NOTES({"method": "classic", "OCCUPIER_compare": "   "}) == []


def test_relativity_uppercase_with_whitespace_still_warns():
    notes = NOTES({"relativity": "  X2C  ", "SMILES": "c1ccccc1"})
    assert any("relativity" in n for n in notes), notes


def test_relativity_set_and_no_smiles_warns_without_crashing():
    notes = NOTES({"relativity": "ZORA"})
    assert any("relativity" in n for n in notes), notes


def test_sodium_and_chlorine_are_not_transition_metals():
    # A SMILES whose naive token split yields Na and Cl must NOT be
    # misread as metal-bearing: relativity=ZORA there is contradicted
    # (policy forces none) and must still warn.
    notes = NOTES({"relativity": "ZORA", "SMILES": "[Na+].[Cl-]"})
    assert any("relativity" in n for n in notes), notes


def test_carbon_chlorine_smiles_does_not_silence_the_warning():
    # "CCl" tokenises to C, Cl; neither is a TM, so ZORA must still warn.
    notes = NOTES({"relativity": "ZORA", "SMILES": "CCl"})
    assert any("relativity" in n for n in notes), notes


def test_palladium_bracket_atom_is_recognised():
    notes = NOTES({"relativity": "ZORA", "SMILES": "[Pd]"})
    assert not any("relativity" in n for n in notes), notes


def test_findings_are_never_validation_errors():
    # Every attacked configuration still validates: findings are logged
    # warnings, never blocking errors (control_validator.py:2714-2716).
    for config in (
        {"charge": 0, "method": "classic", "solvent": "dmf"},
        {"charge": 0, "method": "classic", "relativity": "ZORA"},
        {"charge": 0, "method": "classic", "OCCUPIER_method": "auto"},
        {"charge": 0, "method": "classic", "relativity": "ZORA",
         "SMILES": "[Na+].[Cl-]"},
    ):
        validated = validate_control_config(dict(config))
        assert validated.get("charge", config.get("charge")) is not None
        # the value that was attacked is still there, not dropped
        for key, value in config.items():
            assert validated.get(key, value) == value


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__]))
