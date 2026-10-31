"""Cross-key consistency findings in the CONTROL validator.

Phase-3 of CT1. Each case: red first, warning-level (never blocking),
fired ONLY when the key is explicitly present and its meaning is
CONTRADICTED by the chosen mode/companions -- per the QS ruling
(nacht-s5): a merely-unused key must NOT warn.

Three findings:

1. solvent set, implicit_solvation_model empty: the run goes gas phase
   with no trace at validation time (delfin/common/solvation.py:31-32
   returns "" silently; the reverse case is warned at :34-44). The user
   asked for a solvent and silently does not get one -- that is the
   asymmetry.
2. relativity set to a relativistic method on a LIGHT-ELEMENT system:
   delfin/cli.py:1989-1994 silently overwrites it to 'none'. The value in
   the file is a lie about what runs. Only this direction warns; a
   missing relativity on a heavy-metal system is the auto policy's job
   (delfin/utils.py:203-217) and must NOT warn (reviewer caveat).
3. OCCUPIER-specific keys present while method is classic/manually: the
   keys are inert. Warning only, and only when method is explicit.
"""
import logging

import pytest

from delfin.common.control_validator import (
    cross_key_consistency_notes,
    validate_control_config,
)


def _notes(config):
    return cross_key_consistency_notes(config)


# ---------------------------------------------------------------- case 1

def test_solvent_without_a_solvation_model_warns():
    notes = _notes({"solvent": "dmf"})
    assert any("implicit_solvation_model" in n and "dmf" in n for n in notes), notes


def test_the_reverse_case_does_not_double_warn():
    """Model without solvent is already warned at runtime (solvation.py:34-44)."""
    assert not _notes({"implicit_solvation_model": "CPCM"})


def test_both_set_is_consistent_and_warns_nothing():
    assert _notes({"solvent": "dmf", "implicit_solvation_model": "CPCM"}) == []


def test_neither_set_is_gas_phase_and_warns_nothing():
    assert _notes({}) == []


# ---------------------------------------------------------------- case 2

def test_relativistic_setting_on_a_light_element_system_warns():
    notes = _notes({"relativity": "ZORA", "SMILES": "c1ccccc1"})
    assert any("relativity" in n for n in notes), notes


def test_relativity_none_on_a_heavy_metal_system_does_not_warn():
    """The auto policy (utils.py:203-217) handles this -- never warn on missing."""
    assert _notes({"relativity": "none", "SMILES": "[Pd](Cl)(Cl)(N)N"}) == []


def test_relativistic_setting_on_a_heavy_metal_system_does_not_warn():
    assert _notes({"relativity": "ZORA", "SMILES": "[Pd](Cl)(Cl)(N)N"}) == []


def test_relativity_warning_names_what_would_actually_run():
    notes = _notes({"relativity": "X2C", "SMILES": "c1ccccc1"})
    assert any("none" in n for n in notes), notes


# ---------------------------------------------------------------- case 3

def test_occupier_keys_with_method_classic_warn():
    """A deliberately chosen OCCUPIER setting under method=classic is inert."""
    notes = _notes({"method": "classic", "OCCUPIER_method": "manually"})
    assert any("OCCUPIER" in n for n in notes), notes


def test_a_chosen_occupier_compare_is_left_alone():
    """OCCUPIER_compare is NOT watched: the old test pins a valid value as
    left alone (test_occupier_compare_names_what_it_compares.py:201-210) --
    its value decision is always the user's, never an inconsistency."""
    assert _notes({"method": "classic", "OCCUPIER_compare": "G"}) == []


def test_template_default_occupier_keys_with_method_classic_do_not_warn():
    """The template ships OCCUPIER_compare=FSPE / OCCUPIER_method=auto; that
    shape under method=classic is a normal file, not a contradiction (pinned
    by test_occupier_compare_names_what_it_compares.py:201-210)."""
    assert _notes({"method": "classic", "OCCUPIER_method": "auto",
                   "OCCUPIER_compare": "FSPE"}) == []


def test_an_invalid_occupier_method_value_never_warns():
    """OCCUPIER_method='no' is invalid (coercer _as_occupier_method raises);
    the consistency note only speaks about values that would run."""
    notes = _notes({"method": "classic", "OCCUPIER_method": "no"})
    assert notes == [], notes


def test_geom_opt_default_is_not_warned():
    """geom_opt_OCCUPIER's coercer is _as_geom_opt (default 'OPT'), not
    yes/no -- a raw default-valued key must not trip the finding."""
    assert _notes({"method": "classic", "geom_opt_OCCUPIER": "OPT"}) == []


def test_off_spellings_map_to_the_default_through_the_coercer():
    """Reviewer finding on 64a090dfb: raw OFF-spellings of the _as_yes_no
    key coerce to 'no' (= default) and must not warn."""
    for spelling in ("false", "0", "off", "False", " OFF "):
        notes = _notes({
            "method": "classic",
            "frequency_calculation_OCCUPIER": spelling,
        })
        assert notes == [], (spelling, notes)


def test_a_real_yes_still_warns():
    notes = _notes({
        "method": "classic",
        "frequency_calculation_OCCUPIER": "yes",
    })
    assert any("OCCUPIER" in n for n in notes), notes


def test_an_odd_seq_under_method_classic_warns():
    notes = _notes({"method": "classic", "even_seq": [{"index": 1, "m": 2}]})
    assert any("OCCUPIER" in n for n in notes), notes


def test_occupier_keys_with_method_occupier_do_not_warn():
    assert _notes({"method": "OCCUPIER", "OCCUPIER_method": "auto"}) == []


def test_a_merely_unused_key_does_not_warn():
    """QS ruling: keys a workflow simply does not read are NOT findings."""
    assert _notes({"method": "classic", "TDDFT_nroots": 5}) == []


def test_method_absent_warns_nothing_case3():
    assert _notes({"OCCUPIER_method": "auto"}) == []


# ------------------------------------------------- wired into validation

def test_validate_control_config_reports_the_findings_non_blocking():
    """Findings are logged, not errors: validate_control_config returns the
    validated config (control_validator.py:2607-2608 raises only on errors)."""
    config = {"charge": 0, "method": "classic", "solvent": "dmf"}
    validated = validate_control_config(dict(config))
    assert validated["solvent"] == "dmf"  # nothing dropped, nothing refused


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__]))
