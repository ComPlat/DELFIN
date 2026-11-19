"""The relativity cross-key note must see every metal, in every SMILES spelling.

Regression: a CONTROL with SMILES=[Pd](Cl)(Cl)([NH3])[NH3] and
relativity=ZORA printed "relativity=zora on a light-element system (no
transition metals in the SMILES)".  The note came from the define.py
TEMPLATE (SMILES= empty, relativity=ZORA), which config._load_template_defaults
validates on every run: an absent SMILES was read as "no metals".  Metals
are now read only from the bracket atoms of a SMILES that is actually there.
"""
import logging

import pytest

from delfin.common.control_validator import cross_key_consistency_notes
from delfin.utils import smiles_transition_metals


def _rel_notes(smiles, relativity="ZORA"):
    notes = cross_key_consistency_notes({"relativity": relativity, "SMILES": smiles})
    return [n for n in notes if n.startswith("relativity=")]


@pytest.mark.parametrize("smiles, metal", [
    ("[Pd](Cl)(Cl)([NH3])[NH3]", "Pd"),
    ("[Pd+2].[Cl-].[Cl-]", "Pd"),
    ("[Pd-2](Cl)(Cl)(Cl)Cl", "Pd"),
    ("[PdH](Cl)(P(C)(C)C)P(C)(C)C", "Pd"),
    ("[106Pd]", "Pd"),
    ("Cl[Pt@SP1](Cl)([NH3])[NH3]", "Pt"),
    ("[Ru+2]", "Ru"),
    ("[Ir]", "Ir"),
    ("[Au+]", "Au"),
    ("[W](=O)(=O)(=O)=O", "W"),
    ("[Gd+3]", "Gd"),
    ("[Eu]", "Eu"),
    ("[U+4]", "U"),
    ("O=[U]=O", "U"),
    ("[Th]", "Th"),
])
def test_heavy_metal_in_any_spelling_is_seen_and_does_not_warn(smiles, metal):
    assert metal in smiles_transition_metals(smiles)
    assert _rel_notes(smiles) == []
    assert _rel_notes(smiles, "X2C") == []


@pytest.mark.parametrize("smiles", [
    "CCO",
    "c1ccccc1",
    "CSc1ccccc1",      # thioanisole: S + aromatic c, not scandium
    "Cn1ccnc1",        # N-methylimidazole: C + aromatic n, not copernicium
    "COs",             # no osmium outside a bracket
    "[Na+].[Cl-]",
    "ClC(Cl)Cl",
])
def test_organic_or_main_group_smiles_still_warns(smiles):
    assert smiles_transition_metals(smiles) == []
    notes = _rel_notes(smiles)
    assert len(notes) == 1 and "no transition metals in the SMILES" in notes[0], notes


def test_a_3d_metal_is_named_and_still_warns():
    notes = _rel_notes("[Fe+2]")
    assert len(notes) == 1 and "metals: Fe, class 3d" in notes[0], notes


def test_no_smiles_is_no_evidence_and_never_warns():
    assert _rel_notes("") == []


def test_the_template_defaults_never_print_the_relativity_note(caplog):
    import delfin.config as config_mod

    config_mod._TEMPLATE_DEFAULTS_CACHE = None
    with caplog.at_level(logging.WARNING):
        config_mod._load_template_defaults()
    assert not any("light-element" in r.getMessage() for r in caplog.records), \
        [r.getMessage() for r in caplog.records]
