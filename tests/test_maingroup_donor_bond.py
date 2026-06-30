"""DELFIN_FFFREE_MAINGROUP_DONOR_BOND: CSD-calibrated M-(heavy main-group donor)
bond length.  The generic (cov_sum + offset) estimate over-stretches stannyl /
stibine M-Sn/M-Sb bonds to ~3.1-3.3 A vs the CCDC ~2.6 A; this flag substitutes
the calibrated value for transition/f-block metals binding a heavy main-group
donor.  Default OFF -> byte-identical."""
import os

import delfin.smiles_converter as sc


def _len(m, d):
    return round(sc._get_ml_bond_length(m, d), 2)


def test_off_is_unchanged(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_MAINGROUP_DONOR_BOND", raising=False)
    # generic over-stretched estimate (regression guard on the OFF baseline)
    assert _len("Ir", "Sn") > 3.0
    assert _len("Ru", "Sb") > 3.0


def test_on_uses_csd_value(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_MAINGROUP_DONOR_BOND", "1")
    assert _len("Ir", "Sn") == 2.62
    assert _len("Ru", "Sb") == 2.62
    assert _len("Rh", "Sn") == 2.62
    assert _len("Au", "Sb") == 2.62
    assert _len("Os", "Bi") == 2.78


def test_on_leaves_normal_donors_and_nonTM_alone(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_MAINGROUP_DONOR_BOND", "1")
    # ordinary donors (C/N/O/P) unaffected
    assert _len("Ru", "C") < 2.2
    assert _len("Ir", "P") < 2.5
    # a main-group metal AS THE METAL (Sn-C inside SnMe3) is NOT a corrected M-E bond
    assert _len("Sn", "C") < 2.5
    assert _len("Sn", "C") == round(sc._get_ml_bond_length("Sn", "C"), 2)


def test_table_values_are_physical():
    # every calibrated M-E bond is in a sane covalent range
    for e, d in sc._MAINGROUP_DONOR_BOND.items():
        assert 2.3 <= d <= 3.0, f"{e}: {d}"
