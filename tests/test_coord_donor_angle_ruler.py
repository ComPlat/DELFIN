"""Tests for the donor-internal M-D-X coordination-angle RULER.

scripts/coord_donor_angle_ruler.py closes a metric blindspot: the coordination
ruler measured only the D-M-D angle AT the metal and was BLIND to the donor-
internal M-D-X angle.  A 2-coordinate chalcogen / pnictogen donor placed LINEAR
(M-D-X ~180 deg) must be flagged; a genuinely-linear donor (sp nitrile, aromatic
N, multiply-bonded oxo) must NOT be (CCDC-sane).
"""
from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import delfin.fffree.assemble_complex as AC
from scripts.coord_donor_angle_ruler import measure_donor_angles


def _frag(smiles, donor_sym):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=AC.SEED)
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    di = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == donor_sym)
    return m, di


def _build(donor_smiles, donor_sym, bend, metal="Pt", spectator="[Cl-]",
           spectator_sym="Cl"):
    import os
    os.environ["DELFIN_FFFREE_DONOR_BEND"] = "1" if bend else "0"
    try:
        vs = [_frag(spectator, spectator_sym), _frag(donor_smiles, donor_sym)]
        syms, P = AC.assemble_heteroleptic_from_mols(metal, "L-2 linear", vs,
                                                     refine=True)
        return syms, np.asarray(P)
    finally:
        os.environ.pop("DELFIN_FFFREE_DONOR_BEND", None)


# 1. a linear (un-bent) chalcogen/pnictogen donor is FLAGGED.
@pytest.mark.parametrize("smiles,sym", [
    ("[Se-]c1ccccc1", "Se"),
    ("[S-]c1ccccc1", "S"),
    ("[O-]c1ccccc1", "O"),
])
def test_linearised_donor_flagged(smiles, sym):
    syms, P = _build(smiles, sym, bend=False)
    res = measure_donor_angles(syms, P)
    assert res["n_flagged"] >= 1
    f = res["flagged"][0]
    assert f["D_sym"] == sym
    assert f["observed_angle"] > 150.0
    assert f["linearised"] is True


# 2. the bent (fixed) donor is CLEAN.
@pytest.mark.parametrize("smiles,sym", [
    ("[Se-]c1ccccc1", "Se"),
    ("[S-]c1ccccc1", "S"),
    ("[O-]c1ccccc1", "O"),
])
def test_bent_donor_clean(smiles, sym):
    syms, P = _build(smiles, sym, bend=True)
    res = measure_donor_angles(syms, P)
    assert res["n_donors"] >= 1
    assert res["n_flagged"] == 0


# 3. genuinely-linear donors are NEVER flagged (CCDC-sane).
#    nitrile sp N (N#C) and aromatic in-ring N keep the metal antiperiplanar /
#    on the ring lone-pair axis -> linear is correct, must not flag.
#    (NB: an S-bonded thiocyanate M-S-C IS bent ~100 deg in real crystals, so a
#    LINEARISED S there is correctly a FLAG, not a false positive -> not listed
#    here.)
@pytest.mark.parametrize("smiles,sym", [
    ("N#CC", "N"),            # nitrile sp N
    ("c1ccncc1", "N"),        # aromatic pyridine N
])
def test_linear_donor_not_flagged(smiles, sym):
    for bend in (False, True):
        syms, P = _build(smiles, sym, bend=bend)
        res = measure_donor_angles(syms, P)
        assert res["n_flagged"] == 0, f"{sym} wrongly flagged (bend={bend})"


# 4. YULFAO: the eye-validated case.
def test_yulfao_flagged_then_clean():
    import os
    os.environ["DELFIN_FFFREE_DONOR_BEND"] = "0"
    vs = [_frag("[CH2-]C", "C"), _frag("[Se-]c1ccccc1", "Se")]
    syms, P = AC.assemble_heteroleptic_from_mols("Hg", "L-2 linear", vs, refine=True)
    off = measure_donor_angles(syms, np.asarray(P))
    os.environ["DELFIN_FFFREE_DONOR_BEND"] = "1"
    vs = [_frag("[CH2-]C", "C"), _frag("[Se-]c1ccccc1", "Se")]
    syms, P = AC.assemble_heteroleptic_from_mols("Hg", "L-2 linear", vs, refine=True)
    on = measure_donor_angles(syms, np.asarray(P))
    os.environ.pop("DELFIN_FFFREE_DONOR_BEND", None)
    assert off["n_flagged"] == 1 and off["flagged"][0]["D_sym"] == "Se"
    assert on["n_flagged"] == 0


# 5. determinism: same input => same result.
def test_determinism():
    syms, P = _build("[Se-]c1ccccc1", "Se", bend=False)
    a = measure_donor_angles(syms, P)
    b = measure_donor_angles(syms, P)
    assert a == b
