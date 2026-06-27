"""Tests for the under-coordinated bent-donor VSEPR fix (DELFIN_FFFREE_DONOR_BEND).

A 2-coordinate chalcogen / pyramidal-pnictogen donor (M + ONE substituent) must
be placed BENT at its VSEPR ideal M-D-X angle (chalcogen ~95-109 deg, pyramidal
pnictogen ~107 deg), NOT colinear with the M-D axis (180 deg).  The eye-validated
defect: YULFAO = EtHg-Se-C6H5 was built with Hg-Se-C(phenyl) = 180 deg (linear),
which must become ~96-104 deg while the Hg centre stays linear (180 deg, d10 CN2).

Acceptance contract
-------------------
1. Flag OFF (unset / "0") => byte-identical placement (every donor stays at its
   historic linear 180 deg).
2. Flag ON => single-substituent (k=1) chalcogen / pyramidal-pnictogen donors bend
   to their VSEPR ideal; the metal centre geometry is unchanged.
3. Genuinely-linear donors (sp nitrile N, double-bonded sp chalcogen) stay 180 deg
   even with the flag ON (no over-bending).
4. Deterministic: two ON builds are bit-identical.
"""
from __future__ import annotations

import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import delfin.manta.assemble_complex as AC


def _frag(smiles: str, donor_sym: str):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=AC.SEED)
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    di = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == donor_sym)
    return m, di


def _angle(P, i, j, k):
    a = P[i] - P[j]
    b = P[k] - P[j]
    return float(np.degrees(np.arccos(np.clip(
        np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1.0, 1.0))))


def _build_cn2(donor_smiles, donor_sym, metal="Pt", spectator="[Cl-]",
               spectator_sym="Cl", flag="0", refine=False):
    """Build a CN2 complex (donor ligand + spectator halide) and return the
    M-D-X angle at the donor (X = donor's first heavy substituent)."""
    os.environ["DELFIN_FFFREE_DONOR_BEND"] = flag
    vs = [_frag(spectator, spectator_sym), _frag(donor_smiles, donor_sym)]
    out = AC.assemble_heteroleptic_from_mols(metal, "L-2 linear", vs, refine=refine)
    assert out is not None
    syms, P = out
    D = [i for i, s in enumerate(syms) if s == donor_sym]
    d = min(D, key=lambda i: np.linalg.norm(P[i] - P[0]))
    heavy = [i for i in range(len(syms)) if i not in (0, d) and syms[i] != "H"]
    x = min(heavy, key=lambda i: np.linalg.norm(P[i] - P[d]))
    return _angle(P, 0, d, x), syms, P


@pytest.fixture(autouse=True)
def _clean_env():
    old = os.environ.get("DELFIN_FFFREE_DONOR_BEND")
    yield
    if old is None:
        os.environ.pop("DELFIN_FFFREE_DONOR_BEND", None)
    else:
        os.environ["DELFIN_FFFREE_DONOR_BEND"] = old


# --------------------------------------------------------------------------- #
# 1. bent-capable k=1 donors bend ON, stay linear OFF
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("smiles,sym,ideal", [
    ("[Se-]c1ccccc1", "Se", 98.0),    # selenolate (the YULFAO Se)
    ("[S-]c1ccccc1", "S", 100.0),     # thiolate
    ("[O-]c1ccccc1", "O", 109.0),     # phenoxide
    ("[O-]C", "O", 109.0),            # alkoxide
    ("[Te-]c1ccccc1", "Te", 95.0),    # tellurolate
])
def test_bent_capable_donor_bends_on(smiles, sym, ideal):
    off, _, _ = _build_cn2(smiles, sym, flag="0")
    on, _, _ = _build_cn2(smiles, sym, flag="1")
    assert abs(off - 180.0) < 1.0, f"{sym}: OFF must be linear, got {off:.1f}"
    assert abs(on - ideal) < 2.0, f"{sym}: ON must be ~{ideal}, got {on:.1f}"


# --------------------------------------------------------------------------- #
# 2. genuinely-linear donors stay linear even ON (no over-bending)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("smiles,sym", [
    ("N#CC", "N"),          # nitrile sp N
    ("S=C=NC", "S"),        # isothiocyanate, double-bonded sp chalcogen
])
def test_linear_donor_stays_linear(smiles, sym):
    off, _, _ = _build_cn2(smiles, sym, flag="0")
    on, _, _ = _build_cn2(smiles, sym, flag="1")
    assert abs(off - 180.0) < 1.0
    assert abs(on - 180.0) < 1.0, f"{sym}: must NOT be bent, got {on:.1f}"


# --------------------------------------------------------------------------- #
# 3. aromatic in-ring N is not bent by the flag (ring guard / aromatic guard)
# --------------------------------------------------------------------------- #
def test_aromatic_n_unchanged():
    off, _, _ = _build_cn2("c1ccncc1", "N", flag="0")
    on, _, _ = _build_cn2("c1ccncc1", "N", flag="1")
    assert abs(off - on) < 1e-6


# --------------------------------------------------------------------------- #
# 4. YULFAO: Hg-Se-C bends; Hg centre stays linear
# --------------------------------------------------------------------------- #
def test_yulfao_se_bends_hg_linear():
    for flag, expect_se in (("0", 180.0), ("1", 98.0)):
        os.environ["DELFIN_FFFREE_DONOR_BEND"] = flag
        vs = [_frag("[CH2-]C", "C"), _frag("[Se-]c1ccccc1", "Se")]
        out = AC.assemble_heteroleptic_from_mols("Hg", "L-2 linear", vs, refine=True)
        syms, P = out
        se = [i for i, s in enumerate(syms) if s == "Se"][0]
        cs = [i for i, s in enumerate(syms) if s == "C"]
        cph = min(cs, key=lambda i: np.linalg.norm(P[i] - P[se]))
        cet = min(cs, key=lambda i: np.linalg.norm(P[i] - P[0]))
        assert abs(_angle(P, 0, se, cph) - expect_se) < 2.0
        assert abs(_angle(P, cet, 0, se) - 180.0) < 1.0   # Hg centre linear


# --------------------------------------------------------------------------- #
# 5. determinism ON
# --------------------------------------------------------------------------- #
def test_determinism_on():
    os.environ["DELFIN_FFFREE_DONOR_BEND"] = "1"
    vs1 = [_frag("[CH2-]C", "C"), _frag("[Se-]c1ccccc1", "Se")]
    vs2 = [_frag("[CH2-]C", "C"), _frag("[Se-]c1ccccc1", "Se")]
    a = AC.assemble_heteroleptic_from_mols("Hg", "L-2 linear", vs1, refine=True)[1]
    b = AC.assemble_heteroleptic_from_mols("Hg", "L-2 linear", vs2, refine=True)[1]
    assert np.array_equal(a, b)
