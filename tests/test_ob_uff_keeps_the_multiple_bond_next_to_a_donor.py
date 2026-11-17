"""Open Babel UFF keeps the SMILES bond orders next to a metal donor (DELFIN_FFFREE_OB_KEKULE_ORDERS).

The UFF step reads a bare xyz.  Open Babel then guesses the bond orders with the metal in the
molecule, the metal-donor bond fills the donor's valence, and the multiple bond next to the donor
is guessed single: a coordinated acetone C=O is minimised to single-bond length, and a coordinated
pyridine to a Kekule ring of alternating short and long bonds.  With the switch on, the template's
Kekule orders and aromaticity are used instead.  With the switch off, nothing changes.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("openbabel")
Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta import openbabel_optimize as obo  # noqa: E402

# trans-dichloro(acetone)(pyridine)zinc, RDKit ETKDG geometry (seed 7), template atom order
_SMILES = "CC(C)=[O+][Zn-2]([N+]1=CC=CC=C1)(Cl)Cl"
_XYZ = (
    "C -4.1815 -2.0130 -0.3960\n"
    "C -2.8920 -1.2639 -0.6103\n"
    "C -2.2128 -0.7081 0.5871\n"
    "O -2.4690 -1.1496 -1.7475\n"
    "Zn -0.9306 -0.2675 -2.0969\n"
    "N -0.8836 1.3768 -1.1393\n"
    "C 0.2040 1.7466 -0.4344\n"
    "C 0.2554 2.9215 0.2578\n"
    "C -0.8618 3.7572 0.2257\n"
    "C -1.9892 3.4008 -0.4896\n"
    "C -1.9753 2.1957 -1.1693\n"
    "Cl -0.6946 0.1017 -4.2333\n"
    "Cl 0.7087 -1.5408 -1.3592\n"
    "H -4.8782 -1.7121 -1.1976\n"
    "H -3.9371 -3.0891 -0.5353\n"
    "H -4.5751 -1.8630 0.6126\n"
    "H -1.1290 -0.6947 0.3840\n"
    "H -2.3909 -1.4280 1.4196\n"
    "H -2.5683 0.2682 0.9140\n"
    "H 1.0769 1.1064 -0.4021\n"
    "H 1.1113 3.2535 0.8357\n"
    "H -0.8767 4.6944 0.7524\n"
    "H -2.8738 4.0477 -0.5234\n"
    "H -2.8397 1.8844 -1.7394\n"
)
_CO = (1, 3)
_RING = [(5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 5)]


def _template():
    return Chem.AddHs(Chem.MolFromSmiles(_SMILES))


def _coords(xyz):
    return np.array([[float(v) for v in ln.split()[1:4]] for ln in xyz.strip().splitlines()])


def _d(X, i, j):
    return float(np.linalg.norm(X[i] - X[j]))


def _minimise(monkeypatch, on):
    if on:
        monkeypatch.setenv("DELFIN_FFFREE_OB_KEKULE_ORDERS", "1")
    else:
        monkeypatch.delenv("DELFIN_FFFREE_OB_KEKULE_ORDERS", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_OB_BOND_ORDERS", raising=False)
    return obo._optimize_xyz_openbabel(_XYZ, 500, None, False, _template())


def test_the_template_atom_order_matches_the_geometry():
    t = _template()
    assert [a.GetSymbol() for a in t.GetAtoms()] == [ln.split()[0] for ln in _XYZ.strip().splitlines()]


def test_off_the_template_is_not_read(monkeypatch):
    off = _minimise(monkeypatch, on=False)
    monkeypatch.delenv("DELFIN_FFFREE_OB_KEKULE_ORDERS", raising=False)
    assert off == obo._optimize_xyz_openbabel(_XYZ, 500, None, False, None)


def test_off_the_guess_stretches_the_coordinated_carbonyl(monkeypatch):
    """Documents the defect the switch repairs: the guessed C-O single bond."""
    X = _coords(_minimise(monkeypatch, on=False))
    assert _d(X, *_CO) > 1.33


def test_on_the_coordinated_carbonyl_keeps_its_double_bond(monkeypatch):
    X = _coords(_minimise(monkeypatch, on=True))
    assert 1.18 < _d(X, *_CO) < 1.27


def test_on_the_coordinated_pyridine_stays_aromatic(monkeypatch):
    X = _coords(_minimise(monkeypatch, on=True))
    ring = [_d(X, i, j) for i, j in _RING]
    assert max(ring) - min(ring) < 0.05
    assert max(ring) < 1.43


def test_a_template_in_another_atom_order_is_not_used():
    from openbabel import pybel
    ob = pybel.readstring("xyz", f"{len(_XYZ.splitlines())}\n\n{_XYZ}").OBMol
    t = Chem.AddHs(Chem.MolFromSmiles("Cl[Zn-2](Cl)([N+]1=CC=CC=C1)[O+]=C(C)C"))
    assert obo._apply_template_kekule_orders(ob, t) == 0
