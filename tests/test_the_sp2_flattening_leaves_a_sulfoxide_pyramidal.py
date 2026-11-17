"""The post-UFF sp2 flattening leaves a sulfoxide S pyramidal (DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL).

The flattening projects every atom that carries a double bond and has three heavy non-metal
neighbours into the plane of those neighbours.  A sulfoxide S (O-bound DMSO) is such an atom, but
it is pyramidal: flattening it shortens its S=O and S-C bonds.  With the switch on, a period-3+
atom with three sigma partners and no metal neighbour is left alone.  Off, nothing changes.
"""

from __future__ import annotations

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta import topology_checks as tc  # noqa: E402

_SMILES = "CS(C)=[O+][Zn-2](Cl)(Cl)Cl"
_XYZ = (
    "C -1.9045 -1.1101 -2.3139\n"
    "S -0.3801 -0.4497 -1.6460\n"
    "C 0.0852 -1.3845 -0.1932\n"
    "O -0.5288 1.0190 -1.3854\n"
    "Zn -1.9772 2.0028 -1.6776\n"
    "Cl -1.9381 3.6432 -0.1797\n"
    "Cl -1.8838 2.9365 -3.6805\n"
    "Cl -3.8475 0.8902 -1.4033\n"
    "H -2.5790 -1.3212 -1.4578\n"
    "H -2.3221 -0.3775 -3.0076\n"
    "H -1.7401 -2.0794 -2.8291\n"
    "H -0.8260 -1.5159 0.4160\n"
    "H 0.8994 -0.8803 0.3794\n"
    "H 0.3923 -2.3938 -0.5461\n"
)


def _coords(xyz):
    return np.array([[float(v) for v in ln.split()[1:4]] for ln in xyz.strip().splitlines()])


def _run_flattening(monkeypatch, on):
    if on:
        monkeypatch.setenv("DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL", "1")
    else:
        monkeypatch.delenv("DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL", raising=False)
    return tc._flatten_sp2_atoms_xyz(_XYZ, Chem.AddHs(Chem.MolFromSmiles(_SMILES)))


@pytest.mark.parametrize("keep_length", ["0", "1"])
def test_off_the_sulfoxide_s_is_flattened_and_its_s_o_shortened(monkeypatch, keep_length):
    """Documents the defect the switch repairs."""
    monkeypatch.setenv("DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH", keep_length)
    X0, X = _coords(_XYZ), _coords(_run_flattening(monkeypatch, on=False))
    assert np.linalg.norm(X[1] - X[3]) < np.linalg.norm(X0[1] - X0[3]) - 0.05


@pytest.mark.parametrize("keep_length", ["0", "1"])
def test_on_the_sulfoxide_s_keeps_its_place(monkeypatch, keep_length):
    monkeypatch.setenv("DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH", keep_length)
    assert np.allclose(_coords(_run_flattening(monkeypatch, on=True)), _coords(_XYZ), atol=1e-5)


def test_on_a_carbonyl_carbon_is_still_flattened(monkeypatch):
    """Period 2 keeps the flattening: acetone C pushed out of plane comes back."""
    xyz = ("C 0.0 0.0 0.0\nC 1.3 0.0 0.0\nC 1.95 1.12 0.30\nO 1.95 -1.06 0.30\n"
           "H -0.4 1.0 0.0\nH -0.4 -0.5 0.9\nH -0.4 -0.5 -0.9\n"
           "H 3.0 1.0 0.3\nH 1.6 1.7 1.1\nH 1.6 1.7 -0.5\n")
    monkeypatch.setenv("DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL", "1")
    monkeypatch.delenv("DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH", raising=False)
    out = _coords(tc._flatten_sp2_atoms_xyz(xyz, Chem.AddHs(Chem.MolFromSmiles("CC(C)=O"))))
    assert np.linalg.norm(out[1] - _coords(xyz)[1]) > 0.05


# S-bound DMSO: C, C, O and the metal around the S
_SMILES_S = "C[S+](C)(=O)[Zn-2](Cl)(Cl)Cl"
_XYZ_S = (
    "C -2.1859 -1.4724 -1.9950\n"
    "S -0.9249 -0.2247 -1.7556\n"
    "C -0.0703 -0.6456 -0.2664\n"
    "O 0.0647 -0.2956 -2.8818\n"
    "Zn -1.9084 1.7274 -1.6507\n"
    "Cl -1.6119 2.6825 0.3039\n"
    "Cl -0.9857 3.0439 -3.1698\n"
    "Cl -4.0766 1.5397 -2.0710\n"
    "H -2.8669 -1.3861 -1.1330\n"
    "H -2.7441 -1.1798 -2.9118\n"
    "H -1.6716 -2.4307 -2.0872\n"
    "H -0.7625 -0.8845 0.5669\n"
    "H 0.7078 0.0941 0.0114\n"
    "H 0.4864 -1.5886 -0.4847\n"
)


@pytest.mark.parametrize("on", [False, True])
def test_an_s_bound_sulfoxide_is_flattened_only_off(monkeypatch, on):
    monkeypatch.delenv("DELFIN_FFFREE_DONOR_SIGMA_COUNT", raising=False)
    if on:
        monkeypatch.setenv("DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL", "1")
    else:
        monkeypatch.delenv("DELFIN_FFFREE_FLATTEN_P3_PYRAMIDAL", raising=False)
    X0 = _coords(_XYZ_S)
    X = _coords(tc._flatten_sp2_atoms_xyz(_XYZ_S, Chem.AddHs(Chem.MolFromSmiles(_SMILES_S))))
    moved = float(np.linalg.norm(X[1] - X0[1]))
    assert (moved < 1e-5) if on else (moved > 0.05)
