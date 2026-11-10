"""The template's own metal fixes the spin of a two-donor fragment (DELFIN_FFFREE_PROCRUSTES_METAL_SPIN)."""
import math

import numpy as np
import pytest

pytest.importorskip("rdkit.Chem")

from delfin.manta import mol_prep as MP  # noqa: E402
from delfin.manta.ml_tables import _METAL_SET  # noqa: E402
from delfin.manta.topo_isomers import _build_topology_xyz_from_template  # noqa: E402

SMI = "Cl[Zn-2]1(Cl)[NH2+]CC[NH2+]1"


def _setup():
    mol = MP._prepare_mol_for_embedding(SMI)
    mol.RemoveAllConformers()
    ids = MP._embed_multiple_confs_robust(mol, 1, 42)
    assert ids
    metal = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET)
    donors = [n.GetIdx() for n in mol.GetAtomWithIdx(metal).GetNeighbors()]
    return mol, metal, donors, int(ids[0])


def _parse(xyz):
    return np.array([[float(v) for v in ln.split()[1:4]] for ln in xyz.strip().splitlines()])


def _n_h_to_metal(mol, X, metal, donors):
    """Smallest metal...H distance over the H on the two N donors."""
    out = 9.0
    for d in donors:
        a = mol.GetAtomWithIdx(d)
        if a.GetSymbol() != "N":
            continue
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 1:
                out = min(out, float(np.linalg.norm(X[nb.GetIdx()] - X[metal])))
    return out


def test_metal_spin_off_is_byte_identical(monkeypatch):
    mol, metal, donors, cid = _setup()
    monkeypatch.delenv("DELFIN_FFFREE_PROCRUSTES_METAL_SPIN", raising=False)
    a = _build_topology_xyz_from_template(mol, metal, donors, list(range(len(donors))), "TH", False, conf_id=cid)
    monkeypatch.setenv("DELFIN_FFFREE_PROCRUSTES_METAL_SPIN", "0")
    b = _build_topology_xyz_from_template(mol, metal, donors, list(range(len(donors))), "TH", False, conf_id=cid)
    assert a == b


def test_metal_spin_keeps_donors_and_amine_h_off_the_metal(monkeypatch):
    mol, metal, donors, cid = _setup()
    monkeypatch.setenv("DELFIN_FFFREE_PROCRUSTES_METAL_SPIN", "1")
    xyz = _build_topology_xyz_from_template(mol, metal, donors, list(range(len(donors))), "TH", False, conf_id=cid)
    X = _parse(xyz)
    for d in donors:
        assert 1.5 < float(np.linalg.norm(X[d] - X[metal])) < 3.0
    assert _n_h_to_metal(mol, X, metal, donors) > 2.0
