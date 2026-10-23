"""A compressed whole-complex template fragment is re-embedded alone
(DELFIN_FFFREE_TEMPLATE_COMPRESSED_RESEAT)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from delfin.manta.conformer_pools import (  # noqa: E402
    _frag_xyz_compressed,
    _isolated_fragment_xyz,
    template_compressed_reseat_enabled,
)


def _pyridine_with_metal():
    """4-methylpyridine N-bound to a Zn placeholder, explicit H, ETKDG geometry."""
    m = Chem.AddHs(Chem.MolFromSmiles("Cc1cc[n+](cc1)[Zn]"))
    assert AllChem.EmbedMolecule(m, randomSeed=11) == 0
    frag = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() != "Zn"]
    return m, frag, m.GetConformer().GetPositions()[frag]


def test_switch_is_off_by_default(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_TEMPLATE_COMPRESSED_RESEAT", raising=False)
    assert not template_compressed_reseat_enabled()


def test_a_shrunken_ring_counts_as_compressed_and_the_isolated_embed_does_not():
    mol, frag, X = _pyridine_with_metal()
    assert not _frag_xyz_compressed(mol, frag, X)
    ring = [k for k, a in enumerate(frag) if mol.GetAtomWithIdx(a).GetIsAromatic()]
    Y = X.copy()
    c = Y[ring].mean(axis=0)
    Y[ring] = c + 0.85 * (Y[ring] - c)              # the ring at 85 % of its size (C-N ~1.14 A)
    assert _frag_xyz_compressed(mol, frag, Y)
    Z = _isolated_fragment_xyz(mol, frag)
    assert Z is not None and Z.shape == X.shape
    assert not _frag_xyz_compressed(mol, frag, Z)
