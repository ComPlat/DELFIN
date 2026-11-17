"""A Kekule single bond of a conjugated ring is read against the aromatic reference
(DELFIN_FFFREE_TEMPLATE_COMPRESSED_CONJ, 28.10.2026)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta.conformer_pools import _conjugated_ring_single_bonds, _frag_xyz_compressed  # noqa: E402


def _pyrrolide():
    m = Chem.MolFromSmiles("C1=CC=C[N-]1", sanitize=False)
    m.UpdatePropertyCache(strict=False)
    n = m.GetNumAtoms()
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    X = np.c_[1.15 * np.cos(ang), 1.15 * np.sin(ang), np.zeros(n)]   # ring bond ~1.35 A, C-N(-) 1.35
    return m, X


def test_conjugated_single_bonds_found():
    m, _X = _pyrrolide()
    s = _conjugated_ring_single_bonds(m)
    singles = {b.GetIdx() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE}
    assert singles and singles <= s


def test_flag_changes_only_the_reference(monkeypatch):
    m, X = _pyrrolide()
    X = X * (1.25 / 1.35)          # C-N 1.25 A: below 0.88 x 1.47, above 0.88 x 0.91 x 1.47
    fl = list(range(m.GetNumAtoms()))
    monkeypatch.delenv("DELFIN_FFFREE_TEMPLATE_COMPRESSED_CONJ", raising=False)
    assert _frag_xyz_compressed(m, fl, X)
    monkeypatch.setenv("DELFIN_FFFREE_TEMPLATE_COMPRESSED_CONJ", "1")
    assert not _frag_xyz_compressed(m, fl, X)
    assert _frag_xyz_compressed(m, fl, X * 0.8)


def test_stretched_flag(monkeypatch):
    m, X = _pyrrolide()
    fl = list(range(m.GetNumAtoms()))
    monkeypatch.delenv("DELFIN_FFFREE_TEMPLATE_COMPRESSED_CONJ", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_TEMPLATE_STRETCHED", raising=False)
    Y = X * 1.35                      # ring bonds ~1.82 A
    assert not _frag_xyz_compressed(m, fl, Y)
    monkeypatch.setenv("DELFIN_FFFREE_TEMPLATE_STRETCHED", "1")
    assert _frag_xyz_compressed(m, fl, Y)
    assert not _frag_xyz_compressed(m, fl, X)
