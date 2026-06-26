"""Welle-3 T2.2: pytest for the two distortion-mandatory selector extensions.

Patch file:
  welle3_T2.2_patches/distortion_mandatory_cn4_jt_and_graph_aromatic.patch

Both extensions are env-gated and default OFF.  These tests cover:

  E1. DELFIN_JT_CN4_DISTORTION=0 (default) → CN=4 Cu(II) → None  (legacy).
  E2. DELFIN_JT_CN4_DISTORTION=1 → CN=4 Cu(II) → "jahn_teller_cu2".
  E3. DELFIN_CYCLOMETAL_GRAPH_FALLBACK=0 (default) → unsanitised mol → None.
  E4. DELFIN_CYCLOMETAL_GRAPH_FALLBACK=1 → unsanitised mol → "cyclometal_5ring".
  E5. Determinism (3 runs, same input → same label).
  E6. Sanitised aromatic case still works regardless of either flag.

Skips cleanly when RDKit / _post_optimizer are missing.
"""
from __future__ import annotations

import os
import numpy as np
import pytest

_post_optimizer = pytest.importorskip("delfin._post_optimizer")
Chem = pytest.importorskip("rdkit.Chem")

_classify = getattr(_post_optimizer, "_classify_distortion_mandatory", None)
if _classify is None:
    pytest.skip("classifier missing", allow_module_level=True)


def _cu2_cn4_planar():
    mol = Chem.RWMol()
    cu = Chem.Atom("Cu"); cu.SetFormalCharge(2); ci = mol.AddAtom(cu)
    coords = [np.zeros(3)]
    for p in [(2.0, 0, 0), (-2.0, 0, 0), (0, 2.0, 0), (0, -2.0, 0)]:
        n = mol.AddAtom(Chem.Atom("N")); coords.append(np.array(p, float))
        mol.AddBond(ci, n, Chem.BondType.SINGLE)
    return mol, np.array(coords), ci, [1, 2, 3, 4]


def _ir_cyclometal_5ring(aromatic_flag: bool):
    """5-ring Ir-C(aryl)-C-C-N-Ir with optional aromatic flag."""
    mol = Chem.RWMol()
    ir = mol.AddAtom(Chem.Atom("Ir"))
    c1 = mol.AddAtom(Chem.Atom("C"))
    c2 = mol.AddAtom(Chem.Atom("C"))
    c3 = mol.AddAtom(Chem.Atom("C"))
    n  = mol.AddAtom(Chem.Atom("N"))
    bt = Chem.BondType.AROMATIC if aromatic_flag else Chem.BondType.SINGLE
    mol.AddBond(ir, c1, Chem.BondType.SINGLE)
    mol.AddBond(c1, c2, bt)
    mol.AddBond(c2, c3, bt)
    mol.AddBond(c3, n,  bt)
    mol.AddBond(n, ir, Chem.BondType.SINGLE)
    if aromatic_flag:
        for a in (c1, c2, c3, n):
            mol.GetAtomWithIdx(a).SetIsAromatic(True)
    coords = np.array([
        [0, 0, 0],          # Ir
        [2.0, 0, 0],        # c1
        [2.5, 1.5, 0],      # c2
        [1.5, 2.7, 0],      # c3
        [0, 2.0, 0],        # n
    ], dtype=float)
    return mol, coords, ir, [c1, n]


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("DELFIN_JT_CN4_DISTORTION", raising=False)
    monkeypatch.delenv("DELFIN_CYCLOMETAL_GRAPH_FALLBACK", raising=False)
    yield


def test_E1_cu2_cn4_default_off():
    mol, c, m, d = _cu2_cn4_planar()
    assert _classify(mol, c, m, d) is None


def test_E2_cu2_cn4_env_on(monkeypatch):
    monkeypatch.setenv("DELFIN_JT_CN4_DISTORTION", "1")
    mol, c, m, d = _cu2_cn4_planar()
    assert _classify(mol, c, m, d) == "jahn_teller_cu2"


def test_E3_cyclometal_graph_fallback_default_off():
    mol, c, m, d = _ir_cyclometal_5ring(aromatic_flag=False)
    assert _classify(mol, c, m, d) is None


def test_E4_cyclometal_graph_fallback_env_on(monkeypatch):
    monkeypatch.setenv("DELFIN_CYCLOMETAL_GRAPH_FALLBACK", "1")
    mol, c, m, d = _ir_cyclometal_5ring(aromatic_flag=False)
    assert _classify(mol, c, m, d) == "cyclometal_5ring"


def test_E5_determinism(monkeypatch):
    """Three identical calls must return the same label across all flag combos."""
    monkeypatch.setenv("DELFIN_JT_CN4_DISTORTION", "1")
    mol, c, m, d = _cu2_cn4_planar()
    labels = {_classify(mol, c, m, d) for _ in range(3)}
    assert labels == {"jahn_teller_cu2"}

    monkeypatch.setenv("DELFIN_CYCLOMETAL_GRAPH_FALLBACK", "1")
    mol2, c2, m2, d2 = _ir_cyclometal_5ring(aromatic_flag=False)
    labels2 = {_classify(mol2, c2, m2, d2) for _ in range(3)}
    assert labels2 == {"cyclometal_5ring"}


def test_E6_sanitised_aromatic_invariant():
    """Aromatic-flag path must work regardless of either env flag."""
    mol, c, m, d = _ir_cyclometal_5ring(aromatic_flag=True)
    assert _classify(mol, c, m, d) == "cyclometal_5ring"
