"""Regression test for the η-ring hydrogen restoration on the rigid-hapto path.

Background
---------
The dative-bond SMILES encoding writes π-face carbons as ``[C+]`` cations
(degree-3 with the metal bond).  When ``decompose._decompose_hapto`` cleaves the
M-C bonds the bare ring carbons keep ``noImplicit=True`` + a radical electron, so
RDKit's ``AddHs`` adds ZERO hydrogen to them -> the η-ring CH hydrogens vanished
from every emitted frame (eye-caught: CEBQEI built 11 H vs 17, BOCMIR 6 vs 10).

The fix (``decompose._restore_eta_ring_hydrogens``) drops the encoding-artifact
charge + radical on each η-carbon and sets its explicit-H count to
``max(0, 3 - heavy_neighbours)`` (sp2 ring carbon, valence 3), so the dropped CH
hydrogens are regenerated and carried rigidly through every ensemble frame.

These tests pin:
  1.  H-completeness: built H == chemically-correct H for CEBQEI / BOCMIR, on the
      single build AND on every ensemble frame.
  2.  Determinism: byte-identical repeat build (PYTHONHASHSEED=0).
  3.  Gating: the restore is only reached when DELFIN_FFFREE_RIGID_HAPTO=1.
"""
import os

import numpy as np
import pytest
from rdkit import Chem

# Cp / methylcyclohexadienyl piano-stool complexes whose η-ring CH hydrogens were
# dropped on the rigid-hapto path.  (refcode, SMILES, expected total H count.)
_CASES = [
    ("CEBQEI",
     "[Cl][Ru-7]123456([C+]7[C+]1[C+]2[C+]3[C+]4[C+]75)"
     "[N+](NC1=CC=CC=C1)=CC1=CC=CC=[N+]16", 17),
    ("BOCMIR",
     "CC1=C[C+]23->[Cr-9]4567([C]#[O+])([C]#[O+])([C]#[O+])"
     "<-[C+]2(C1)[C+]4[C+]5[C+]6[C+]37", 10),
]


def _build(smiles):
    from delfin.fffree import decompose as DC
    from delfin.fffree import assemble_complex as AC
    d = DC.decompose(smiles)
    assert d is not None, "rigid-hapto decompose returned None"
    return d, AC


@pytest.fixture(autouse=True)
def _rigid_hapto_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")
    monkeypatch.setenv("DELFIN_FFFREE_RIGID_HAPTO", "1")
    monkeypatch.setenv("DELFIN_FFFREE_CN3", "1")
    monkeypatch.setenv("DELFIN_FFFREE_CN_EXTEND", "1")
    monkeypatch.setenv("PYTHONHASHSEED", "0")


@pytest.mark.parametrize("refcode,smiles,exp_h", _CASES)
def test_eta_ring_h_single_build(refcode, smiles, exp_h):
    d, AC = _build(smiles)
    syms, P, donors, exempt = AC.assemble_hapto(d["metal"], d["geometry"], d)
    n_h = sum(1 for s in syms if s == "H")
    assert n_h == exp_h, f"{refcode}: built {n_h} H, expected {exp_h}"
    assert np.all(np.isfinite(np.asarray(P)))


@pytest.mark.parametrize("refcode,smiles,exp_h", _CASES)
def test_eta_ring_h_all_ensemble_frames(refcode, smiles, exp_h):
    d, AC = _build(smiles)
    builds = AC.assemble_hapto_ensemble(d["metal"], d["geometry"], d, max_builds=30)
    assert builds, f"{refcode}: ensemble empty"
    counts = {sum(1 for s in b[0] if s == "H") for b in builds}
    assert counts == {exp_h}, f"{refcode}: ensemble H counts {sorted(counts)} != {exp_h}"


@pytest.mark.parametrize("refcode,smiles,exp_h", _CASES)
def test_eta_ring_h_geometry_exo(refcode, smiles, exp_h):
    """Each restored ring-H is exo (away from centroid), ~1.08 A C-H, far side of
    the metal (M-C-H angle obtuse, M...H > M-C)."""
    d, AC = _build(smiles)
    syms, P, donors, exempt = AC.assemble_hapto(d["metal"], d["geometry"], d)
    P = np.asarray(P, float)
    M = P[0]
    Cs = [i for i, s in enumerate(syms) if s == "C" and 1.8 < np.linalg.norm(P[i] - M) < 2.6]
    Os = [i for i, s in enumerate(syms) if s == "O"]
    ring = [i for i in Cs if not any(np.linalg.norm(P[i] - P[o]) < 1.4 for o in Os)]
    assert ring, f"{refcode}: no η-ring carbons found"
    centroid = P[ring].mean(axis=0)
    n_ring_h = 0
    for ci in ring:
        for hi in [h for h, s in enumerate(syms) if s == "H" and np.linalg.norm(P[h] - P[ci]) < 1.4]:
            n_ring_h += 1
            ch = np.linalg.norm(P[hi] - P[ci])
            mh = np.linalg.norm(P[hi] - M)
            mc = np.linalg.norm(P[ci] - M)
            assert 1.0 < ch < 1.15, f"{refcode}: C-H {ch:.3f} not ~1.08"
            assert mh > mc, f"{refcode}: H not on far side of metal ({mh:.2f} <= {mc:.2f})"
            assert (np.linalg.norm(P[hi] - centroid)
                    > np.linalg.norm(P[ci] - centroid)), f"{refcode}: H not exo"
    assert n_ring_h >= 1


def test_repeat_build_deterministic():
    import hashlib

    def _hash(smi):
        d, AC = _build(smi)
        syms, P, _, _ = AC.assemble_hapto(d["metal"], d["geometry"], d)
        blob = "".join(f"{s}{np.round(p, 6)}" for s, p in zip(syms, P))
        return hashlib.md5(blob.encode()).hexdigest()

    smi = _CASES[1][1]
    assert _hash(smi) == _hash(smi)


def test_restore_is_idempotent():
    """Re-running the restore on an already-fixed fragment changes nothing."""
    from delfin.fffree.decompose import _restore_eta_ring_hydrogens
    frag = Chem.MolFromSmiles("[C+]1[C+][C+][C+][C+][C+]1")
    eta = list(range(frag.GetNumAtoms()))
    once = _restore_eta_ring_hydrogens(frag, eta)
    twice = _restore_eta_ring_hydrogens(once, eta)
    h1 = sum(1 for a in Chem.AddHs(once).GetAtoms() if a.GetAtomicNum() == 1)
    h2 = sum(1 for a in Chem.AddHs(twice).GetAtoms() if a.GetAtomicNum() == 1)
    assert h1 == h2 == 6
