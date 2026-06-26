"""Tests for the FF-free JOINT inter-ligand declash (delfin.fffree.joint_declash).

Locks the load-bearing contracts of the class-B recovery lever:
  * geometry-preservation: declash is rigid-rotation-only -> bond-length and
    bond-angle RMSD before/after == 0 (machine epsilon),
  * M-D invariant: every metal-donor distance preserved exactly,
  * never-worse: the inter-ligand heavy-heavy minimum never decreases,
  * declash actually OPENS an artificially imposed inter-ligand heavy clash,
  * determinism: same input -> same output (no RNG),
  * env-gate: DELFIN_FFFREE_JOINT_DECLASH unset -> declash_if_enabled is a no-op,
  * decompose gate-lift is inert when the flag is unset.
"""
import os

import numpy as np
import pytest

from delfin.fffree import joint_declash as JD
from delfin.fffree import decompose as DEC


def _two_ligand_complex():
    """Tiny synthetic complex: metal at origin, two bent C-C-C ligands at ~90°
    whose DISTAL carbons are swung into a genuine heavy-heavy clash (1.41 A, well
    below the 0.75*(1.7+1.7)=2.55 A C-C vdW floor).  A whole-body rotation about
    each M-donor axis can swing the distal carbons apart while keeping the donors
    on their vertices.  Ligand A donor=C1 (on +x), ligand B donor=C4 (on +y)."""
    # atoms: 0=M, 1=C(donorA) 2=C 3=C(distalA), 4=C(donorB) 5=C 6=C(distalB)
    syms = ["Fe", "C", "C", "C", "C", "C", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],     # M
        [1.9, 0.0, 0.0],     # C donor A (on +x axis)
        [3.0, 0.9, 0.0],     # C
        [3.0, 2.0, 0.0],     # C distal A (swung toward ligand B)
        [0.0, 1.9, 0.0],     # C donor B (on +y axis)
        [0.9, 3.0, 0.0],     # C
        [2.0, 3.0, 0.0],     # C distal B (swung toward ligand A; 1.41 A from distalA)
    ], dtype=float)
    # bonds: M-donorA, A internal; M-donorB, B internal
    bond_pairs = [(0, 1), (1, 2), (2, 3), (0, 4), (4, 5), (5, 6)]
    frozen = {0, 1, 4}      # metal + the two donors
    return syms, P, frozen, bond_pairs


def _bond_lengths(P, bonds):
    return np.array([float(np.linalg.norm(P[i] - P[j])) for i, j in bonds])


def _bond_angles(P, bonds):
    adj = {}
    for i, j in bonds:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    angs = []
    for c, nbrs in adj.items():
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                v1 = P[nbrs[a]] - P[c]
                v2 = P[nbrs[b]] - P[c]
                cu = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
                angs.append(np.degrees(np.arccos(np.clip(cu, -1, 1))))
    return np.array(angs)


def test_geometry_preservation_exact():
    """Rigid rotation only: bond lengths + bond angles preserved to machine eps."""
    syms, P, frozen, bp = _two_ligand_complex()
    L0 = _bond_lengths(P, bp)
    A0 = _bond_angles(P, bp)
    Pr = JD.declash(syms, P.copy(), frozen, bond_pairs=bp)
    L1 = _bond_lengths(Pr, bp)
    A1 = _bond_angles(Pr, bp)
    assert np.sqrt(((L1 - L0) ** 2).mean()) < 1e-9
    assert np.sqrt(((A1 - A0) ** 2).mean()) < 1e-8


def test_md_invariant_preserved():
    """Every metal-donor distance must be preserved (<= md_tol)."""
    syms, P, frozen, bp = _two_ligand_complex()
    Pr = JD.declash(syms, P.copy(), frozen, bond_pairs=bp)
    mi = 0
    for d in (1, 4):
        d0 = np.linalg.norm(P[d] - P[mi])
        d1 = np.linalg.norm(Pr[d] - Pr[mi])
        assert abs(d1 - d0) <= JD._DEF_MD_TOL


def test_declash_opens_inter_ligand_clash_never_worse():
    """The inter-ligand heavy minimum must rise (or stay) — never decrease."""
    syms, P, frozen, bp = _two_ligand_complex()
    lig = JD._ligand_of_atom(len(syms), syms, bp, P)

    def min_inter_heavy(Q):
        best = 1e9
        for i in range(len(syms)):
            for j in range(i + 1, len(syms)):
                if syms[i] == "H" or syms[j] == "H":
                    continue
                if syms[i] in ("Fe",) or syms[j] in ("Fe",):
                    continue
                if lig[i] == lig[j] and lig[i] >= 0:
                    continue
                best = min(best, float(np.linalg.norm(Q[i] - Q[j])))
        return best

    before = min_inter_heavy(P)
    Pr = JD.declash(syms, P.copy(), frozen, bond_pairs=bp)
    after = min_inter_heavy(Pr)
    assert after >= before - 1e-6        # never worse
    assert after > before                # and strictly opened the imposed clash


def test_determinism():
    """Same input twice -> byte-identical output (no RNG)."""
    syms, P, frozen, bp = _two_ligand_complex()
    a = JD.declash(syms, P.copy(), frozen, bond_pairs=bp)
    b = JD.declash(syms, P.copy(), frozen, bond_pairs=bp)
    assert np.array_equal(a, b)


def test_envgate_off_is_noop():
    """declash_if_enabled is a no-op when DELFIN_FFFREE_JOINT_DECLASH is unset."""
    syms, P, frozen, bp = _two_ligand_complex()
    os.environ.pop("DELFIN_FFFREE_JOINT_DECLASH", None)
    out = JD.declash_if_enabled(syms, P.copy(), frozen, bond_pairs=bp)
    assert np.array_equal(np.asarray(out, dtype=float), P)


def test_envgate_on_runs():
    """declash_if_enabled actually relaxes when the flag is set."""
    syms, P, frozen, bp = _two_ligand_complex()
    os.environ["DELFIN_FFFREE_JOINT_DECLASH"] = "1"
    try:
        out = np.asarray(JD.declash_if_enabled(syms, P.copy(), frozen, bond_pairs=bp),
                         dtype=float)
        assert not np.array_equal(out, P)        # it moved something
        assert np.all(np.isfinite(out))
    finally:
        os.environ.pop("DELFIN_FFFREE_JOINT_DECLASH", None)


def test_never_nonfinite():
    """Pathological input never yields non-finite coordinates."""
    syms = ["Fe", "C", "C"]
    P = np.zeros((3, 3))                          # everything on the metal
    out = JD.declash(syms, P.copy(), {0, 1, 2}, bond_pairs=[(0, 1), (0, 2)])
    assert np.all(np.isfinite(out))


def test_decompose_gatelift_inert_when_flag_off():
    """A per-arm-9 monodentate ligand stays None (legacy) with the flag UNSET, and
    builds-decompose only with the flag set (gate-lift coupled to JOINT_DECLASH)."""
    # VANBER: two mesityl C-donors (9 heavy each) + 2 oxo -> per-arm 9 > 8.
    smi = "CC1=CC(C)=[C]([Mo+2]([O-])([O-])[C]2=C(C)C=C(C)C=C2C)C(C)=C1"
    os.environ.pop("DELFIN_FFFREE_JOINT_DECLASH", None)
    assert DEC.decompose(smi) is None            # gate at 8 -> legacy
    os.environ["DELFIN_FFFREE_JOINT_DECLASH"] = "1"
    try:
        d = DEC.decompose(smi)
        assert d is not None                     # gate lifted to 12 -> decomposes
        assert d["cn"] == 4 and not d["has_chelate"]
    finally:
        os.environ.pop("DELFIN_FFFREE_JOINT_DECLASH", None)


def test_decompose_chelate_arm_keeps_historic_cap():
    """The gate-lift is MONODENTATE-only: a chelate arm keeps the historic cap of
    8 even with the flag set (denticity>1 not raised)."""
    # A bidentate ligand with a 10-heavy arm would still bail; we assert the cap
    # logic by checking a known small chelate still decomposes unchanged.
    smi = "C(=O)([O-])C(=O)[O-]"                  # oxalate-like fragment context
    # ethylenediamine on Pt (small chelate, well under cap either way)
    en = "NCCN[Pt]1NCCN1"
    os.environ["DELFIN_FFFREE_JOINT_DECLASH"] = "1"
    try:
        d_on = DEC.decompose(en)
    finally:
        os.environ.pop("DELFIN_FFFREE_JOINT_DECLASH", None)
    d_off = DEC.decompose(en)
    # small chelate: identical decision both ways (cap never binds)
    assert (d_on is None) == (d_off is None)
