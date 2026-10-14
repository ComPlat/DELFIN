"""DELFIN_FFFREE_LONE_PAIR_SPACE / _LONE_PAIR_MD: lone pairs on free tetrahedral vertices."""
import math

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from delfin.manta import _lone_pair_space as L
from delfin.manta import assemble_orient as AO
from delfin.manta import polyhedra as PH


def _lps_angle(a, b):
    c = float(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _lps_free_ligand(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=42)
    return m, np.array(m.GetConformer().GetPositions())


def test_two_substituents_free_vertices_are_tetrahedral():
    t = math.radians(109.4712)
    u0 = np.array([1.0, 0.0, 0.0])
    u1 = np.array([math.cos(t), math.sin(t), 0.0])
    fv = L.free_vertices([u0, u1])
    assert len(fv) == 2
    for v in fv:
        assert abs(_lps_angle(v, u0) - 109.47) < 0.05 and abs(_lps_angle(v, u1) - 109.47) < 0.05
    assert abs(_lps_angle(fv[0], fv[1]) - 109.47) < 0.05
    assert abs(float(np.dot(fv[0], np.cross(u0, u1)))) > 0.5      # out of the substituent plane


def test_water_metal_direction_is_pyramidal_and_switch_gates_it(monkeypatch):
    m, P = _lps_free_ligand("O")
    o = 0
    hs = [a.GetIdx() for a in m.GetAtomWithIdx(o).GetNeighbors()]
    n = np.cross(P[hs[0]] - P[o], P[hs[1]] - P[o])
    monkeypatch.delenv(L.EARLY_FLAG, raising=False)
    old = AO._donor_and_lp(None, P, m, o)
    assert abs(float(np.dot(old, n))) < 1e-6          # old rule: in the H-O-H plane
    monkeypatch.setenv(L.EARLY_FLAG, "1")
    new = AO._donor_and_lp(None, P, m, o)
    s = (_lps_angle(new, P[hs[0]] - P[o]) + _lps_angle(new, P[hs[1]] - P[o])
         + _lps_angle(P[hs[0]] - P[o], P[hs[1]] - P[o]))
    assert s < 345.0                                    # pyramidal at O


def test_scope():
    for smi, idx, want in [("O", 0, True), ("CO", 1, True), ("CSC", 1, True), ("COC", 1, False),
                           ("Oc1ccccc1", 0, False), ("CC(=O)O", 3, False), ("C1CCOC1", 4, False)]:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        assert L.is_lone_pair_donor(m.GetAtomWithIdx(idx)) is want, smi


def test_ionic_md_only_with_switch(monkeypatch):
    m = Chem.AddHs(Chem.MolFromSmiles("O"))
    o = m.GetAtomWithIdx(0)
    monkeypatch.delenv(L.LONE_PAIR_MD_FLAG + "_EARLY", raising=False)
    base = PH.md_distance("Zn", "O", atom=o, mol=m)
    L.set_current_os(2)
    assert PH.md_distance("Zn", "O", atom=o, mol=m) == base
    monkeypatch.setenv(L.LONE_PAIR_MD_FLAG + "_EARLY", "1")
    assert abs(PH.md_distance("Zn", "O", atom=o, mol=m, cn=5) - 2.04) < 1e-9
    L.set_current_os(None)
    assert PH.md_distance("Zn", "O", atom=o, mol=m, cn=5) == base   # unknown state -> old sum


def _lps_planar_aqua_frame():
    # Zn with a planar water (angle sum 360) and an ammonia-like N donor trans to it
    syms = ["Zn", "O", "H", "H", "N", "H", "H", "H"]
    a = math.radians(104.5 / 2)
    P = np.array([[0, 0, 0], [1.88, 0, 0],
                  [1.88 + 0.97 * math.cos(a), 0.97 * math.sin(a), 0],
                  [1.88 + 0.97 * math.cos(a), -0.97 * math.sin(a), 0],
                  [-2.0, 0, 0], [-2.37, 0.95, 0], [-2.37, -0.47, 0.82], [-2.37, -0.47, -0.82]])
    return syms, P


def test_seat_frame_pyramidalises_planar_water(monkeypatch):
    syms, P = _lps_planar_aqua_frame()
    Q, k, kept = L.seat_frame(syms, P)
    assert k == 1 and kept == 0
    s = (_lps_angle(Q[0] - Q[1], Q[2] - Q[1]) + _lps_angle(Q[0] - Q[1], Q[3] - Q[1])
         + _lps_angle(Q[2] - Q[1], Q[3] - Q[1]))
    assert s < 345.0
    assert abs(float(np.linalg.norm(Q[1] - Q[0])) - 1.88) < 1e-9       # M-D unchanged
    assert abs(_lps_angle(Q[2] - Q[1], Q[3] - Q[1]) - 109.47) < 1.5       # H on its own vertex
    assert np.allclose(Q[4:], P[4:])                                     # others untouched
    Q2, _, _ = L.seat_frame(syms, P, ox=2)
    assert np.allclose(Q2, Q)                                            # MD switch off
    monkeypatch.setenv(L.LONE_PAIR_MD_FLAG, "1")
    Q3, _, _ = L.seat_frame(syms, P, ox=2)
    assert abs(float(np.linalg.norm(Q3[1] - Q3[0])) - (0.60 + 1.36)) < 1e-6   # CN 2 -> nearest tabulated CN 4


def test_seat_results_off_returns_same_object(monkeypatch):
    monkeypatch.delenv(L.FLAG, raising=False)
    r = [("x", "y")]
    assert L.seat_results(r, "O") is r


def test_seat_never_trades_one_short_contact_for_another():
    """#525 VI (EZIBOF): a better minimum may not hide a NEW contact below the floor."""
    syms = ["H", "C", "C"]
    # before: H-C1 1.70 A (ratio 0.59, the old minimum), H-C2 2.10 A (0.72, just clear)
    P0 = np.array([[0.0, 0.0, 0.0], [1.70, 0.0, 0.0], [-2.10, 0.0, 0.0]])
    # after: H-C1 1.75 A (minimum improves), H-C2 1.90 A (0.655 -> a new short contact)
    P1 = np.array([[0.0, 0.0, 0.0], [1.75, 0.0, 0.0], [-1.90, 0.0, 0.0]])
    assert L._contact_worsened(P0, P1, syms, [0], [1, 2])
    # a move that only lengthens contacts is allowed
    P2 = np.array([[0.0, 0.0, 0.0], [1.75, 0.0, 0.0], [-2.20, 0.0, 0.0]])
    assert not L._contact_worsened(P0, P2, syms, [0], [1, 2])
    # a polar H against an acceptor is a hydrogen bond, guarded elsewhere
    assert not L._contact_worsened(P0, P1, ["H", "C", "O"], [0], [1, 2], polar_h={0})
