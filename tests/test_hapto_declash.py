"""Safety + efficacy tests for delfin.manta.hapto_declash (eta-ring centroid spin).

Proves the invariants the morning eye-gate relies on:
  * every M-(ring atom) distance is preserved (coordination polyhedron invariant),
  * the inter-ligand overlap is never increased (never-worse),
  * the relaxation is deterministic (identical output on repeat),
  * a genuine eclipsing clash IS relieved (efficacy).
"""
import math

import numpy as np

from delfin.manta import hapto_declash as HD


def _total(syms, P):
    """The module's partition-free guarded measure (heavy-heavy, H zeroed)."""
    P = np.asarray(P, float)
    vdw = np.array([(HD._VDW.get(s, HD._VDW_D) if s != "H" else 0.0) for s in syms])
    cov = np.array([HD._COV.get(s, HD._COV_D) for s in syms])
    return HD._total_overlap(P, vdw, cov, HD._DEF_CLASH_F)


def _ferrocene_like_piano_stool():
    """A synthetic (eta6-C6H5Cl)Cr(CO)3 piano stool with the ring Cl ECLIPSING a
    carbonyl leg.

    Cr at origin; benzene ring (6 C, radius 1.42) at z = +1.6 (centroid on the +z
    axis); a Cl substituent on ring-C0 pointing radially out + down so it sits
    directly over the CO leg at the same azimuth (angle 0); three CO legs below at
    0/120/240 degrees.  A 60-degree ring spin (one benzene step) carries the Cl to
    the gap between legs and relieves the clash -- the staggered conformer the
    analytical seat does not choose.
    """
    syms = ["Cr"]
    P = [[0.0, 0.0, 0.0]]
    bonds = []
    ring_idx = []
    rC = 1.42
    zC = 1.60
    for k in range(6):
        ang = 2 * math.pi * k / 6.0
        idx = len(syms)
        syms.append("C")
        P.append([rC * math.cos(ang), rC * math.sin(ang), zC])
        ring_idx.append(idx)
        bonds.append((0, idx))  # M-C(ring) eta bonds
    for k in range(6):
        bonds.append((ring_idx[k], ring_idx[(k + 1) % 6]))
    # H on ring C1..C5 (C0 carries Cl)
    for k in range(1, 6):
        ang = 2 * math.pi * k / 6.0
        idx = len(syms)
        syms.append("H")
        P.append([2.1 * math.cos(ang), 2.1 * math.sin(ang), zC + 0.05])
        bonds.append((ring_idx[k], idx))
    # Cl substituent on C0, pointing radially out (angle 0) and DOWN toward a leg
    cl = len(syms)
    syms.append("Cl")
    P.append([2.60, 0.0, 0.90])
    bonds.append((ring_idx[0], cl))
    # three CO legs at 0/120/240 deg below the metal; the one at angle 0 is
    # eclipsed by the Cl
    co_c = []
    for kk in range(3):
        ang = 2 * math.pi * kk / 3.0
        ci = len(syms)
        syms.append("C")
        P.append([1.5 * math.cos(ang), 1.5 * math.sin(ang), -1.10])
        bonds.append((0, ci))
        oi = len(syms)
        syms.append("O")
        P.append([2.4 * math.cos(ang), 2.4 * math.sin(ang), -1.75])
        bonds.append((ci, oi))
        co_c.append(ci)
    metal_indices = [0]
    by_metal = {0: [ring_idx]}
    all_hapto = list(ring_idx)
    frozen_donors = list(co_c)
    return syms, np.array(P, float), metal_indices, by_metal, all_hapto, frozen_donors, bonds


def test_mc_distance_invariant_and_never_worse():
    syms, P, mi, bm, hap, fd, bonds = _ferrocene_like_piano_stool()
    md0 = {c: float(np.linalg.norm(P[c] - P[0])) for c in hap}
    Pn = HD.declash_hapto(P, syms, mi, bm, hap, fd, bonds)
    # M-C(ring) distances preserved to < 1e-3 A
    for c in hap:
        md1 = float(np.linalg.norm(Pn[c] - Pn[0]))
        assert abs(md1 - md0[c]) < 1e-3, f"M-C{c} moved {abs(md1-md0[c]):.2e}"
    # metal itself fixed
    assert np.allclose(Pn[0], P[0], atol=1e-9)
    # partition-free total overlap (the module's guarded measure) never increased
    l0, w0, h0 = _total(syms, P)
    l1, w1, h1 = _total(syms, Pn)
    assert l1 <= l0 + 1e-9 and w1 >= w0 - 1e-6 and h1 <= h0


def test_deterministic():
    syms, P, mi, bm, hap, fd, bonds = _ferrocene_like_piano_stool()
    a = HD.declash_hapto(P, syms, mi, bm, hap, fd, bonds)
    b = HD.declash_hapto(P, syms, mi, bm, hap, fd, bonds)
    assert np.array_equal(a, b)


def test_relieves_a_real_eclipse():
    """With a constructed eclipsing clash present, the spin must strictly lower the
    total overlap."""
    syms, P, mi, bm, hap, fd, bonds = _ferrocene_like_piano_stool()
    l0, _, _ = _total(syms, P)
    Pn = HD.declash_hapto(P, syms, mi, bm, hap, fd, bonds)
    l1, _, _ = _total(syms, Pn)
    assert l0 > 0.0, "fixture has no clash to relieve"
    assert l1 < l0, f"overlap not relieved: {l0:.3f} -> {l1:.3f}"


def test_geometry_detection_and_frame_declash():
    """declash_frame must detect the eta-ring from geometry alone, preserve the
    M-C distances, and never increase the inter-ligand overlap."""
    syms, P, mi, bm, hap, fd, bonds = _ferrocene_like_piano_stool()
    det = HD.detect_hapto(syms, P)
    assert det is not None, "eta-ring not detected from geometry"
    metals, by_metal, all_hapto, frozen_donors, _ = det
    assert metals == [0]
    assert len(all_hapto) >= 4
    md0 = {c: float(np.linalg.norm(P[c] - P[0])) for c in all_hapto}
    Pn = HD.declash_frame(syms, P)
    for c in all_hapto:
        assert abs(float(np.linalg.norm(Pn[c] - Pn[0])) - md0[c]) < 1e-3
    # never-worse on the partition-free measure
    l0, w0, h0 = _total(syms, P)
    l1, w1, h1 = _total(syms, Pn)
    assert l1 <= l0 + 1e-9 and w1 >= w0 - 1e-6 and h1 <= h0


def test_declash_frame_no_metal_is_identity():
    """A pure-organic frame (no metal) is returned unchanged."""
    syms = ["C", "C", "O", "H", "H"]
    P = np.array([[0, 0, 0], [1.5, 0, 0], [2.2, 0.5, 0], [-0.5, 0.9, 0], [-0.5, -0.9, 0]], float)
    Pn = HD.declash_frame(syms, P)
    assert np.array_equal(Pn, P)


def test_carbonyl_contraction():
    """A stretched M-C#O (C-O ~1.43) is contracted toward 1.15; only the terminal
    O moves (inward along C->O), C/metal fixed; a good carbonyl is left alone."""
    # Cr-C-O linear along x; C at 1.9 (M-C), O stretched to 1.9+1.43
    syms = ["Cr", "C", "O"]
    P = np.array([[0, 0, 0], [1.9, 0, 0], [1.9 + 1.43, 0, 0]], float)
    Pn = HD.correct_carbonyls(syms, P)
    co = float(np.linalg.norm(Pn[2] - Pn[1]))
    assert abs(co - HD._CO_TRIPLE) < 1e-6, f"C-O not contracted: {co:.3f}"
    assert np.allclose(Pn[0], P[0]) and np.allclose(Pn[1], P[1]), "metal/C moved"
    # O moved inward (closer to C), not outward
    assert np.linalg.norm(Pn[2] - Pn[0]) < np.linalg.norm(P[2] - P[0])
    # a correct carbonyl (already 1.15) is untouched
    P2 = np.array([[0, 0, 0], [1.9, 0, 0], [1.9 + 1.15, 0, 0]], float)
    assert np.array_equal(HD.correct_carbonyls(syms, P2), P2)


def test_carbonyl_no_metal_or_nonterminal_is_identity():
    # no metal -> identity
    syms = ["C", "C", "O"]
    P = np.array([[0, 0, 0], [1.5, 0, 0], [2.9, 0, 0]], float)
    assert np.array_equal(HD.correct_carbonyls(syms, P), P)
    # bridging O (2 heavy neighbours) -> not a terminal carbonyl -> identity
    syms2 = ["Cr", "C", "O", "C"]
    P2 = np.array([[0, 0, 0], [1.9, 0, 0], [3.3, 0, 0], [4.7, 0, 0]], float)
    assert np.array_equal(HD.correct_carbonyls(syms2, P2), P2)


def test_chelate_through_face_and_sidearm_is_skipped():
    """A ring whose ligand component also carries a sigma-donor is NOT spun
    (spinning would break the sidearm M-D bond)."""
    syms, P, mi, bm, hap, fd, bonds = _ferrocene_like_piano_stool()
    # tether the aryl sigma-donor (a frozen donor) into the Cp ring's component
    bonds = list(bonds) + [(hap[0], fd[-1])]  # ring C -- aryl ipso C
    Pn = HD.declash_hapto(P, syms, mi, bm, hap, fd, bonds)
    # ring component now carries a frozen donor -> skipped -> identity
    assert np.array_equal(Pn, np.array(P, float))
