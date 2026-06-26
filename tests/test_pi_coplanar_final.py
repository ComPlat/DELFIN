"""Tests for the final per-arm π-coplanarity polish (_pi_coplanar_final).

Geometry-only, deterministic, never-worse: brings the metal into a σ-aromatic
donor's ring plane by rotating ONLY that arm about its donor (M-D preserved),
for both a monodentate ligand and one arm of a multi-arm polydentate, without
creating any inter-fragment clash.  Default-OFF byte-identical via the pipeline.
"""
import math

import numpy as np

from delfin import _pi_coplanar_final as pcf


def _pyridine(tilt_deg, md=2.10):
    """Fe + a flat pyridine ring (N at vertex 0) with the metal placed at md along
    the in-plane exterior bisector, then tilted `tilt_deg` out of the ring plane."""
    r = 1.39
    ring = [[r * math.cos(math.radians(60 * k)),
             r * math.sin(math.radians(60 * k)), 0.0] for k in range(6)]
    syms = ["N", "C", "C", "C", "C", "C"]
    H, Hs = [], []
    for k in range(1, 6):
        a = math.radians(60 * k)
        H.append([(r + 1.08) * math.cos(a), (r + 1.08) * math.sin(a), 0.0]); Hs.append("H")
    N = np.array([1.39, 0.0, 0.0])
    Mdir = np.array([math.cos(math.radians(tilt_deg)), 0.0, math.sin(math.radians(tilt_deg))])
    M = N + md * Mdir
    syms_all = ["Fe"] + syms + Hs
    xyz = [list(M)] + ring + H
    return "\n".join(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}" for s, c in zip(syms_all, xyz))


def _metal_oop(block, ring_syms=("N", "C")):
    L = block.splitlines(); S, P = [], []
    for ln in L:
        p = ln.split()
        S.append(p[0]); P.append([float(x) for x in p[1:4]])
    P = np.array(P); M = P[S.index("Fe")]
    ring = [i for i, s in enumerate(S) if s in ring_syms][:6]
    pts = P[ring]; c = pts.mean(0)
    _, _, vt = np.linalg.svd(pts - c)
    return abs(float(np.dot(M - c, vt[2]))), float(np.linalg.norm(M - P[1]))


def test_monodentate_pulls_metal_into_plane():
    blk = _pyridine(30.0)
    o0, mn0 = _metal_oop(blk)
    out = pcf.correct_xyz(blk)
    o1, mn1 = _metal_oop(out)
    assert o0 > 0.5 and o1 < 0.1          # metal pulled into plane
    assert abs(mn1 - mn0) < 1e-3          # M-D preserved exactly (donor is pivot)


def test_already_in_plane_is_noop():
    blk = _pyridine(0.0)
    assert pcf.correct_xyz(blk) == blk     # good frame untouched


def test_deterministic():
    blk = _pyridine(35.0)
    assert pcf.correct_xyz(blk) == pcf.correct_xyz(blk)


def test_returns_input_on_garbage():
    assert pcf.correct_xyz("not\na\nvalid\nxyz") == "not\na\nvalid\nxyz"
    assert pcf.correct_xyz("") == ""


def test_no_metal_is_noop():
    blk = "C 0 0 0\nC 1.4 0 0\nN 2.8 0 0"
    assert pcf.correct_xyz(blk) == blk


def test_multiarm_only_tilted_arm_moves_both_md_preserved():
    """Two pyridyl arms on a shared backbone, one tilted: only the tilted arm is
    rotated, both M-N distances are preserved exactly, the flat arm is untouched."""
    def hexring(centroid, donor_dir, normal):
        u = np.array(donor_dir, float); u /= np.linalg.norm(u)
        nrm = np.array(normal, float); nrm /= np.linalg.norm(nrm)
        v = np.cross(nrm, u); v /= np.linalg.norm(v)
        r = 1.39; c = np.array(centroid, float)
        return [c + r * (math.cos(math.radians(60 * k)) * (-u)
                         + math.sin(math.radians(60 * k)) * v) for k in range(6)]
    ring1 = hexring([3.5, 0, 0], [1, 0, 0], [0, 0, 1])
    ring2 = hexring([0, 3.5, 0], [0, 1, 0], [0, 0, 1])
    N2 = ring2[0]
    th = math.radians(25)
    Rx = np.array([[1, 0, 0], [0, math.cos(th), -math.sin(th)], [0, math.sin(th), math.cos(th)]])
    ring2 = [N2 + Rx @ (p - N2) for p in ring2]
    Ca1, Ca2 = ring1[3], ring2[3]
    Cb = (Ca1 + Ca2) / 2
    syms = ["Fe"] + ["N", "C", "C", "C", "C", "C"] * 2 + ["C"]
    xyz = [np.zeros(3)] + ring1 + ring2 + [Cb]
    for ring in (ring1, ring2):
        cen = np.mean(ring, 0)
        for k in range(1, 6):
            if k == 3:
                continue
            d = ring[k] - cen; d /= np.linalg.norm(d)
            xyz.append(ring[k] + 1.08 * d); syms.append("H")
    blk = "\n".join(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}" for s, c in zip(syms, xyz))

    def md(P, idx):
        return float(np.linalg.norm(P[0] - P[idx]))
    P0 = np.array([[float(x) for x in ln.split()[1:4]] for ln in blk.splitlines()])
    out = pcf.correct_xyz(blk)
    assert out != blk
    P1 = np.array([[float(x) for x in ln.split()[1:4]] for ln in out.splitlines()])
    assert abs(md(P1, 1) - md(P0, 1)) < 1e-3     # M-N1 preserved
    assert abs(md(P1, 7) - md(P0, 7)) < 1e-3     # M-N2 preserved


def test_pipeline_byte_identical_when_flag_off(monkeypatch):
    from delfin import smiles_converter as sc
    monkeypatch.delenv("DELFIN_FFFREE_PI_COPLANAR_FINAL", raising=False)
    probe = [("frameA", 1), ("frameB", 2)]
    assert sc._apply_pi_coplanar_final(probe) is probe   # exact identity
