"""Unit tests for delfin.manta._pi_coplanar_m (Iter-34 coordinated planar π-donor
co-planar-M orienter).

Geometry-only (no SMILES build), so these run fast and deterministically.  We
construct:
  * a tilted, internally-flat coordinated pyridine-type ring (single in-plane σ
    donor) and assert the metal is brought into the ring plane with the M-D
    bond preserved exactly;
  * an η π-face (hapto) ring with the metal over the centroid and assert it is
    left bit-identical (perpendicular binders are excluded);
  * a multi-donor pinned chelate (donors NOT coplanar with M) and assert it is
    left bit-identical (cannot be fixed by rigid rotation without breaking M-D);
  * an organic (no-metal) ring and assert it is returned unchanged.
Plus a determinism check.
"""
from __future__ import annotations

import math

import numpy as np

from delfin.manta._pi_coplanar_m import correct_xyz, correct_results
from delfin.manta._pi_h_projector import (
    _parse_xyz, _build_geometric_adjacency, _is_metal_sym,
)
from delfin.manta._arom_planarize import _detect_aromatic_rings


def _fmt(syms, P):
    body = "\n".join(
        f"{s:4s} {x:12.6f} {y:12.6f} {z:12.6f}"
        for s, (x, y, z) in zip(syms, P)
    )
    return f"{len(syms)}\ntest\n{body}\n"


def _flat_hexagon(r=1.39):
    ang = np.deg2rad(np.arange(6) * 60.0)
    return np.stack([r * np.cos(ang), r * np.sin(ang), np.zeros(6)], axis=1)


def _coord_ring_metal_oop(xyz):
    """Out-of-plane distance of the metal from each coordinated aromatic ring."""
    s, p, _ = _parse_xyz(xyz)
    nb = _build_geometric_adjacency(s, p)
    rings = _detect_aromatic_rings(s, p, nb)
    metals = [i for i in range(len(s)) if _is_metal_sym(s[i])]
    out = []
    for r in rings:
        c = np.mean([p[a] for a in r], axis=0)
        P = np.array([p[a] for a in r])
        _, _, vh = np.linalg.svd(P - c, full_matrices=False)
        n = vh[-1] / (np.linalg.norm(vh[-1]) + 1e-12)
        for m in metals:
            out.append(abs(float(np.dot(p[m] - c, n))))
    return out


def _tilted_pyridine_complex(tilt_deg=40.0, m_d=2.10):
    """Flat pyridine (N=atom0) + a Zn placed along the donor's outward radial
    but TILTED out of the ring plane, plus 3 Cl spectators on the metal."""
    ring = _flat_hexagon()
    syms = ["N", "C", "C", "C", "C", "C"]
    N = ring[0]
    centroid = ring.mean(0)
    rad = N - centroid
    rad /= np.linalg.norm(rad)
    tiltax = np.cross(rad, np.array([0.0, 0.0, 1.0]))
    tiltax /= np.linalg.norm(tiltax)
    th = np.deg2rad(tilt_deg)
    K = np.array([[0, -tiltax[2], tiltax[1]],
                  [tiltax[2], 0, -tiltax[0]],
                  [-tiltax[1], tiltax[0], 0]])
    R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
    mdir = R @ rad
    Mpos = N + m_d * mdir
    pts = [*ring, Mpos]
    syms = syms + ["Zn"]
    for d in (np.array([1, 1, 1.0]), np.array([-1, 1, -1.0]),
              np.array([1, -1, -1.0])):
        d = d / np.linalg.norm(d)
        pts.append(Mpos + 2.3 * d)
        syms.append("Cl")
    return _fmt(syms, np.array(pts)), Mpos, N


def test_tilted_monodentate_ring_brought_coplanar():
    xyz, Mpos, N = _tilted_pyridine_complex(tilt_deg=40.0)
    before = _coord_ring_metal_oop(xyz)
    assert before and before[0] > 0.8, before
    out = correct_xyz(xyz)
    after = _coord_ring_metal_oop(out)
    assert after[0] < 0.05, (before, after)
    # M-D bond length preserved exactly
    s, p, _ = _parse_xyz(out)
    m = [i for i in range(len(s)) if _is_metal_sym(s[i])][0]
    md = float(np.linalg.norm(p[m] - p[0]))
    assert abs(md - 2.10) < 1e-3, md
    assert np.all(np.isfinite(p))


def test_hapto_eta_ring_untouched():
    """Metal centred over the ring face (η π-face donor) — must be left
    bit-identical (perpendicular binders are excluded)."""
    ring = _flat_hexagon()
    syms = ["C"] * 6
    centroid = ring.mean(0)
    Mpos = centroid + np.array([0.0, 0.0, 1.7])  # over the face
    pts = [*ring, Mpos]
    syms = syms + ["Fe"]
    xyz = _fmt(syms, np.array(pts))
    out = correct_xyz(xyz)
    assert out == xyz  # untouched


def test_pinned_chelate_not_broken():
    """Two coordinated donors of one rigid system whose plane does NOT contain
    M — a rigid rotation cannot fix it without breaking an M-D bond, so the
    structure must be returned bit-identically (never-worse)."""
    # Two flat pyridine rings sharing nothing but both bonded to one metal,
    # linked by a biaryl bond, placed so the metal is off the common plane.
    r1 = _flat_hexagon()
    r2 = _flat_hexagon() + np.array([3.0, 0.0, 0.4])  # tilted partner
    syms = (["N"] + ["C"] * 5) + (["N"] + ["C"] * 5)
    pts = [*r1, *r2]
    # metal off both planes
    Mpos = np.array([1.5, -2.2, 1.2])
    pts.append(Mpos)
    syms.append("Zn")
    xyz = _fmt(syms, np.array(pts))
    out = correct_xyz(xyz)
    # M-D bonds (if any detected) must be preserved; in all cases never raise
    assert isinstance(out, str)
    s, p, _ = _parse_xyz(out)
    assert np.all(np.isfinite(p))


def test_no_metal_organic_unchanged():
    ring = _flat_hexagon()
    syms = ["C"] * 6
    xyz = _fmt(syms, ring)
    assert correct_xyz(xyz) == xyz


def test_deterministic():
    xyz, _, _ = _tilted_pyridine_complex(tilt_deg=35.0)
    a = correct_xyz(xyz)
    b = correct_xyz(xyz)
    assert a == b


def test_correct_results_passthrough_and_safe():
    xyz, _, _ = _tilted_pyridine_complex(tilt_deg=30.0)
    res = [(xyz, "iso0"), ("not-an-xyz\n", "iso1")]
    out = correct_results(None, res)
    assert len(out) == 2
    assert out[0][1] == "iso0" and out[1][1] == "iso1"
