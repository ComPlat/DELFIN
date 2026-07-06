#!/usr/bin/env python3
"""_d8_square_planar.py — additive d8-CN4 square-planar SEATING enforcement (FF-free, geometric).

Root cause (eye-flagged, 2026-07-06): for a d8 metal (Pt/Pd/Ni/Au/Rh/Ir) at CN4 the crystal is
ALWAYS square-planar (SP-4), never tetrahedral.  `decompose._default_geometry` already targets SP-4
for d8, and `metal_sphere_builder.place_sphere` seats monodentate donors on SP-4 vertices — BUT for
CHELATE complexes the coordination geometry comes from the chelate's ETKDG embed (tetrahedral bite),
so the donors end up TETRAHEDRAL and the SP-4 preference is silently lost (HUPJUY Pd-Se2P2, VIHFAT
Pt-P2 built T-4; the whole manifold was ALARM).  The existing DELFIN_FFFREE_CN4_BOTH flag only
enumerates geometry LABELS — it does NOT enforce donor COPLANARITY — so it does not fix this.

This module takes a BUILT frame and, for each d8 metal at CN4 that is not already planar, produces an
ADDITIONAL frame where the four donors are moved into a common plane through the metal by a RIGID-BODY
rotation of each coordinated ligand fragment (ligand-internal geometry is preserved exactly; only the
donor's DIRECTION from the metal changes).  Chelate fragments (one fragment carrying >=2 donors) are
rotated as ONE unit by a Kabsch fit of all their donors to their in-plane targets, so the bite is
kept.  Strictly ADDITIVE (the caller keeps the original T-4 frame too) and env-gated default-OFF.

License-clean: covalent radii + numpy only; no CCDC/COD, no force field.
"""
from __future__ import annotations
import os
from typing import List, Optional, Tuple
import numpy as np

# d8 square-planar-preferring metals (matches decompose._D8)
_D8 = {"Pt", "Pd", "Ni", "Au", "Rh", "Ir"}
# Cordero covalent radii (A) for the M-D coordination-distance cutoff (metal-aware; license-clean)
_COV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20, "Br": 1.20, "Te": 1.38, "I": 1.39,
    "Ni": 1.24, "Pd": 1.39, "Pt": 1.36, "Au": 1.36, "Rh": 1.42, "Ir": 1.41,
}


def d8_sp4_enabled() -> bool:
    """Additive d8-CN4 square-planar seating (default OFF -> byte-identical when unset)."""
    return os.environ.get("DELFIN_D8_SP4_SEAT", "0") == "1"


def _cov(e: str) -> float:
    return _COV.get(e, 0.95)


def _donor_cn(syms: List[str], P: np.ndarray, m: int, f: float = 1.30
              ) -> List[int]:
    """Heavy donors coordinated to metal m: within f*(r_cov(M)+r_cov(D)); C only if short (organo-
    metallic), so backbone C at 1.3x+ is not mistaken for a donor (matches metric_coord_shape sense)."""
    out = []
    for j in range(len(syms)):
        if j == m or syms[j] == "H":
            continue
        d = float(np.linalg.norm(P[j] - P[m]))
        cut = f * (_cov(syms[m]) + _cov(syms[j]))
        if syms[j] == "C":
            cut = 1.15 * (_cov(syms[m]) + _cov(syms[j]))   # only a real short M-C bond
        if d < cut:
            out.append(j)
    return out


def _fragments(syms: List[str], P: np.ndarray, m: int, donors: List[int]):
    """Ligand fragment (set of atom indices) attached to each donor, by BFS over geometric bonds
    WITHOUT crossing the metal.  Two donors of a chelate share one fragment."""
    n = len(syms)
    # geometric bonds among non-metal atoms (metal excluded so fragments don't merge through it)
    nbr = {i: [] for i in range(n)}
    for i in range(n):
        if i == m:
            continue
        ri = _cov(syms[i])
        for j in range(i + 1, n):
            if j == m:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * (ri + _cov(syms[j])):
                nbr[i].append(j); nbr[j].append(i)
    frag_of = {}
    frags = []
    for d in donors:
        if d in frag_of:
            continue
        seen = set([d]); stack = [d]
        while stack:
            x = stack.pop()
            for y in nbr[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        fid = len(frags); frags.append(seen)
        for a in seen:
            frag_of[a] = fid
    return frags, frag_of


def _sp4_targets(dirs: np.ndarray) -> np.ndarray:
    """Given 4 current donor unit directions, return 4 COPLANAR target directions.  We only remove
    the OUT-OF-PLANE tilt (project each direction onto the donors' best-fit plane through the metal)
    and keep the IN-PLANE angular structure — including each chelate's natural bite.  Forcing an exact
    90-deg spacing squeezed wide-bite chelates (HUPJUY Se2 116deg) and swung bulky ligands into each
    other; a pure de-tilt is the minimal move that reaches SP-4 CShM without body overlap."""
    U, S, Vt = np.linalg.svd(dirs - dirs.mean(0))
    n = Vt[2] / (np.linalg.norm(Vt[2]) + 1e-12)
    proj = dirs - (dirs @ n)[:, None] * n              # drop the component along the plane normal
    proj = proj / (np.linalg.norm(proj, axis=1, keepdims=True) + 1e-12)
    return proj


def _kabsch(Aw: np.ndarray, Bw: np.ndarray) -> np.ndarray:
    """Proper rotation R minimizing |R·A - B| for point sets through the origin (donors vs targets)."""
    H = Aw.T @ Bw
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def square_planarize_frame(syms: List[str], P: np.ndarray) -> Optional[np.ndarray]:
    """Return a NEW coordinate array where every d8 CN4 metal's donors are coplanar (SP-4), via a
    rigid rotation of each coordinated ligand fragment about the metal.  None if nothing to do."""
    P = np.asarray(P, float)
    metals = [i for i in range(len(syms)) if syms[i] in _D8]
    if not metals:
        return None
    out = P.copy()
    changed = False
    for m in metals:
        donors = _donor_cn(syms, out, m)
        if len(donors) != 4:
            continue
        dirs = np.array([out[d] - out[m] for d in donors])
        dn = dirs / (np.linalg.norm(dirs, axis=1, keepdims=True) + 1e-12)
        # already planar? (max |out-of-plane| small) -> skip
        _U, _S, _Vt = np.linalg.svd(dn - dn.mean(0))
        nrm = _Vt[2]
        if np.max(np.abs(dn @ nrm)) < 0.25:      # ~<14 deg out of plane -> already ~square-planar
            continue
        targ = _sp4_targets(dn)
        frags, frag_of = _fragments(syms, out, m, donors)
        # per fragment: Kabsch-rotate the whole fragment (about the metal) so its donor DIRECTIONS
        # move to their SP-4 targets; chelate fragment fits all its donors at once (keeps the bite).
        for fid, atoms in enumerate(frags):
            fdon = [k for k, d in enumerate(donors) if frag_of.get(d) == fid]
            if not fdon:
                continue
            A = np.array([dn[k] for k in fdon])
            B = np.array([targ[k] for k in fdon])
            R = _kabsch(A, B)
            for a in atoms:
                out[a] = m_center_rotate(out[a], out[m], R)
        changed = True
    return out if changed else None


def m_center_rotate(p: np.ndarray, center: np.ndarray, R: np.ndarray) -> np.ndarray:
    return center + R @ (p - center)


if __name__ == "__main__":   # tiny self-test on a synthetic tetrahedral d8-CN4
    import numpy as _np
    td = _np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    td = td / _np.linalg.norm(td, axis=1, keepdims=True) * 2.0
    syms = ["Pd", "N", "N", "N", "N"]
    P = _np.vstack([[0, 0, 0], td])
    out = square_planarize_frame(syms, P)
    dirs = out[1:] - out[0]
    dn = dirs / _np.linalg.norm(dirs, axis=1, keepdims=True)
    _U, _S, _Vt = _np.linalg.svd(dn - dn.mean(0))
    print("max out-of-plane after:", round(float(_np.max(_np.abs(dn @ _Vt[2]))), 3), "(expect ~0 = planar)")
