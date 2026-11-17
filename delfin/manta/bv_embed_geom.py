"""Chemically right chelate embedding at bond-valence M-D lengths (DELFIN_FFFREE_BV_EMBED_GEOM, default OFF).

THE DEFECT (geometry audit 29.10.2026, root 1).  DELFIN_FFFREE_ML_TARGETS_BV_EMBED embeds every
chelate ring at its bond-valence M-D radii.  To make those radii feasible it (i) swaps every ligand
pair to the bounds of the ligand ALONE, which frees the torsion that closes the chelate ring (the
bipyridine N-C-C'-N'), (ii) resets every metal...non-donor bound to [1.2, 100] A, which erases the
1-3 and 1-4 bounds through the metal (in-plane lone pair, chelate-ring planarity), (iii) keeps only a
95 deg cone at the donor, and (iv) embeds with plain ``EmbedParameters()``, i.e. without ETKDG's
basic knowledge (flat aromatic rings, sp2 impropers) and torsion preferences.  Nothing restores the
geometry afterwards.  Measured: 85 % of the FF-free chelate frames of lc18 carry a crystal-calibrated
defect (pi 72 %, metal out of the donor plane 40 %, donor angle 62 %, hybridisation 66 %; main 46 %);
AFODIJ bipy twisted 79 deg, Ru 58 deg out of the pyridine plane.

THE FIX keeps the BV radii and restores what the ring closure means chemically:

1. The metal row is written from the donor's hybridisation instead of being erased.  For every
   substituent X of a donor D (heavy atoms and H), the 1-3 bound M...X follows from the M-D radius,
   the D-X bond and a window on the angle M-D-X:
     sp2 donor, two substituents (pyridine, imine N): 180 - alpha/2 +- 10 deg, alpha = X1-D-X2 of
       the ligand (lone pair along the external bisector);
     sp3 donor, three substituents: N/O 97.5-121.5 deg, heavier donors (P, As, S) 100-125 deg;
     sp3 donor, two substituents (ether/thioether, water): 95-135 deg;
     sp2 donor, one substituent (C=O, N=): 100-180 deg;  sp donor (nitrile, NCS-N, CO): 150-180 deg;
     one-substituent sp3 donor (alkoxide, thiolate): the historic 95 deg cone.
2. The embed uses ETKDGv3 (experimental torsions + basic knowledge) on that bounds matrix.
3. Every conformer is then refined with MMFF on the LIGAND alone (MMFF has no metal parameters), with
   the metal as a fixed extra point: M-D held at the BV radius, the M-D-X windows above as angle
   constraints, and for sp2 donors with two substituents the metal held in the donor plane by
   torsion constraints M-D-X-Y (0 or 180 deg, whichever side the embed put it on).  Aromatic rings,
   sp2 centres including the donors, and the chelate ring come out planar.  A conformer whose
   refinement fails or moves a bond by more than 10 % keeps its embedded coordinates (never dropped).

Universal: every input is graph-derived (hybridisation from bond orders and aromaticity, angles from
RDKit's own ligand bounds); no element, ligand or system list.  Deterministic (fixed seed, one thread).
Off -> nothing here is called (byte-identical).
"""
from __future__ import annotations

import math
import os
import threading
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

ENV = "DELFIN_FFFREE_BV_EMBED_GEOM"


_SUSPEND = threading.local()


def geom_enabled() -> bool:
    if getattr(_SUSPEND, "on", False):
        return False
    return os.environ.get(ENV, "0") == "1"


def without_geom(fn):
    """fn() with this switch suspended (the historic BV embed, for the fallback of a ring that cannot close)."""
    prev = getattr(_SUSPEND, "on", False)
    _SUSPEND.on = True
    try:
        return fn()
    finally:
        _SUSPEND.on = prev


def conformer_ok(mh, mi: int, P: np.ndarray, donor_idxs: Sequence[int], radii: Sequence[float], bm) -> bool:
    """The embedded ring closed as asked: every M-D within 0.06 A of its target and every ligand bond within
    0.88-1.12 of its bounds-matrix length (KEXHOL: radii 0.19 A off, C-N at 0.86)."""
    try:
        for d, r in zip(donor_idxs, radii):
            if abs(float(np.linalg.norm(P[int(d)] - P[mi])) - float(r)) > 0.06:
                return False
        for b in mh.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if mi in (i, j):
                continue
            d0 = _bm_mid(bm, i, j)
            if d0 <= 0:
                continue
            q = float(np.linalg.norm(P[i] - P[j])) / d0
            if q < 0.88 or q > 1.12:
                return False
        return True
    except Exception:
        return False


def _geom_parts() -> set:
    # Ablation only (DELFIN_FFFREE_BV_EMBED_GEOM_PARTS=mrow,etkdg,refine; unset = all three).  Probe
    # 30.10.2026 on bipy/phen at an octahedral cis pair: plain BV embed metal 2.0-2.2 A out of the
    # pyridine plane, lone pair 76 deg off; ETKDG alone flattens (0.00 A) but leaves the lone pair 13-18
    # deg off; + metal row 10 deg; + refine M-D within 0.01 A of the target.
    v = os.environ.get(ENV + "_PARTS", "")
    return set(v.split(",")) if v else {"mrow", "etkdg", "refine"}


def use_mrow() -> bool:
    return geom_enabled() and "mrow" in _geom_parts()


def use_etkdg() -> bool:
    return geom_enabled() and "etkdg" in _geom_parts()


def use_refine() -> bool:
    return geom_enabled() and "refine" in _geom_parts()


def _bm_mid(bm, i, j) -> float:
    lo, hi = (i, j) if i < j else (j, i)
    return 0.5 * (float(bm[lo][hi]) + float(bm[hi][lo]))


def donor_kind(mh, d: int, mi: int) -> Tuple[str, List[int]]:
    """('sp' | 'sp2' | 'sp3', substituents) of donor ``d`` in ``mh`` (metal ``mi`` excluded).

    Read from the bonds: a triple bond or two double bonds -> sp, an aromatic atom or one double bond
    -> sp2, otherwise sp3.  The dative M-D bond is ignored."""
    a = mh.GetAtomWithIdx(int(d))
    subs = [n.GetIdx() for n in a.GetNeighbors() if n.GetIdx() != mi]
    n_double = n_triple = 0
    for b in a.GetBonds():
        if b.GetOtherAtomIdx(int(d)) == mi:
            continue
        bt = b.GetBondType()
        if bt == Chem.BondType.TRIPLE:
            n_triple += 1
        elif bt == Chem.BondType.DOUBLE:
            n_double += 1
    if n_triple or n_double >= 2:
        return "sp", subs
    if a.GetIsAromatic() or n_double == 1:
        return "sp2", subs
    return "sp3", subs


def angle_windows(mh, bm, d: int, mi: int) -> Dict[int, Tuple[float, float]]:
    """{substituent X: (lo, hi) of the angle M-D-X in degrees} for donor ``d`` (see module docstring)."""
    kind, subs = donor_kind(mh, d, mi)
    if not subs:
        return {}
    el = mh.GetAtomWithIdx(int(d)).GetAtomicNum()
    if kind == "sp":
        return {x: (150.0, 180.0) for x in subs}
    if kind == "sp2" and len(subs) == 2:
        x1, x2 = subs
        d1, d2, d12 = _bm_mid(bm, d, x1), _bm_mid(bm, d, x2), _bm_mid(bm, x1, x2)
        c = (d1 * d1 + d2 * d2 - d12 * d12) / max(2.0 * d1 * d2, 1e-6)
        alpha = math.degrees(math.acos(max(-1.0, min(1.0, c))))
        ideal = 180.0 - 0.5 * alpha
        return {x: (ideal - 10.0, min(ideal + 10.0, 180.0)) for x in subs}
    if kind == "sp2" and len(subs) == 1:
        return {subs[0]: (100.0, 180.0)}
    if kind == "sp2":
        return {x: (95.0, 180.0) for x in subs}
    if len(subs) >= 3:
        if el in (7, 8):
            return {x: (97.5, 121.5) for x in subs}
        return {x: (100.0, 125.0) for x in subs}
    if len(subs) == 2:
        return {x: (95.0, 135.0) for x in subs}
    return {subs[0]: (95.0, 180.0)}


def write_metal_row(bm, mh, mi: int, donor_idxs: Sequence[int], radii: Sequence[float]) -> None:
    """Set the M...X 1-3 bounds of every donor substituent from its angle window (in place)."""
    dset = {int(x) for x in donor_idxs}
    for d, r in zip(donor_idxs, radii):
        d = int(d)
        for x, (lo_a, hi_a) in angle_windows(mh, bm, d, mi).items():
            if x in dset:
                continue
            dx = _bm_mid(bm, d, x)
            lo = math.sqrt(max(r * r + dx * dx - 2.0 * r * dx * math.cos(math.radians(lo_a)), 0.0))
            hi = math.sqrt(max(r * r + dx * dx - 2.0 * r * dx * math.cos(math.radians(hi_a)), 0.0))
            a, b = min(x, mi), max(x, mi)
            bm[b][a] = float(lo)
            bm[a][b] = float(max(hi, lo + 0.02))


def embed_params(bm, seed: int):
    from rdkit.Chem import rdDistGeom as _DG
    ep = _DG.ETKDGv3()
    ep.randomSeed = int(seed)
    ep.useRandomCoords = True
    ep.numThreads = 1
    ep.SetBoundsMat(bm)
    return ep


def _torsion_deg(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - (b0 @ b1n) * b1n
    w = b2 - (b2 @ b1n) * b1n
    return float(np.degrees(np.arctan2(np.cross(b1n, v) @ w, v @ w)))


def refine_conformer(lmol, mh, mi: int, P: np.ndarray, donor_idxs: Sequence[int],
                     radii: Sequence[float], bm) -> Optional[np.ndarray]:
    """MMFF refinement of one embedded conformer of ``mh`` (metal at index ``mi``), metal as a fixed
    extra point.  Returns the refined positions of ALL atoms of ``mh`` (metal unchanged) or None."""
    try:
        lh = Chem.AddHs(lmol)
        nl = lh.GetNumAtoms()
        if mh.GetNumAtoms() != nl + 1:
            return None
        idx = [i if i < mi else i + 1 for i in range(nl)]           # ligand atom -> mh atom
        for i in range(nl):
            if lh.GetAtomWithIdx(i).GetAtomicNum() != mh.GetAtomWithIdx(idx[i]).GetAtomicNum():
                return None
        inv = {m: i for i, m in enumerate(idx)}
        conf = Chem.Conformer(nl)
        for i in range(nl):
            conf.SetAtomPosition(i, P[idx[i]].tolist())
        lh.RemoveAllConformers()
        lh.AddConformer(conf, assignId=True)
        mp = AllChem.MMFFGetMoleculeProperties(lh)
        if mp is None:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(lh, mp, ignoreInterfragInteractions=False)
        if ff is None:
            return None
        pm = P[mi]
        ff.AddExtraPoint(float(pm[0]), float(pm[1]), float(pm[2]), True)
        ff.Initialize()
        m_ff = nl                                                   # index of the metal point
        dset = {int(x) for x in donor_idxs}
        for d, r in zip(donor_idxs, radii):
            ld = inv[int(d)]
            ff.MMFFAddDistanceConstraint(ld, m_ff, False, float(r) - 0.01, float(r) + 0.01, 500.0)
            for x, (lo_a, hi_a) in angle_windows(mh, bm, int(d), mi).items():
                if x in dset:
                    continue
                ff.MMFFAddAngleConstraint(m_ff, ld, inv[x], False, float(lo_a), float(hi_a), 200.0)
            kind, subs = donor_kind(mh, int(d), mi)
            if kind == "sp2" and len(subs) == 2:
                for x in subs:
                    for y in mh.GetAtomWithIdx(x).GetNeighbors():
                        yi = y.GetIdx()
                        if yi in (int(d), mi):
                            continue
                        t = _torsion_deg(P[mi], P[int(d)], P[x], P[yi])
                        t0 = 0.0 if abs(t) < 90.0 else 180.0
                        ff.MMFFAddTorsionConstraint(m_ff, ld, inv[x], inv[yi], False,
                                                    t0 - 2.0, t0 + 2.0, 100.0)
        ff.Minimize(maxIts=2000)
        pos = np.asarray(ff.Positions(), float).reshape(-1, 3)
        Q = np.array(P, float)
        for i in range(nl):
            Q[idx[i]] = pos[i]
        if not np.all(np.isfinite(Q)):
            return None
        for b in lh.GetBonds():                 # never worse on the bonds
            i, j = idx[b.GetBeginAtomIdx()], idx[b.GetEndAtomIdx()]
            d0 = np.linalg.norm(P[i] - P[j])
            d1 = np.linalg.norm(Q[i] - Q[j])
            if abs(d1 - d0) > 0.10 * max(d0, 1e-6) + 0.05:
                return None
        return Q
    except Exception:
        return None
