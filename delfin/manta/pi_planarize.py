"""Conjugated-pi-ligand planarizer (DELFIN_FFFREE_PI_PLANARIZE, default-OFF).

RECOVERS a champion advantage main lost: champion ddde3273 / V2R sent conjugated
metallacycles through the legacy Open-Babel-UFF whole-complex relax, whose aromatic
planarity terms FLATTENED the conjugated pi-macrocycle AND pulled the metal into that
plane (DECTEN: V2R macroRMS 0.25 / metal-in-plane 0.00 vs the native MANTA builder's
0.98 / 0.21).  The native builder (rigid seat + frozen-donor relaxers) never planarises
the fused pi-system, so the ETKDG embed's random pucker/saddle survives into the emitted
frame (user-eye 2026-07-02: "komplett verzogen, nicht planar; Metall muss in die pi-Ebene").

Local sp2-flattening (_flatten_sp2_atoms_xyz) makes each sp2 atom flat on its OWN
neighbour plane but does NOT force the WHOLE conjugated macrocycle into ONE plane.  This
does: detect the conjugated core (sp2 C/N of the coordination macrocycle, near the metal),
fit its best plane, and gently pull the core + the metal toward it while 1-2 (bond) and
1-3 (angle) distance springs referenced to the CURRENT geometry hold the in-plane shape
so no bond collapses.  Deterministic pure-numpy (reproduces the OB-UFF planarity effect
license-clean, no OpenBabel/CCDC).

NEVER-WORSE by construction when wired ADDITIVELY (the caller adds the planarised frame as
an EXTRA manifold entry, never replacing the original): a genuinely-ruffled macrocycle
(CAZQED beta-octabromo: crystal macroRMS 1.55) keeps its ruffled original as the best
isomer; a system that should be planar (DECTEN/GISHOF/BEGNOT) gets the flat variant that
wins.  best-isomer RMSD (min over frames) can only improve.
"""
import os
import numpy as np

_COV = {"H":0.31,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"Na":1.66,"Mg":1.41,
        "Al":1.21,"Si":1.11,"P":1.07,"S":1.05,"Cl":1.02,"Se":1.20,"Br":1.20,"I":1.39}
_TM = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg La".split())


def _enabled():
    return os.environ.get("DELFIN_FFFREE_PI_PLANARIZE", "0") == "1"


def _f(name, d):
    try:
        return float(os.environ.get(name, ""))
    except Exception:
        return d


def _bonds(syms, P):
    n = len(syms)
    out = []
    for i in range(n):
        ri = _COV.get(syms[i], 0.9)
        for j in range(i + 1, n):
            if float(np.linalg.norm(P[i] - P[j])) < 1.15 * (ri + _COV.get(syms[j], 0.9)):
                out.append((i, j))
    return out


def _angle_pairs(bond_pairs, n):
    adj = [[] for _ in range(n)]
    for i, j in bond_pairs:
        adj[i].append(j); adj[j].append(i)
    out = set()
    for j in range(n):
        nb = adj[j]
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                out.add((min(nb[a], nb[b]), max(nb[a], nb[b])))
    return list(out)


def planarize_frame(syms, P, metal_idx=0, core_radius=None, bond_pairs=None):
    """Return a planarised copy of P, or P unchanged.  Pure numpy, deterministic."""
    if not _enabled():
        return P
    try:
        P = np.asarray(P, float)
        n = len(syms)
        if n < 6 or P.shape[0] != n:
            return P
        if syms[metal_idx] not in _TM:
            metal_idx = next((i for i, s in enumerate(syms) if s in _TM), None)
            if metal_idx is None:
                return P
        m0 = P[metal_idx]
        rad = core_radius if core_radius is not None else _f("DELFIN_PIPL_RADIUS", 3.7)
        # conjugated core = sp2-capable C/N of the coordination macrocycle near the metal
        core = [i for i in range(n) if syms[i] in ("C", "N") and i != metal_idx
                and float(np.linalg.norm(P[i] - m0)) < rad]
        if len(core) < 6:
            return P                                  # not a conjugated coordination core
        bp = bond_pairs if bond_pairs is not None else _bonds(syms, P)
        bp = [(i, j) for (i, j) in bp if 0 <= i < n and 0 <= j < n]
        if not bp:
            return P
        tp = _angle_pairs(bp, n)
        L12 = [float(np.linalg.norm(P[i] - P[j])) for (i, j) in bp]
        L13 = [float(np.linalg.norm(P[i] - P[j])) for (i, j) in tp]
        grp = core + [metal_idx]
        kpl = _f("DELFIN_PIPL_KPL", 0.5)
        k12 = _f("DELFIN_PIPL_K12", 1.0)
        k13 = _f("DELFIN_PIPL_K13", 0.7)
        cap = _f("DELFIN_PIPL_CAP", 0.12)
        iters = int(_f("DELFIN_PIPL_ITERS", 250))
        X = P.copy()
        for _ in range(iters):
            cc = X[core]; c = cc.mean(0)
            _, _, vv = np.linalg.svd(cc - c)
            nrm = vv[2]
            F = np.zeros_like(X)
            for i in grp:                              # pull core + metal toward the plane
                F[i] -= kpl * float(np.dot(X[i] - c, nrm)) * nrm
            for (i, j), L0 in zip(bp, L12):            # bond springs (hold shape, no collapse)
                d = X[i] - X[j]; r = float(np.linalg.norm(d))
                if r > 1e-9:
                    f = -k12 * (r - L0) * d / r; F[i] += f; F[j] -= f
            for (i, j), L0 in zip(tp, L13):            # 1-3 springs (hold angles)
                d = X[i] - X[j]; r = float(np.linalg.norm(d))
                if r > 1e-9:
                    f = -k13 * (r - L0) * d / r; F[i] += f; F[j] -= f
            step = 0.10 * F
            mag = np.linalg.norm(step, axis=1)
            over = mag > cap
            if np.any(over):
                step[over] *= (cap / mag[over])[:, None]
            X = X + step
        # safety: never emit a collapsed bond
        for (i, j) in bp:
            lim = 0.60 * (_COV.get(syms[i], 0.9) + _COV.get(syms[j], 0.9))
            if float(np.linalg.norm(X[i] - X[j])) < lim:
                return P
        return X
    except Exception:
        return np.asarray(P, float)


def macro_rms(syms, P, metal_idx=0, core_radius=3.7):
    """Diagnostic metric: RMS out-of-plane of the conjugated coordination core [A]."""
    try:
        P = np.asarray(P, float)
        mi = metal_idx if syms[metal_idx] in _TM else next(
            (i for i, s in enumerate(syms) if s in _TM), None)
        if mi is None:
            return None
        m = P[mi]
        core = [i for i in range(len(syms)) if syms[i] in ("C", "N")
                and float(np.linalg.norm(P[i] - m)) < core_radius]
        if len(core) < 6:
            return None
        cc = P[core]; c = cc.mean(0)
        _, _, vv = np.linalg.svd(cc - c)
        return float(np.sqrt(np.mean((np.dot(cc - c, vv[2])) ** 2)))
    except Exception:
        return None
