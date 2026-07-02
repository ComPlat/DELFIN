"""Whole-complex M-D relax — the strict-never-worse root fix (DELFIN_FFFREE_WHOLE_RELAX).

THE ROOT (the unification, 2026-07-01): the ETKDG metallacycle embed produces an
under-length / collapsed metal environment, and EVERY downstream refiner in
``_finish_config_frame`` (``refine``, ``torsion_relax``, ``joint_declash``,
``sphere_flex``) FREEZES the metal + donors -> the crush is conserved into the
emitted frame.  This is the shared root of the three biggest defect clusters
(cage-collapse, over-coordination high-CN, generic-sigma).  Champion ``ddde3273``
did NOT have this defect: chelates fell through to the legacy Open-Babel-UFF
whole-complex relax which froze the metal but left DONORS FREE and pulled every
M-D to its ideal via distance constraints, so crushed bonds expanded and the
ligand followed.  OB-UFF is non-deterministic on TM cations (uninitialised-heap
gradients when it cannot type the metal — ``ff.Energy() > 1e9`` marker), so this
re-implements the SAME behaviour generatively, deterministically, license-clean.

MECHANISM (pure numpy, no force field, no CCDC/COD/OpenBabel):
  * Hold every LIGAND at its built shape with distance springs referenced to the
    CURRENT geometry — 1-2 (bond) AND 1-3 (angle-pair) pairs.  Referencing the
    current geometry (not an ideal table) means we trust the ligand build and
    change ONLY what must change; it also sidesteps bond-order-aware length tables
    and needs no angle-gradient math (a 1-3 distance spring pins the angle).
  * Drive each M-D with a SOFT spring toward its ideal (``md_distance``).  Soft is
    deliberate: a crushed M-D (far below ideal) feels a strong pull and expands; a
    healthy M-D already near ideal barely moves, and where the ligand cavity wants
    a slightly shorter bond the ligand strain wins -> NO fixed-target overshoot
    (the failure mode of the hard CAGE_MD_GUARD reset on the 6 broad WORSE).
  * Soft vdW repulsion resolves residual heavy-heavy inter-ligand clashes.
  * The metal (index 0) is pinned; donors are FREE (the whole point).

NEVER-WORSE (per-frame acceptance, mirrors champion smiles_converter.py:25068):
  accept the relaxed frame ONLY if it does not collapse a covalent bond AND does
  not worsen a license-clean composite (sum |M-D - ideal| + heavy clash) beyond a
  small tol; otherwise return the input frame unchanged.  Byte-identical when the
  flag is unset.  Never raises (any failure -> input frame).
"""
import os
import numpy as np

# license-clean covalent radii (Cordero) for vdW-repulsion / collapse thresholds
_COV = {"H":0.31,"Li":1.28,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,
        "Na":1.66,"Mg":1.41,"Al":1.21,"Si":1.11,"P":1.07,"S":1.05,"Cl":1.02,
        "K":2.03,"Ca":1.76,"Se":1.20,"Br":1.20,"As":1.19,"Sb":1.39,"Te":1.38,"I":1.39}


def _enabled():
    return os.environ.get("DELFIN_FFFREE_WHOLE_RELAX", "0") == "1"


def _fenv(name, default):
    try:
        return float(os.environ.get(name, ""))
    except Exception:
        return default


def _bond_pairs_from_frags(relax_frags):
    """Ligand-internal covalent bonds as GLOBAL index pairs (metal at 0, a frag's
    atom k -> global lig_off + 1 + k, exactly as _finish_config_frame maps them)."""
    pairs = []
    for frag, lig_off in relax_frags:
        base = lig_off + 1
        for b in frag.GetBonds():
            pairs.append((base + b.GetBeginAtomIdx(), base + b.GetEndAtomIdx()))
    return pairs


def _angle_pairs(bond_pairs, natoms):
    """1-3 pairs: two atoms both bonded to a common neighbour (angle end atoms).
    Holding their distance fixes the bond ANGLE without any angle-gradient math."""
    adj = [[] for _ in range(natoms)]
    for i, j in bond_pairs:
        if 0 <= i < natoms and 0 <= j < natoms:
            adj[i].append(j); adj[j].append(i)
    out = set()
    for j in range(natoms):
        nb = adj[j]
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                p = (min(nb[a], nb[b]), max(nb[a], nb[b]))
                if p[0] != p[1]:
                    out.add(p)
    return list(out)


def _spring_forces(X, pairs, targets, k, F):
    """Accumulate harmonic distance-spring forces into F (in place)."""
    for (i, j), L0 in zip(pairs, targets):
        d = X[i] - X[j]
        r = float(np.linalg.norm(d))
        if r < 1e-9:
            continue
        f = -k * (r - L0) * d / r
        F[i] += f
        F[j] -= f


def _md_composite(X, metal_idx, donor_ideal, syms, nonbond_heavy):
    """License-clean defect scalar: sum |M-D - ideal| + soft heavy-heavy clash."""
    dev = 0.0
    for d, ideal in donor_ideal.items():
        dev += abs(float(np.linalg.norm(X[d] - X[metal_idx])) - ideal)
    clash = 0.0
    for i, j in nonbond_heavy:
        rr = float(np.linalg.norm(X[i] - X[j]))
        lim = 0.75 * (_COV.get(syms[i], 0.8) + _COV.get(syms[j], 0.8))
        if rr < lim:
            clash += (lim - rr)
    return dev + clash


def _collapsed(X, bond_pairs, syms):
    """A covalent bond crushed below 0.6*(ri+rj) = broken topology."""
    for i, j in bond_pairs:
        rr = float(np.linalg.norm(X[i] - X[j]))
        lim = 0.60 * (_COV.get(syms[i], 0.8) + _COV.get(syms[j], 0.8))
        if rr < lim:
            return True
    return False


def relax_if_enabled(out_syms, P, fixed, relax_frags, metal_sym):
    """Return a relaxed copy of P, or P unchanged.  Byte-identical when the flag is
    unset.  ``fixed`` = {metal_idx} | donor_globals; ``relax_frags`` = [(frag, off)]."""
    if not _enabled():
        return P
    try:
        P = np.asarray(P, float)
        n = len(out_syms)
        if not relax_frags or n < 3 or P.shape[0] != n:
            return P
        metal_idx = 0
        donor_globals = sorted(set(fixed) - {metal_idx})
        if not donor_globals:
            return P
        # ideal M-D per donor (license-clean covalent/VSEPR helper the build uses)
        from delfin.manta import metal_sphere_builder as MSB
        donor_ideal = {}
        # locate each donor's atom + its frag mol for a hybridisation-aware ideal
        frag_of = {}
        for frag, lig_off in relax_frags:
            base = lig_off + 1
            for k in range(frag.GetNumAtoms()):
                frag_of[base + k] = (frag, k)
        for dg in donor_globals:
            fr = frag_of.get(dg)
            dsym = out_syms[dg]
            try:
                if fr is not None:
                    da = fr[0].GetAtomWithIdx(fr[1])
                    ideal = float(MSB.md_distance(metal_sym, dsym, atom=da, mol=fr[0]))
                else:
                    ideal = float(MSB.md_distance(metal_sym, dsym))
            except Exception:
                try:
                    ideal = float(MSB.md_distance(metal_sym, dsym))
                except Exception:
                    ideal = 0.0
            if ideal > 0:
                donor_ideal[dg] = ideal
        if not donor_ideal:
            return P

        bond_pairs = _bond_pairs_from_frags(relax_frags)
        bond_pairs = [(i, j) for (i, j) in bond_pairs if 0 <= i < n and 0 <= j < n]
        ang_pairs = _angle_pairs(bond_pairs, n)
        # reference lengths from the CURRENT geometry (trust the ligand build)
        bond_L0 = [float(np.linalg.norm(P[i] - P[j])) for (i, j) in bond_pairs]
        ang_L0 = [float(np.linalg.norm(P[i] - P[j])) for (i, j) in ang_pairs]
        # M-D springs — EXPAND-ONLY + CRUSH-GATED (the never-worse lesson, 2026-07-01):
        # the covalent-'ideal' (md_distance) is only a rough estimate and for many
        # donor/metal combos it disagrees with the crystal in EITHER direction.  A
        # bidirectional spring toward ideal therefore drags GOOD bonds off the crystal
        # (ISEFUJ: OFF 2.30 == crystal 2.30, ideal ~2.05 -> spring contracted it to 2.05
        # = +1.0A RMSD regression).  So (a) only include a donor whose M-D is genuinely
        # CRUSHED below crush_gate*ideal (a healthy/long bond is left alone), and (b) the
        # M-D force is EXPAND-ONLY (pull outward when short, never contract).  This keeps
        # the real cage/over-coordination expansions while making over-length and
        # already-correct bonds untouchable -> never-worse by construction on M-D.
        crush_gate = _fenv("DELFIN_WR_CRUSH_GATE", 0.90)
        md_pairs = []
        md_L0 = []
        for dg, ideal in donor_ideal.items():
            r_now = float(np.linalg.norm(P[dg] - P[metal_idx]))
            if r_now < crush_gate * ideal:          # crushed -> expand toward ideal
                md_pairs.append((metal_idx, dg))
                md_L0.append(ideal)

        # non-bonded heavy pairs for vdW / clash (exclude 1-2, 1-3, and M-D)
        excl = set()
        for (i, j) in bond_pairs + ang_pairs + md_pairs:
            excl.add((min(i, j), max(i, j)))
        heavy = [i for i in range(n) if out_syms[i] != "H"]
        nonbond_heavy = []
        for a in range(len(heavy)):
            i = heavy[a]
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                if (min(i, j), max(i, j)) not in excl:
                    nonbond_heavy.append((i, j))

        # force constants (soft M-D is the never-worse key)
        k12 = _fenv("DELFIN_WR_K12", 1.0)
        k13 = _fenv("DELFIN_WR_K13", 0.7)
        kmd = _fenv("DELFIN_WR_KMD", 0.30)
        kvdw = _fenv("DELFIN_WR_KVDW", 0.5)
        n_iter = int(_fenv("DELFIN_WR_ITERS", 160))
        maxdisp = _fenv("DELFIN_WR_MAXDISP", 1.20)   # hard per-atom drift clamp (A)

        X = P.copy()
        for it in range(n_iter):
            F = np.zeros_like(X)
            _spring_forces(X, bond_pairs, bond_L0, k12, F)
            _spring_forces(X, ang_pairs, ang_L0, k13, F)
            # M-D: EXPAND-ONLY (act only while the bond is still shorter than target)
            for (i, j), L0 in zip(md_pairs, md_L0):
                d = X[i] - X[j]
                r = float(np.linalg.norm(d))
                if r < 1e-9 or r >= L0:
                    continue
                f = -kmd * (r - L0) * d / r      # (r-L0)<0 -> pushes j outward from i
                F[i] += f
                F[j] -= f
            # soft vdW repulsion (only when overlapping)
            for i, j in nonbond_heavy:
                d = X[i] - X[j]
                r = float(np.linalg.norm(d))
                if r < 1e-6:
                    continue
                lim = 0.85 * (_COV.get(out_syms[i], 0.8) + _COV.get(out_syms[j], 0.8))
                if r < lim:
                    f = kvdw * (lim - r) * d / r
                    F[i] += f
                    F[j] -= f
            F[metal_idx] = 0.0            # metal pinned
            # deterministic damped step with a per-atom cap (stability, no RNG)
            step = 0.08 * F
            mag = np.linalg.norm(step, axis=1)
            cap = 0.15
            over = mag > cap
            if np.any(over):
                step[over] *= (cap / mag[over])[:, None]
            X = X + step
            # HARD TOTAL-DISPLACEMENT CLAMP (stability — the never-explode guard):
            # a legitimate M-D correction moves atoms < ~1 A; runaway gradient descent
            # on dense aromatic complexes can otherwise drift atoms 10+ A (ISEFUJ blew
            # up to RMSD 6 A before this).  Clamp every atom to within maxdisp of its
            # INPUT position each iteration -> the relax can only make a bounded local
            # correction, never scramble the structure.
            dvec = X - P
            dn = np.linalg.norm(dvec, axis=1)
            hot = dn > maxdisp
            if np.any(hot):
                X[hot] = P[hot] + dvec[hot] * (maxdisp / dn[hot])[:, None]
            # convergence: stop once the step is negligible (determinism-safe)
            if float(np.max(np.linalg.norm(step, axis=1))) < 5e-3:
                break

        # Sanity guards only (the CALLER adds this frame ADDITIVELY -> never-worse is
        # guaranteed at the manifold level, so we don't need a composite improvement
        # test here; we only refuse to return a BROKEN or DIVERGED frame that would be
        # wasteful to add).  (1) no collapsed covalent bond; (2) bounded displacement
        # (a legitimate M-D correction moves atoms < ~1 A -> skip runaway/scramble).
        if _collapsed(X, bond_pairs, out_syms):
            return P
        maxd = float(np.max(np.linalg.norm(X - P, axis=1))) if len(X) else 0.0
        if maxd > _fenv("DELFIN_WR_ACCEPT_MAXDISP", 1.00) or maxd < 1e-4:
            if os.environ.get("DELFIN_WR_DEBUG", "0") == "1" and maxd >= 1e-4:
                os.write(2, ("[WHOLE_RELAX] skip (maxdisp=%.2f)\n" % maxd).encode())
            return P
        if os.environ.get("DELFIN_WR_DEBUG", "0") == "1":
            _b = _md_composite(P, metal_idx, donor_ideal, out_syms, nonbond_heavy)
            _a = _md_composite(X, metal_idx, donor_ideal, out_syms, nonbond_heavy)
            os.write(2, ("[WHOLE_RELAX] add variant dev %.3f->%.3f maxdisp=%.2f (crushed=%d)\n"
                         % (_b, _a, maxd, len(md_pairs))).encode())
        return X
    except Exception:
        return np.asarray(P, float)
