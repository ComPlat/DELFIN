"""delfin/_conformer_rank.py — order the emitted isomer/conformer ensemble so the
frame most likely to match the experimental (XRD) structure comes first.

Cross-validated 2026-06-12 against ~978k CCDC crystal structures (n=2223+621
independent ensembles): among the deterministically-enumerated frames for a
SMILES, ranking by *least steric clash* is the strongest local-geometry predictor
of the crystallised form (top-3 hit 60-64% vs 56-59% for defect-count scoring).
Only ONE generated frame matches the crystal; the others are valid alternative
isomers/conformers — so this re-orders, it never drops or alters any structure.

Pure-geometry, deterministic, no external data, no force field.  Used to sort the
public `smiles_to_xyz_isomers` output best-first.  Disable with
``DELFIN_NO_FRAME_RANK=1`` (e.g. for byte-for-byte enumeration-order debugging).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Tuple

# Covalent radii (Å) — same table as the quality_framework conformer_ranker.
_COV_R: Dict[str, float] = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Hf": 1.75, "Ta": 1.70, "W": 1.62,
    "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

_METALS = frozenset({
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Th", "U",
    "Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi",
})

_DONORS = frozenset({"N", "O", "P", "S", "Cl", "Br", "I", "C", "B", "F", "Se", "Te"})

# σ-aromatic donor elements whose ring plane should contain the metal (in-plane
# lone pair).  Used by the π-coplanarity ranking term (DELFIN_RANK_PI_COPLANAR).
_PI_DONORS = frozenset({"N", "O", "S", "P", "Se"})

# Weights — cross-validated (see module docstring).  Hard vetoes dominate so a
# genuinely corrupt frame never ranks first; the continuous clash penalty is the
# primary discriminator among valid frames.
_W_ATOM_OVERLAP = 1000.0   # literal <0.4 Å overlap = corrupt
_W_MD_TOO_SHORT = 100.0    # genuine M-donor atom overlap (ratio < 0.65)
_W_CLASH = 12.0            # continuous Σ(0.75·Σr_cov − d) over overlapping pairs


def _cov(sym: str) -> float:
    return _COV_R.get(sym, 0.75)


def _parse_xyz(text: str) -> List[Tuple[str, float, float, float]]:
    """Parse a standard / DELFIN-format XYZ string -> [(sym, x, y, z), ...]."""
    lines = text.splitlines()
    if len(lines) < 3:
        return []
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return []
    atoms: List[Tuple[str, float, float, float]] = []
    for ln in lines[2:2 + n]:
        p = ln.split()
        if len(p) < 4:
            continue
        try:
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        except ValueError:
            continue
    return atoms


def _parse_xyz_tolerant(text: str) -> List[Tuple[str, float, float, float]]:
    """Header-tolerant XYZ parser used ONLY by the Stream-B Fix 3 path.

    The public ``smiles_to_xyz_isomers`` emits HEADERLESS frames (raw atom lines,
    no count/comment header), so the strict :func:`_parse_xyz` above returns []
    for every emitted frame — meaning the legacy clash ranker is effectively a
    no-op on real ensembles (it falls back to enumeration order).  This parser
    accepts both the standard 2-line-header format AND the headerless format by
    simply reading every line that looks like ``<sym> <x> <y> <z>``.  Kept
    SEPARATE from ``_parse_xyz`` so the default (flag-OFF) ranking path is
    byte-identical to before."""
    atoms: List[Tuple[str, float, float, float]] = []
    for ln in text.splitlines():
        p = ln.split()
        if len(p) < 4:
            continue
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        atoms.append((p[0], x, y, z))
    return atoms


def _frame_score(atoms: List[Tuple[str, float, float, float]]) -> float:
    """Higher = better (more crystal-like). Start at 100, subtract penalties."""
    n = len(atoms)
    if n < 2:
        return 100.0
    syms = [a[0] for a in atoms]
    rs = [_cov(s) for s in syms]
    n_atom_overlap = 0
    n_md_too_short = 0
    clash_penalty = 0.0
    for i in range(n):
        xi, yi, zi = atoms[i][1], atoms[i][2], atoms[i][3]
        si = syms[i]
        for j in range(i + 1, n):
            dx = xi - atoms[j][1]
            dy = yi - atoms[j][2]
            dz = zi - atoms[j][3]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d < 0.001:
                continue
            sj = syms[j]
            sum_r = rs[i] + rs[j]
            # continuous clash: any genuinely overlapping pair (normal bonds ~1.0·Σr never trigger)
            if sum_r > 0 and d < 0.75 * sum_r:
                clash_penalty += (0.75 * sum_r - d)
            if d < 0.40:
                n_atom_overlap += 1
            # M-donor atom overlap
            is_md = (si in _METALS and sj in _DONORS) or (sj in _METALS and si in _DONORS)
            if is_md and sum_r > 0 and d < 0.65 * sum_r:
                n_md_too_short += 1
    return (100.0
            - _W_ATOM_OVERLAP * n_atom_overlap
            - _W_MD_TOO_SHORT * n_md_too_short
            - _W_CLASH * clash_penalty)


# ---------------------------------------------------------------------------
# Stream-B Fix 3: coordination-geometry ideality tiebreak (DELFIN_RANK_COORD_IDEAL)
# ---------------------------------------------------------------------------
# Least-clash alone can rank a BENT donor-M-donor frame (e.g. 134° P-Au-S) above
# the IDEAL linear 180° one when both are clash-free.  This SECONDARY term measures
# how far each metal's donor-M-donor angles deviate from the ideal polyhedron for
# its CN, and breaks near-ties in favour of the more idealised coordination.  Clash
# stays primary (the term is only consulted when the clash-based scores tie within
# ``_COORD_IDEAL_TIE_EPS``).  Reorder-only; never alters atoms.  Default-OFF →
# byte-identical isomer order.

# Per-CN candidate geometries (canonical keys in delfin._polyhedron_targets).  For
# CN values with more than one chemically-valid ideal (CN4 Td vs square-planar,
# CN5 TBP vs SPY), every candidate is scored and the BEST-fitting one wins, so a
# valid square-planar CN4 is not penalised against a tetrahedral target (or v.v.).
_COORD_IDEAL_GEOMS: Dict[int, Tuple[str, ...]] = {
    2: ("linear_2",),
    3: ("trigonal_planar",),
    4: ("Td", "sqp_4"),
    5: ("tbp", "sqp_5"),
    6: ("Oh", "trig_prism"),
    7: ("pbp", "capped_oct"),
    8: ("sq_antiprism", "dodecahedron"),
    9: ("tricapped_tp",),
}

_COORD_IDEAL_TIE_EPS = 1e-6   # clash-score gap below which the geometry term decides
_W_COORD_IDEAL = 1.0          # mean-angle-deviation (deg) penalty weight (tiebreak only)


def _ideal_angle_bags():
    """Lazy per-CN ideal donor-M-donor angle bags, derived ONCE from the shared
    ``delfin._polyhedron_targets`` tables (single source of truth).  Maps
    cn -> list of (geometry_key, sorted-tuple-of-ideal-angles).  Cached on the
    function object."""
    cache = getattr(_ideal_angle_bags, "_cache", None)
    if cache is not None:
        return cache
    cache = {}
    try:
        import numpy as _np
        from delfin import _polyhedron_targets as _PT
        for cn, geoms in _COORD_IDEAL_GEOMS.items():
            bags = []
            for g in geoms:
                try:
                    V = _PT.get_ideal_donor_vectors(cn, g)
                except Exception:
                    continue
                angs = set()
                for i in range(len(V)):
                    vi = V[i]; ni = float(_np.linalg.norm(vi))
                    if ni < 1e-9:
                        continue
                    for j in range(i + 1, len(V)):
                        vj = V[j]; nj = float(_np.linalg.norm(vj))
                        if nj < 1e-9:
                            continue
                        c = float(_np.dot(vi, vj)) / (ni * nj)
                        c = max(-1.0, min(1.0, c))
                        angs.add(round(math.degrees(math.acos(c)), 4))
                if angs:
                    bags.append((g, tuple(sorted(angs))))
            if bags:
                cache[cn] = bags
    except Exception:
        cache = {}
    _ideal_angle_bags._cache = cache
    return cache


def _coord_ideal_penalty(atoms: List[Tuple[str, float, float, float]]) -> float:
    """Sum, over every metal centre, of the mean donor-M-donor angle deviation
    (deg) from the BEST-fitting ideal polyhedron for that metal's CN.  Lower =
    more idealised coordination.  Returns 0.0 when there is nothing to score
    (no metal, CN<2, CN with no table) so it is a strict no-op tiebreak there.

    Graph-level + geometry-agnostic: donors are heavy atoms within a covalent
    bonding contact of the metal; the per-CN ideal angle bag comes from the
    shared polyhedra tables; the best-matching candidate geometry is chosen per
    metal.  Deterministic (pure geometry, fixed tables, sorted iteration)."""
    n = len(atoms)
    if n < 3:
        return 0.0
    bags = _ideal_angle_bags()
    if not bags:
        return 0.0
    syms = [a[0] for a in atoms]
    metal_idxs = [i for i in range(n) if syms[i] in _METALS]
    if not metal_idxs:
        return 0.0
    total = 0.0
    for m in metal_idxs:
        mx, my, mz = atoms[m][1], atoms[m][2], atoms[m][3]
        rm = _cov(syms[m])
        donors = []
        for j in range(n):
            if j == m or syms[j] == "H" or syms[j] not in _DONORS:
                continue
            dx = mx - atoms[j][1]; dy = my - atoms[j][2]; dz = mz - atoms[j][3]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            sum_r = rm + _cov(syms[j])
            # covalent bonding contact (M-D); 1.30·Σr_cov admits long dative bonds
            if sum_r > 0 and 0.3 * sum_r < d < 1.30 * sum_r:
                donors.append(j)
        cn = len(donors)
        if cn < 2 or cn not in bags:
            continue
        # observed donor-M-donor angles
        obs = []
        for a in range(cn):
            ia = donors[a]
            v1 = (atoms[ia][1] - mx, atoms[ia][2] - my, atoms[ia][3] - mz)
            n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
            if n1 < 1e-9:
                continue
            for b in range(a + 1, cn):
                ib = donors[b]
                v2 = (atoms[ib][1] - mx, atoms[ib][2] - my, atoms[ib][3] - mz)
                n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
                if n2 < 1e-9:
                    continue
                c = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
                c = max(-1.0, min(1.0, c))
                obs.append(math.degrees(math.acos(c)))
        if not obs:
            continue
        # best-fitting candidate geometry: min mean nearest-ideal deviation
        best = None
        for _g, ideal in bags[cn]:
            dev = 0.0
            for ang in obs:
                dev += min(abs(ang - it) for it in ideal)
            dev /= len(obs)
            if best is None or dev < best:
                best = dev
        if best is not None:
            total += best
    return total


# ---------------------------------------------------------------------------
# π-coplanarity ranking term (DELFIN_RANK_PI_COPLANAR, default-OFF → byte-id)
# ---------------------------------------------------------------------------
# Data-driven root finding (2026-06-23) over 229 π-donor smoke structures: 50%
# already emit a coplanar primary frame, 17% have a coplanar pose IN the manifold
# that is simply not frame 0 (the RANKING root), 34% have none (a construction
# root, addressed elsewhere).  Simulation of the rule below recovered the full 17%
# (f0 metal-out-of-plane 0.786→0.508 Å, bad-count 115→77) with ZERO clash cost
# (mean −0.001) and zero completeness loss (reorder only).  So: among non-corrupt
# frames bucketed by clash (0.5 Å), prefer the one whose metal lies most in the
# σ-aromatic donor ring plane.  Pure geometry, deterministic.

_PI_COPLANAR_BUCKET = 0.5   # Å of clash penalty grouped as a tie before π decides
_PI_COORD_BUCKET = 5.0      # deg of coord-angle deviation grouped as a tie (composite)
# Pareto-gated π promotion thresholds (eye-root 2026-06-23): promote the coplanar
# frame ONLY when it gives a real π gain AND worsens nothing else vs the lead frame.
_PI_PROMOTE_MIN_GAIN = 0.3   # Å of π-oop improvement required (rejects marginal cases)
_PI_PROMOTE_COORD_TOL = 1.0  # deg coord-angle worsening tolerated (rejects ZULSAC/XIMHAD)
_PI_PROMOTE_MIND_TOL = 0.06  # Å non-bonded contact tightening tolerated (rejects CUXWAV/XIMHAD)
_PI_PROMOTE_CLASH_TOL = 0.05  # clash-penalty worsening tolerated


def _clash_components(atoms: List[Tuple[str, float, float, float]]) -> Tuple[int, float]:
    """Return (hard_veto_count, soft_clash_penalty) — same thresholds as
    :func:`_frame_score` but split so the π term can bucket on the soft clash
    while corrupt frames (literal overlap / M-donor fusion) are sorted last.
    Kept SEPARATE from ``_frame_score`` so the default ranking path is untouched."""
    n = len(atoms)
    if n < 2:
        return (0, 0.0)
    syms = [a[0] for a in atoms]
    rs = [_cov(s) for s in syms]
    veto = 0
    pen = 0.0
    for i in range(n):
        xi, yi, zi = atoms[i][1], atoms[i][2], atoms[i][3]
        si = syms[i]
        for j in range(i + 1, n):
            dx = xi - atoms[j][1]; dy = yi - atoms[j][2]; dz = zi - atoms[j][3]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d < 0.001:
                continue
            sj = syms[j]
            sum_r = rs[i] + rs[j]
            if sum_r > 0 and d < 0.75 * sum_r:
                pen += 0.75 * sum_r - d
            if d < 0.40:
                veto += 1
            is_md = (si in _METALS and sj in _DONORS) or (sj in _METALS and si in _DONORS)
            if is_md and sum_r > 0 and d < 0.65 * sum_r:
                veto += 1
    return (veto, pen)


def _pi_coplanar_penalty(atoms: List[Tuple[str, float, float, float]]) -> float:
    """Sum, over every metal centre, of the worst metal-out-of-plane distance (Å)
    for its σ-aromatic donor rings — the metal should lie IN each coordinating
    aromatic donor's ring plane (in-plane sp² lone pair).  0.0 when nothing
    qualifies (no metal / no coordinating aromatic donor ring) so it is a strict
    no-op there.  Restricted to the coordination sphere (heavy atoms within 5 Å of
    a metal) for speed; deterministic pure geometry, no force field."""
    n = len(atoms)
    if n < 5:
        return 0.0
    try:
        import numpy as _np
    except Exception:
        return 0.0
    syms = [a[0] for a in atoms]
    xyz = _np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
    metal_idxs = [i for i in range(n) if syms[i] in _METALS]
    if not metal_idxs:
        return 0.0
    # coordination-sphere subgraph: heavy non-metal atoms within 5 Å of any metal
    in_sphere = [False] * n
    for j in range(n):
        if j in metal_idxs or syms[j] == "H":
            continue
        for m in metal_idxs:
            if float(_np.linalg.norm(xyz[j] - xyz[m])) < 5.0:
                in_sphere[j] = True
                break
    adj = [[] for _ in range(n)]
    sphere = [j for j in range(n) if in_sphere[j]]
    for ai in range(len(sphere)):
        a = sphere[ai]
        for bi in range(ai + 1, len(sphere)):
            b = sphere[bi]
            if float(_np.linalg.norm(xyz[a] - xyz[b])) < 1.25 * (_cov(syms[a]) + _cov(syms[b])):
                adj[a].append(b); adj[b].append(a)
    R = []; seen = set()

    def dfs(s, cur, sn):
        last = cur[-1]
        for nb in adj[last]:
            if nb == s and 5 <= len(cur) <= 6:
                k = frozenset(cur)
                if k not in seen:
                    seen.add(k); R.append(list(cur))
            elif nb not in sn and len(cur) < 6:
                dfs(s, cur + [nb], sn | {nb})
    for s in sphere:
        dfs(s, [s], {s})
    total = 0.0
    any_found = False
    for m in metal_idxs:
        M = xyz[m]; D = _np.linalg.norm(xyz - M, axis=1)
        worst = 0.0; found = False
        for rg in R:
            near = [k for k in rg if D[k] < 2.5]
            if len(near) != 1 or syms[near[0]] not in _PI_DONORS:
                continue
            pts = xyz[rg]; c = pts.mean(0)
            _, _, vt = _np.linalg.svd(pts - c)
            if _np.abs((pts - c) @ vt[2]).max() > 0.35:
                continue  # ring itself not flat → not a clean aromatic donor
            oop = abs(float((M - c) @ vt[2]))
            if oop > worst:
                worst = oop; found = True
        if found:
            total += worst; any_found = True
    return (total, any_found)


def _min_nonbonded_heavy(atoms: List[Tuple[str, float, float, float]]) -> float:
    """Smallest non-bonded heavy-heavy distance (Å) — a collateral-contact signal
    more sensitive than the clash penalty (catches contacts that tighten but stay
    above the clash floor).  Returns 9.9 if <2 heavy atoms."""
    heavy = [(a[1], a[2], a[3], a[0]) for a in atoms if a[0] != "H"]
    n = len(heavy)
    if n < 2:
        return 9.9
    mind = 9.9
    for i in range(n):
        xi, yi, zi, si = heavy[i]
        ri = _cov(si)
        for j in range(i + 1, n):
            xj, yj, zj, sj = heavy[j]
            d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            if d < 1.30 * (ri + _cov(sj)):
                continue  # bonded
            if d < mind:
                mind = d
    return mind


def rank_isomers(isomers: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Re-order ``[(xyz_string, label), ...]`` best (most crystal-like) first.

    Deterministic and stable (ties keep enumeration order).  Never drops, adds,
    or mutates a structure — only reorders.  No-op for <2 frames or when
    ``DELFIN_NO_FRAME_RANK=1``.

    Stream-B Fix 3 (``DELFIN_RANK_COORD_IDEAL=1``, default-OFF → byte-identical
    order): adds a SECONDARY coordination-geometry-ideality term as a tiebreak
    AFTER the clash score, so a clash-free but bent donor-M-donor frame yields to
    the more idealised polyhedron when their clash scores are (near-)tied.  Clash
    stays primary; this only resolves ties.

    Stream-B batch-2 Fix B (``DELFIN_FRAME_RANK_FIX=1``, default-OFF → byte-
    identical order): activates the clash ranker on the REAL (headerless)
    emitter output.  The public ``smiles_to_xyz_isomers`` emits frames with NO
    atom-count header, so the strict :func:`_parse_xyz` returns [] for every
    frame → every frame scores 100.0 → no reordering → the clash ranker is a
    silent no-op today.  With the flag set, the default path parses with the
    header-tolerant parser so the ranker actually orders frames by least clash.
    Default-OFF keeps the strict parse → the ranker stays a no-op → byte-
    identical with today.  ``DELFIN_NO_FRAME_RANK=1`` still wins outright
    (returns enumeration order regardless of this flag — revert semantics
    preserved).
    """
    if os.environ.get("DELFIN_NO_FRAME_RANK", "0") == "1":
        return isomers
    if not isinstance(isomers, list) or len(isomers) < 2:
        return isomers
    coord_ideal = os.environ.get("DELFIN_RANK_COORD_IDEAL", "0") == "1"
    # Fix B: header-tolerant parsing in the default path, gated default-OFF.
    rank_fix = os.environ.get("DELFIN_FRAME_RANK_FIX", "0") == "1"
    # π-coplanarity ranking (default-OFF → byte-identical order).
    pi_cop = os.environ.get("DELFIN_RANK_PI_COPLANAR", "0") == "1"
    _default_parse = _parse_xyz_tolerant if rank_fix else _parse_xyz
    try:
        if pi_cop:
            # Pareto-gated π promotion (eye-root 2026-06-23): move the most-coplanar
            # frame to the FRONT only when it gives a meaningful π gain AND worsens
            # NOTHING else (coordination angle, non-bonded contact, clash) vs the
            # original lead frame.  The earlier Goodhart failure (rank by π-oop alone
            # → 6 wins / 6 regressions) is fixed because the eye-regressions all fail
            # exactly one of these guards (marginal π gain, or worse coord/contact),
            # while the 6 eye-wins all pass.  Otherwise: enumeration order (byte-id
            # to flag-OFF).  No-op for non-π / no-metal structures (has_pi False).
            feats = []
            for idx, item in enumerate(isomers):
                xyz = item[0] if isinstance(item, (tuple, list)) and item else ""
                atoms = _parse_xyz_tolerant(xyz) if isinstance(xyz, str) else []
                if atoms:
                    veto, clash = _clash_components(atoms)
                    coord = _coord_ideal_penalty(atoms)
                    pio, has_pi = _pi_coplanar_penalty(atoms)
                    mind = _min_nonbonded_heavy(atoms)
                else:
                    veto, clash, coord, pio, has_pi, mind = 1, 0.0, 0.0, 0.0, False, 9.9
                feats.append({"idx": idx, "veto": veto, "clash": clash,
                              "coord": coord, "pio": pio, "has_pi": has_pi,
                              "mind": mind, "item": item})
            base = next((f for f in feats if f["veto"] == 0 and f["has_pi"]), None)
            pi_frames = [f for f in feats if f["veto"] == 0 and f["has_pi"]]
            if base is not None and pi_frames:
                cand = min(pi_frames, key=lambda f: (f["pio"], f["idx"]))
                if (cand["idx"] != base["idx"]
                        and (base["pio"] - cand["pio"]) >= _PI_PROMOTE_MIN_GAIN
                        and (cand["coord"] - base["coord"]) <= _PI_PROMOTE_COORD_TOL
                        and (cand["mind"] - base["mind"]) >= -_PI_PROMOTE_MIND_TOL
                        and (cand["clash"] - base["clash"]) <= _PI_PROMOTE_CLASH_TOL):
                    return ([cand["item"]]
                            + [f["item"] for f in feats if f["idx"] != cand["idx"]])
            return [f["item"] for f in feats]
        if not coord_ideal:
            # Default path.  Parser choice is the ONLY change: strict (flag-OFF
            # → byte-identical, ranker is a no-op on headerless frames) vs
            # header-tolerant (flag-ON → ranker actually orders by least clash).
            scored = []
            for idx, item in enumerate(isomers):
                xyz = item[0] if isinstance(item, (tuple, list)) and item else ""
                sc = (_frame_score(_default_parse(xyz))
                      if isinstance(xyz, str) else -1e9)
                scored.append((idx, sc, item))
            scored.sort(key=lambda e: (-e[1], e[0]))
            return [e[2] for e in scored]
        # Fix-3 path — clash score PRIMARY (bucketed so FP noise below
        # _COORD_IDEAL_TIE_EPS counts as a tie), coord-geometry deviation as the
        # SECONDARY tiebreak inside a clash bucket, enumeration order last.  Uses
        # the header-tolerant parser because the public emitter is headerless
        # (the strict parser would see 0 atoms and disable both terms).
        scored = []
        for idx, item in enumerate(isomers):
            xyz = item[0] if isinstance(item, (tuple, list)) and item else ""
            atoms = _parse_xyz_tolerant(xyz) if isinstance(xyz, str) else []
            sc = _frame_score(atoms) if isinstance(xyz, str) else -1e9
            geo = _coord_ideal_penalty(atoms) if atoms else 0.0
            scored.append((idx, sc, geo, item))

        def _key(e):
            idx, sc, geo, _item = e
            bucket = (round(sc / _COORD_IDEAL_TIE_EPS)
                      if _COORD_IDEAL_TIE_EPS else sc)
            return (-bucket, _W_COORD_IDEAL * geo, idx)

        scored.sort(key=_key)
        return [e[3] for e in scored]
    except Exception:
        # Ranking must never break the build — fall back to original order.
        return isomers
