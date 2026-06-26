"""Phase 1.5 v2 ensemble per-frame composite scorer + polyhedron fingerprint.

Self-contained — NO runtime dependency on the agent_workspace quality_framework
detector tree.  All detector logic that the ensemble dispatcher relies on has
been copied / re-derived here so the scoring module is importable from a clean
DELFIN install (production-grade requirement: HEAD must not break when the
research workspace is absent).

Detectors implemented (all return a *cost* — lower is better):

    1. extras            — mismatched heavy-bonds vs SMILES graph (mirrors
                           ``smiles_converter._count_extra_heavy_bonds``)
    2. md_invariant      — max |‖M-D‖ - 2.0 Å| deviation across all metal-donor
                           pairs detected via covalent-radius cutoffs
    3. heavy_clash       — count of heavy-heavy pairs at d < 0.85 · Σr_vdw
    4. h_clash           — count of H-H pairs at d < 1.6 Å
    5. tridentate_oop    — count of tridentate donor sets with metal-OOP > 0.30 Å
    6. bond_outliers     — count of heavy-heavy bonds at d > 1.3 · (r_i + r_j)

Composite cost: ``Σ (weight_k · raw_metric_k)`` with default weights selected
in the dispatcher (5 / 3 / 1 / 0.5 / 2 / 0.5).  Lower = better.

Fingerprint:
    ``fingerprint(xyz)`` returns a string-stable canonical-form key (geometry
    polyhedron + sorted donor-element tuple) that groups isomer-equivalent
    frames across paths.  When a polyhedron cannot be classified (no metal,
    low-cn, etc.) we fall back to a heavy-atom-signature hash so that fingerprint
    is *always* defined and the dispatcher's grouping logic is total.

Performance target: <50 ms per frame on a 60-atom complex.  Achieved by
    - pre-computing covalent-radius lookups once per call
    - O(N²) heavy-atom pair scan with squared-distance early-exit
    - lazy-loading scipy only inside polyhedron classifier (a no-metal frame
      never pays the import cost)
"""
from __future__ import annotations

import hashlib
import math
import re
from collections import Counter, defaultdict
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------- #
# Constants  (mirrors find_chelate_planarity / find_isomer_coverage_named)
# ---------------------------------------------------------------------- #

_METALS: FrozenSet[str] = frozenset(
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "Hf Ta W  Re Os Ir Pt Au Hg "
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split()
)

# Cordero 2008 covalent radii (Å); fallback 1.50 for unknowns.
_COV: Dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Se": 1.20, "As": 1.19, "Br": 1.20, "I": 1.39, "Te": 1.38,
    # Common metals (Cordero / Pyykkö averages)
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47,
    "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
}

# van der Waals radii (Bondi) for clash detection.
_VDW: Dict[str, float] = {
    "H": 1.20, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75,
    "Se": 1.90, "As": 1.85, "Br": 1.85, "I": 1.98, "Te": 2.06,
}
_VDW_DEFAULT_HEAVY: float = 2.00
_VDW_DEFAULT_METAL: float = 2.00  # Iter-7 lesson: M-D ~2.0 Å, do NOT clash-flag

_DONOR_CUT: Dict[str, float] = {
    "N": 2.40, "O": 2.30, "P": 2.70, "S": 2.70, "C": 2.55,
    "Cl": 2.80, "Br": 2.95, "I": 3.10, "F": 2.20, "H": 1.85,
}
_BOND_CUT: float = 1.85  # generic non-metal-non-metal cutoff

# Bond-tolerance for the extras detector (mirrors smiles_converter).
_TOL_ORG: float = 0.45
_TOL_METAL: float = 0.60

# Idealised polyhedron vertex sets (mirror find_isomer_coverage_named.py).
_GEOM_VECTORS_RAW: Dict[str, List[Tuple[float, float, float]]] = {
    'LIN': [(2, 0, 0), (-2, 0, 0)],
    'TP':  [(2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'TS':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0)],
    'OH':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)],
    'SQ':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0)],
    'TH':  [(1.155, 1.155, 1.155), (-1.155, -1.155, 1.155),
            (-1.155, 1.155, -1.155), (1.155, -1.155, -1.155)],
    'TBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'SP':  [(0, 0, 2.2), (2, 0, 0.4), (0, 2, 0.4), (-2, 0, 0.4), (0, -2, 0.4)],
    'PBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (0.618, 1.902, 0),
            (-1.618, 1.176, 0), (-1.618, -1.176, 0), (0.618, -1.902, 0)],
    'SAP': [(1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
            (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)],
}

# Default scoring weights — every entry mappable via score_weights kwarg.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "extras":          5.0,
    "md_invariant":    3.0,
    "heavy_clash":     1.0,
    "h_clash":         0.5,
    "tridentate_oop":  2.0,
    "bond_outliers":   0.5,
}

# Canonical M-D bond length (Å) for invariant deviation scoring.
_MD_CANONICAL: float = 2.00


# ---------------------------------------------------------------------- #
# Lightweight XYZ parser (no RDKit needed)
# ---------------------------------------------------------------------- #

def _parse_xyz(xyz_text: str) -> Optional[List[Tuple[str, float, float, float]]]:
    """Parse an XYZ string into ``[(symbol, x, y, z), ...]``.

    Handles three flavours:
        1. DELFIN ``smiles_to_xyz_isomers`` output: pure atom rows (no count
           line, no comment line) — the most common form.
        2. Standard XYZ with leading ``n`` and comment header.
        3. First frame of a multi-frame XYZ.

    Detection: parse line 0 — if it is a single integer in [1, 1e5] AND there
    are at least ``n + 2`` remaining lines that look like atom rows, treat it
    as a header.  Otherwise treat the whole input as atom rows.

    Returns ``None`` on any parse failure.
    """
    if not xyz_text:
        return None
    try:
        lines = [ln for ln in xyz_text.splitlines() if ln.strip()]
        if not lines:
            return None
        # Header sniff
        header_offset = 0
        try:
            n = int(lines[0].strip())
            if 1 <= n <= 100000 and len(lines) >= n + 2:
                # Likely standard XYZ.  Validate by checking that lines[2..n+1]
                # parse as atom rows.
                for ln in lines[2 : 2 + n]:
                    parts = ln.split()
                    if len(parts) < 4:
                        raise ValueError
                    float(parts[1]); float(parts[2]); float(parts[3])
                header_offset = 2
                lines = lines[2 : 2 + n]
        except (ValueError, IndexError):
            pass
        atoms: List[Tuple[str, float, float, float]] = []
        for ln in lines:
            parts = ln.split()
            if len(parts) < 4:
                continue
            sym = parts[0]
            try:
                x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
            except ValueError:
                continue
            atoms.append((sym, x, y, z))
        if not atoms:
            return None
        return atoms
    except Exception:
        return None


def _dist2(a: Tuple[str, float, float, float],
           b: Tuple[str, float, float, float]) -> float:
    dx = a[1] - b[1]; dy = a[2] - b[2]; dz = a[3] - b[3]
    return dx * dx + dy * dy + dz * dz


def _dist(a: Tuple[str, float, float, float],
          b: Tuple[str, float, float, float]) -> float:
    return math.sqrt(_dist2(a, b))


# ---------------------------------------------------------------------- #
# Detector 1: extra heavy-bonds vs SMILES (mirrors _count_extra_heavy_bonds)
# ---------------------------------------------------------------------- #

def _count_extras(atoms: Sequence[Tuple[str, float, float, float]],
                  smi_pairs: Counter) -> int:
    """XYZ heavy-pair multiset minus SMILES heavy-pair multiset (positive part).

    ``smi_pairs`` is pre-computed once per (smiles, mol) and reused across all
    frames of that SMILES — saves the RDKit parse cost per frame.
    """
    if not atoms:
        return 0
    heavy = [a for a in atoms if a[0] != "H"]
    if not heavy:
        return 0
    n = len(heavy)
    xyz_pairs: Counter = Counter()
    # Pre-cache covalent radii and metal flags.
    radii = [_COV.get(a[0], 1.50) for a in heavy]
    is_metal = [a[0] in _METALS for a in heavy]
    syms = [a[0] for a in heavy]
    for i in range(n):
        si = syms[i]; ri = radii[i]; mi = is_metal[i]
        ai = heavy[i]
        for j in range(i + 1, n):
            sj = syms[j]; rj = radii[j]; mj = is_metal[j]
            tol = _TOL_METAL if (mi or mj) else _TOL_ORG
            cutoff = ri + rj + tol
            # Squared-distance early-exit
            d2 = _dist2(ai, heavy[j])
            if d2 <= cutoff * cutoff:
                xyz_pairs[tuple(sorted([si, sj]))] += 1
    n_extra = 0
    for key in set(xyz_pairs) | set(smi_pairs):
        d = xyz_pairs.get(key, 0) - smi_pairs.get(key, 0)
        if d > 0:
            n_extra += d
    return n_extra


def build_smiles_pair_counter(mol) -> Counter:
    """Pre-compute the SMILES heavy-bond pair multiset for a mol.

    Caller is expected to invoke this *once* per SMILES and pass the result
    to ``score_frame`` — turns a per-frame O(bonds) re-parse into a per-SMILES
    O(bonds) one-shot.
    """
    smi_pairs: Counter = Counter()
    if mol is None:
        return smi_pairs
    try:
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom(); b = bond.GetEndAtom()
            if a.GetAtomicNum() <= 1 or b.GetAtomicNum() <= 1:
                continue
            smi_pairs[tuple(sorted([a.GetSymbol(), b.GetSymbol()]))] += 1
    except Exception:
        pass
    return smi_pairs


# ---------------------------------------------------------------------- #
# Detector 2: M-D invariant deviation
# ---------------------------------------------------------------------- #

def _md_invariant_deviation(
    atoms: Sequence[Tuple[str, float, float, float]],
) -> float:
    """Maximum |‖M-D‖ - canonical| across all M-D pairs in the frame.

    Returns 0.0 when no metal is present (no-metal SMILES → no penalty).
    Returns the max-deviation in Å otherwise.  Only donor atoms (N/O/P/S/C/
    halogens within element-specific cutoff) are counted as bonded to the
    metal; that mirrors find_chelate_planarity's donor identification.
    """
    metals_idx = [i for i, a in enumerate(atoms) if a[0] in _METALS]
    if not metals_idx:
        return 0.0
    max_dev = 0.0
    for mi in metals_idx:
        m = atoms[mi]
        for j, a in enumerate(atoms):
            if j == mi or a[0] == "H" or a[0] in _METALS:
                continue
            cut = _DONOR_CUT.get(a[0], 2.30)
            d = _dist(m, a)
            if d < cut:
                dev = abs(d - _MD_CANONICAL)
                if dev > max_dev:
                    max_dev = dev
    return max_dev


# ---------------------------------------------------------------------- #
# Detector 3 & 4: clashes (heavy-heavy, H-H)
# ---------------------------------------------------------------------- #

def _count_clashes(
    atoms: Sequence[Tuple[str, float, float, float]],
    heavy_factor: float = 0.85,
    h_h_threshold: float = 1.6,
) -> Tuple[int, int]:
    """Return ``(n_heavy_clashes, n_h_clashes)``.

    Heavy clashes: heavy-heavy pair at d < ``heavy_factor`` · (r_vdw_i + r_vdw_j),
        excluding M-D pairs (those are coordination bonds, not clashes).
    H clashes:     H-H pair at d < ``h_h_threshold`` Å.

    Bonded heavy-heavy pairs (covalent + metal-donor) are NOT excluded by an
    explicit graph here; instead we use the vdw-factor 0.85 cutoff which sits
    well below typical bond lengths (a 1.5 Å C-C bond at 0.85·(1.7+1.7)=2.89 Å
    threshold would be a clash by definition; correct behaviour: real bonds
    are always shorter than ``Σr_vdw`` so the factor cuts them).  Net: only
    "true" non-bonded close contacts ping the counter.

    The exclusion that DOES matter: M-D pairs sit at ~2.0 Å vs typical
    Σr_vdw(M, donor) ≈ 3.5 Å → 0.85·3.5 = 2.97 Å; a 2.0 Å M-D would
    spuriously count.  We exclude any pair where one atom is a metal and the
    other is a donor (covered by ``_DONOR_CUT``).
    """
    heavy: List[Tuple[int, str, float, float, float]] = []
    hs: List[Tuple[int, str, float, float, float]] = []
    for i, a in enumerate(atoms):
        if a[0] == "H":
            hs.append((i, *a))
        else:
            heavy.append((i, *a))
    n_heavy = len(heavy)
    n_heavy_clash = 0
    for i in range(n_heavy):
        _, si, xi, yi, zi = heavy[i]
        ri = _VDW.get(si, _VDW_DEFAULT_HEAVY)
        mi = si in _METALS
        for j in range(i + 1, n_heavy):
            _, sj, xj, yj, zj = heavy[j]
            rj = _VDW.get(sj, _VDW_DEFAULT_HEAVY)
            mj = sj in _METALS
            # Skip coordination M-D pairs.
            if mi != mj:  # exactly one is metal
                donor_sym = sj if mi else si
                cut = _DONOR_CUT.get(donor_sym, 2.30)
                d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                if d2 < cut * cut:
                    continue
            else:
                d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            cutoff = heavy_factor * (ri + rj)
            if d2 < cutoff * cutoff:
                n_heavy_clash += 1
    # H-H pairs
    n_h = len(hs)
    n_h_clash = 0
    h_thresh_sq = h_h_threshold * h_h_threshold
    for i in range(n_h):
        _, _, xi, yi, zi = hs[i]
        for j in range(i + 1, n_h):
            _, _, xj, yj, zj = hs[j]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if d2 < h_thresh_sq:
                n_h_clash += 1
    return n_heavy_clash, n_h_clash


# ---------------------------------------------------------------------- #
# Detector 5: tridentate planarity
# ---------------------------------------------------------------------- #

def _build_nonmetal_graph(
    atoms: Sequence[Tuple[str, float, float, float]],
) -> Dict[int, set]:
    n = len(atoms)
    nbrs: Dict[int, set] = {
        i: set() for i in range(n) if atoms[i][0] not in _METALS
    }
    keys = list(nbrs.keys())
    for ii, i in enumerate(keys):
        ai = atoms[i]
        for j in keys[ii + 1:]:
            d = _dist(ai, atoms[j])
            if d >= _BOND_CUT:
                continue
            cut_pair = (
                1.30 if (atoms[i][0] == "H" or atoms[j][0] == "H") else _BOND_CUT
            )
            if d < cut_pair:
                nbrs[i].add(j); nbrs[j].add(i)
    return nbrs


def _find_donors(
    atoms: Sequence[Tuple[str, float, float, float]],
    metal_idx: int,
) -> List[int]:
    m = atoms[metal_idx]
    donors: List[int] = []
    for i, a in enumerate(atoms):
        if i == metal_idx or a[0] == "H" or a[0] in _METALS:
            continue
        cut = _DONOR_CUT.get(a[0], 2.30)
        if _dist(m, a) < cut:
            donors.append(i)
    return donors


def _cluster_donors(donor_idxs: List[int],
                    nbrs: Dict[int, set]) -> List[List[int]]:
    seen: set = set()
    groups: List[List[int]] = []
    donor_set = set(donor_idxs)
    for d in donor_idxs:
        if d in seen:
            continue
        stack = [d]; comp_donors: List[int] = []
        local_seen = {d}
        while stack:
            cur = stack.pop()
            if cur in donor_set:
                comp_donors.append(cur)
            for nb in nbrs.get(cur, ()):
                if nb not in local_seen:
                    local_seen.add(nb); stack.append(nb)
        for c in comp_donors:
            seen.add(c)
        groups.append(sorted(comp_donors))
    return groups


def _plane_oop(metal_xyz: Tuple[float, float, float],
               donors_xyz: List[Tuple[float, float, float]]) -> Optional[float]:
    if len(donors_xyz) < 3:
        return None
    d0 = donors_xyz[0]
    v1 = [donors_xyz[1][k] - d0[k] for k in range(3)]
    v2 = [donors_xyz[2][k] - d0[k] for k in range(3)]
    n = [v1[1] * v2[2] - v1[2] * v2[1],
         v1[2] * v2[0] - v1[0] * v2[2],
         v1[0] * v2[1] - v1[1] * v2[0]]
    norm = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    if norm < 1e-9:
        return None
    n = [c / norm for c in n]
    v_m = [metal_xyz[k] - d0[k] for k in range(3)]
    return abs(v_m[0] * n[0] + v_m[1] * n[1] + v_m[2] * n[2])


def _count_tridentate_oop_fails(
    atoms: Sequence[Tuple[str, float, float, float]],
    threshold: float = 0.30,
) -> int:
    metals = [i for i, a in enumerate(atoms) if a[0] in _METALS]
    if not metals:
        return 0
    nbrs = _build_nonmetal_graph(atoms)
    n_fail = 0
    for mi in metals:
        donors = _find_donors(atoms, mi)
        if len(donors) < 3:
            continue
        groups = _cluster_donors(donors, nbrs)
        m_xyz = (atoms[mi][1], atoms[mi][2], atoms[mi][3])
        for grp in groups:
            if len(grp) != 3:
                continue
            d_xyz = [(atoms[d][1], atoms[d][2], atoms[d][3]) for d in grp]
            oop = _plane_oop(m_xyz, d_xyz)
            if oop is not None and oop > threshold:
                n_fail += 1
    return n_fail


# ---------------------------------------------------------------------- #
# Detector 6: bond-length outliers
# ---------------------------------------------------------------------- #

def _count_bond_outliers(
    atoms: Sequence[Tuple[str, float, float, float]],
    long_factor: float = 1.30,
) -> int:
    """Heavy-heavy bonds whose length exceeds ``long_factor · (r_i + r_j + tol)``.

    "Bond" is defined here as: a heavy-heavy pair at d ≤ Σr_cov + tol (the
    extras-detector criterion).  An outlier is a bond above 1.30× that cutoff
    — typically a stretched / partially-broken bond from a UFF blow-up.
    """
    heavy = [(i, a) for i, a in enumerate(atoms) if a[0] != "H"]
    n = len(heavy)
    n_out = 0
    for ii in range(n):
        i, ai = heavy[ii]
        ri = _COV.get(ai[0], 1.50)
        mi = ai[0] in _METALS
        for jj in range(ii + 1, n):
            j, aj = heavy[jj]
            rj = _COV.get(aj[0], 1.50)
            mj = aj[0] in _METALS
            tol = _TOL_METAL if (mi or mj) else _TOL_ORG
            cutoff = ri + rj + tol
            d = _dist(ai, aj)
            if d <= cutoff and d > long_factor * cutoff:
                n_out += 1
    return n_out


# ---------------------------------------------------------------------- #
# Public scoring entry point
# ---------------------------------------------------------------------- #

def score_frame(
    xyz_text: str,
    smiles: str,
    mol=None,
    weights: Optional[Dict[str, float]] = None,
    smi_pair_counter: Optional[Counter] = None,
) -> Tuple[float, Dict[str, float]]:
    """Composite score for a single XYZ frame.  Lower = better.

    Returns ``(composite, breakdown)`` where ``breakdown`` is the per-metric
    raw count (NOT weight-multiplied), useful for adaptive-memory logging
    + per-class diagnostics.

    ``mol`` and ``smi_pair_counter`` are optional pre-computed inputs:
        - If ``smi_pair_counter`` is provided, it overrides any ``mol`` re-parse.
        - If ``mol`` is provided but no counter, build counter from mol.
        - If neither is provided, parse SMILES via RDKit (if available); on
          failure the extras detector is skipped (counter = empty Counter →
          all pairs in XYZ count as "extra" — this is correct fail-safe
          behaviour: a SMILES we can't parse cannot validate against XYZ).

    Total cost on a 60-atom complex with weights = DEFAULT_WEIGHTS: <50 ms.
    """
    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update(weights)

    atoms = _parse_xyz(xyz_text)
    if not atoms:
        # Unparseable XYZ → infinite cost so any successful frame wins.
        return float("inf"), {
            "extras": float("inf"), "md_invariant": float("inf"),
            "heavy_clash": float("inf"), "h_clash": float("inf"),
            "tridentate_oop": float("inf"), "bond_outliers": float("inf"),
        }

    # Pair counter resolution
    if smi_pair_counter is None:
        if mol is None:
            try:
                from rdkit import Chem  # type: ignore
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is not None:
                    try:
                        mol.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
            except Exception:
                mol = None
        smi_pair_counter = build_smiles_pair_counter(mol)

    extras = _count_extras(atoms, smi_pair_counter)
    md_dev = _md_invariant_deviation(atoms)
    heavy_clash, h_clash = _count_clashes(atoms)
    tri_oop = _count_tridentate_oop_fails(atoms)
    bond_out = _count_bond_outliers(atoms)

    breakdown: Dict[str, float] = {
        "extras":         float(extras),
        "md_invariant":   float(md_dev),
        "heavy_clash":    float(heavy_clash),
        "h_clash":        float(h_clash),
        "tridentate_oop": float(tri_oop),
        "bond_outliers":  float(bond_out),
    }
    composite = sum(w.get(k, 0.0) * v for k, v in breakdown.items())
    return composite, breakdown


# ---------------------------------------------------------------------- #
# Polyhedron fingerprint (canonical-form, scipy-Hungarian)
# ---------------------------------------------------------------------- #

def _kabsch(P, Q):
    import numpy as np  # local import — only when classifier runs
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def _best_perm(donor_unit, ref_unit):
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    n = len(ref_unit)
    cost = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            d = ref_unit[i] - donor_unit[j]
            cost[i, j] = float(np.dot(d, d))
    row, col = linear_sum_assignment(cost)
    perm = [0] * n
    for r, c in zip(row, col):
        perm[r] = int(c)
    err = float(sum(cost[r, c] for r, c in zip(row, col)) / max(1, n))
    return perm, err


def _normalise(vs):
    import numpy as np
    arr = np.asarray(vs, dtype=float)
    norms = np.linalg.norm(arr, axis=1, keepdims=True)
    norms[norms < 1e-9] = 1.0
    return arr / norms


# Lazy-init geom unit-vertex dict.  Avoids numpy import at module load.
_GEOM_UNIT_CACHE: Optional[Dict[str, "object"]] = None
_CN_TO_GEOMS_CACHE: Optional[Dict[int, List[str]]] = None


def _geom_unit():
    global _GEOM_UNIT_CACHE, _CN_TO_GEOMS_CACHE
    if _GEOM_UNIT_CACHE is None:
        _GEOM_UNIT_CACHE = {
            g: _normalise(v) for g, v in _GEOM_VECTORS_RAW.items()
        }
        cn_map: Dict[int, List[str]] = defaultdict(list)
        for g, arr in _GEOM_UNIT_CACHE.items():
            cn_map[len(arr)].append(g)
        _CN_TO_GEOMS_CACHE = dict(cn_map)
    return _GEOM_UNIT_CACHE, _CN_TO_GEOMS_CACHE


def _classify_polyhedron(donor_xyz, metal_xyz, cn_hint: int):
    """Return ``(best_geom, mean_sq_err, perm)``.

    Classifier with single random-start (faster than the 6-start version in
    find_isomer_coverage_named — for fingerprinting we do NOT need the
    absolute best alignment, only a stable canonical-form across small noise).
    """
    import numpy as np
    geom_unit, cn_to_geoms = _geom_unit()
    candidates = cn_to_geoms.get(cn_hint, [])
    if not candidates:
        return ("?", float("inf"), list(range(len(donor_xyz))))

    centred = donor_xyz - metal_xyz
    norms = np.linalg.norm(centred, axis=1, keepdims=True)
    norms[norms < 1e-9] = 1.0
    donor_unit = centred / norms

    best_geom = candidates[0]
    best_err = float("inf")
    best_perm = list(range(len(donor_unit)))
    for g in candidates:
        ref = geom_unit[g]
        if len(ref) != len(donor_unit):
            continue
        # Single Kabsch-Hungarian iteration with two refinement rounds.
        Pr = donor_unit
        for _ in range(2):
            perm, _err = _best_perm(Pr, ref)
            paired = donor_unit[perm]
            R = _kabsch(paired, ref)
            Pr = donor_unit @ R.T
        perm, err = _best_perm(Pr, ref)
        if err < best_err:
            best_err = err; best_geom = g; best_perm = perm
    return best_geom, best_err, best_perm


def _build_canonical_form(geom: str, types: List[str]) -> tuple:
    """Mirror find_isomer_coverage_named._build_canonical_form (subset of
    geometries that the classifier emits).  Donor symbols are bare elements.
    """
    t = list(types)
    if geom == 'OH':
        pairs = tuple(sorted([
            tuple(sorted([t[0], t[1]])),
            tuple(sorted([t[2], t[3]])),
            tuple(sorted([t[4], t[5]])),
        ]))
        return ('OH', pairs)
    if geom == 'SQ':
        pairs = tuple(sorted([
            tuple(sorted([t[0], t[1]])),
            tuple(sorted([t[2], t[3]])),
        ]))
        return ('SQ', pairs)
    if geom == 'TBP':
        axial = tuple(sorted([t[0], t[1]]))
        eq = tuple(sorted([t[2], t[3], t[4]]))
        return ('TBP', axial, eq)
    if geom == 'SP':
        basal = tuple(sorted([
            tuple(sorted([t[1], t[3]])),
            tuple(sorted([t[2], t[4]])),
        ]))
        return ('SP', t[0], basal)
    if geom == 'PBP':
        axial = tuple(sorted([t[0], t[1]]))
        eq = tuple(sorted([t[2], t[3], t[4], t[5], t[6]]))
        return ('PBP', axial, eq)
    if geom == 'TH':
        return ('TH', tuple(sorted(t)))
    if geom == 'LIN':
        return ('LIN', tuple(sorted(t)))
    if geom == 'TP':
        return ('TP', tuple(sorted(t)))
    if geom == 'TS':
        trans = tuple(sorted([t[0], t[1]]))
        return ('TS', trans, t[2])
    if geom == 'SAP':
        trans_pairs = tuple(sorted([
            tuple(sorted([t[0], t[6]])),
            tuple(sorted([t[1], t[7]])),
            tuple(sorted([t[2], t[4]])),
            tuple(sorted([t[3], t[5]])),
        ]))
        return ('SAP', trans_pairs)
    return (geom, tuple(sorted(t)))


def fingerprint(xyz_text: str) -> str:
    """Canonical-form fingerprint string for an XYZ frame.

    Strategy:
        1. Parse XYZ.
        2. For each metal atom, identify donors + classify polyhedron.
        3. Combine the per-metal canonical-form tuples (sorted by metal symbol)
           into a single string → final fingerprint.
        4. Fall back to a heavy-atom signature hash if no metal exists OR
           classification fails (low-cn etc.).

    Output is *always* a non-empty string so dispatcher grouping is total.
    Output is stable across identical structures up to small numeric noise
    (donors are unit-normalised before classification; Kabsch is rotation-
    invariant).
    """
    atoms = _parse_xyz(xyz_text)
    if not atoms:
        return "fp:empty"

    metal_idxs = [i for i, a in enumerate(atoms) if a[0] in _METALS]
    if not metal_idxs:
        return "fp:no_metal:" + _heavy_signature(atoms)

    try:
        import numpy as np
    except Exception:
        return "fp:nonumpy:" + _heavy_signature(atoms)

    parts: List[Tuple[str, tuple]] = []
    classified_any = False
    for mi in metal_idxs:
        donors = _find_donors(atoms, mi)
        cn = len(donors)
        if cn < 2:
            parts.append((atoms[mi][0], ("low-cn", cn)))
            continue
        donor_xyz = np.array([(atoms[d][1], atoms[d][2], atoms[d][3])
                              for d in donors])
        metal_xyz = np.array([atoms[mi][1], atoms[mi][2], atoms[mi][3]])
        donor_elems = [atoms[d][0] for d in donors]
        try:
            geom, _err, perm = _classify_polyhedron(donor_xyz, metal_xyz, cn)
            ordered_types = [donor_elems[perm[k]] for k in range(cn)]
            cf = _build_canonical_form(geom, ordered_types)
            parts.append((atoms[mi][0], cf))
            classified_any = True
        except Exception:
            parts.append((atoms[mi][0], ("err", cn)))

    if not classified_any:
        return "fp:cnonly:" + _heavy_signature(atoms)

    parts.sort(key=lambda t: (t[0], repr(t[1])))
    return "fp:" + "|".join(f"{m}={cf!r}" for m, cf in parts)


def _heavy_signature(atoms: Sequence[Tuple[str, float, float, float]]) -> str:
    """SHA-1 of the rounded heavy-atom positions — fallback fingerprint."""
    body: List[str] = []
    for sym, x, y, z in atoms:
        if sym == "H":
            continue
        body.append(f"{sym} {round(x, 2):.2f} {round(y, 2):.2f} {round(z, 2):.2f}")
    body.sort()
    return hashlib.sha1("\n".join(body).encode("utf-8")).hexdigest()[:16]
