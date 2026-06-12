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


def rank_isomers(isomers: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Re-order ``[(xyz_string, label), ...]`` best (most crystal-like) first.

    Deterministic and stable (ties keep enumeration order).  Never drops, adds,
    or mutates a structure — only reorders.  No-op for <2 frames or when
    ``DELFIN_NO_FRAME_RANK=1``.
    """
    if os.environ.get("DELFIN_NO_FRAME_RANK", "0") == "1":
        return isomers
    if not isinstance(isomers, list) or len(isomers) < 2:
        return isomers
    try:
        scored = []
        for idx, item in enumerate(isomers):
            xyz = item[0] if isinstance(item, (tuple, list)) and item else ""
            sc = _frame_score(_parse_xyz(xyz)) if isinstance(xyz, str) else -1e9
            scored.append((idx, sc, item))
        # sort by descending score; idx as stable tiebreaker (enumeration order)
        scored.sort(key=lambda e: (-e[1], e[0]))
        return [e[2] for e in scored]
    except Exception:
        # Ranking must never break the build — fall back to original order.
        return isomers
