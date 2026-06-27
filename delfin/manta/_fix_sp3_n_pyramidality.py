"""F25 — Targeted sp3-N pyramidality fixer.

For sp3 nitrogen atoms whose three substituents are flattened (sum of three
N-X-N' angles around N > planarity_threshold, default 348°, i.e. less than
12° pyramidalization), bend an H substituent out of the local trigonal plane
toward the ideal pyramidal apex (target sum ~328° = 3 × 109.47°).

Common root cause: UFF flattens sp3 N adjacent to aromatic rings or carbonyls
via spurious π-conjugation, but the chemistry (free amine, sp3 hybridization
asserted by RDKit) demands a pyramidal geometry.

Touch-rules (chemical safety):
  - Only sp3 N (RDKit hybridization).
  - Exactly 3 heavy/H neighbours (proper 3-coord N — discards 4-coord
    ammonium / 2-coord imine).
  - Skip if N is bonded to a metal (M-N geometry obeys coordination, not
    free-amine geometry).
  - Skip if N is part of an amide (one neighbour is C=O — sp2 by resonance
    even though RDKit may mark sp3 / amide-N).
  - Only H substituents are displaced (heavy chain remains rigid).  If the
    sp3 N has no H, the fix is a no-op (e.g. trimethylamine).

Bend amount is computed analytically: for an ideal pyramidal apex at OOP
distance ``h`` from the trigonal plane through three neighbours, each at
distance ``r`` from N, the sum of the three N-X angles equals
    ``S = 3 * arccos((r^2 cos(120°) - h^2) / (r^2 + h^2))``.
We solve numerically for the apex height that yields S = target.  The
H atoms are then translated along the plane normal by exactly that amount
(N and heavy neighbours stay put — only H moves).  This is conservative:
the heavy-atom geometry is unchanged, so the rest of the molecule is
untouched.

Doctrine (consistent with B3/B4 sibling fixers):
  - Universal: covalent radii + atomic number range for metal detection.
  - Per-violation rollback if a NEW non-bonded clash appears.
  - Co-evolution safe: post-UFF, no embedding disturbance.
  - Diversity safe: per-conformer; n_frames invariant.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Vendored helpers (kept in sync with delfin/_coord_angle_corrector.py)
# ---------------------------------------------------------------------------

_COV_RADII: Dict[str, float] = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
}

_METAL_Z_RANGES = (
    set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
    | set(range(89, 104))
    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51, 55, 56, 81, 82, 83}
)
_Z_BY_SYMBOL: Dict[str, int] = {
    "H": 1, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
    "As": 33, "Se": 34, "Br": 35, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47,
    "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Cs": 55,
    "Ba": 56, "La": 57, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76,
    "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}


def _is_metal_sym(sym: str) -> bool:
    z = _Z_BY_SYMBOL.get(sym)
    return z is not None and z in _METAL_Z_RANGES


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse XYZ → (symbols, Nx3 positions, raw_lines).  Header / comment
    lines are preserved in raw_lines for the round-trip writer."""
    syms: List[str] = []
    pts: List[np.ndarray] = []
    lines = xyz_str.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append(np.array([float(m.group(2)), float(m.group(3)),
                                 float(m.group(4))]))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(orig_lines: List[str], syms: List[str],
                positions: np.ndarray) -> str:
    """Reformat XYZ preserving header but rewriting atom coordinates."""
    out: List[str] = []
    atom_i = 0
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(
                f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}"
            )
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if orig_lines else "")


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _angle_deg(p_center: np.ndarray, p_a: np.ndarray,
               p_b: np.ndarray) -> Optional[float]:
    v1 = p_a - p_center
    v2 = p_b - p_center
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return None
    cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos)))


def _sum_three_angles_at(p_center: np.ndarray, p_a: np.ndarray,
                         p_b: np.ndarray, p_c: np.ndarray) -> Optional[float]:
    a_ab = _angle_deg(p_center, p_a, p_b)
    a_ac = _angle_deg(p_center, p_a, p_c)
    a_bc = _angle_deg(p_center, p_b, p_c)
    if a_ab is None or a_ac is None or a_bc is None:
        return None
    return a_ab + a_ac + a_bc


def _is_amide_n(mol, n_atom) -> bool:
    """Return True if N is directly bonded to a carbonyl C (C=O / C=S),
    i.e. classic amide / thioamide / carbamate / urea nitrogen.  These
    are sp2 by π-conjugation even when RDKit's hybridization is sp3."""
    try:
        from rdkit import Chem
    except Exception:
        return False
    for nb in n_atom.GetNeighbors():
        if nb.GetSymbol() != "C":
            continue
        for nb2 in nb.GetNeighbors():
            if nb2.GetIdx() == n_atom.GetIdx():
                continue
            bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
            if bond is None:
                continue
            if bond.GetBondType() == Chem.BondType.DOUBLE and \
               nb2.GetSymbol() in ("O", "S", "N"):
                return True
    return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def detect_sp3_n_planar_violations(xyz: str, mol,
                                   threshold_deg: float = 348.0) -> List[Dict]:
    """Return list of sp3-N atoms whose three-angle sum exceeds threshold.

    Each entry: {"n_idx": int, "neighbors": [a, b, c],
                 "current_sum_deg": float, "angles": (a_ab, a_ac, a_bc)}.

    Filters: requires RDKit sp3 hybridization, exactly 3 neighbours,
    no metal neighbour, not amide-N.
    """
    try:
        from rdkit import Chem
    except Exception:
        return []
    syms, pts, _ = _parse_xyz(xyz)
    if pts.shape[0] == 0:
        return []
    if mol is None or mol.GetNumAtoms() != len(syms):
        return []

    violations: List[Dict] = []
    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "N":
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        # 3-coordinate sp3 N only — exclude 4-coord ammonium etc.
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        # Exclude metal-coordinated N (different chemistry).
        if any(nb.GetSymbol() in _Z_BY_SYMBOL and _is_metal_sym(nb.GetSymbol())
               for nb in nbrs):
            continue
        # Exclude amide-N (sp2 by resonance).
        if _is_amide_n(mol, atom):
            continue
        nbr_idxs = [nb.GetIdx() for nb in nbrs]
        sum_deg = _sum_three_angles_at(
            pts[atom_idx], pts[nbr_idxs[0]], pts[nbr_idxs[1]], pts[nbr_idxs[2]],
        )
        if sum_deg is None:
            continue
        if sum_deg > threshold_deg:
            a_ab = _angle_deg(pts[atom_idx], pts[nbr_idxs[0]], pts[nbr_idxs[1]])
            a_ac = _angle_deg(pts[atom_idx], pts[nbr_idxs[0]], pts[nbr_idxs[2]])
            a_bc = _angle_deg(pts[atom_idx], pts[nbr_idxs[1]], pts[nbr_idxs[2]])
            violations.append({
                "n_idx": atom_idx,
                "neighbors": nbr_idxs,
                "current_sum_deg": sum_deg,
                "angles": (a_ab, a_ac, a_bc),
            })
    return violations


def _solve_apex_height(r: float, target_sum_deg: float) -> float:
    """Return apex height h such that for three N-X bonds of length ``r``
    arranged symmetrically with apex (N) at height h above the centroid of
    the X-triangle (X-X-X equilateral edge), the sum of the three N-X-N
    angles equals ``target_sum_deg``.

    Geometry: place 3 X atoms in a plane at distance ``d_xy`` from centroid.
    Then ``r^2 = d_xy^2 + h^2``.  Pairwise X-X distance = ``sqrt(3) * d_xy``.
    The angle X-N-X' is then
        ``alpha = arccos((r^2 - 1.5 * d_xy^2) / r^2)
                = arccos(1 - 1.5 * d_xy^2 / r^2)``.
    Sum = 3 * alpha.  Solve numerically by bisection in h ∈ [0, r].

    For h = 0 (planar) → d_xy = r → cos(alpha) = -0.5 → alpha = 120° →
    sum = 360°.  For h = r (collapsed) → d_xy = 0 → cos(alpha) = 1 →
    alpha = 0° → sum = 0°.  Monotonic in h, so bisection converges.
    """
    target_rad = math.radians(target_sum_deg / 3.0)  # per-angle target
    cos_target = math.cos(target_rad)
    # cos(alpha) = 1 - 1.5 * d_xy^2 / r^2 = 1 - 1.5 * (r^2 - h^2) / r^2
    #            = 1 - 1.5 + 1.5 * h^2 / r^2 = -0.5 + 1.5 * (h/r)^2
    # So (h/r)^2 = (cos_target + 0.5) / 1.5
    val = (cos_target + 0.5) / 1.5
    if val < 0.0:
        return 0.0
    if val > 1.0:
        return r
    return r * math.sqrt(val)


def fix_sp3_n_pyramidality(xyz: str, mol,
                           planarity_threshold_deg: float = 348.0,
                           target_sum_deg: float = 328.0,
                           ) -> Tuple[str, Dict]:
    """Detect flat-too-much sp3 N atoms and pyramidalize by bending H
    substituents out of the local trigonal plane.

    For each violating sp3 N:
        1. Identify H substituents (only these will move; heavy chain rigid).
        2. If no H substituent → skip (e.g. trimethylamine).
        3. Compute mean N-H bond length r_NH among the H substituents.
        4. Solve for apex height h_target such that the symmetric pyramid
           with r = r_NH yields ``target_sum_deg``.
        5. Compute the current apex height of N above the three-neighbour
           plane (signed along the plane normal).  Determine displacement
           ``delta = h_target - h_current`` along the normal.  Choose the
           sign that pushes H atoms TOWARD the side opposite from N's
           current apex (lone-pair side reinforced).
        6. Translate each H substituent by ``-delta * n_plane`` (i.e. push
           it through the plane to deepen the pyramid).
        7. Rollback the whole N if any moved H gains a NEW non-bonded
           clash with non-H heavy atoms within 0.85 · Σr_vdw.

    Args:
        xyz: DELFIN XYZ string.
        mol: RDKit Mol matching the XYZ atom order.
        planarity_threshold_deg: trigger; sum > this means "too flat".
        target_sum_deg: target sum after fix (default 328° ≈ 3·109.47°).

    Returns:
        (new_xyz, report).  report keys: ``n_violations``, ``n_fixed``,
        ``max_correction_oop_A``, ``per_n`` (list per violation).
    """
    syms, pts, orig_lines = _parse_xyz(xyz)
    report: Dict = {"n_violations": 0, "n_fixed": 0,
                    "max_correction_oop_A": 0.0, "per_n": []}
    if pts.shape[0] == 0 or mol is None:
        return xyz, report
    if mol.GetNumAtoms() != len(syms):
        return xyz, report

    violations = detect_sp3_n_planar_violations(
        xyz, mol, threshold_deg=planarity_threshold_deg,
    )
    report["n_violations"] = len(violations)
    if not violations:
        return xyz, report

    new_pts = pts.copy()
    max_disp = 0.0

    for viol in violations:
        n_idx = viol["n_idx"]
        nbr_idxs = viol["neighbors"]

        # H substituents only — heavy chain stays rigid.
        h_nbrs = [i for i in nbr_idxs if syms[i] == "H"]
        if not h_nbrs:
            report["per_n"].append({
                "n_idx": n_idx, "fixed": False,
                "reason": "no_H_substituent",
                "current_sum_deg": viol["current_sum_deg"],
            })
            continue

        # Mean N-H bond length (target apex height parameter).
        r_NH = float(np.mean([
            np.linalg.norm(new_pts[h] - new_pts[n_idx]) for h in h_nbrs
        ]))
        if r_NH < 0.5 or r_NH > 1.5:
            # pathological — skip
            report["per_n"].append({
                "n_idx": n_idx, "fixed": False,
                "reason": "pathological_NH_length",
                "current_sum_deg": viol["current_sum_deg"],
            })
            continue

        # Plane through the three neighbours.
        p_n = new_pts[n_idx]
        p_a = new_pts[nbr_idxs[0]]
        p_b = new_pts[nbr_idxs[1]]
        p_c = new_pts[nbr_idxs[2]]
        normal_raw = np.cross(p_b - p_a, p_c - p_a)
        nn = float(np.linalg.norm(normal_raw))
        if nn < 1e-9:
            report["per_n"].append({
                "n_idx": n_idx, "fixed": False,
                "reason": "degenerate_plane",
                "current_sum_deg": viol["current_sum_deg"],
            })
            continue
        n_plane = normal_raw / nn

        # Apex height of N above the neighbour plane (signed).
        centroid = (p_a + p_b + p_c) / 3.0
        h_current = float(np.dot(p_n - centroid, n_plane))

        # Target apex height for symmetric pyramid with r = r_NH.
        h_target = _solve_apex_height(r_NH, target_sum_deg)

        # If N is already on one side of the plane, we want to DEEPEN that
        # side (i.e. drive H atoms through to the opposite side).  If N
        # straddles the plane (|h_current| < 1e-3), pick the +n_plane side.
        if abs(h_current) < 1e-3:
            sign_h = +1.0
        else:
            sign_h = math.copysign(1.0, h_current)

        # The H atoms should be on the OPPOSITE side of the plane from N
        # (lone-pair-N apex, H going down).  Required H apex height (signed):
        #   h_H_target = -sign_h * h_target.
        # Current H apex height (per H) = projection onto signed normal.
        candidate_pts = new_pts.copy()
        for h_idx in h_nbrs:
            h_cur_signed = float(np.dot(new_pts[h_idx] - centroid, n_plane))
            target_signed = -sign_h * h_target
            # Translate H along plane normal so its signed projection becomes
            # target_signed.  This preserves the H's in-plane position
            # (so N-H bond length changes slightly but stays within tolerance).
            delta = target_signed - h_cur_signed
            candidate_pts[h_idx] = new_pts[h_idx] + delta * n_plane
            max_disp = max(max_disp, abs(delta))

        # Clash-rollback gate: any NEW heavy contact for a moved H?
        if _introduces_new_clash(new_pts, candidate_pts, syms, h_nbrs):
            report["per_n"].append({
                "n_idx": n_idx, "fixed": False,
                "reason": "rollback_clash",
                "current_sum_deg": viol["current_sum_deg"],
            })
            continue

        # Verify sum-of-angles actually improved.
        new_sum = _sum_three_angles_at(
            candidate_pts[n_idx], candidate_pts[nbr_idxs[0]],
            candidate_pts[nbr_idxs[1]], candidate_pts[nbr_idxs[2]],
        )
        if new_sum is None or new_sum >= viol["current_sum_deg"]:
            report["per_n"].append({
                "n_idx": n_idx, "fixed": False,
                "reason": "no_improvement",
                "current_sum_deg": viol["current_sum_deg"],
                "new_sum_deg": new_sum,
            })
            continue

        new_pts = candidate_pts
        report["n_fixed"] += 1
        report["per_n"].append({
            "n_idx": n_idx, "fixed": True,
            "current_sum_deg": viol["current_sum_deg"],
            "new_sum_deg": new_sum,
            "h_moved": h_nbrs,
        })

    report["max_correction_oop_A"] = max_disp
    if report["n_fixed"] == 0:
        return xyz, report
    return _format_xyz(orig_lines, syms, new_pts), report


def _introduces_new_clash(old_pts: np.ndarray, new_pts: np.ndarray,
                          syms: List[str], moved_idxs: List[int],
                          factor: float = 0.85) -> bool:
    """Return True if any moved atom in new_pts gains a non-bonded contact
    with a stationary heavy atom at d < factor · Σr_vdw that wasn't present
    in old_pts.

    Bonded pairs (d < 1.10 · Σr_cov in old_pts) are excluded.
    """
    n = len(syms)
    for i in moved_idxs:
        ri_v = _VDW_RADII.get(syms[i], 1.20 if syms[i] == "H" else 1.70)
        ri_c = _COV_RADII.get(syms[i], 0.31 if syms[i] == "H" else 1.50)
        for j in range(n):
            if j == i or j in moved_idxs:
                continue
            sj = syms[j]
            if sj == "H":
                continue
            rj_v = _VDW_RADII.get(sj, 1.70)
            rj_c = _COV_RADII.get(sj, 1.50)
            d_old = float(np.linalg.norm(old_pts[i] - old_pts[j]))
            if d_old < 1.10 * (ri_c + rj_c):
                continue  # bonded
            d_new = float(np.linalg.norm(new_pts[i] - new_pts[j]))
            thr = factor * (ri_v + rj_v)
            if d_new < thr <= d_old:
                return True
    return False
