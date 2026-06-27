"""SP2N-PLANARIZE — universal sp2-nitrogen planarisation fixer.

Sibling of ``_fix_sp3_n_pyramidality`` (F25) but for the OPPOSITE failure
mode: an sp2 nitrogen that the embed / force-field has **pyramidalised** and
desymmetrised (the classic broken nitro group: ``N`` pushed ~1 Å out of the
``C-O-O`` plane, one ``N-O`` stretched to ~1.43 Å while the other collapses to
~1.12 Å, and ``O-N-O`` opening to ~107° instead of the planar-sp2 ~125°).

The primary target is the **C-nitro group** (organic ``-NO2`` substituent on a
carbon): a nitrogen bonded to exactly two oxygens and one carbon, total degree
three.  This group is planar sp2 (resonance-delocalised, D-bond order 1.5 to
each O) and must be coplanar with its three neighbours.  ``conformer_enum`` /
UFF / MMFF leave the isolated ligand fine, but the **complex-assembly** path
(metal-sphere construction + global heavy-atom relax) can distort it.

Optionally (``include_amide_imine=True``) the same planarisation is applied to
other genuinely-planar acyclic sp2 N that has been pyramidalised:
  - amide / carbamate / urea N (one neighbour is a carbonyl C=O / C=S);
  - imine / Schiff-base N (N=C, two-coordinate is skipped — needs 3 neighbours).
These are *geometry-only* (positions are projected onto the best-fit plane);
no bond-length retargeting is applied to them (only nitro gets N-O = 1.22 Å).

Touch-rules (chemical safety, mirrors the F25 / sibling doctrine):
  - Universal: pure graph + element rules, never SMILES-specific.
  - Skip if N is bonded to a metal (M-N geometry obeys coordination).
  - Skip ring N (aromatic / heterocyclic planarity is handled elsewhere and
    its in-plane geometry is fixed by the ring, not by this fixer).
  - Per-violation rollback if the projection introduces a NEW non-bonded
    clash with a stationary heavy atom (never-worse).
  - Deterministic: closed-form projection, no random state, no FF call.
  - Never emits a non-finite coordinate (guarded; rollback to input on NaN).

Construction (per nitro N with neighbours C, O1, O2):
  1. Best-fit plane through {C, O1, O2}: normal = (O1-C) x (O2-C).
  2. Move N into that plane (drop its out-of-plane component), then place it
     so that it sits at the apex of a symmetric, planar NO2 with the C-N
     direction fixed: keep the C-N bond direction (in-plane) and the C-N bond
     length, set the two N-O bonds to ``no_target`` (1.22 Å) at ±(O-N-O)/2
     about the C->N axis, all in the C-O-O plane.
  3. The two O atoms are repositioned to the symmetric planar geometry; the
     N stays bonded to C at its original C-N length and direction.

This keeps the carbon skeleton rigid and only moves {N, O1, O2}.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Vendored helpers (kept in sync with delfin/_fix_sp3_n_pyramidality.py)
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
    out: List[str] = []
    atom_i = 0
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if orig_lines else "")


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


def _oop_distance(p_n: np.ndarray, p_a: np.ndarray, p_b: np.ndarray,
                  p_c: np.ndarray) -> Optional[float]:
    """Signed-magnitude out-of-plane distance of ``p_n`` from the plane
    through ``p_a, p_b, p_c``."""
    nrm = np.cross(p_b - p_a, p_c - p_a)
    nn = float(np.linalg.norm(nrm))
    if nn < 1e-9:
        return None
    return abs(float(np.dot(p_n - p_a, nrm / nn)))


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def detect_nitro_groups(mol) -> List[Dict]:
    """Return every C-nitro nitrogen (N bonded to exactly 2 O + 1 C, degree 3).

    Each entry: {"n_idx", "o_idxs": [o1, o2], "c_idx"}.

    Filters: not metal-bonded, not in a ring (a ring nitro is exotic and its
    in-plane geometry is constrained by the ring; we leave it alone).
    """
    try:
        from rdkit import Chem  # noqa: F401
    except Exception:
        return []
    out: List[Dict] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        if atom.IsInRing():
            continue
        o_idxs = [nb.GetIdx() for nb in nbrs if nb.GetSymbol() == "O"]
        c_idxs = [nb.GetIdx() for nb in nbrs if nb.GetSymbol() == "C"]
        if len(o_idxs) != 2 or len(c_idxs) != 1:
            continue
        # Skip if any neighbour (or N itself) is metal-bonded.
        if any(_is_metal_sym(nb.GetSymbol()) for nb in nbrs):
            continue
        # An O bonded to a metal => coordinated nitro/nitrite, leave to coord.
        metal_touch = False
        for o_idx in o_idxs:
            for onb in mol.GetAtomWithIdx(o_idx).GetNeighbors():
                if _is_metal_sym(onb.GetSymbol()):
                    metal_touch = True
                    break
            if metal_touch:
                break
        if metal_touch:
            continue
        out.append({"n_idx": atom.GetIdx(), "o_idxs": o_idxs,
                    "c_idx": c_idxs[0]})
    return out


def detect_planar_sp2n_groups(mol, include_amide_imine: bool) -> List[Dict]:
    """Return acyclic, non-metal-bonded sp2 N that should be planar but is not
    constrained by a ring.

    Always includes nitro groups (kind="nitro").  When ``include_amide_imine``
    is True, also includes amide/carbamate/urea N (kind="amide") and 3-coord
    imine N (kind="imine") — geometry-only planarisation, no bond retarget.
    """
    groups = [dict(g, kind="nitro") for g in detect_nitro_groups(mol)]
    if not include_amide_imine:
        return groups
    try:
        from rdkit import Chem
    except Exception:
        return groups
    nitro_n = {g["n_idx"] for g in groups}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        if atom.GetIdx() in nitro_n:
            continue
        if atom.IsInRing():
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        if any(_is_metal_sym(nb.GetSymbol()) for nb in nbrs):
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP2:
            continue
        kind = None
        # amide / carbamate / urea: a neighbour C bears a =O / =S / =N.
        for nb in nbrs:
            if nb.GetSymbol() != "C":
                continue
            for nb2 in nb.GetNeighbors():
                if nb2.GetIdx() == atom.GetIdx():
                    continue
                bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE and \
                        nb2.GetSymbol() in ("O", "S", "N"):
                    kind = "amide"
                    break
            if kind:
                break
        # imine: N itself is double-bonded to a C.
        if kind is None:
            for b in atom.GetBonds():
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    other = b.GetOtherAtom(atom)
                    if other.GetSymbol() == "C":
                        kind = "imine"
                        break
        if kind is None:
            continue
        groups.append({"n_idx": atom.GetIdx(),
                       "nbr_idxs": [nb.GetIdx() for nb in nbrs],
                       "kind": kind})
    return groups


# ---------------------------------------------------------------------------
# Clash gate
# ---------------------------------------------------------------------------

def _introduces_new_clash(old_pts: np.ndarray, new_pts: np.ndarray,
                          syms: List[str], moved_idxs: List[int],
                          factor: float = 0.80) -> bool:
    """True if a moved atom gains a NEW non-bonded heavy contact
    (d < factor·Σr_vdw) that wasn't present in old_pts.  Bonded pairs and
    pairs among the moved set are exempt."""
    n = len(syms)
    moved = set(moved_idxs)
    for i in moved_idxs:
        ri_v = _VDW_RADII.get(syms[i], 1.55)
        ri_c = _COV_RADII.get(syms[i], 0.71)
        for j in range(n):
            if j == i or j in moved:
                continue
            sj = syms[j]
            rj_v = _VDW_RADII.get(sj, 1.70)
            rj_c = _COV_RADII.get(sj, 0.76)
            d_old = float(np.linalg.norm(old_pts[i] - old_pts[j]))
            if d_old < 1.15 * (ri_c + rj_c):
                continue  # bonded / near-bonded
            d_new = float(np.linalg.norm(new_pts[i] - new_pts[j]))
            thr = factor * (ri_v + rj_v)
            if d_new < thr <= d_old:
                return True
    return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def planarize_sp2_nitrogen(xyz: str, mol,
                           oop_threshold_A: float = 0.20,
                           no_target_A: float = 1.22,
                           ono_target_deg: float = 125.0,
                           include_amide_imine: bool = False,
                           ) -> Tuple[str, Dict]:
    """Flatten pyramidalised / desymmetrised acyclic sp2 N.

    For each nitro group (N + 2 O + 1 C) whose N is more than
    ``oop_threshold_A`` out of the C-O-O plane OR whose N-O lengths differ by
    more than 0.12 Å OR whose O-N-O deviates from 125° by more than 12°,
    rebuild a symmetric planar NO2:
        - keep C and the C->N bond direction & length;
        - place the two O atoms in the C-O-O plane at ``no_target_A`` from N,
          symmetric about the C->N axis, separated by ``ono_target_deg``.

    For amide/imine N (only when ``include_amide_imine``), project N and its
    out-of-plane component back into the plane of its three neighbours
    (geometry-only; no bond retarget).

    Per-group rollback on new clash or non-finite output (never-worse).

    Returns ``(new_xyz, report)``.  report keys: ``n_candidates``,
    ``n_fixed``, ``max_oop_before_A``, ``per_n``.
    """
    syms, pts, orig_lines = _parse_xyz(xyz)
    report: Dict = {"n_candidates": 0, "n_fixed": 0,
                    "max_oop_before_A": 0.0, "per_n": []}
    if pts.shape[0] == 0 or mol is None:
        return xyz, report
    if mol.GetNumAtoms() != len(syms):
        return xyz, report

    try:
        groups = detect_planar_sp2n_groups(mol, include_amide_imine)
    except Exception:
        return xyz, report
    report["n_candidates"] = len(groups)
    if not groups:
        return xyz, report

    new_pts = pts.copy()

    for g in groups:
        kind = g["kind"]
        n_idx = g["n_idx"]
        try:
            if kind == "nitro":
                o1, o2 = g["o_idxs"]
                c = g["c_idx"]
                p_n = new_pts[n_idx]
                p_o1 = new_pts[o1]
                p_o2 = new_pts[o2]
                p_c = new_pts[c]

                # Current metrics (trigger gate).
                oop = _oop_distance(p_n, p_o1, p_o2, p_c)
                no1 = float(np.linalg.norm(p_n - p_o1))
                no2 = float(np.linalg.norm(p_n - p_o2))
                ono = _angle_deg(p_n, p_o1, p_o2)
                if oop is None or ono is None:
                    continue
                report["max_oop_before_A"] = max(
                    report["max_oop_before_A"], float(oop))
                triggered = (oop > oop_threshold_A
                             or abs(no1 - no2) > 0.12
                             or abs(ono - 125.0) > 12.0)
                if not triggered:
                    report["per_n"].append({
                        "n_idx": n_idx, "kind": kind, "fixed": False,
                        "reason": "within_tolerance",
                        "oop_A": round(float(oop), 3),
                        "no_A": (round(no1, 3), round(no2, 3)),
                        "ono_deg": round(float(ono), 1)})
                    continue

                # Plane through {C, O1, O2}; if degenerate, use {C, O1, N}.
                nrm = np.cross(p_o1 - p_c, p_o2 - p_c)
                if float(np.linalg.norm(nrm)) < 1e-6:
                    nrm = np.cross(p_o1 - p_c, p_n - p_c)
                nn = float(np.linalg.norm(nrm))
                if nn < 1e-6:
                    continue
                normal = nrm / nn

                # C->N axis (in-plane); keep current C-N length & direction.
                cn = p_n - p_c
                cn_len = float(np.linalg.norm(cn))
                if cn_len < 0.5:
                    continue
                # Project the C->N direction onto the plane so N stays planar.
                cn_dir = cn - np.dot(cn, normal) * normal
                cn_dir_n = float(np.linalg.norm(cn_dir))
                if cn_dir_n < 1e-6:
                    continue
                cn_dir = cn_dir / cn_dir_n
                # N sits at C + cn_len along the in-plane C->N direction.
                new_n = p_c + cn_len * cn_dir

                # In-plane axis perpendicular to C->N for placing the two O.
                perp = np.cross(normal, cn_dir)
                perp_n = float(np.linalg.norm(perp))
                if perp_n < 1e-6:
                    continue
                perp = perp / perp_n

                # The N-O bonds point "away" from C (outward), each at half the
                # O-N-O angle off the (N->away-from-C) axis = +cn_dir from N.
                half = math.radians(ono_target_deg / 2.0)
                out_axis = cn_dir  # direction from C to N continued = outward
                # Decide perp sign so O1 keeps its current side (less motion).
                side1 = float(np.dot(p_o1 - new_n, perp))
                s1 = 1.0 if side1 >= 0 else -1.0
                dir_o1 = math.cos(half) * out_axis + math.sin(half) * s1 * perp
                dir_o2 = math.cos(half) * out_axis - math.sin(half) * s1 * perp
                new_o1 = new_n + no_target_A * dir_o1
                new_o2 = new_n + no_target_A * dir_o2

                cand = new_pts.copy()
                cand[n_idx] = new_n
                cand[o1] = new_o1
                cand[o2] = new_o2
                moved = [n_idx, o1, o2]

                if not np.all(np.isfinite(cand[moved])):
                    continue
                if _introduces_new_clash(new_pts, cand, syms, moved):
                    report["per_n"].append({
                        "n_idx": n_idx, "kind": kind, "fixed": False,
                        "reason": "rollback_clash"})
                    continue

                # Verify it actually improved planarity (never-worse).
                new_oop = _oop_distance(cand[n_idx], cand[o1], cand[o2],
                                        cand[c])
                if new_oop is not None and new_oop > float(oop) + 1e-6:
                    report["per_n"].append({
                        "n_idx": n_idx, "kind": kind, "fixed": False,
                        "reason": "no_improvement"})
                    continue

                new_pts = cand
                report["n_fixed"] += 1
                report["per_n"].append({
                    "n_idx": n_idx, "kind": kind, "fixed": True,
                    "oop_before_A": round(float(oop), 3),
                    "oop_after_A": round(float(new_oop), 3)
                    if new_oop is not None else None,
                    "no_before_A": (round(no1, 3), round(no2, 3)),
                    "ono_before_deg": round(float(ono), 1)})

            else:  # amide / imine — geometry-only plane projection
                nbr = g["nbr_idxs"]
                if len(nbr) != 3:
                    continue
                p_n = new_pts[n_idx]
                p_a, p_b, p_c = (new_pts[nbr[0]], new_pts[nbr[1]],
                                 new_pts[nbr[2]])
                oop = _oop_distance(p_n, p_a, p_b, p_c)
                if oop is None or oop <= oop_threshold_A:
                    continue
                nrm = np.cross(p_b - p_a, p_c - p_a)
                nn = float(np.linalg.norm(nrm))
                if nn < 1e-6:
                    continue
                normal = nrm / nn
                # Drop N's out-of-plane component (project onto the plane).
                disp = np.dot(p_n - p_a, normal)
                new_n = p_n - disp * normal
                cand = new_pts.copy()
                cand[n_idx] = new_n
                if not np.all(np.isfinite(cand[n_idx])):
                    continue
                if _introduces_new_clash(new_pts, cand, syms, [n_idx]):
                    report["per_n"].append({
                        "n_idx": n_idx, "kind": kind, "fixed": False,
                        "reason": "rollback_clash"})
                    continue
                new_pts = cand
                report["n_fixed"] += 1
                report["per_n"].append({
                    "n_idx": n_idx, "kind": kind, "fixed": True,
                    "oop_before_A": round(float(oop), 3)})
        except Exception:
            # Per-group fallback: leave this group untouched.
            continue

    if report["n_fixed"] == 0:
        return xyz, report
    if not np.all(np.isfinite(new_pts)):
        return xyz, report
    return _format_xyz(orig_lines, syms, new_pts), report
