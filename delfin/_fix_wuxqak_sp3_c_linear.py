"""Targeted WUXQAK / sp3-C linear-collapse fixer.

Pattern: M-CH2-R or M-CHR-R' complexes where an sp3 C bonded to a metal
collapses to a linear M-C-X angle (~180 degrees) instead of the expected
tetrahedral 109.5 degrees.  Memory `feedback_sp3_c_donor_linear.md` ties this
to WUXQAK and confirms HEAD has 10.6% files affected vs 5.2% in e6761e4.

This module exposes a small, surgical fixer that:
    1. Detects all sp3 C atoms that are (a) bonded to a metal and (b) carry
       at least one substituent X with M-C-X angle > ``angle_threshold_deg``.
    2. Rotates the rigid BFS sub-fragment of X around C (M and C fixed) by
       just enough rotation to bring M-C-X closer to ``target_angle_deg``.
    3. Reverts any rotation that introduces a new inter-fragment clash.

Design constraints (per mission spec):
    - Only the X sub-fragment moves; M and C are immovable pivots.
    - Rigid-body rotation preserves all bond lengths and internal angles
      inside the X side.
    - Conservative trigger: > 150 degrees default (only the truly-linear
      collapses, not mild distortions which UFF/B3 handle).
    - Topology safety: per-violation rollback if a new clash appears or
      any pre-existing intact M-D bond changes by more than 0.05 A.
    - No env-flag gate inside; caller decides when to apply.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Shared geometry helpers (kept local — module is self-contained).
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
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
    "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00, "Fe": 2.00,
    "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87, "Ge": 2.11,
    "As": 1.85, "Se": 1.90, "Y": 2.20, "Zr": 2.10, "Nb": 2.10, "Mo": 2.10,
    "Tc": 2.10, "Ru": 2.10, "Rh": 2.10, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "La": 2.30, "Hf": 2.10,
    "Ta": 2.10, "W": 2.10, "Re": 2.10, "Os": 2.10, "Ir": 2.10, "Pt": 1.75,
    "Au": 1.66, "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
}

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


def _metal_set() -> Set[str]:
    """Return the canonical metal-symbol set used by the converter."""
    try:
        from delfin.smiles_converter import _METAL_SET
        return set(_METAL_SET)
    except Exception:  # pragma: no cover - graceful fallback
        # Universal atomic-number fallback identical to _coord_angle_corrector.
        z_ranges = (set(range(21, 31)) | set(range(39, 49))
                    | set(range(57, 81)) | set(range(89, 104))
                    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51,
                       55, 56, 81, 82, 83})
        z_to_sym = {
            1: "H", 3: "Li", 4: "Be", 6: "C", 7: "N", 8: "O", 11: "Na",
            12: "Mg", 13: "Al", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
            23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni",
            29: "Cu", 30: "Zn", 31: "Ga", 37: "Rb", 38: "Sr", 39: "Y",
            40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh",
            46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb",
            55: "Cs", 56: "Ba", 57: "La", 72: "Hf", 73: "Ta", 74: "W",
            75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
            81: "Tl", 82: "Pb", 83: "Bi",
        }
        return {z_to_sym[z] for z in z_ranges if z in z_to_sym}


def _bfs_fragment(mol, start_idx: int, blocked: Set[int]) -> Set[int]:
    """BFS over mol bonds from ``start_idx``, never traversing through any
    atom in ``blocked``.  Returns the set of reachable atom indices
    (including ``start_idx`` itself when not blocked)."""
    if start_idx in blocked:
        return set()
    visited: Set[int] = {start_idx}
    queue: List[int] = [start_idx]
    while queue:
        cur = queue.pop()
        try:
            atom = mol.GetAtomWithIdx(cur)
        except Exception:
            continue
        for nb in atom.GetNeighbors():
            j = nb.GetIdx()
            if j in blocked or j in visited:
                continue
            visited.add(j)
            queue.append(j)
    return visited


def _rodrigues(axis: np.ndarray, theta_rad: float) -> np.ndarray:
    """Rodrigues rotation matrix; ``axis`` must be unit-length."""
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(axis, axis)


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def detect_wuxqak_linear_sp3_c(
    xyz: str, mol, threshold_deg: float = 150.0,
) -> List[Dict]:
    """Detect M-C-X violations on sp3 C atoms bonded to a metal.

    Returns a list of dicts with keys: ``c_idx``, ``metal_idx``, ``x_idx``,
    ``angle_deg``, ``x_symbol``.  Empty list when nothing offends.
    """
    syms, coords, _ = _parse_xyz(xyz)
    if coords.shape[0] == 0:
        return []
    metals = _metal_set()
    n_atoms = mol.GetNumAtoms()
    if n_atoms != len(syms):
        # XYZ and mol disagree on atom count — refuse to act.
        return []

    try:
        from rdkit import Chem
        SP3 = Chem.HybridizationType.SP3
    except Exception:  # pragma: no cover - rdkit always present in DELFIN
        return []

    violations: List[Dict] = []
    for c_idx in range(n_atoms):
        atom = mol.GetAtomWithIdx(c_idx)
        if atom.GetSymbol() != "C":
            continue
        try:
            if atom.GetHybridization() != SP3:
                continue
        except Exception:
            continue

        # First metal neighbour (M-CH... has typically a single metal bond).
        metal_idx: Optional[int] = None
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() in metals:
                metal_idx = nbr.GetIdx()
                break
        if metal_idx is None:
            continue

        v_cm = coords[metal_idx] - coords[c_idx]
        n_cm = float(np.linalg.norm(v_cm))
        if n_cm < 1e-9:
            continue
        u_cm = v_cm / n_cm

        for x_nbr in atom.GetNeighbors():
            x_idx = x_nbr.GetIdx()
            if x_idx == metal_idx:
                continue
            v_cx = coords[x_idx] - coords[c_idx]
            n_cx = float(np.linalg.norm(v_cx))
            if n_cx < 1e-9:
                continue
            u_cx = v_cx / n_cx
            cos_a = float(np.clip(np.dot(u_cm, u_cx), -1.0, 1.0))
            angle_deg = float(np.degrees(math.acos(cos_a)))

            if angle_deg > threshold_deg:
                violations.append({
                    "c_idx": c_idx,
                    "metal_idx": metal_idx,
                    "x_idx": x_idx,
                    "angle_deg": angle_deg,
                    "x_symbol": x_nbr.GetSymbol(),
                })
    return violations


# ---------------------------------------------------------------------------
# Topology / clash safety
# ---------------------------------------------------------------------------

def _intact_md_snapshot(coords: np.ndarray, syms: List[str],
                        metals: Set[str]) -> Dict[Tuple[int, int], float]:
    """Snapshot every (metal, donor) bond currently inside the intact range
    (d <= 1.30 * (r_M + r_D))."""
    metal_idxs = [i for i, s in enumerate(syms) if s in metals]
    snap: Dict[Tuple[int, int], float] = {}
    for m_idx in metal_idxs:
        rm = _COV_RADII.get(syms[m_idx], 1.5)
        for j in range(len(syms)):
            if j == m_idx:
                continue
            sj = syms[j]
            if sj == "H" or sj in metals:
                continue
            rj = _COV_RADII.get(sj, 1.5)
            d = float(np.linalg.norm(coords[m_idx] - coords[j]))
            if d <= 1.30 * (rm + rj):
                snap[(m_idx, j)] = d
    return snap


def _md_invariant_ok(pre: Dict[Tuple[int, int], float],
                     new_coords: np.ndarray, tol: float = 0.05) -> bool:
    for (m_idx, d_idx), pre_d in pre.items():
        new_d = float(np.linalg.norm(new_coords[m_idx] - new_coords[d_idx]))
        if abs(new_d - pre_d) > tol:
            return False
    return True


def _new_clash_pair_appeared(
    pre_coords: np.ndarray,
    new_coords: np.ndarray,
    syms: List[str],
    moving: Set[int],
    pivot_idx: int,
    metal_idx: int,
    threshold_factor_heavy: float = 0.85,
    threshold_factor_h: float = 0.80,
) -> bool:
    """Return True if rotating ``moving`` introduced a new (moving, rest)
    clash pair that did not exist in ``pre_coords``."""
    n = len(syms)

    def _clash_set(coords: np.ndarray) -> Set[Tuple[int, int]]:
        pairs: Set[Tuple[int, int]] = set()
        for i in moving:
            si = syms[i]
            ri_v = _VDW_RADII.get(si, 1.7 if si != "H" else 1.2)
            ri_c = _COV_RADII.get(si, 1.5 if si != "H" else 0.31)
            for j in range(n):
                if j == i or j in moving:
                    continue
                if j == pivot_idx or j == metal_idx:
                    continue
                sj = syms[j]
                rj_v = _VDW_RADII.get(sj, 1.7 if sj != "H" else 1.2)
                rj_c = _COV_RADII.get(sj, 1.5 if sj != "H" else 0.31)
                d = float(np.linalg.norm(coords[i] - coords[j]))
                if d < (ri_c + rj_c) * 1.10:
                    continue  # bonded
                thr = (threshold_factor_h
                       if (si == "H" or sj == "H")
                       else threshold_factor_heavy)
                if d < thr * (ri_v + rj_v):
                    pairs.add((min(i, j), max(i, j)))
        return pairs

    pre_pairs = _clash_set(pre_coords)
    new_pairs = _clash_set(new_coords)
    return bool(new_pairs - pre_pairs)


# ---------------------------------------------------------------------------
# Single-violation bend
# ---------------------------------------------------------------------------

def _bend_one(
    coords: np.ndarray,
    syms: List[str],
    mol,
    c_idx: int,
    metal_idx: int,
    x_idx: int,
    target_deg: float,
    metals: Set[str],
    md_snapshot: Dict[Tuple[int, int], float],
) -> Tuple[bool, float]:
    """Try a rigid-body rotation of X's BFS fragment about C so that M-C-X
    approaches ``target_deg``.  Returns ``(applied, correction_deg)``.
    Mutates ``coords`` in place on success."""
    pivot = coords[c_idx]
    v_cm = coords[metal_idx] - pivot
    v_cx = coords[x_idx] - pivot
    n_cm = float(np.linalg.norm(v_cm))
    n_cx = float(np.linalg.norm(v_cx))
    if n_cm < 1e-9 or n_cx < 1e-9:
        return False, 0.0

    u_cm = v_cm / n_cm
    u_cx = v_cx / n_cx

    cos_now = float(np.clip(np.dot(u_cm, u_cx), -1.0, 1.0))
    current_deg = float(np.degrees(math.acos(cos_now)))
    delta = target_deg - current_deg  # negative — we are too linear

    # Rotation axis perpendicular to the (M, C, X) plane.
    axis = np.cross(v_cm, v_cx)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-9:
        # Collinear: pick an arbitrary perpendicular to u_cm.  This is the
        # near-180 degree case — without a defined plane, choose any direction.
        ref = (np.array([1.0, 0.0, 0.0])
               if abs(u_cm[0]) < 0.9 else np.array([0.0, 1.0, 0.0]))
        axis = np.cross(u_cm, ref)
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm < 1e-9:
            return False, 0.0
    axis = axis / axis_norm

    # BFS fragment from X with both C and the metal blocked.  Chelate guard:
    # if the BFS reaches another metal, the fragment is part of a chelate
    # ring; rotating it would tear that ring — skip.
    blocked: Set[int] = {c_idx, metal_idx}
    fragment = _bfs_fragment(mol, x_idx, blocked)
    if not fragment:
        return False, 0.0
    # Reject pathological fragments that contain almost everything else.
    if len(fragment) >= len(syms) - 2:
        return False, 0.0
    for f_idx in fragment:
        if f_idx == metal_idx:
            return False, 0.0
        if syms[f_idx] in metals:
            return False, 0.0  # chelate — abort
    fragment_list = sorted(fragment)

    # Candidate corrections — favour full delta, fall back to fractions when
    # the gates reject the full rotation.  Both signs are tried since the
    # geometric axis can point either way relative to the desired bend.
    candidate_thetas: List[float] = []
    for sign in (+1, -1):
        for fac in (1.0, 0.75, 0.5, 0.25):
            candidate_thetas.append(sign * delta * fac)

    best_applied: Optional[np.ndarray] = None
    best_score = abs(current_deg - target_deg)  # must strictly improve
    best_correction = 0.0

    for theta_deg in candidate_thetas:
        if abs(theta_deg) < 1.0:
            continue  # noise-level; skip
        R = _rodrigues(axis, math.radians(theta_deg))
        new_coords = coords.copy()
        for i in fragment_list:
            new_coords[i] = pivot + R @ (coords[i] - pivot)

        new_v_cx = new_coords[x_idx] - pivot
        new_n_cx = float(np.linalg.norm(new_v_cx))
        if new_n_cx < 1e-9:
            continue
        new_cos = float(np.clip(np.dot(u_cm, new_v_cx / new_n_cx),
                                -1.0, 1.0))
        new_angle = float(np.degrees(math.acos(new_cos)))
        score = abs(new_angle - target_deg)
        if score >= best_score:
            continue  # not an improvement

        # Topology gate: M-D invariant on every pre-existing bond.
        if not _md_invariant_ok(md_snapshot, new_coords):
            continue

        # Clash identity gate.
        if _new_clash_pair_appeared(
            coords, new_coords, syms, fragment, c_idx, metal_idx,
        ):
            continue

        best_applied = new_coords
        best_score = score
        best_correction = abs(theta_deg)

    if best_applied is None:
        return False, 0.0

    coords[:] = best_applied
    return True, best_correction


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fix_wuxqak_sp3_c_linear(
    xyz: str,
    mol,
    angle_threshold_deg: float = 150.0,
    target_angle_deg: float = 109.5,
) -> Tuple[str, Dict]:
    """Detect and bend M-C-X linear collapses on sp3 C donors.

    Parameters
    ----------
    xyz : str
        DELFIN-format XYZ string.
    mol : rdkit Mol
        Molecule supplying bond connectivity and hybridization.
    angle_threshold_deg : float, default 150.0
        Only M-C-X angles strictly greater than this trigger a fix.
    target_angle_deg : float, default 109.5
        Target M-C-X angle after the bend.

    Returns
    -------
    new_xyz : str
        Corrected XYZ.  Equals input if nothing was applied.
    report : dict
        Keys ``n_violations``, ``n_fixed``, ``max_correction_deg``,
        ``topology_preserved``.
    """
    report: Dict = {
        "n_violations": 0,
        "n_fixed": 0,
        "max_correction_deg": 0.0,
        "topology_preserved": True,
    }

    violations = detect_wuxqak_linear_sp3_c(xyz, mol,
                                            threshold_deg=angle_threshold_deg)
    report["n_violations"] = len(violations)
    if not violations:
        return xyz, report

    syms, coords, orig_lines = _parse_xyz(xyz)
    if coords.shape[0] == 0 or len(syms) != mol.GetNumAtoms():
        return xyz, report

    metals = _metal_set()
    md_snapshot = _intact_md_snapshot(coords, syms, metals)

    # Sort worst-deviation first so the most catastrophic case gets the first
    # rotation budget; later passes re-detect because earlier moves may have
    # changed downstream angles.
    violations.sort(key=lambda v: -v["angle_deg"])

    coords = coords.copy()
    fixed = 0
    max_corr = 0.0
    for v in violations:
        applied, corr = _bend_one(
            coords, syms, mol,
            c_idx=v["c_idx"],
            metal_idx=v["metal_idx"],
            x_idx=v["x_idx"],
            target_deg=target_angle_deg,
            metals=metals,
            md_snapshot=md_snapshot,
        )
        if applied:
            fixed += 1
            if corr > max_corr:
                max_corr = corr

    report["n_fixed"] = fixed
    report["max_correction_deg"] = max_corr
    # Topology preserved iff every snapshotted M-D bond is still within
    # tolerance — _bend_one rolls back individually, but we re-verify the
    # cumulative state defensively.
    report["topology_preserved"] = _md_invariant_ok(md_snapshot, coords)

    if fixed == 0:
        return xyz, report

    new_xyz = _format_xyz(orig_lines, syms, coords)
    return new_xyz, report


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _self_test() -> None:  # pragma: no cover - exercised manually
    """Synthetic Pt-CH3 with one M-C-H at 180 degrees must bend to ~109.5."""
    print("[_fix_wuxqak_sp3_c_linear self-test]")

    try:
        from rdkit import Chem
        SP3 = Chem.HybridizationType.SP3
    except Exception:
        class _HT:
            SP3 = "SP3"
        SP3 = _HT.SP3  # type: ignore[assignment]

    class _MockAtom:
        def __init__(self, idx, sym, hyb):
            self.idx = idx
            self.sym = sym
            self.hyb = hyb
            self._nbrs: List["_MockAtom"] = []

        def GetIdx(self):
            return self.idx

        def GetSymbol(self):
            return self.sym

        def GetHybridization(self):
            return self.hyb

        def GetNeighbors(self):
            return list(self._nbrs)

    class _MockMol:
        def __init__(self, atoms):
            self._atoms = atoms

        def GetNumAtoms(self):
            return len(self._atoms)

        def GetAtomWithIdx(self, i):
            return self._atoms[i]

        def GetAtoms(self):
            return list(self._atoms)

    # Build Pt(0) - C(1) - H(2), H(3), H(4).  C is sp3.
    atoms = [
        _MockAtom(0, "Pt", SP3),
        _MockAtom(1, "C", SP3),
        _MockAtom(2, "H", SP3),
        _MockAtom(3, "H", SP3),
        _MockAtom(4, "H", SP3),
    ]
    atoms[0]._nbrs = [atoms[1]]
    atoms[1]._nbrs = [atoms[0], atoms[2], atoms[3], atoms[4]]
    atoms[2]._nbrs = [atoms[1]]
    atoms[3]._nbrs = [atoms[1]]
    atoms[4]._nbrs = [atoms[1]]
    mol = _MockMol(atoms)

    # Scenario 1 — H at 2 is collinear (180 deg) with Pt-C.
    coords = np.array([
        [0.00, 0.00, 0.00],   # Pt
        [2.00, 0.00, 0.00],   # C
        [3.10, 0.00, 0.00],   # H, M-C-H = 180 (BAD)
        [1.65, 0.95, 0.00],   # H tetrahedral
        [1.65, -0.45, 0.82],  # H tetrahedral
    ], dtype=float)
    xyz_in = "5\nPt-CH3 linear-H\n" + "\n".join(
        f"{s:4s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}"
        for s, c in zip(["Pt", "C", "H", "H", "H"], coords)
    ) + "\n"

    viols = detect_wuxqak_linear_sp3_c(xyz_in, mol, threshold_deg=150.0)
    print(f"  scenario-1 detected: {len(viols)} violation(s)")
    assert len(viols) == 1, viols
    assert viols[0]["c_idx"] == 1 and viols[0]["metal_idx"] == 0
    assert abs(viols[0]["angle_deg"] - 180.0) < 1e-3

    new_xyz, rep = fix_wuxqak_sp3_c_linear(xyz_in, mol)
    print(f"  scenario-1 report : {rep}")
    assert rep["n_violations"] == 1
    assert rep["n_fixed"] == 1
    assert rep["topology_preserved"] is True

    # Recheck the angle in the corrected XYZ.
    _, new_coords, _ = _parse_xyz(new_xyz)
    v_cm = new_coords[0] - new_coords[1]
    v_cx = new_coords[2] - new_coords[1]
    cos_a = float(np.dot(v_cm, v_cx)
                  / (np.linalg.norm(v_cm) * np.linalg.norm(v_cx)))
    new_ang = float(np.degrees(math.acos(max(-1.0, min(1.0, cos_a)))))
    print(f"  scenario-1 new M-C-H: {new_ang:.2f} deg")
    assert abs(new_ang - 109.5) < 5.0, new_ang

    # M-C bond length unchanged.
    d_mc = float(np.linalg.norm(new_coords[0] - new_coords[1]))
    print(f"  scenario-1 M-C bond : {d_mc:.4f} A (was 2.0000)")
    assert abs(d_mc - 2.0) < 1e-3

    # Scenario 2 — already-correct tetrahedral CH3 (noop expected).
    coords2 = np.array([
        [0.00, 0.00, 0.00],   # Pt
        [2.00, 0.00, 0.00],   # C
        [2.36, 1.02, 0.00],   # H ~ 109 deg
        [2.36, -0.51, 0.88],
        [2.36, -0.51, -0.88],
    ], dtype=float)
    xyz_in2 = "5\nPt-CH3 tetrahedral\n" + "\n".join(
        f"{s:4s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}"
        for s, c in zip(["Pt", "C", "H", "H", "H"], coords2)
    ) + "\n"
    viols2 = detect_wuxqak_linear_sp3_c(xyz_in2, mol, threshold_deg=150.0)
    print(f"  scenario-2 detected: {len(viols2)} violation(s)")
    assert len(viols2) == 0
    new_xyz2, rep2 = fix_wuxqak_sp3_c_linear(xyz_in2, mol)
    print(f"  scenario-2 report : {rep2}")
    assert rep2["n_violations"] == 0
    assert rep2["n_fixed"] == 0
    assert new_xyz2 == xyz_in2  # bit-exact passthrough

    # Scenario 3 — two H atoms above threshold simultaneously.
    coords3 = np.array([
        [0.00, 0.00, 0.00],   # Pt
        [2.00, 0.00, 0.00],   # C
        [3.10, 0.05, 0.00],   # H ~ 178 deg
        [3.05, -0.05, 0.30],  # H ~ 175 deg
        [1.65, -0.50, -0.82], # H tetrahedral
    ], dtype=float)
    xyz_in3 = "5\nPt-CH3 dual-linear\n" + "\n".join(
        f"{s:4s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}"
        for s, c in zip(["Pt", "C", "H", "H", "H"], coords3)
    ) + "\n"
    viols3 = detect_wuxqak_linear_sp3_c(xyz_in3, mol, threshold_deg=150.0)
    print(f"  scenario-3 detected: {len(viols3)} violation(s)")
    assert len(viols3) == 2
    _, rep3 = fix_wuxqak_sp3_c_linear(xyz_in3, mol)
    print(f"  scenario-3 report : {rep3}")
    assert rep3["n_violations"] == 2
    assert rep3["n_fixed"] >= 1
    assert rep3["topology_preserved"] is True

    print("  result             : PASS")


if __name__ == "__main__":
    _self_test()
