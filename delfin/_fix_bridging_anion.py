"""Targeted bridging-anion (μ-X) M-X-M angle fixer.

Pattern: bimetallic complexes with bridging anionic donors X (Cl, Br, OH,
OR, NR2, SR, CR3) where the M-X-M angle deviates from chemistry-realistic
values, typically because ETKDG/UFF puts the anion on the wrong side of
the metal-metal axis or stretches a 4-ring [M-X-M-X'] to near-linearity.

This module is a sibling of ``_fix_wuxqak_sp3_c_linear`` and uses the same
architecture (BFS rigid-body rotation + per-violation revert).

Universal design (no SMILES/refcode):
  - μ-X detection is purely graph-driven: any non-metal atom bonded to at
    least 2 metals (via RDKit GetNeighbors / mol bond graph).  No element
    whitelist.
  - Target M-X-M angle is inferred from the local graph topology:
      • If X is part of a 4-membered metallacycle [M_a - X - M_b - X'] (i.e.
        the two metals are also bridged by a second μ-X'): target 90-100°
        (classic dimer geometry, e.g. Pt2(μ-Cl)2, Cu2(μ-OR)2).
      • Else, single-bridge [M - X - M]:
          - X is sp2/oxo (no H/R substituents apart from M's): target 145-170°
            (oxo bridge, [Fe2(μ-O)] type).
          - X is sp3 (has at least one non-metal substituent or H): target
            105-125° (μ-OH, μ-OR, μ-NR2, μ-CR3).
  - Geometry is purely positional (heavy-atom angle); RDKit is used only
    for the bond graph + hybridization hint, fallback to graph-only via
    covalent radii when RDKit hybridization is unavailable.

Operation:
    1. Detect μ-X bridges in the input mol (graph-only).
    2. For each, compute current M-X-M angle from xyz.
    3. If angle outside target window, attempt a rigid-body rotation that
       bends X (and its rigid substituent BFS fragment) toward target.
       Both metals stay fixed; X moves on a circle on the perpendicular
       bisector plane of the M-M axis.
    4. Per-violation rollback if:
        a) any intact M-D bond drifts > 0.05 Å, or
        b) a new heavy-atom clash appears.

Doctrine consistent with B3/B4/F19/F25/WUXQAK:
    - Universal: covalent radii + atomic-number range for metal detection.
    - Per-violation rollback; cumulative state defensively re-verified.
    - Post-ETKDG/UFF: bit-exact when ``DELFIN_FIX_BRIDGING_ANION=0`` (default).
    - Class-conditional via ``DELFIN_FIX_BRIDGING_ANION_CLASSES`` (allow-list).
    - Per-conformer; ``n_frames`` invariant.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Shared geometry helpers (kept local — module is self-contained).
# Tables match _fix_wuxqak_sp3_c_linear / _coord_angle_corrector.
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
    """Parse XYZ → (symbols, Nx3 coords, raw_lines for round-trip writer)."""
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
    """Round-trip XYZ writer (preserves header / comment lines)."""
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
    """Canonical metal-symbol set used by the converter."""
    try:
        from delfin.smiles_converter import _METAL_SET
        return set(_METAL_SET)
    except Exception:  # pragma: no cover - graceful fallback
        z_ranges = (set(range(21, 31)) | set(range(39, 49))
                    | set(range(57, 81)) | set(range(89, 104))
                    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51,
                       55, 56, 81, 82, 83})
        z_to_sym = {
            3: "Li", 4: "Be", 11: "Na", 12: "Mg", 13: "Al", 19: "K",
            20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn",
            26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga",
            37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
            43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
            49: "In", 50: "Sn", 51: "Sb", 55: "Cs", 56: "Ba", 57: "La",
            72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
            78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi",
        }
        return {z_to_sym[z] for z in z_ranges if z in z_to_sym}


def _bfs_fragment(mol, start_idx: int, blocked: Set[int]) -> Set[int]:
    """BFS over mol bonds from start_idx, never traversing blocked atoms."""
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
    """Rodrigues rotation matrix; axis must be unit-length."""
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(axis, axis)


# ---------------------------------------------------------------------------
# Detection + target-angle inference
# ---------------------------------------------------------------------------

def _classify_bridge_target_deg(
    mol, x_idx: int, metal_idxs: List[int],
) -> Tuple[float, float, float, str]:
    """Infer target M-X-M angle window for a μ-X bridge.

    Returns ``(target_deg, lo_ok, hi_ok, motif)``.  ``lo_ok``/``hi_ok``
    delimit the in-range angle window (no fix if current angle is inside).
    ``motif`` is a short label used in the report.

    Universal classification (no element/SMILES list):
      • 4-ring motif: X is in a 4-membered metallacycle with a second μ-X'
        bridge between the same two metals → ``target=92°`` (typical
        [M2(μ-X)2] dimer).
      • sp2/oxo motif: X has no non-metal neighbours apart from the two
        metals (i.e. it is truly μ-X with no H/C substituents) → linear-
        leaning, ``target=145°``.  This covers μ-O / μ-N=M' systems.
      • sp3/single motif (default): single-bridge with one or more non-
        metal substituents (μ-OH, μ-OR, μ-NR2, μ-CR3) → ``target=110°``.
    """
    # Default windows (chemistry-realistic ranges).
    DEFAULT_TARGET = 110.0
    DEFAULT_LO = 95.0
    DEFAULT_HI = 130.0

    # 4-ring motif: is x_idx in a 4-membered ring containing at least 2
    # metals?  RDKit's RingInfo may be empty for mols built via RWMol
    # without an explicit GetSymmSSSR call (e.g. unit-test fixtures), so
    # we fall back to a graph-only 4-cycle search via the bond graph.
    in_4ring_with_two_metals = False
    try:
        ring_info = mol.GetRingInfo()
        atom_rings = list(ring_info.AtomRings()) if ring_info is not None else []
    except Exception:
        atom_rings = []
    if not atom_rings:
        try:
            from rdkit import Chem as _Chem
            _Chem.GetSymmSSSR(mol)
            ring_info = mol.GetRingInfo()
            atom_rings = list(ring_info.AtomRings()) if ring_info is not None else []
        except Exception:
            atom_rings = []
    for ring in atom_rings:
        if len(ring) != 4:
            continue
        if x_idx not in ring:
            continue
        metals_in_ring = sum(
            1 for r in ring if r in metal_idxs
        )
        if metals_in_ring >= 2:
            in_4ring_with_two_metals = True
            break
    if not in_4ring_with_two_metals:
        # Graph-only fallback: does X share a *different* μ-X' partner
        # with both metal_idxs?  i.e. is there a second non-metal atom
        # Y, Y != X, such that Y is bonded to both metal_idxs[0] AND
        # metal_idxs[1]?  This is the defining feature of a [M2(μ-X)2]
        # 4-ring without depending on RDKit ring detection.
        try:
            atom = mol.GetAtomWithIdx(x_idx)
            atoms_bonded_to_metal_a = {
                nb.GetIdx() for nb in mol.GetAtomWithIdx(
                    metal_idxs[0]).GetNeighbors()
            }
            atoms_bonded_to_metal_b = {
                nb.GetIdx() for nb in mol.GetAtomWithIdx(
                    metal_idxs[1]).GetNeighbors()
            }
            shared = atoms_bonded_to_metal_a & atoms_bonded_to_metal_b
            # shared contains x_idx itself; we need a second non-metal partner.
            shared.discard(x_idx)
            metals_set = set(metal_idxs)
            for y in shared:
                if y in metals_set:
                    continue
                if mol.GetAtomWithIdx(y).GetSymbol() in _metal_set():
                    continue
                in_4ring_with_two_metals = True
                break
        except Exception:
            pass

    if in_4ring_with_two_metals:
        return 92.0, 80.0, 105.0, "4-ring"

    # Substituent count (non-metal heavy/H neighbours other than the two
    # bridging metals).  If zero → oxo/nitrido-style, sp2-leaning.
    try:
        atom = mol.GetAtomWithIdx(x_idx)
        non_metal_subs = [
            nb.GetIdx() for nb in atom.GetNeighbors()
            if nb.GetIdx() not in metal_idxs
        ]
    except Exception:
        non_metal_subs = []

    if not non_metal_subs:
        # μ-oxo / μ-nitrido style: bent but more open than sp3.
        return 145.0, 125.0, 175.0, "oxo"

    # sp3 single-bridge (μ-OH, μ-OR, μ-NR2, μ-CR3).
    return DEFAULT_TARGET, DEFAULT_LO, DEFAULT_HI, "sp3"


def detect_mu_x_bridges(
    xyz: str, mol, *, min_metals: int = 2,
) -> List[Dict]:
    """Detect μ-X bridging anions on the mol graph.

    Returns a list of dicts with keys ``x_idx``, ``metal_idxs``,
    ``x_symbol``, ``current_angle_deg`` (worst angle if >2 metals),
    ``target_deg``, ``lo_ok``, ``hi_ok``, ``motif``.

    A μ-X is any non-metal atom bonded to ``>=min_metals`` metals (graph,
    not distance — uses RDKit mol bonds).  Empty list when nothing
    detected (no metal pair, no shared anion, etc.).
    """
    if mol is None:
        return []

    syms, coords, _ = _parse_xyz(xyz)
    if coords.shape[0] == 0:
        return []
    n_atoms = mol.GetNumAtoms()
    if n_atoms != len(syms):
        # XYZ and mol disagree on atom count — refuse to act (defensive).
        return []

    metals = _metal_set()
    bridges: List[Dict] = []

    for x_idx in range(n_atoms):
        atom = mol.GetAtomWithIdx(x_idx)
        if atom.GetSymbol() in metals:
            continue
        # H is not chemistry-meaningful as a μ-X carrier here (hydride
        # bridges exist but are rare and have their own coord chemistry —
        # excluded to avoid false positives on stretched H atoms).
        if atom.GetSymbol() == "H":
            continue

        metal_nbrs = [
            nb.GetIdx() for nb in atom.GetNeighbors()
            if nb.GetSymbol() in metals
        ]
        if len(metal_nbrs) < min_metals:
            continue

        # Compute the worst (largest |Δ-from-target|) M-X-M angle across
        # all metal pairs at this X.  We need the worst pair for both
        # the classify call (target depends on motif) and the trigger
        # decision.
        target_deg, lo_ok, hi_ok, motif = _classify_bridge_target_deg(
            mol, x_idx, metal_nbrs,
        )

        worst_pair: Optional[Tuple[int, int]] = None
        worst_angle: Optional[float] = None
        worst_dev = -1.0
        for i_m in range(len(metal_nbrs)):
            for j_m in range(i_m + 1, len(metal_nbrs)):
                ma = metal_nbrs[i_m]
                mb = metal_nbrs[j_m]
                v1 = coords[ma] - coords[x_idx]
                v2 = coords[mb] - coords[x_idx]
                n1 = float(np.linalg.norm(v1))
                n2 = float(np.linalg.norm(v2))
                if n1 < 1e-9 or n2 < 1e-9:
                    continue
                cos_a = float(np.clip(np.dot(v1, v2) / (n1 * n2),
                                      -1.0, 1.0))
                angle = float(math.degrees(math.acos(cos_a)))
                # Deviation = how far outside the window [lo_ok, hi_ok]
                if angle < lo_ok:
                    dev = lo_ok - angle
                elif angle > hi_ok:
                    dev = angle - hi_ok
                else:
                    dev = 0.0
                if dev > worst_dev:
                    worst_dev = dev
                    worst_pair = (ma, mb)
                    worst_angle = angle

        if worst_pair is None or worst_angle is None:
            continue

        bridges.append({
            "x_idx": x_idx,
            "x_symbol": atom.GetSymbol(),
            "metal_idxs": metal_nbrs,
            "worst_pair": worst_pair,
            "current_angle_deg": worst_angle,
            "target_deg": target_deg,
            "lo_ok": lo_ok,
            "hi_ok": hi_ok,
            "motif": motif,
            "dev_deg": worst_dev,
        })
    return bridges


# ---------------------------------------------------------------------------
# Topology / clash safety
# ---------------------------------------------------------------------------

def _intact_md_snapshot(
    coords: np.ndarray, syms: List[str], metals: Set[str],
) -> Dict[Tuple[int, int], float]:
    """Snapshot every (metal, donor) bond inside the intact range
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


def _md_invariant_ok(
    pre: Dict[Tuple[int, int], float], new_coords: np.ndarray,
    tol: float = 0.05,
) -> bool:
    for (m_idx, d_idx), pre_d in pre.items():
        new_d = float(np.linalg.norm(new_coords[m_idx] - new_coords[d_idx]))
        if abs(new_d - pre_d) > tol:
            return False
    return True


def _new_clash_pair_appeared(
    pre_coords: np.ndarray, new_coords: np.ndarray, syms: List[str],
    moving: Set[int], frozen: Set[int],
    threshold_factor_heavy: float = 0.85,
    threshold_factor_h: float = 0.80,
) -> bool:
    """True if a new (moving, rest) clash pair appeared after the move."""
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
                if j in frozen:
                    # Frozen atoms (e.g. the bridging metals themselves)
                    # are part of the natural geometry — don't count them
                    # as "new" clashes since they are intrinsic.
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
    coords: np.ndarray, syms: List[str], mol,
    bridge: Dict, metals: Set[str],
    md_snapshot: Dict[Tuple[int, int], float],
) -> Tuple[bool, float, Set[int]]:
    """Try a rigid translation of X (and its non-metal BFS substituent
    fragment) so that the worst M-X-M angle approaches ``target_deg``.

    Both bridging metals stay fixed.  X is moved onto the perpendicular
    bisector plane of the M_a-M_b axis at a distance that makes M-X-M
    equal target_deg.  The M-X bonds are allowed to change length (the
    fix is FOR the angle — bond lengths follow); other bonds keep their
    pre-fix length via the M-D snapshot filter.

    Returns ``(applied, correction_deg, moved_fragment)``.  Mutates
    ``coords`` in place on success.  ``moved_fragment`` is the set of
    atom indices that were translated (empty when not applied).
    """
    x_idx = bridge["x_idx"]
    ma, mb = bridge["worst_pair"]
    target_deg = bridge["target_deg"]

    p_ma = coords[ma]
    p_mb = coords[mb]
    p_x = coords[x_idx]
    midpoint = 0.5 * (p_ma + p_mb)
    mm_axis = p_mb - p_ma
    mm_norm = float(np.linalg.norm(mm_axis))
    if mm_norm < 1e-6:
        return False, 0.0, set()
    u_mm = mm_axis / mm_norm

    # Current vector from midpoint to X.  Decompose into component along
    # M-M axis (parallel, must be preserved to keep |M-X| invariant) and
    # perpendicular component (which is what we rotate).
    v_mid_x = p_x - midpoint
    v_par = float(np.dot(v_mid_x, u_mm)) * u_mm  # along axis
    v_perp = v_mid_x - v_par
    perp_norm = float(np.linalg.norm(v_perp))
    if perp_norm < 1e-6:
        # X lies on M-M line (collinear) — can't define rotation plane.
        # Pick an arbitrary perpendicular and lift X off the axis a bit.
        ref = (np.array([1.0, 0.0, 0.0])
               if abs(u_mm[0]) < 0.9 else np.array([0.0, 1.0, 0.0]))
        v_perp = ref - float(np.dot(ref, u_mm)) * u_mm
        perp_norm = float(np.linalg.norm(v_perp))
        if perp_norm < 1e-6:
            return False, 0.0, set()
        v_perp = v_perp / perp_norm  # unit vector for direction only
        perp_norm = 0.01  # tiny offset to bootstrap

    # Compute the perp_distance and parallel offset that yield the target
    # M-X-M angle, while preserving the geometric mean of the two M-X
    # distances.  For symmetric placement (X on perpendicular bisector,
    # v_par == 0), the M-X-M angle is 2 * atan(d_mm/2 / perp_dist).  We
    # solve for the new perp_dist that gives target_deg.
    d_mm = mm_norm
    half_d_mm = 0.5 * d_mm

    # Current perp_dist (signed via the unit perp direction).
    u_perp = v_perp / perp_norm

    # Target perp distance: perp_target = (d_mm/2) / tan(target/2)
    target_rad = math.radians(target_deg)
    tan_half = math.tan(0.5 * target_rad)
    if tan_half < 1e-6:
        return False, 0.0, set()
    perp_target = half_d_mm / tan_half

    # The current preserves the parallel offset v_par (otherwise we would
    # alter |M_a-X| and |M_b-X| unequally).  We just move X along the
    # u_perp direction.
    new_x = midpoint + v_par + perp_target * u_perp

    # Sanity: should not be NaN, should not be coincident with metal.
    if not np.all(np.isfinite(new_x)):
        return False, 0.0, set()
    if (float(np.linalg.norm(new_x - p_ma)) < 0.3
            or float(np.linalg.norm(new_x - p_mb)) < 0.3):
        return False, 0.0, set()

    # Determine BFS fragment to move along with X: X itself + any non-
    # metal substituent BFS, with both bridging metals blocked.  Chelate
    # guard: if BFS reaches another metal not in {ma, mb}, the fragment
    # is part of a chelate ring → skip.
    blocked: Set[int] = set(bridge["metal_idxs"])
    fragment = _bfs_fragment(mol, x_idx, blocked)
    if not fragment:
        return False, 0.0, set()
    # Reject pathological fragments that contain almost everything else.
    if len(fragment) >= len(syms) - len(bridge["metal_idxs"]):
        return False, 0.0, set()
    for f_idx in fragment:
        if syms[f_idx] in metals and f_idx not in bridge["metal_idxs"]:
            return False, 0.0, set()  # chelate spanning a third metal — abort

    # Translation applied to X is also applied to the rigid fragment.
    delta = new_x - p_x

    candidate_coords = coords.copy()
    fragment_list = sorted(fragment)
    for i in fragment_list:
        candidate_coords[i] = coords[i] + delta

    # Topology gate: M-D invariant on every pre-existing bond EXCEPT the
    # two M-X bonds we are actively correcting (those distances are
    # expected to change — we are moving X precisely to fix the angle).
    snapshot_filtered = {
        k: v for k, v in md_snapshot.items()
        if not (k[1] in fragment)  # any bond into the moving fragment is exempt
    }
    if not _md_invariant_ok(snapshot_filtered, candidate_coords):
        return False, 0.0, set()

    # Clash gate.  Freeze both bridging metals so the intrinsic M-X
    # contact is not flagged.
    frozen: Set[int] = set(bridge["metal_idxs"])
    if _new_clash_pair_appeared(
        coords, candidate_coords, syms, fragment, frozen,
    ):
        return False, 0.0, set()

    # Verify the angle actually improved.
    new_v1 = candidate_coords[ma] - candidate_coords[x_idx]
    new_v2 = candidate_coords[mb] - candidate_coords[x_idx]
    n1 = float(np.linalg.norm(new_v1))
    n2 = float(np.linalg.norm(new_v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return False, 0.0, set()

    # Chemistry-realistic M-X bond-length gate: each post-move |M-X| must
    # lie within [0.75, 1.40] × (r_M + r_X) covalent-radius sum so we
    # never produce a Pt-Cl of 4 Å or a fused 1 Å contact.
    r_x = _COV_RADII.get(syms[x_idx], 1.5)
    r_ma = _COV_RADII.get(syms[ma], 1.5)
    r_mb = _COV_RADII.get(syms[mb], 1.5)
    target_a_lo = 0.75 * (r_ma + r_x)
    target_a_hi = 1.40 * (r_ma + r_x)
    target_b_lo = 0.75 * (r_mb + r_x)
    target_b_hi = 1.40 * (r_mb + r_x)
    if not (target_a_lo <= n1 <= target_a_hi):
        return False, 0.0, set()
    if not (target_b_lo <= n2 <= target_b_hi):
        return False, 0.0, set()

    cos_new = float(np.clip(np.dot(new_v1, new_v2) / (n1 * n2), -1.0, 1.0))
    new_angle = float(math.degrees(math.acos(cos_new)))

    pre_dev = abs(bridge["current_angle_deg"] - target_deg)
    new_dev = abs(new_angle - target_deg)
    if new_dev >= pre_dev - 1.0:  # require >1° improvement; else no-op
        return False, 0.0, set()

    coords[:] = candidate_coords
    return True, float(abs(new_angle - bridge["current_angle_deg"])), set(fragment)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fix_bridging_anion_angles(
    xyz: str, mol,
) -> Tuple[str, Dict]:
    """Detect and bend μ-X bridges so that M-X-M angles fall into a
    chemistry-realistic window for their motif.

    Parameters
    ----------
    xyz : str
        DELFIN XYZ string.
    mol : rdkit.Chem.Mol
        Molecule with the same atom order as ``xyz``.

    Returns
    -------
    new_xyz : str
        Corrected XYZ (equals input if nothing fired or nothing was
        successfully applied).
    report : dict
        Keys: ``n_bridges`` (detected), ``n_violations`` (current angle
        outside window), ``n_fixed`` (successful bends),
        ``max_correction_deg``, ``topology_preserved``, ``per_bridge``.
    """
    report: Dict = {
        "n_bridges": 0,
        "n_violations": 0,
        "n_fixed": 0,
        "max_correction_deg": 0.0,
        "topology_preserved": True,
        "per_bridge": [],
    }

    bridges = detect_mu_x_bridges(xyz, mol)
    report["n_bridges"] = len(bridges)
    if not bridges:
        return xyz, report

    violations = [b for b in bridges if b["dev_deg"] > 0.0]
    report["n_violations"] = len(violations)
    if not violations:
        return xyz, report

    syms, coords, orig_lines = _parse_xyz(xyz)
    if coords.shape[0] == 0 or len(syms) != mol.GetNumAtoms():
        return xyz, report

    metals = _metal_set()
    md_snapshot = _intact_md_snapshot(coords, syms, metals)

    # Sort worst-deviation first.
    violations.sort(key=lambda v: -v["dev_deg"])

    coords = coords.copy()
    fixed = 0
    max_corr = 0.0
    moved_atoms: Set[int] = set()
    for v in violations:
        applied, corr, frag_moved = _bend_one(
            coords, syms, mol, v, metals, md_snapshot,
        )
        report["per_bridge"].append({
            "x_idx": v["x_idx"],
            "x_symbol": v["x_symbol"],
            "metals": list(v["worst_pair"]),
            "motif": v["motif"],
            "pre_angle_deg": v["current_angle_deg"],
            "target_deg": v["target_deg"],
            "applied": bool(applied),
            "correction_deg": corr,
        })
        if applied:
            fixed += 1
            if corr > max_corr:
                max_corr = corr
            moved_atoms |= frag_moved

    report["n_fixed"] = fixed
    report["max_correction_deg"] = max_corr
    # Cumulative M-D invariant sanity-check.  Bonds whose donor atom was
    # part of a moved fragment are intentionally exempt — those distances
    # are the ones we corrected.
    cum_snapshot = {
        k: v for k, v in md_snapshot.items() if k[1] not in moved_atoms
    }
    report["topology_preserved"] = _md_invariant_ok(cum_snapshot, coords)

    if fixed == 0 or not report["topology_preserved"]:
        # If overall topology slipped (shouldn't, given per-bend gates),
        # bail out conservatively — return original XYZ.
        return xyz, report

    new_xyz = _format_xyz(orig_lines, syms, coords)
    return new_xyz, report
