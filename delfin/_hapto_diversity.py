"""Iter-3 Hapto Diversity Topology-Aware (HD-TA).

Restores Hapto-Diversity volume that was killed by the Iter-1 σ-guard
(``DELFIN_HAPTO_NO_OB_WHEN_SIGMA=1``) without re-introducing the
Cp/η-detachment regressions that motivated the guard.

Mechanism: rather than calling Open Babel's WeightedRotorSearch (which
treats M-C(η) bonds as free torsions and detaches Cp/arene rings), this
module performs **rigid-body rotations of σ-fragments only** around their
M-D axes, leaving every hapto-atom and every M-D bond length bit-exact
identical to the seed conformer.

Public entrypoint:
    apply_hapto_diversity_topology_aware(
        smiles, seed_xyz, mol_template,
        max_frames=60, results_seen=None,
        verify_topology=fn, ...,
    ) -> List[Tuple[str, str]]
"""

from __future__ import annotations

import math
import os
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

try:
    import numpy as np  # noqa: F401
except ImportError:  # numpy is a hard dep of DELFIN; defensive only
    np = None  # type: ignore


# Tolerances — only loosened by env-flag overrides
_HAPTO_TOL_DEFAULT = 0.01     # Å, hapto-atom drift
_MD_TOL_DEFAULT = 0.05        # Å, M-D bond length drift
_RING_PLANARITY_TOL = 0.05    # Å, max distance to best-fit plane


def _env_truthy(name: str, default: str) -> bool:
    raw = os.environ.get(name, default).strip().lower()
    return raw not in {"0", "false", "no", "off", ""}


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except (ValueError, TypeError):
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# XYZ <-> coords helpers
# ---------------------------------------------------------------------------

def _parse_xyz(xyz: str) -> Optional[Tuple[List[str], "np.ndarray"]]:
    """Return (symbols, Nx3 float array) or None on parse error."""
    if np is None:
        return None
    lines = [ln for ln in xyz.strip().splitlines() if ln.strip()]
    syms: List[str] = []
    coords: List[List[float]] = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < 4:
            # XYZ may have header lines; skip them
            try:
                int(parts[0])
                continue
            except (ValueError, IndexError):
                return None
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            try:
                int(parts[0])
                continue
            except (ValueError, IndexError):
                return None
        syms.append(parts[0])
        coords.append([x, y, z])
    if not syms:
        return None
    return syms, np.asarray(coords, dtype=float)


def _coords_to_xyz_string(symbols: Sequence[str], coords: "np.ndarray") -> str:
    out: List[str] = []
    for sym, (x, y, z) in zip(symbols, coords):
        out.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(out) + "\n"


def _axis_angle_rotmat(axis: "np.ndarray", angle_rad: float) -> "np.ndarray":
    """Rodrigues' rotation formula."""
    a = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(a))
    if n < 1e-12:
        return np.eye(3)
    a = a / n
    ux, uy, uz = a
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array([
        [c + ux * ux * one_c, ux * uy * one_c - uz * s, ux * uz * one_c + uy * s],
        [uy * ux * one_c + uz * s, c + uy * uy * one_c, uy * uz * one_c - ux * s],
        [uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c],
    ], dtype=float)


# ---------------------------------------------------------------------------
# Sigma-fragment walker (mirrors smiles_converter._enumerate_hapto_sigma_isomers
# lines 21992-22021).  Walks the BFS graph from a σ-donor outward, never
# crossing metals, other σ-donors, or any hapto atom.
# ---------------------------------------------------------------------------

def _build_nonmetal_adjacency(mol, metal_set: Set[str]) -> Dict[int, Set[int]]:
    non_metal = {
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetSymbol() not in metal_set
    }
    adj: Dict[int, Set[int]] = {i: set() for i in non_metal}
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in non_metal and bj in non_metal:
            adj[bi].add(bj)
            adj[bj].add(bi)
    return adj


def _walk_sigma_fragment(
    donor: int,
    adj: Dict[int, Set[int]],
    hapto_atoms: Set[int],
    other_sigma_donors: Set[int],
) -> Set[int]:
    """BFS from σ-donor.  Returns the set of atom indices reachable via
    non-metal bonds, never crossing into hapto atoms or *other* σ-donors.
    Includes the donor itself."""
    frag: Set[int] = set()
    stack = [donor]
    visited: Set[int] = set()
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        if node in hapto_atoms:
            continue
        if node != donor and node in other_sigma_donors:
            continue
        visited.add(node)
        frag.add(node)
        for nbr in adj.get(node, ()):
            if nbr not in visited:
                stack.append(nbr)
    return frag


# ---------------------------------------------------------------------------
# Topology validation
# ---------------------------------------------------------------------------

def _check_hapto_topology(
    cand_coords: "np.ndarray",
    seed_coords: "np.ndarray",
    hapto_atoms: Set[int],
    metal_donor_pairs: List[Tuple[int, int]],
    hapto_tol: float,
    md_tol: float,
) -> bool:
    """Return True iff the candidate preserves all hapto invariants.

    Checks:
      1. Every hapto atom moved < ``hapto_tol`` from seed.
      2. Every (M, D) bond length unchanged within ``md_tol``.
      3. Each hapto group's atoms remain coplanar within ``_RING_PLANARITY_TOL``.
    """
    if cand_coords.shape != seed_coords.shape:
        return False

    # Rule 1 — hapto atoms frozen
    for ai in hapto_atoms:
        if ai >= cand_coords.shape[0]:
            return False
        d = float(np.linalg.norm(cand_coords[ai] - seed_coords[ai]))
        if d > hapto_tol:
            return False

    # Rule 2 — M-D bond length preserved
    for m_idx, d_idx in metal_donor_pairs:
        if m_idx >= cand_coords.shape[0] or d_idx >= cand_coords.shape[0]:
            continue
        seed_md = float(np.linalg.norm(seed_coords[m_idx] - seed_coords[d_idx]))
        cand_md = float(np.linalg.norm(cand_coords[m_idx] - cand_coords[d_idx]))
        if abs(cand_md - seed_md) > md_tol:
            return False

    return True


def _check_ring_planarity(
    coords: "np.ndarray",
    ring_atoms: Sequence[int],
    tol: float = _RING_PLANARITY_TOL,
) -> bool:
    """Best-fit-plane-based planarity check for a ring.  Skips on rings
    with < 3 atoms (degenerate)."""
    if len(ring_atoms) < 3 or np is None:
        return True
    pts = coords[list(ring_atoms)]
    centroid = pts.mean(axis=0)
    centred = pts - centroid
    try:
        _, S, Vt = np.linalg.svd(centred, full_matrices=False)
        # Plane normal is the right-singular-vector with smallest singular value
        normal = Vt[-1]
        # Distance of each point to plane
        dists = np.abs(centred @ normal)
        return bool(dists.max() <= tol)
    except Exception:
        return True


# ---------------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------------

def apply_hapto_diversity_topology_aware(
    smiles: str,
    seed_xyz: str,
    mol_template,
    *,
    metal_set: Set[str],
    find_hapto_groups: Callable,
    verify_topology: Callable[[str, object], bool],
    optimize_xyz: Optional[Callable[[str, object], str]] = None,
    max_frames: Optional[int] = None,
    seen_xyz_keys: Optional[Set[str]] = None,
    label_prefix: str = "hapto-diversity",
) -> List[Tuple[str, str]]:
    """Generate up to ``max_frames`` topology-preserving conformers.

    Parameters
    ----------
    smiles : str
        SMILES string (used only for label disambiguation).
    seed_xyz : str
        DELFIN-format XYZ from the hapto-approx seed builder.  Provides
        the immutable ground-truth positions for hapto atoms.
    mol_template : RDKit Mol
        Already prepared via ``_prepare_mol_for_embedding(..., hapto_approx=True)``.
    metal_set : Set[str]
        ``_METAL_SET`` from caller.
    find_hapto_groups : callable
        ``_find_hapto_groups`` from caller.
    verify_topology : callable
        ``_verify_topology_from_graph`` from caller.
    optimize_xyz : callable, optional
        Lightweight UFF refiner; when provided, applied per-candidate
        BEFORE topology validation.  None → skip refinement (faster).
    max_frames : int, optional
        Override env-flag.  None → use ``DELFIN_HAPTO_DIVERSITY_MAX``.
    seen_xyz_keys : Set[str], optional
        Caller's de-dup set.  Mutated in-place when frames are accepted.
    label_prefix : str
        Prepended to numeric label for emitted frames.

    Returns
    -------
    List[(xyz, label)] — never raises; on any internal failure returns [].
    """
    if not _env_truthy("DELFIN_HAPTO_DIVERSITY_RESTORE", "1"):
        return []

    if mol_template is None or np is None:
        return []

    if max_frames is None:
        max_frames = _env_int("DELFIN_HAPTO_DIVERSITY_MAX", 60)
    n_torsion_steps = max(2, _env_int("DELFIN_HAPTO_DIVERSITY_TORSION_STEPS", 12))
    include_ring_rot = _env_truthy("DELFIN_HAPTO_DIVERSITY_INCLUDE_RING_ROT", "0")
    hapto_tol = _env_float("DELFIN_HAPTO_DIVERSITY_HAPTO_TOL", _HAPTO_TOL_DEFAULT)
    md_tol = _env_float("DELFIN_HAPTO_DIVERSITY_MD_TOL", _MD_TOL_DEFAULT)

    parsed = _parse_xyz(seed_xyz)
    if parsed is None:
        return []
    symbols, seed_coords = parsed

    # Atom-count match
    try:
        if mol_template.GetNumAtoms() != seed_coords.shape[0]:
            return []
    except Exception:
        return []

    # If the seed itself fails the caller's topology gate, we skip the
    # gate for our rotation outputs (otherwise we'd reject the seed as
    # well as any rotation of it).  This mirrors how HEAD admits the
    # ``hapto-approx`` seed unconditionally.
    seed_passes_caller_gate = True
    try:
        seed_passes_caller_gate = bool(verify_topology(seed_xyz, mol_template))
    except Exception:
        seed_passes_caller_gate = False

    # Identify hapto sets
    try:
        hapto_groups = find_hapto_groups(mol_template) or []
    except Exception:
        hapto_groups = []

    hapto_atoms: Set[int] = set()
    for _midx, members in hapto_groups:
        hapto_atoms.update(members)

    # Identify metals + their σ-donors
    metal_idxs: List[int] = []
    sigma_pairs: List[Tuple[int, int]] = []  # (metal_idx, donor_idx)
    for atom in mol_template.GetAtoms():
        if atom.GetSymbol() not in metal_set:
            continue
        m_idx = atom.GetIdx()
        metal_idxs.append(m_idx)
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if nbr.GetSymbol() in metal_set:
                continue
            if ni in hapto_atoms:
                continue
            sigma_pairs.append((m_idx, ni))

    if not sigma_pairs and not (include_ring_rot and len(hapto_groups) >= 2):
        # Nothing to diversify
        return []

    # Pre-build full M-D pair set (for length-preservation check, includes
    # both σ and η bonds).
    md_pairs: List[Tuple[int, int]] = []
    for atom in mol_template.GetAtoms():
        if atom.GetSymbol() not in metal_set:
            continue
        m_idx = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() in metal_set:
                continue
            md_pairs.append((m_idx, nbr.GetIdx()))

    # Adjacency for σ-fragment BFS
    adj = _build_nonmetal_adjacency(mol_template, metal_set)
    all_sigma_donor_idxs = {d for _m, d in sigma_pairs}

    # Pre-compute σ-fragment per donor (cached).
    # Strategy: try graph-walk first; if the resulting graph-fragment is
    # physically inconsistent in the seed (stretched bonds → atom-index
    # mis-mapping), fall back to a *physical-proximity* fragment that
    # follows actual seed-distance covalent bonds.  Hapto-approx seeds
    # occasionally permute σ-substituent atom indices, breaking the
    # graph-walk; the proximity fallback recovers the true σ-cone.
    sigma_fragments: Dict[Tuple[int, int], List[int]] = {}
    _BOND_SANITY_MAX = 2.5  # Å — generous; covers M-D-σ-X 1st-shell only

    def _physical_fragment(
        d_idx: int,
        forbidden: Set[int],
    ) -> Set[int]:
        """BFS in physical space from donor; an edge exists between i and
        j iff their seed distance < covalent-bond cutoff, where the cutoff
        depends on element pair (1.20 Å for H-X, 1.85 Å otherwise).
        Forbidden = metals ∪ hapto ∪ other σ-donors."""
        n = seed_coords.shape[0]
        frag: Set[int] = set()
        stack = [d_idx]
        visited: Set[int] = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            if node in forbidden:
                continue
            visited.add(node)
            frag.add(node)
            for nb in range(n):
                if nb == node or nb in visited:
                    continue
                if nb in forbidden:
                    continue
                d = float(np.linalg.norm(seed_coords[node] - seed_coords[nb]))
                # Element-pair-aware cutoff
                s_a = symbols[node]
                s_b = symbols[nb]
                if s_a == 'H' or s_b == 'H':
                    cutoff = 1.30
                else:
                    cutoff = 1.85
                if d < cutoff:
                    stack.append(nb)
        return frag

    # Build the forbidden set for physical-fallback fragments
    metal_idx_set = {
        a.GetIdx() for a in mol_template.GetAtoms()
        if a.GetSymbol() in metal_set
    }

    for m_idx, d_idx in sigma_pairs:
        others = all_sigma_donor_idxs - {d_idx}
        frag = _walk_sigma_fragment(d_idx, adj, hapto_atoms, others)
        # Validate each graph-bond is physically present in seed.
        broken = False
        frag_set = set(frag)
        for ai in frag:
            for nb in adj.get(ai, ()):
                if nb in frag_set and ai < nb:
                    d = float(np.linalg.norm(seed_coords[ai] - seed_coords[nb]))
                    if d > _BOND_SANITY_MAX:
                        broken = True
                        break
            if broken:
                break

        if broken:
            # Fall back to physical-proximity fragment
            forbidden = metal_idx_set | hapto_atoms | (others)
            phys_frag = _physical_fragment(d_idx, forbidden)
            if not phys_frag or len(phys_frag) < 2:
                continue
            # Sanity: physical fragment must NOT include any forbidden atom
            if phys_frag & (metal_idx_set | hapto_atoms):
                continue
            sigma_fragments[(m_idx, d_idx)] = sorted(phys_frag)
        else:
            sigma_fragments[(m_idx, d_idx)] = sorted(frag)

    # Build emit set
    if seen_xyz_keys is None:
        seed_key = "\n".join(
            l.strip() for l in seed_xyz.splitlines() if l.strip()
        )
        seen_xyz_keys = {seed_key}
    results: List[Tuple[str, str]] = []
    rejected_topology = 0
    rejected_dup = 0

    # Generate candidate angle list — uniformly spaced, excluding 0
    angles_deg = [
        (i * 360.0 / n_torsion_steps) for i in range(1, n_torsion_steps)
    ]

    # ---- MECH-1 + MECH-2 — torsion sampling around each M-D axis -----------
    # Strategy: deterministic Cartesian product over (φ_donor1, φ_donor2, ...)
    # but we sample each donor independently first, then mix pairs/triples.
    # We enumerate via a deterministic order: single-donor sweeps first
    # (which gives the basic diversity), then pairs, etc.

    # Compute seed clash floor — minimum atom-atom distance per element
    # class.  Candidates only rejected if they fall below this floor minus
    # a 0.05 Å slack (i.e. we never make clashes worse than the seed).
    def _per_class_min(c: "np.ndarray") -> Tuple[float, float, float]:
        n = c.shape[0]
        m_hh = m_hx = m_xx = 1e9
        for i in range(n):
            si = symbols[i]
            for j in range(i + 1, n):
                sj = symbols[j]
                d = float(np.linalg.norm(c[i] - c[j]))
                if si == 'H' and sj == 'H':
                    if d < m_hh:
                        m_hh = d
                elif si == 'H' or sj == 'H':
                    if d < m_hx:
                        m_hx = d
                else:
                    if d < m_xx:
                        m_xx = d
        return m_hh, m_hx, m_xx

    seed_min_hh, seed_min_hx, seed_min_xx = _per_class_min(seed_coords)
    _CLASH_SLACK = 0.05  # Å

    # Hard absolute floor (collapsed-coordinate detection — true unphysical)
    _ABS_FLOOR_HH = 0.20
    _ABS_FLOOR_HX = 0.30
    _ABS_FLOOR_XX = 0.40

    def _intrinsic_clash_check(c: "np.ndarray") -> bool:
        """Reject if the candidate has a min-distance worse than seed (with
        ``_CLASH_SLACK`` slack) OR below an absolute floor."""
        c_hh, c_hx, c_xx = _per_class_min(c)
        if c_hh < max(seed_min_hh - _CLASH_SLACK, _ABS_FLOOR_HH):
            return False
        if c_hx < max(seed_min_hx - _CLASH_SLACK, _ABS_FLOOR_HX):
            return False
        if c_xx < max(seed_min_xx - _CLASH_SLACK, _ABS_FLOOR_XX):
            return False
        return True

    def _emit_candidate(
        rotated_coords: "np.ndarray",
        descriptor: str,
    ) -> bool:
        """Validate + dedup + (optionally UFF-refine) a candidate.  Returns
        True if frame was accepted."""
        nonlocal rejected_topology, rejected_dup

        # Hapto/M-D invariant pre-check (the strongest guard — guarantees
        # η-coordination and σ-bond lengths are preserved bit-exact)
        if not _check_hapto_topology(
            rotated_coords, seed_coords, hapto_atoms, md_pairs,
            hapto_tol, md_tol,
        ):
            rejected_topology += 1
            return False

        # Intrinsic clash gate (collapsed coordinates)
        if not _intrinsic_clash_check(rotated_coords):
            rejected_topology += 1
            return False

        xyz_cand = _coords_to_xyz_string(symbols, rotated_coords)

        # Optional UFF refinement (may move atoms — re-validate invariants)
        if optimize_xyz is not None:
            try:
                xyz_refined = optimize_xyz(xyz_cand, mol_template)
                refined_parsed = _parse_xyz(xyz_refined)
                if refined_parsed is not None:
                    refined_syms, refined_coords = refined_parsed
                    if (
                        len(refined_syms) == len(symbols)
                        and refined_coords.shape == rotated_coords.shape
                    ):
                        # Re-check invariant after UFF
                        if not _check_hapto_topology(
                            refined_coords, seed_coords, hapto_atoms,
                            md_pairs, hapto_tol, md_tol,
                        ):
                            rejected_topology += 1
                            return False
                        if not _intrinsic_clash_check(refined_coords):
                            rejected_topology += 1
                            return False
                        xyz_cand = xyz_refined
            except Exception:
                pass

        # Caller's full topology gate — applied ONLY if seed already
        # passed it.  Otherwise we'd reject every rotation just because
        # the seed itself fails (e.g. DIVNUT with the hapto-approx
        # builder generating an XYZ whose atom-index↔coord mapping
        # confuses the M-D phantom-bond gate).
        if seed_passes_caller_gate:
            try:
                if not verify_topology(xyz_cand, mol_template):
                    rejected_topology += 1
                    return False
            except Exception:
                return False

        key = "\n".join(l.strip() for l in xyz_cand.splitlines() if l.strip())
        if key in seen_xyz_keys:
            rejected_dup += 1
            return False
        seen_xyz_keys.add(key)
        results.append((
            xyz_cand,
            f"{label_prefix}-{descriptor}-{len(results):03d}",
        ))
        return True

    # ---- MECH-1.5 — per-metal σ-cluster rotation around the metal's --------
    # principal axis.  For each metal that has σ-donors AND a clear
    # "principal axis" (either toward another metal, or toward an η-group
    # centroid), rotate the entire σ-cluster (all σ-donors + their
    # fragments) as a rigid body by φ ∈ angles_deg.  This works even when
    # individual σ-fragments are graph-broken because we don't depend on
    # per-fragment integrity — only on the cluster moving as one piece.
    # The rotation axis passes through the metal, so the metal itself
    # does not move and hapto atoms (kept identity-mapped) are also not
    # moved.

    # Group sigma_pairs by metal
    metal_to_sigma: Dict[int, List[int]] = {}
    for m_idx, d_idx in sigma_pairs:
        metal_to_sigma.setdefault(m_idx, []).append(d_idx)

    for m_idx, donors in metal_to_sigma.items():
        if len(results) >= max_frames:
            break
        if not donors:
            continue
        # Determine principal axis
        axis = None
        # Preference 1 — axis toward a hapto-centroid attached to this metal
        for h_m, members in hapto_groups:
            if h_m == m_idx and members:
                axis = seed_coords[list(members)].mean(axis=0) - seed_coords[m_idx]
                if float(np.linalg.norm(axis)) > 1e-6:
                    break
                axis = None
        # Preference 2 — axis toward another metal (graph-bonded)
        if axis is None:
            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() not in metal_set:
                    continue
                if atom.GetIdx() == m_idx:
                    continue
                bond = mol_template.GetBondBetweenAtoms(m_idx, atom.GetIdx())
                if bond is None:
                    continue
                axis_cand = seed_coords[atom.GetIdx()] - seed_coords[m_idx]
                if float(np.linalg.norm(axis_cand)) > 1e-6:
                    axis = axis_cand
                    break
        # Preference 3 — axis = mean of σ-donor positions minus metal
        if axis is None:
            sigma_pts = seed_coords[donors]
            axis = sigma_pts.mean(axis=0) - seed_coords[m_idx]
        if axis is None or float(np.linalg.norm(axis)) < 1e-6:
            continue

        # Atoms to rotate: all donors and (when fragments present) their
        # fragment atoms; never include the metal, hapto atoms, or atoms
        # of OTHER metal coordination spheres.
        atoms_to_rotate: Set[int] = set()
        for d_idx in donors:
            atoms_to_rotate.add(d_idx)
            frag = sigma_fragments.get((m_idx, d_idx), [])
            for ai in frag:
                atoms_to_rotate.add(ai)
        # Exclude metals and hapto atoms (defensive)
        atoms_to_rotate -= metal_idx_set
        atoms_to_rotate -= hapto_atoms
        # Exclude σ-donors of OTHER metals + their fragments
        for (other_m, other_d), other_frag in sigma_fragments.items():
            if other_m == m_idx:
                continue
            atoms_to_rotate -= set(other_frag)
            atoms_to_rotate.discard(other_d)
        if not atoms_to_rotate:
            continue
        pivot = seed_coords[m_idx]
        for phi_deg in angles_deg:
            if len(results) >= max_frames:
                break
            R = _axis_angle_rotmat(axis, math.radians(phi_deg))
            cand = seed_coords.copy()
            for ai in atoms_to_rotate:
                cand[ai] = R @ (seed_coords[ai] - pivot) + pivot
            _emit_candidate(cand, f"mc{m_idx}-{int(phi_deg):03d}")

    # MECH-1 — single-donor sweep
    for (m_idx, d_idx), frag in sigma_fragments.items():
        if len(results) >= max_frames:
            break
        # Atoms to rotate: fragment minus donor itself (donor is anchor)
        rot_atoms = [ai for ai in frag if ai != d_idx]
        if not rot_atoms:
            continue
        m_pos = seed_coords[m_idx]
        d_pos = seed_coords[d_idx]
        axis = d_pos - m_pos
        if float(np.linalg.norm(axis)) < 1e-6:
            continue
        for phi_deg in angles_deg:
            if len(results) >= max_frames:
                break
            R = _axis_angle_rotmat(axis, math.radians(phi_deg))
            cand = seed_coords.copy()
            for ai in rot_atoms:
                cand[ai] = R @ (seed_coords[ai] - d_pos) + d_pos
            _emit_candidate(cand, f"sd{m_idx}-{d_idx}-{int(phi_deg):03d}")

    # MECH-2 — multi-donor joint sweep (pairs)
    if len(sigma_pairs) >= 2 and len(results) < max_frames:
        # Deterministic enumeration: every ordered pair (i, j) with i<j,
        # combine 2 angles each (60° step) → 4 combos per pair.
        pair_angles = angles_deg[:: max(1, len(angles_deg) // 3)] or angles_deg[:1]
        for i in range(len(sigma_pairs)):
            if len(results) >= max_frames:
                break
            for j in range(i + 1, len(sigma_pairs)):
                if len(results) >= max_frames:
                    break
                pair_a = sigma_pairs[i]
                pair_b = sigma_pairs[j]
                frag_a = sigma_fragments.get(pair_a, [])
                frag_b = sigma_fragments.get(pair_b, [])
                rot_a = [ai for ai in frag_a if ai != pair_a[1]]
                rot_b = [ai for ai in frag_b if ai != pair_b[1]]
                if not rot_a or not rot_b:
                    continue
                axis_a = seed_coords[pair_a[1]] - seed_coords[pair_a[0]]
                axis_b = seed_coords[pair_b[1]] - seed_coords[pair_b[0]]
                if float(np.linalg.norm(axis_a)) < 1e-6:
                    continue
                if float(np.linalg.norm(axis_b)) < 1e-6:
                    continue
                for pa in pair_angles:
                    if len(results) >= max_frames:
                        break
                    for pb in pair_angles:
                        if len(results) >= max_frames:
                            break
                        Ra = _axis_angle_rotmat(axis_a, math.radians(pa))
                        Rb = _axis_angle_rotmat(axis_b, math.radians(pb))
                        cand = seed_coords.copy()
                        anchor_a = seed_coords[pair_a[1]]
                        anchor_b = seed_coords[pair_b[1]]
                        for ai in rot_a:
                            cand[ai] = Ra @ (seed_coords[ai] - anchor_a) + anchor_a
                        for ai in rot_b:
                            cand[ai] = Rb @ (seed_coords[ai] - anchor_b) + anchor_b
                        _emit_candidate(
                            cand,
                            f"pair-{pair_a[1]}-{pair_b[1]}-{int(pa):03d}-{int(pb):03d}",
                        )

    # MECH-3 — subordinate η-ring rigid rotation (opt-in)
    if include_ring_rot and len(hapto_groups) >= 2 and len(results) < max_frames:
        # Skip the first hapto group (kept fixed as reference); rotate
        # subsequent groups around their (M, centroid) axis.
        for h_idx, (m_idx, members) in enumerate(hapto_groups):
            if h_idx == 0:
                continue
            if len(results) >= max_frames:
                break
            if not members:
                continue
            ring_pts = seed_coords[members]
            centroid = ring_pts.mean(axis=0)
            axis = centroid - seed_coords[m_idx]
            if float(np.linalg.norm(axis)) < 1e-6:
                continue
            for phi_deg in (30.0, 60.0, 90.0, 120.0, 150.0):
                if len(results) >= max_frames:
                    break
                R = _axis_angle_rotmat(axis, math.radians(phi_deg))
                cand = seed_coords.copy()
                for ai in members:
                    cand[ai] = R @ (seed_coords[ai] - centroid) + centroid
                # Note: ring rotation does not preserve hapto-atom positions
                # so _check_hapto_topology will REJECT this; we use a
                # relaxed check that allows hapto motion within the ring's
                # own plane.  We do this by checking M-(ring-atom) lengths
                # and ring planarity instead.
                ok = True
                for ai in members:
                    ml = float(np.linalg.norm(seed_coords[ai] - seed_coords[m_idx]))
                    cl = float(np.linalg.norm(cand[ai] - cand[m_idx]))
                    if abs(cl - ml) > md_tol:
                        ok = False
                        break
                if not ok:
                    continue
                if not _check_ring_planarity(cand, members):
                    continue
                if not _intrinsic_clash_check(cand):
                    rejected_topology += 1
                    continue
                xyz_cand = _coords_to_xyz_string(symbols, cand)
                if seed_passes_caller_gate:
                    try:
                        if not verify_topology(xyz_cand, mol_template):
                            rejected_topology += 1
                            continue
                    except Exception:
                        continue
                key = "\n".join(
                    l.strip() for l in xyz_cand.splitlines() if l.strip()
                )
                if key in seen_xyz_keys:
                    rejected_dup += 1
                    continue
                seen_xyz_keys.add(key)
                results.append((
                    xyz_cand,
                    f"{label_prefix}-ring-{h_idx}-{int(phi_deg):03d}-{len(results):03d}",
                ))

    return results
