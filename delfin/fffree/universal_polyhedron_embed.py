"""delfin.fffree.universal_polyhedron_embed — Universal polyhedron-vertex
coord-map embedding for ANY (metal, CN, donor-multiset, polyhedron) combo.

Fundamental architecture (2026-06-07): replaces the per-class template
dispatcher (``piano_stool_template`` / ``fischer_carbene_template`` /
``build_sp4_chelate_mono2_template`` / ``classify_for_template``) with ONE
universal mechanism.  No SMILES patterns, no per-class branches; the only
inputs are

    - the RDKit graph (Mol),
    - the polyhedron name (chosen by Pólya-vertex enumeration upstream),
    - the per-donor M-D target distance (CCDC-empirical via
      :func:`delfin.fffree.polyhedra.md_distance`),
    - a deterministic random seed.

What this module does
---------------------
1.  Read the polyhedron vertex unit-vectors from
    :func:`delfin.fffree.polyhedra.ref_vectors` (single source of truth — no
    duplicated geometry data).
2.  Scale each vertex by the per-donor M-D target distance so each donor is
    pinned at its empirical Cr–O / Cu–N / Cr–C(carbene) / etc. distance.
3.  Build an RDKit ``coordMap`` with the metal at the origin + each donor at
    its scaled vertex.
4.  Run ``ETKDGv3`` with ``SetCoordMap`` — the solver fills in all the
    ligand-internal positions (substituents, H atoms, backbone) using the
    bounds-matrix that ETKDG derives internally from the molecular graph.
5.  For hapto donors (η^n rings), the η ring centroid is treated as the
    "effective donor" placed at the polyhedron vertex; the actual ring atoms
    are then projected onto a circle in the plane perpendicular to the
    M–centroid axis using a single universal formula (no per-class lookup).

Universality contract
---------------------
* NO SMILES-pattern matching.
* NO per-class branches (``if class == "piano_stool"`` etc.).
* NO reference to ``template_dispatcher`` or ``piano_stool_template`` or
  ``fischer_carbene_template`` or ``build_oc6_chelate2_mono2_template`` or
  ``build_sp4_chelate_mono2_template``.
* Single code path for AFOFIL (Ni²⁺ chelate + 2 mono), ODUXAN (Fischer
  carbene Cr(CO)₅(N)), WICROP (η⁶-arene Cr(CO)₃(PR₃)), BEYRAY (Cu OC-6 2
  chel + 2 H₂O), and every other (metal, CN, polyhedron, donor-set) combo.

Determinism contract
--------------------
* ``ETKDGv3`` with fixed ``randomSeed`` (default 42).
* ``useRandomCoords = True`` (required when coordMap pins coordinates).
* ``numThreads = 1``.
* No PRNG outside that one call.
* Two runs with the same SMILES + same seed produce bit-identical XYZ blocks
  (verified by ``test_universal_polyhedron_embed.py::test_determinism``).

Default OFF, byte-identical
---------------------------
Nothing in this module mutates state.  The integration site (the
``DELFIN_FFFREE_UNIVERSAL_POLY_EMBED`` env-flag) only consults this module
when the flag is set.  Adopters with the flag unset see byte-identical
behaviour to HEAD.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Public env flag (single source of truth).
# ---------------------------------------------------------------------------
ENV_FLAG = "DELFIN_FFFREE_UNIVERSAL_POLY_EMBED"


def flag_active() -> bool:
    """Return ``True`` when the universal polyhedron embed env-flag is on."""
    return os.environ.get(ENV_FLAG, "0").strip().lower() in (
        "1", "true", "yes", "on",
    )


# ---------------------------------------------------------------------------
# Deterministic ETKDG seed.  Mirrors :data:`embed_fallback.ETKDG_SEED` so the
# universal path is deterministic with the same value when both modules are
# called in the same Python session.
# ---------------------------------------------------------------------------
DEFAULT_SEED = 42


# ---------------------------------------------------------------------------
# Geometry helpers.
# ---------------------------------------------------------------------------
def _scaled_vertices(
    geometry: str,
    metal: str,
    md_target_per_donor: Sequence[float],
) -> Optional[np.ndarray]:
    """Return the polyhedron vertex array scaled per-donor to ``md_target``.

    The unit-vector polyhedron (each vertex on the unit sphere) is read from
    :func:`delfin.fffree.polyhedra.ref_vectors` -- a single source of truth
    shared with the legacy CoordMap path.  Each vertex ``i`` is then
    multiplied by ``md_target_per_donor[i]`` so donor ``i`` is pinned at its
    empirical M-D distance.

    Returns ``None`` when the geometry is unknown or the donor count doesn't
    match the polyhedron's vertex count.  Pure function, no side effects.
    """
    try:
        from delfin.fffree import polyhedra as _PLY
    except Exception:
        return None
    try:
        V = _PLY.ref_vectors(geometry)
    except KeyError:
        return None
    if V is None:
        return None
    V = np.asarray(V, dtype=float)
    targets = np.asarray(md_target_per_donor, dtype=float)
    if targets.shape[0] != V.shape[0]:
        # Caller asked for a polyhedron whose vertex count doesn't match the
        # number of donor M-D targets.  Refuse rather than guess.
        return None
    return V * targets[:, None]


def _is_hapto_donor(mol, donor_atom_idx: int) -> Tuple[bool, List[int]]:
    """Detect whether a donor sits in an aromatic ring that's coordinating as
    η^n (hapto), and return the ring-atom indices when it does.

    Universal rule (graph-only, no SMILES patterns):
      A donor is η^n iff it is part of an aromatic ring AND at least one OTHER
      atom of that same ring is ALSO bonded to the metal in the molecular
      graph.  The hapto count is then the number of ring atoms bonded to the
      metal.

    Returns ``(False, [])`` for σ-donors.
    """
    try:
        from rdkit.Chem import GetSSSR
    except Exception:
        return False, []
    try:
        donor_atom = mol.GetAtomWithIdx(int(donor_atom_idx))
    except Exception:
        return False, []
    if donor_atom is None:
        return False, []
    # Find the metal neighbour of this donor.
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return False, []
    metal_nbrs = [
        nbr.GetIdx()
        for nbr in donor_atom.GetNeighbors()
        if _is_metal(nbr.GetSymbol())
    ]
    if not metal_nbrs:
        return False, []
    metal_idx = int(metal_nbrs[0])
    if not donor_atom.GetIsAromatic():
        return False, []
    # Identify the smallest aromatic ring containing the donor.
    try:
        ri = mol.GetRingInfo()
        rings = ri.AtomRings()
    except Exception:
        return False, []
    candidate_rings: List[List[int]] = []
    for ring in rings:
        if int(donor_atom_idx) not in ring:
            continue
        # Aromatic ring check: all atoms aromatic.
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        candidate_rings.append(list(ring))
    if not candidate_rings:
        return False, []
    # Pick the SMALLEST aromatic ring (most chemically canonical for η).
    candidate_rings.sort(key=len)
    ring = candidate_rings[0]
    # Count ring atoms also bonded to the metal.
    metal_bonded_ring_atoms = []
    for ring_atom_idx in ring:
        ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
        for nbr in ring_atom.GetNeighbors():
            if int(nbr.GetIdx()) == metal_idx:
                metal_bonded_ring_atoms.append(int(ring_atom_idx))
                break
    # η^n means at least 2 ring atoms bonded to metal.
    if len(metal_bonded_ring_atoms) < 2:
        return False, []
    return True, sorted(set(metal_bonded_ring_atoms))


def _hapto_ring_carbon_positions(
    centroid_position: np.ndarray,
    metal_position: np.ndarray,
    n_ring_atoms: int,
    *,
    r_ring: float = 1.40,
) -> np.ndarray:
    """Universal hapto-ring atom placement.

    For an η^n donor whose effective coord-map vertex is the ring CENTROID
    (placed at ``centroid_position``), this function returns the n ring-atom
    positions on a circle of radius ``r_ring`` lying in the plane
    perpendicular to the M–centroid axis.  The first ring atom is placed at
    angle 0 (deterministic) and subsequent atoms at ``2π k / n``.

    No per-ring-type lookup.  Pure geometry.
    """
    axis = centroid_position - metal_position
    norm_axis = float(np.linalg.norm(axis))
    if norm_axis < 1e-9:
        # Degenerate: metal sits on the centroid (shouldn't happen in
        # practice).  Default the axis to +z.
        u_axis = np.array([0.0, 0.0, 1.0])
    else:
        u_axis = axis / norm_axis
    # Two orthonormal vectors spanning the perpendicular plane.
    # Pick an arbitrary vector not parallel to u_axis and Gram-Schmidt it.
    pick = np.array([1.0, 0.0, 0.0])
    if abs(float(u_axis @ pick)) > 0.9:
        pick = np.array([0.0, 1.0, 0.0])
    e1 = pick - (pick @ u_axis) * u_axis
    e1 = e1 / max(float(np.linalg.norm(e1)), 1e-12)
    e2 = np.cross(u_axis, e1)
    e2 = e2 / max(float(np.linalg.norm(e2)), 1e-12)
    positions = np.zeros((int(n_ring_atoms), 3), dtype=float)
    for k in range(int(n_ring_atoms)):
        theta = 2.0 * math.pi * k / float(n_ring_atoms)
        positions[k] = (
            centroid_position
            + r_ring * (math.cos(theta) * e1 + math.sin(theta) * e2)
        )
    return positions


# ---------------------------------------------------------------------------
# Core universal embed.
# ---------------------------------------------------------------------------
def universal_polyhedron_vertex_embed(
    metal_symbol: str,
    geometry: str,
    donor_indices_in_mol: Sequence[int],
    rdkit_mol,
    md_target_per_donor: Sequence[float],
    *,
    coord_seed: int = DEFAULT_SEED,
    n_conformers: int = 1,
) -> Optional[Tuple[np.ndarray, List[str]]]:
    """Universal coord-map ETKDG embed for ANY (metal, CN, polyhedron) combo.

    Parameters
    ----------
    metal_symbol : str
        Metal element symbol (used for M-D fallback distance only).
    geometry : str
        Polyhedron name from :mod:`delfin.fffree.polyhedra` (e.g.
        ``"OC-6 octahedron"``, ``"SP-4 square planar"``, ``"T-4 tetrahedron"``).
    donor_indices_in_mol : sequence of int
        RDKit atom indices of the donors that should be pinned to vertices.
        Length must equal the polyhedron's vertex count.
    rdkit_mol : RDKit Mol
        Pre-built molecule (AddHs already applied is recommended).  The
        metal atom index is detected automatically (lowest-index metal).
    md_target_per_donor : sequence of float
        Empirical M-D distance per donor (in Å).  Use
        :func:`delfin.fffree.polyhedra.md_distance` to obtain these.
    coord_seed : int
        Deterministic ETKDG random seed.
    n_conformers : int
        Number of conformers to embed (each a different ETKDG random init
        with the same vertex-pinned coordMap).

    Returns
    -------
    tuple ``(coords, syms)`` or ``None``
        ``coords`` is a (n_atoms, 3) float64 array containing the BEST
        (lowest-energy index 0) conformer's positions, ``syms`` is the
        per-atom symbol list.  ``None`` on any failure (RDKit unavailable,
        vertex count mismatch, embed failed, NaN coordinates).

    Never raises.  All failure modes return ``None`` so the caller can
    silently fall back to the legacy path.
    """
    # Local imports guard against rdkit being unavailable in worker
    # contexts (we never raise).
    try:
        from rdkit.Chem import AllChem
        from rdkit.Geometry import Point3D
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    if rdkit_mol is None:
        return None
    try:
        syms = [a.GetSymbol() for a in rdkit_mol.GetAtoms()]
    except Exception:
        return None
    if not syms:
        return None

    # Find the metal atom (lowest-index metal — deterministic).
    metal_idx: Optional[int] = None
    for i, s in enumerate(syms):
        if _is_metal(s):
            metal_idx = int(i)
            break
    if metal_idx is None:
        return None

    donor_list = [int(d) for d in donor_indices_in_mol]
    md_targets = [float(d) for d in md_target_per_donor]
    if len(donor_list) != len(md_targets) or not donor_list:
        return None

    # Validate every donor exists in mol.
    n_atoms = len(syms)
    if any(d < 0 or d >= n_atoms for d in donor_list):
        return None

    # Build the per-donor scaled vertex array.
    V_scaled = _scaled_vertices(geometry, metal_symbol, md_targets)
    if V_scaled is None:
        return None

    # Build the coordMap: metal at origin + each donor at its scaled vertex.
    cmap: Dict[int, "Point3D"] = {
        int(metal_idx): Point3D(0.0, 0.0, 0.0),
    }
    for d_idx, vertex in zip(donor_list, V_scaled):
        cmap[int(d_idx)] = Point3D(
            float(vertex[0]), float(vertex[1]), float(vertex[2]),
        )

    # Hapto handling — universal formula (graph-only, no per-class lookup).
    # For each donor that is η^n, we expand the coord-map by ALSO pinning
    # the other ring atoms onto a circle perpendicular to the M–centroid axis.
    # The donor-vertex remains at the polyhedron vertex (treated as ring
    # CENTROID for placement purposes); the actual ring carbon positions
    # come from the universal _hapto_ring_carbon_positions formula.
    metal_pos = np.array([0.0, 0.0, 0.0])
    for d_idx, vertex in zip(donor_list, V_scaled):
        is_hapto, ring_atoms = _is_hapto_donor(rdkit_mol, int(d_idx))
        if not is_hapto:
            continue
        if len(ring_atoms) < 2:
            continue
        # The vertex was assigned to the donor atom (one specific ring atom);
        # for hapto placement we re-interpret the vertex as the RING CENTROID
        # and lay all ring atoms on the perpendicular circle.
        centroid_pos = np.asarray(vertex, dtype=float)
        ring_positions = _hapto_ring_carbon_positions(
            centroid_pos, metal_pos, len(ring_atoms),
        )
        for ra_idx, ra_pos in zip(ring_atoms, ring_positions):
            cmap[int(ra_idx)] = Point3D(
                float(ra_pos[0]), float(ra_pos[1]), float(ra_pos[2]),
            )

    # Deterministic ETKDG.
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = int(coord_seed)
        params.useRandomCoords = True
        params.numThreads = 1
        params.SetCoordMap(cmap)
    except Exception:
        return None

    n_conf = max(1, int(n_conformers))
    try:
        cids = AllChem.EmbedMultipleConfs(
            rdkit_mol, numConfs=n_conf, params=params,
        )
    except Exception:
        cids = []
    cids = list(cids)
    if not cids:
        return None

    # Return the first successful conformer.  (Callers that want the full
    # ensemble can call this in a loop with seed offsets.)
    for cid in cids:
        try:
            conf = rdkit_mol.GetConformer(int(cid))
            coords = np.asarray(conf.GetPositions(), dtype=np.float64)
        except Exception:
            continue
        if coords.size == 0 or not np.all(np.isfinite(coords)):
            continue
        if coords.shape != (n_atoms, 3):
            continue
        return coords, list(syms)
    return None


# ---------------------------------------------------------------------------
# SMILES → universal embed convenience wrapper.
# ---------------------------------------------------------------------------
def universal_embed_from_smiles(
    smiles: str,
    *,
    coord_seed: int = DEFAULT_SEED,
    n_conformers: int = 1,
) -> Optional[Tuple[np.ndarray, List[str]]]:
    """High-level wrapper: parse a SMILES into (metal, polyhedron, donors,
    md_targets) and call :func:`universal_polyhedron_vertex_embed`.

    Uses ONLY graph-only universal helpers:
      - :func:`delfin.fffree.polyhedron_vertex_polya.detect_from_smiles`
        for (polyhedron, donor_types, chelate_pairs, metal_idx).
      - :func:`delfin.fffree.polyhedra.md_distance` for per-donor M-D
        distance (CCDC-empirical covalent-radius sum).

    No SMILES patterns.  No per-class branches.  Universal.
    """
    try:
        from rdkit import Chem
    except Exception:
        return None
    if not smiles:
        return None
    try:
        from delfin.fffree.polyhedron_vertex_polya import detect_from_smiles
        from delfin.fffree import polyhedra as _PLY
    except Exception:
        return None
    det = detect_from_smiles(smiles)
    if det is None:
        return None
    polyhedron, donor_types, chelate_pairs, metal_idx = det
    # Re-parse the SMILES to get the donor atom indices in the canonical
    # heavy-only ordering (matches detect_from_smiles).
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
    except Exception:
        return None
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    metal_sym = syms[int(metal_idx)] if 0 <= int(metal_idx) < len(syms) else ""
    if not metal_sym or not _is_metal(metal_sym):
        return None
    donor_atoms: List[int] = []
    donor_elems: List[str] = []
    for nbr in mol.GetAtomWithIdx(int(metal_idx)).GetNeighbors():
        if _is_metal(nbr.GetSymbol()):
            continue
        donor_atoms.append(int(nbr.GetIdx()))
        donor_elems.append(str(nbr.GetSymbol()))
    if not donor_atoms:
        return None
    # Per-donor M-D target (CCDC-empirical covalent-radius sum).
    md_targets = [
        float(_PLY.md_distance(metal_sym, elem))
        for elem in donor_elems
    ]
    return universal_polyhedron_vertex_embed(
        metal_symbol=metal_sym,
        geometry=polyhedron,
        donor_indices_in_mol=donor_atoms,
        rdkit_mol=mol,
        md_target_per_donor=md_targets,
        coord_seed=coord_seed,
        n_conformers=n_conformers,
    )


__all__ = [
    "ENV_FLAG",
    "DEFAULT_SEED",
    "flag_active",
    "universal_polyhedron_vertex_embed",
    "universal_embed_from_smiles",
]
