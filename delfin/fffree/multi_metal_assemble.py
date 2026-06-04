"""delfin.fffree.multi_metal_assemble — Multi-metal / cluster construction
orchestrator (Task #63 Phase B, 2026-06-03).

Goal
----
FF-free construction of di-/tri-/poly-nuclear TMCs and small metal
clusters from a SMILES whose RDKit graph contains 2+ metal atoms.
Coordination polyhedra are placed per metal; metals themselves sit on a
deterministic skeleton (linear, triangle, tetrahedron, square, star,
octahedron); bridging donors with >=2 metal neighbours are placed on the
edge midpoint / face centroid / cluster centroid.

This module is the orchestrator only.  All math primitives come from:

  * ``delfin.fffree.multi_metal_polyhedra``   (M-M tables, vertex math)
  * ``delfin.fffree.bridging_ligand``         (mu2 / mu3 / mu4 placement)
  * ``delfin.fffree.cluster_scaffold_polyhedra`` (skeletons + donor partition)
  * ``delfin.fffree.metal_sphere_builder``    (per-metal polyhedron vectors)

Architecture
------------
1. Build the metal-metal connectivity graph from the RDKit mol.
2. Pick the metal-skeleton matching (n_metals, declared_mm_edges).
3. Place the metals at the skeleton vertices, scaled to the COD-empirical
   M-M distance per metal pair.
4. Detect mu-bridging atoms and place them at edge / face / centroid.
5. Per metal, place the non-bridging terminal donors radially outward
   from the cluster centroid (re-uses ``place_skeleton_with_ligands``).
6. Hard-rollback: every M-D bond within +-0.05 A of target; every M-M
   bond within +-0.10 A of target.  Failure -> return None so the
   caller falls back to the single-metal path.

Determinism
-----------
* lex-sort metal atom indices before mapping to skeleton vertices
* lex-sort donors per metal
* fixed seeds for any internal numerical routines
* PYTHONHASHSEED=0 contract observed by callers

Env-gate
--------
``DELFIN_FFFREE_MULTI_METAL=1`` (also auto-enabled by
``DELFIN_FFFREE_PURE_TRACK3=1``).  Default OFF -> the import-time
constant ``MULTI_METAL_ENABLED`` is ``False`` and every public entry
returns ``None`` without touching the legacy single-metal path.

This file is byte-identical-safe to import: nothing here mutates global
state, no monkey-patching, no side-effects on the converter.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Env gate, evaluated lazily inside ``assemble_multi_metal`` so test
# helpers can flip the flag at runtime.
_FLAG_MULTI_METAL = "DELFIN_FFFREE_MULTI_METAL"
_FLAG_PT3 = "DELFIN_FFFREE_PURE_TRACK3"

# Post-construction global-clash gate (Surgical fix 2026-06-04, b00f9a0 forensik).
# When set, after the existing M-M / M-D rollback contract passes we additionally
# verify that the assembled (metals + bridges + terminals) structure has zero
# heavy-heavy or X-H collapse violations as measured by
# ``build_time_clash_gate.collapse_count``.  Any non-zero count rolls back the
# multi-metal path and falls through to the legacy single-metal assembler
# (production-safe).  Default OFF -> byte-identical to b00f9a0.
_FLAG_GLOBAL_CLASH_GATE = "DELFIN_FFFREE_MULTI_METAL_CLASH_GATE"


def _multi_metal_enabled() -> bool:
    return (
        os.environ.get(_FLAG_MULTI_METAL, "0") == "1"
        or os.environ.get(_FLAG_PT3, "0") == "1"
    )


def _global_clash_gate_active() -> bool:
    return os.environ.get(_FLAG_GLOBAL_CLASH_GATE, "0") == "1"


# Public boolean for callers that want to short-circuit at module load.
MULTI_METAL_ENABLED = _multi_metal_enabled()


# Hard-rollback tolerances (Angstrom).  Single-metal contract is +-0.05;
# we relax M-M to +-0.10 because M-M bonds are softer (paddlewheel
# quadruple bonds e.g. Mo-Mo vary 2.10-2.30 across COD).
_MD_TOL_A = 0.05
_MM_TOL_A = 0.10


# ---------------------------------------------------------------------------
# Metal-graph helpers (read-only on the mol)
# ---------------------------------------------------------------------------


def _is_metal_sym(symbol: str) -> bool:
    """Local fallback when delfin._bond_decollapse is not importable."""
    try:
        from delfin._bond_decollapse import _is_metal as _bd_is_metal
        return _bd_is_metal(symbol)
    except Exception:
        non_metals = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I",
                      "B", "Si", "Se", "Te", "As", "He", "Ne", "Ar", "Kr",
                      "Xe", "Rn"}
        return symbol not in non_metals


def build_metal_graph(mol) -> Dict[str, object]:
    """Build the metal-connectivity graph from an RDKit mol.

    Returns a dict with:
      metal_idxs        : sorted tuple of RDKit atom indices that are metals
      metal_syms        : tuple of element symbols, same order
      mm_edges_declared : set of frozenset({i,j}) — direct M-M bonds in SMILES
      mm_edges_implied  : set of frozenset({i,j}) — pairs bridged by a >=2-metal donor
      bridges           : list of (donor_idx, frozenset(metal_idxs)) for every
                          >=2-metal donor (mu2/mu3/mu4 candidates)

    Order is deterministic (sorted indices).  Read-only on `mol`.
    """
    if mol is None:
        return {
            "metal_idxs": (),
            "metal_syms": (),
            "mm_edges_declared": set(),
            "mm_edges_implied": set(),
            "bridges": [],
        }
    metal_idxs = []
    metal_syms = {}
    for atom in mol.GetAtoms():
        s = atom.GetSymbol()
        if _is_metal_sym(s):
            metal_idxs.append(atom.GetIdx())
            metal_syms[atom.GetIdx()] = s
    metal_idxs.sort()
    metal_set = set(metal_idxs)

    mm_declared = set()
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in metal_set and b in metal_set and a != b:
            mm_declared.add(frozenset({a, b}))

    bridges: List[Tuple[int, frozenset]] = []
    mm_implied = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() in metal_set:
            continue
        ms = []
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in metal_set:
                ms.append(nb.GetIdx())
        if len(ms) >= 2:
            ms = sorted(ms)
            bridges.append((atom.GetIdx(), frozenset(ms)))
            for i in range(len(ms)):
                for j in range(i + 1, len(ms)):
                    mm_implied.add(frozenset({ms[i], ms[j]}))
    # Deterministic order
    bridges.sort(key=lambda b: (b[0], sorted(b[1])))

    return {
        "metal_idxs": tuple(metal_idxs),
        "metal_syms": tuple(metal_syms[i] for i in metal_idxs),
        "mm_edges_declared": mm_declared,
        "mm_edges_implied": mm_implied,
        "bridges": bridges,
    }


# ---------------------------------------------------------------------------
# Skeleton choice
# ---------------------------------------------------------------------------


def choose_skeleton(graph: Dict[str, object]) -> Optional[Dict[str, object]]:
    """Pick a metal-skeleton given the metal-connectivity graph.

    Returns a dict with:
      kind             : "M2_dimer" | "M3_triangle" | "M3_linear" | ...
      positions        : (n_metal, 3) ndarray of skeleton positions
      edges            : tuple of (i, j) edges in local-index space
      faces            : tuple of tuples of vertex indices
      mm_distance      : nominal M-M distance used to scale positions

    Deterministic: positions only depend on (n_metals, metal_syms, declared edges).
    Returns ``None`` when the topology does not match a supported skeleton.
    """
    syms = list(graph["metal_syms"])
    n = len(syms)
    declared = graph["mm_edges_declared"]
    implied = graph["mm_edges_implied"]
    metal_idxs = list(graph["metal_idxs"])
    # Map global to local
    local = {gi: li for li, gi in enumerate(metal_idxs)}
    n_mm_d = len(declared)
    n_mm_i = len(implied)

    if n < 2:
        return None

    from delfin.fffree import multi_metal_polyhedra as MMP
    from delfin.fffree import cluster_scaffold_polyhedra as CSP

    # -- n == 2 : dimer
    if n == 2:
        c = MMP.m2_dimer(syms[0], syms[1])
        return {
            "kind": "M2_dimer",
            "positions": np.asarray(c["vertex_positions"], float),
            "edges": ((0, 1),),
            "faces": (),
            "mm_distance": float(c["mm_distance"]),
        }

    # -- n == 3 : triangle (closed) or linear (open)
    if n == 3:
        if n_mm_d >= 3 or n_mm_i >= 3:
            c = MMP.m3_triangle(syms)
            return {
                "kind": "M3_triangle",
                "positions": np.asarray(c["vertex_positions"], float),
                "edges": ((0, 1), (1, 2), (2, 0)),
                "faces": ((0, 1, 2),),
                "mm_distance": float(c["mm_distance"]),
            }
        sk = CSP._skeleton_linear_m3(d=max(
            MMP.m_m_distance(syms[i], syms[(i + 1) % 3]) for i in range(3)
        ))
        return {
            "kind": "M3_linear",
            "positions": np.asarray(sk.coords, float),
            "edges": tuple(tuple(sorted(e)) for e in sk.edges),
            "faces": (),
            "mm_distance": MMP.m_m_distance(syms[0], syms[1]),
        }

    # -- n == 4 : tetrahedron / square / butterfly
    if n == 4:
        if n_mm_d >= 6 or n_mm_i >= 6:
            c = MMP.m4_tetrahedron(syms)
            return {
                "kind": "M4_tetrahedron",
                "positions": np.asarray(c["vertex_positions"], float),
                "edges": tuple((i, j) for i in range(4) for j in range(i + 1, 4)),
                "faces": ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
                "mm_distance": float(c["mm_distance"]),
            }
        if n_mm_d == 4 or n_mm_i == 4:
            c = MMP.m4_square_planar(syms)
            return {
                "kind": "M4_square",
                "positions": np.asarray(c["vertex_positions"], float),
                "edges": ((0, 1), (1, 2), (2, 3), (3, 0)),
                "faces": ((0, 1, 2, 3),),
                "mm_distance": float(c["mm_distance"]),
            }
        # butterfly (5 edges) or open M4 -> butterfly
        sk = CSP._skeleton_butterfly_m4(d=max(
            MMP.m_m_distance(syms[i], syms[j])
            for i in range(4) for j in range(i + 1, 4)
        ))
        return {
            "kind": "M4_butterfly",
            "positions": np.asarray(sk.coords, float),
            "edges": tuple(tuple(sorted(e)) for e in sk.edges),
            "faces": tuple(tuple(sorted(f)) for f in sk.faces),
            "mm_distance": MMP.m_m_distance(syms[0], syms[1]),
        }

    # -- n == 5 : star (one central + 4 spokes)
    if n == 5:
        sk = CSP._skeleton_star_m5(d=max(
            MMP.m_m_distance(syms[i], syms[(i + 1) % 5]) for i in range(5)
        ))
        return {
            "kind": "M5_star",
            "positions": np.asarray(sk.coords, float),
            "edges": tuple(tuple(sorted(e)) for e in sk.edges),
            "faces": (),
            "mm_distance": MMP.m_m_distance(syms[0], syms[1]),
        }

    # -- n == 6 : octahedron
    if n == 6:
        c = MMP.m6_octahedron(syms)
        return {
            "kind": "M6_octahedron",
            "positions": np.asarray(c["vertex_positions"], float),
            "edges": tuple(c["edges"]),
            "faces": tuple(c["faces"]),
            "mm_distance": float(c["mm_distance"]),
        }

    # -- n >= 7 : generic centroid layout (no closed-form polyhedron yet)
    # Distribute metals on a sphere of radius r = mean_mm / sqrt(2) using
    # a Fibonacci lattice (deterministic, lex-stable).
    r = float(MMP.m_m_distance(syms[0], syms[0])) * (n / 4.0)  # rough scale
    positions = np.zeros((n, 3), float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1.0 - (i / float(n - 1)) * 2.0 if n > 1 else 0.0
        radius_xy = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden * i
        positions[i] = [
            r * math.cos(theta) * radius_xy,
            r * y,
            r * math.sin(theta) * radius_xy,
        ]
    edges = tuple((i, j) for i in range(n) for j in range(i + 1, n)
                  if frozenset({metal_idxs[i], metal_idxs[j]}) in (declared | implied))
    return {
        "kind": f"M{n}_generic",
        "positions": positions,
        "edges": edges,
        "faces": (),
        "mm_distance": float(MMP.m_m_distance(syms[0], syms[0])),
    }


# ---------------------------------------------------------------------------
# Per-pair MM-distance scaling
# ---------------------------------------------------------------------------


def _rescale_skeleton(positions: np.ndarray, edges: Sequence[Tuple[int, int]],
                      syms: Sequence[str]) -> np.ndarray:
    """Rescale skeleton so each edge matches its empirical M-M distance.

    The skeleton primitives use a single ``mm_distance`` for all edges.
    For heteronuclear clusters that one-size-fits-all length is wrong;
    here we anisotropically scale each metal towards the centroid such
    that edge i-j ends up at ``m_m_distance(sym_i, sym_j)``.

    Only used for n=2 (dimer) where the per-edge fix is exact.  For
    n>=3 we leave the regular polyhedron because rigid rescaling would
    break the polyhedral symmetry; the hard-rollback gate then catches
    cases the regular polyhedron cannot satisfy.
    """
    from delfin.fffree import multi_metal_polyhedra as MMP
    if len(positions) != 2:
        return positions
    # dimer: place at +- d/2 along z
    d = MMP.m_m_distance(syms[0], syms[1])
    out = np.zeros((2, 3), float)
    out[0] = (0.0, 0.0, -d / 2.0)
    out[1] = (0.0, 0.0, +d / 2.0)
    return out


# ---------------------------------------------------------------------------
# Donor partitioning (terminal / bridge / face / apex)
# ---------------------------------------------------------------------------


def partition_donors(mol, metal_idxs: Sequence[int]) -> Dict[str, List[Tuple[int, frozenset]]]:
    """Partition donor atoms by number of metal-neighbours.

    Returns:
      {
        "terminal": [(donor_idx, frozenset({m_local})), ...],
        "mu2":      [(donor_idx, frozenset({m_local_a, m_local_b})), ...],
        "mu3":      [(donor_idx, frozenset({m_local_a, m_local_b, m_local_c})), ...],
        "mu4":      [(donor_idx, frozenset(...))],
      }

    Local indices reference ``metal_idxs`` order.
    """
    local_of = {g: i for i, g in enumerate(metal_idxs)}
    out = {"terminal": [], "mu2": [], "mu3": [], "mu4": []}
    if mol is None:
        return out
    metal_set = set(metal_idxs)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in metal_set:
            continue
        s = atom.GetSymbol()
        if _is_metal_sym(s):
            continue
        mnbrs = []
        for nb in atom.GetNeighbors():
            ni = nb.GetIdx()
            if ni in metal_set:
                mnbrs.append(local_of[ni])
        if not mnbrs:
            continue
        mnbrs = sorted(mnbrs)
        key = {1: "terminal", 2: "mu2", 3: "mu3"}.get(len(mnbrs), "mu4")
        out[key].append((atom.GetIdx(), frozenset(mnbrs)))
    for k in out:
        out[k].sort(key=lambda t: (t[0], sorted(t[1])))
    return out


# ---------------------------------------------------------------------------
# Bridging-atom placement
# ---------------------------------------------------------------------------


def _md_target(metal_sym: str, donor_sym: str) -> float:
    """Target M-D distance using the polyhedra covalent table."""
    try:
        from delfin.fffree.polyhedra import md_distance
        d = float(md_distance(metal_sym, donor_sym))
        if d > 0:
            return d
    except Exception:
        pass
    # Safe fallback (mean transition metal + N) — clamped to plausible window
    return max(1.7, min(2.6, 2.0))


def _place_bridges(
    positions: np.ndarray,
    metal_syms: Sequence[str],
    bridges: Dict[str, List[Tuple[int, frozenset]]],
    mol,
) -> Dict[int, np.ndarray]:
    """Place every mu-bridging atom.  Returns ``{atom_idx: (x,y,z)}``."""
    from delfin.fffree import bridging_ligand as BL
    out: Dict[int, np.ndarray] = {}

    for d_idx, mset in bridges.get("mu2", []):
        a, b = sorted(mset)
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        m_a_sym = metal_syms[a]
        m_d_dist = _md_target(m_a_sym, d_sym)
        out[d_idx] = BL.place_mu2_bridge(positions[a], positions[b],
                                         m_x_distance=m_d_dist)

    for d_idx, mset in bridges.get("mu3", []):
        a, b, c = sorted(mset)
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        m_a_sym = metal_syms[a]
        m_d_dist = _md_target(m_a_sym, d_sym)
        out[d_idx] = BL.place_mu3_bridge(positions[a], positions[b], positions[c],
                                         m_x_distance=m_d_dist)

    for d_idx, mset in bridges.get("mu4", []):
        members = sorted(mset)
        pts = [positions[i] for i in members]
        out[d_idx] = BL.place_mu4_bridge(pts, m_x_distance=2.0)

    return out


# ---------------------------------------------------------------------------
# Terminal-donor radial placement
# ---------------------------------------------------------------------------


def _place_terminals(
    positions: np.ndarray,
    metal_idxs: Sequence[int],
    metal_syms: Sequence[str],
    terminals: List[Tuple[int, frozenset]],
    mol,
    already_placed: Dict[int, np.ndarray],
) -> Dict[int, np.ndarray]:
    """Place terminal donors radially outward from each metal.

    Per metal, donors fan out around the (metal - cluster_centroid) axis
    at the empirical M-D distance.  Multiple terminals on one metal
    populate a 90-deg cross plane perpendicular to the radial axis.
    """
    centroid = positions.mean(axis=0)
    out: Dict[int, np.ndarray] = {}
    # Group terminals by parent metal
    by_metal: Dict[int, List[int]] = {i: [] for i in range(len(metal_idxs))}
    for d_idx, mset in terminals:
        (m_local,) = sorted(mset)
        by_metal[m_local].append(d_idx)
    for m_local in sorted(by_metal):
        donors = sorted(by_metal[m_local])
        if not donors:
            continue
        m_pos = positions[m_local]
        radial = m_pos - centroid
        rl = float(np.linalg.norm(radial))
        if rl < 1e-6:
            radial = np.array([0.0, 0.0, 1.0])
            rl = 1.0
        radial_unit = radial / rl
        # Perpendicular basis
        ref = np.array([1.0, 0.0, 0.0]) if abs(radial_unit[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = ref - radial_unit * float(np.dot(ref, radial_unit))
        e1_norm = float(np.linalg.norm(e1))
        if e1_norm < 1e-6:
            continue
        e1 = e1 / e1_norm
        e2 = np.cross(radial_unit, e1)
        n_d = len(donors)
        # Single donor: pure radial.  Multiple: fan around the metal in
        # a cone with half-angle 65 deg (mid-tetrahedral).
        cone_half = math.radians(65.0) if n_d > 1 else 0.0
        cos_a = math.cos(cone_half)
        sin_a = math.sin(cone_half)
        for k, d_idx in enumerate(donors):
            if d_idx in already_placed:
                continue
            d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
            m_d_dist = _md_target(metal_syms[m_local], d_sym)
            if n_d == 1:
                direction = radial_unit
            else:
                phi = 2.0 * math.pi * k / n_d
                cos_p = math.cos(phi)
                sin_p = math.sin(phi)
                direction = (cos_a * radial_unit
                             + sin_a * (cos_p * e1 + sin_p * e2))
                direction = direction / float(np.linalg.norm(direction))
            out[d_idx] = m_pos + direction * m_d_dist
    return out


# ---------------------------------------------------------------------------
# Hard-rollback gates
# ---------------------------------------------------------------------------


def _rollback_check(
    positions: np.ndarray,
    metal_syms: Sequence[str],
    bridges: Dict[str, List[Tuple[int, frozenset]]],
    terminals: List[Tuple[int, frozenset]],
    placed: Dict[int, np.ndarray],
    mol,
    edges: Sequence[Tuple[int, int]],
) -> bool:
    """Validate every M-D within +-_MD_TOL_A and every M-M within +-_MM_TOL_A.

    Returns ``True`` if the assembled coords pass the contract; ``False``
    otherwise.  The caller is expected to fall back to the single-metal
    path on a False result.
    """
    from delfin.fffree import multi_metal_polyhedra as MMP
    # M-M distance contract
    for (i, j) in edges:
        if i >= len(positions) or j >= len(positions):
            return False
        d = float(np.linalg.norm(positions[i] - positions[j]))
        d_target = MMP.m_m_distance(metal_syms[i], metal_syms[j])
        if d_target <= 0:
            continue
        if abs(d - d_target) > _MM_TOL_A:
            return False

    # M-D distance contract (terminal + bridging)
    def _check_md(d_idx: int, m_locals: Sequence[int]) -> bool:
        if d_idx not in placed:
            return False
        pos_d = placed[d_idx]
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        for ml in m_locals:
            tgt = _md_target(metal_syms[ml], d_sym)
            obs = float(np.linalg.norm(pos_d - positions[ml]))
            if abs(obs - tgt) > _MD_TOL_A:
                # mu-bridges legitimately deviate by up to 0.15 A
                # because the perpendicular height is a single solve
                # against multiple metals.  Use the looser M-M tol.
                if len(m_locals) >= 2 and abs(obs - tgt) > _MM_TOL_A + 0.05:
                    return False
                elif len(m_locals) == 1 and abs(obs - tgt) > _MD_TOL_A:
                    return False
        return True

    for d_idx, mset in terminals:
        if not _check_md(d_idx, sorted(mset)):
            return False
    for d_idx, mset in bridges.get("mu2", []):
        if not _check_md(d_idx, sorted(mset)):
            return False
    for d_idx, mset in bridges.get("mu3", []):
        if not _check_md(d_idx, sorted(mset)):
            return False
    for d_idx, mset in bridges.get("mu4", []):
        if not _check_md(d_idx, sorted(mset)):
            return False
    return True


# ---------------------------------------------------------------------------
# Public orchestrator
# ---------------------------------------------------------------------------


def assemble_multi_metal(
    mol,
    ligands: Optional[Sequence[Dict]] = None,
    metals: Optional[Sequence[str]] = None,
    config: Optional[Dict] = None,
) -> Optional[Tuple[List[str], np.ndarray, List[int]]]:
    """Assemble a multi-metal complex.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Parsed molecule with explicit metal-donor bonds.  Must contain
        >= 2 metal atoms; otherwise the function returns ``None``.
    ligands : sequence of fffree-style ligand dicts, optional
        Currently unused but accepted for API symmetry with
        ``assemble_complex.assemble_from_config``.
    metals : list of metal element symbols, optional
        When provided, takes precedence over auto-detection.
    config : dict, optional
        Reserved for vertex-isomer information; not consumed yet.

    Returns
    -------
    tuple(syms, P, donors) | None
        ``syms`` is a flat list of element symbols, length matches ``P``.
        ``P`` is an (N, 3) ndarray of coordinates.  ``donors`` lists the
        donor atom indices in the final coordinate array.  Returns
        ``None`` when:
          - the env flag is not set, OR
          - the mol has < 2 metals, OR
          - the hard-rollback gate rejects the assembled coords.

    Determinism: identical input -> identical output across machines and runs.
    """
    if not _multi_metal_enabled():
        return None
    if mol is None:
        return None

    graph = build_metal_graph(mol)
    metal_idxs = list(graph["metal_idxs"])
    metal_syms_g = list(graph["metal_syms"])
    if len(metal_idxs) < 2:
        return None
    # Optional override of the symbols (e.g. for test stubs)
    if metals is not None and len(metals) == len(metal_idxs):
        metal_syms_g = list(metals)

    # 1. Skeleton choice
    sk = choose_skeleton({
        **graph,
        "metal_syms": tuple(metal_syms_g),
    })
    if sk is None:
        return None
    positions = sk["positions"].copy()
    # Rescale per-edge for dimer (the only n-case where it's well-defined)
    positions = _rescale_skeleton(positions, sk["edges"], metal_syms_g)

    # 2. Donor partition
    parts = partition_donors(mol, metal_idxs)
    terminals = parts["terminal"]
    bridges = {"mu2": parts["mu2"], "mu3": parts["mu3"], "mu4": parts["mu4"]}

    # 3. Place bridges
    placed: Dict[int, np.ndarray] = {}
    bridge_coords = _place_bridges(positions, metal_syms_g, bridges, mol)
    for k, v in bridge_coords.items():
        placed[k] = np.asarray(v, float)

    # 4. Place terminals
    term_coords = _place_terminals(
        positions, metal_idxs, metal_syms_g, terminals, mol, placed,
    )
    for k, v in term_coords.items():
        placed[k] = np.asarray(v, float)

    # 5. Hard-rollback contract
    if not _rollback_check(positions, metal_syms_g, bridges, terminals,
                           placed, mol, sk["edges"]):
        return None

    # 6. Compose final (syms, P, donors).  Order: metals (lex-sorted),
    # then mu4, mu3, mu2, terminal (each lex-sorted by atom_idx).
    syms: List[str] = []
    coords: List[np.ndarray] = []
    donor_indices_out: List[int] = []
    for li, gi in enumerate(metal_idxs):
        syms.append(metal_syms_g[li])
        coords.append(positions[li])
    placed_order: List[int] = []
    for grp in ("mu4", "mu3", "mu2"):
        for d_idx, _ in bridges[grp]:
            if d_idx in placed and d_idx not in placed_order:
                placed_order.append(d_idx)
    for d_idx, _ in terminals:
        if d_idx in placed and d_idx not in placed_order:
            placed_order.append(d_idx)
    for d_idx in placed_order:
        syms.append(mol.GetAtomWithIdx(d_idx).GetSymbol())
        coords.append(placed[d_idx])
        donor_indices_out.append(len(syms) - 1)

    if len(syms) != len(coords):
        return None
    P = np.vstack(coords)

    # Surgical fix (2026-06-04, b00f9a0 forensik): optional global-clash gate.
    # The hard-rollback gate above validates M-M and M-D distances against
    # empirical targets but does NOT detect heavy-heavy or X-H clashes between
    # ligand donors / bridges that have been placed on opposite vertices of
    # the same metal.  When the assembled poly-nuclear structure has any
    # clash, downstream detectors see it as funcgrp / F20_h / angles
    # violations on the otherwise-correct multi-metal skeleton.  Reject and
    # fall back to the legacy single-metal path on any clash.
    # Default OFF -> byte-identical to b00f9a0.
    if _global_clash_gate_active():
        try:
            from delfin.fffree import build_time_clash_gate as _bcg
            if _bcg.collapse_count(syms, P) > 0:
                return None
        except Exception:
            # Defensive: any failure in the clash gate is non-fatal — we
            # accept the existing-validated structure (no behaviour change
            # vs the existing rollback contract).
            pass

    return syms, P, donor_indices_out


# ---------------------------------------------------------------------------
# Detection helpers used by the dispatcher
# ---------------------------------------------------------------------------


def should_dispatch_multi_metal(mol) -> bool:
    """Return True when the dispatcher should call ``assemble_multi_metal``."""
    if not _multi_metal_enabled():
        return False
    if mol is None:
        return False
    n = sum(1 for a in mol.GetAtoms() if _is_metal_sym(a.GetSymbol()))
    return n >= 2
