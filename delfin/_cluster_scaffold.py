"""Two-level scaffold construction for ``cluster`` class TMC SMILES.

Status: PROTOTYPE (T7.3 design — Welle-3, 2026-05-15).  NOT wired into
``smiles_converter.smiles_to_xyz_isomers`` yet.  This module documents
how a future ``cluster`` class would be detected, how its metal-skeleton
would be enumerated, and how ligands would be placed relative to the
skeleton.  Designed to live alongside the existing
``_build_multimetal_scaffold`` (dimer-only) and the per-metal Pólya path.

The module is environment-flag-gated when wired (``DELFIN_CLUSTER_PATH=1``);
no runtime behaviour is changed by *importing* this file.

Author: T7.3 Welle-3 agent (chirality of the metal core, Polya-on-skeleton,
ligand placement relative to face-centroids).  See companion report
``iters/welle3_T7.3_cluster_class_scaffold_20260515.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Classifier extension — call from smiles_converter._classify_complex_class
# ---------------------------------------------------------------------------

def is_cluster_complex(mol, metal_set: Iterable[str]) -> bool:
    """Return True when *mol* qualifies as a metal cluster.

    Heuristic (chemistry-feature-based, no SMILES literals):
      1. n_metals >= 3                   (M3/M4/M6 cores)
      2. n_mm_edges >= 2                 (cluster connectivity)
      3. n_mm_edges / max(1, n_ml_edges) >= 0.10
         (M-M backbone is non-trivial relative to M-L decoration)

    Tier-relaxed: also accepts (n_metals == 2 AND n_mm_edges >= 1) when
    the caller signals dimer-clusters should fold into the cluster path.
    Default below is *strict* (≥3 metals); dimers continue through the
    existing multi_sigma / multi_hapto pipelines.
    """
    if mol is None:
        return False
    metal_set = set(metal_set)
    n_metals = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in metal_set)
    if n_metals < 3:
        return False
    mm_edges = 0
    ml_edges = 0
    for bond in mol.GetBonds():
        s1 = bond.GetBeginAtom().GetSymbol()
        s2 = bond.GetEndAtom().GetSymbol()
        m1 = s1 in metal_set
        m2 = s2 in metal_set
        if m1 and m2:
            mm_edges += 1
        elif m1 or m2:
            ml_edges += 1
    if mm_edges < 2:
        return False
    ratio = mm_edges / max(1, ml_edges)
    return ratio >= 0.10


# ---------------------------------------------------------------------------
# Metal-skeleton enumeration (level 1)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MetalSkeleton:
    """A 3D arrangement of metal atoms forming the cluster backbone.

    Coordinates are in a *local* skeleton frame; ligands are placed in
    that frame, then the whole assembly is rigid-Procrustes-aligned to
    the ETKDG template (if available).
    """
    name: str
    coords: Tuple[Tuple[float, float, float], ...]      # one per metal
    edges: FrozenSet[FrozenSet[int]]                     # M-M edges (local indices)
    faces: Tuple[FrozenSet[int], ...]                    # 3+ atom faces

    @property
    def n_metals(self) -> int:
        return len(self.coords)


def _skeleton_triangle(d: float = 2.7) -> MetalSkeleton:
    """Equilateral M3 triangle (e.g. Os3(CO)12, Ru3, Rh3)."""
    h = d * math.sqrt(3.0) / 2.0
    return MetalSkeleton(
        name="triangle",
        coords=((0.0, 0.0, 0.0), (d, 0.0, 0.0), (d / 2.0, h, 0.0)),
        edges=frozenset({frozenset({0, 1}), frozenset({1, 2}), frozenset({0, 2})}),
        faces=(frozenset({0, 1, 2}),),
    )


def _skeleton_linear_m3(d: float = 2.7) -> MetalSkeleton:
    """Open M3 linear chain (M-M-M, fewer M-M bonds, e.g. mixed-valence chain)."""
    return MetalSkeleton(
        name="linear_m3",
        coords=((0.0, 0.0, 0.0), (d, 0.0, 0.0), (2.0 * d, 0.0, 0.0)),
        edges=frozenset({frozenset({0, 1}), frozenset({1, 2})}),
        faces=(),
    )


def _skeleton_tetrahedron(d: float = 2.7) -> MetalSkeleton:
    """Regular M4 tetrahedron (e.g. Ir4(CO)12, Rh4, Co4)."""
    # Standard tetrahedron with edge length d.
    a = d / math.sqrt(2.0)
    coords = (
        (a, a, a), (a, -a, -a), (-a, a, -a), (-a, -a, a),
    )
    edges = frozenset(
        frozenset({i, j}) for i in range(4) for j in range(i + 1, 4)
    )
    faces = tuple(
        frozenset({i, j, k})
        for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)
    )
    return MetalSkeleton(name="tetrahedron", coords=coords, edges=edges, faces=faces)


def _skeleton_butterfly_m4(d: float = 2.7) -> MetalSkeleton:
    """M4 "butterfly" — tetrahedron minus one edge (5 M-M edges)."""
    a = d / math.sqrt(2.0)
    coords = ((a, a, a), (a, -a, -a), (-a, a, -a), (-a, -a, a))
    # Remove the (2,3) edge to open the butterfly hinge.
    edges = frozenset({
        frozenset({0, 1}), frozenset({0, 2}), frozenset({0, 3}),
        frozenset({1, 2}), frozenset({1, 3}),
    })
    faces = (frozenset({0, 1, 2}), frozenset({0, 1, 3}))
    return MetalSkeleton(name="butterfly_m4", coords=coords, edges=edges, faces=faces)


def _skeleton_octahedron_m6(d: float = 2.7) -> MetalSkeleton:
    """Regular M6 octahedron (e.g. Mo6, [Mo6Cl14]2-, Zr6 cluster cores)."""
    r = d / math.sqrt(2.0)
    coords = (
        (r, 0, 0), (-r, 0, 0), (0, r, 0), (0, -r, 0), (0, 0, r), (0, 0, -r),
    )
    edges = frozenset(
        frozenset({i, j})
        for i in range(6) for j in range(i + 1, 6)
        if abs(sum((coords[i][k] - coords[j][k]) ** 2 for k in range(3)) - d * d) < 1e-3
    )
    faces = tuple(
        frozenset({i, j, k})
        for i in range(6) for j in range(i + 1, 6) for k in range(j + 1, 6)
        if frozenset({i, j}) in edges
        and frozenset({i, k}) in edges
        and frozenset({j, k}) in edges
    )
    return MetalSkeleton(name="octahedron_m6", coords=coords, edges=edges, faces=faces)


def _skeleton_star_m5(d: float = 2.7) -> MetalSkeleton:
    """Centro-star M5 (one central metal + 4 peripheral spokes).

    Models the BOQYIR pattern: Co at center with three Sn(Me)3 spokes + Cp anchor.
    Generalises Os(Sn)5(Cl) (MASNCL) when n_peripheral=5 (octahedral spokes).
    """
    coords = (
        (0.0, 0.0, 0.0),
        (d, 0.0, 0.0), (-d, 0.0, 0.0),
        (0.0, d, 0.0), (0.0, -d, 0.0),
    )
    edges = frozenset({frozenset({0, i}) for i in range(1, 5)})
    return MetalSkeleton(name="star_m5", coords=coords, edges=edges, faces=())


def enumerate_skeletons_for_mol(mol, metal_set: Iterable[str]) -> List[MetalSkeleton]:
    """Choose plausible skeletons based on (n_metals, n_mm_edges) signature.

    Strategy: pick skeletons whose connectivity matches the SMILES-declared
    M-M-edge count exactly (canonical match) OR matches within ±1 (relaxed,
    for SMILES that under-specifies metal-metal bonds).
    """
    metal_set = set(metal_set)
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    n_metals = len(metals)
    metal_idx_to_local = {midx: i for i, midx in enumerate(metals)}
    declared_mm_edges = 0
    declared_pairs: set = set()
    for bond in mol.GetBonds():
        s1 = bond.GetBeginAtom().GetSymbol()
        s2 = bond.GetEndAtom().GetSymbol()
        if s1 in metal_set and s2 in metal_set:
            declared_mm_edges += 1
            declared_pairs.add(frozenset({
                metal_idx_to_local[bond.GetBeginAtomIdx()],
                metal_idx_to_local[bond.GetEndAtomIdx()],
            }))

    candidates: List[MetalSkeleton] = []
    if n_metals == 3:
        if declared_mm_edges >= 3:
            candidates.append(_skeleton_triangle())
        if declared_mm_edges <= 2:
            candidates.append(_skeleton_linear_m3())
        # Star-M3 = one center + 2 spokes (peripheral connectivity).
        # Reuse star_m5 truncated would need its own helper — omitted in PoC.
    elif n_metals == 4:
        if declared_mm_edges >= 6:
            candidates.append(_skeleton_tetrahedron())
        if 4 <= declared_mm_edges <= 5:
            candidates.append(_skeleton_butterfly_m4())
        if declared_mm_edges == 3:
            # M4 with central spoke = star-M4 — model as butterfly opened.
            candidates.append(_skeleton_butterfly_m4())
        if declared_mm_edges <= 3:
            # Central-metal star (BOQYIR-like Co + 3 Sn): use star_m5 truncated;
            # fall back to butterfly for prototype simplicity.
            candidates.append(_skeleton_butterfly_m4())
    elif n_metals == 5:
        candidates.append(_skeleton_star_m5())
    elif n_metals == 6:
        if declared_mm_edges >= 6:
            candidates.append(_skeleton_octahedron_m6())
        else:
            # Os(SnR3)5 = central Os + 5 Sn-spokes. Star with extra Cl/anchor:
            # use star_m5 + 1 extra peripheral metal (PoC: still star_m5).
            candidates.append(_skeleton_star_m5())

    # Deduplicate by name.
    seen: set = set()
    unique: List[MetalSkeleton] = []
    for sk in candidates:
        if sk.name not in seen:
            seen.add(sk.name)
            unique.append(sk)
    return unique


# ---------------------------------------------------------------------------
# Two-level construction — wire-in stub (NOT YET CALLED from converter)
# ---------------------------------------------------------------------------

@dataclass
class ClusterBuildPlan:
    """One concrete attempt to build a cluster: skeleton + ligand placement."""
    skeleton: MetalSkeleton
    metal_indices: Tuple[int, ...]            # SMILES atom indices for skeleton metals
    terminal_donors: Tuple[Tuple[int, int], ...]      # (donor_idx, parent_metal_local)
    edge_bridges: Tuple[Tuple[int, FrozenSet[int]], ...]   # (donor_idx, {m_local_a, m_local_b})
    face_caps: Tuple[Tuple[int, FrozenSet[int]], ...]      # (donor_idx, {m1,m2,m3})


def classify_donors_for_skeleton(
    mol,
    metal_indices: Sequence[int],
    skeleton: MetalSkeleton,
    metal_set: Iterable[str],
) -> ClusterBuildPlan:
    """Partition non-metal donor atoms into terminal/edge-bridge/face-cap.

    A donor d is:
      - face-cap if bonded to ≥3 of the skeleton metals
      - edge-bridge if bonded to exactly 2 skeleton metals
      - terminal if bonded to exactly 1 skeleton metal
    """
    metal_set = set(metal_set)
    local_of = {midx: i for i, midx in enumerate(metal_indices)}
    terminals: List[Tuple[int, int]] = []
    edge_bridges: List[Tuple[int, FrozenSet[int]]] = []
    face_caps: List[Tuple[int, FrozenSet[int]]] = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metal_set:
            continue
        parent_locals = []
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() in metal_set and nbr.GetIdx() in local_of:
                parent_locals.append(local_of[nbr.GetIdx()])
        if not parent_locals:
            continue
        if len(parent_locals) == 1:
            terminals.append((atom.GetIdx(), parent_locals[0]))
        elif len(parent_locals) == 2:
            edge_bridges.append((atom.GetIdx(), frozenset(parent_locals)))
        else:
            face_caps.append((atom.GetIdx(), frozenset(parent_locals)))

    return ClusterBuildPlan(
        skeleton=skeleton,
        metal_indices=tuple(metal_indices),
        terminal_donors=tuple(terminals),
        edge_bridges=tuple(edge_bridges),
        face_caps=tuple(face_caps),
    )


def place_skeleton_with_ligands(
    plan: ClusterBuildPlan,
    mol,
    metal_set: Iterable[str],
    ml_bond_lookup,
) -> Dict[int, Tuple[float, float, float]]:
    """Produce initial coordinates {atom_idx: (x,y,z)} for skeleton metals
    and their first-shell donors.  Heavy ligand atoms (beyond the donor)
    are *not* placed here — the downstream per-metal fragment-docker
    handles them, exactly like the existing dimer pathway.
    """
    coords: Dict[int, Tuple[float, float, float]] = {}
    sk = plan.skeleton

    # 1) Metals
    for m_local, m_idx in enumerate(plan.metal_indices):
        coords[m_idx] = sk.coords[m_local]

    # 2) Edge bridges: midpoint of the two metals, lifted off the M-M line
    #    by a fraction of an ideal M-D bond so the bridge angle is ~120-130°.
    for d_idx, m_set in plan.edge_bridges:
        m_a, m_b = sorted(m_set)
        pa = sk.coords[m_a]
        pb = sk.coords[m_b]
        mid = tuple((pa[k] + pb[k]) / 2.0 for k in range(3))
        # Lift in +z (skeleton local frame) — alternating sign per bridge
        # to spread bridges around the polyhedron.
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        m_sym = mol.GetAtomWithIdx(plan.metal_indices[m_a]).GetSymbol()
        d_md = float(ml_bond_lookup(m_sym, d_sym))
        # Cosine of bridge angle 130° = -0.643.  M-D-M -> D sits at
        # height h such that |D - pa| = d_md and (D - pa).(pb - pa) is
        # consistent with the chosen angle.  Simple lift: h = d_md * 0.55.
        h = d_md * 0.55
        # Pick a stable lift axis: cross of (pb - pa) with z if available.
        ab = tuple(pb[k] - pa[k] for k in range(3))
        ab_len = math.sqrt(sum(c * c for c in ab))
        # default lift = global +z; if ab is parallel to z, use +y.
        if ab_len > 1e-6 and abs(ab[2]) / ab_len < 0.9:
            lift = (0.0, 0.0, h)
        else:
            lift = (0.0, h, 0.0)
        coords[d_idx] = tuple(mid[k] + lift[k] for k in range(3))

    # 3) Face caps: centroid of the face, lifted off the face plane.
    for d_idx, m_set in plan.face_caps:
        members = sorted(m_set)
        pts = [sk.coords[m] for m in members]
        centroid = tuple(sum(p[k] for p in pts) / len(pts) for k in range(3))
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        m_sym = mol.GetAtomWithIdx(plan.metal_indices[members[0]]).GetSymbol()
        d_md = float(ml_bond_lookup(m_sym, d_sym))
        # Lift along face normal: cross of (p1-p0) x (p2-p0).
        if len(pts) >= 3:
            a = tuple(pts[1][k] - pts[0][k] for k in range(3))
            b = tuple(pts[2][k] - pts[0][k] for k in range(3))
            n = (
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            )
            nl = math.sqrt(sum(c * c for c in n))
            if nl > 1e-6:
                lift = tuple(c / nl * d_md * 0.85 for c in n)
            else:
                lift = (0.0, 0.0, d_md * 0.85)
        else:
            lift = (0.0, 0.0, d_md * 0.85)
        coords[d_idx] = tuple(centroid[k] + lift[k] for k in range(3))

    # 4) Terminal donors: place radially outward from parent metal.
    #    Direction = (parent - cluster_centroid), normalised.  This
    #    ensures terminals point AWAY from the M-M backbone — the
    #    most important rule for clusters (cf. Os3(CO)12 axial vs
    #    radial CO orientation).  Multiple terminals on the same
    #    metal are tilted by ±15° around the radial axis (PoC keeps
    #    them coincident; the downstream per-metal Pólya then enumerates
    #    the orientation degree of freedom).
    centroid_all = tuple(
        sum(sk.coords[i][k] for i in range(sk.n_metals)) / sk.n_metals
        for k in range(3)
    )
    for d_idx, m_local in plan.terminal_donors:
        pm = sk.coords[m_local]
        radial = tuple(pm[k] - centroid_all[k] for k in range(3))
        rl = math.sqrt(sum(c * c for c in radial))
        if rl < 1e-6:
            # metal AT centroid (e.g. star_m5 center) — pick +z
            radial = (0.0, 0.0, 1.0)
            rl = 1.0
        d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
        m_sym = mol.GetAtomWithIdx(plan.metal_indices[m_local]).GetSymbol()
        d_md = float(ml_bond_lookup(m_sym, d_sym))
        scale = d_md / rl
        coords[d_idx] = tuple(pm[k] + radial[k] * scale for k in range(3))

    return coords


# ---------------------------------------------------------------------------
# Self-test / PoC entry point
# ---------------------------------------------------------------------------

def proof_of_concept_run(mol, metal_set: Iterable[str], ml_bond_lookup) -> Optional[
    Tuple[ClusterBuildPlan, Dict[int, Tuple[float, float, float]]]
]:
    """One-shot PoC: classify, pick first skeleton, build coords.

    Returns ``None`` if *mol* is not a cluster, otherwise the plan +
    initial coordinates.  Used only by the determinism/PoC test in
    the agent report — *NOT* called from production code.
    """
    if not is_cluster_complex(mol, metal_set):
        return None
    skeletons = enumerate_skeletons_for_mol(mol, metal_set)
    if not skeletons:
        return None
    metal_set = set(metal_set)
    metal_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    skeleton = skeletons[0]
    # Match skeleton size to actual metal count (PoC slicing for partial fits)
    if len(metal_indices) < skeleton.n_metals:
        return None
    metals_used = metal_indices[: skeleton.n_metals]
    plan = classify_donors_for_skeleton(mol, metals_used, skeleton, metal_set)
    coords = place_skeleton_with_ligands(plan, mol, metal_set, ml_bond_lookup)
    return plan, coords


# ---------------------------------------------------------------------------
# Welle-5f-G — multi_hapto skeleton-first scaffold (2-level construction)
# ---------------------------------------------------------------------------

def _hapto_groups_for_metal(mol, metal_idx: int, metal_set: Iterable[str]):
    """Return list[list[int]] of hapto-carbon groups bonded to ``metal_idx``.

    A hapto group = set of carbon neighbours of ``metal_idx`` that are
    (a) mutually bonded OR (b) co-members of the same SSSR ring.
    This handles η4/η5/η6 ligands where the metal-bound carbons span
    a CC=CC pattern with two interior sp3 carbons (e.g. cyclohexadienyl)
    that breaks naive mutual-bond BFS into multiple sub-fragments.
    Mirrors ``_find_hapto_groups`` intent but constrained to one metal.
    """
    metal_set = set(metal_set)
    try:
        atom = mol.GetAtomWithIdx(metal_idx)
    except Exception:
        return []
    if atom.GetSymbol() not in metal_set:
        return []
    c_neighbours = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(c_neighbours) < 2:
        return []
    c_set = set(c_neighbours)
    # Step 1: BFS-by-direct-bond.
    raw_groups: List[List[int]] = []
    seen: set = set()
    for start in c_neighbours:
        if start in seen:
            continue
        comp: List[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            try:
                cur_atom = mol.GetAtomWithIdx(cur)
            except Exception:
                continue
            for nb in cur_atom.GetNeighbors():
                ni = nb.GetIdx()
                if ni in c_set and ni not in seen:
                    seen.add(ni)
                    stack.append(ni)
        if len(comp) >= 1:
            raw_groups.append(sorted(comp))
    if not raw_groups:
        return []
    # Step 2: merge groups whose atoms share an SSSR ring (cyclopentadienyl
    # with two interior sp3 carbons gets split by BFS into 2 sub-groups,
    # but both belong to the same 5-ring → merge into one η-ligand).
    # GetSymmSSSR works on unsanitized RDKit mols where AtomRings() would
    # otherwise be empty.
    try:
        from rdkit import Chem as _Chem
        _Chem.GetSymmSSSR(mol)
        ring_info = mol.GetRingInfo()
        rings = [set(r) for r in ring_info.AtomRings()]
    except Exception:
        rings = []
    # Union-find over raw groups using "share-a-ring" relation.
    parent = list(range(len(raw_groups)))
    def _find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def _union(i, j):
        ri, rj = _find(i), _find(j)
        if ri != rj:
            parent[ri] = rj
    for ring in rings:
        members = [
            i for i, g in enumerate(raw_groups)
            if any(c in ring for c in g)
        ]
        for k in range(1, len(members)):
            _union(members[0], members[k])
    merged: Dict[int, List[int]] = {}
    for i, g in enumerate(raw_groups):
        r = _find(i)
        merged.setdefault(r, []).extend(g)
    groups = [sorted(set(g)) for g in merged.values() if len(set(g)) >= 2]
    # Deterministic order: by sorted first-atom index.
    groups.sort(key=lambda g: g[0])
    return groups


def _ring_centroid_distance(metal_sym: str, eta: int) -> float:
    """Approximate metal-centroid distance for an eta-n ring.

    Empirical defaults: η5 ~ 1.80 A (Cp), η6 ~ 1.65 A (arene), η4 ~ 1.95 A (diene),
    η3 ~ 2.00 A (allyl), η2 ~ 2.05 A (alkene).  Heavier metals get +0.10 A.
    Used only when no chemistry-specific lookup table is available.
    """
    base = {2: 2.05, 3: 2.00, 4: 1.95, 5: 1.80, 6: 1.65}.get(eta, 1.80)
    heavy_d = {"Ru", "Rh", "Os", "Ir", "Pt", "Re", "W", "Mo", "Ta", "Hf"}
    if metal_sym in heavy_d:
        return base + 0.10
    return base


def is_multihapto_skeleton_candidate(mol, metal_set: Iterable[str]) -> bool:
    """Return True when *mol* is a multi-metal complex with at least one
    hapto-coordinated metal — the target chemistry for the multi_hapto
    skeleton-first scaffold (Welle-5f-G).
    """
    if mol is None:
        return False
    metal_set = set(metal_set)
    metal_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    if len(metal_idxs) < 2:
        return False
    for mi in metal_idxs:
        if _hapto_groups_for_metal(mol, mi, metal_set):
            return True
    return False


def build_multihapto_scaffold(
    mol,
    metal_set: Iterable[str],
    ml_bond_lookup,
    mm_bond_lookup=None,
) -> Optional[Dict[int, Tuple[float, float, float]]]:
    """Two-level scaffold for multi_hapto: metal-skeleton first, hapto-rings
    + σ-donors placed relative to the skeleton frame.

    Strategy
    --------
    1. **Metal skeleton** (level 1): place hapto-metal M_h at origin, every
       other metal at ideal M-M distance (``mm_bond_lookup`` lookup, fall
       back to ``r_cov(M1) + r_cov(M2) + 0.4``).  Layout: M_h at origin,
       subsequent metals along the +x axis spaced by ideal distance.
       For 3+ metals connected to M_h, place them on a circle in the xy
       plane around M_h at ideal distance, evenly distributed.  This
       avoids the dimer-only ``_build_multimetal_scaffold`` hard-gate
       and the cluster-only ``place_skeleton_with_ligands`` gate.

    2. **Hapto rings** (level 2): for each hapto-coordinated metal, place
       the ring centroid along -z (away from skeleton plane), with the
       ring atoms in a regular polygon perpendicular to -z at
       ``_ring_centroid_distance(M, eta)`` from the metal.  Ring atoms are
       returned in ``coords`` keyed by their RDKit atom indices so the
       downstream rigid-Procrustes can dock the full ligand using the
       ring as anchor.

    3. **σ-donors** (level 2): non-hapto, non-metal neighbours of each
       metal are placed radially outward — direction = (metal_pos −
       cluster_centroid) normalised, distance = ``_get_ml_bond_length``.
       This is exactly the cluster-class terminal placement (so the M-D
       σ-bonds land within the post-build validator's [0.85, 1.10] × ideal
       window).  When multiple σ-donors share a metal, they fan out
       around the radial axis (90° increments).

    Returns ``{atom_idx: (x, y, z)}`` for the metal-skeleton + first donor
    shell (hapto rings + σ-donors).  Heavy ligand-arm atoms are *not*
    placed here — the downstream fragment-Procrustes handles them.

    Default-OFF wire-in: the caller gates this behind
    ``DELFIN_5F_G_MULTIHAPTO_3PLUS_METALS=1``.  When unset the legacy
    pipeline (dimer scaffold for n_metals==2, per-metal Pólya for >=3)
    runs unchanged.
    """
    if mol is None:
        return None
    metal_set = set(metal_set)
    metal_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    n_m = len(metal_idxs)
    if n_m < 2:
        return None

    # Pick the hapto-coordinated metal as the skeleton origin.  When more
    # than one metal carries a hapto group, prefer the one with the LARGEST
    # ring (Cp >> arene by donor-richness) — empirically the hapto-metal
    # frame is the most rigid and least likely to drift in UFF.
    hapto_info: Dict[int, List[List[int]]] = {}
    for mi in metal_idxs:
        gs = _hapto_groups_for_metal(mol, mi, metal_set)
        if gs:
            hapto_info[mi] = gs
    if not hapto_info:
        return None
    # Sort by (max group size desc, atom idx asc) for determinism.
    sorted_hapto = sorted(
        hapto_info.items(),
        key=lambda kv: (-max(len(g) for g in kv[1]), kv[0]),
    )
    m_h = sorted_hapto[0][0]
    other_metals = [mi for mi in metal_idxs if mi != m_h]

    coords: Dict[int, Tuple[float, float, float]] = {}
    coords[m_h] = (0.0, 0.0, 0.0)
    m_h_sym = mol.GetAtomWithIdx(m_h).GetSymbol()

    # --- Level 1: secondary metals on a +xy circle at ideal M-M distance.
    # For n_metals==2: just place M2 along +x.
    # For n_metals>=3: evenly distribute on a circle in the xy-plane around M_h.
    if len(other_metals) == 1:
        m_o = other_metals[0]
        m_o_sym = mol.GetAtomWithIdx(m_o).GetSymbol()
        target_mm = None
        if mm_bond_lookup is not None:
            try:
                target_mm = float(mm_bond_lookup(m_h_sym, m_o_sym))
            except Exception:
                target_mm = None
        if target_mm is None or target_mm <= 0:
            # Heuristic fallback: 2 × covalent-radius sum + 0.4 A buffer.
            # Use ml_bond_lookup as a proxy when available (M-D bond length
            # roughly tracks atomic size).
            try:
                d_mh_h = float(ml_bond_lookup(m_h_sym, "H"))
                d_mo_h = float(ml_bond_lookup(m_o_sym, "H"))
                target_mm = (d_mh_h + d_mo_h) + 0.4
            except Exception:
                target_mm = 3.0
        coords[m_o] = (target_mm, 0.0, 0.0)
    else:
        # 3+ metals: M_h at origin, others on circle in xy-plane.
        # Radius = ideal M_h-M_other distance (use first as anchor; assume
        # all secondary metals are roughly chemically equivalent — a
        # reasonable approximation for star-cluster patterns).
        n_o = len(other_metals)
        for k, m_o in enumerate(other_metals):
            m_o_sym = mol.GetAtomWithIdx(m_o).GetSymbol()
            target_mm = None
            if mm_bond_lookup is not None:
                try:
                    target_mm = float(mm_bond_lookup(m_h_sym, m_o_sym))
                except Exception:
                    target_mm = None
            if target_mm is None or target_mm <= 0:
                try:
                    d_mh_h = float(ml_bond_lookup(m_h_sym, "H"))
                    d_mo_h = float(ml_bond_lookup(m_o_sym, "H"))
                    target_mm = (d_mh_h + d_mo_h) + 0.4
                except Exception:
                    target_mm = 3.0
            angle = 2.0 * math.pi * k / n_o
            coords[m_o] = (
                target_mm * math.cos(angle),
                target_mm * math.sin(angle),
                0.0,
            )

    # --- Level 2a: hapto rings.  For each hapto-coordinated metal, place
    # the ring centroid along -z (out of the metal-skeleton plane) at
    # ``_ring_centroid_distance(M, eta)``, and the ring atoms in a regular
    # polygon perpendicular to -z.  This guarantees the η-cage is NOT
    # coplanar with another η-cage (the welle3 T7.2 Root Cause #2 = Cp
    # ring collapse).
    cluster_centroid = (
        sum(coords[mi][0] for mi in metal_idxs) / n_m,
        sum(coords[mi][1] for mi in metal_idxs) / n_m,
        sum(coords[mi][2] for mi in metal_idxs) / n_m,
    )
    hapto_atoms_placed: set = set()
    for mi, groups in hapto_info.items():
        m_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        m_pos = coords[mi]
        # Multiple hapto groups per metal: alternate +z and -z so two
        # rings on the same metal end up on opposite faces (sandwich
        # complex pattern, e.g. ferrocene-like Cp-M-Cp).
        for g_idx, ring in enumerate(groups):
            eta = len(ring)
            d_centroid = _ring_centroid_distance(m_sym, eta)
            # Ring plane normal: by default -z (so the ring sits below
            # the metal-skeleton plane).  When the metal is at origin,
            # also tilt slightly outward (radial component) so two
            # adjacent hapto-metals don't sandwich into the SAME ring.
            radial = (
                m_pos[0] - cluster_centroid[0],
                m_pos[1] - cluster_centroid[1],
                m_pos[2] - cluster_centroid[2],
            )
            r_len = math.sqrt(sum(c * c for c in radial))
            if r_len < 1e-6:
                radial_unit = (0.0, 0.0, 0.0)
            else:
                radial_unit = (radial[0] / r_len, radial[1] / r_len, radial[2] / r_len)
            # Centroid direction: dominantly -z, with a 30% radial tilt
            # away from cluster centroid.  Alternate sign per group so
            # two rings on the same metal end on opposite faces.
            z_sign = -1.0 if g_idx % 2 == 0 else 1.0
            tilt = 0.3
            cn_dir = (
                tilt * radial_unit[0],
                tilt * radial_unit[1],
                z_sign * (1.0 - tilt),
            )
            cn_norm = math.sqrt(sum(c * c for c in cn_dir))
            if cn_norm < 1e-6:
                cn_unit = (0.0, 0.0, z_sign)
            else:
                cn_unit = (cn_dir[0] / cn_norm, cn_dir[1] / cn_norm, cn_dir[2] / cn_norm)
            centroid = (
                m_pos[0] + cn_unit[0] * d_centroid,
                m_pos[1] + cn_unit[1] * d_centroid,
                m_pos[2] + cn_unit[2] * d_centroid,
            )
            # Build two orthonormal vectors spanning the plane perpendicular
            # to cn_unit, then place ring atoms in a regular polygon.
            # Vector e1: any vector not parallel to cn_unit, projected onto
            # the perpendicular plane.
            ref = (1.0, 0.0, 0.0) if abs(cn_unit[0]) < 0.9 else (0.0, 1.0, 0.0)
            # e1 = ref - (ref . cn) * cn
            dot = ref[0] * cn_unit[0] + ref[1] * cn_unit[1] + ref[2] * cn_unit[2]
            e1 = (
                ref[0] - dot * cn_unit[0],
                ref[1] - dot * cn_unit[1],
                ref[2] - dot * cn_unit[2],
            )
            e1n = math.sqrt(sum(c * c for c in e1))
            if e1n < 1e-6:
                continue
            e1 = (e1[0] / e1n, e1[1] / e1n, e1[2] / e1n)
            # e2 = cn x e1
            e2 = (
                cn_unit[1] * e1[2] - cn_unit[2] * e1[1],
                cn_unit[2] * e1[0] - cn_unit[0] * e1[2],
                cn_unit[0] * e1[1] - cn_unit[1] * e1[0],
            )
            # Ring radius from centroid: r = d_centroid * tan(half-angle of
            # M-centroid-C cone) ~ d_centroid * 0.65 for η5 (matches Cp
            # geometry).  Empirical heuristic per eta.
            ring_r = {2: 0.70, 3: 0.85, 4: 1.00, 5: 1.20, 6: 1.40}.get(eta, 1.20)
            for k, atom_idx in enumerate(sorted(ring)):
                if atom_idx in hapto_atoms_placed:
                    continue
                theta = 2.0 * math.pi * k / eta
                cos_t = math.cos(theta)
                sin_t = math.sin(theta)
                pos = (
                    centroid[0] + ring_r * (cos_t * e1[0] + sin_t * e2[0]),
                    centroid[1] + ring_r * (cos_t * e1[1] + sin_t * e2[1]),
                    centroid[2] + ring_r * (cos_t * e1[2] + sin_t * e2[2]),
                )
                coords[atom_idx] = pos
                hapto_atoms_placed.add(atom_idx)

    # --- Level 2b: σ-donors (non-hapto, non-metal neighbours).
    # Place radially outward from metal through (metal - cluster_centroid)
    # at ``_get_ml_bond_length(M, D)``.  Multiple σ-donors per metal: fan
    # out at 90° increments in the plane perpendicular to radial.
    for mi in metal_idxs:
        m_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        m_pos = coords[mi]
        sigma_donors = []
        for nb in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nb.GetIdx()
            if nb.GetSymbol() in metal_set:
                continue
            if ni in hapto_atoms_placed:
                continue
            if nb.GetAtomicNum() == 1:
                continue  # H placed downstream
            sigma_donors.append(ni)
        if not sigma_donors:
            continue
        radial = (
            m_pos[0] - cluster_centroid[0],
            m_pos[1] - cluster_centroid[1],
            m_pos[2] - cluster_centroid[2],
        )
        r_len = math.sqrt(sum(c * c for c in radial))
        if r_len < 1e-6:
            radial_unit = (0.0, 0.0, 1.0)  # M at centroid → +z fallback
        else:
            radial_unit = (radial[0] / r_len, radial[1] / r_len, radial[2] / r_len)
        # Build perpendicular basis for the fan.
        ref = (0.0, 0.0, 1.0) if abs(radial_unit[2]) < 0.9 else (1.0, 0.0, 0.0)
        dot = ref[0] * radial_unit[0] + ref[1] * radial_unit[1] + ref[2] * radial_unit[2]
        e1 = (
            ref[0] - dot * radial_unit[0],
            ref[1] - dot * radial_unit[1],
            ref[2] - dot * radial_unit[2],
        )
        e1n = math.sqrt(sum(c * c for c in e1))
        if e1n < 1e-6:
            continue
        e1 = (e1[0] / e1n, e1[1] / e1n, e1[2] / e1n)
        e2 = (
            radial_unit[1] * e1[2] - radial_unit[2] * e1[1],
            radial_unit[2] * e1[0] - radial_unit[0] * e1[2],
            radial_unit[0] * e1[1] - radial_unit[1] * e1[0],
        )
        n_d = len(sigma_donors)
        # Distribute donors: single donor purely radial; multiple donors
        # at 70° cone half-angle from radial, evenly around.  Matches
        # piano-stool / half-sandwich M-X axial geometries.
        cone_half_angle = math.radians(70.0) if n_d > 1 else 0.0
        cos_a = math.cos(cone_half_angle)
        sin_a = math.sin(cone_half_angle)
        for k, d_idx in enumerate(sorted(sigma_donors)):
            d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
            try:
                bond_len = float(ml_bond_lookup(m_sym, d_sym))
            except Exception:
                bond_len = 2.0
            if n_d == 1:
                direction = radial_unit
            else:
                phi = 2.0 * math.pi * k / n_d
                cos_p = math.cos(phi)
                sin_p = math.sin(phi)
                direction = (
                    cos_a * radial_unit[0] + sin_a * (cos_p * e1[0] + sin_p * e2[0]),
                    cos_a * radial_unit[1] + sin_a * (cos_p * e1[1] + sin_p * e2[1]),
                    cos_a * radial_unit[2] + sin_a * (cos_p * e1[2] + sin_p * e2[2]),
                )
            coords[d_idx] = (
                m_pos[0] + direction[0] * bond_len,
                m_pos[1] + direction[1] * bond_len,
                m_pos[2] + direction[2] * bond_len,
            )

    return coords
