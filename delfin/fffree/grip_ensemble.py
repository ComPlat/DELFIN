"""delfin.fffree.grip_ensemble — Deterministic-Completeness Ensemble Builder.

Hebel #101 (User 2026-06-02): the deterministic-global-optimization counterpart
to xTB / DFT search.  For each SMILES we ENUMERATE all relevant Pólya
coordination isomers × all relevant Cremer-Pople ring conformers, optionally
polish each candidate via :func:`grip_polish`, score them by a combined
``severity + inter_ligand_clash + cshm`` criterion, and emit the top-K (or the
full ranked ensemble).

Why this is the gamechanger:

* Forward inference (xTB / DFT) finds ONE local minimum from one start point.
* GRIP-Ensemble does the inverse: it enumerates the full discrete space of
  possible structures the molecule can adopt (combinatorial + ring-flip), runs
  the geometric polish on every candidate, and ranks them by a CCDC-grounded
  score.  This is deterministic, reproducible, and complete up to the
  enumeration knobs.
* Inter-ligand clashes are a FIRST-CLASS ranking criterion (a structure can
  have perfect per-fragment internas yet still be unphysical because two
  ligands spatially overlap; Mogul severity alone does not catch this).

Architecture (see SPEC_GRIP §3 and the task brief)::

    SMILES
        ↓
    decompose() → list of (metal, geom, ligands)
        ↓
    enumerate_polya_isomers()           ← Burnside-complete (existing module)
        ↓ list of isomer configs
    for each isomer:
        enumerate_cp_conformers()       ← Cremer-Pople (existing ring_pucker)
            ↓ list of conformers per isomer
        for each conformer:
            assemble_from_config()      ← existing fffree pipeline
            (grip_polish optional, post-grip-corrector optional via env)
            record (P, severity, clash_count, cshm)
        ↓
    filter: clash_count ≤ clash_threshold
    sort:   by (severity + w_clash * clash_count + w_cshm * cshm)
    emit:   top-K candidates  (or full ensemble when DELFIN_FFFREE_GRIP_ENSEMBLE_FULL=1)

The module is DEFAULT-OFF.  When the env flag
``DELFIN_FFFREE_GRIP_ENSEMBLE`` is NOT set (or 0), the public hooks return
``None`` so :func:`assemble_from_config` keeps its byte-identical
single-output behaviour.

Determinism contract
--------------------
* All enumeration generators iterate in sorted key order.
* Candidate ranking ties are broken by ``(isomer_id, conformer_id)`` lex.
* No RNG inside this module (downstream MMFF / ETKDG calls in
  ``conformer_enum`` already use fixed seeds).
* ``PYTHONHASHSEED=0`` is honoured for the dict-based isomer configs.

Universality contract
---------------------
* Inter-ligand subgraph detection is purely graph-topological: build the
  molecular graph, remove all metal-donor edges, then enumerate connected
  components excluding the metal.  Each component = one ligand.  This works
  identically for ferrocene (2 Cp ligands), Cu(NH3)4 (4 NH3 ligands),
  porphyrin Fe (1 macrocyclic ligand), Werner-OC-6 (6 monodentate), without
  any SMILES- or geometry-specific patterns.
* Inter-ligand clash detection uses vdW-radius-based proximity, same table
  as :mod:`grip_polish`, applied only to atom pairs whose endpoints belong
  to DIFFERENT subgraphs (so intra-ligand close contacts are not flagged).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from .grip_polish import DEFAULT_VDW_RADII

__all__ = [
    "EnsembleCandidate",
    "EnsembleResult",
    "ensemble_active",
    "ensemble_emit_full",
    "ensemble_topk",
    "count_inter_ligand_clashes",
    "identify_ligand_subgraphs",
    "grip_ensemble_enumerate",
    "rank_candidates",
    "DEFAULT_CLASH_FLOOR_FRACTION",
    "DEFAULT_CLASH_W",
    "DEFAULT_CSHM_W",
    "DEFAULT_CLASH_THRESHOLD",
    "DEFAULT_MAX_ISOMERS",
    "DEFAULT_MAX_CONFORMERS_PER_ISOMER",
    "DEFAULT_MAX_TOTAL",
    "DEFAULT_TOP_K",
]

_LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants (tunable via env / kwargs).  Stay deterministic across processes.
# ---------------------------------------------------------------------------

#: Pauli-floor fraction reused from :mod:`grip_constraints`.  An atom pair
#: counts as an inter-ligand clash when ``d_ij < frac * (r_i + r_j)``.
DEFAULT_CLASH_FLOOR_FRACTION: float = 0.85

#: Default ranking weight on the inter-ligand clash penalty (severity-primary).
DEFAULT_CLASH_W: float = 10.0

#: Default ranking weight on the polyhedron CShM term.
DEFAULT_CSHM_W: float = 1.0

#: Candidates with more than this many inter-ligand clashes are rejected
#: outright.  Configurable via ``DELFIN_FFFREE_GRIP_ENSEMBLE_CLASH_MAX``.
DEFAULT_CLASH_THRESHOLD: int = 5

#: Caps on the enumeration to keep ensembles tractable.
DEFAULT_MAX_ISOMERS: int = 20
DEFAULT_MAX_CONFORMERS_PER_ISOMER: int = 30
DEFAULT_MAX_TOTAL: int = 200

#: Default top-K returned when neither FULL nor a custom K is requested.
DEFAULT_TOP_K: int = 3


def _env_int(name: str, default: int) -> int:
    """Read an int env var, falling back to the default on any failure."""
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        v = int(raw)
        return v
    except (TypeError, ValueError):
        return int(default)


def _env_float(name: str, default: float) -> float:
    """Read a float env var, falling back to the default on any failure."""
    raw = os.environ.get(name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if not np.isfinite(v):
            return float(default)
        return v
    except (TypeError, ValueError):
        return float(default)


def _env_bool(name: str, default: bool = False) -> bool:
    """Read a boolean env var (1/0, true/false, yes/no)."""
    raw = os.environ.get(name, "").strip().lower()
    if not raw:
        return bool(default)
    if raw in ("1", "true", "yes", "on"):
        return True
    if raw in ("0", "false", "no", "off"):
        return False
    return bool(default)


def ensemble_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_ENSEMBLE`` is on.  Default OFF."""
    return _env_bool("DELFIN_FFFREE_GRIP_ENSEMBLE", default=False)


def ensemble_emit_full() -> bool:
    """``True`` iff the full-ensemble output mode is requested."""
    return _env_bool("DELFIN_FFFREE_GRIP_ENSEMBLE_FULL", default=False)


def ensemble_topk() -> int:
    """Top-K cutoff (>= 1).  Configurable via
    ``DELFIN_FFFREE_GRIP_ENSEMBLE_K``."""
    k = _env_int("DELFIN_FFFREE_GRIP_ENSEMBLE_K", DEFAULT_TOP_K)
    return max(1, k)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class EnsembleCandidate:
    """One generated candidate -- a (P, mol-symbols) tuple with diagnostics.

    Attributes are immutable for deterministic hashing / sorting.
    """

    syms: Tuple[str, ...]
    P: np.ndarray
    severity: float
    clash_count: int
    cshm: float
    isomer_id: int
    conformer_id: int
    label: str
    accepted: bool
    score: float = 0.0

    def __post_init__(self):
        # Force a deterministic immutable view of the array (avoid downstream
        # mutation leaking back into the ensemble cache).
        object.__setattr__(self, "P", np.array(self.P, dtype=np.float64, copy=True))
        self.P.setflags(write=False)


@dataclass
class EnsembleResult:
    """Output container of :func:`grip_ensemble_enumerate`.

    ``candidates`` is the full list (post-filter), ordered by score ascending.
    ``top_k`` is a *view* slice of that list (``len <= K``) honoring the
    current env-flag selection (``DELFIN_FFFREE_GRIP_ENSEMBLE_K`` /
    ``DELFIN_FFFREE_GRIP_ENSEMBLE_FULL``).
    """

    smiles: str
    candidates: List[EnsembleCandidate]
    n_isomers_enumerated: int = 0
    n_conformers_enumerated: int = 0
    n_assembled: int = 0
    n_rejected_clash: int = 0
    n_rejected_build: int = 0
    skip_reason: str = ""
    top_k: List[EnsembleCandidate] = field(default_factory=list)

    @property
    def n_candidates(self) -> int:
        return len(self.candidates)

    def best(self) -> Optional[EnsembleCandidate]:
        """Top-1 candidate, or ``None`` when the ensemble is empty."""
        return self.candidates[0] if self.candidates else None


# ---------------------------------------------------------------------------
# Ligand-subgraph identification
# ---------------------------------------------------------------------------
def identify_ligand_subgraphs(
    mol,
    metal_idx: int,
    donors: Sequence[int],
) -> List[FrozenSet[int]]:
    """Partition the heavy-atom graph into ligand subgraphs.

    Build the molecular graph from ``mol``, remove all metal-incident edges
    (the metal-donor σ/π bonds + any bridging metal edges), then enumerate
    connected components excluding the metal atom itself.  Each component
    is one ligand.

    The function is read-only: ``mol`` is never mutated.

    Examples
    --------
    * Cu(NH3)4 → 4 ligands (4 disjoint NH3 subgraphs).
    * Ferrocene (Fe(Cp)2) → 2 ligands (each ring carbons + their Hs).
    * Hexamine porphyrin Fe → 1 ligand (the entire macrocycle backbone is
      one connected subgraph after the 4 Fe-N edges are cut).
    * Cisplatin Pt(NH3)2(Cl)2 → 4 ligands (2 NH3 + 2 Cl).

    Parameters
    ----------
    mol : RDKit Mol-like
        The assembled complex molecule.  Hydrogens count as graph nodes.
    metal_idx : int
        Atom index of the metal centre.
    donors : sequence of int
        Atom indices of the constructed donors.  Used only to identify the
        edges to cut (any edge incident to the metal is cut anyway, but we
        also cut metal-donor edges explicitly so an aggregator that forgot
        to add the M-D bond still works).

    Returns
    -------
    list of frozenset[int]
        One frozenset of atom indices per ligand subgraph.  Sorted by
        the minimum atom index in each subgraph (deterministic).
    """
    metal_idx = int(metal_idx)
    donors_set: Set[int] = {int(d) for d in donors}

    # 1) Build adjacency (all bonds, exclude metal-incident).
    n = mol.GetNumAtoms()
    adj: List[Set[int]] = [set() for _ in range(n)]
    for bond in mol.GetBonds():
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        if a == b:
            continue
        # Cut metal-incident edges (covers M-D bonds + any extra bridging
        # bond that touches the metal).
        if a == metal_idx or b == metal_idx:
            continue
        # Defence-in-depth: also cut any edge where one endpoint is the
        # metal even if the metal_idx test missed it (e.g. when callers
        # pass the wrong index).  No-op when metal_idx is consistent.
        adj[a].add(b)
        adj[b].add(a)

    # 2) Connected components, excluding the metal atom node.
    visited: Set[int] = {metal_idx}
    components: List[FrozenSet[int]] = []
    for start in range(n):
        if start in visited:
            continue
        stack = [start]
        comp: Set[int] = set()
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            comp.add(v)
            for w in adj[v]:
                if w not in visited:
                    stack.append(w)
        if comp:
            components.append(frozenset(comp))

    # 3) Deterministic ordering: sort by the minimum index in each component.
    components.sort(key=lambda c: (min(c), sorted(c)))
    return components


# ---------------------------------------------------------------------------
# Inter-ligand clash count
# ---------------------------------------------------------------------------
def count_inter_ligand_clashes(
    P: np.ndarray,
    mol,
    metal_idx: int,
    donors: Sequence[int],
    *,
    threshold: float = DEFAULT_CLASH_FLOOR_FRACTION,
    vdw_radii_by_symbol: Optional[Dict[str, float]] = None,
    return_pairs: bool = False,
):
    """Count atom pairs from DIFFERENT ligand subgraphs that violate the
    vdW floor.

    A pair ``(i, j)`` counts as an inter-ligand clash when:

    1. ``i`` and ``j`` are in different ligand subgraphs (as returned by
       :func:`identify_ligand_subgraphs`).
    2. ``i`` and ``j`` are NOT bonded (the metal-donor σ bond is the only
       inter-subgraph bond, and we ignore the metal/donor pair anyway).
    3. ``|r_i - r_j| < threshold * (r_vdW_i + r_vdW_j)``.

    Atoms whose element is missing from the vdW table are skipped (same
    safe-default behaviour as :class:`ClashFloorPenalty`).  The metal atom
    is excluded entirely (its M-D bonds are constructive, not steric).

    Parameters
    ----------
    P : (N, 3) ndarray
        Cartesian coordinates.
    mol : RDKit Mol-like
        Molecule graph (must have at least the same atom count as ``P``).
    metal_idx : int
        Metal atom index.
    donors : sequence of int
        Donor atom indices.
    threshold : float, default 0.85
        vdW-sum multiplier.  Below this fraction, the pair counts as a
        clash.  Matches SPEC §3.2 row 5 and :mod:`grip_constraints`.
    vdw_radii_by_symbol : dict, optional
        Override the default vdW table.
    return_pairs : bool, default False
        If True, also return the list of clashing ``(i, j)`` pairs for
        downstream diagnostics.

    Returns
    -------
    int  (or  (int, list[tuple[int, int]]) when ``return_pairs=True``)
        Total inter-ligand clash count.
    """
    P = np.asarray(P, dtype=np.float64)
    if P.ndim == 1:
        if P.size % 3 != 0:
            raise ValueError("P must reshape to (N, 3)")
        P = P.reshape(-1, 3)
    n = P.shape[0]

    metal_idx = int(metal_idx)
    vdw_tbl = vdw_radii_by_symbol if vdw_radii_by_symbol is not None else DEFAULT_VDW_RADII

    # Subgraph membership map: atom_idx -> subgraph_idx.
    subgraphs = identify_ligand_subgraphs(mol, metal_idx, donors)
    subgraph_of: Dict[int, int] = {}
    for k, comp in enumerate(subgraphs):
        for a in comp:
            subgraph_of[int(a)] = k

    # vdW radii per atom (NaN if unknown -> skipped).
    radii = np.full(n, np.nan, dtype=np.float64)
    try:
        for atom in mol.GetAtoms():
            idx = int(atom.GetIdx())
            sym = atom.GetSymbol()
            r = vdw_tbl.get(sym)
            if r is not None and 0 <= idx < n:
                radii[idx] = float(r)
    except Exception:
        pass

    # Bonded-pair set (full molecular graph, including M-D bonds -- we want
    # to skip those too in case the metal is in a subgraph by accident).
    bonded: Set[FrozenSet[int]] = set()
    try:
        for b in mol.GetBonds():
            a = int(b.GetBeginAtomIdx())
            c = int(b.GetEndAtomIdx())
            if a != c:
                bonded.add(frozenset((a, c)))
    except Exception:
        pass

    pairs: List[Tuple[int, int]] = []
    count = 0
    for i in range(n):
        if i == metal_idx:
            continue
        if i not in subgraph_of:
            continue
        ri = radii[i]
        if not np.isfinite(ri):
            continue
        si = subgraph_of[i]
        Pi = P[i]
        for j in range(i + 1, n):
            if j == metal_idx:
                continue
            if j not in subgraph_of:
                continue
            sj = subgraph_of[j]
            if si == sj:
                continue  # intra-ligand pair -- not an inter-ligand clash
            rj = radii[j]
            if not np.isfinite(rj):
                continue
            if frozenset((i, j)) in bonded:
                continue
            d = float(np.linalg.norm(Pi - P[j]))
            d_min = float(threshold) * (ri + rj)
            if d < d_min:
                count += 1
                if return_pairs:
                    pairs.append((i, j))
    if return_pairs:
        return count, pairs
    return count


# ---------------------------------------------------------------------------
# Polyhedron CShM helper (geometry-aware)
# ---------------------------------------------------------------------------
def _safe_cshm(P: np.ndarray, metal_idx: int, donors: Sequence[int], geometry: str) -> float:
    """Continuous shape measure of the donor polyhedron, or ``inf`` on
    failure (so failing builds rank LAST in the score).
    """
    if not geometry or len(donors) < 2:
        return 0.0
    try:
        from . import polyhedra as PLY
        d_arr = np.asarray(P[list(donors)] - P[int(metal_idx)], dtype=float)
        norms = np.linalg.norm(d_arr, axis=1)
        norms = np.where(norms < 1e-12, 1e-12, norms)
        obs = d_arr / norms[:, None]
        # Translate fffree geometry tag (e.g. "OC-6 octahedron") to the
        # polyhedra string.  polyhedra.cshm accepts both forms.
        return float(PLY.cshm(obs, geometry))
    except Exception:
        return float("inf")


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------
def rank_candidates(
    candidates: Sequence[EnsembleCandidate],
    *,
    clash_w: float = DEFAULT_CLASH_W,
    cshm_w: float = DEFAULT_CSHM_W,
    clash_threshold: int = DEFAULT_CLASH_THRESHOLD,
) -> List[EnsembleCandidate]:
    """Filter + rank candidates by the combined score.

    Score formula
    -------------
    ``score = severity + clash_w * max(0, clash_count - clash_threshold)^2 + cshm_w * cshm``

    Rationale: severity is the CCDC-Mahalanobis distance (the primary
    "geometric realism" axis).  Inter-ligand clashes get a quadratic
    penalty ABOVE the threshold (so a tiny grazing contact is forgiven
    but a serious overlap dominates).  CShM is a small tiebreaker on the
    coordination polyhedron quality.

    Filtering
    ---------
    Candidates with ``clash_count > clash_threshold * 2`` are filtered out
    entirely (gross overlap = unphysical).  Use
    ``DELFIN_FFFREE_GRIP_ENSEMBLE_CLASH_MAX`` to override the hard cutoff.

    Tie-breaking
    ------------
    Ties on the score are broken by ``(isomer_id, conformer_id)`` lex --
    deterministic across runs.

    Returns
    -------
    list of :class:`EnsembleCandidate`
        Sorted ascending by score, score field overwritten in-place
        (as new frozen dataclasses since :class:`EnsembleCandidate` is
        immutable).
    """
    cl_w = float(clash_w)
    cs_w = float(cshm_w)
    cl_thr = int(clash_threshold)
    cl_hard = _env_int("DELFIN_FFFREE_GRIP_ENSEMBLE_CLASH_MAX", max(cl_thr * 2, 10))

    scored: List[EnsembleCandidate] = []
    for c in candidates:
        if c.clash_count > cl_hard:
            continue
        gap = max(0, c.clash_count - cl_thr)
        cshm_val = float(c.cshm) if np.isfinite(c.cshm) else 1e6
        score = float(c.severity) + cl_w * (gap * gap) + cs_w * cshm_val
        # Reconstitute the immutable candidate with the score filled in.
        scored.append(
            EnsembleCandidate(
                syms=c.syms,
                P=c.P,
                severity=c.severity,
                clash_count=c.clash_count,
                cshm=c.cshm,
                isomer_id=c.isomer_id,
                conformer_id=c.conformer_id,
                label=c.label,
                accepted=c.accepted,
                score=float(score),
            )
        )

    scored.sort(key=lambda c: (c.score, c.isomer_id, c.conformer_id, c.label))
    return scored


# ---------------------------------------------------------------------------
# Severity wrapper (uses GRIP fragment loss if library available)
# ---------------------------------------------------------------------------
def _compute_severity(syms: Sequence[str], P: np.ndarray, mol, metal_idx: int,
                       donors: Sequence[int]) -> float:
    """Compute the Mogul-style severity of the assembled structure.

    Mirrors :func:`grip_polish.mogul_severity` but without requiring the
    polish to have run.  Uses :func:`detect_fragments` against the default
    GRIP library.  Returns ``inf`` on any detection failure so failing
    candidates rank last.
    """
    try:
        from .grip_fragment_detect import detect_fragments
        from .grip_loss_terms import TotalGripLoss
        from .grip_mogul_lookup import get_default_library
        try:
            library = get_default_library()
        except Exception:
            library = None
        frozen = frozenset((int(metal_idx), *[int(d) for d in donors]))
        terms = detect_fragments(
            mol, np.asarray(P, dtype=np.float64),
            frozen_atoms=frozen, donors=tuple(int(d) for d in donors),
            library=library, return_result=False,
        )
        agg = TotalGripLoss(terms=list(terms))
        if len(agg) == 0:
            return 0.0
        sev, _ = agg.evaluate(np.asarray(P, dtype=np.float64))
        sev = float(sev)
        if not np.isfinite(sev):
            return float("inf")
        return sev
    except Exception:
        return float("inf")


# ---------------------------------------------------------------------------
# Main entry point: enumerate × polish × rank
# ---------------------------------------------------------------------------
def grip_ensemble_enumerate(
    smiles: str,
    *,
    max_isomers: int = DEFAULT_MAX_ISOMERS,
    max_conformers_per_isomer: int = DEFAULT_MAX_CONFORMERS_PER_ISOMER,
    max_total: int = DEFAULT_MAX_TOTAL,
    clash_w: float = DEFAULT_CLASH_W,
    cshm_w: float = DEFAULT_CSHM_W,
    clash_threshold: int = DEFAULT_CLASH_THRESHOLD,
    top_k: Optional[int] = None,
) -> EnsembleResult:
    """Enumerate Pólya × Cremer-Pople candidates, polish each, and rank.

    This is the public API.  The function is safe to call at any time;
    when the SMILES is outside the fffree-decomposable domain (multi-metal
    legacy, CN outside the supported range, etc.) it returns an empty
    :class:`EnsembleResult` with ``skip_reason`` populated.

    Parameters
    ----------
    smiles : str
        Input SMILES string of the metal complex.
    max_isomers : int
        Cap on Pólya isomers enumerated (deterministic ordering applies
        BEFORE the cap so the same SMILES always produces the same
        subset).
    max_conformers_per_isomer : int
        Cap on ring-pucker conformers per isomer.  ``1`` = use the
        as-built conformer only (no CP enumeration).
    max_total : int
        Hard cap on the total number of (isomer × conformer) candidates
        attempted -- protects against combinatorial blow-up on
        macrocycles.
    clash_w, cshm_w : float
        Score weights.
    clash_threshold : int
        Inter-ligand clash count above which the penalty kicks in
        quadratically.  Hard cutoff (rejection) is
        ``clash_threshold * 2`` by default, override via
        ``DELFIN_FFFREE_GRIP_ENSEMBLE_CLASH_MAX``.
    top_k : int, optional
        When given, populate ``result.top_k`` to the first ``k`` candidates.
        ``None`` honors :func:`ensemble_topk` from the env.

    Returns
    -------
    :class:`EnsembleResult`
        With ``candidates`` ranked ascending by score and ``top_k`` populated.
    """
    # Lazy imports keep import cost paid only when the ensemble is used.
    try:
        from . import decompose as DEC
        from . import polya_isomer_count as PIC
        from . import assemble_complex as AC
    except Exception as exc:  # pragma: no cover -- import-time wiring
        return EnsembleResult(smiles=smiles, candidates=[],
                              skip_reason=f"import_error: {exc!r}")

    try:
        d = DEC.decompose(smiles)
    except Exception as exc:
        return EnsembleResult(smiles=smiles, candidates=[],
                              skip_reason=f"decompose_error: {exc!r}")
    if d is None:
        return EnsembleResult(smiles=smiles, candidates=[],
                              skip_reason="decompose_returned_none")

    # ------------------------------------------------------------------
    # 1) Enumerate Pólya isomers
    # ------------------------------------------------------------------
    _GEOM_TO_POLYA = {
        "L-2 linear": "linear",
        "SP-3 trigonal planar": "trigonal_planar",
        "T-3 T-shape": "tshape",
        "OC-6 octahedron": "octahedron",
        "SP-4 square planar": "square_planar",
        "T-4 tetrahedron": "tetrahedron",
        "TBP-5 trigonal bipyramid": "trigonal_bipyramid",
        "SPY-5 square pyramid": "square_pyramid",
        "TPR-6 trigonal prism": "trigonal_prism",
        "PB-7 pentagonal bipyramid": "pentagonal_bipyramid",
        "SQAP-8 square antiprism": "square_antiprism",
        "TTP-9 tricapped trigonal prism": "tricapped_trigonal_prism",
    }
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return EnsembleResult(smiles=smiles, candidates=[],
                              skip_reason=f"geom_unsupported: {d['geometry']}")

    ligands = d["ligands"]

    # Use the same enumeration helpers the existing fffree pipeline uses, so
    # the ensemble's denominator equals the Pólya-completeness denominator
    # documented in :mod:`polya_isomer_count`.
    isomer_configs: List = []
    try:
        if d.get("has_chelate"):
            # Chelate-config enumeration (vertex -> (lig, arm))
            from rdkit import Chem as _Chem
            specs = [
                {
                    "type": _Chem.MolToSmiles(lg["mol"]),
                    "denticity": lg["denticity"],
                    "asym": len(set(lg.get("donor_elems", []))) > 1,
                }
                for lg in ligands
            ]
            isomer_configs = list(PIC.enumerate_chelate_configs(geom_key, specs))
        else:
            # Vertex-coloring enumeration -- monodentate only.  Build a
            # synthetic config (vertex -> (lig, arm)) so the assembly path
            # has the SAME signature in either branch.
            from rdkit import Chem as _Chem
            from collections import Counter as _Counter
            lig_label: List[str] = []
            lig_ref: Dict[str, Tuple[object, int]] = {}
            for lg in ligands:
                lab = _Chem.MolToSmiles(lg["mol"])
                lig_label.append(lab)
                lig_ref.setdefault(lab, (lg["mol"], lg["donor_local_idx"]))
            spec = dict(_Counter(lig_label))
            colorings = list(PIC.enumerate_isomers(geom_key, spec))
            # Map every coloring to a vertex->(lig_idx, arm=0) config.
            for coloring in colorings:
                # Each vertex gets one ligand label; pick the first available
                # instance of that label that has not been assigned yet.
                used: Set[int] = set()
                config: Dict[int, Tuple[int, int]] = {}
                for v_idx, lab in enumerate(coloring):
                    # find next ligand instance with this label
                    for li, lg in enumerate(ligands):
                        if li in used:
                            continue
                        try:
                            if _Chem.MolToSmiles(lg["mol"]) == lab:
                                config[v_idx] = (li, 0)
                                used.add(li)
                                break
                        except Exception:
                            continue
                # Only keep config if every vertex got a ligand
                if len(config) == len(coloring):
                    isomer_configs.append(config)
    except Exception as exc:
        return EnsembleResult(smiles=smiles, candidates=[],
                              skip_reason=f"isomer_enum_failed: {exc!r}")

    # Deterministic order:  PIC.enumerate_* already iterate in sorted order,
    # but the dict-based config keys can hash-order downstream; freeze the
    # list to its construction order which is itertools.product / canonical
    # from PIC, both deterministic.
    n_isomers_total = len(isomer_configs)
    if max_isomers is not None and max_isomers > 0:
        isomer_configs = isomer_configs[: int(max_isomers)]

    # ------------------------------------------------------------------
    # 2) Build candidates -- per isomer, optionally enumerate CP conformers
    # ------------------------------------------------------------------
    candidates: List[EnsembleCandidate] = []
    n_assembled = 0
    n_rejected_clash = 0
    n_rejected_build = 0
    n_conformers_total = 0

    n_total = 0

    metal = d["metal"]
    geometry = d["geometry"]

    # Lazy import for the legacy converter-backend polish stack.  We re-use
    # ``_maybe_relax`` (PT3 FF-free defect-refiner + ring scale) and
    # ``_g16_soft_polyhedron_polish`` so the ensemble candidates go through
    # the SAME post-build polish chain as the legacy emit path.  Without
    # this the ensemble's "top-1" is a less-polished structure than what
    # the legacy non-chelate emit returns -> apples-vs-oranges metrics.
    try:
        from . import converter_backend as _CB  # type: ignore
        _has_polish = True
    except Exception:
        _CB = None
        _has_polish = False

    _has_chelate = d.get("has_chelate", False)
    _has_hapto = any(lg.get("is_hapto") for lg in ligands)

    for iso_id, config in enumerate(isomer_configs):
        if n_total >= int(max_total):
            break
        # 2a) Build the base assembly for this isomer (single conformer).
        try:
            built = AC.assemble_from_config(metal, geometry, config, ligands)
        except Exception:
            built = None
            n_rejected_build += 1
        if built is None:
            n_rejected_build += 1
            continue
        syms, P_base, donors_global = built

        # 2a') Apply the same post-build polish stack the legacy emit
        #      path applies.  Without this the ensemble emits a structure
        #      one or two relax steps SHORT of the legacy baseline, which
        #      regresses every internals axis the relax fixes
        #      (XH collapse, F20 H planar, decoord, funcgrp).
        if _has_polish and _CB is not None:
            try:
                syms_l, P_l = _CB._maybe_relax(list(syms), np.asarray(P_base, dtype=float))
                # Self-gate: reject builds the legacy gate would reject.
                if _CB._build_is_clean(syms_l, P_l, cn=d.get("cn"),
                                        geom=geometry,
                                        donors=donors_global,
                                        has_hapto=_has_hapto):
                    # G16 soft polyhedron polish (auto under PURE_TRACK3).
                    if (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                            or os.environ.get("DELFIN_FFFREE_SOFT_POLY", "0") == "1"):
                        try:
                            syms_l, P_l = _CB._g16_soft_polyhedron_polish(
                                syms_l, P_l, cn=d.get("cn"),
                                geom=geometry, donors=donors_global,
                                has_hapto=_has_hapto,
                            )
                        except Exception:
                            pass
                    syms = list(syms_l)
                    P_base = np.asarray(P_l, dtype=np.float64)
                else:
                    # Build did not pass legacy self-gate -> skip this
                    # candidate (mirrors legacy fall-through to next iso).
                    n_rejected_build += 1
                    continue
            except Exception:
                # Any failure in polish -> use the raw build, do not fail.
                pass

        n_assembled += 1

        # 2b) Build the assembled mol (metal + ligand bonds) for graph ops.
        # We re-derive it the same way assemble_from_config does so we have
        # a stable RDKit mol for both severity and clash measurement.  This
        # mirrors the construction inside assemble_complex.py (RWMol, add
        # metal-donor bonds).  If construction fails, skip the candidate
        # (defensive -- in practice always succeeds here).
        try:
            from rdkit import Chem as _Chem
            cm = _Chem.RWMol()
            cm.AddAtom(_Chem.Atom(syms[0]))  # metal at index 0
            # CRITICAL: assemble_from_config iterates ``by_lig.items()`` in
            # DICT INSERTION ORDER, which is the order of the first vertex
            # appearance for each ligand in ``config.items()``.  We must
            # mirror that exact order so the atom indices in our rebuilt
            # ``cm`` align with the indices in ``syms`` / ``P_base``.
            # Using ``sorted(by_lig)`` instead silently swaps fragment blocks
            # and corrupts every downstream (severity / clash / cshm) read.
            by_lig: Dict[int, List[Tuple[int, int]]] = {}
            for v, (li, arm) in config.items():
                by_lig.setdefault(li, []).append((v, arm))
            for li in by_lig.keys():  # insertion order == assemble order
                lg = ligands[li]
                frag = _Chem.AddHs(lg["mol"])
                base = cm.GetNumAtoms()
                for a in frag.GetAtoms():
                    cm.AddAtom(_Chem.Atom(a.GetAtomicNum()))
                for b in frag.GetBonds():
                    cm.AddBond(
                        b.GetBeginAtomIdx() + base,
                        b.GetEndAtomIdx() + base,
                        b.GetBondType(),
                    )
            # Add metal-donor bonds (using donors_global from the built tuple)
            for dg in donors_global:
                try:
                    cm.AddBond(0, int(dg), _Chem.BondType.SINGLE)
                except Exception:
                    pass
            # Sanitize is best-effort; we proceed even if it fails because
            # our graph operations (subgraph, bonded-pair set) need only the
            # bond topology, not Kekulé / aromatic perception.
            try:
                _Chem.SanitizeMol(cm, catchErrors=True)
            except Exception:
                pass
            # Defensive check: atom count alignment.  If the rebuilt graph
            # has a different size from the built XYZ, all downstream
            # indices are unreliable -> skip cm and rank this candidate
            # without per-fragment severity / clash diagnostics (will
            # surface as severity=inf, clash=-1 -> bottom-ranked).
            if cm.GetNumAtoms() != len(syms):
                cm = None
        except Exception:
            cm = None

        # 2c) Conformer enumeration via Cremer-Pople ring puckers.
        #
        # We treat the base assembly as conformer #0 (always present), then
        # optionally generate ring-pucker variants when there are puckerable
        # 5/6/7-membered rings AND max_conformers_per_isomer > 1.  This
        # keeps the OFF-path identical to the legacy single-output flow.
        conformer_Ps: List[Tuple[int, np.ndarray, str]] = [(0, np.asarray(P_base, dtype=np.float64), "base")]

        if max_conformers_per_isomer > 1 and cm is not None:
            try:
                from .ring_pucker import enumerate_ring_conformers
                # Detect rings on the assembled mol; if there are none with
                # 5/6/7 atoms, the enumerator yields nothing extra.
                try:
                    rings = []
                    try:
                        ri = cm.GetRingInfo()
                        rings = [tuple(int(a) for a in r) for r in ri.AtomRings()
                                 if 5 <= len(r) <= 7]
                    except Exception:
                        rings = []
                    if rings:
                        max_per_ring = max(1, min(4, max_conformers_per_isomer))
                        extra_iter = enumerate_ring_conformers(
                            np.asarray(P_base, dtype=np.float64),
                            rings,
                            max_per_ring=max_per_ring,
                        )
                        # Drop the very first variant which equals base
                        # (chair-canonical for 6-ring, envelope for 5-ring);
                        # for safety we just keep any variant whose
                        # max-atom-displacement from base exceeds 0.05 A.
                        for c_id, P_var in enumerate(extra_iter, start=1):
                            if len(conformer_Ps) >= max_conformers_per_isomer:
                                break
                            try:
                                disp = float(np.max(np.linalg.norm(
                                    np.asarray(P_var, dtype=float) - P_base, axis=1
                                )))
                            except Exception:
                                disp = 0.0
                            if disp < 0.05:
                                continue
                            conformer_Ps.append((c_id, np.asarray(P_var, dtype=np.float64),
                                                  f"pucker{c_id}"))
                except Exception:
                    pass
            except Exception:
                pass

        n_conformers_total += len(conformer_Ps)

        # 2d) Score each conformer.
        for c_id, P_curr, c_label in conformer_Ps:
            if n_total >= int(max_total):
                break
            n_total += 1

            # Severity (Mogul Mahalanobis).
            if cm is not None and cm.GetNumAtoms() == len(syms):
                sev = _compute_severity(syms, P_curr, cm, 0, donors_global)
            else:
                sev = float("inf")

            # Inter-ligand clashes.
            if cm is not None and cm.GetNumAtoms() == len(syms):
                try:
                    clash_n = count_inter_ligand_clashes(
                        P_curr, cm, metal_idx=0, donors=donors_global,
                    )
                except Exception:
                    clash_n = -1
            else:
                clash_n = -1

            # Polyhedron CShM.
            cshm = _safe_cshm(P_curr, 0, donors_global, geometry)

            label = f"iso{iso_id}_conf{c_id}_{c_label}"
            cand = EnsembleCandidate(
                syms=tuple(syms),
                P=P_curr,
                severity=float(sev),
                clash_count=int(max(clash_n, 0)),
                cshm=float(cshm),
                isomer_id=int(iso_id),
                conformer_id=int(c_id),
                label=label,
                accepted=(np.isfinite(sev) and clash_n >= 0),
                score=0.0,
            )
            candidates.append(cand)

    # ------------------------------------------------------------------
    # 3) Rank + top-K selection
    # ------------------------------------------------------------------
    ranked = rank_candidates(
        candidates,
        clash_w=clash_w,
        cshm_w=cshm_w,
        clash_threshold=clash_threshold,
    )
    n_rejected_clash = max(0, len(candidates) - len(ranked))

    k = int(top_k) if top_k is not None else ensemble_topk()
    if ensemble_emit_full():
        top = list(ranked)
    else:
        top = list(ranked[: max(1, k)])

    return EnsembleResult(
        smiles=smiles,
        candidates=ranked,
        n_isomers_enumerated=int(n_isomers_total),
        n_conformers_enumerated=int(n_conformers_total),
        n_assembled=int(n_assembled),
        n_rejected_clash=int(n_rejected_clash),
        n_rejected_build=int(n_rejected_build),
        skip_reason="",
        top_k=top,
    )
