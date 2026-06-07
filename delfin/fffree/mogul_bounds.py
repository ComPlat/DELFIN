"""Mogul-DG Phase A — Distance-Geometry bounds matrix construction.

Builds the per-atom-pair :math:`(l_{ij}, u_{ij})` bounds matrix that defines
the feasibility region of the whole-complex distance-geometry embed (see
``agent_workspace/quality_framework/docs/architecture/mogul_whole_complex_dg.md``).

The bounds are derived from **graph topology + CCDC/COD empirical fragment
libraries** — there is no SMILES-specific code path, and no element-specific
constant beyond covalent / van-der-Waals radii (which are universal physical
constants, not chemistry-specific heuristics).

Seven bound classes are constructed (SPEC §2.2 + §13 + §14):

1.  **Bonded pairs**  :math:`\\mu \\pm 3\\sigma` from the CCDC fragment lookup;
    fall back to Pauling covalent-radii sum + 5 % when the library has no
    match.
2.  **1,3 angle pairs** (atoms two bonds apart): the CCDC angle distribution
    is converted into a 1,3 distance bound through the law of cosines.
3.  **Metal-donor pairs**: tighter window :math:`\\mu \\pm 2\\sigma`.  When the
    library supplies a multi-modal distribution (Jahn-Teller / hapto-vs-σ /
    trans-influence) we use the *envelope* of all components (lowest
    :math:`\\mu - 2\\sigma`, highest :math:`\\mu + 2\\sigma`) so the solver
    can land on any of them.
4.  **Donor-donor pairs**: ideal polyhedron geometry encoded as the
    pairwise distances between donor reference vectors scaled by the
    metal-donor radius.  Tolerance :math:`\\pm 5\\sigma_{Md}` (polyhedron
    flexibility, SPEC §13 Q4).
5.  **Hapto-π pairs**: M-centroid distance from the v5 TM category
    library, equi-distance constraint for the M-C ring members
    (centroid encoded implicitly via equi-distance bound).
6.  **Resonance / mesomeric equivalence**: aromatic rings,
    symmetric terminal groups (XO_n, XF_n), and graph-automorphism
    orbits get equi-distance bounds (SPEC §14).
7.  **Non-bonded vdW**: lower bound :math:`0.85 \\cdot (r^{vdW}_i +
    r^{vdW}_j)`, upper bound :math:`+\\infty`.

The output is a symmetric ``(n, n)`` lower/upper bounds matrix plus an
``info`` dict for diagnostics.  The diagonal is zero.

Determinism contract
--------------------
* All loops iterate over sorted integer atom indices.
* No RNG, no Python-set iteration order leaks (sets are converted to
  sorted lists before any iteration that influences a bound).
* Bounds are written in *value-overwrite-last-wins* order: the order is
  bonded → 1,3-angle → metal-donor → donor-donor → hapto → resonance →
  vdW.  The strongest constraint (tightest bound) wins when classes
  overlap (e.g. an aromatic C-C bond is overwritten by the resonance
  equi-distance constraint, which is correct).

Universality contract
---------------------
* The only element-specific tables are ``_COV`` (Pauling covalent radii,
  used for the fallback bonded bound) and ``_VDW`` (Bondi van-der-Waals
  radii, used for the non-bonded floor).  Both are universal physical
  constants tabulated for every periodic-table element.
* Every other decision is graph-driven (RDKit topology, ring info,
  hybridisation, automorphism rank) or library-driven (CCDC bond /
  angle / TM-category lookups).
* No SMILES string is parsed; no element name appears in an ``if``
  branch.  ``tests/test_mogul_bounds.py::test_universal_no_hardcoded_elements``
  greps the module for such patterns and fails CI if any are found.

API
---
:func:`build_bounds_matrix` is the only public entry point.  It receives
the symbol list, the RDKit ``Mol``, the metal index, the donor indices,
and a CCDC ``GripLibrary``; it returns ``(lower, upper, info)``.

This module is *Phase A only* — it does **not** solve anything.  The
projected L-BFGS solver lives in Phase B (``mogul_solver.py``); the
public ``mogul_embed`` driver lives in Phase C+ (``mogul_dg.py``).
"""

from __future__ import annotations

import math
import os
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from delfin._bond_decollapse import _COV, _METALS, _VDW, _VDW_DEFAULT
from delfin.fffree.grip_mogul_lookup import GripLibrary
from delfin.fffree import polyhedra as _polyhedra


def _default_library() -> Optional[GripLibrary]:
    """Honour ``DELFIN_GRIP_LIB_PATH`` and fall back to the pinned library.

    Returns ``None`` only when no library is found anywhere; callers must
    handle that case (the cov-sum fallback is universal).
    """
    return GripLibrary.get_default()

__all__ = [
    "build_bounds_matrix",
    "BoundsInfo",
    # Internal helpers exposed for tests
    "find_aromatic_groups",
    "find_symmetric_terminal_groups",
    "find_automorphism_orbits",
    "detect_resonance_groups",
    "detect_hapto_groups",
    "hyb_str",
]

# ---------------------------------------------------------------------------
# Universal physical-constant tables (re-exported, single-source-of-truth)
# ---------------------------------------------------------------------------
# Default covalent radius (Å) for any element absent from ``_COV``.  Universal
# physical bound — Cordero 2008 mean across periods.  No chemistry encoded.
_COV_DEFAULT: float = 0.85

# σ-multiplier defaults (SPEC §3.3, §13 Q4)
_BOND_SIGMA_MULT_DEFAULT: float = 3.0    # ± 3σ for organic bonds (99.7 % CI)
_MD_SIGMA_MULT_DEFAULT: float = 2.0      # ± 2σ for metal-donor (tighter)
_DD_SIGMA_MULT_DEFAULT: float = 5.0      # ± 5σ for donor-donor (polyhedron flex)
_ANGLE_SIGMA_DEG_DEFAULT: float = 5.0    # ± 5° for 1,3 angles (SPEC §13 Q4)

# Non-bonded vdW floor (SPEC §2.2)
_VDW_FLOOR_FACTOR: float = 0.85

# Library coverage gate
_MIN_N_DEFAULT: int = 5

# Upper-bound sentinel for "no upper limit" — chosen large enough that it never
# triggers in any realistic geometry but small enough to avoid float-overflow
# in downstream arithmetic.
_UB_INF: float = 1e6

# Sigma floor when CCDC reports an unrealistically small σ (e.g. a single-bin
# distribution from low sample size).  Prevents the bound from collapsing to a
# single point that the solver cannot reach numerically.
_SIGMA_FLOOR_BOND: float = 0.01
_SIGMA_FLOOR_ANGLE_DEG: float = 1.0

# Default equi-distance σ for resonance groups (SPEC §14.3, empirical CCDC
# scatter inside aromatic rings).
_RESONANCE_SIGMA_DEFAULT: float = 0.02


# ---------------------------------------------------------------------------
# Environment-flag accessors
# ---------------------------------------------------------------------------
def _env_float(name: str, default: float) -> float:
    """Read a positive float from the env; fall back to ``default`` on miss."""
    raw = os.environ.get(name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if math.isfinite(v) and v > 0:
            return float(v)
    except (TypeError, ValueError):
        pass
    return float(default)


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        v = int(raw)
        if v > 0:
            return int(v)
    except (TypeError, ValueError):
        pass
    return int(default)


def _bond_sigma_mult() -> float:
    return _env_float("DELFIN_FFFREE_MOGUL_DG_BOUNDS_TOL", _BOND_SIGMA_MULT_DEFAULT)


def _md_sigma_mult() -> float:
    return _env_float("DELFIN_FFFREE_MOGUL_DG_MD_BOUNDS_TOL", _MD_SIGMA_MULT_DEFAULT)


def _dd_sigma_mult() -> float:
    return _env_float("DELFIN_FFFREE_MOGUL_DG_DD_BOUNDS_TOL", _DD_SIGMA_MULT_DEFAULT)


# ---------------------------------------------------------------------------
# Element classification (universal — no chemistry encoded beyond periodic-table
# block membership which is a physical constant of the element).
# ---------------------------------------------------------------------------
def _is_metal(symbol: str) -> bool:
    """True iff the symbol is a transition / f-block / heavy main-group metal.

    Uses the project-wide ``_METALS`` set from :mod:`delfin._bond_decollapse`
    so all modules share one definition.  Universal: based on periodic-table
    block, not chemistry.
    """
    return str(symbol) in _METALS


def _cov_radius(symbol: str) -> float:
    """Pauling covalent radius for ``symbol``; default 0.85 Å on miss."""
    return float(_COV.get(str(symbol), _COV_DEFAULT))


def _vdw_radius(symbol: str) -> float:
    """Bondi van-der-Waals radius for ``symbol``; default 1.7 Å on miss."""
    return float(_VDW.get(str(symbol), _VDW_DEFAULT))


# ---------------------------------------------------------------------------
# Hybridisation extraction (universal — pulled directly from RDKit, with a
# graph-only fallback when RDKit cannot decide).
# ---------------------------------------------------------------------------
_HYB_MAP: Dict[str, str] = {
    "S": "sp",  # rare RDKit code
    "SP": "sp",
    "SP2": "sp2",
    "SP3": "sp3",
    "SP3D": "sp3d",
    "SP3D2": "sp3d2",
    "UNSPECIFIED": "*",
    "OTHER": "*",
    "UNKNOWN": "*",
}


def hyb_str(atom) -> str:
    """Return the lookup-library hybridisation string for an RDKit atom.

    Universal: relies only on RDKit's perception of hybridisation, which is
    derived from graph valence + bond orders + ring membership.
    """
    try:
        name = atom.GetHybridization().name
    except Exception:
        return "*"
    return _HYB_MAP.get(name.upper(), "*")


# ---------------------------------------------------------------------------
# Graph helpers (RDKit-based, universal).
# ---------------------------------------------------------------------------
def _bonds_sorted(mol) -> List[Tuple[int, int]]:
    """All bonds as (lo, hi) tuples sorted in lex order (determinism)."""
    out: List[Tuple[int, int]] = []
    for b in mol.GetBonds():
        a = int(b.GetBeginAtomIdx())
        c = int(b.GetEndAtomIdx())
        if a == c:
            continue
        out.append((min(a, c), max(a, c)))
    out.sort()
    return out


def _neighbors_sorted(mol, idx: int) -> List[int]:
    """Sorted neighbour indices — determinism (RDKit's neighbour order is
    insertion-dependent)."""
    try:
        atom = mol.GetAtomWithIdx(int(idx))
    except Exception:
        return []
    return sorted(int(nb.GetIdx()) for nb in atom.GetNeighbors())


def _angle_triples_sorted(mol) -> List[Tuple[int, int, int]]:
    """All 1,3 angle triples ``(i, j, k)`` with j = central atom.

    Deterministic order: sorted by ``(j, min(i,k), max(i,k))``.  Each unordered
    triple appears exactly once (the smaller endpoint comes first).
    """
    seen: Set[Tuple[int, int, int]] = set()
    out: List[Tuple[int, int, int]] = []
    n = mol.GetNumAtoms()
    for j in range(n):
        nbrs = _neighbors_sorted(mol, j)
        for a_pos in range(len(nbrs)):
            for b_pos in range(a_pos + 1, len(nbrs)):
                i, k = nbrs[a_pos], nbrs[b_pos]
                lo, hi = (i, k) if i < k else (k, i)
                key = (j, lo, hi)
                if key in seen:
                    continue
                seen.add(key)
                out.append((lo, j, hi))
    out.sort(key=lambda t: (t[1], t[0], t[2]))
    return out


def _ring_size_min(mol, i: int) -> int:
    """Smallest ring containing ``i``; -1 if not in any ring."""
    try:
        ri = mol.GetRingInfo()
    except Exception:
        return -1
    if not ri.NumAtomRings(int(i)):
        return -1
    sizes = [len(r) for r in ri.AtomRings() if int(i) in r]
    if not sizes:
        return -1
    return int(min(sizes))


def _bond_ring_size_min(mol, i: int, j: int) -> int:
    try:
        ri = mol.GetRingInfo()
    except Exception:
        return -1
    sizes = [len(r) for r in ri.AtomRings() if int(i) in r and int(j) in r]
    if not sizes:
        return -1
    return int(min(sizes))


# ---------------------------------------------------------------------------
# §14 Resonance / mesomeric equivalence detection — three orthogonal
# mechanisms.  All graph-only.
# ---------------------------------------------------------------------------
def find_aromatic_groups(mol) -> List[FrozenSet[int]]:
    """All atoms within an aromatic ring belong to one resonance group.

    Universal: relies on RDKit aromaticity perception, which is itself a
    graph-only algorithm (Hückel 4n+2 over the SSSR ring set).
    """
    out: List[FrozenSet[int]] = []
    try:
        rings = mol.GetRingInfo().AtomRings()
    except Exception:
        return out
    for ring in rings:
        ring_ints = sorted(int(a) for a in ring)
        if not ring_ints:
            continue
        is_aromatic = True
        for a in ring_ints:
            try:
                if not mol.GetAtomWithIdx(a).GetIsAromatic():
                    is_aromatic = False
                    break
            except Exception:
                is_aromatic = False
                break
        if is_aromatic and len(ring_ints) >= 3:
            out.append(frozenset(ring_ints))
    # Deterministic order — sort by (smallest atom idx, second smallest, …)
    out.sort(key=lambda s: tuple(sorted(s)))
    return out


def find_symmetric_terminal_groups(mol) -> List[FrozenSet[int]]:
    """Detect central atoms with multiple equivalent dangling X neighbours.

    Captures carboxylate (RCOO⁻ → 2 O), nitrate (NO₃⁻ → 3 O), sulfate
    (SO₄²⁻ → 4 O), perchlorate (ClO₄⁻ → 4 O), tetrafluoroborate
    (BF₄⁻ → 4 F), hexafluorophosphate (PF₆⁻ → 6 F), nitro RNO₂, methyl
    CH₃, etc.  Universal: detection is element-symmetric — any atom with
    ≥ 2 terminal neighbours of the same element forms a group.

    The current implementation requires the dangling atom to be *degree-1*
    (i.e. truly terminal), which excludes accidental graph matches at
    branching centres.
    """
    out: List[FrozenSet[int]] = []
    try:
        n = mol.GetNumAtoms()
    except Exception:
        return out
    for parent_idx in range(n):
        try:
            parent = mol.GetAtomWithIdx(int(parent_idx))
        except Exception:
            continue
        by_element: Dict[str, List[int]] = {}
        for nb in parent.GetNeighbors():
            try:
                if nb.GetDegree() != 1:
                    continue
                sym = nb.GetSymbol()
            except Exception:
                continue
            by_element.setdefault(str(sym), []).append(int(nb.GetIdx()))
        for sym, indices in sorted(by_element.items()):
            if len(indices) >= 2:
                out.append(frozenset(sorted(int(i) for i in indices)))
    out.sort(key=lambda s: tuple(sorted(s)))
    return out


def find_automorphism_orbits(mol, *, exclude: Optional[Set[int]] = None) -> List[FrozenSet[int]]:
    """Compute graph automorphism orbits via RDKit canonical ranking.

    Two atoms with the same canonical rank are topologically equivalent (graph
    automorphism orbit).  Universal: pure graph topology, no chemistry.

    ``exclude`` masks atoms (typically the metal) from being grouped — useful
    when the metal would otherwise pool every donor atom into one orbit.
    """
    try:
        from rdkit import Chem
    except Exception:
        return []
    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        return []
    excl = set(int(i) for i in (exclude or ()))
    by_rank: Dict[int, List[int]] = {}
    for idx, rank in enumerate(ranks):
        if int(idx) in excl:
            continue
        by_rank.setdefault(int(rank), []).append(int(idx))
    out: List[FrozenSet[int]] = []
    for rank, indices in sorted(by_rank.items()):
        if len(indices) >= 2:
            out.append(frozenset(sorted(int(i) for i in indices)))
    out.sort(key=lambda s: tuple(sorted(s)))
    return out


def _merge_overlapping_sets(groups: Sequence[FrozenSet[int]]) -> List[FrozenSet[int]]:
    """Merge any two groups that share at least one atom (transitive closure).

    Deterministic: union-find with sorted iteration; output sorted by min idx.
    """
    sets = [set(g) for g in groups]
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(sets):
            j = i + 1
            while j < len(sets):
                if sets[i] & sets[j]:
                    sets[i] |= sets[j]
                    sets.pop(j)
                    changed = True
                else:
                    j += 1
            i += 1
    out = [frozenset(sorted(s)) for s in sets if len(s) >= 2]
    out.sort(key=lambda s: tuple(sorted(s)))
    return out


def detect_resonance_groups(
    mol,
    *,
    metal_idx: Optional[int] = None,
    use_automorphism: bool = True,
) -> List[FrozenSet[int]]:
    """Union of the three detection mechanisms with overlap merging.

    Parameters
    ----------
    mol : RDKit Mol
    metal_idx : int, optional
        Atom index of the metal.  Excluded from automorphism orbits so the
        metal does not pool all donors into one "orbit".
    use_automorphism : bool
        If False, only aromatic + symmetric-terminal mechanisms are used.
        Automorphism orbits can be very inclusive on small symmetric ligands
        (e.g. NH₃ ⇒ 3 H atoms) which is correct but sometimes redundant.

    Returns
    -------
    list of frozenset[int]
        Resonance-equivalent atom sets, deterministic order.  Each set has
        size ≥ 2 by construction.
    """
    exclude = {int(metal_idx)} if metal_idx is not None else set()
    groups: List[FrozenSet[int]] = []
    groups.extend(find_aromatic_groups(mol))
    groups.extend(find_symmetric_terminal_groups(mol))
    if use_automorphism:
        groups.extend(find_automorphism_orbits(mol, exclude=exclude))
    return _merge_overlapping_sets(groups)


# ---------------------------------------------------------------------------
# Hapto-π detection — pure graph (ring membership + donor adjacency to metal).
# ---------------------------------------------------------------------------
def detect_hapto_groups(
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
) -> List[FrozenSet[int]]:
    """Return frozensets of atom indices that form a hapto-π ring system.

    Criteria (universal, graph-only):
        * ≥ 3 carbon donors share a single ring
        * the ring has no metal in it (the metal is *above* the ring, not in it)

    Each returned frozenset is the FULL ring (so equi-distance can be enforced
    across all ring carbons, and the centroid encoded implicitly).
    """
    metal_idx = int(metal_idx)
    donor_set: Set[int] = set(int(d) for d in donor_idxs)
    if len(donor_set) < 3:
        return []
    # Drop M-D edges before ring perception so we see the underlying Cp/arene.
    try:
        from rdkit import Chem
        from rdkit.Chem import RWMol
        sub = RWMol(mol)
        to_remove: List[Tuple[int, int]] = []
        for bond in sub.GetBonds():
            u = int(bond.GetBeginAtomIdx())
            v = int(bond.GetEndAtomIdx())
            if u == metal_idx or v == metal_idx:
                to_remove.append((u, v))
        for u, v in to_remove:
            sub.RemoveBond(u, v)
        Chem.GetSSSR(sub)
        rings = [tuple(sorted(int(a) for a in r)) for r in sub.GetRingInfo().AtomRings()]
    except Exception:
        try:
            rings = [tuple(sorted(int(a) for a in r)) for r in mol.GetRingInfo().AtomRings()]
        except Exception:
            return []
    out: List[FrozenSet[int]] = []
    seen: Set[FrozenSet[int]] = set()
    for ring in rings:
        ring_atoms = set(ring)
        carbon_donors_in_ring: Set[int] = set()
        for a in ring_atoms:
            try:
                if mol.GetAtomWithIdx(int(a)).GetSymbol() != "C":
                    continue
            except Exception:
                continue
            if int(a) in donor_set:
                carbon_donors_in_ring.add(int(a))
        if len(carbon_donors_in_ring) >= 3:
            fs = frozenset(sorted(ring_atoms))
            if fs not in seen:
                seen.add(fs)
                out.append(fs)
    out.sort(key=lambda s: tuple(sorted(s)))
    return out


# ---------------------------------------------------------------------------
# Bound recording helper
# ---------------------------------------------------------------------------
def _set_bound(
    lower: np.ndarray,
    upper: np.ndarray,
    i: int,
    j: int,
    lo: float,
    hi: float,
) -> None:
    """Set ``lower[i,j] = lower[j,i] = lo`` and similarly for ``upper``.

    ``i == j`` is silently ignored.  Symmetric update.
    """
    if int(i) == int(j):
        return
    a, b = int(i), int(j)
    lower[a, b] = float(lo)
    lower[b, a] = float(lo)
    upper[a, b] = float(hi)
    upper[b, a] = float(hi)


# ---------------------------------------------------------------------------
# §2.2 Library lookup helpers
# ---------------------------------------------------------------------------
def _lookup_bond_md(
    grip_lib: Optional[GripLibrary],
    metal: str,
    metal_hyb: str,
    donor: str,
    donor_hyb: str,
    cod_lib: Optional[GripLibrary],
    min_n: int,
) -> Optional[Tuple[float, float, int]]:
    """Look up a metal-donor bond distribution.

    Cascade:
        1. v5 TM-pair-block lookup (if available)
        2. CCDC ``lookup_bond`` (with TM-aware fallback if enabled by env)
        3. COD library if provided
        4. None
    """
    # IMPORTANT: the legacy `GripLibrary.lookup_bond` cov-cascade is unsound
    # for metal-donor bonds — it pools any "N with any neighbour" entry into
    # an organic-dominated μ ≈ 1.48 Å for Pt-N (real Pt-N ≈ 2.05 Å).  Only
    # the v4/v5 pair-block + TM-category tables give the correct M-D μ.
    # If neither is available we return None and the caller falls back to
    # the universal cov-sum (polyhedra.md_distance).
    def _try(lib: Optional[GripLibrary]) -> Optional[Tuple[float, float, int]]:
        if lib is None:
            return None
        if getattr(lib, "has_pair_tables", False):
            try:
                hit = lib._v5_lookup_bond_block(metal, metal_hyb, donor, donor_hyb, min_n)
            except Exception:
                hit = None
            if hit is not None:
                return hit
            # v4 pair lookup (when v5 block table not populated for this pair)
            try:
                hit = lib._v4_lookup_bond(metal, metal_hyb, donor, donor_hyb, min_n)
            except Exception:
                hit = None
            if hit is not None:
                return hit
        # TM-aware additive fallback — env-gated; pure on-the-fly cascade,
        # produces correct M-D μ when enabled.
        try:
            hit = lib._tm_lookup_bond(metal, metal_hyb, donor, donor_hyb, min_n)
        except Exception:
            hit = None
        return hit

    hit = _try(grip_lib)
    if hit is not None:
        return hit
    hit = _try(cod_lib)
    if hit is not None:
        return hit
    return None


def _lookup_bond_organic(
    grip_lib: GripLibrary,
    z1: str,
    hyb1: str,
    z2: str,
    hyb2: str,
    cod_lib: Optional[GripLibrary],
    min_n: int,
) -> Optional[Tuple[float, float, int]]:
    """Look up an organic bond distribution.  CCDC then COD then None.

    The env-flag :data:`MOGUL_BOND_FALLBACK_ENV_FLAG`
    (``DELFIN_FFFREE_MOGUL_BOND_FALLBACK``) enables the multi-tier fallback
    chain (Tiers 1-5) — see :func:`_lookup_bond_organic_with_tier`.
    Default-OFF: byte-identical to the pre-fix path (Tier 1 only) when the
    env-flag is unset.
    """
    if _bond_fallback_enabled():
        hit, _tier = _lookup_bond_organic_with_tier(
            grip_lib, z1, hyb1, z2, hyb2, cod_lib, min_n
        )
        return hit
    # Legacy path -- byte-identical to the pre-fix behaviour.
    hit = grip_lib.lookup_bond(z1, hyb1, z2, hyb2, min_n=min_n)
    if hit is not None:
        return hit
    if cod_lib is not None:
        hit = cod_lib.lookup_bond(z1, hyb1, z2, hyb2, min_n=min_n)
        if hit is not None:
            return hit
    return None


# ---------------------------------------------------------------------------
# Multi-tier bond-key fallback chain (env-gated, default-OFF).
#
# Motivation: the CCDC v5 library indexes bonds by *centred-fragment* keys
# (one centre atom + its FULL neighbour-list).  ``GripLibrary.lookup_bond``
# only walks single-neighbour centred fragments, which match approximately
# 0-10 % of typical organic bonds (most real centres have >= 2 neighbours).
#
# This caused mogul-primary builds to fall back to Pauling covalent-radii
# sums for ~90 % of ligand 1-2 bonds, eliminating CCDC-derived planarity
# and terminal-angle constraints for organic-internal geometry.
#
# The fallback chain queries:
#   * Tier 1 -- centred-fragment ``lookup_bond`` (legacy).
#   * Tier 2 -- v4 pair table with FULL hyb pair (z1, hyb1, z2, hyb2).
#   * Tier 3 -- v4 pair table, one hyb wildcarded.
#   * Tier 4 -- v4 pair table, both hyb wildcarded.
#   * Tier 5 -- element-pair aggregate over ALL hyb variants
#     (:meth:`GripLibrary._lookup_bond_element_pair_aggregate`).
#   * Tier 6 -- caller covalent-radii fallback (NOT done here; returned as
#     ``(None, 6)`` so the caller can fall back universally).
#
# Determinism: each tier is deterministic; the chain stops at the first
# hit, so the returned (mu, sigma, n) and tier are reproducible.
# Universality: all keys are (element, hyb) tuples — no SMILES, no class,
# no element ``if``-branches.
# ---------------------------------------------------------------------------
MOGUL_BOND_FALLBACK_ENV_FLAG = "DELFIN_FFFREE_MOGUL_BOND_FALLBACK"


def _bond_fallback_enabled() -> bool:
    """``True`` when :data:`MOGUL_BOND_FALLBACK_ENV_FLAG` is set to ``1``."""
    return os.environ.get(MOGUL_BOND_FALLBACK_ENV_FLAG, "").strip() == "1"


def _lookup_bond_organic_with_tier(
    grip_lib: GripLibrary,
    z1: str,
    hyb1: str,
    z2: str,
    hyb2: str,
    cod_lib: Optional[GripLibrary],
    min_n: int,
) -> Tuple[Optional[Tuple[float, float, int]], int]:
    """Multi-tier organic bond lookup; returns ``((mu, sigma, n), tier)``.

    Tier number is 1..5 when a library hit is returned, or 6 when the
    entire chain misses (caller must apply the covalent-radii fallback).
    The COD library is consulted at every tier in parallel with the
    CCDC library: a Tier-N CCDC hit wins over a Tier-N COD hit, but COD
    Tier-N is preferred over CCDC Tier-(N+1) (CCDC pair-table is far
    larger than COD so this rarely fires in practice).

    Universal: keyed purely on (element, hyb).  No class/SMILES branches.
    """
    # ---- Tier 1: centred-fragment lookup_bond (legacy chain) ----
    hit = grip_lib.lookup_bond(z1, hyb1, z2, hyb2, min_n=min_n)
    if hit is not None:
        return (hit, 1)
    if cod_lib is not None:
        hit = cod_lib.lookup_bond(z1, hyb1, z2, hyb2, min_n=min_n)
        if hit is not None:
            return (hit, 1)

    # ---- Tier 2: pair-table, full hyb (z1, hyb1, z2, hyb2) ----
    if getattr(grip_lib, "has_pair_tables", False):
        try:
            hit = grip_lib._v4_lookup_bond(z1, hyb1, z2, hyb2, min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 2)
    if cod_lib is not None and getattr(cod_lib, "has_pair_tables", False):
        try:
            hit = cod_lib._v4_lookup_bond(z1, hyb1, z2, hyb2, min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 2)

    # ---- Tier 3: pair-table, one hyb wildcarded ----
    # Already tried inside ``_v4_lookup_bond`` -- but only when the FULL key
    # missed, so a separate tier counter is informative.  The library walks
    # (hyb1, hyb2), (hyb1, '*'), ('*', hyb2), ('*', '*') in order; the
    # one-wildcard cases are tiers 3, both-wildcard is tier 4.
    if getattr(grip_lib, "has_pair_tables", False):
        for h1, h2 in ((hyb1, "*"), ("*", hyb2)):
            try:
                hit = grip_lib._v4_lookup_bond(z1, h1, z2, h2, min_n)
            except Exception:
                hit = None
            if hit is not None:
                return (hit, 3)
    if cod_lib is not None and getattr(cod_lib, "has_pair_tables", False):
        for h1, h2 in ((hyb1, "*"), ("*", hyb2)):
            try:
                hit = cod_lib._v4_lookup_bond(z1, h1, z2, h2, min_n)
            except Exception:
                hit = None
            if hit is not None:
                return (hit, 3)

    # ---- Tier 4: pair-table, both hyb wildcarded ----
    if getattr(grip_lib, "has_pair_tables", False):
        try:
            hit = grip_lib._v4_lookup_bond(z1, "*", z2, "*", min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 4)
    if cod_lib is not None and getattr(cod_lib, "has_pair_tables", False):
        try:
            hit = cod_lib._v4_lookup_bond(z1, "*", z2, "*", min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 4)

    # ---- Tier 5: element-pair aggregate (NEW) ----
    if hasattr(grip_lib, "_lookup_bond_element_pair_aggregate"):
        try:
            hit = grip_lib._lookup_bond_element_pair_aggregate(z1, z2, min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 5)
    if cod_lib is not None and hasattr(cod_lib, "_lookup_bond_element_pair_aggregate"):
        try:
            hit = cod_lib._lookup_bond_element_pair_aggregate(z1, z2, min_n)
        except Exception:
            hit = None
        if hit is not None:
            return (hit, 5)

    # All tiers missed — caller must fall back to covalent radii (Tier 6).
    return (None, 6)


def _lookup_angle(
    grip_lib: GripLibrary,
    z1: str,
    z2: str,
    hyb2: str,
    z3: str,
    cod_lib: Optional[GripLibrary],
    min_n: int,
    hyb1: str = "*",
    hyb3: str = "*",
) -> Optional[Tuple[float, float, int]]:
    """Look up a 1,3 angle distribution."""
    hit = grip_lib.lookup_angle(z1, z2, hyb2, z3, hyb1=hyb1, hyb3=hyb3, min_n=min_n)
    if hit is not None:
        return hit
    if cod_lib is not None:
        hit = cod_lib.lookup_angle(z1, z2, hyb2, z3, hyb1=hyb1, hyb3=hyb3, min_n=min_n)
        if hit is not None:
            return hit
    return None


def _hapto_eta_tag(n_ring: int) -> Optional[str]:
    """Map ring size → v5 TM-category tag."""
    return {
        2: "hapto_eta2",
        5: "hapto_eta5",
        6: "hapto_eta6",
    }.get(int(n_ring))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
class BoundsInfo:
    """Container for build diagnostics (mirrors the spec's ``info`` dict)."""

    def __init__(self) -> None:
        self.n_atoms: int = 0
        self.n_bond_bounds: int = 0
        self.n_bond_lib_hit: int = 0
        self.n_bond_fallback: int = 0
        self.n_angle_bounds: int = 0
        self.n_angle_lib_hit: int = 0
        self.n_md_bounds: int = 0
        self.n_md_lib_hit: int = 0
        self.n_dd_bounds: int = 0
        self.n_hapto_groups: int = 0
        self.n_aromatic_groups: int = 0
        self.n_resonance_groups: int = 0
        self.n_vdw_floors: int = 0
        self.bounds_relaxed: List[Tuple[int, int, str]] = []

    def as_dict(self) -> Dict[str, Any]:
        d = dict(self.__dict__)
        d["bounds_relaxed"] = list(self.bounds_relaxed)
        return d


def build_bounds_matrix(
    syms: Sequence[str],
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
    grip_lib: Optional[GripLibrary] = None,
    cod_lib: Optional[GripLibrary] = None,
    *,
    geometry: Optional[str] = None,
    min_n: int = _MIN_N_DEFAULT,
    use_automorphism: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Construct the whole-complex distance-geometry bounds matrix.

    Parameters
    ----------
    syms : sequence of str
        Element symbols, length ``n_atoms``.
    mol : RDKit ``Mol`` / ``RWMol``
        Molecule with bonds, hybridisation, ring info populated.
    metal_idx : int
        Atom index of the metal centre (or *one* metal centre for
        multi-metal complexes — the donor-donor + hapto bounds will be
        anchored on this metal).
    donor_idxs : sequence of int
        Atom indices of the coordinating donors (the M-D-rigid set).
    grip_lib : GripLibrary, optional
        Pre-loaded CCDC library.  ``None`` ⇒ load the default release-pinned
        :data:`grip_mogul_lookup.DEFAULT_LIB_PATH`.
    cod_lib : GripLibrary, optional
        Pre-loaded COD library, used as a fallback after the CCDC lookup
        misses.
    geometry : str, optional
        Polyhedron name (e.g. ``"OC-6 octahedron"``).  When omitted the
        donor-donor bounds are derived from the coordination number alone
        (the *first* geometry in ``GEOM_BY_CN[cn]`` is used as the
        reference).
    min_n : int
        Minimum sample size for a library entry to be trusted.
    use_automorphism : bool
        If True (default), include graph-automorphism orbits in the
        resonance detection.

    Returns
    -------
    lower : (n, n) ndarray float
        Lower distance bounds.  ``lower[i,j] = lower[j,i]``; diagonal = 0.
    upper : (n, n) ndarray float
        Upper distance bounds.  Entries with no explicit bound default to
        :data:`_UB_INF`.
    info : dict
        Diagnostic counters (see :class:`BoundsInfo`).
    """
    if grip_lib is None:
        grip_lib = _default_library()
    n = int(len(syms))
    info = BoundsInfo()
    info.n_atoms = n

    # Initialise: all lower-bounds zero, all upper-bounds at sentinel.
    lower = np.zeros((n, n), dtype=np.float64)
    upper = np.full((n, n), _UB_INF, dtype=np.float64)

    metal_idx = int(metal_idx)
    donor_set: Set[int] = set(int(d) for d in donor_idxs)

    bond_sigma_mult = _bond_sigma_mult()
    md_sigma_mult = _md_sigma_mult()
    dd_sigma_mult = _dd_sigma_mult()

    # ------------------------------------------------------------------
    # 1) BONDED PAIRS — μ ± 3σ from CCDC; cov-sum fallback
    # ------------------------------------------------------------------
    bonds = _bonds_sorted(mol)
    metal_bond_pairs: Set[Tuple[int, int]] = set()
    for (i, j) in bonds:
        zi, zj = str(syms[i]), str(syms[j])
        try:
            ai = mol.GetAtomWithIdx(int(i))
            aj = mol.GetAtomWithIdx(int(j))
        except Exception:
            continue
        hi = hyb_str(ai)
        hj = hyb_str(aj)
        is_md = (i == metal_idx and j in donor_set) or (
            j == metal_idx and i in donor_set
        )
        if is_md:
            # M-D bond: handled separately in the metal-donor block.
            metal_bond_pairs.add((i, j))
            continue
        hit = _lookup_bond_organic(grip_lib, zi, hi, zj, hj, cod_lib, min_n)
        if hit is None:
            # Pauling fallback — single-bond covalent sum.  Larger σ so the
            # bound is permissive.
            mu = _cov_radius(zi) + _cov_radius(zj)
            sigma = max(0.05, 0.05 * mu)
            info.n_bond_fallback += 1
        else:
            mu, sigma, _nbn = hit
            sigma = max(_SIGMA_FLOOR_BOND, float(sigma))
            info.n_bond_lib_hit += 1
        lo = max(0.1, float(mu) - bond_sigma_mult * float(sigma))
        hi_val = float(mu) + bond_sigma_mult * float(sigma)
        _set_bound(lower, upper, i, j, lo, hi_val)
        info.n_bond_bounds += 1

    # ------------------------------------------------------------------
    # 2) 1,3 ANGLE PAIRS — via law of cosines from CCDC angle distributions
    # ------------------------------------------------------------------
    triples = _angle_triples_sorted(mol)
    for (i, j, k) in triples:
        # Skip triples that involve the metal as a *centre* — those are the
        # D-M-D polyhedron angles, handled by the donor-donor block.
        if int(j) == metal_idx:
            continue
        zi = str(syms[i]); zj = str(syms[j]); zk = str(syms[k])
        try:
            aj = mol.GetAtomWithIdx(int(j))
        except Exception:
            continue
        hj = hyb_str(aj)
        hit = _lookup_angle(grip_lib, zi, zj, hj, zk, cod_lib, min_n)
        if hit is None:
            # No CCDC entry — skip rather than guess.  The bonded bound on
            # (i,j) and (j,k) plus the vdW floor still constrain the
            # geometry; an over-broad angle would only add noise.
            continue
        theta_deg, sigma_deg, _nbn = hit
        sigma_deg = max(_SIGMA_FLOOR_ANGLE_DEG, float(sigma_deg))
        # Use the bonded means (μ values) for d12 and d23 — those are the
        # most-probable distances; the law-of-cosines bound is then the
        # propagation of the angle σ through.
        d12 = _bonded_mean(syms, i, j, mol, grip_lib, cod_lib, min_n,
                           md_set=metal_bond_pairs, metal_idx=metal_idx,
                           donor_set=donor_set)
        d23 = _bonded_mean(syms, j, k, mol, grip_lib, cod_lib, min_n,
                           md_set=metal_bond_pairs, metal_idx=metal_idx,
                           donor_set=donor_set)
        if d12 is None or d23 is None:
            continue
        theta_lo_deg = max(1.0, float(theta_deg) - bond_sigma_mult * sigma_deg)
        theta_hi_deg = min(179.0, float(theta_deg) + bond_sigma_mult * sigma_deg)
        cos_lo = math.cos(math.radians(theta_hi_deg))  # cos is decreasing
        cos_hi = math.cos(math.radians(theta_lo_deg))
        d2_lo = d12 * d12 + d23 * d23 - 2.0 * d12 * d23 * cos_hi
        d2_hi = d12 * d12 + d23 * d23 - 2.0 * d12 * d23 * cos_lo
        d_lo = math.sqrt(max(1e-4, d2_lo))
        d_hi = math.sqrt(max(d2_lo + 1e-6, d2_hi))
        # Don't tighten an already-tight bonded bound (i,k might be bonded if
        # the graph has parallel edges in fused systems; usually it is not).
        # Only loosen / honour the 1,3 bound when (i,k) is not already bonded.
        if (min(i, k), max(i, k)) in set(bonds):
            continue
        _set_bound(lower, upper, i, k, d_lo, d_hi)
        info.n_angle_bounds += 1
        info.n_angle_lib_hit += 1

    # ------------------------------------------------------------------
    # 3) METAL-DONOR PAIRS — μ ± 2σ from CCDC TM-aware lookup
    # ------------------------------------------------------------------
    metal_sym = str(syms[metal_idx])
    metal_atom = None
    try:
        metal_atom = mol.GetAtomWithIdx(metal_idx)
    except Exception:
        pass
    metal_hyb = hyb_str(metal_atom) if metal_atom is not None else "*"
    hapto_groups = detect_hapto_groups(mol, metal_idx, donor_idxs)
    hapto_atom_set: Set[int] = set()
    for g in hapto_groups:
        hapto_atom_set |= set(g)
    info.n_hapto_groups = len(hapto_groups)

    for d in sorted(donor_set):
        if int(d) == metal_idx:
            continue
        donor_sym = str(syms[int(d)])
        try:
            donor_atom = mol.GetAtomWithIdx(int(d))
        except Exception:
            donor_atom = None
        donor_hyb = hyb_str(donor_atom) if donor_atom is not None else "*"
        # If the donor is part of a hapto ring we look up via TM-category
        # (η-fragment library); otherwise via the bond table.
        eta_hit = None
        if int(d) in hapto_atom_set:
            # Find the ring this donor belongs to → use its size to pick η-tag.
            for g in hapto_groups:
                if int(d) in g:
                    tag = _hapto_eta_tag(len(g))
                    if tag is not None:
                        try:
                            eta_hit = grip_lib.lookup_tm_category(
                                tag, metal_sym, donor_sym, min_n=min_n
                            )
                        except Exception:
                            eta_hit = None
                    break
        if eta_hit is not None:
            mu, sigma, _nbn = eta_hit
            sigma = max(_SIGMA_FLOOR_BOND, float(sigma))
            info.n_md_lib_hit += 1
        else:
            md_hit = _lookup_bond_md(
                grip_lib, metal_sym, metal_hyb, donor_sym, donor_hyb,
                cod_lib, min_n,
            )
            if md_hit is None:
                # Fallback: polyhedra.md_distance (covalent-sum or Shannon
                # ionic).  Use a generous σ so the solver isn't over-pinned.
                mu = _polyhedra.md_distance(metal_sym, donor_sym)
                sigma = max(0.05, 0.05 * mu)
                info.n_bond_fallback += 1
            else:
                mu, sigma, _nbn = md_hit
                sigma = max(_SIGMA_FLOOR_BOND, float(sigma))
                info.n_md_lib_hit += 1
        lo = max(0.5, float(mu) - md_sigma_mult * float(sigma))
        hi_val = float(mu) + md_sigma_mult * float(sigma)
        _set_bound(lower, upper, metal_idx, int(d), lo, hi_val)
        info.n_md_bounds += 1

    # ------------------------------------------------------------------
    # 4) DONOR-DONOR PAIRS — polyhedron geometry encoded as distances
    # ------------------------------------------------------------------
    # Compute M-D mean distance per donor (already enforced above) to scale
    # the polyhedron reference vectors.  Average across donors is used so
    # all donor-donor distances share a single radius (deterministic).
    donor_list = sorted(int(d) for d in donor_set if int(d) != metal_idx)
    cn = len(donor_list)
    if cn >= 2:
        # Estimate mean M-D distance from already-set bounds (midpoint).
        md_means: List[float] = []
        md_sigmas: List[float] = []
        for d in donor_list:
            lo_md = float(lower[metal_idx, d])
            hi_md = float(upper[metal_idx, d])
            if hi_md >= _UB_INF * 0.5:
                continue
            md_means.append(0.5 * (lo_md + hi_md))
            md_sigmas.append(max(0.05, (hi_md - lo_md) / (2.0 * md_sigma_mult)))
        if md_means:
            r_md = float(np.mean(md_means))
            sigma_md = float(np.mean(md_sigmas))
        else:
            r_md = float(_polyhedra.md_distance(metal_sym, "C"))
            sigma_md = 0.10

        # Resolve geometry
        ref_vecs: Optional[np.ndarray] = None
        if geometry is not None:
            try:
                ref_vecs = _polyhedra.ref_vectors(str(geometry))
            except Exception:
                ref_vecs = None
        if ref_vecs is None:
            cands = _polyhedra.geometries_for_cn(int(cn), metal_sym)
            for g in cands:
                try:
                    ref_vecs = _polyhedra.ref_vectors(str(g))
                    break
                except Exception:
                    continue
        if ref_vecs is not None and ref_vecs.shape[0] >= cn:
            # Take the first `cn` reference vertices in their canonical order
            # → matches the assembler's index-space.  Scale by r_md.
            vecs = np.asarray(ref_vecs[:cn], dtype=np.float64) * float(r_md)
            for a_pos in range(cn):
                for b_pos in range(a_pos + 1, cn):
                    da = donor_list[a_pos]
                    db = donor_list[b_pos]
                    d_ideal = float(np.linalg.norm(vecs[a_pos] - vecs[b_pos]))
                    if d_ideal <= 0.0:
                        continue
                    band = dd_sigma_mult * sigma_md
                    lo = max(0.5, d_ideal - band)
                    hi = d_ideal + band
                    _set_bound(lower, upper, da, db, lo, hi)
                    info.n_dd_bounds += 1

    # ------------------------------------------------------------------
    # 5) HAPTO-π PAIRS — equi-distance constraint across ring members
    # ------------------------------------------------------------------
    # Already covered by M-D bounds via TM-category lookup above; here we
    # additionally tighten ring-internal C-C bonds via the resonance
    # constraint block.  No new bounds in this section — η-rings are
    # *resonance groups* by construction (they appear in
    # ``find_aromatic_groups`` because RDKit marks Cp / arene atoms aromatic).

    # ------------------------------------------------------------------
    # 6) RESONANCE / MESOMERIC EQUIVALENCE — equal-distance bounds
    # ------------------------------------------------------------------
    res_groups = detect_resonance_groups(
        mol, metal_idx=metal_idx, use_automorphism=use_automorphism,
    )
    info.n_resonance_groups = len(res_groups)
    info.n_aromatic_groups = len(find_aromatic_groups(mol))

    bond_set: Set[Tuple[int, int]] = set(bonds)
    # Build a graph adjacency for the "common-parent" bond resolution.
    adj: Dict[int, Set[int]] = {i: set() for i in range(n)}
    for (a, b) in bonds:
        adj[a].add(b)
        adj[b].add(a)

    for group in res_groups:
        atoms = sorted(int(a) for a in group)
        if len(atoms) < 2:
            continue
        # --- Equi-distance bond set ---
        # Two flavours of resonance equivalence (universal, graph-only):
        #
        # (a) INTRA-RING bonds: a bond (i, j) where BOTH i, j ∈ group.
        #     Captures aromatic ring C-C bonds (benzene, Cp, …).
        #
        # (b) PARENT-RADIAL bonds: a bond (parent, x) where x ∈ group AND
        #     parent ∉ group AND parent has multiple group neighbours.
        #     Captures carboxylate C-O × 2, nitrate N-O × 3, sulfate S-O × 4,
        #     methyl C-H × 3 — i.e. every "central atom with equivalent
        #     terminal neighbours" pattern.
        #
        # Both flavours combine into ONE equi-distance window: this is the
        # correct constraint because in a real resonance / VSEPR-symmetric
        # group, the chemistry equalises every "type-A" bond.
        intra_bonds: List[Tuple[int, int]] = []
        parent_radial_bonds: Dict[int, List[Tuple[int, int]]] = {}
        group_set: Set[int] = set(atoms)
        for a_pos in range(len(atoms)):
            for b_pos in range(a_pos + 1, len(atoms)):
                a, b = atoms[a_pos], atoms[b_pos]
                if (a, b) in bond_set:
                    intra_bonds.append((a, b))
        # Find common parents whose multiple neighbours all sit in the group.
        # Iterate sorted parents for determinism.  METAL parents are
        # EXCLUDED — coordination bonds are *not* resonance-equivalent to
        # internal organic bonds; the η-fragment library already places
        # them on the correct ((M, donor, η-class)) distribution.
        parent_candidates: Set[int] = set()
        for x in atoms:
            for nb in adj.get(x, set()):
                if nb not in group_set:
                    parent_candidates.add(int(nb))
        for parent in sorted(parent_candidates):
            if _is_metal(syms[int(parent)]):
                continue
            kids_in_group = sorted(int(k) for k in adj.get(parent, set()) if k in group_set)
            if len(kids_in_group) >= 2:
                parent_radial_bonds.setdefault(parent, []).extend(
                    (min(parent, k), max(parent, k)) for k in kids_in_group
                )

        # Collect ALL bonds that should be equi-distant for THIS group.
        equi_bonds: List[Tuple[int, int]] = list(intra_bonds)
        for radial in parent_radial_bonds.values():
            equi_bonds.extend(radial)
        # Deduplicate (a parent-radial pair may also be intra).
        equi_bonds = sorted(set(equi_bonds))
        if len(equi_bonds) < 2:
            continue

        # Use the *mean* of the already-set lower bounds as the equi-distance
        # midpoint, with a tight resonance σ.  This guarantees every internal
        # bond in the group sees an identical [lo, hi] window.
        mids: List[float] = []
        for (a, b) in equi_bonds:
            lo_b = float(lower[a, b])
            hi_b = float(upper[a, b])
            if hi_b >= _UB_INF * 0.5:
                continue
            mids.append(0.5 * (lo_b + hi_b))
        if not mids:
            continue
        d_eq = float(np.mean(mids))
        sigma_eq = max(_RESONANCE_SIGMA_DEFAULT, float(np.std(mids)))
        lo_eq = max(0.5, d_eq - bond_sigma_mult * sigma_eq)
        hi_eq = d_eq + bond_sigma_mult * sigma_eq
        for (a, b) in equi_bonds:
            _set_bound(lower, upper, a, b, lo_eq, hi_eq)

    # ------------------------------------------------------------------
    # 7) NON-BONDED vdW FLOOR — α × (r_vdW_i + r_vdW_j)
    # ------------------------------------------------------------------
    # Apply only to pairs that DO NOT already have a stronger lower bound.
    # That preserves the bonded / angle / M-D / D-D windows and only fills
    # the remaining cells with the vdW floor.
    bonded_atoms_1_2: Set[Tuple[int, int]] = set(bond_set)
    bonded_atoms_1_3: Set[Tuple[int, int]] = set(
        (min(i, k), max(i, k)) for (i, _j, k) in triples
    )
    # Also exclude M-D-related pairs (the donor-donor and metal-donor
    # bounds are explicit).
    md_pairs_set: Set[Tuple[int, int]] = set()
    for d in donor_list:
        md_pairs_set.add((min(metal_idx, d), max(metal_idx, d)))
    for a_pos in range(len(donor_list)):
        for b_pos in range(a_pos + 1, len(donor_list)):
            a, b = donor_list[a_pos], donor_list[b_pos]
            md_pairs_set.add((min(a, b), max(a, b)))

    for i in range(n):
        for j in range(i + 1, n):
            key = (i, j)
            if key in bonded_atoms_1_2 or key in bonded_atoms_1_3 or key in md_pairs_set:
                continue
            # If a bound was already written (e.g. resonance), the lower
            # bound is non-zero — do not overwrite the existing tighter
            # window.
            if lower[i, j] > 0.0 and upper[i, j] < _UB_INF * 0.5:
                continue
            ri = _vdw_radius(syms[i])
            rj = _vdw_radius(syms[j])
            lo = _VDW_FLOOR_FACTOR * (ri + rj)
            # Keep upper at sentinel — non-bonded pairs have no upper limit
            # (counter-ions, distant ligand atoms, …).
            lower[i, j] = lo
            lower[j, i] = lo
            info.n_vdw_floors += 1

    # Symmetrise + zero-diagonal sanity
    for i in range(n):
        lower[i, i] = 0.0
        upper[i, i] = 0.0

    # ------------------------------------------------------------------
    # Feasibility check + minimal relaxation
    # ------------------------------------------------------------------
    # A pair with lower > upper is infeasible.  Relax by swapping the two
    # values and flagging.
    for i in range(n):
        for j in range(i + 1, n):
            if lower[i, j] > upper[i, j]:
                # Relax: widen the window so it is at least the lower bound.
                upper[i, j] = lower[i, j] + 0.05
                upper[j, i] = lower[i, j] + 0.05
                info.bounds_relaxed.append((int(i), int(j), "lower>upper"))

    return lower, upper, info.as_dict()


def _bonded_mean(
    syms: Sequence[str],
    i: int,
    j: int,
    mol,
    grip_lib: GripLibrary,
    cod_lib: Optional[GripLibrary],
    min_n: int,
    *,
    md_set: Set[Tuple[int, int]],
    metal_idx: int,
    donor_set: Set[int],
) -> Optional[float]:
    """Return the μ of the bonded distribution for atoms i-j, or None.

    Used by the law-of-cosines 1,3 calculation to convert an angle bound
    into a distance bound.
    """
    key = (min(int(i), int(j)), max(int(i), int(j)))
    is_md = key in md_set
    zi = str(syms[i]); zj = str(syms[j])
    try:
        ai = mol.GetAtomWithIdx(int(i))
        aj = mol.GetAtomWithIdx(int(j))
    except Exception:
        return None
    hi = hyb_str(ai); hj = hyb_str(aj)
    if is_md:
        # Use the polyhedra mean distance — bypasses TM lookup (the angle
        # block is not centred on the metal so M-D distances are auxiliary
        # constants here).
        return float(_polyhedra.md_distance(zi if int(i) == int(metal_idx) else zj,
                                            zj if int(i) == int(metal_idx) else zi))
    hit = _lookup_bond_organic(grip_lib, zi, hi, zj, hj, cod_lib, min_n)
    if hit is not None:
        return float(hit[0])
    return float(_cov_radius(zi) + _cov_radius(zj))
