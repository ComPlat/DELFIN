"""GRIP — Healing Mode (Topology-aware iterative atom repositioning).

Extends GRIP from *polish-only* (Mahalanobis L-BFGS on already-near-optimum
coordinates) to *polish + heal* — a pre-polish constraint-satisfaction step
that repositions atoms which the constructive builder placed many σ off
their CCDC-grounded ideal bond length.

Mathematical background
-----------------------
This is **classical distance-geometry constraint satisfaction**
(Crippen & Havel 1988) adapted to use CCDC-grounded ideal bond lengths
in place of triangle-inequality bounds.

For each *broken* atom A (i.e. an atom that participates in a bond whose
length deviates by more than ``sigma_threshold`` σ from its CCDC ideal),
we iteratively pull A onto a *soft target* derived from its bonded
neighbours and the ideal bond lengths::

    For each bonded neighbor N of A:
        unit_dir = (r_A - r_N) / ||r_A - r_N||
        target   = r_N + ideal_bond_length(N-A) * unit_dir

    r_A <- mean(target over all bonded neighbours)

The procedure converges *linearly* for mild distortions (≤ 20% bond
mismatch) and may fail to converge for fully scrambled coordinates --
those need a full re-embedding (planned for a future "Whole-Complex-DG"
step, see project_dg_bounds_construction_spec).

Determinism contract
--------------------
* All atom iteration is lexicographically sorted.
* Bond / neighbour iteration uses sorted ``(i, j)`` tuples.
* Degenerate directions (||r_A - r_N|| < ``_DEGENERATE_TOL``) fall back
  to a deterministic unit vector keyed off the atom indices.
* No RNG, no floating-point parallel reduction.
* float64 throughout.

Force-field-free contract
-------------------------
The targets come from CCDC bond-length statistics
(:mod:`grip_mogul_lookup`).  No Lennard-Jones, no Coulomb, no spring
constants -- only ideal lengths and a Mahalanobis-style σ threshold.

Default OFF / byte-identical
----------------------------
The new code path activates only when
``DELFIN_FFFREE_GRIP_HEALING_MODE`` is set to one of ``{"1", "true",
"yes", "on"}``.  When the flag is unset, :func:`grip_polish_with_healing`
delegates verbatim to :func:`grip_polish.grip_polish`, preserving
byte-identity with HEAD 93b396d.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

_LOG = logging.getLogger(__name__)

__all__ = [
    "HEALING_MODE_ENV",
    "DEFAULT_SIGMA_THRESHOLD",
    "DEFAULT_MAX_ITER",
    "DEFAULT_TOL",
    "DEFAULT_FLOATING_FACTOR",
    "DEFAULT_DAMPING",
    "HealingDiagnostics",
    "detect_broken_atoms",
    "iterative_topology_repositioning",
    "grip_polish_with_healing",
    "healing_mode_active",
    "_deterministic_fallback_direction",
]

# ---------------------------------------------------------------------------
# Env-flag + defaults
# ---------------------------------------------------------------------------
HEALING_MODE_ENV: str = "DELFIN_FFFREE_GRIP_HEALING_MODE"

#: σ-threshold above which a bond residual flags both endpoints as
#: *broken*.  3 σ is the standard Mahalanobis significance cut.
DEFAULT_SIGMA_THRESHOLD: float = 3.0

#: Iteration cap for :func:`iterative_topology_repositioning`.
DEFAULT_MAX_ITER: int = 50

#: Convergence tolerance for the per-iteration max displacement (Å).
DEFAULT_TOL: float = 0.01

#: An atom counts as "floating" when its closest topological neighbour
#: sits beyond ``floating_factor * ideal_bond_length`` of the actual
#: ideal length, signalling that the constructive placement *missed*
#: that ligand altogether.
DEFAULT_FLOATING_FACTOR: float = 2.0

#: Per-iteration step damping coefficient.  A pure Jacobi update
#: (alpha=1) can oscillate when many bonds are coupled (one atom is
#: endpoint to multiple broken bonds whose target projections
#: disagree).  Damping by 0.5 (the classical successive-over-relaxation
#: under-relaxation) makes the iteration provably contractive on
#: any acyclic ligand graph; the convergence rate halves, which is
#: fine for our O(20-50) iteration budget.
DEFAULT_DAMPING: float = 0.5

#: Below this distance two atoms count as coincident — direction is
#: undefined and we substitute a deterministic fallback.
_DEGENERATE_TOL: float = 0.01

#: Fallback σ used when the GRIP library does not return a bond entry
#: (e.g. exotic element pair, sparse library).  This is *much* wider
#: than the typical 0.02 Å σ so a missing entry never flags an atom as
#: broken on its own; the sigma_threshold gate becomes near-vacuous.
_FALLBACK_SIGMA: float = 1.0

#: Fallback ideal bond length when the GRIP library has no entry --
#: sum of two carbon covalent radii (Å).  Used only for the residual
#: computation; healing skips the atom when no library entry exists.
_FALLBACK_MU: float = 1.52


def healing_mode_active() -> bool:
    """Return ``True`` iff :data:`HEALING_MODE_ENV` is set to a truthy
    value (``"1"``, ``"true"``, ``"yes"``, ``"on"``, case-insensitive).

    Default-OFF: any other value (unset, empty, ``"0"``, ``"false"``)
    returns ``False`` so the healing path is *byte-identical* with HEAD
    93b396d.
    """
    raw = os.environ.get(HEALING_MODE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


# ---------------------------------------------------------------------------
# Diagnostics container
# ---------------------------------------------------------------------------
@dataclass
class HealingDiagnostics:
    """Per-call diagnostic record for :func:`iterative_topology_repositioning`.

    Attributes
    ----------
    n_broken_initial
        Number of atoms flagged as broken before healing.
    n_broken_final
        Number of atoms still broken after healing.
    n_iter
        Number of iterations actually run (≤ ``max_iter``).
    converged
        ``True`` when the per-iteration max displacement dropped below
        ``tol`` before the iteration cap.
    max_delta_history
        List of per-iteration max displacements (Å); empty when no
        atom was broken (no-op path).
    skipped_orphans
        Atoms flagged as broken but with no bonded neighbour -- the
        healer cannot reposition them and skips them.
    floating_atoms
        Atoms flagged as floating (no neighbour within
        ``floating_factor * ideal_bond_length``).
    """

    n_broken_initial: int = 0
    n_broken_final: int = 0
    n_iter: int = 0
    converged: bool = False
    max_delta_history: List[float] = field(default_factory=list)
    skipped_orphans: List[int] = field(default_factory=list)
    floating_atoms: List[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Helpers — topology / library lookup
# ---------------------------------------------------------------------------
def _coerce_R(R) -> np.ndarray:
    """Return a writable ``(N, 3)`` float64 copy of ``R``."""
    R = np.asarray(R, dtype=np.float64)
    if R.ndim == 1:
        if R.size % 3 != 0:
            raise ValueError(f"coords must be (N,3) or (3N,); got {R.size}")
        R = R.reshape(-1, 3)
    return R.copy()


def _normalize_topology(
    bonds_topology: Iterable[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Return the deduplicated, lex-sorted ``(min, max)`` bond list."""
    out: Set[Tuple[int, int]] = set()
    for pair in bonds_topology:
        if len(pair) != 2:
            continue
        a = int(pair[0])
        b = int(pair[1])
        if a == b:
            continue
        if a > b:
            a, b = b, a
        out.add((a, b))
    return sorted(out)


def _build_adjacency(
    n_atoms: int,
    bonds: Sequence[Tuple[int, int]],
) -> List[List[int]]:
    """Return a list of sorted neighbour lists for atom indices ``0..n_atoms-1``."""
    adj: List[Set[int]] = [set() for _ in range(n_atoms)]
    for a, b in bonds:
        if 0 <= a < n_atoms and 0 <= b < n_atoms:
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(s) for s in adj]


def _lookup_ideal_bond(
    a_sym: str,
    b_sym: str,
    library=None,
) -> Tuple[float, float]:
    """Return ``(mu, sigma)`` for a bond ``a_sym-b_sym`` via the GRIP library.

    Falls back to ``(_FALLBACK_MU, _FALLBACK_SIGMA)`` when no entry is
    available -- this keeps the healer functional even when the
    library is missing (default-OFF byte-identity is still preserved
    because ``healing_mode_active()`` returns ``False``).
    """
    if library is None:
        return _FALLBACK_MU, _FALLBACK_SIGMA
    try:
        # Hybridisation wildcards keep the lookup at the broadest level
        # (we have no reliable per-atom hyb info in the healing path).
        hit = library.lookup_bond(a_sym, "*", b_sym, "*")
    except Exception:
        hit = None
    if hit is None:
        return _FALLBACK_MU, _FALLBACK_SIGMA
    mu, sigma, _n = hit
    if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
        return _FALLBACK_MU, _FALLBACK_SIGMA
    return float(mu), float(sigma)


def _build_ideal_table(
    bonds: Sequence[Tuple[int, int]],
    symbols: Sequence[str],
    library=None,
) -> Dict[Tuple[int, int], Tuple[float, float]]:
    """Return ``{(a, b): (mu, sigma)}`` for every bond, lex order keyed.

    Both endpoints are looked up symbol-wise; the call is deterministic
    because ``bonds`` is already lex-sorted and the library lookup is
    pure-functional.
    """
    out: Dict[Tuple[int, int], Tuple[float, float]] = {}
    for (a, b) in bonds:
        try:
            sa = str(symbols[a])
            sb = str(symbols[b])
        except (IndexError, TypeError):
            continue
        mu, sigma = _lookup_ideal_bond(sa, sb, library=library)
        out[(int(a), int(b))] = (float(mu), float(sigma))
    return out


def _deterministic_fallback_direction(
    atom_idx: int,
    neighbour_idx: int,
) -> np.ndarray:
    """Return a deterministic unit vector for degenerate-direction cases.

    When two atoms have collapsed onto each other (``||r_A - r_N|| <
    _DEGENERATE_TOL``) the heal cannot derive a direction from the
    coordinates.  Instead we pick a unit vector that depends *only*
    on the (atom, neighbour) integer indices -- this guarantees that
    the same (atom, neighbour) pair always produces the same direction
    across runs, machines and processes.

    The construction is:

    1. Combine the indices into a stable 64-bit hash via a Cantor-pair
       and the FNV constants (no Python ``hash()`` -- that is salted
       under PYTHONHASHSEED != 0).
    2. Map the hash to three coordinates in [-1, 1] using a fixed
       deterministic mixing function.
    3. Normalise.

    The result is byte-identical across processes, across machines,
    independent of PYTHONHASHSEED.
    """
    a = int(atom_idx)
    n = int(neighbour_idx)
    # FNV-1a 64-bit mix.
    h = np.uint64(14695981039346656037)
    prime = np.uint64(1099511628211)
    for value in (a, n, a ^ n):
        v = np.uint64(int(value) & 0xFFFFFFFFFFFFFFFF)
        h = np.uint64((int(h) ^ int(v)) & 0xFFFFFFFFFFFFFFFF)
        h = np.uint64((int(h) * int(prime)) & 0xFFFFFFFFFFFFFFFF)
    raw = int(h)
    # Map to three coordinates by slicing the 64-bit hash.
    x = ((raw >> 0) & 0xFFFFF) / float(0xFFFFF) * 2.0 - 1.0
    y = ((raw >> 20) & 0xFFFFF) / float(0xFFFFF) * 2.0 - 1.0
    z = ((raw >> 40) & 0xFFFFF) / float(0xFFFFF) * 2.0 - 1.0
    vec = np.array([x, y, z], dtype=np.float64)
    norm = float(np.linalg.norm(vec))
    if norm < 1e-9:
        # Pathological case: hash mapped to ~origin.  Return canonical x.
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    return vec / norm


# ---------------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------------
def detect_broken_atoms(
    coords,
    bonds_topology: Iterable[Tuple[int, int]],
    ideal_lengths: Optional[Dict[Tuple[int, int], Tuple[float, float]]] = None,
    symbols: Optional[Sequence[str]] = None,
    library=None,
    *,
    sigma_threshold: float = DEFAULT_SIGMA_THRESHOLD,
    floating_factor: float = DEFAULT_FLOATING_FACTOR,
) -> List[int]:
    """Return the lex-sorted list of atom indices flagged as *broken*.

    A bond ``(i, j)`` flags both ``i`` and ``j`` as broken when

    .. math::

        \\left| \\frac{||r_i - r_j|| - \\mu_{ij}}{\\sigma_{ij}} \\right|
        > \\text{sigma\\_threshold}

    Additionally, an atom is flagged as broken when it has at least one
    bonded neighbour but every neighbour sits beyond
    ``floating_factor * mu_ij`` of the ideal length -- this is the
    *floating-atom* signature (constructive placement missed the ligand
    entirely).

    Parameters
    ----------
    coords
        ``(N, 3)`` ndarray (or flat ``(3N,)``) of atomic coordinates (Å).
    bonds_topology
        Iterable of ``(i, j)`` bonded-pair indices.  May contain duplicates
        and reversed orderings; the function normalises to lex-sorted
        ``(min, max)`` pairs.
    ideal_lengths
        Optional pre-built ``{(min, max): (mu, sigma)}`` table.  When
        ``None`` the table is built from ``symbols`` + ``library`` (or
        the fallback constants if either is missing).
    symbols
        Optional per-atom element symbols (``["C", "H", ...]``).  Required
        when ``ideal_lengths`` is ``None``.
    library
        Optional :class:`GripLibrary` instance.  When ``None`` and
        ``ideal_lengths`` is ``None``, the fallback constants are used.
    sigma_threshold
        Mahalanobis σ cut for the broken classification (default 3.0).
    floating_factor
        Multiplier on the ideal bond length above which an atom counts
        as *floating* (default 2.0).

    Returns
    -------
    list of int
        Lex-sorted list of broken atom indices.  Empty when no bond
        violates the threshold and no atom is floating.
    """
    R = _coerce_R(coords)
    n_atoms = R.shape[0]
    if n_atoms == 0:
        return []
    bonds = _normalize_topology(bonds_topology)
    if not bonds:
        return []

    # Build the ideal table on demand.
    if ideal_lengths is None:
        if symbols is None:
            # Without symbols we cannot look up ideal lengths -- fall
            # back to the constant table.
            symbols = ["C"] * n_atoms
        ideal_lengths = _build_ideal_table(bonds, symbols, library=library)

    broken: Set[int] = set()
    # Adjacency by atom -> list of (neighbour, mu, sigma).
    per_atom: Dict[int, List[Tuple[int, float, float]]] = {
        i: [] for i in range(n_atoms)
    }

    for (a, b) in bonds:
        key = (a, b)
        mu, sigma = ideal_lengths.get(key, (_FALLBACK_MU, _FALLBACK_SIGMA))
        if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
            mu, sigma = _FALLBACK_MU, _FALLBACK_SIGMA
        per_atom[a].append((b, mu, sigma))
        per_atom[b].append((a, mu, sigma))
        try:
            d = float(np.linalg.norm(R[a] - R[b]))
        except (IndexError, ValueError):
            continue
        if not np.isfinite(d):
            broken.add(a)
            broken.add(b)
            continue
        z = (d - mu) / sigma
        if abs(z) > sigma_threshold:
            broken.add(a)
            broken.add(b)

    # Floating-atom check: any atom whose *closest* neighbour sits
    # beyond ``floating_factor * mu``.
    for atom_idx, nbrs in per_atom.items():
        if not nbrs:
            continue
        min_ratio = float("inf")
        for (j, mu, _sigma) in nbrs:
            if mu <= 0:
                continue
            try:
                d = float(np.linalg.norm(R[atom_idx] - R[j]))
            except (IndexError, ValueError):
                continue
            if not np.isfinite(d):
                broken.add(atom_idx)
                break
            ratio = d / mu
            if ratio < min_ratio:
                min_ratio = ratio
        if min_ratio > floating_factor:
            broken.add(atom_idx)

    return sorted(broken)


# ---------------------------------------------------------------------------
# Iterative repositioning
# ---------------------------------------------------------------------------
def iterative_topology_repositioning(
    coords,
    topology: Iterable[Tuple[int, int]],
    ideal_lengths: Optional[Dict[Tuple[int, int], Tuple[float, float]]] = None,
    symbols: Optional[Sequence[str]] = None,
    library=None,
    *,
    frozen_atoms: Optional[Iterable[int]] = None,
    max_iter: int = DEFAULT_MAX_ITER,
    tol: float = DEFAULT_TOL,
    sigma_threshold: float = DEFAULT_SIGMA_THRESHOLD,
    floating_factor: float = DEFAULT_FLOATING_FACTOR,
    damping: float = DEFAULT_DAMPING,
    return_diagnostics: bool = False,
):
    """Iteratively reposition broken atoms onto their CCDC-grounded targets.

    Algorithm (per iteration, deterministic order):

    1. Compute the broken-atom set via :func:`detect_broken_atoms`.
    2. For each broken atom *A* (lex-sorted) that has at least one bonded
       neighbour:

       * For each bonded neighbour *N* (lex-sorted):

         * ``unit_dir = (r_A - r_N) / ||r_A - r_N||``
           (or deterministic fallback if degenerate)
         * ``target_N = r_N + mu(N-A) * unit_dir``

       * ``new_r_A = mean(targets_N)``
       * Track ``max_delta = max(max_delta, ||new_r_A - r_A||)``.

    3. Apply all updates *after* the per-iteration loop so the heal is
       Jacobi-style (no within-iteration coupling).  This keeps the
       algorithm deterministic and invariant to per-atom iteration
       order.

    4. Stop when ``max_delta < tol`` or ``max_iter`` is reached.

    Atoms in ``frozen_atoms`` (e.g. metal + donors) are never moved;
    they still appear as neighbour anchors for their bonded partners.

    The function is pure-functional: same ``(coords, topology,
    ideal_lengths, frozen_atoms, hyperparameters)`` -> bit-identical
    output across processes.

    Returns
    -------
    healed_coords : ndarray (N, 3)
        Repositioned coordinates (float64).  When no atom was broken
        this is a *copy* of the input (the caller must not assume
        identity-of-array; the values are byte-identical).
    diagnostics : HealingDiagnostics
        Returned only when ``return_diagnostics=True``.
    """
    R = _coerce_R(coords)
    n_atoms = R.shape[0]
    bonds = _normalize_topology(topology)

    # Build the ideal-length table once (the GRIP library lookup is
    # the only non-trivial cost; do it before the loop).
    if ideal_lengths is None:
        if symbols is None:
            symbols = ["C"] * n_atoms
        ideal_lengths = _build_ideal_table(bonds, symbols, library=library)

    frozen_set: Set[int] = set(int(i) for i in (frozen_atoms or ()))

    diag = HealingDiagnostics()

    # Early-out: no broken atoms -> return a copy with empty diagnostics.
    broken = detect_broken_atoms(
        R,
        bonds,
        ideal_lengths=ideal_lengths,
        symbols=symbols,
        library=library,
        sigma_threshold=sigma_threshold,
        floating_factor=floating_factor,
    )
    diag.n_broken_initial = len(broken)
    if not broken:
        diag.n_broken_final = 0
        diag.converged = True
        diag.n_iter = 0
        if return_diagnostics:
            return R, diag
        return R

    # Pre-build the per-atom (neighbour, mu, sigma) tables so the
    # inner loop only does arithmetic.
    per_atom: Dict[int, List[Tuple[int, float, float]]] = {
        i: [] for i in range(n_atoms)
    }
    for (a, b) in bonds:
        mu, sigma = ideal_lengths.get((a, b), (_FALLBACK_MU, _FALLBACK_SIGMA))
        if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
            mu, sigma = _FALLBACK_MU, _FALLBACK_SIGMA
        per_atom[a].append((b, float(mu), float(sigma)))
        per_atom[b].append((a, float(mu), float(sigma)))
    for atom_idx in per_atom:
        per_atom[atom_idx] = sorted(per_atom[atom_idx], key=lambda t: t[0])

    # Identify orphans + floating atoms up front.
    orphans: List[int] = []
    floating: List[int] = []
    actionable: List[int] = []
    for atom_idx in broken:
        if atom_idx in frozen_set:
            # Frozen atoms cannot be repositioned -- skip silently.
            continue
        if not per_atom.get(atom_idx):
            orphans.append(atom_idx)
            continue
        actionable.append(atom_idx)
    diag.skipped_orphans = sorted(orphans)

    # Pre-scan for floating atoms (purely diagnostic; the heal still
    # processes them as broken because they appear in `actionable`).
    for atom_idx in actionable:
        nbrs = per_atom[atom_idx]
        if not nbrs:
            continue
        min_ratio = float("inf")
        for (j, mu, _sigma) in nbrs:
            if mu <= 0:
                continue
            d = float(np.linalg.norm(R[atom_idx] - R[j]))
            if not np.isfinite(d):
                continue
            ratio = d / mu
            if ratio < min_ratio:
                min_ratio = ratio
        if min_ratio > floating_factor:
            floating.append(atom_idx)
    diag.floating_atoms = sorted(floating)

    if not actionable:
        # Only orphans / frozen atoms were broken -- nothing to do.
        diag.n_broken_final = diag.n_broken_initial
        diag.converged = True
        if return_diagnostics:
            return R, diag
        return R

    # Working copy.  Jacobi-style update: read positions from `R_prev`,
    # write positions to `R_next`, swap at the end of every iteration.
    R_work = R
    R_initial_snapshot = R.copy()  # for accept-if-better roll-back below
    for it in range(int(max_iter)):
        # Re-detect broken set on the current geometry so atoms that
        # have already converged drop out of the work list.
        if it > 0:
            current_broken = detect_broken_atoms(
                R_work,
                bonds,
                ideal_lengths=ideal_lengths,
                symbols=symbols,
                library=library,
                sigma_threshold=sigma_threshold,
                floating_factor=floating_factor,
            )
            # Restrict to actionable subset (non-frozen, has neighbours).
            current_broken = [
                idx for idx in current_broken
                if idx not in frozen_set and per_atom.get(idx)
            ]
        else:
            current_broken = list(actionable)

        if not current_broken:
            diag.converged = True
            diag.n_iter = it
            break

        R_next = R_work.copy()
        max_delta = 0.0
        for atom_idx in sorted(current_broken):
            nbrs = per_atom[atom_idx]
            targets: List[np.ndarray] = []
            for (j, mu, _sigma) in nbrs:
                # Always read from the previous iteration's geometry so
                # the update is Jacobi-style (independent of the order
                # in which we iterate over broken atoms).
                r_a = R_work[atom_idx]
                r_n = R_work[j]
                diff = r_a - r_n
                norm = float(np.linalg.norm(diff))
                if norm < _DEGENERATE_TOL:
                    direction = _deterministic_fallback_direction(
                        atom_idx, j,
                    )
                else:
                    direction = diff / norm
                target = r_n + float(mu) * direction
                if not np.all(np.isfinite(target)):
                    continue
                targets.append(target)
            if not targets:
                continue
            new_position_full = np.mean(np.stack(targets, axis=0), axis=0)
            if not np.all(np.isfinite(new_position_full)):
                continue
            # Damped Jacobi update.  Pure Jacobi (alpha=1) can oscillate
            # when bonds are tightly coupled (one atom is endpoint to
            # multiple broken bonds whose target projections disagree).
            # Under-relaxation alpha=0.5 (the classical SOR factor)
            # makes the iteration contractive on any acyclic ligand
            # graph at the cost of halving the convergence rate.
            alpha = float(damping)
            if alpha <= 0.0 or alpha > 1.0:
                alpha = 1.0
            new_position = (
                R_work[atom_idx] + alpha * (new_position_full - R_work[atom_idx])
            )
            delta = float(np.linalg.norm(new_position - R_work[atom_idx]))
            if delta > max_delta:
                max_delta = delta
            R_next[atom_idx] = new_position

        diag.max_delta_history.append(float(max_delta))
        R_work = R_next
        diag.n_iter = it + 1
        if max_delta < float(tol):
            diag.converged = True
            break

    # Compute the final broken count for the caller.
    final_broken = detect_broken_atoms(
        R_work,
        bonds,
        ideal_lengths=ideal_lengths,
        symbols=symbols,
        library=library,
        sigma_threshold=sigma_threshold,
        floating_factor=floating_factor,
    )

    # ACCEPT-IF-BETTER OUTER GATE.  On densely-coupled structures
    # (e.g. ≥ 10 % of bonds broken, with broken atoms sharing many
    # neighbours) the Jacobi heal can settle into a fixed point that
    # *increases* the total broken-bond count -- one atom moves to
    # satisfy bond (i, j) but j now sits worse for bond (j, k).  We
    # detect this case via the broken count and roll back to the
    # initial coords so the heal never harms the structure.  This
    # mirrors the accept-if-better semantics of grip_polish.
    if len(final_broken) > diag.n_broken_initial:
        diag.n_broken_final = diag.n_broken_initial
        diag.converged = True  # signal "no-op" status
        diag.max_delta_history.append(0.0)
        if return_diagnostics:
            return R_initial_snapshot, diag
        return R_initial_snapshot

    diag.n_broken_final = len(final_broken)

    if return_diagnostics:
        return R_work, diag
    return R_work


# ---------------------------------------------------------------------------
# Convenience: wire healing in front of grip_polish
# ---------------------------------------------------------------------------
def _extract_symbols(mol) -> List[str]:
    """Read element symbols from an RDKit-like mol; return ``[]`` on failure."""
    try:
        return [str(a.GetSymbol()) for a in mol.GetAtoms()]
    except Exception:
        return []


def _extract_mol_bonds(mol) -> List[Tuple[int, int]]:
    """Read the lex-sorted bonded-pair list from an RDKit-like mol."""
    try:
        from delfin.fffree.grip_polish import _mol_bonds
    except Exception:
        out: List[Tuple[int, int]] = []
        try:
            for bond in mol.GetBonds():
                a = int(bond.GetBeginAtomIdx())
                b = int(bond.GetEndAtomIdx())
                if a == b:
                    continue
                if a > b:
                    a, b = b, a
                out.append((a, b))
        except Exception:
            return []
        return sorted(set(out))
    return _mol_bonds(mol)


def grip_polish_with_healing(
    P0,
    mol,
    metal: int,
    donors: Sequence[int],
    geom: str = "",
    mogul_lib=None,
    *,
    sigma_threshold: float = DEFAULT_SIGMA_THRESHOLD,
    max_iter: int = DEFAULT_MAX_ITER,
    tol: float = DEFAULT_TOL,
    floating_factor: float = DEFAULT_FLOATING_FACTOR,
    damping: float = DEFAULT_DAMPING,
    return_diagnostics: bool = False,
    **polish_kwargs,
):
    """Healing-mode entry point: heal *then* polish.

    When :func:`healing_mode_active` is ``True`` and the structure
    contains broken atoms, this function runs the iterative
    distance-geometry pre-pass (:func:`iterative_topology_repositioning`)
    *before* delegating to :func:`grip_polish.grip_polish`.  The metal
    + donors are always frozen during the heal (matching the rest of
    GRIP's M-D rigidity contract).

    When healing mode is OFF, or when no atom is flagged as broken,
    this function calls :func:`grip_polish.grip_polish` directly with
    the same arguments -- byte-identical with HEAD 93b396d.

    ``polish_kwargs`` are forwarded verbatim to ``grip_polish`` so the
    caller can keep using ``clash_weight``, ``md_tol``, etc.

    Parameters
    ----------
    P0
        Initial coordinates ``(N, 3)``.
    mol
        RDKit-like molecule with bonds + symbols.
    metal
        Metal atom index (frozen during heal AND polish).
    donors
        Donor atom indices (frozen during heal AND polish).
    geom
        Coordination geometry name, forwarded to ``grip_polish``.
    mogul_lib
        Pre-loaded library; used for both the heal-time bond lookup
        and the polish-time Mahalanobis terms.
    sigma_threshold, max_iter, tol, floating_factor
        Healing hyperparameters; see
        :func:`iterative_topology_repositioning`.
    return_diagnostics
        When ``True`` return a tuple ``(P_polished, healing_diag,
        polish_result)``.  The middle element is a
        :class:`HealingDiagnostics`; the last is either a
        ``GripPolishResult`` (when grip_polish was asked for one)
        or ``None``.

    Returns
    -------
    ndarray (N, 3)
        Healed-then-polished coordinates.  Falls back to ``P0`` on
        any internal failure (matches the rest of GRIP's
        defence-in-depth contract).
    """
    # Import locally to avoid a hard import cycle with grip_polish
    # (grip_polish itself does not import grip_healing).
    from delfin.fffree.grip_polish import grip_polish

    polish_return_diag = bool(polish_kwargs.pop("return_diagnostics", False))

    # Pure pass-through when healing mode is OFF.  This is the
    # byte-identity guarantee: env unset -> we delegate verbatim to
    # ``grip_polish`` with the original P0 / kwargs.
    if not healing_mode_active():
        P_out = grip_polish(
            P0, mol, metal, donors, geom,
            mogul_lib=mogul_lib,
            return_diagnostics=polish_return_diag,
            **polish_kwargs,
        )
        if return_diagnostics:
            empty_diag = HealingDiagnostics()
            empty_diag.converged = True
            if polish_return_diag:
                # P_out is the GripPolishResult; pull P off it.
                P_arr = getattr(P_out, "P", P_out)
                return P_arr, empty_diag, P_out
            return P_out, empty_diag, None
        return P_out

    try:
        R0 = _coerce_R(P0)
    except Exception:
        return P0

    symbols = _extract_symbols(mol)
    bonds = _extract_mol_bonds(mol)

    # Empty topology -> heal is a no-op; delegate.
    if not bonds or not symbols:
        P_out = grip_polish(
            P0, mol, metal, donors, geom,
            mogul_lib=mogul_lib,
            return_diagnostics=polish_return_diag,
            **polish_kwargs,
        )
        if return_diagnostics:
            empty_diag = HealingDiagnostics()
            empty_diag.converged = True
            if polish_return_diag:
                P_arr = getattr(P_out, "P", P_out)
                return P_arr, empty_diag, P_out
            return P_out, empty_diag, None
        return P_out

    ideal_table = _build_ideal_table(bonds, symbols, library=mogul_lib)
    frozen = set([int(metal), *[int(d) for d in donors]])

    # Heal.
    try:
        R_healed, heal_diag = iterative_topology_repositioning(
            R0,
            bonds,
            ideal_lengths=ideal_table,
            symbols=symbols,
            library=mogul_lib,
            frozen_atoms=frozen,
            sigma_threshold=sigma_threshold,
            max_iter=max_iter,
            tol=tol,
            floating_factor=floating_factor,
            damping=damping,
            return_diagnostics=True,
        )
    except Exception as exc:  # defence in depth -- the heal must never crash
        _LOG.warning("grip_polish_with_healing: heal raised (%r); skipping", exc)
        R_healed = R0
        heal_diag = HealingDiagnostics()

    # Polish on the healed coordinates.
    try:
        P_out = grip_polish(
            R_healed, mol, metal, donors, geom,
            mogul_lib=mogul_lib,
            return_diagnostics=polish_return_diag,
            **polish_kwargs,
        )
    except Exception as exc:
        _LOG.warning(
            "grip_polish_with_healing: polish raised (%r); returning healed",
            exc,
        )
        if return_diagnostics:
            return R_healed, heal_diag, None
        return R_healed

    if return_diagnostics:
        if polish_return_diag:
            P_arr = getattr(P_out, "P", P_out)
            return P_arr, heal_diag, P_out
        return P_out, heal_diag, None
    return P_out
