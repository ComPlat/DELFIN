"""GRIP — L-BFGS Polish + Constraint Stack (Phase 3, v1).

Operational heart of GRIP.  Turns the per-fragment Mahalanobis loss
(:mod:`grip_loss_terms`, Phase 2) and the constraint stack
(:mod:`grip_constraints`, Phase 3) into a single

.. code:: python

    P_refined = grip_polish(P0, mol, metal, donors, geom, ...)

call.  The function is deterministic, FF-free, env-agnostic and rolls back to
``P0`` (returning the original array) whenever a hard constraint would be
violated by the polished geometry or whenever the polish degrades the total
mogul severity (accept-if-better gate).

Algorithm (SPEC §3.4):

1. Detect bond / angle / improper fragments via
   :func:`grip_fragment_detect.detect_fragments` -- frozen set is
   ``{metal, *donors}``.
2. Read the M-D rigidity targets, topology bond lengths, and chiral signs
   off ``P0`` (so the constraints lock in the *initial* state, not a target
   pulled from elsewhere).
3. Build the combined objective ``L = L_GRIP + clash_weight × L_clash``.
4. Project the gradient: zero on every frozen atom -- M-D + polyhedron are
   rigid by construction; the optimiser literally cannot move them.
5. Run :func:`scipy.optimize.minimize` ``method='L-BFGS-B'`` with the SPEC
   §3.3 settings (``maxiter=200``, ``gtol=1e-4``, ``ftol=1e-7``).
6. Validate the result against M-D, topology and chirality constraints --
   any violation returns the original ``P0`` (rollback).
7. Accept-if-better gate on the mogul severity (sum of per-fragment
   ``z²`` terms): if the polished severity is ``≥`` the initial severity,
   return ``P0``.

The function is *pure*: same input -> bit-identical output across runs (no
RNG, sorted iteration order, float64 throughout, PYTHONHASHSEED-respecting).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

_LOG = logging.getLogger(__name__)

# Default Pauli-floor multiplier (Phase 3 calibration).  Overridable per call
# via the ``clash_weight`` kwarg, OR globally via the env var
# ``DELFIN_FFFREE_GRIP_CLASH_WEIGHT`` (read every call so it is honoured under
# subprocess fork pools).  Default-OFF (env unset) keeps Phase 4 byte-identity.
DEFAULT_CLASH_WEIGHT: float = 5.0
_CLASH_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_CLASH_WEIGHT"

# Optimiser-method dispatcher (2026-06-04).  Default ``lbfgs`` preserves
# byte-identity with HEAD bcf56f8 (L-BFGS-B via scipy.minimize).  Set to
# ``lm`` to dispatch to :func:`grip_polish_lm` (TRF, Trust-Region-Reflective
# Levenberg-Marquardt-like solver with bounds support).  Unknown values
# fall back to the L-BFGS default with a single warning -- the polish must
# never crash on a misconfigured env-flag.
_GRIP_METHOD_ENV: str = "DELFIN_FFFREE_GRIP_METHOD"
_GRIP_METHOD_LBFGS: str = "lbfgs"
_GRIP_METHOD_LM: str = "lm"
_GRIP_METHOD_ALIASES: Dict[str, str] = {
    "": _GRIP_METHOD_LBFGS,
    "lbfgs": _GRIP_METHOD_LBFGS,
    "l-bfgs": _GRIP_METHOD_LBFGS,
    "l-bfgs-b": _GRIP_METHOD_LBFGS,
    "default": _GRIP_METHOD_LBFGS,
    "lm": _GRIP_METHOD_LM,
    "trf": _GRIP_METHOD_LM,
    "levenberg-marquardt": _GRIP_METHOD_LM,
    "levenberg_marquardt": _GRIP_METHOD_LM,
}


def _resolve_grip_method() -> str:
    """Return the canonical method tag (``"lbfgs"`` or ``"lm"``).

    Read ``DELFIN_FFFREE_GRIP_METHOD`` (case-insensitive).  Unknown values
    fall back to ``"lbfgs"`` (default-OFF byte-identity) and emit a single
    warning so the misconfiguration is visible.
    """
    raw = os.environ.get(_GRIP_METHOD_ENV, "").strip().lower()
    tag = _GRIP_METHOD_ALIASES.get(raw)
    if tag is None:
        _LOG.warning(
            "grip_polish: unknown %s=%r; falling back to %s",
            _GRIP_METHOD_ENV, raw, _GRIP_METHOD_LBFGS,
        )
        return _GRIP_METHOD_LBFGS
    return tag

# Fix A (2026-06-02 User-Direktive): inter-ligand clash boost weight + the
# acceptance-gate extension that consumes it.  Default values keep the
# operator byte-identical with the legacy path -- only when the L-BFGS
# wrapper passes a ``ligand_atom_id`` map (or the new env-flag is set) does
# the boost kick in.
DEFAULT_INTER_LIGAND_CLASH_WEIGHT: float = 15.0
_INTER_LIGAND_CLASH_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_INTER_LIGAND_CLASH_WEIGHT"
_ACCEPT_WITH_CLASH_ENV: str = "DELFIN_FFFREE_GRIP_ACCEPT_WITH_CLASH"
_ACCEPT_WITH_CLASH_ALPHA_ENV: str = "DELFIN_FFFREE_GRIP_ACCEPT_CLASH_ALPHA"
DEFAULT_ACCEPT_CLASH_ALPHA: float = 1.0


def _resolve_inter_ligand_clash_weight(arg_value) -> float:
    """Pick the effective inter-ligand clash multiplier (fix A).

    Resolution order mirrors :func:`_resolve_clash_weight`:
    explicit arg -> env -> DEFAULT_INTER_LIGAND_CLASH_WEIGHT.
    """
    if arg_value is not None:
        try:
            v = float(arg_value)
            if np.isfinite(v):
                return v
        except (TypeError, ValueError):
            pass
    raw = os.environ.get(_INTER_LIGAND_CLASH_WEIGHT_ENV, "")
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v):
                return v
        except (TypeError, ValueError):
            pass
    return DEFAULT_INTER_LIGAND_CLASH_WEIGHT


def _accept_with_clash_active() -> bool:
    """``True`` iff the acceptance-gate clash extension is enabled.

    Reads ``DELFIN_FFFREE_GRIP_ACCEPT_WITH_CLASH``.  Default OFF (the
    severity-only gate preserves byte-identity with the legacy path); set
    to ``1``/``true`` to enable.  When the GRIP-Ensemble is active the
    caller typically pairs this flag with
    ``DELFIN_FFFREE_GRIP_INTER_LIGAND_CLASH_WEIGHT`` to also bias the
    L-BFGS gradient towards relieving inter-ligand contacts.
    """
    raw = os.environ.get(_ACCEPT_WITH_CLASH_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _accept_with_clash_alpha() -> float:
    """Coefficient for the inter-ligand clash component in the acceptance
    score (fix A).  Defaults to :data:`DEFAULT_ACCEPT_CLASH_ALPHA` (1.0).
    """
    raw = os.environ.get(_ACCEPT_WITH_CLASH_ALPHA_ENV, "").strip()
    if not raw:
        return DEFAULT_ACCEPT_CLASH_ALPHA
    try:
        v = float(raw)
        if np.isfinite(v):
            return v
    except (TypeError, ValueError):
        pass
    return DEFAULT_ACCEPT_CLASH_ALPHA


def _resolve_clash_weight(arg_value) -> float:
    """Pick the effective ClashFloorPenalty multiplier.

    Resolution order:

    1. If the caller passed an explicit numeric value (not ``None``), honour
       it.  This preserves backwards compatibility for every internal test
       that already pins ``clash_weight=5.0`` (or any other number).
    2. Otherwise consult the env var ``DELFIN_FFFREE_GRIP_CLASH_WEIGHT``.
       If it parses as a finite float, that value wins.
    3. Otherwise (unset, empty, non-numeric, NaN / inf) fall back to
       ``DEFAULT_CLASH_WEIGHT`` (5.0); a non-numeric / non-finite value is
       reported via a single warning so misconfiguration is visible without
       spamming.

    Pure-functional: same arg + same env -> same output.
    """
    if arg_value is not None:
        try:
            v = float(arg_value)
            if np.isfinite(v):
                return v
        except (TypeError, ValueError):
            pass
        # Caller passed something non-numeric / non-finite -- fall through to
        # env / default rather than crash mid-build.
        _LOG.warning(
            "grip_polish: ignoring non-numeric clash_weight=%r; using env/default",
            arg_value,
        )
    raw = os.environ.get(_CLASH_WEIGHT_ENV, "")
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v):
                return v
            _LOG.warning(
                "grip_polish: %s=%r is not finite; falling back to %.3f",
                _CLASH_WEIGHT_ENV, raw, DEFAULT_CLASH_WEIGHT,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r is not numeric; falling back to %.3f",
                _CLASH_WEIGHT_ENV, raw, DEFAULT_CLASH_WEIGHT,
            )
    return DEFAULT_CLASH_WEIGHT


from .grip_constraints import (
    ChiralVolumeConstraint,
    ClashFloorPenalty,
    DonorPolyhedronConstraint,
    MDInvariantConstraint,
    TopologyConstraint,
)
from .grip_fragment_detect import detect_fragments
from .grip_loss_terms import TotalGripLoss
from .grip_mogul_lookup import GripLibrary, get_default_library

__all__ = [
    "grip_polish",
    "mogul_severity",
    "GripPolishResult",
    "DEFAULT_VDW_RADII",
    "DEFAULT_CLASH_WEIGHT",
    "DEFAULT_INTER_LIGAND_CLASH_WEIGHT",
    "DEFAULT_ACCEPT_CLASH_ALPHA",
    "detect_hapto_atoms",
    "build_ligand_atom_id_map",
    "expand_hapto_for_sigma_only",
    "_sigma_only_mode_active",
    "_resolve_grip_method",
    "_GRIP_METHOD_ENV",
    "_GRIP_METHOD_LBFGS",
    "_GRIP_METHOD_LM",
    # GRACE env-gated dispatcher (default OFF, byte-identical to HEAD 00f1a5b).
    "grace_dispatch_active",
    "grace_polish_or_enumerate",
]


def grace_dispatch_active() -> bool:
    """``True`` iff GRACE-Ensemble dispatch is enabled.

    Reads the GRACE master env flag without importing grace_ensemble
    eagerly (so the OFF path keeps the import surface byte-identical
    to HEAD 00f1a5b).
    """
    raw = os.environ.get("DELFIN_FFFREE_GRACE_ENABLE", "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def grace_polish_or_enumerate(smiles: Optional[str] = None, **kwargs):
    """Public dispatcher: when GRACE is active, run
    :func:`delfin.fffree.grace_ensemble.grace_enumerate` and return the
    :class:`GraceResult`; otherwise return ``None`` (caller falls back
    to the single-shot :func:`grip_polish`).

    This hook is intended to be called from
    :func:`delfin.fffree.assemble_complex.assemble_from_config` (or any
    other top-level orchestrator) when a deterministic global ensemble
    is preferred over a single local polish.  The default-OFF check is
    cheap (one env-var read) so the dispatcher is safe to call from any
    hot path.
    """
    if not grace_dispatch_active():
        return None
    if not smiles:
        return None
    try:
        from .grace_ensemble import grace_enumerate
    except Exception as exc:  # pragma: no cover -- import-time wiring
        _LOG.warning("grace_polish_or_enumerate: grace_ensemble import failed: %r", exc)
        return None
    return grace_enumerate(smiles, **kwargs)


def build_ligand_atom_id_map(
    mol,
    metal_idx: int,
    donors: Sequence[int],
) -> Dict[int, int]:
    """Build the ``{atom_idx -> ligand_id}`` map used by the inter-ligand
    clash boost (fix A).

    Thin wrapper that defers to :func:`grip_ensemble.identify_ligand_subgraphs`
    so the ligand partitioning is identical between the L-BFGS loss term
    and the discrete clash filter.  Returns ``{}`` on any failure (the
    caller falls back to the legacy single-weight path).
    """
    try:
        from .grip_ensemble import identify_ligand_subgraphs
    except Exception:
        return {}
    try:
        comps = identify_ligand_subgraphs(mol, int(metal_idx), tuple(int(d) for d in donors))
    except Exception:
        return {}
    out: Dict[int, int] = {}
    for lig_id, comp in enumerate(comps):
        for atom in comp:
            out[int(atom)] = int(lig_id)
    return out


# A small, well-tested vdW table -- same numbers as :mod:`refine` plus a
# few extras for transition metals (Pauling/Bondi).  Atoms not in this
# table are skipped by the clash floor (they get no Pauli contribution),
# which is the safe behaviour (don't penalise what you can't size).
DEFAULT_VDW_RADII = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92,
    "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.11, "Ti": 1.95, "V": 1.06,
    "Cr": 1.89, "Mn": 1.97, "Fe": 1.94, "Co": 1.92, "Ni": 1.84,
    "Cu": 1.86, "Zn": 2.10, "Ga": 1.87, "Ge": 2.11, "As": 1.85,
    "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.11, "Zr": 1.86, "Nb": 2.07,
    "Mo": 2.09, "Ru": 1.97, "Rh": 1.95, "Pd": 2.02, "Ag": 2.03,
    "Cd": 2.30, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06,
    "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43,
    "Hf": 2.12, "Ta": 2.17, "W": 2.10, "Re": 2.17, "Os": 2.16,
    "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23, "Tl": 1.96,
    "Pb": 2.02, "Bi": 2.07,
}


# ---------------------------------------------------------------------------
# Helpers (no RDKit hard import; only needed when mol is real)
# ---------------------------------------------------------------------------
def _coerce_R(R: np.ndarray) -> np.ndarray:
    R = np.asarray(R, dtype=np.float64)
    if R.ndim == 1:
        if R.size % 3 != 0:
            raise ValueError(f"P0 must be (N,3) or (3N,), got length {R.size}")
        R = R.reshape(-1, 3)
    return R.copy()


def _mol_bonds(mol) -> List[Tuple[int, int]]:
    """Return sorted ``(a, b)`` pairs from a RDKit-like mol (deterministic)."""
    out: List[Tuple[int, int]] = []
    for bond in mol.GetBonds():
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        if a == b:
            continue
        if a > b:
            a, b = b, a
        out.append((a, b))
    return sorted(set(out))


def _mol_symbol_table(mol) -> Dict[int, str]:
    return {int(a.GetIdx()): str(a.GetSymbol()) for a in mol.GetAtoms()}


def _build_13_exclusions(
    mol_bonds: Sequence[Tuple[int, int]],
    n_atoms: int,
) -> Set[FrozenSet[int]]:
    """Build the set of bonded + 1,3 atom-index pairs (frozensets)."""
    adj: List[Set[int]] = [set() for _ in range(n_atoms)]
    excl: Set[FrozenSet[int]] = set()
    for (a, b) in mol_bonds:
        if 0 <= a < n_atoms and 0 <= b < n_atoms and a != b:
            adj[a].add(b)
            adj[b].add(a)
            excl.add(frozenset((a, b)))
    for b in range(n_atoms):
        nbrs = sorted(adj[b])
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                excl.add(frozenset((nbrs[i], nbrs[j])))
    return excl


def _vdw_table_for_mol(mol, vdw_radii_by_symbol: Dict[str, float]) -> Dict[int, float]:
    """Map atom indices to vdW radii using the symbol table."""
    syms = _mol_symbol_table(mol)
    out: Dict[int, float] = {}
    for idx, sym in syms.items():
        r = vdw_radii_by_symbol.get(sym)
        if r is None:
            # Unknown element -- skip; clash floor will ignore it.
            continue
        out[idx] = float(r)
    return out


def _stereocenter_quadruples(mol) -> List[Tuple[int, int, int, int]]:
    """Return ``(center, a, b, c)`` quadruples for every heavy atom with
    exactly three heavy neighbours (deterministic, sorted)."""
    quads: List[Tuple[int, int, int, int]] = []
    for atom in mol.GetAtoms():
        c = int(atom.GetIdx())
        try:
            sym = atom.GetSymbol()
        except Exception:
            continue
        if sym == "H":
            continue
        nbrs = sorted(int(n.GetIdx()) for n in atom.GetNeighbors())
        # Sterochirality lives at sp3 (or square-planar metal-free) centres
        # with at least three neighbours.  We use the FIRST THREE in sorted
        # order to define a determinant; the constraint freezes only the
        # sign of that determinant so it is invariant to which triple was
        # picked as long as the choice is deterministic and consistent.
        if len(nbrs) < 3:
            continue
        quads.append((c, nbrs[0], nbrs[1], nbrs[2]))
    return sorted(set(quads))


# ---------------------------------------------------------------------------
# Hapto-π atom-set detector (2026-06-02, voll-pool regression fix)
# ---------------------------------------------------------------------------
def detect_hapto_atoms(
    mol,
    metal_idx: int,
    donors: Sequence[int],
) -> Set[int]:
    """Return the set of atom indices that belong to a hapto-coordinated
    π-system in ``mol``.

    A donor cluster qualifies as hapto when ALL of the following hold
    (universal graph criteria, no SMILES patterns):

    * ≥ 3 of the metal's bonded donors are carbon
    * those carbons sit in the same ring (closed hapto: η³-η⁸)
      OR they form an unbroken chain in the molecular graph
      (open-chain hapto: η³-allyl, η⁴-diene)
    * the η-count (3, 4, 5, 6, 7, 8) matches the constructive
      :mod:`hapto_modes` enumeration used by ``assemble_complex.py``

    When a ring matches, the FULL ring is treated as the hapto-π system
    (so a single π-coordinated atom does not stay as a discrete σ donor
    flagged downstream); for open-chain hapto only the donor chain is
    used.

    The detector is read-only and pure: same ``(mol, metal, donors)``
    yields the same set every call (sorted iteration; no RNG).

    Parameters
    ----------
    mol : RDKit ``Mol`` / ``RWMol``
        Assembled metal-plus-ligand mol (as built by
        :func:`assemble_complex`).
    metal_idx : int
        Atom index of the metal centre.
    donors : sequence of int
        Atom indices of the constructed donors (M-D-rigid).  For a
        hapto-η⁵ Cp this contains all five ring carbons; the detector
        clusters them back into a hapto group from the molecular graph.

    Returns
    -------
    set of int
        Hapto-π atom indices.  Empty when no hapto cluster is found
        (purely σ-coordinated complex, monatomic donors, organic-only).
    """
    metal_idx = int(metal_idx)
    donor_set: Set[int] = set(int(d) for d in donors)
    if len(donor_set) < 3:
        return set()

    # The assembled ``cm`` carries metal-donor bonds (added by
    # assemble_complex.py).  Vanilla SSSR on it picks small M-C-C
    # triangles instead of the underlying Cp/arene ring, so we ring-
    # detect on a CLONE that has the metal-donor edges removed (and any
    # other metal edge to be safe).  Read-only on the caller's mol.
    rings: List[Tuple[int, ...]] = []
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
        rings = [tuple(int(a) for a in r) for r in sub.GetRingInfo().AtomRings()]
    except Exception:
        # Fall back to the caller's ring info (still works for ligands
        # whose ring is not collapsed by the M-edge).
        try:
            ri = mol.GetRingInfo()
            rings = [tuple(int(a) for a in r) for r in ri.AtomRings()]
        except Exception:
            rings = []

    # Pre-compute donor symbols + atom objects (read-only).
    donor_carbon: Set[int] = set()
    for d in donor_set:
        try:
            atom = mol.GetAtomWithIdx(int(d))
        except Exception:
            continue
        try:
            if atom.GetSymbol() == "C":
                donor_carbon.add(int(d))
        except Exception:
            continue
    if len(donor_carbon) < 3:
        return set()

    hapto: Set[int] = set()

    # --- 1) Closed-ring hapto (η³-η⁸) -----------------------------------
    # A ring qualifies when ≥ 3 of its atoms are carbon donors AND those
    # donors are all in the same ring.  The full ring atoms enter the
    # hapto set so neighbouring impropers / angles on ring carbons are
    # also protected.
    for ring in rings:
        ring_set = set(int(a) for a in ring)
        carbon_donors_in_ring = ring_set & donor_carbon
        n_eta = len(carbon_donors_in_ring)
        if n_eta in (3, 4, 5, 6, 7, 8):
            # Confirm the ring carbons are all C (otherwise this is a
            # mixed heteroaromatic ring with several donors — not the
            # ferrocene-Cp pattern we want to protect).
            try:
                ring_all_c = all(
                    mol.GetAtomWithIdx(int(a)).GetSymbol() == "C"
                    for a in ring
                )
            except Exception:
                ring_all_c = False
            if not ring_all_c:
                continue
            hapto.update(ring_set)

    # --- 2) Open-chain hapto (η³-allyl, η⁴-diene) ----------------------
    # If we did not already cover the carbon donors via a ring, look for
    # a chain of length 3 or 4 in the molecular graph.  Build a tiny
    # adjacency graph restricted to donor carbons + their bonds, then
    # check if it is a path (2 endpoints with deg 1, middle atoms deg 2).
    uncovered = donor_carbon - hapto
    if 3 <= len(uncovered) <= 4:
        adj = {a: set() for a in uncovered}
        try:
            for bond in mol.GetBonds():
                u = int(bond.GetBeginAtomIdx())
                v = int(bond.GetEndAtomIdx())
                if u in adj and v in adj:
                    adj[u].add(v)
                    adj[v].add(u)
        except Exception:
            adj = {}
        if adj:
            degs = sorted(len(adj[a]) for a in adj)
            n = len(adj)
            if n in (3, 4) and degs == [1, 1] + [2] * (n - 2):
                # Walk the chain — connectedness (single component).
                start = next(a for a in adj if len(adj[a]) == 1)
                seen = {start}
                stack = [start]
                while stack:
                    cur = stack.pop()
                    for nb in adj[cur]:
                        if nb not in seen:
                            seen.add(nb)
                            stack.append(nb)
                if seen == set(adj.keys()):
                    hapto.update(seen)

    # The metal itself is M-D-rigid (frozen) so it does not need to be in
    # the hapto set — but the hapto-π atoms are typically a SUPERSET of
    # the carbon donors, so we leave the donors in the set and rely on
    # the detect_fragments frozen/donor union to keep their priors
    # excluded.  Including them in hapto_atoms costs nothing.
    return hapto


# ---------------------------------------------------------------------------
# Class-conditional σ-only mode (2026-06-03):
# DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE expands the hapto-protected atom set so
# the L-BFGS polish only acts on σ-bonded ligand internals.  The π-system
# (every ring carrying a C-donor and every H attached to such a carbon) is
# excluded from the loss entirely, which:
#   (a) addresses the hapto_geom_per_frame regression observed in the
#       race-stack smoke (+2862 % vs f8c9905) — the Burnside-Konformer /
#       Symmetry-Priority touches piano-stool / sandwich geometries the
#       hapto-protection (Heal Option 1) handled at the η3-η8 ring level
#       but not at the η1/η2 / aromatic-σ-N ring level;
#   (b) supports the "fffree+GRIP achieves σ-only-cshm < UFF on the σ
#       sub-polyhedron" publication claim — the polish only optimises σ
#       internals; the π contribution is left to the placement step.
# Default OFF -> byte-identical to HEAD c03a550.
# ---------------------------------------------------------------------------
_SIGMA_ONLY_MODE_ENV: str = "DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE"


def _sigma_only_mode_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE`` is on (default OFF)."""
    raw = os.environ.get(_SIGMA_ONLY_MODE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def expand_hapto_for_sigma_only(
    mol,
    metal_idx: int,
    donors: Sequence[int],
    base_hapto: Set[int],
    md_cut: float = 2.6,
) -> Set[int]:
    """Extend ``base_hapto`` to cover the WHOLE π-system of every
    C-donor in the M-D-cut shell, plus every H attached to a π atom.

    The σ-only mode protects more than the η³-η⁸ rings detect_hapto_atoms
    catches: any C-donor in an aromatic ring (even η¹ phenyl, η² alkene)
    is treated as π and excluded together with the ring carbons + the H
    atoms riding on the ring (the H–C bond/angle priors otherwise drag
    the ring into a sp² geometry that fights the M-η placement).

    Pure-functional, deterministic; same ``(mol, metal_idx, donors,
    base_hapto)`` -> same returned set.  ``base_hapto`` is NOT mutated.

    Parameters
    ----------
    mol : RDKit Mol-like
        Read-only molecular graph.
    metal_idx : int
        Atom index of the metal.
    donors : sequence of int
        Atom indices of the M-D-rigid donors (the constructed shell).
    base_hapto : set of int
        Starting set from :func:`detect_hapto_atoms`.  May be empty.
    md_cut : float
        M-C distance cutoff (Å) above which a C-donor is not eligible for
        π-protection in σ-only mode.  Default 2.6 Å — matches the SPEC
        and the standard η M-C distance (typical η⁵-Cp M-C ≈ 2.0-2.2 Å).

    Returns
    -------
    set of int
        Expanded π-set.  Always a SUPERSET of ``base_hapto``.
    """
    out: Set[int] = set(int(i) for i in base_hapto)
    if mol is None:
        return out
    try:
        donor_set: Set[int] = set(int(d) for d in donors)
    except Exception:
        return out
    metal_idx = int(metal_idx)

    # --- 1) Identify C-donors that are within md_cut of the metal AND in a
    #        ring (aromatic OR any ring).  Use mol's built-in geometry/ring
    #        info; the M-D rigidity has already been enforced by the
    #        builder, so the M-C distance comes from the bond table when no
    #        coordinates are available.
    try:
        from rdkit import Chem
        from rdkit.Chem import RWMol
    except Exception:
        # No RDKit available -> we can only return the base set.
        return out

    # Build a ring view that skips the metal-donor edges so the underlying
    # Cp/arene ring is not collapsed into M-C-C triangles.
    try:
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
        rings = [tuple(int(a) for a in r) for r in sub.GetRingInfo().AtomRings()]
    except Exception:
        try:
            ri = mol.GetRingInfo()
            rings = [tuple(int(a) for a in r) for r in ri.AtomRings()]
        except Exception:
            rings = []

    # --- 2) For each donor C, find the rings it sits in.  When the ring
    #        is all-C OR any-aromatic -> include the WHOLE ring.  This
    #        catches η¹ phenyl (1 C donor in a phenyl ring) too, which
    #        detect_hapto_atoms (which requires >= 3 C donors per ring)
    #        does not.
    donor_C: Set[int] = set()
    for d in donor_set:
        try:
            atom = mol.GetAtomWithIdx(int(d))
        except Exception:
            continue
        try:
            if atom.GetSymbol() == "C":
                donor_C.add(int(d))
        except Exception:
            continue

    for ring in rings:
        ring_set = set(int(a) for a in ring)
        if not (ring_set & donor_C):
            continue
        # Aromatic ring OR all-C ring -> include in π-protection.
        try:
            ring_aromatic = all(
                mol.GetAtomWithIdx(int(a)).GetIsAromatic() for a in ring
            )
        except Exception:
            ring_aromatic = False
        try:
            ring_all_c = all(
                mol.GetAtomWithIdx(int(a)).GetSymbol() == "C" for a in ring
            )
        except Exception:
            ring_all_c = False
        if ring_aromatic or ring_all_c:
            out.update(ring_set)

    # --- 3) Include H atoms attached to any π atom.  H-C ring-bond
    #        priors otherwise hold the ring geometry to its free-ligand
    #        sp² template, fighting the M-η placement.
    pi_atoms_snapshot = set(out)
    for a in pi_atoms_snapshot:
        try:
            atom = mol.GetAtomWithIdx(int(a))
        except Exception:
            continue
        try:
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() == "H":
                    out.add(int(nb.GetIdx()))
        except Exception:
            continue

    return out


# ---------------------------------------------------------------------------
# Mogul severity (the iter_gate-facing metric)
# ---------------------------------------------------------------------------
def mogul_severity(
    R: np.ndarray,
    fragments: TotalGripLoss,
    mogul_lib: Optional[GripLibrary] = None,  # noqa: ARG001 -- API parity
) -> float:
    """Total weighted Mahalanobis severity of ``R`` under ``fragments``.

    This is the same quantity the loss minimises -- pulled out as a standalone
    function so the accept-if-better gate compares apples-to-apples.

    ``mogul_lib`` is accepted for API symmetry with :func:`grip_polish` but
    not used here: the ``fragments`` aggregator already carries its
    ``(mu, sigma)`` per term.
    """
    if fragments is None or len(fragments) == 0:
        return 0.0
    R = np.asarray(R, dtype=np.float64)
    if R.ndim == 1:
        R = R.reshape(-1, 3)
    sev, _ = fragments.evaluate(R)
    return float(sev)


# ---------------------------------------------------------------------------
# Result container (debug-friendly; main API still returns ndarray)
# ---------------------------------------------------------------------------
@dataclass
class GripPolishResult:
    """Extended return value for callers that want the diagnostic detail."""

    P: Optional[np.ndarray]
    accepted: bool
    severity_before: float
    severity_after: float
    n_iter: int
    n_terms: int
    rollback_reason: str = ""


# ---------------------------------------------------------------------------
# Main polish
# ---------------------------------------------------------------------------
def grip_polish(
    P0: np.ndarray,
    mol,
    metal: int,
    donors: Sequence[int],
    geom: str = "",
    mogul_lib: Optional[GripLibrary] = None,
    *,
    max_iter: int = 200,
    gtol: float = 1e-4,
    ftol: float = 1e-7,
    clash_weight: Optional[float] = None,
    md_tol: float = 0.05,
    topo_max_multiplier: float = 1.5,
    cshm_tol: float = 5.0,
    vdw_radii_by_symbol: Optional[Dict[str, float]] = None,
    return_diagnostics: bool = False,
):
    """Polish ``P0`` by minimising the GRIP loss under hard constraints.

    Parameters
    ----------
    P0 : ndarray (N, 3) or (3N,)
        Initial atomic coordinates from fffree's constructive builder.
    mol : RDKit Mol-like
        Molecule with bonds, hybridisation, neighbours populated.
    metal : int
        Atom index of the metal centre.
    donors : sequence of int
        Atom indices of the donor atoms (M-D rigid).
    geom : str
        Coordination geometry name ('OC-6', 'TBP-5', ...).  Used by the
        donor-polyhedron safety net only -- empty string disables it.
    mogul_lib : GripLibrary, optional
        Pre-loaded library; defaults to the singleton at
        :data:`grip_mogul_lookup.DEFAULT_LIB_PATH`.
    max_iter, gtol, ftol : float
        L-BFGS-B hyperparameters (SPEC §3.3).
    clash_weight : float, optional
        Multiplier on the Pauli-floor penalty.  ``None`` (default) =
        consult ``$DELFIN_FFFREE_GRIP_CLASH_WEIGHT``; if that is unset or
        not a finite float, fall back to :data:`DEFAULT_CLASH_WEIGHT`
        (5.0).  An explicit numeric value always wins.  5.0 is the
        Phase-3 calibration where the soft repulsion is strong enough to
        keep non-bonded pairs apart but not so strong that it dominates
        the Mahalanobis bond/angle pull.
    md_tol : float
        Half-width tolerance for the M-D-invariant validator (Å).
    topo_max_multiplier : float
        Bond stretch multiplier above which topology is considered broken.
    cshm_tol : float
        CShM tolerance for the donor-polyhedron safety net.
    vdw_radii_by_symbol : dict, optional
        Override the default vdW radii table.  Atoms missing from the
        table get no clash contribution.
    return_diagnostics : bool
        If True returns a :class:`GripPolishResult`; otherwise returns the
        polished ndarray (or ``P0`` on failure -- never ``None`` to keep
        the function safe to drop in).

    Returns
    -------
    ndarray (N, 3) : the polished coordinates, or the original ``P0`` if the
        polish failed any constraint or did not improve severity.

    When called with ``return_diagnostics=True``, returns
    :class:`GripPolishResult`.
    """
    # ------------------------------------------------------------------
    # 0. Coerce inputs.  Always operate on a local float64 copy so the
    # caller's array is never mutated.
    # ------------------------------------------------------------------
    P_init = _coerce_R(P0)
    n_atoms = P_init.shape[0]
    if not np.all(np.isfinite(P_init)):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False, severity_before=float("nan"),
                severity_after=float("nan"), n_iter=0, n_terms=0,
                rollback_reason="P0 contains non-finite values",
            )
        return P_init

    metal = int(metal)
    donors_t = tuple(int(d) for d in donors)
    vdw_table = vdw_radii_by_symbol if vdw_radii_by_symbol is not None else DEFAULT_VDW_RADII

    # ------------------------------------------------------------------
    # 1. Detect fragments.  Frozen = {metal, *donors}.
    # ------------------------------------------------------------------
    library = mogul_lib if mogul_lib is not None else None
    if library is None:
        try:
            library = get_default_library()
        except Exception:
            library = None
    frozen: FrozenSet[int] = frozenset((metal, *donors_t))
    # Option-B (2026-06-01): adaptive shell-1 protection.  Env-gated so
    # the polish operator can be A/B-tested in equal-n smokes without a
    # code change.  Default ON because Option-B is the current best-known
    # heuristic for preserving CCDC pressure on funcgrp internals; pass
    # ``DELFIN_FFFREE_GRIP_ADAPTIVE_SHELL1=0`` to fall back to pure
    # Heal-1/1b for forensic comparison.
    _adaptive_env = os.environ.get(
        "DELFIN_FFFREE_GRIP_ADAPTIVE_SHELL1", "1"
    ).strip()
    _adaptive_shell1 = _adaptive_env not in ("0", "false", "False", "no", "")
    # Hapto-class protection (2026-06-02): detect the hapto-π atom set
    # from the molecular graph and pass it to the fragment detector.
    # ``DELFIN_FFFREE_GRIP_SKIP_HAPTO`` (default "1") toggles the skip
    # inside detect_fragments — we always compute the set so the
    # forensic A/B remains a pure env-flag operation (no code path
    # divergence).  When the molecule has no hapto cluster the detector
    # returns an empty set and the call is byte-identical to the pre-fix
    # path.
    try:
        hapto_atoms = detect_hapto_atoms(mol, metal, donors_t)
    except Exception:
        hapto_atoms = set()
    # Class-conditional σ-only mode (2026-06-03): when active, expand the
    # hapto-protected atom set to cover the WHOLE π-system (every aromatic /
    # all-C ring containing a C-donor + the H atoms on those rings).  This
    # is a SUPERSET of the η³-η⁸ detection, so byte-identity with the env
    # OFF code-path is preserved (default-OFF returns the original set).
    if _sigma_only_mode_active():
        try:
            hapto_atoms = expand_hapto_for_sigma_only(
                mol, metal, donors_t, hapto_atoms,
            )
        except Exception:
            # On any failure fall back to the base set; the polish then
            # behaves exactly as if σ-only-mode was off, which is the
            # safe default.
            pass
    try:
        # Heal-1 (2026-06-01, mddir-fix): pass donors explicitly so the
        # detector can additionally protect the donor's first-shell
        # neighbours from being central-atom in angle / improper terms.
        # Legacy frozen_atoms still also passed for byte-compat with
        # other callers and for safety.
        # Option-B (2026-06-01): the library is also used inside
        # detect_fragments to query per-class coverage; passing it
        # explicitly keeps the assignment deterministic across forks.
        # Hapto (2026-06-02): pass the graph-derived hapto-π atom set
        # so all terms touching it are skipped — env-flag controlled.
        terms = detect_fragments(
            mol, P_init, frozen_atoms=frozen, donors=donors_t,
            hapto_atoms=hapto_atoms,
            library=library, return_result=False,
            adaptive_shell1=_adaptive_shell1,
        )
    except Exception:
        terms = []

    fragments = TotalGripLoss(terms=list(terms))

    # ------------------------------------------------------------------
    # 2. Build the constraint stack from the INITIAL state.
    # ------------------------------------------------------------------
    md_constraint = MDInvariantConstraint(
        metal_idx=metal,
        donor_idxs=donors_t,
        target_distances=tuple(
            float(np.linalg.norm(P_init[d] - P_init[metal])) for d in donors_t
        ),
        tol=md_tol,
    )

    mol_bonds = _mol_bonds(mol)
    topo_constraint = TopologyConstraint.from_initial(
        mol_bonds, P_init, max_distance_multiplier=topo_max_multiplier,
    )

    chiral_constraint = ChiralVolumeConstraint.from_initial(
        _stereocenter_quadruples(mol), P_init,
    )

    # Donor polyhedron: only if geom + ≥ 2 donors known.  Falls back to a
    # benign no-op (CShM 0.0) when geom == "".
    poly_constraint: Optional[DonorPolyhedronConstraint] = None
    if geom and len(donors_t) >= 2:
        try:
            from .polyhedra import ref_vectors
            poly_constraint = DonorPolyhedronConstraint(
                metal_idx=metal,
                donor_idxs=donors_t,
                ref_vectors_normalized=ref_vectors(geom),
                target_md_distances=md_constraint.target_distances,
                geometry=geom,
                cshm_tol=cshm_tol,
            )
        except Exception:
            poly_constraint = None

    # ------------------------------------------------------------------
    # 3. Clash floor (Pauli).
    # ------------------------------------------------------------------
    vdw_idx_table = _vdw_table_for_mol(mol, vdw_table)
    excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
    effective_clash_weight = _resolve_clash_weight(clash_weight)
    # Fix A: inter-ligand clash boost (User-Direktive 2026-06-02).  Only
    # activated when the env-flag is set to a finite numeric value (the
    # default behaviour resolves w_inter to the intra-weight, leaving the
    # gradient byte-identical with the legacy path).  Build the ligand
    # partition lazily so we don't pay the cost when the boost is OFF.
    _inter_env_raw = os.environ.get(_INTER_LIGAND_CLASH_WEIGHT_ENV, "").strip()
    ligand_map: Optional[Dict[int, int]] = None
    if _inter_env_raw:
        try:
            ligand_map = build_ligand_atom_id_map(mol, metal, donors_t)
            if not ligand_map:
                ligand_map = None
        except Exception:
            ligand_map = None
    w_inter_resolved = (
        _resolve_inter_ligand_clash_weight(None)
        if (ligand_map is not None) else effective_clash_weight
    )
    clash = ClashFloorPenalty(
        vdw_radii=vdw_idx_table,
        exclude_13_pairs=excl_13,
        floor_fraction=0.85,
        weight=float(effective_clash_weight),
        ligand_atom_id=ligand_map,
        w_inter=float(w_inter_resolved),
    )

    sev_before = mogul_severity(P_init, fragments)

    # ------------------------------------------------------------------
    # 4. Combined loss + frozen-gradient projection.
    # ------------------------------------------------------------------
    frozen_indices = np.array(sorted(frozen), dtype=np.int64)

    def loss_and_grad(x_flat: np.ndarray) -> Tuple[float, np.ndarray]:
        R = np.asarray(x_flat, dtype=np.float64).reshape(n_atoms, 3)
        L_grip, G_grip = fragments.evaluate(R) if len(fragments) > 0 else (0.0, np.zeros_like(R))
        L_clash, G_clash = clash.value_and_grad(R)
        L = float(L_grip) + float(L_clash)
        G = G_grip + G_clash
        # Zero the gradient on the frozen sphere -- L-BFGS-B will then leave
        # those atoms in place (the metal + donors stay locked).
        if frozen_indices.size:
            G[frozen_indices] = 0.0
        # Guard against NaNs creeping in (e.g. coincident atoms in a step).
        if not np.isfinite(L):
            L = 0.0
        np.nan_to_num(G, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        return L, G.reshape(-1)

    # ------------------------------------------------------------------
    # 5. L-BFGS-B.
    # ------------------------------------------------------------------
    # Special case: no terms AND no clash pairs -> nothing to do.  Return
    # the original.  Also: if the optimiser would receive a constant-zero
    # objective it can churn for max_iter steps and return spurious noise.
    L0, _ = loss_and_grad(P_init.reshape(-1))
    if len(fragments) == 0 and L0 == 0.0:
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=0.0, severity_after=0.0,
                n_iter=0, n_terms=0,
                rollback_reason="no fragments + no clashes",
            )
        return P_init

    # ------------------------------------------------------------------
    # 5a. Method dispatcher (2026-06-04).
    #
    # The default ``lbfgs`` path runs scipy.minimize(L-BFGS-B) exactly as
    # before -- byte-identical with HEAD bcf56f8 when the env-flag is unset
    # (resolver returns "lbfgs").
    #
    # The ``lm`` path runs scipy.optimize.least_squares(method='trf') via
    # :func:`grip_lm.grip_polish_lm` on the SAME prepared inputs (frozen
    # set, fragments, clash data, ligand-aware weights).  When LM diverges
    # or scipy is missing we fall through to L-BFGS so the polish never
    # crashes on a misconfigured env-flag (defence in depth).
    # ------------------------------------------------------------------
    grip_method = _resolve_grip_method()

    if grip_method == _GRIP_METHOD_LM:
        try:
            from .grip_lm import grip_polish_lm
        except Exception as exc:
            _LOG.warning(
                "grip_polish: failed to import grip_lm (%r); falling back to L-BFGS",
                exc,
            )
            grip_method = _GRIP_METHOD_LBFGS

    if grip_method == _GRIP_METHOD_LM:
        # Build the pre-screened clash pair list (deterministic, sorted by
        # (i, j)) and the radii ndarray expected by the LM residual path.
        from .grip_lm import (
            _enumerate_clash_pairs as _lm_enumerate_clash_pairs,
            grip_polish_lm,
        )
        radii_arr = np.full(n_atoms, np.nan, dtype=np.float64)
        for idx, r in vdw_idx_table.items():
            if 0 <= idx < n_atoms:
                radii_arr[idx] = float(r)
        try:
            lm_clash_pairs = _lm_enumerate_clash_pairs(
                radii_arr, excl_13, n_atoms,
            )
        except Exception as exc:
            _LOG.warning(
                "grip_polish: clash-pair enumeration failed (%r); LM falling back to L-BFGS",
                exc,
            )
            lm_clash_pairs = None

        # Optional inter-ligand weight override: when the L-BFGS path
        # would boost inter-ligand pairs (ligand_map is not None), expose
        # the same boosted weight to the LM residual via a per-pair map.
        inter_lig_weight_map: Optional[Dict[Tuple[int, int], float]] = None
        if ligand_map is not None and lm_clash_pairs is not None:
            inter_lig_weight_map = {}
            for (i, j) in lm_clash_pairs:
                li = ligand_map.get(int(i))
                lj = ligand_map.get(int(j))
                if li is not None and lj is not None and li != lj:
                    inter_lig_weight_map[(int(i), int(j))] = float(w_inter_resolved)
                else:
                    inter_lig_weight_map[(int(i), int(j))] = float(effective_clash_weight)

        if lm_clash_pairs is not None:
            try:
                P_refined, lm_diag = grip_polish_lm(
                    P_init,
                    fragments=fragments,
                    clash_pairs=lm_clash_pairs,
                    radii=radii_arr,
                    clash_weight=float(effective_clash_weight),
                    floor_fraction=0.85,
                    frozen_atoms=frozen,
                    n_atoms=n_atoms,
                    max_nfev=int(max_iter),
                    ftol=float(ftol),
                    xtol=float(gtol),
                    gtol=float(gtol),
                    inter_ligand_weight_map=inter_lig_weight_map,
                )
                n_iter = int(lm_diag.get("n_nfev", 0))
                if not lm_diag.get("ok", False):
                    _LOG.info(
                        "grip_polish[lm]: solver bailed (%s); rolling back to P0",
                        lm_diag.get("reason", "unknown"),
                    )
                    if return_diagnostics:
                        return GripPolishResult(
                            P=P_init, accepted=False,
                            severity_before=sev_before, severity_after=sev_before,
                            n_iter=n_iter, n_terms=len(fragments),
                            rollback_reason=f"lm: {lm_diag.get('reason', 'unknown')}",
                        )
                    return P_init
            except Exception as exc:
                _LOG.warning(
                    "grip_polish: LM path raised (%r); falling back to L-BFGS",
                    exc,
                )
                grip_method = _GRIP_METHOD_LBFGS

    if grip_method == _GRIP_METHOD_LBFGS:
        try:
            from scipy.optimize import minimize
        except Exception as exc:  # pragma: no cover -- scipy is in the env
            if return_diagnostics:
                return GripPolishResult(
                    P=P_init, accepted=False,
                    severity_before=sev_before, severity_after=sev_before,
                    n_iter=0, n_terms=len(fragments),
                    rollback_reason=f"scipy unavailable: {exc!r}",
                )
            return P_init

        res = minimize(
            loss_and_grad,
            P_init.reshape(-1).copy(),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": int(max_iter),
                "gtol": float(gtol),
                "ftol": float(ftol),
            },
        )
        P_refined = np.asarray(res.x, dtype=np.float64).reshape(n_atoms, 3)
        n_iter = int(getattr(res, "nit", 0))

    # ------------------------------------------------------------------
    # 6. Hard-constraint validation.  Any failure -> rollback to P0.
    # ------------------------------------------------------------------
    if not np.all(np.isfinite(P_refined)):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before, severity_after=float("nan"),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="non-finite polish result",
            )
        return P_init

    if not md_constraint.validate(P_refined):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before,
                severity_after=mogul_severity(P_refined, fragments),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="M-D invariant violated",
            )
        return P_init

    if not topo_constraint.validate(P_refined):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before,
                severity_after=mogul_severity(P_refined, fragments),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="topology bond stretched past multiplier",
            )
        return P_init

    if not chiral_constraint.validate(P_refined):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before,
                severity_after=mogul_severity(P_refined, fragments),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="chiral volume sign flipped",
            )
        return P_init

    if poly_constraint is not None and not poly_constraint.validate(P_refined):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before,
                severity_after=mogul_severity(P_refined, fragments),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="donor polyhedron CShM drifted past tolerance",
            )
        return P_init

    # ------------------------------------------------------------------
    # 7. Accept-if-better gate.  Default = pure mogul severity (legacy).
    #    Fix A (2026-06-02): when DELFIN_FFFREE_GRIP_ACCEPT_WITH_CLASH=1
    #    is set, extend the score with the inter-ligand clash count so a
    #    step that relieves an inter-ligand clash but marginally raises
    #    severity is still accepted -- this matches the loss the L-BFGS
    #    is now minimising (severity + boosted clash floor).
    # ------------------------------------------------------------------
    sev_after = mogul_severity(P_refined, fragments)
    use_clash_gate = _accept_with_clash_active()
    score_before = float(sev_before)
    score_after = float(sev_after)
    if use_clash_gate:
        try:
            from .grip_ensemble import count_inter_ligand_clashes
            n_before = int(count_inter_ligand_clashes(
                P_init, mol, metal_idx=metal, donors=donors_t,
            ))
            n_after = int(count_inter_ligand_clashes(
                P_refined, mol, metal_idx=metal, donors=donors_t,
            ))
            alpha = _accept_with_clash_alpha()
            score_before = float(sev_before) + alpha * float(n_before)
            score_after = float(sev_after) + alpha * float(n_after)
        except Exception:
            # Clash-augmented gate failed -> fall back to severity-only.
            score_before = float(sev_before)
            score_after = float(sev_after)
    if not np.isfinite(score_after) or score_after >= score_before:
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before, severity_after=sev_after,
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="accept-if-better gate (no improvement)",
            )
        return P_init

    if return_diagnostics:
        return GripPolishResult(
            P=P_refined, accepted=True,
            severity_before=sev_before, severity_after=sev_after,
            n_iter=n_iter, n_terms=len(fragments),
            rollback_reason="",
        )
    return P_refined
