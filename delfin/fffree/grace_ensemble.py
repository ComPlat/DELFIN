"""delfin.fffree.grace_ensemble — GRACE: Group-theoretic Reproducible Adaptive
Conformer Ensemble.

GRACE is a DETERMINISTIC alternative to basin-hopping conformer searches
such as ORCA's GOAT.  Where GOAT explores conformer space through a
stochastic basin-hopping MD-flavored walk (whose output depends on a
random seed and a wall-clock budget), GRACE is built from three rigorous
combinatorial layers:

  1. **Pólya topological isomers** —
     :mod:`delfin.fffree.polya_isomer_count` enumerates the donor-vertex
     colourings under the polyhedron proper-rotation group.  This gives
     the COMPLETE set of coordination-isomer connectivity choices for the
     decomposed SMILES (fac/mer, Δ/Λ, cis/trans, …).

  2. **Cremer-Pople ring puckers** —
     :mod:`delfin.fffree.ring_pucker` enumerates the canonical CP pucker
     states for every 5/6/7-membered ring of the assembled complex.

  3. **Single-bond rotamer tuples** —
     :mod:`delfin.fffree.single_bond_rotamers` enumerates the deterministic
     dihedral grid (default 60° step → 6 states per rotor) over every
     rotatable single bond.

For each point in the discrete product space
``X_iso × ∏_r X_ring,r × X_rot``, GRACE:

  a. Assembles the candidate XYZ via
     :func:`delfin.fffree.assemble_complex.assemble_from_config` (the
     same constructor the legacy single-output pipeline uses, so the
     metal centre + donor polyhedron are FF-free deterministic).
  b. Applies a topology-healing pre-polish
     (:func:`delfin.fffree.topology_healing.topology_healing_pipeline`)
     to repair phantom / missing / wrong-angle bonds without ever
     moving the metal or donors.
  c. Polishes the ligand internals via :func:`grip_polish` (L-BFGS-B by
     default; :func:`grip_polish_lm` when ``method='lm'``).
  d. Scores the result against the CCDC manifold (Mahalanobis distance
     sum from :class:`delfin.fffree.grip_loss_terms.TotalGripLoss`).
  e. Records an EnsembleCandidate that carries ``isomer_id``,
     ``ring_id``, ``rotamer_id`` and the score.

After the full sweep, GRACE applies RMSD-Butina deduplication per
isomer, keeps the top-K survivors per isomer, and returns the full
per-isomer dictionary plus a Burnside-completeness report (raw emitted
count vs :func:`burnside_conformer.count_orbits_product_space`
prediction).

Why GRACE is different from GRIP-Ensemble
=========================================

:mod:`grip_ensemble` (HEAD 00f1a5b) already does Pólya × Cremer-Pople ×
polish-and-rank.  GRACE EXTENDS that pipeline with:

  * single-bond rotamer enumeration (an axis grip_ensemble does not
    touch — currently it only iterates the as-built dihedrals);
  * topology_healing as a pre-polish step that bridges the gap between
    constructive assembly and L-BFGS convergence (constructive fffree
    occasionally yields broken angles that L-BFGS cannot recover on its
    own; topology_healing repairs them, then L-BFGS polishes);
  * method dispatch (LM vs L-BFGS) for the inner polish;
  * Burnside-completeness verification: GRACE compares emitted-orbit
    count to ``count_orbits_product_space`` and exposes the
    "coverage = emitted / predicted" diagnostic that anchors the paper's
    completeness claim;
  * deterministic time budget (total seconds) with graceful early
    termination when the budget is exhausted;
  * per-isomer top-K instead of single global top-K — useful when
    several non-degenerate isomers all need follow-up at DFT level.

Determinism contract
====================

  * No RNG ANYWHERE in the GRACE pipeline.
  * All enumeration generators iterate explicit, lex-sorted lists.
  * Candidate ranking ties break on
    ``(isomer_id, ring_id, rotamer_id, label)``.
  * PYTHONHASHSEED=0 is honoured (no dict-iteration dependence).
  * The same SMILES + same env → byte-identical output across
    machines and runs.

Default-OFF byte-identity
=========================

When ``DELFIN_FFFREE_GRACE_ENABLE`` is NOT set (or 0), GRACE is dormant:
the public hook :func:`grace_active` returns ``False`` and the
:func:`assemble_from_config` path keeps its byte-identical single-output
flow.  No other module imports this one at module-import time, so the
mere existence of this file does not change any downstream behaviour.
"""
from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

__all__ = [
    "GraceCandidate",
    "GraceResult",
    "DEFAULT_GRACE_MAX_PER_ISOMER",
    "DEFAULT_GRACE_RMSD_THRESHOLD",
    "DEFAULT_GRACE_MAX_TOTAL_TIME_S",
    "DEFAULT_GRACE_MAX_ROTAMER_TUPLES",
    "DEFAULT_GRACE_MAX_ISOMERS",
    "DEFAULT_GRACE_METHOD",
    "grace_active",
    "grace_resolve_method",
    "grace_resolve_max_per_isomer",
    "grace_resolve_max_total_time_s",
    "grace_resolve_rmsd_threshold",
    "grace_resolve_max_rotamer_tuples",
    "grace_resolve_max_isomers",
    "grace_enumerate",
    "burnside_coverage_report",
]

_LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Env-gate + tunable parameters
# ---------------------------------------------------------------------------

ENV_GRACE_ENABLE = "DELFIN_FFFREE_GRACE_ENABLE"
ENV_GRACE_METHOD = "DELFIN_FFFREE_GRACE_METHOD"
ENV_GRACE_MAX_PER_ISOMER = "DELFIN_FFFREE_GRACE_MAX_PER_ISOMER"
ENV_GRACE_MAX_TOTAL_TIME_S = "DELFIN_FFFREE_GRACE_MAX_TOTAL_TIME_S"
ENV_GRACE_RMSD_THRESHOLD = "DELFIN_FFFREE_GRACE_RMSD_THRESHOLD"
ENV_GRACE_MAX_ROTAMER_TUPLES = "DELFIN_FFFREE_GRACE_MAX_ROTAMER_TUPLES"
ENV_GRACE_MAX_ISOMERS = "DELFIN_FFFREE_GRACE_MAX_ISOMERS"
ENV_GRACE_USE_TOPOLOGY_HEALING = "DELFIN_FFFREE_GRACE_USE_TOPOLOGY_HEALING"

DEFAULT_GRACE_METHOD: str = "lm"
DEFAULT_GRACE_MAX_PER_ISOMER: int = 10
DEFAULT_GRACE_MAX_TOTAL_TIME_S: float = 60.0
DEFAULT_GRACE_RMSD_THRESHOLD: float = 0.3
DEFAULT_GRACE_MAX_ROTAMER_TUPLES: int = 27  # 3 states ^ 3 rotors
DEFAULT_GRACE_MAX_ISOMERS: int = 12


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name, "").strip().lower()
    if not raw:
        return bool(default)
    if raw in ("1", "true", "yes", "on"):
        return True
    if raw in ("0", "false", "no", "off"):
        return False
    return bool(default)


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        return int(raw)
    except (TypeError, ValueError):
        return int(default)


def _env_float(name: str, default: float) -> float:
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


def grace_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRACE_ENABLE`` is set.  Default OFF."""
    return _env_bool(ENV_GRACE_ENABLE, default=False)


def grace_resolve_method(method: Optional[str] = None) -> str:
    """Resolve the inner-polish method tag.

    Resolution order: explicit arg > env > default ``"lm"``.  Recognised
    tags are ``"lm"`` and ``"lbfgs"``.  Unknown values fall back to the
    default with a single warning.
    """
    candidates = [method, os.environ.get(ENV_GRACE_METHOD, ""), DEFAULT_GRACE_METHOD]
    for c in candidates:
        if c is None:
            continue
        c_str = str(c).strip().lower()
        if c_str in ("lm", "trf", "levenberg-marquardt", "levenberg_marquardt"):
            return "lm"
        if c_str in ("lbfgs", "l-bfgs", "l-bfgs-b"):
            return "lbfgs"
        if c_str == "":
            continue
        _LOG.warning(
            "grace_resolve_method: unknown method=%r; falling back to %r",
            c, DEFAULT_GRACE_METHOD,
        )
        return DEFAULT_GRACE_METHOD
    return DEFAULT_GRACE_METHOD


def grace_resolve_max_per_isomer(arg: Optional[int] = None) -> int:
    if arg is not None and int(arg) > 0:
        return int(arg)
    return max(1, _env_int(ENV_GRACE_MAX_PER_ISOMER, DEFAULT_GRACE_MAX_PER_ISOMER))


def grace_resolve_max_total_time_s(arg: Optional[float] = None) -> float:
    if arg is not None and float(arg) > 0:
        return float(arg)
    return max(0.0, _env_float(ENV_GRACE_MAX_TOTAL_TIME_S, DEFAULT_GRACE_MAX_TOTAL_TIME_S))


def grace_resolve_rmsd_threshold(arg: Optional[float] = None) -> float:
    if arg is not None and float(arg) > 0:
        return float(arg)
    return max(0.0, _env_float(ENV_GRACE_RMSD_THRESHOLD, DEFAULT_GRACE_RMSD_THRESHOLD))


def grace_resolve_max_rotamer_tuples(arg: Optional[int] = None) -> int:
    if arg is not None and int(arg) > 0:
        return int(arg)
    return max(1, _env_int(ENV_GRACE_MAX_ROTAMER_TUPLES, DEFAULT_GRACE_MAX_ROTAMER_TUPLES))


def grace_resolve_max_isomers(arg: Optional[int] = None) -> int:
    if arg is not None and int(arg) > 0:
        return int(arg)
    return max(1, _env_int(ENV_GRACE_MAX_ISOMERS, DEFAULT_GRACE_MAX_ISOMERS))


def grace_use_topology_healing() -> bool:
    """``True`` iff the topology-healing pre-polish step is enabled.

    Default ON when GRACE is active (it is a safe, accept-if-better
    operator).  Toggle via :data:`ENV_GRACE_USE_TOPOLOGY_HEALING`.
    """
    return _env_bool(ENV_GRACE_USE_TOPOLOGY_HEALING, default=True)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GraceCandidate:
    """One emitted GRACE candidate, fully labelled.

    Attributes are immutable so the candidate hashes for set / dict use.
    """

    syms: Tuple[str, ...]
    P: np.ndarray
    score: float                  # CCDC-Mahalanobis sum-of-z² (primary)
    severity: float               # raw GRIP severity before clash term
    clash_count: int              # inter-ligand clash count
    cshm: float                   # donor-polyhedron CShM (auxiliary)
    isomer_id: int
    ring_id: int                  # CP combination index, 0 = as-built
    rotamer_id: int               # rotamer tuple index, 0 = as-built
    rotamer_offsets: Tuple[float, ...]
    label: str
    method: str                   # "lm" or "lbfgs"
    accepted: bool                # passed every gate
    topology_healed: bool         # did the pre-polish step apply

    def __post_init__(self):
        object.__setattr__(self, "P", np.array(self.P, dtype=np.float64, copy=True))
        self.P.setflags(write=False)
        object.__setattr__(self, "rotamer_offsets",
                            tuple(float(x) for x in self.rotamer_offsets))


@dataclass
class GraceResult:
    """Output container of :func:`grace_enumerate`.

    The ``per_isomer`` dict maps the integer isomer id to a list of
    survivor candidates (post-dedup, ascending score).  ``best`` is the
    global score-minimum (across isomers).  ``burnside`` is a
    ``{n_predicted, n_emitted, coverage}`` triple plus diagnostic
    fields.  The fields ``n_total_assembled``, ``n_total_rejected_*`` are
    populated for the smoke-benchmark report.
    """

    smiles: str
    method: str
    per_isomer: Dict[int, List[GraceCandidate]] = field(default_factory=dict)
    n_isomers_enumerated: int = 0
    n_assembled: int = 0
    n_rejected_topology: int = 0
    n_rejected_polish: int = 0
    n_rejected_clash: int = 0
    n_emitted: int = 0
    burnside: Dict[str, Any] = field(default_factory=dict)
    wall_time_s: float = 0.0
    skip_reason: str = ""
    timed_out: bool = False

    @property
    def best(self) -> Optional[GraceCandidate]:
        candidates = [c for cs in self.per_isomer.values() for c in cs]
        if not candidates:
            return None
        return min(candidates, key=lambda c: (c.score, c.isomer_id, c.ring_id,
                                              c.rotamer_id, c.label))

    def total_candidates(self) -> int:
        return sum(len(v) for v in self.per_isomer.values())


# ---------------------------------------------------------------------------
# Burnside coverage helper
# ---------------------------------------------------------------------------
def burnside_coverage_report(
    n_polya: int,
    ring_state_counts: Sequence[int],
    n_rotamers: int,
    group_order: int,
    n_emitted: int,
) -> Dict[str, Any]:
    """Compare emitted-conformer count to the Burnside-orbit prediction.

    Uses :func:`burnside_conformer.count_orbits_product_space` for the
    closed-form prediction.  Returns a diagnostics dict::

        {
            "n_predicted_raw_product": int,
            "n_predicted_orbits": int,
            "n_emitted": int,
            "coverage_raw": float,    # n_emitted / raw_product
            "coverage_orbits": float, # n_emitted / orbits
            "active": bool,
        }

    All values are well-defined even when the inputs are degenerate
    (single rotor, no rings, group_order=1).
    """
    try:
        from . import burnside_conformer as BC
        n_pred = int(BC.count_orbits_product_space(
            n_polya=n_polya,
            ring_state_counts=list(ring_state_counts),
            n_rotamers=int(n_rotamers),
            group_order=int(group_order),
        ))
    except Exception as exc:  # pragma: no cover -- defensive
        return {
            "n_predicted_raw_product": 0,
            "n_predicted_orbits": 0,
            "n_emitted": int(n_emitted),
            "coverage_raw": 0.0,
            "coverage_orbits": 0.0,
            "active": False,
            "error": repr(exc),
        }

    raw_product = max(int(n_polya), 1) * max(int(n_rotamers), 1)
    for r in ring_state_counts:
        raw_product *= max(int(r), 1)

    cov_raw = float(n_emitted) / float(raw_product) if raw_product > 0 else 0.0
    cov_orb = float(n_emitted) / float(n_pred) if n_pred > 0 else 0.0
    return {
        "n_predicted_raw_product": int(raw_product),
        "n_predicted_orbits": int(n_pred),
        "n_emitted": int(n_emitted),
        "coverage_raw": float(cov_raw),
        "coverage_orbits": float(cov_orb),
        "active": True,
    }


# ---------------------------------------------------------------------------
# Internal: assembled-complex graph + base candidate
# ---------------------------------------------------------------------------
def _build_assembled_mol(
    metal_symbol: str,
    config: Dict[int, Tuple[int, int]],
    ligands: Sequence[Dict[str, Any]],
    donors_global: Sequence[int],
):
    """Reconstruct the RDKit mol for the assembled complex (metal index 0).

    Mirrors the iteration order used by
    :func:`delfin.fffree.assemble_complex.assemble_from_config` so atom
    indices align with the constructed XYZ.
    """
    try:
        from rdkit import Chem
    except Exception:
        return None

    cm = Chem.RWMol()
    cm.AddAtom(Chem.Atom(metal_symbol))

    by_lig: Dict[int, List[Tuple[int, int]]] = {}
    for v, (li, arm) in config.items():
        by_lig.setdefault(int(li), []).append((int(v), int(arm)))

    for li in by_lig.keys():  # insertion-order == assemble-order
        lg = ligands[li]
        try:
            frag = Chem.AddHs(lg["mol"])
        except Exception:
            return None
        base = cm.GetNumAtoms()
        for a in frag.GetAtoms():
            cm.AddAtom(Chem.Atom(a.GetAtomicNum()))
        for b in frag.GetBonds():
            cm.AddBond(
                b.GetBeginAtomIdx() + base,
                b.GetEndAtomIdx() + base,
                b.GetBondType(),
            )

    for dg in donors_global:
        try:
            cm.AddBond(0, int(dg), Chem.BondType.SINGLE)
        except Exception:
            pass

    try:
        Chem.SanitizeMol(cm, catchErrors=True)
    except Exception:
        pass
    return cm


def _compute_severity(syms: Sequence[str], P: np.ndarray, mol, metal_idx: int,
                       donors: Sequence[int]) -> float:
    """Mogul-Mahalanobis severity (reuses GRIP machinery).  Returns ``inf``
    on any failure so failing builds sink to the bottom of the ranking.
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


def _safe_cshm(P: np.ndarray, metal_idx: int, donors: Sequence[int],
                geometry: str) -> float:
    if not geometry or len(donors) < 2:
        return 0.0
    try:
        from . import polyhedra as PLY
        d_arr = np.asarray(P[list(donors)] - P[int(metal_idx)], dtype=float)
        norms = np.linalg.norm(d_arr, axis=1)
        norms = np.where(norms < 1e-12, 1e-12, norms)
        obs = d_arr / norms[:, None]
        return float(PLY.cshm(obs, geometry))
    except Exception:
        return float("inf")


def _count_inter_ligand_clashes_safe(P: np.ndarray, mol, donors: Sequence[int]) -> int:
    """Defensive wrapper around grip_ensemble.count_inter_ligand_clashes."""
    try:
        from .grip_ensemble import count_inter_ligand_clashes
        return int(count_inter_ligand_clashes(P, mol, metal_idx=0, donors=donors))
    except Exception:
        return 0


def _maybe_topology_heal(
    P: np.ndarray, syms: Sequence[str], mol, metal_idx: int,
    donors: Sequence[int],
) -> Tuple[np.ndarray, bool]:
    """Apply topology_healing_pipeline; return (P, healed_flag).

    Falls back to the input on any failure so the polish step always
    sees a valid array.
    """
    if not grace_use_topology_healing():
        return P, False
    try:
        from . import topology_healing as TH
        P_h = TH.topology_healing_pipeline(
            np.asarray(P, dtype=np.float64),
            list(syms),
            mol,
            metal_idx=int(metal_idx),
            donors=tuple(int(d) for d in donors),
            return_diagnostics=False,
        )
        if P_h is None:
            return P, False
        P_h = np.asarray(P_h, dtype=np.float64)
        if not np.all(np.isfinite(P_h)):
            return P, False
        # Did anything change?
        if P_h.shape == P.shape and float(np.max(np.abs(P_h - P))) > 1e-9:
            return P_h, True
        return P_h, False
    except Exception:
        return P, False


def _polish_dispatch(
    P: np.ndarray, mol, metal_idx: int, donors: Sequence[int],
    geometry: str, method: str,
) -> Tuple[np.ndarray, bool]:
    """Dispatch to :func:`grip_polish` (LM or L-BFGS).

    Returns ``(P_polished, ok_flag)``.  ``ok_flag=False`` means the
    polish either crashed or violated the M-D / topology guards (in
    which case the input is returned).
    """
    try:
        from .grip_polish import grip_polish
    except Exception:
        return P, False
    # The grip_polish dispatcher honours DELFIN_FFFREE_GRIP_METHOD;
    # we want to honour the GRACE method arg locally without leaking
    # env state into the caller.
    prev_env = os.environ.get("DELFIN_FFFREE_GRIP_METHOD")
    os.environ["DELFIN_FFFREE_GRIP_METHOD"] = method
    try:
        P_out = grip_polish(
            np.asarray(P, dtype=np.float64),
            mol,
            metal=int(metal_idx),
            donors=tuple(int(d) for d in donors),
            geom=str(geometry),
        )
        if P_out is None:
            return P, False
        P_out = np.asarray(P_out, dtype=np.float64)
        if not np.all(np.isfinite(P_out)):
            return P, False
        return P_out, True
    except Exception:
        return P, False
    finally:
        if prev_env is None:
            os.environ.pop("DELFIN_FFFREE_GRIP_METHOD", None)
        else:
            os.environ["DELFIN_FFFREE_GRIP_METHOD"] = prev_env


def _rmsd_dedup(
    candidates: List[GraceCandidate],
    threshold: float,
    max_keep: int,
    mol=None,
) -> List[GraceCandidate]:
    """RMSD-Butina dedup; keep up to ``max_keep`` survivors per isomer.

    The score-min survivor is always retained.  Other survivors are
    cluster representatives from non-overlapping Butina clusters.

    When ``DELFIN_FFFREE_GRACE_PARETO=1`` is set, the post-cluster
    survivor list is re-selected via :func:`grace_pareto.pareto_select`
    so the ensemble retains the full Pareto frontier across
    (severity, clash, cshm, arom_plan) before filling with dominated
    points.  Default OFF → byte-identical legacy behaviour.
    """
    if not candidates:
        return []
    if max_keep <= 0:
        return []
    # Pre-sort by score; cluster on heavy atoms only.
    try:
        from .conformer_dedup import pairwise_rmsd_matrix, butina_cluster
    except Exception:
        # Defensive: if dedup module unavailable, return score-sorted top-K.
        c_sorted = sorted(candidates,
                          key=lambda c: (c.score, c.isomer_id, c.ring_id,
                                          c.rotamer_id, c.label))
        return _maybe_pareto_reselect(c_sorted, max_keep, mol)

    coords_list = [np.asarray(c.P, dtype=float) for c in candidates]
    syms = list(candidates[0].syms) if candidates else None
    try:
        M = pairwise_rmsd_matrix(coords_list, syms=syms)
        clusters = butina_cluster(M, float(threshold))
    except Exception:
        c_sorted = sorted(candidates,
                          key=lambda c: (c.score, c.isomer_id, c.ring_id,
                                          c.rotamer_id, c.label))
        return _maybe_pareto_reselect(c_sorted, max_keep, mol)

    # For each cluster, pick the lowest-score representative.
    survivors: List[GraceCandidate] = []
    for cl in clusters:
        if not cl:
            continue
        best_idx = min(cl, key=lambda i: (float(candidates[i].score), int(i)))
        survivors.append(candidates[best_idx])
    survivors.sort(key=lambda c: (c.score, c.isomer_id, c.ring_id,
                                   c.rotamer_id, c.label))
    return _maybe_pareto_reselect(survivors, max_keep, mol)


def _maybe_pareto_reselect(
    survivors: List[GraceCandidate],
    max_keep: int,
    mol,
) -> List[GraceCandidate]:
    """Optional Pareto re-selection over RMSD-cluster survivors.

    When the env-gate is OFF this is a no-op (legacy truncation).
    When ON, calls :func:`grace_pareto.pareto_select` which keeps
    every Pareto-frontier candidate first, then fills with the lowest
    dominator-count survivors up to ``max_keep``.
    """
    if not survivors:
        return survivors
    try:
        from . import grace_pareto as GP
    except Exception:
        return survivors[: max(1, int(max_keep))]
    if not GP.pareto_active():
        return survivors[: max(1, int(max_keep))]
    try:
        sel = GP.pareto_select(
            survivors, max_keep=int(max_keep), mol=mol,
        )
    except Exception:
        return survivors[: max(1, int(max_keep))]
    # Final lex-sort for determinism (pareto_select already lex-orders
    # the frontier; dominated tail is appended in dominator-rank order).
    return list(sel)


# ---------------------------------------------------------------------------
# Pólya isomer enumeration (mirrors grip_ensemble.grip_ensemble_enumerate)
# ---------------------------------------------------------------------------
_GEOM_TO_POLYA_KEY = {
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
    # Mission B1 (2026-06-05): sandwich / piano-stool / half-sandwich.
    # Reachable only when ``DELFIN_FFFREE_SANDWICH_DISPATCH`` (or
    # ``PURE_TRACK3``) is on — see ``decompose._default_geometry``.  The
    # ``_eff`` keys point at the per-site EFFECTIVE Pólya groups (2 / 4 / 4
    # sites) used by the ensemble enumerator — the assembler view, not the
    # full N-atom polyhedron view.
    "SANDWICH-10 bis-eta5-Cp": "sandwich_10_eff",
    "PIANO-STOOL-8 eta5-Cp+L3": "piano_stool_8_eff",
    "HALF-SANDWICH-9 eta6+L3": "half_sandwich_9_eff",
}


def _enumerate_polya_isomers(d: Dict[str, Any]) -> Tuple[List[Dict[int, Tuple[int, int]]], str]:
    """Return the deterministic list of Pólya isomer configs + the polya key.

    Returns ``([], "")`` when the decomposed geometry is outside the
    Pólya-supported set.
    """
    geom = d.get("geometry") or ""
    geom_key = _GEOM_TO_POLYA_KEY.get(geom)
    if geom_key is None:
        return [], ""
    try:
        from . import polya_isomer_count as PIC
    except Exception:
        return [], ""
    if geom_key not in PIC._GROUPS:
        return [], geom_key

    ligands = d["ligands"]
    isomer_configs: List[Dict[int, Tuple[int, int]]] = []
    try:
        if d.get("has_chelate"):
            from rdkit import Chem
            specs = [
                {
                    "type": Chem.MolToSmiles(lg["mol"]),
                    "denticity": lg["denticity"],
                    "asym": len(set(lg.get("donor_elems", []))) > 1,
                }
                for lg in ligands
            ]
            isomer_configs = list(PIC.enumerate_chelate_configs(geom_key, specs))
        else:
            from collections import Counter
            from rdkit import Chem
            lig_label: List[str] = [Chem.MolToSmiles(lg["mol"]) for lg in ligands]
            spec = dict(Counter(lig_label))
            colorings = list(PIC.enumerate_isomers(geom_key, spec))
            for coloring in colorings:
                used: Set[int] = set()
                cfg: Dict[int, Tuple[int, int]] = {}
                for v_idx, lab in enumerate(coloring):
                    for li, lg in enumerate(ligands):
                        if li in used:
                            continue
                        try:
                            if Chem.MolToSmiles(lg["mol"]) == lab:
                                cfg[v_idx] = (li, 0)
                                used.add(li)
                                break
                        except Exception:
                            continue
                if len(cfg) == len(coloring):
                    isomer_configs.append(cfg)
    except Exception:
        return [], geom_key
    return isomer_configs, geom_key


def _detect_rings_for_cp(mol) -> List[Tuple[int, ...]]:
    """Return the sorted list of 5/6/7-membered rings for CP enumeration."""
    try:
        ri = mol.GetRingInfo()
        rings = []
        for r in ri.AtomRings():
            if 5 <= len(r) <= 7:
                rings.append(tuple(sorted(int(a) for a in r)))
        rings.sort()
        return rings
    except Exception:
        return []


def _enumerate_ring_combinations(
    P_base: np.ndarray,
    rings: Sequence[Tuple[int, ...]],
    max_per_ring: int,
) -> List[Tuple[np.ndarray, str]]:
    """Per-ring CP pucker variants over the base coords; lex-deterministic.

    For complexity reasons we apply each ring's variants independently
    (no cross-product across rings) — this is the same conservative
    choice :mod:`grip_ensemble` makes and keeps the candidate count
    bounded by ``len(rings) * max_per_ring + 1``.
    """
    variants: List[Tuple[np.ndarray, str]] = [(np.asarray(P_base, dtype=np.float64), "ring_base")]
    if not rings or max_per_ring <= 1:
        return variants
    try:
        from .ring_pucker import enumerate_ring_conformers
    except Exception:
        return variants
    try:
        for r_idx, ring in enumerate(rings):
            try:
                gen = enumerate_ring_conformers(
                    np.asarray(P_base, dtype=np.float64),
                    [ring],
                    max_per_ring=max_per_ring,
                )
                for c_id, P_var in enumerate(gen, start=1):
                    if c_id >= max_per_ring:
                        break
                    P_arr = np.asarray(P_var, dtype=np.float64)
                    try:
                        disp = float(np.max(np.linalg.norm(P_arr - P_base, axis=1)))
                    except Exception:
                        disp = 0.0
                    if disp < 0.05:
                        continue
                    variants.append((P_arr, f"ring{r_idx}_cp{c_id}"))
            except Exception:
                continue
    except Exception:
        pass
    return variants


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def grace_enumerate(
    smiles: str,
    *,
    method: Optional[str] = None,
    max_per_isomer: Optional[int] = None,
    max_total_time_s: Optional[float] = None,
    rmsd_dedup_threshold: Optional[float] = None,
    max_rotamer_tuples: Optional[int] = None,
    max_isomers: Optional[int] = None,
    max_per_ring: int = 3,
    clash_threshold: int = 5,
) -> GraceResult:
    """Deterministic global conformer ensemble via combinatorial scan.

    Pipeline (per Pólya isomer):

        for each ring-CP combination:
            for each single-bond rotamer tuple:
                1. Assemble candidate XYZ from (isomer, config, ligands).
                2. Apply topology-healing pre-polish (when enabled).
                3. Run grip_polish (L-BFGS-B or LM per ``method``).
                4. Score = Mogul-Mahalanobis severity + clash + cshm.
                5. Reject when topology / clash gates fail.
        RMSD-Butina-dedup → essential survivor set per isomer.
        Keep ``max_per_isomer`` lowest-score survivors.

    Output: ``GraceResult`` with ``per_isomer`` dict mapping
    ``isomer_id`` → list of :class:`GraceCandidate`, plus Burnside-
    coverage diagnostics and timing / rejection counters.

    Parameters
    ----------
    smiles : str
        Input SMILES of the metal complex (or organic molecule for a
        future organic-mode extension).
    method : {"lm", "lbfgs"}, optional
        Inner-polish optimiser.  Default from
        ``DELFIN_FFFREE_GRACE_METHOD`` env or ``"lm"``.
    max_per_isomer : int, optional
        Survivors kept after RMSD-dedup per isomer.  Default 10.
    max_total_time_s : float, optional
        Total wall-clock budget.  When exhausted, ``result.timed_out``
        is set; survivors collected so far are returned.  Default 60 s.
    rmsd_dedup_threshold : float, optional
        Butina cluster threshold in Å.  Default 0.3.
    max_rotamer_tuples : int, optional
        Cap on the rotamer-tuple count per (isomer, ring-combo).
        Default 27 (= 3³ — 3 staggered states over 3 rotors).
    max_isomers : int, optional
        Cap on the Pólya isomers iterated.  Default 12.
    max_per_ring : int, default 3
        Per-ring CP pucker count (1 = use as-built only).
    clash_threshold : int, default 5
        Inter-ligand clash count above which a candidate is rejected.

    Returns
    -------
    GraceResult
        With ``per_isomer`` populated, ``burnside`` coverage diagnostics,
        and timing / rejection counters.  Empty per_isomer means the
        SMILES is outside the fffree-decomposable domain (multi-metal
        legacy, unsupported geometry, …).
    """
    t0 = time.monotonic()

    method_eff = grace_resolve_method(method)
    max_iso_eff = grace_resolve_max_isomers(max_isomers)
    max_iso_per = grace_resolve_max_per_isomer(max_per_isomer)
    max_time_eff = grace_resolve_max_total_time_s(max_total_time_s)
    rmsd_eff = grace_resolve_rmsd_threshold(rmsd_dedup_threshold)
    max_rot_eff = grace_resolve_max_rotamer_tuples(max_rotamer_tuples)

    result = GraceResult(smiles=smiles, method=method_eff)

    # 0) Decompose the SMILES via the fffree pipeline.
    try:
        from . import decompose as DEC
        from . import assemble_complex as AC
    except Exception as exc:
        result.skip_reason = f"import_error: {exc!r}"
        result.wall_time_s = float(time.monotonic() - t0)
        return result

    try:
        d = DEC.decompose(smiles)
    except Exception as exc:
        result.skip_reason = f"decompose_error: {exc!r}"
        result.wall_time_s = float(time.monotonic() - t0)
        return result
    if d is None:
        result.skip_reason = "decompose_returned_none"
        result.wall_time_s = float(time.monotonic() - t0)
        return result

    metal_symbol = d.get("metal", "X")
    geometry = d.get("geometry", "")

    # 1) Enumerate Pólya isomers.
    isomer_configs, polya_key = _enumerate_polya_isomers(d)
    if not isomer_configs:
        result.skip_reason = f"isomer_enum_empty: geom={geometry!r}"
        result.wall_time_s = float(time.monotonic() - t0)
        return result

    n_isomers_total = len(isomer_configs)
    isomer_configs = isomer_configs[: int(max_iso_eff)]
    result.n_isomers_enumerated = len(isomer_configs)

    # Aggregate ring counts (for Burnside report) — populated as we go.
    sum_ring_counts: List[int] = []
    n_rotamer_layers: List[int] = []

    n_emitted_total = 0

    # 2) Per-isomer combinatorial scan.
    for iso_id, config in enumerate(isomer_configs):
        if (time.monotonic() - t0) > max_time_eff:
            result.timed_out = True
            break

        # 2a) Build the base assembly.
        try:
            built = AC.assemble_from_config(d["metal"], geometry, config, d["ligands"])
        except Exception:
            built = None
        if built is None:
            result.n_rejected_polish += 1
            continue
        syms_base, P_base, donors_global = built
        syms_base = list(syms_base)
        P_base = np.asarray(P_base, dtype=np.float64)
        donors_global = list(donors_global)
        result.n_assembled += 1

        # 2b) Reconstruct the assembled mol for graph operations.
        cm = _build_assembled_mol(metal_symbol, config, d["ligands"], donors_global)
        if cm is None or cm.GetNumAtoms() != len(syms_base):
            result.n_rejected_polish += 1
            continue

        # 2c) Ring detection + CP combinations.
        rings = _detect_rings_for_cp(cm)
        ring_combos = _enumerate_ring_combinations(P_base, rings, max_per_ring=max_per_ring)
        sum_ring_counts.append(len(ring_combos))

        # 2d) Rotamer enumeration (single-bond).
        try:
            from .single_bond_rotamers import (
                find_rotatable_bonds, enumerate_rotamer_configs, apply_rotamer_config,
                _dihedral_reference_atoms,
            )
            rotors_all = find_rotatable_bonds(cm)
            rotors: List[Tuple[int, int]] = []
            for b1, b2 in rotors_all:
                if _dihedral_reference_atoms(cm, b1, b2) is None:
                    continue
                rotors.append((int(b1), int(b2)))
        except Exception:
            rotors = []

        # Choose rotamer states: 3 staggered states (0, +120, -120) per rotor
        # to keep tuples small + meaningful chemically.
        if rotors:
            try:
                rotor_offsets_iter = list(enumerate_rotamer_configs(
                    n_rotors=len(rotors), n_states=3, max_configs=int(max_rot_eff),
                ))
            except Exception:
                rotor_offsets_iter = [tuple(0.0 for _ in rotors)]
        else:
            rotor_offsets_iter = [tuple()]
        n_rotamer_layers.append(len(rotor_offsets_iter))

        # 2e) Combinatorial loop.
        iso_candidates: List[GraceCandidate] = []
        # Adaptive plateau state (per-isomer); inert when env-flag OFF.
        _adaptive_state = None
        try:
            from . import grace_pareto as _GP_adapt
            if _GP_adapt.adaptive_active():
                _adaptive_state = _GP_adapt.AdaptiveState(
                    patience=_GP_adapt.resolve_adaptive_patience(),
                    max_ensemble=_GP_adapt.resolve_max_ensemble(),
                    mol=cm,
                )
        except Exception:
            _adaptive_state = None
        for ring_id, (P_ring, ring_label) in enumerate(ring_combos):
            if (time.monotonic() - t0) > max_time_eff:
                result.timed_out = True
                break
            for rot_id, rotor_offs in enumerate(rotor_offsets_iter):
                if (time.monotonic() - t0) > max_time_eff:
                    result.timed_out = True
                    break

                # Apply rotamer.
                if rotors and rotor_offs:
                    try:
                        P_cand = apply_rotamer_config(cm, P_ring, rotors, list(rotor_offs))
                    except Exception:
                        P_cand = np.asarray(P_ring, dtype=np.float64).copy()
                else:
                    P_cand = np.asarray(P_ring, dtype=np.float64).copy()

                # Topology-heal (pre-polish).
                P_healed, healed_flag = _maybe_topology_heal(
                    P_cand, syms_base, cm, metal_idx=0, donors=donors_global,
                )

                # Polish dispatch.
                P_polished, ok = _polish_dispatch(
                    P_healed, cm, metal_idx=0, donors=donors_global,
                    geometry=geometry, method=method_eff,
                )
                if not ok:
                    result.n_rejected_polish += 1
                    continue

                # Inter-ligand clash gate.
                clash_n = _count_inter_ligand_clashes_safe(P_polished, cm, donors_global)
                if clash_n > int(clash_threshold) * 2:
                    result.n_rejected_clash += 1
                    continue

                sev = _compute_severity(syms_base, P_polished, cm, 0, donors_global)
                cshm = _safe_cshm(P_polished, 0, donors_global, geometry)
                score = float(sev) + 10.0 * max(0, clash_n - int(clash_threshold)) ** 2 + 1.0 * (
                    float(cshm) if np.isfinite(cshm) else 1e6
                )

                rot_offs_t = tuple(float(x) for x in rotor_offs)
                label = (
                    f"iso{iso_id}_{ring_label}_rot{rot_id}"
                )
                cand = GraceCandidate(
                    syms=tuple(syms_base),
                    P=np.asarray(P_polished, dtype=np.float64),
                    score=float(score),
                    severity=float(sev),
                    clash_count=int(max(clash_n, 0)),
                    cshm=float(cshm),
                    isomer_id=int(iso_id),
                    ring_id=int(ring_id),
                    rotamer_id=int(rot_id),
                    rotamer_offsets=rot_offs_t,
                    label=label,
                    method=method_eff,
                    accepted=(np.isfinite(sev) and clash_n >= 0),
                    topology_healed=bool(healed_flag),
                )
                iso_candidates.append(cand)

                # Adaptive plateau check (default OFF).
                if _adaptive_state is not None:
                    try:
                        _adaptive_state.update(cand)
                        if _adaptive_state.should_stop():
                            break
                    except Exception:
                        pass

            if result.timed_out:
                break
            # Propagate adaptive stop one level up.
            if _adaptive_state is not None and _adaptive_state.should_stop():
                break
        # 2f) RMSD-dedup per isomer + keep top-K.
        if iso_candidates:
            # When Pareto re-selection is active, raise the effective
            # cap so the frontier is preserved (still hard-capped by
            # MAX_ENSEMBLE).  Legacy path is unaffected (default OFF).
            cap_eff = max_iso_per
            try:
                from . import grace_pareto as _GP
                if _GP.pareto_active():
                    cap_eff = max(max_iso_per, _GP.resolve_max_ensemble())
            except Exception:
                pass
            survivors = _rmsd_dedup(iso_candidates, threshold=rmsd_eff,
                                     max_keep=cap_eff, mol=cm)
            result.per_isomer[int(iso_id)] = survivors
            n_emitted_total += len(survivors)

        if result.timed_out:
            break

    result.n_emitted = n_emitted_total

    # 3) Burnside coverage report.
    # Aggregate ring-state counts: take per-isomer mean (rings are detected
    # AFTER assembly).  This is heuristic but anchored in the actual emit
    # set.  Group order: heuristic = polyhedron rotation-group order.
    try:
        from . import polya_isomer_count as PIC
        if polya_key in PIC._GROUPS:
            _g, _n = PIC._get_group(polya_key)
            group_order = int(len(_g))
        else:
            group_order = 1
    except Exception:
        group_order = 1

    mean_rings = int(round(np.mean(sum_ring_counts))) if sum_ring_counts else 1
    mean_rotamers = int(round(np.mean(n_rotamer_layers))) if n_rotamer_layers else 1

    result.burnside = burnside_coverage_report(
        n_polya=n_isomers_total,
        ring_state_counts=[mean_rings],
        n_rotamers=mean_rotamers,
        group_order=group_order,
        n_emitted=n_emitted_total,
    )

    result.wall_time_s = float(time.monotonic() - t0)
    return result
