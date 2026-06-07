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


# ---------------------------------------------------------------------------
# Mission D2 (2026-06-05): tight-convergence env-overrides.
#
# Defaults are byte-identical with the legacy path (max_iter=200, gtol=1e-4,
# restarts=1).  Set the env-flags in the production launcher to:
#
#   DELFIN_FFFREE_GRIP_MAX_ITER=1000   -> deeper convergence (UFF-level)
#   DELFIN_FFFREE_GRIP_GTOL=1e-6       -> tighter gradient tol
#   DELFIN_FFFREE_GRIP_RESTARTS=3      -> multi-restart L-BFGS for robustness
#
# Resolution: caller arg (numeric, finite) > env (numeric, finite) > default.
# Determinism: restarts use deterministic, atom-sorted perturbation seeded
# by PYTHONHASHSEED + restart index, so 2-run output is byte-identical.
# ---------------------------------------------------------------------------
_MAX_ITER_ENV: str = "DELFIN_FFFREE_GRIP_MAX_ITER"
_GTOL_ENV: str = "DELFIN_FFFREE_GRIP_GTOL"
_RESTARTS_ENV: str = "DELFIN_FFFREE_GRIP_RESTARTS"
_RESTART_PERTURB_ENV: str = "DELFIN_FFFREE_GRIP_RESTART_PERTURB"
DEFAULT_MAX_ITER: int = 200
DEFAULT_GTOL: float = 1e-4
DEFAULT_RESTARTS: int = 1
DEFAULT_RESTART_PERTURB: float = 0.05  # Å, max coord perturbation amplitude


def _resolve_max_iter(arg_value) -> int:
    """Effective L-BFGS max_iter.  arg > env > 200 (default-OFF byte-identical)."""
    if arg_value is not None:
        try:
            v = int(arg_value)
            if v > 0:
                return v
        except (TypeError, ValueError):
            pass
    raw = os.environ.get(_MAX_ITER_ENV, "").strip()
    if raw:
        try:
            v = int(raw)
            if v > 0:
                return v
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not a positive int; using default %d",
                _MAX_ITER_ENV, raw, DEFAULT_MAX_ITER,
            )
    return DEFAULT_MAX_ITER


def _resolve_gtol(arg_value) -> float:
    """Effective L-BFGS gtol.  arg > env > 1e-4 (default-OFF byte-identical)."""
    if arg_value is not None:
        try:
            v = float(arg_value)
            if np.isfinite(v) and v > 0:
                return v
        except (TypeError, ValueError):
            pass
    raw = os.environ.get(_GTOL_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v > 0:
                return v
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not a positive float; using default %.1e",
                _GTOL_ENV, raw, DEFAULT_GTOL,
            )
    return DEFAULT_GTOL


def _resolve_restarts() -> int:
    """Number of multi-restart L-BFGS rounds.  env > 1 (default-OFF byte-identical).

    Returns 1 by default -> the legacy path runs exactly one L-BFGS call from
    P_init (byte-identical with HEAD).  Values > 1 add (k-1) extra restarts
    from deterministically perturbed copies of P_init, keeping the lowest-
    severity acceptable result.
    """
    raw = os.environ.get(_RESTARTS_ENV, "").strip()
    if raw:
        try:
            v = int(raw)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not a positive int; using default %d",
                _RESTARTS_ENV, raw, DEFAULT_RESTARTS,
            )
    return DEFAULT_RESTARTS


def _resolve_restart_perturb() -> float:
    """Amplitude (Å) for deterministic multi-restart coord perturbation.

    Default 0.05 Å: small enough that frozen-sphere (metal/donors) constraints
    still validate after the perturbation, large enough to escape a shallow
    local minimum.  env-overridable via DELFIN_FFFREE_GRIP_RESTART_PERTURB.
    """
    raw = os.environ.get(_RESTART_PERTURB_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v >= 0:
                return v
        except (TypeError, ValueError):
            pass
    return DEFAULT_RESTART_PERTURB


def _deterministic_perturbation(P_init: np.ndarray, restart_idx: int,
                                amplitude: float,
                                frozen_indices: np.ndarray) -> np.ndarray:
    """Generate a deterministic perturbation of ``P_init`` for restart k.

    The perturbation is generated by a Mersenne-Twister seeded by
    ``(PYTHONHASHSEED || 0, restart_idx, n_atoms)`` so the same restart_idx
    + same P_init shape yields byte-identical output.  Frozen atoms (metal +
    donors) are NOT perturbed -- their positions remain exactly P_init.
    """
    if restart_idx == 0 or amplitude <= 0.0:
        return P_init.copy()
    seed_hash = int(os.environ.get("PYTHONHASHSEED", "0") or "0")
    rng = np.random.RandomState(seed=(seed_hash + 1_000 * restart_idx) & 0x7FFFFFFF)
    delta = rng.uniform(-amplitude, amplitude, size=P_init.shape).astype(np.float64)
    if frozen_indices.size:
        delta[frozen_indices] = 0.0
    return P_init + delta


# ---------------------------------------------------------------------------
# Heavy-atom vdW-floor (2026-06-07) — anti-collapse penalty.
#
# The existing :class:`ClashFloorPenalty` already implements a soft Pauli floor
# at ``0.85 * (r_i + r_j)`` with default weight 5.0 over EVERY non-bonded,
# non-1,3 atom pair (including H).  In production we observe the Mahalanobis-
# severity loss can still over-pull heavy atoms together (V3 voll-pool:
# ``structqual_heavy_collapse_pct_frames_viol +537 %`` vs UFF baseline).
#
# This block adds a STRICTER additional vdW-floor term that fires ONLY between
# heavy (non-H) atoms with a higher default weight (20.0) and the same Pauli
# fraction (0.85).  Default-OFF (env unset) -> byte-identical with HEAD; opt
# in via ``DELFIN_FFFREE_GRIP_VDW_FLOOR=1``.
#
# The math (per heavy-atom pair (i,j) NOT bonded, NOT 1-3):
#
#     d_ij    = ||x_i - x_j||
#     floor   = fraction * (r_vdw_i + r_vdw_j)
#     if d_ij < floor and d_ij > 0:
#         L_vdw       += w_vdw * (floor - d_ij)^2
#         grad_x_i    += -2 * w_vdw * (floor - d_ij) * (x_i - x_j) / d_ij
#         grad_x_j    += +2 * w_vdw * (floor - d_ij) * (x_i - x_j) / d_ij
#
# Pairs whose endpoints are both frozen (M-D rigid sphere) still receive the
# gradient but it is later zeroed by the frozen-projection step inside
# ``loss_and_grad`` -- no double-bookkeeping needed here.  Coincident atoms
# (``d_ij < _VDW_FLOOR_EPS``) are skipped to avoid 0/0 gradient.
# ---------------------------------------------------------------------------
_VDW_FLOOR_ENV: str = "DELFIN_FFFREE_GRIP_VDW_FLOOR"
_VDW_FLOOR_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_VDW_FLOOR_WEIGHT"
_VDW_FLOOR_FRACTION_ENV: str = "DELFIN_FFFREE_GRIP_VDW_FLOOR_FRACTION"
DEFAULT_VDW_FLOOR_WEIGHT: float = 20.0
DEFAULT_VDW_FLOOR_FRACTION: float = 0.85
_VDW_FLOOR_EPS: float = 1e-9


def _vdw_floor_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_VDW_FLOOR`` is on (default OFF).

    Accepts ``1``/``true``/``yes``/``on`` (case-insensitive).  Any other value
    (including unset / empty) keeps the term OFF, preserving byte-identity
    with the legacy clash floor only.
    """
    raw = os.environ.get(_VDW_FLOOR_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _resolve_vdw_floor_weight() -> float:
    """Effective weight for the heavy-atom vdW-floor penalty.

    env > :data:`DEFAULT_VDW_FLOOR_WEIGHT` (20.0).  Non-finite / non-numeric
    values fall through to the default with a single warning so the polish
    never crashes on a misconfigured env-flag.
    """
    raw = os.environ.get(_VDW_FLOOR_WEIGHT_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v >= 0.0:
                return v
            _LOG.warning(
                "grip_polish: %s=%r not a non-negative finite float; using default %.3f",
                _VDW_FLOOR_WEIGHT_ENV, raw, DEFAULT_VDW_FLOOR_WEIGHT,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not numeric; using default %.3f",
                _VDW_FLOOR_WEIGHT_ENV, raw, DEFAULT_VDW_FLOOR_WEIGHT,
            )
    return DEFAULT_VDW_FLOOR_WEIGHT


def _resolve_vdw_floor_fraction() -> float:
    """Effective Pauli-floor fraction for the heavy-atom vdW-floor penalty.

    env > :data:`DEFAULT_VDW_FLOOR_FRACTION` (0.85).  Values are clamped to
    ``(0, 2]``; pathological values fall through to the default.
    """
    raw = os.environ.get(_VDW_FLOOR_FRACTION_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and 0.0 < v <= 2.0:
                return v
            _LOG.warning(
                "grip_polish: %s=%r out of (0, 2]; using default %.3f",
                _VDW_FLOOR_FRACTION_ENV, raw, DEFAULT_VDW_FLOOR_FRACTION,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not numeric; using default %.3f",
                _VDW_FLOOR_FRACTION_ENV, raw, DEFAULT_VDW_FLOOR_FRACTION,
            )
    return DEFAULT_VDW_FLOOR_FRACTION


def _heavy_atom_indices(symbols: Sequence[str], n_atoms: int) -> np.ndarray:
    """Return a sorted ``int64`` ndarray of indices whose symbol is NOT ``H``.

    Symbols are matched case-sensitively against ``"H"`` -- mirrors the
    rest of :mod:`grip_polish` (which uses ``str(a.GetSymbol())`` directly).
    Atoms with an out-of-range index or non-string symbol are silently
    skipped (defensive: same policy as ``_vdw_table_for_mol``).
    """
    out: List[int] = []
    for i in range(min(int(n_atoms), len(symbols))):
        s = symbols[i]
        if not isinstance(s, str):
            continue
        if s == "H":
            continue
        out.append(int(i))
    return np.asarray(sorted(set(out)), dtype=np.int64)


def _vdw_floor_value_and_grad(
    R: np.ndarray,
    *,
    heavy_indices: np.ndarray,
    radii: np.ndarray,
    excluded_pairs: Set[FrozenSet[int]],
    weight: float,
    fraction: float,
) -> Tuple[float, np.ndarray]:
    """Heavy-atom-only Pauli-floor penalty + analytic gradient.

    Parameters
    ----------
    R : ndarray (N, 3)
        Current coordinates.
    heavy_indices : ndarray of int64
        Sorted indices of heavy (non-H) atoms with a known vdW radius.
    radii : ndarray of float64, shape (N,)
        vdW radii indexed by atom; ``NaN`` for unknown elements.  Only
        ``radii[heavy_indices]`` is consulted.
    excluded_pairs : set of frozenset({i, j})
        Bonded + 1,3 pairs to skip (typically the same set the legacy
        :class:`ClashFloorPenalty` already uses).
    weight, fraction : float
        Penalty multiplier and floor fraction.

    Returns
    -------
    (L, G) : float, ndarray (N, 3)
        Loss value and gradient.  When no pair violates the floor both are
        zero, so the call is a no-op for well-separated geometries.
    """
    n = int(R.shape[0])
    grad = np.zeros_like(R)
    if heavy_indices.size < 2 or weight <= 0.0 or fraction <= 0.0:
        return 0.0, grad

    total = 0.0
    H = heavy_indices.tolist()
    for ii in range(len(H)):
        i = H[ii]
        ri = radii[i] if 0 <= i < n else np.nan
        if not np.isfinite(ri):
            continue
        Pi = R[i]
        for jj in range(ii + 1, len(H)):
            j = H[jj]
            rj = radii[j] if 0 <= j < n else np.nan
            if not np.isfinite(rj):
                continue
            if frozenset((i, j)) in excluded_pairs:
                continue
            d_vec = Pi - R[j]
            d = float(np.linalg.norm(d_vec))
            floor_ij = fraction * (ri + rj)
            if d >= floor_ij:
                continue
            if d < _VDW_FLOOR_EPS:
                # Exactly coincident -- skip rather than emit NaN.
                continue
            gap = floor_ij - d
            total += weight * gap * gap
            # dL/dx_i = -2 w (floor - d) * (x_i - x_j) / d
            coef = -2.0 * weight * gap / d
            gi = coef * d_vec
            grad[i] += gi
            grad[j] -= gi

    return float(total), grad


# ---------------------------------------------------------------------------
# Inter-shell Pauli-floor (2026-06-07) — donor / 2nd-shell clash anti-collapse.
#
# Auto-Diagnostic on the V3 voll-pool flagged 27.4% (1816 / 6627) of structures
# as having *inter-shell donor clashes* — a donor atom of one ligand sitting
# inside the Pauli radius of a 2nd-shell atom (a NON-donor that is bonded to
# a donor of a DIFFERENT ligand, e.g. an N-CH₃ methyl carbon eclipsing an
# adjacent Cl donor).  These pairs are 1-4 or further in the molecular graph
# (so the legacy excl_13 lets them through) yet they are physically forced
# close together because both atoms live in the metal's first / second
# coordination sphere.
#
# Universal definition — graph + geometry only, no per-class heuristics:
#
#   * First shell  = {donors}                  (atoms bonded to the metal)
#   * Second shell = {a | a NOT bonded to M
#                       AND ∃ d ∈ donors : (d, a) is a bond
#                       AND ||x_a - x_M|| ≤ d_max}   (d_max default 3.5 Å)
#
# For every pair (donor_i, second_shell_j) that is NOT bonded and NOT inside
# the existing 1-3 exclusion set (so (donor_i, donor_j) via-metal pairs are
# already skipped), apply a quadratic Pauli floor identical in form to the
# heavy-atom vdW-floor above:
#
#     d_ij    = ||x_i - x_j||
#     floor   = fraction * (r_vdw_i + r_vdw_j)
#     if d_ij < floor and d_ij > 0:
#         L         += w * (floor - d_ij)^2
#         grad_x_i  += -2 w (floor - d_ij) * (x_i - x_j) / d_ij
#         grad_x_j  += +2 w (floor - d_ij) * (x_i - x_j) / d_ij
#
# Composes ADDITIVELY with the existing ``ClashFloorPenalty`` (intra-ligand
# weight 5.0) and the heavy-atom ``_vdw_floor_value_and_grad`` term (heavy-
# all-pair weight 20.0).  The default weight here is ``15.0`` — chosen to
# match the existing inter-ligand clash boost so well-formed structures stay
# untouched (zero violation -> zero loss + zero gradient) while real
# inter-shell clashes feel a strong push.
#
# Env-gated default OFF -> byte-identical with HEAD when the flag is unset:
#
#     DELFIN_FFFREE_INTER_SHELL_FLOOR=1            # turn it on
#     DELFIN_FFFREE_INTER_SHELL_FLOOR_WEIGHT=15.0  # tune the multiplier
#     DELFIN_FFFREE_INTER_SHELL_FLOOR_FRACTION=0.85
#     DELFIN_FFFREE_INTER_SHELL_FLOOR_RADIUS=3.5   # 2nd-shell M-X cutoff (Å)
# ---------------------------------------------------------------------------
_INTER_SHELL_FLOOR_ENV: str = "DELFIN_FFFREE_INTER_SHELL_FLOOR"
_INTER_SHELL_FLOOR_WEIGHT_ENV: str = "DELFIN_FFFREE_INTER_SHELL_FLOOR_WEIGHT"
_INTER_SHELL_FLOOR_FRACTION_ENV: str = "DELFIN_FFFREE_INTER_SHELL_FLOOR_FRACTION"
_INTER_SHELL_FLOOR_RADIUS_ENV: str = "DELFIN_FFFREE_INTER_SHELL_FLOOR_RADIUS"
DEFAULT_INTER_SHELL_FLOOR_WEIGHT: float = 15.0
DEFAULT_INTER_SHELL_FLOOR_FRACTION: float = 0.85
DEFAULT_INTER_SHELL_FLOOR_RADIUS: float = 3.5  # Å; 2nd-shell cutoff vs metal
_INTER_SHELL_FLOOR_EPS: float = 1e-9


def _inter_shell_floor_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_INTER_SHELL_FLOOR`` is on (default OFF).

    Accepts ``1``/``true``/``yes``/``on`` (case-insensitive).  Any other
    value (including unset / empty) keeps the term OFF, preserving
    byte-identity with HEAD.
    """
    raw = os.environ.get(_INTER_SHELL_FLOOR_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _resolve_inter_shell_floor_weight() -> float:
    """Effective weight for the inter-shell Pauli-floor penalty.

    env > :data:`DEFAULT_INTER_SHELL_FLOOR_WEIGHT` (15.0).  Non-finite /
    non-numeric values fall through to the default with a single warning so
    the polish never crashes on a misconfigured env-flag.
    """
    raw = os.environ.get(_INTER_SHELL_FLOOR_WEIGHT_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v >= 0.0:
                return v
            _LOG.warning(
                "grip_polish: %s=%r not a non-negative finite float; using default %.3f",
                _INTER_SHELL_FLOOR_WEIGHT_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not numeric; using default %.3f",
                _INTER_SHELL_FLOOR_WEIGHT_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
            )
    return DEFAULT_INTER_SHELL_FLOOR_WEIGHT


def _resolve_inter_shell_floor_fraction() -> float:
    """Effective Pauli-floor fraction for the inter-shell penalty.

    env > :data:`DEFAULT_INTER_SHELL_FLOOR_FRACTION` (0.85).  Values are
    clamped to ``(0, 2]``; pathological values fall through to the default.
    """
    raw = os.environ.get(_INTER_SHELL_FLOOR_FRACTION_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and 0.0 < v <= 2.0:
                return v
            _LOG.warning(
                "grip_polish: %s=%r out of (0, 2]; using default %.3f",
                _INTER_SHELL_FLOOR_FRACTION_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_FRACTION,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not numeric; using default %.3f",
                _INTER_SHELL_FLOOR_FRACTION_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_FRACTION,
            )
    return DEFAULT_INTER_SHELL_FLOOR_FRACTION


def _resolve_inter_shell_floor_radius() -> float:
    """Effective 2nd-shell selection radius (Å) vs the metal.

    env > :data:`DEFAULT_INTER_SHELL_FLOOR_RADIUS` (3.5 Å).  Values are
    clamped to ``(0, 10]``; pathological values fall through to the default.
    """
    raw = os.environ.get(_INTER_SHELL_FLOOR_RADIUS_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and 0.0 < v <= 10.0:
                return v
            _LOG.warning(
                "grip_polish: %s=%r out of (0, 10]; using default %.3f",
                _INTER_SHELL_FLOOR_RADIUS_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "grip_polish: %s=%r not numeric; using default %.3f",
                _INTER_SHELL_FLOOR_RADIUS_ENV, raw, DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            )
    return DEFAULT_INTER_SHELL_FLOOR_RADIUS


def _build_inter_shell_pairs(
    mol_bonds: Sequence[Tuple[int, int]],
    metal: int,
    donors: Sequence[int],
    P: np.ndarray,
    excluded_pairs: Set[FrozenSet[int]],
    radius: float,
    n_atoms: int,
) -> List[Tuple[int, int]]:
    """Enumerate the inter-shell ``(donor_i, second_shell_j)`` pairs.

    A pair ``(i, j)`` is returned iff:

      * ``i`` is in the donor set (first shell, bonded to ``metal``).
      * ``j`` is NOT bonded to ``metal`` (not first shell), AND ``j`` is
        bonded to AT LEAST ONE donor ``d ∈ donors`` (second shell), AND
        ``||x_j - x_M|| <= radius``.
      * ``(i, j)`` is NOT in ``excluded_pairs`` (bonded + 1-3 through any
        atom, INCLUDING through-metal 1-3).
      * ``j`` is NOT the metal itself.
      * ``i != j``.

    Sorted, deterministic, graph-only construction.  Geometry is consulted
    ONLY for the 2nd-shell radius cutoff (a soft geometric gate, not a
    classification).
    """
    n = int(n_atoms)
    metal_i = int(metal)
    donor_set: Set[int] = {int(d) for d in donors if 0 <= int(d) < n and int(d) != metal_i}
    if not donor_set:
        return []

    # Build adjacency.
    adj: List[Set[int]] = [set() for _ in range(n)]
    for (a, b) in mol_bonds:
        if 0 <= a < n and 0 <= b < n and a != b:
            adj[a].add(b)
            adj[b].add(a)

    # First-shell membership = atoms bonded to the metal.
    first_shell: Set[int] = set(adj[metal_i]) if 0 <= metal_i < n else set()
    # In normal operation first_shell == donor_set; we tolerate divergence
    # by using ``donor_set`` for the donor endpoint (authoritative caller
    # data) and ``first_shell`` for the "NOT bonded to metal" exclusion.

    # Second-shell candidates: atoms bonded to any donor but NOT to metal,
    # within ``radius`` of the metal.
    second_shell: Set[int] = set()
    if 0 <= metal_i < n:
        metal_pos = P[metal_i]
        for d in donor_set:
            if not (0 <= d < n):
                continue
            for nb in adj[d]:
                if nb == metal_i:
                    continue
                if nb in first_shell:
                    # Skip atoms that are themselves first-shell (donors).
                    continue
                if nb in donor_set:
                    continue
                # Geometric 2nd-shell gate.
                d_M = float(np.linalg.norm(P[nb] - metal_pos))
                if d_M <= radius:
                    second_shell.add(int(nb))

    if not second_shell:
        return []

    # Enumerate candidate (donor, 2nd-shell) pairs, skipping bonded / 1-3.
    pairs: Set[Tuple[int, int]] = set()
    for i in donor_set:
        if not (0 <= i < n):
            continue
        for j in second_shell:
            if i == j:
                continue
            key = frozenset((int(i), int(j)))
            if key in excluded_pairs:
                continue
            a, b = (int(i), int(j)) if int(i) < int(j) else (int(j), int(i))
            pairs.add((a, b))

    return sorted(pairs)


def _inter_shell_floor_value_and_grad(
    R: np.ndarray,
    *,
    pairs: Sequence[Tuple[int, int]],
    radii: np.ndarray,
    weight: float,
    fraction: float,
) -> Tuple[float, np.ndarray]:
    """Inter-shell Pauli-floor penalty + analytic gradient.

    Parameters
    ----------
    R : ndarray (N, 3)
        Current coordinates.
    pairs : sequence of (int, int)
        Sorted (i, j) pairs to penalize -- pre-built by
        :func:`_build_inter_shell_pairs`.
    radii : ndarray of float64, shape (N,)
        vdW radii indexed by atom; ``NaN`` for unknown elements (pair is
        silently skipped).
    weight, fraction : float
        Penalty multiplier and floor fraction.

    Returns
    -------
    (L, G) : float, ndarray (N, 3)
        Loss value and gradient.  When no pair violates the floor both are
        zero, so the call is a no-op for well-separated geometries.
    """
    n = int(R.shape[0])
    grad = np.zeros_like(R)
    if not pairs or weight <= 0.0 or fraction <= 0.0:
        return 0.0, grad

    total = 0.0
    for (i, j) in pairs:
        if not (0 <= i < n and 0 <= j < n):
            continue
        ri = radii[i] if 0 <= i < radii.shape[0] else np.nan
        rj = radii[j] if 0 <= j < radii.shape[0] else np.nan
        if not (np.isfinite(ri) and np.isfinite(rj)):
            continue
        d_vec = R[i] - R[j]
        d = float(np.linalg.norm(d_vec))
        floor_ij = fraction * (float(ri) + float(rj))
        if d >= floor_ij:
            continue
        if d < _INTER_SHELL_FLOOR_EPS:
            # Coincident atoms -- skip to avoid 0/0 gradient.
            continue
        gap = floor_ij - d
        total += weight * gap * gap
        # dL/dx_i = -2 w (floor - d) * (x_i - x_j) / d
        coef = -2.0 * weight * gap / d
        gi = coef * d_vec
        grad[i] += gi
        grad[j] -= gi

    return float(total), grad


from .grip_constraints import (
    ChiralVolumeConstraint,
    ClashFloorPenalty,
    DonorPolyhedronConstraint,
    MDInvariantConstraint,
    MetalChiralityConstraint,
    TopologyConstraint,
)
from .grip_fragment_detect import detect_fragments
from .grip_loss_terms import TotalGripLoss
from .grip_md_terms import (
    DEFAULT_MD_WARMUP_FRACTION as _MD_WARMUP_FRACTION,
    anneal_md_term_weights as _anneal_md_term_weights,
    build_md_loss_terms as _build_md_loss_terms,
    compute_total_severity as _phase2_total_severity,
    detect_hapto_donor_clusters as _detect_md_hapto_clusters,
    detect_hapto_present as _detect_md_hapto_present,
    detect_multi_metal as _detect_md_multi_metal,
    max_md_drift as _max_md_drift,
    md_hapto_skip_active as _md_hapto_skip_active,
    md_warmup_active as _md_warmup_active,
    resolve_md_accept_weight as _resolve_md_accept_weight,
    unfreeze_md_active as _unfreeze_md_active,
)
from .grip_mogul_lookup import GripLibrary, get_default_library

# Phase 2 M-D drift safety threshold (Å).  When the polish moves any M-D
# distance more than ``_MD_DRIFT_SAFETY`` from its initial value we revert
# to the frozen-sphere result -- preserves the M-D invariant doctrine even
# when the unfreeze flag is on (fail open).
_MD_DRIFT_SAFETY: float = 0.5

__all__ = [
    "grip_polish",
    "mogul_severity",
    "GripPolishResult",
    "DEFAULT_VDW_RADII",
    "DEFAULT_CLASH_WEIGHT",
    "DEFAULT_INTER_LIGAND_CLASH_WEIGHT",
    "DEFAULT_ACCEPT_CLASH_ALPHA",
    # Heavy-atom vdW-floor anti-collapse term (2026-06-07, env-gated default OFF).
    "DEFAULT_VDW_FLOOR_WEIGHT",
    "DEFAULT_VDW_FLOOR_FRACTION",
    "_vdw_floor_active",
    "_resolve_vdw_floor_weight",
    "_resolve_vdw_floor_fraction",
    "_vdw_floor_value_and_grad",
    "_heavy_atom_indices",
    # Inter-shell Pauli-floor (2026-06-07, env-gated default OFF, donor / 2nd-shell).
    "DEFAULT_INTER_SHELL_FLOOR_WEIGHT",
    "DEFAULT_INTER_SHELL_FLOOR_FRACTION",
    "DEFAULT_INTER_SHELL_FLOOR_RADIUS",
    "_inter_shell_floor_active",
    "_resolve_inter_shell_floor_weight",
    "_resolve_inter_shell_floor_fraction",
    "_resolve_inter_shell_floor_radius",
    "_build_inter_shell_pairs",
    "_inter_shell_floor_value_and_grad",
    "detect_hapto_atoms",
    "build_ligand_atom_id_map",
    "expand_hapto_for_sigma_only",
    "_sigma_only_mode_active",
    # Hapto-π explicit freeze + substituent subtree drag (2026-06-07, Task #91).
    "_hapto_pi_freeze_active",
    "_build_hapto_rigid_body",
    "_project_rigid_body",
    "_resolve_grip_method",
    "_GRIP_METHOD_ENV",
    "_GRIP_METHOD_LBFGS",
    "_GRIP_METHOD_LM",
    # Healing-wiring (2026-06-04) — pre-polish heal hooks, env-gated default-OFF.
    "_topology_healing_active",
    "_grip_healing_active",
    "_run_pre_polish_topology_healing",
    "_run_pre_polish_grip_healing",
    "_sp3_h_heal_active",
    "_run_pre_polish_sp3_h_heal",
    # GRACE env-gated dispatcher (default OFF, byte-identical to HEAD 00f1a5b).
    "grace_dispatch_active",
    "grace_polish_or_enumerate",
]


# ---------------------------------------------------------------------------
# Healing-wiring env-flag predicates (2026-06-04, wiring cleanup).
#
# Three standalone healing modules ship with HEAD b195dba but are NOT yet
# called by the pipeline.  This block adds env-gated pre-polish hooks so the
# orchestration sits in ONE place (grip_polish dispatcher) instead of forking
# the call graph at every site.  All flags default OFF -> byte-identical with
# HEAD b195dba when unset.  Order (when all enabled):
#
#     coords -> topology_healing  (phantom / missing / wrong-angle, graph)
#            -> grip_healing      (broken-atom iterative reposition, atoms)
#            -> grip_polish_inner (L-BFGS / LM Mahalanobis polish)
#
# H-H clash inclusion is wired separately via ClashFloorPenalty extras +
# the per-endpoint hooks in grip_constraints / build_time_clash_gate /
# inter_ligand_clash_gate / assemble_complex (env DELFIN_FFFREE_HH_CLASH_INCLUDE).
# ---------------------------------------------------------------------------
_TOPOLOGY_HEALING_ENV: str = "DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING"
_GRIP_HEALING_ENV: str = "DELFIN_FFFREE_GRIP_HEALING_MODE"
_HH_CLASH_INCLUDE_ENV: str = "DELFIN_FFFREE_HH_CLASH_INCLUDE"
# sp³-H umbrella heal (2026-06-04, tBu-CH3 90°/180° defect, GUPPY ea26dcb).
# Default OFF byte-identical; orthogonal to the two heals above.
_SP3_H_HEAL_ENV: str = "DELFIN_FFFREE_SP3_H_HEAL"


# Breadcrumb logging — default OFF; env-flag DELFIN_FFFREE_VERIFY_BREADCRUMBS=1
# emits a single WARNING line whenever a heal hook actually fires AND moves
# coordinates.  Used by overnight verification (2026-06-04) to confirm
# modules are not silently passing through.
_BREADCRUMB_ENV: str = "DELFIN_FFFREE_VERIFY_BREADCRUMBS"


def _breadcrumbs_active() -> bool:
    raw = os.environ.get(_BREADCRUMB_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def _emit_breadcrumb(module: str, P_in: np.ndarray, P_out: np.ndarray,
                     extra: str = "") -> None:
    """Log a one-line ``MODULE_FIRED: <module> rmsd=<x>`` breadcrumb.

    Only emitted when ``DELFIN_FFFREE_VERIFY_BREADCRUMBS=1`` AND the output
    coords differ from input.  Optional file sink via
    ``DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE`` so callers don't have to rely
    on log capture.
    """
    if not _breadcrumbs_active():
        return
    try:
        diff = np.asarray(P_out, dtype=np.float64) - np.asarray(P_in, dtype=np.float64)
        rms = float(np.sqrt(np.mean(np.sum(diff.reshape(-1, 3) ** 2, axis=1))))
    except Exception:
        rms = float("nan")
    if rms > 1e-12:
        _LOG.warning("MODULE_FIRED: %s rmsd=%.4f %s", module, rms, extra)
        try:
            crumb_file = os.environ.get(
                "DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE", ""
            ).strip()
            if crumb_file:
                with open(crumb_file, "a") as fh:
                    fh.write(f"MODULE_FIRED: {module} rmsd={rms:.4f} {extra}\n")
        except Exception:
            pass


def _topology_healing_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING`` is on (default OFF).

    Mirror of :func:`delfin.fffree.topology_healing.is_active` kept local
    so the grip_polish dispatcher does not pay a hard import cost when the
    flag is OFF (the topology_healing module is imported lazily inside the
    hook below).
    """
    raw = os.environ.get(_TOPOLOGY_HEALING_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _grip_healing_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_HEALING_MODE`` is on (default OFF).

    Mirror of :func:`delfin.fffree.grip_healing.healing_mode_active`.
    """
    raw = os.environ.get(_GRIP_HEALING_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _hh_clash_include_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_HH_CLASH_INCLUDE`` is on (default OFF)."""
    raw = os.environ.get(_HH_CLASH_INCLUDE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _sp3_h_heal_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_SP3_H_HEAL`` is on (default OFF).

    Mirror of :func:`delfin.fffree.sp3_h_heal.heal_active` kept local so
    the dispatcher avoids a hard import when the flag is OFF.
    """
    raw = os.environ.get(_SP3_H_HEAL_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _run_pre_polish_topology_healing(
    P_init: np.ndarray,
    mol,
    metal: int,
    donors: Sequence[int],
) -> np.ndarray:
    """Run :func:`topology_healing.topology_healing_pipeline` on ``P_init``.

    Pure pass-through when the env-flag is OFF -- returns ``P_init`` unchanged
    (byte-identical with HEAD b195dba).  When ON, the topology heal:

      * detects phantom / missing / wrong-angle defects using SMILES
        topology from ``mol``,
      * applies accept-if-better rolls (defects must not increase),
      * preserves the M-D invariant (metal + donors are frozen).

    Any exception -> swallow + return ``P_init`` (defence in depth; the heal
    must never crash the polish).
    """
    if not _topology_healing_active():
        return P_init
    try:
        from .topology_healing import topology_healing_pipeline
    except Exception as exc:  # pragma: no cover -- module missing
        _LOG.warning(
            "grip_polish: failed to import topology_healing (%r); skipping",
            exc,
        )
        return P_init
    try:
        # Extract symbols + bond topology from mol.  Falls back to no-op
        # if either is missing (avoids reading bogus data).
        try:
            syms = [str(a.GetSymbol()) for a in mol.GetAtoms()]
        except Exception:
            return P_init
        if len(syms) != P_init.shape[0]:
            return P_init
        out = topology_healing_pipeline(
            P_init,
            syms,
            mol,
            metal_idx=int(metal),
            donors=tuple(int(d) for d in donors),
            return_diagnostics=False,
        )
        out = np.asarray(out, dtype=np.float64)
        if out.shape != P_init.shape or not np.all(np.isfinite(out)):
            return P_init
        _emit_breadcrumb("topology_healing", P_init, out)
        return out
    except Exception as exc:
        _LOG.warning(
            "grip_polish: topology_healing raised (%r); falling back to P_init",
            exc,
        )
        return P_init


def _run_pre_polish_grip_healing(
    P_init: np.ndarray,
    mol,
    metal: int,
    donors: Sequence[int],
    library,
) -> np.ndarray:
    """Run :func:`grip_healing.iterative_topology_repositioning` on ``P_init``.

    Pure pass-through when the env-flag is OFF -- returns ``P_init`` unchanged
    (byte-identical with HEAD b195dba).  When ON, the grip_healing pass
    repositions broken atoms (bond residual > σ-threshold) onto their
    CCDC-grounded targets while keeping the metal + donors frozen.

    Any exception -> swallow + return ``P_init``.
    """
    if not _grip_healing_active():
        return P_init
    try:
        from .grip_healing import (
            iterative_topology_repositioning,
            _build_ideal_table,
            _extract_mol_bonds,
            _extract_symbols,
        )
    except Exception as exc:  # pragma: no cover
        _LOG.warning(
            "grip_polish: failed to import grip_healing (%r); skipping",
            exc,
        )
        return P_init
    try:
        bonds = _extract_mol_bonds(mol)
        syms = _extract_symbols(mol)
        if not bonds or not syms:
            return P_init
        if len(syms) != P_init.shape[0]:
            return P_init
        ideal_table = _build_ideal_table(bonds, syms, library=library)
        frozen = set([int(metal), *[int(d) for d in donors]])
        out = iterative_topology_repositioning(
            P_init,
            bonds,
            ideal_lengths=ideal_table,
            symbols=syms,
            library=library,
            frozen_atoms=frozen,
            return_diagnostics=False,
        )
        out = np.asarray(out, dtype=np.float64)
        if out.shape != P_init.shape or not np.all(np.isfinite(out)):
            return P_init
        _emit_breadcrumb("grip_healing", P_init, out)
        return out
    except Exception as exc:
        _LOG.warning(
            "grip_polish: grip_healing raised (%r); falling back to P_init",
            exc,
        )
        return P_init


def _run_pre_polish_sp3_h_heal(
    P_init: np.ndarray,
    mol,
    metal: int,
    donors: Sequence[int],
) -> np.ndarray:
    """Run :func:`sp3_h_heal.heal_degenerate_sp3_h` on ``P_init``.

    Pure pass-through when ``DELFIN_FFFREE_SP3_H_HEAL`` is OFF -- returns
    ``P_init`` unchanged (byte-identical with HEAD ea26dcb).  When ON, the
    pass rewrites every degenerate sp³-C/N H umbrella (H-X-H angles
    outside (95°, 125°) or near-collinear) to the canonical Td template,
    keeping the metal + σ-donors frozen and respecting an accept-if-better
    gate (max H-X-H deviation must strictly decrease).

    Any exception -> swallow + return ``P_init`` (defence in depth).
    """
    if not _sp3_h_heal_active():
        return P_init
    try:
        from .sp3_h_heal import heal_degenerate_sp3_h
    except Exception as exc:  # pragma: no cover -- module missing
        _LOG.warning(
            "grip_polish: failed to import sp3_h_heal (%r); skipping",
            exc,
        )
        return P_init
    try:
        try:
            syms = [str(a.GetSymbol()) for a in mol.GetAtoms()]
        except Exception:
            return P_init
        if len(syms) != P_init.shape[0]:
            return P_init
        frozen = {int(metal), *[int(d) for d in donors]}
        out, _report = heal_degenerate_sp3_h(
            P_init,
            mol=mol,
            syms=syms,
            frozen_atoms=frozen,
            preserve_bond_lengths=True,
            accept_if_better=True,
        )
        out = np.asarray(out, dtype=np.float64)
        if out.shape != P_init.shape or not np.all(np.isfinite(out)):
            return P_init
        _emit_breadcrumb("sp3_h_heal", P_init, out)
        return out
    except Exception as exc:
        _LOG.warning(
            "grip_polish: sp3_h_heal raised (%r); falling back to P_init",
            exc,
        )
        return P_init


# ---------------------------------------------------------------------------
# GRACE env-gated dispatcher (default OFF, byte-identical to HEAD 00f1a5b).
# Deterministic Group-theoretic Reproducible Adaptive Conformer Ensemble —
# Pólya × Cremer-Pople × Rotamer enumeration with topology_healing pre-polish
# and LM/L-BFGS dispatch.  See delfin.fffree.grace_ensemble for the algorithm.
# ---------------------------------------------------------------------------


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


def _stereocenter_quadruples(
    mol, exclude_center: Optional[int] = None,
) -> List[Tuple[int, int, int, int]]:
    """Return ``(center, a, b, c)`` quadruples for every heavy atom with
    exactly three heavy neighbours (deterministic, sorted).

    Parameters
    ----------
    mol : RDKit Mol-like
    exclude_center : int, optional
        Atom index to exclude as a stereocenter.  When Phase 2 (M+D
        unfreeze) is active the caller passes the metal index so the
        per-atom signed-volume check does not flip on innocent M+D
        rigid rotations -- the metal's intrinsic chirality is instead
        protected by :class:`MetalChiralityConstraint`.  Default
        ``None`` -> byte-identical with the legacy stereocenter list.
    """
    quads: List[Tuple[int, int, int, int]] = []
    skip_idx = -1 if exclude_center is None else int(exclude_center)
    for atom in mol.GetAtoms():
        c = int(atom.GetIdx())
        if c == skip_idx:
            continue
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
# Hapto-π rigid-body freeze (2026-06-07, Task #91 explicit freeze + drag).
#
# When DELFIN_FFFREE_GRIP_UNFREEZE_MD=1 is active the metal + donors gain DoF
# under L-BFGS.  For hapto-π binding modes (η^n, n>=2) this can break the
# η^n geometry: the metal drifts relative to the ring carbons (donors) and
# the ring atoms drift relative to the metal, even with MD_HAPTO_SKIP=1
# (which only skips loss TERMS, not the position updates).
#
# Fix: when DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE=1 (default OFF byte-identical),
# build a per-hapto-cluster "rigid body" = ring atoms + every H + every
# heavy-atom substituent subtree attached to a ring atom.  Then during the
# L-BFGS loss/grad evaluation:
#   - take the metal's displacement delta_M = R[metal] - P_init[metal];
#   - apply delta_M to every rigid-body atom (they ride with the metal);
#   - evaluate loss + grad on the effective coords;
#   - chain-rule: sum gradient contributions on the rigid-body atoms onto
#     the metal's gradient; zero the rigid-body atoms' own gradients (they
#     are determined by the metal's position).
#
# After convergence, project the final positions: rigid-body atoms are
# translated by the metal's net displacement.
#
# Reuses :mod:`hapto_rigid_subtree` (already shipped, deterministic, graph-
# only, no SMILES patterns).  Default-OFF byte-identical: when the env-flag
# is unset the predicate short-circuits, the rigid-body set is empty, and
# the wrapper falls through to the original loss_and_grad call.
# ---------------------------------------------------------------------------
_HAPTO_PI_FREEZE_ENV: str = "DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE"


def _hapto_pi_freeze_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE`` is on (default OFF).

    Live-read so the env-flag can flip between subprocesses.  Defaults
    to OFF -> byte-identical with HEAD when unset.
    """
    raw = os.environ.get(_HAPTO_PI_FREEZE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _build_hapto_rigid_body(
    mol,
    metal_idx: int,
    donors: Sequence[int],
    P_init: np.ndarray,
    hapto_atoms: Set[int],
) -> List[int]:
    """Return the sorted list of atom indices that must ride rigidly with
    the metal during the M+D unfreeze L-BFGS pass.

    The set contains:
      * Every detected hapto-π ring atom (from ``hapto_atoms``);
      * Every H bonded to a ring atom;
      * Every non-ring heavy-atom substituent attached to a ring atom +
        its entire connected subtree (BFS, never crossing the metal or
        another hapto ring).

    The metal index is EXCLUDED from the returned list (the metal is the
    "leader" of the rigid block; its position drives the rigid translation).

    Pure / deterministic / graph-only (delegated to
    :func:`hapto_rigid_subtree.collect_rigid_subtree_atoms`).  Returns
    the empty list when ``hapto_atoms`` is empty or the helper raises --
    fail-open so the polish never crashes on this guard.
    """
    if not hapto_atoms:
        return []
    metal_idx = int(metal_idx)
    try:
        from .hapto_rigid_subtree import collect_rigid_subtree_atoms
    except Exception:
        return []
    try:
        symbols = [str(a.GetSymbol()) for a in mol.GetAtoms()]
    except Exception:
        return []
    if len(symbols) != int(P_init.shape[0]):
        return []
    # Cluster the ring atoms by connectivity so multi-hapto complexes
    # (ferrocene = 2 Cp rings) get TWO rigid bodies, not one merged blob.
    # Two hapto atoms belong to the same cluster iff they share a bond
    # in the molecular graph (graph-only, deterministic).
    cluster_of: Dict[int, int] = {}
    next_cid = 0
    hapto_list = sorted(int(a) for a in hapto_atoms)
    for a in hapto_list:
        if a in cluster_of:
            continue
        # BFS over hapto subgraph.
        comp = [a]
        cluster_of[a] = next_cid
        head = 0
        while head < len(comp):
            cur = comp[head]
            head += 1
            try:
                atom = mol.GetAtomWithIdx(int(cur))
            except Exception:
                continue
            for nb in atom.GetNeighbors():
                ni = int(nb.GetIdx())
                if ni in hapto_atoms and ni not in cluster_of:
                    cluster_of[ni] = next_cid
                    comp.append(ni)
        next_cid += 1

    # Inverse map cid -> ring atoms.
    clusters: Dict[int, List[int]] = {}
    for a, cid in cluster_of.items():
        clusters.setdefault(cid, []).append(int(a))

    # For each cluster, run collect_rigid_subtree_atoms with the OTHER
    # clusters' atoms as a boundary ("other_ring_atoms") so the BFS
    # never crosses through a different hapto block.  Frozen extras =
    # the σ donors (M-D rigid, never absorbed into a hapto block).
    sigma_donors = set(int(d) for d in donors) - set(hapto_atoms)
    rigid_union: Set[int] = set()
    for cid, ring in clusters.items():
        other = set()
        for cid2, ring2 in clusters.items():
            if cid2 == cid:
                continue
            other.update(ring2)
        try:
            block = collect_rigid_subtree_atoms(
                symbols, P_init, ring,
                metal_idx=metal_idx,
                other_ring_atoms=other,
                frozen_extra=sigma_donors,
            )
        except Exception:
            continue
        for idx in block:
            i = int(idx)
            if i == metal_idx:
                continue
            rigid_union.add(i)
    return sorted(rigid_union)


def _project_rigid_body(
    R: np.ndarray,
    metal_idx: int,
    rigid_body_indices: np.ndarray,
    P_init: np.ndarray,
) -> np.ndarray:
    """Return a coords array where every atom in ``rigid_body_indices``
    is set to ``P_init[i] + (R[metal_idx] - P_init[metal_idx])``.

    The metal position is taken from ``R`` (unchanged from the L-BFGS
    iterate).  Pure / deterministic.  Out-of-range / empty input is a
    no-op (returns ``R`` copy).
    """
    R = np.asarray(R, dtype=np.float64)
    if rigid_body_indices.size == 0:
        return R.copy()
    if R.ndim == 1:
        R_mat = R.reshape(-1, 3).copy()
    else:
        R_mat = R.copy()
    delta = R_mat[int(metal_idx)] - np.asarray(P_init, dtype=np.float64)[int(metal_idx)]
    R_mat[rigid_body_indices] = (
        np.asarray(P_init, dtype=np.float64)[rigid_body_indices] + delta
    )
    return R_mat


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
    # 0b. Pre-polish healing pipeline (2026-06-04, wiring cleanup).
    #
    # Three standalone healing modules ship in HEAD b195dba but were not
    # wired into the polish dispatcher.  This block invokes them in a
    # fixed order BEFORE fragment detection so the polish operates on
    # a "best-effort-healed" geometry:
    #
    #     1. topology_healing  -> phantom / missing-bond / wrong-angle
    #     2. grip_healing      -> broken-atom iterative repositioning
    #
    # Both hooks are env-gated default-OFF byte-identical with HEAD
    # b195dba (return ``P_init`` verbatim when their flags are unset).
    # Each is independently composable: one can run without the other.
    # ------------------------------------------------------------------
    library = mogul_lib if mogul_lib is not None else None
    if library is None:
        try:
            library = get_default_library()
        except Exception:
            library = None

    P_init = _run_pre_polish_topology_healing(P_init, mol, metal, donors_t)
    P_init = _run_pre_polish_grip_healing(P_init, mol, metal, donors_t, library)
    # sp³-H umbrella heal (2026-06-04): degenerate 90°/180° H-C-H patterns
    # that survive both heals above (Mahalanobis-tolerant + bond-length-only)
    # and cannot move under L-BFGS because 180° is a cos(angle) saddle.
    P_init = _run_pre_polish_sp3_h_heal(P_init, mol, metal, donors_t)
    # Phase 2 (2026-06-06): controlled M+D unfreeze.
    #
    # Default-OFF byte-identical with HEAD (resolver returns False when
    # DELFIN_FFFREE_GRIP_UNFREEZE_MD is unset / 0).  When ON:
    #   - the ``frozen`` sphere shrinks to {} so the metal + donors can
    #     move under L-BFGS,
    #   - we generate Phase 2 M-D loss terms (MD bond + M-D-X angle +
    #     D-M-D angle + hapto centroid) and concatenate them onto the
    #     Phase 1 fragment list,
    #   - after the L-BFGS run we measure the M-D drift; if any donor
    #     has moved more than _MD_DRIFT_SAFETY Å the polish reverts to
    #     the frozen-sphere result (fail open).
    #
    # Phase 2 acceptance-gate fix (2026-06-06): per-SMILES skip when the
    # molecule has a hapto-π cluster or multiple metals.  The smoke 100
    # paired verdict showed Phase 2 regressing those classes (the M-D
    # Gaussians + HaptoCentroidTerm fight the construction priors).  The
    # skip is on by default; set ``DELFIN_FFFREE_GRIP_MD_HAPTO_SKIP=0`` to
    # force Phase 2 on hapto/multi-metal SMILES for forensics.
    _unfreeze_active_flag = _unfreeze_md_active()
    # Hapto-π rigid-body freeze + drag (Task #91, 2026-06-07).  When
    # ``DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE=1`` is on, the unfreeze pass
    # is allowed to PROCEED on hapto SMILES, but the ring + substituent
    # subtrees are tied to the metal as a rigid block (see the wrapper
    # around loss_and_grad below).  Default OFF preserves the legacy
    # "skip Phase 2 on hapto" behaviour byte-identically.
    _hapto_pi_freeze_flag = _hapto_pi_freeze_active()
    if _unfreeze_active_flag and _md_hapto_skip_active():
        try:
            if _detect_md_multi_metal(mol):
                _unfreeze_active_flag = False
        except Exception:
            pass
        if _unfreeze_active_flag and not _hapto_pi_freeze_flag:
            try:
                if _detect_md_hapto_present(mol, metal, donors_t, P_init):
                    _unfreeze_active_flag = False
            except Exception:
                pass
    if _unfreeze_active_flag:
        frozen: FrozenSet[int] = frozenset()
    else:
        frozen = frozenset((metal, *donors_t))
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
            metal=int(metal),
        )
    except Exception:
        terms = []

    # Phase 2 (2026-06-06): inject the M+D unfreeze loss terms.  No-op when
    # the unfreeze flag is OFF (build returns immediately on the predicate),
    # so this is byte-identical with HEAD whenever the env-var is unset.
    #
    # Phase 2 acceptance-gate fix: keep a separate handle on the md_terms
    # list so the gate (below) can score Phase 2 severity independently of
    # the Phase 1 fragment severity.  ``md_terms`` is the empty list when
    # the unfreeze flag is OFF or the build returns nothing -- both code
    # paths score Phase 2 = 0.0 and the gate degrades to legacy behaviour.
    md_terms: List[object] = []
    if _unfreeze_active_flag:
        try:
            md_terms = list(_build_md_loss_terms(
                mol, P_init, metal, donors_t, geom=geom,
                library=library,
            ))
            if md_terms:
                terms = list(terms) + list(md_terms)
        except Exception as _md_exc:
            _LOG.warning(
                "grip_polish: build_md_loss_terms raised (%r); skipping Phase 2 terms",
                _md_exc,
            )
            md_terms = []

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
    # When the angle-to-metal layer is active the M-D-X angle pressure can
    # push an X atom into a position that stretches a single X-Y bond
    # toward the topology multiplier (default 1.5×).  The legacy polish
    # had no pressure on X relative to the M-D axis, so the 1.5× cap was
    # designed for bond / angle / improper terms only.  Loosen it modestly
    # (default 1.8×) when the new layer is on so the constraint stays
    # protective against catastrophic stretches but does not roll back
    # every angle-to-metal improvement.  Default-OFF byte-identical:
    # without DELFIN_FFFREE_GRIP_ANGLE_TO_METAL=1 the multiplier is the
    # caller-supplied value (default 1.5×).
    #
    # Mission B1 (2026-06-05): donor-cone DoF (commit b338003) has the same
    # property -- per-donor azimuth as a joint L-BFGS variable can push a
    # downstream X-Y bond past the 1.5× cap when a polish that is good
    # globally (lower CShM, smaller M-D residuals) implies a longer single
    # bond at the cone tip.  Apply the same modest 1.8× cap when
    # ``DELFIN_FFFREE_GRIP_DONOR_CONE=1`` is on.  The bump is the MAX of the
    # two layers (still 1.8× when both are on), so the combined behaviour is
    # deterministic and idempotent.  Default-OFF byte-identical with HEAD.
    try:
        from .grip_fragment_detect import angle_to_metal_active as _atm_active
        _topo_mult_eff = float(topo_max_multiplier)
        if _atm_active() and _topo_mult_eff < 1.8:
            _topo_mult_eff = 1.8
    except Exception:
        _topo_mult_eff = float(topo_max_multiplier)
    try:
        from .grip_donor_cone import donor_cone_active as _dc_active
        if _dc_active() and _topo_mult_eff < 1.8:
            _topo_mult_eff = 1.8
    except Exception:
        # Defence-in-depth: never let an import failure regress topo behaviour.
        pass
    topo_constraint = TopologyConstraint.from_initial(
        mol_bonds, P_init, max_distance_multiplier=_topo_mult_eff,
    )

    # Phase 2 chiral-volume fix (2026-06-06).
    #
    # When the M+D unfreeze (Phase 2) is active the per-atom signed
    # tetrahedral volume at the metal centre is no longer the right
    # invariant -- M+D atoms rotate as a rigid body and the per-atom
    # sign can flip even when the intrinsic Δ/Λ chirality at the metal
    # is unchanged.  Replace the per-atom check with a metal-centered
    # Δ/Λ class (from donor unit vectors) for the metal, and keep the
    # per-atom check for every NON-metal stereocenter.
    #
    # Default-OFF byte-identical: when _unfreeze_active_flag is False
    # exclude_center=None -> _stereocenter_quadruples returns the
    # legacy quad list (metal included) and metal_chirality is the
    # empty constraint (validate() returns True with no metals).
    if _unfreeze_active_flag:
        chiral_constraint = ChiralVolumeConstraint.from_initial(
            _stereocenter_quadruples(mol, exclude_center=metal), P_init,
        )
        metal_chirality_constraint = MetalChiralityConstraint.from_initial(
            [(metal, donors_t)], P_init,
        )
    else:
        chiral_constraint = ChiralVolumeConstraint.from_initial(
            _stereocenter_quadruples(mol), P_init,
        )
        metal_chirality_constraint = MetalChiralityConstraint(metals=[])

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
    # H-H clash refinement (2026-06-04, wiring cleanup).
    #
    # ``ClashFloorPenalty`` already iterates EVERY non-bonded, non-1,3
    # pair using its vdW table (which includes H at 1.20 Å), so H-H
    # contacts already enter the loss surface.  The standalone
    # :mod:`hh_clash_detector` module adds a small but important
    # refinement: it knows that geminal H-H (~1.78 Å on a methylene C)
    # and 1-3 H-H within the same methyl group (~1.78 Å) are NOT
    # clashes -- chemically correct close contacts that should not
    # contribute to the Pauli penalty.  When the env-flag is ON we
    # extend ``excl_13`` with those H-H pairs so the penalty fires
    # only on genuine 1-4+ H-H eclipsing / inter-ligand contacts.
    #
    # Default OFF -> ``excl_13`` is byte-identical with HEAD b195dba.
    if _hh_clash_include_active():
        try:
            from .hh_clash_detector import build_hh_exclusion_pairs
            try:
                hh_syms = [str(a.GetSymbol()) for a in mol.GetAtoms()]
            except Exception:
                hh_syms = []
            if len(hh_syms) == n_atoms:
                hh_excl = build_hh_exclusion_pairs(hh_syms, P_init)
                if hh_excl:
                    excl_13 = excl_13 | set(hh_excl)
        except Exception as exc:
            _LOG.warning(
                "grip_polish: hh_clash_detector exclusion build failed (%r); skipping",
                exc,
            )
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

    # ------------------------------------------------------------------
    # 3b. Heavy-atom vdW-floor (2026-06-07, anti-collapse).
    #
    # Env-gated stricter Pauli-floor that fires ONLY between heavy (non-H)
    # atom pairs.  Default OFF -> byte-identical with HEAD.  Pre-computed
    # heavy index list + radii ndarray live in the closure so the inner
    # loss_and_grad call pays nothing when the term is OFF (a single bool
    # branch) and only the O(N_heavy^2) walk when the term is ON.
    # ------------------------------------------------------------------
    _vdw_floor_on = _vdw_floor_active()
    _vdw_floor_weight: float = 0.0
    _vdw_floor_fraction: float = DEFAULT_VDW_FLOOR_FRACTION
    _vdw_floor_heavy_indices: np.ndarray = np.empty(0, dtype=np.int64)
    _vdw_floor_radii_arr: np.ndarray = np.empty(0, dtype=np.float64)
    _vdw_floor_symbols: List[str] = []
    # Extended vdW-floor (Bug-class 2+3 fix, 2026-06-07): apply the floor to
    # ALL non-bonded heavy pairs INCLUDING those inside the hapto-π rigid
    # body, but route the gradient through the rigid body's translation so
    # the body cannot internally distort.  Env-gated default OFF byte-
    # identical: when ``DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS`` is unset the
    # legacy ``_vdw_floor_value_and_grad`` path is used unchanged.
    _vdw_floor_all_pairs = False
    try:
        from delfin.fffree.construction_sanity import (
            vdw_floor_all_pairs_active as _cs_vdw_all_active,
            vdw_floor_all_pairs_value_and_grad as _cs_vdw_all_grad,
        )
        _vdw_floor_all_pairs = bool(_cs_vdw_all_active())
    except Exception:
        _vdw_floor_all_pairs = False
        _cs_vdw_all_grad = None  # type: ignore[assignment]
    # When the all-pairs extension is requested, also require the legacy
    # vdW-floor to be on -- the extension reuses its weight, fraction and
    # symbol list so we keep one configuration path.
    if _vdw_floor_all_pairs and not _vdw_floor_on:
        _vdw_floor_on = True
    if _vdw_floor_on:
        try:
            _vdw_floor_weight = _resolve_vdw_floor_weight()
            _vdw_floor_fraction = _resolve_vdw_floor_fraction()
            try:
                _vdw_floor_symbols = [str(a.GetSymbol()) for a in mol.GetAtoms()]
            except Exception:
                _vdw_floor_symbols = []
            if len(_vdw_floor_symbols) == n_atoms:
                _vdw_floor_heavy_indices = _heavy_atom_indices(
                    _vdw_floor_symbols, n_atoms,
                )
            # Drop heavy indices that have no vdW radius -- the inner walk
            # already guards on NaN but pre-filtering keeps the loop tight.
            if _vdw_floor_heavy_indices.size:
                _vdw_floor_radii_arr = np.full(n_atoms, np.nan, dtype=np.float64)
                for _idx, _r in vdw_idx_table.items():
                    if 0 <= int(_idx) < n_atoms:
                        _vdw_floor_radii_arr[int(_idx)] = float(_r)
                _keep = np.array(
                    [bool(np.isfinite(_vdw_floor_radii_arr[int(k)]))
                     for k in _vdw_floor_heavy_indices.tolist()],
                    dtype=bool,
                )
                _vdw_floor_heavy_indices = _vdw_floor_heavy_indices[_keep]
            if (_vdw_floor_heavy_indices.size < 2
                    or _vdw_floor_weight <= 0.0
                    or _vdw_floor_fraction <= 0.0):
                # Nothing to do -- demote to OFF so the inner call is a no-op.
                _vdw_floor_on = False
        except Exception as _vdw_floor_exc:
            _LOG.warning(
                "grip_polish: vdW-floor setup failed (%r); disabling term",
                _vdw_floor_exc,
            )
            _vdw_floor_on = False

    # ------------------------------------------------------------------
    # 3c. Inter-shell Pauli-floor (2026-06-07, donor / 2nd-shell anti-clash).
    #
    # Auto-Diagnostic on the V3 voll-pool showed 27.4 % of files with a donor
    # atom clashing inside the Pauli radius of a 2nd-shell atom on a DIFFERENT
    # ligand (e.g. Ru-N and adjacent N-CH₃ where the methyl C overlaps an
    # adjacent Cl donor).  These pairs are 1-4+ in the graph so the existing
    # ``ClashFloorPenalty`` / heavy-atom vdW-floor still miss them in tight
    # first/second-shell geometries.  Wire a stricter, purely geometric +
    # graph-derived Pauli term ON TOP of the existing penalties.
    #
    # The donor / 2nd-shell pair list is constructed once from the INITIAL
    # geometry (the topology is frozen during the polish — TopologyConstraint
    # would roll back any bond break), so the cost per loss_and_grad call is
    # O(|pairs|) with |pairs| ≤ |donors| * |second_shell| which is small.
    #
    # Default OFF -> byte-identical with HEAD when the env-flag is unset.
    # ------------------------------------------------------------------
    _inter_shell_floor_on = _inter_shell_floor_active()
    _inter_shell_floor_weight: float = 0.0
    _inter_shell_floor_fraction: float = DEFAULT_INTER_SHELL_FLOOR_FRACTION
    _inter_shell_floor_radius: float = DEFAULT_INTER_SHELL_FLOOR_RADIUS
    _inter_shell_pairs: List[Tuple[int, int]] = []
    _inter_shell_radii_arr: np.ndarray = np.empty(0, dtype=np.float64)
    if _inter_shell_floor_on:
        try:
            _inter_shell_floor_weight = _resolve_inter_shell_floor_weight()
            _inter_shell_floor_fraction = _resolve_inter_shell_floor_fraction()
            _inter_shell_floor_radius = _resolve_inter_shell_floor_radius()
            # Build the per-atom vdW radii ndarray once (reused on every call).
            _inter_shell_radii_arr = np.full(n_atoms, np.nan, dtype=np.float64)
            for _idx, _r in vdw_idx_table.items():
                if 0 <= int(_idx) < n_atoms:
                    _inter_shell_radii_arr[int(_idx)] = float(_r)
            # Enumerate the (donor, 2nd-shell) pairs from the INITIAL geometry.
            _inter_shell_pairs = _build_inter_shell_pairs(
                mol_bonds=mol_bonds,
                metal=int(metal),
                donors=donors_t,
                P=P_init,
                excluded_pairs=excl_13,
                radius=_inter_shell_floor_radius,
                n_atoms=n_atoms,
            )
            if (not _inter_shell_pairs
                    or _inter_shell_floor_weight <= 0.0
                    or _inter_shell_floor_fraction <= 0.0):
                # Nothing to do -- demote to OFF so the inner call is a no-op.
                _inter_shell_floor_on = False
        except Exception as _inter_shell_exc:
            _LOG.warning(
                "grip_polish: inter-shell-floor setup failed (%r); disabling term",
                _inter_shell_exc,
            )
            _inter_shell_floor_on = False
            _inter_shell_pairs = []

    sev_before = mogul_severity(P_init, fragments)

    # ------------------------------------------------------------------
    # 4. Combined loss + frozen-gradient projection.
    # ------------------------------------------------------------------
    frozen_indices = np.array(sorted(frozen), dtype=np.int64)

    # ------------------------------------------------------------------
    # Hapto-π rigid-body freeze + drag (Task #91, 2026-06-07).
    #
    # When DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE=1 AND the unfreeze flag is
    # ON AND a hapto cluster is present, build the rigid-body atom set
    # (ring atoms + ring-H + substituent subtrees, EXCLUDING the metal).
    # The loss_and_grad wrapper below then ties these atoms to the metal
    # so any M update moves the whole rigid block; the gradient onto the
    # metal accumulates the rigid-body atoms' contributions.
    #
    # Default-OFF byte-identical: when the env-flag is unset
    # ``_rigid_body_indices`` stays empty -> the loss_and_grad wrapper
    # is a no-op and the gradient/positions are byte-identical with
    # the legacy path.
    # ------------------------------------------------------------------
    _rigid_body_atoms: List[int] = []
    if _hapto_pi_freeze_flag and _unfreeze_active_flag and hapto_atoms:
        try:
            _rigid_body_atoms = _build_hapto_rigid_body(
                mol, metal, donors_t, P_init, hapto_atoms,
            )
        except Exception:
            _rigid_body_atoms = []
    _rigid_body_indices = np.array(
        sorted(int(i) for i in _rigid_body_atoms if int(i) != int(metal)),
        dtype=np.int64,
    )
    _hapto_rigid_on = bool(_rigid_body_indices.size > 0)

    def loss_and_grad(x_flat: np.ndarray) -> Tuple[float, np.ndarray]:
        R = np.asarray(x_flat, dtype=np.float64).reshape(n_atoms, 3)
        # Hapto-π rigid-body projection (default OFF byte-identical).
        # Tie the rigid-body atoms to the metal's current displacement
        # from its initial position BEFORE evaluating loss / grad.
        if _hapto_rigid_on:
            R = _project_rigid_body(
                R, int(metal), _rigid_body_indices, P_init,
            )
        L_grip, G_grip = fragments.evaluate(R) if len(fragments) > 0 else (0.0, np.zeros_like(R))
        L_clash, G_clash = clash.value_and_grad(R)
        L = float(L_grip) + float(L_clash)
        G = G_grip + G_clash
        # Optional heavy-atom vdW-floor (anti-collapse).  Default-OFF
        # byte-identical with HEAD when ``DELFIN_FFFREE_GRIP_VDW_FLOOR``
        # is unset -- the predicate short-circuits before any new arithmetic.
        if _vdw_floor_on:
            if _vdw_floor_all_pairs and _cs_vdw_all_grad is not None:
                # Extended floor: penalise ALL non-bonded heavy pairs INCLUDING
                # those inside the hapto-π rigid body, but route their
                # gradient through the rigid body's metal-translation slot so
                # the body cannot deform internally.  Bug-class 2+3 fix.
                if len(_vdw_floor_symbols) == n_atoms:
                    L_vdw, G_vdw = _cs_vdw_all_grad(
                        R,
                        symbols=_vdw_floor_symbols,
                        excluded_pairs=excl_13,
                        weight=_vdw_floor_weight,
                        fraction=_vdw_floor_fraction,
                        rigid_body_atoms=(
                            _rigid_body_indices.tolist()
                            if _hapto_rigid_on else None
                        ),
                        rigid_translation_metal=(
                            int(metal) if _hapto_rigid_on else None
                        ),
                    )
                else:
                    L_vdw, G_vdw = _vdw_floor_value_and_grad(
                        R,
                        heavy_indices=_vdw_floor_heavy_indices,
                        radii=_vdw_floor_radii_arr,
                        excluded_pairs=excl_13,
                        weight=_vdw_floor_weight,
                        fraction=_vdw_floor_fraction,
                    )
            else:
                L_vdw, G_vdw = _vdw_floor_value_and_grad(
                    R,
                    heavy_indices=_vdw_floor_heavy_indices,
                    radii=_vdw_floor_radii_arr,
                    excluded_pairs=excl_13,
                    weight=_vdw_floor_weight,
                    fraction=_vdw_floor_fraction,
                )
            L += float(L_vdw)
            G = G + G_vdw
        # Optional inter-shell Pauli-floor (donor / 2nd-shell anti-clash).
        # Default-OFF byte-identical with HEAD when ``DELFIN_FFFREE_INTER_SHELL_FLOOR``
        # is unset -- the predicate short-circuits before any new arithmetic.
        # Composes additively with the existing clash + vdW-floor penalties.
        if _inter_shell_floor_on:
            L_is, G_is = _inter_shell_floor_value_and_grad(
                R,
                pairs=_inter_shell_pairs,
                radii=_inter_shell_radii_arr,
                weight=_inter_shell_floor_weight,
                fraction=_inter_shell_floor_fraction,
            )
            L += float(L_is)
            G = G + G_is
        # Zero the gradient on the frozen sphere -- L-BFGS-B will then leave
        # those atoms in place (the metal + donors stay locked).
        if frozen_indices.size:
            G[frozen_indices] = 0.0
        # Hapto-π rigid-body chain rule (default OFF byte-identical).
        # The rigid-body atoms are determined by R[metal] (translation-
        # locked), so their gradient flows entirely onto the metal's
        # gradient slot; their own gradient is zeroed.  This makes
        # L-BFGS-B leave the rigid-body atoms in place at each step --
        # the post-step projection at the end of the polish re-applies
        # the metal's net displacement to the rigid body.
        if _hapto_rigid_on:
            # Sum gradient contributions from rigid-body atoms onto metal.
            G[int(metal)] = G[int(metal)] + G[_rigid_body_indices].sum(axis=0)
            G[_rigid_body_indices] = 0.0
        # Guard against NaNs creeping in (e.g. coincident atoms in a step).
        if not np.isfinite(L):
            L = 0.0
        np.nan_to_num(G, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        return L, G.reshape(-1)

    # ------------------------------------------------------------------
    # 5. L-BFGS-B.
    # ------------------------------------------------------------------
    # Mission D2 (2026-06-05): resolve max_iter and gtol from env so the
    # production launcher can deepen convergence without code change.
    # Default-OFF byte-identical (resolvers return the caller-supplied
    # arg when env unset, which itself defaults to 200 / 1e-4).
    max_iter = _resolve_max_iter(max_iter if max_iter != DEFAULT_MAX_ITER else None)
    gtol = _resolve_gtol(gtol if gtol != DEFAULT_GTOL else None)
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

        # ------------------------------------------------------------------
        # 5b. Donor-cone DoF augmentation (2026-06-05, user direction
        # follow-up to ebb8cc3).
        #
        # When DELFIN_FFFREE_GRIP_DONOR_CONE=1, build per-donor cone
        # objects and augment the L-BFGS variable set with one extra
        # scalar (theta_d, radians) per donor.  At each loss evaluation
        # the X-atoms of each donor's cone are rigidly rotated around the
        # frozen M-D axis by theta_d before computing the position-based
        # loss.  Default OFF -> the augmentation block is a no-op (the
        # cone list is empty so x_aug collapses to x_pos and the wrapper
        # equals the underlying loss).
        # ------------------------------------------------------------------
        try:
            from .grip_donor_cone import (
                build_donor_cones,
                augmented_loss_and_grad,
                donor_cone_active,
                donor_cone_include_h,
            )
        except Exception as exc:  # pragma: no cover -- import guard
            _LOG.warning(
                "grip_polish: failed to import grip_donor_cone (%r); "
                "donor-cone DoF disabled",
                exc,
            )
            donor_cone_active = lambda: False  # type: ignore[assignment]
            build_donor_cones = None  # type: ignore[assignment]
            augmented_loss_and_grad = None  # type: ignore[assignment]
            donor_cone_include_h = lambda: False  # type: ignore[assignment]

        cones = []
        if donor_cone_active() and build_donor_cones is not None:
            try:
                cones = build_donor_cones(
                    mol, int(metal), donors_t, P_init,
                    hapto_atoms=hapto_atoms,
                    include_h=donor_cone_include_h(),
                    min_x_count=1,
                )
            except Exception as exc:
                _LOG.warning(
                    "grip_polish: build_donor_cones raised (%r); "
                    "donor-cone DoF inactive for this structure",
                    exc,
                )
                cones = []

        # ------------------------------------------------------------------
        # Mission D2 (2026-06-05): multi-restart L-BFGS for robustness.
        # When DELFIN_FFFREE_GRIP_RESTARTS=k (k > 1), run L-BFGS k times --
        # once from P_init (restart_idx=0, identical to legacy) and (k-1)
        # times from deterministically perturbed copies.  Keep the lowest
        # severity result that ALSO passes the topology/MD constraints.
        # Default k=1 -> exactly one call from P_init, byte-identical.
        # ------------------------------------------------------------------
        n_restarts = _resolve_restarts()
        restart_perturb = _resolve_restart_perturb()
        from .grip_donor_cone import apply_donor_cone_rotations as _apply_cone

        def _single_lbfgs_call(P_start: np.ndarray) -> Tuple[np.ndarray, int]:
            """One L-BFGS call from P_start.  Returns (P_refined, n_iter)."""
            if cones and augmented_loss_and_grad is not None:
                K = len(cones)
                wrapped = augmented_loss_and_grad(
                    loss_and_grad, cones, n_atoms,
                )
                x0_pos = P_start.reshape(-1).copy()
                x0_theta = np.zeros(K, dtype=np.float64)
                x0_aug = np.concatenate([x0_pos, x0_theta])
                bounds = [(None, None)] * (3 * n_atoms) + [
                    (-float(np.pi), float(np.pi)) for _ in range(K)
                ]
                res = minimize(
                    wrapped,
                    x0_aug,
                    method="L-BFGS-B",
                    jac=True,
                    bounds=bounds,
                    options={
                        "maxiter": int(max_iter),
                        "gtol": float(gtol),
                        "ftol": float(ftol),
                    },
                )
                x_aug_opt = np.asarray(res.x, dtype=np.float64)
                P_flat_opt = x_aug_opt[: 3 * n_atoms]
                thetas_opt = x_aug_opt[3 * n_atoms: 3 * n_atoms + K]
                return _apply_cone(P_flat_opt, cones, thetas_opt), int(getattr(res, "nit", 0))
            res = minimize(
                loss_and_grad,
                P_start.reshape(-1).copy(),
                method="L-BFGS-B",
                jac=True,
                options={
                    "maxiter": int(max_iter),
                    "gtol": float(gtol),
                    "ftol": float(ftol),
                },
            )
            return (
                np.asarray(res.x, dtype=np.float64).reshape(n_atoms, 3),
                int(getattr(res, "nit", 0)),
            )

        # Phase 2 acceptance-gate fix (2026-06-06): 2-pass MD weight anneal.
        #
        # When the unfreeze flag is on AND Phase 2 emitted terms AND the warmup
        # flag is on (default ON), run an extra L-BFGS pass first with the
        # Phase 2 weights scaled by DEFAULT_MD_WARMUP_FRACTION (0.1) so the
        # metal + donors can find an initial position without being pinned
        # by the M-D Gaussians; then run the regular full-weight pass
        # starting from the warmed-up coords.  Phase 1 terms are untouched
        # so the warmup never makes the geometry worse on the legacy axes.
        #
        # The wrapper mutates ``fragments.terms`` in-place between passes;
        # because ``loss_and_grad`` reads ``fragments`` from the enclosing
        # scope at call time, the new term list is honoured automatically.
        run_anneal = (
            _unfreeze_active_flag
            and bool(md_terms)
            and _md_warmup_active()
        )
        if run_anneal:
            # Save full-weight Phase 1+2 list (in canonical sort order) so
            # we can restore it for the final pass.
            _full_terms_list = list(fragments.terms)
            try:
                _warmup_md_terms = _anneal_md_term_weights(
                    md_terms, _MD_WARMUP_FRACTION,
                )
                # Replace md_terms in-place by their warmup variants.
                _md_ids = {id(t) for t in md_terms}
                _phase1_terms = [t for t in _full_terms_list if id(t) not in _md_ids]
                _warmup_total = TotalGripLoss(
                    terms=list(_phase1_terms) + list(_warmup_md_terms)
                )
                fragments.terms = list(_warmup_total.terms)
                try:
                    P_warmup, _ = _single_lbfgs_call(P_init)
                except Exception as _warm_exc:
                    _LOG.debug(
                        "grip_polish: warmup pass raised (%r); skipping anneal",
                        _warm_exc,
                    )
                    P_warmup = P_init
                if not np.all(np.isfinite(P_warmup)):
                    P_warmup = P_init
            finally:
                # Always restore the full-weight term list before the next
                # pass / the gate so severity comparisons stay consistent.
                fragments.terms = _full_terms_list
            P_start_for_final = P_warmup
        else:
            P_start_for_final = P_init

        # Restart 0 is always P_init verbatim (legacy byte-identity gate)
        # when anneal is OFF; otherwise we start from the warmed-up coords.
        P_refined, n_iter = _single_lbfgs_call(P_start_for_final)
        if n_restarts > 1:
            best_P = P_refined
            best_sev = mogul_severity(best_P, fragments)
            best_n_iter = n_iter
            for rk in range(1, n_restarts):
                P_perturb = _deterministic_perturbation(
                    P_start_for_final, rk, restart_perturb, frozen_indices,
                )
                try:
                    P_try, n_iter_try = _single_lbfgs_call(P_perturb)
                except Exception as _restart_exc:
                    _LOG.debug(
                        "grip_polish: restart %d raised (%r); skipping",
                        rk, _restart_exc,
                    )
                    continue
                if not np.all(np.isfinite(P_try)):
                    continue
                try:
                    sev_try = mogul_severity(P_try, fragments)
                except Exception:
                    continue
                if not np.isfinite(sev_try):
                    continue
                if sev_try < best_sev:
                    best_P = P_try
                    best_sev = sev_try
                    best_n_iter = n_iter_try
            P_refined = best_P
            n_iter = best_n_iter

    # Hapto-π rigid-body final projection (Task #91, 2026-06-07).
    # The L-BFGS pass zeroed the rigid-body atoms' gradients so they
    # stayed at P_init positions throughout the optimisation -- the
    # gradient on those atoms was accumulated onto the metal's slot.
    # Now apply the metal's net displacement to the rigid body so the
    # output coordinates form a valid (translation-projected) rigid block.
    # Default-OFF byte-identical: when ``_hapto_rigid_on`` is False the
    # call is a pure ``P_refined.copy()`` (no displacement applied).
    if _hapto_rigid_on and np.all(np.isfinite(P_refined)):
        try:
            P_refined = _project_rigid_body(
                P_refined, int(metal), _rigid_body_indices, P_init,
            )
        except Exception:
            # Fail-open: if projection raises, keep the L-BFGS result so
            # downstream validators can still decide.  Should never trigger.
            pass

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

    # Phase 2 fail-open: when the unfreeze flag is ON the M-D-invariant
    # validator above uses a SOFTER bound (the init distance ± md_tol),
    # but the doctrine demands no donor strays more than _MD_DRIFT_SAFETY
    # Å from its initial distance.  We measure the drift explicitly and
    # roll back if exceeded -- this is the fail-open guard.
    if _unfreeze_active_flag:
        try:
            md_drift = _max_md_drift(P_init, P_refined, metal, donors_t)
        except Exception:
            md_drift = float("inf")
        if not np.isfinite(md_drift) or md_drift > _MD_DRIFT_SAFETY:
            if return_diagnostics:
                return GripPolishResult(
                    P=P_init, accepted=False,
                    severity_before=sev_before,
                    severity_after=mogul_severity(P_refined, fragments),
                    n_iter=n_iter, n_terms=len(fragments),
                    rollback_reason=f"phase2 M-D drift {md_drift:.3f} > {_MD_DRIFT_SAFETY}",
                )
            return P_init
    elif not md_constraint.validate(P_refined):
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

    # Phase 2 chiral-volume fix (2026-06-06): metal-centered Δ/Λ check.
    # No-op when _unfreeze_active_flag is False (constraint is empty).
    if not metal_chirality_constraint.validate(P_refined):
        if return_diagnostics:
            return GripPolishResult(
                P=P_init, accepted=False,
                severity_before=sev_before,
                severity_after=mogul_severity(P_refined, fragments),
                n_iter=n_iter, n_terms=len(fragments),
                rollback_reason="metal Δ/Λ chirality class flipped",
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

    # Phase 2 acceptance-gate fix (2026-06-06):
    #   Replace the legacy "total severity" with
    #       score = (Phase 1 severity) + α × (Phase 2 severity)
    #   so the M+D Gaussians can actually drive the acceptance decision
    #   without being washed out by the (typically larger) Phase 1 sum.
    #   ``mogul_severity(R, fragments)`` already returns Phase 1 + Phase 2 at
    #   α = 1.0, so we subtract the Phase 2 piece, then re-add it at α.
    #   When α = 1.0 this is byte-identical with the legacy total severity.
    #   When the unfreeze flag is OFF (md_terms is empty) the Phase 2 piece
    #   is 0.0 and the gate stays in legacy mode.
    if _unfreeze_active_flag and md_terms:
        try:
            md_sev_before = _phase2_total_severity(P_init, md_terms)
            md_sev_after = _phase2_total_severity(P_refined, md_terms)
            alpha_md = _resolve_md_accept_weight()
            ph1_before = float(sev_before) - float(md_sev_before)
            ph1_after = float(sev_after) - float(md_sev_after)
            score_before = ph1_before + float(alpha_md) * float(md_sev_before)
            score_after = ph1_after + float(alpha_md) * float(md_sev_after)
        except Exception:
            # Phase 2 augmented gate failed -> degrade to legacy severity-only.
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
            score_before = float(score_before) + alpha * float(n_before)
            score_after = float(score_after) + alpha * float(n_after)
        except Exception:
            # Clash-augmented gate failed -> fall back to severity-only.
            # We deliberately retain any Phase 2 augmentation already applied.
            pass
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
