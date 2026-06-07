"""Donor-angle polyhedron floor (2026-06-07).

Auto-Diagnostic finding (V3 voll-pool, 2026-06-07): 890 / 6627 files
(13.4 %) showed the ``metal_axis_linearisation`` anomaly -- three or four
donors around a Pd(SP-4) / Pt(SP-4) centre drifting toward 180 deg, so the
coordination geometry collapses from square-planar (90° / 180°) to a near-
linear arrangement.

Root cause: the Mahalanobis severity loss
(:mod:`delfin.fffree.grip_loss_terms`) plus the clash floor
(:mod:`delfin.fffree.grip_constraints.ClashFloorPenalty`) only see
**bonded** geometric features.  Donor-M-donor angles are not bonded
fragments (M and the two donors live in the **frozen** sphere by default),
so they are NEVER penalised by the polish even when they drift far from
the ideal polyhedron angles.

This module adds a **universal**, **geometry-only** donor-angle floor:
for each donor pair ``(d_i, d_j)`` we compute the actual angle
``theta = angle(d_i, M, d_j)`` and the nearest expected angle for the
declared polyhedron ``theta_nearest``.  A smooth quadratic penalty fires
only when ``|theta - theta_nearest|`` exceeds the threshold (default
15 deg):

    L = w * (|delta| - tau)**2   if |delta| > tau   else 0

where ``delta = theta - theta_nearest`` and ``tau`` is the threshold.  At
the boundary the value and the gradient are both zero, so the term is C^1
continuous.  Analytical gradient -- no autodiff dependency.

The expected angles come from the same canonical vertex vectors that
:func:`delfin.fffree.polyhedra.ref_vectors` returns -- the table is built
by enumerating every pair of vertices and collecting their dot-product
angles.  No per-polyhedron / per-class branch is needed; new geometries
appear automatically as soon as they appear in :mod:`polyhedra`.

Env-flag (default OFF, byte-identical with HEAD c1e0fde when unset):

    DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR=1   -> term ON
    DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_WEIGHT=5.0
    DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG=15.0

The term is wired into :func:`delfin.fffree.grip_polish.grip_polish` --
see the ``loss_and_grad`` closure for the call site.  It is purely
geometric: no SMILES patterns, no per-class branches, no element-aware
heuristics.  The full universality contract (CN2 -> CN10) is covered by
:mod:`tests.test_donor_angle_polyhedron`.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

_LOG = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Env-flag / resolver scaffolding.
# ---------------------------------------------------------------------------
_FLOOR_ENV: str = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR"
_WEIGHT_ENV: str = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_WEIGHT"
_THRESHOLD_ENV: str = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG"

DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT: float = 5.0
DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG: float = 15.0

# Numerical guard rails -- match the clash/vdW floor module conventions so
# the polish does not emit NaNs when two donors land exactly on top of the
# metal (an edge case that should never happen, but defence in depth).
_NORM_EPS: float = 1e-9
_ANGLE_EPS: float = 1e-9


def donor_angle_polyhedron_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR`` is on.

    Default OFF (env unset / empty / "0" / "false" / "no" / "off").  Any
    other value (including "1", "true", "yes", "on", case-insensitive)
    activates the term.  Mirrors the convention used by the heavy-atom
    vdW-floor module in :mod:`grip_polish`.
    """
    raw = os.environ.get(_FLOOR_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _resolve_donor_angle_polyhedron_weight() -> float:
    """Effective penalty weight.  env > 5.0 (byte-identical default-OFF).

    Non-finite / non-numeric values fall through to the default with a
    single warning so the polish never crashes on a misconfigured flag.
    """
    raw = os.environ.get(_WEIGHT_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v >= 0.0:
                return v
            _LOG.warning(
                "donor_angle_polyhedron: %s=%r not a non-negative finite float; "
                "using default %.3f",
                _WEIGHT_ENV, raw, DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "donor_angle_polyhedron: %s=%r not numeric; using default %.3f",
                _WEIGHT_ENV, raw, DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            )
    return DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT


def _resolve_donor_angle_polyhedron_threshold_deg() -> float:
    """Effective deviation threshold (degrees).  env > 15.0.

    Pathological values (non-numeric, negative, infinite) fall back to the
    default.  Values >= 90 deg would silence the penalty for SP-4
    cis/trans confusion (the only mode this module exists to fix), so we
    do not clamp -- the user must opt in.
    """
    raw = os.environ.get(_THRESHOLD_ENV, "").strip()
    if raw:
        try:
            v = float(raw)
            if np.isfinite(v) and v >= 0.0:
                return v
            _LOG.warning(
                "donor_angle_polyhedron: %s=%r not a non-negative finite float; "
                "using default %.3f",
                _THRESHOLD_ENV, raw, DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
            )
        except (TypeError, ValueError):
            _LOG.warning(
                "donor_angle_polyhedron: %s=%r not numeric; using default %.3f",
                _THRESHOLD_ENV, raw, DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
            )
    return DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG


# ---------------------------------------------------------------------------
# Expected angle table -- universal, derived from ref_vectors.
# ---------------------------------------------------------------------------
def expected_angles_for_geometry(geometry: str) -> List[float]:
    """Return the sorted list of unique pairwise donor angles (degrees)
    encoded by the polyhedron ``geometry``.

    Mechanism: pull the unit-vector ideal vertex set from
    :func:`delfin.fffree.polyhedra.ref_vectors`; for every ordered pair
    (i < j) compute ``acos(v_i . v_j)``; cluster duplicates within 0.5 deg
    and return the sorted set.

    Returns ``[]`` for unknown geometries -- the caller short-circuits the
    penalty in that case (rather than picking an arbitrary template).

    Universality contract: NO per-class branches; the function works for
    every name :mod:`polyhedra` defines (CN2 .. CN10) plus any future
    addition (sandwich / f-block / Wells-canonical CN10).  Determinism:
    same geometry name -> same list, every call (no RNG, no env-read).
    """
    if not geometry:
        return []
    try:
        from delfin.fffree.polyhedra import ref_vectors
        V = np.asarray(ref_vectors(geometry), dtype=np.float64)
    except Exception:
        return []
    n = int(V.shape[0]) if V.ndim == 2 else 0
    if n < 2:
        return []
    # Unit-normalise (defensive: the table is already unit-row, but a
    # caller-supplied geometry may not be).
    norms = np.linalg.norm(V, axis=1)
    norms = np.where(norms < _NORM_EPS, 1.0, norms)
    U = V / norms[:, None]
    raw: List[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            c = float(np.dot(U[i], U[j]))
            # clamp into [-1, 1] for numerical safety
            c = max(-1.0, min(1.0, c))
            raw.append(math.degrees(math.acos(c)))
    # Deduplicate within 0.5 deg.  Sort the raw angles, then walk; any
    # angle within 0.5 deg of the last accepted one is dropped.  This
    # keeps the SP-4 set as {90, 180}, the T-4 set as {109.471}, etc.
    raw.sort()
    out: List[float] = []
    for v in raw:
        if not out or abs(v - out[-1]) > 0.5:
            out.append(v)
    return out


# Memoise the expected-angle list per geometry string -- the polish calls
# loss_and_grad O(max_iter) times per build, and the table never changes.
_EXPECTED_CACHE: Dict[str, Tuple[float, ...]] = {}


def get_expected_angles(geometry: str) -> Tuple[float, ...]:
    """Cached read-only version of :func:`expected_angles_for_geometry`."""
    if geometry in _EXPECTED_CACHE:
        return _EXPECTED_CACHE[geometry]
    angles = tuple(expected_angles_for_geometry(geometry))
    _EXPECTED_CACHE[geometry] = angles
    return angles


# ---------------------------------------------------------------------------
# Value + analytical gradient.
# ---------------------------------------------------------------------------
def _nearest_expected(theta_deg: float,
                      expected: Sequence[float]) -> Tuple[float, float]:
    """Return ``(theta_nearest_deg, signed_delta_deg)`` against the
    expected-angle list.  Tie-break: lowest expected angle wins (sorted
    list, first-min stable).
    """
    best_diff = float("inf")
    best_ref = expected[0]
    for ref in expected:
        d = abs(theta_deg - ref)
        if d < best_diff:
            best_diff = d
            best_ref = ref
    return best_ref, theta_deg - best_ref


def donor_angle_polyhedron_value_and_grad(
    R: np.ndarray,
    *,
    metal: int,
    donors: Sequence[int],
    expected_angles_deg: Sequence[float],
    weight: float,
    threshold_deg: float,
) -> Tuple[float, np.ndarray]:
    """Penalise donor-M-donor angles that drift more than ``threshold_deg``
    from the nearest expected polyhedron angle.

    Parameters
    ----------
    R : (N, 3) ndarray
        Current coordinates.
    metal : int
        Atom index of the metal centre.
    donors : sequence of int
        Atom indices of the donor atoms; pairs are formed over (i < j).
    expected_angles_deg : sequence of float
        Sorted, deduplicated list of allowed donor-M-donor angles for the
        declared polyhedron (degrees).  Empty -> no-op.
    weight : float
        Penalty multiplier.  ``<= 0`` -> no-op.
    threshold_deg : float
        Tolerance below which the penalty is zero (degrees).  ``< 0`` ->
        clamped to 0.

    Returns
    -------
    (L, G) : (float, (N, 3) ndarray)
        Loss value and gradient.  Pairs whose deviation is within
        ``threshold_deg`` of an expected angle contribute zero to both.

    Math (per donor pair (i, j), with d_i = R[i] - R[metal],
    d_j = R[j] - R[metal]):

        c = (d_i . d_j) / (|d_i| |d_j|)              # cosine of angle
        theta_rad = acos(clamp(c, -1, 1))
        theta_deg = degrees(theta_rad)
        theta_nearest, delta = nearest_expected(theta_deg, expected)
        gap_deg = max(0, |delta| - threshold_deg)
        L_pair = weight * gap_deg ** 2
        # gradient
        if |delta| > threshold_deg:
            dL/d_theta_deg = 2 weight * gap_deg * sign(delta)
            dL/d_theta_rad = dL/d_theta_deg * (180 / pi)
            dL/dc = -dL/d_theta_rad / sqrt(max(1 - c^2, eps))
            dL/d_d_i = (1/(|d_i||d_j|)) * (d_j - c * d_i / |d_i|^2 * (|d_i||d_j|))
                     = d_j/(|d_i||d_j|) - c * d_i / |d_i|^2
            dL/d_d_j = symmetric
            grad_R[i]      += dL/dc * dL/d_d_i
            grad_R[j]      += dL/dc * dL/d_d_j
            grad_R[metal]  -= dL/dc * (dL/d_d_i + dL/d_d_j)
    """
    n = int(R.shape[0])
    grad = np.zeros_like(R)
    if (
        len(donors) < 2
        or weight <= 0.0
        or not expected_angles_deg
        or not (0 <= int(metal) < n)
    ):
        return 0.0, grad

    tau = max(0.0, float(threshold_deg))
    expected = tuple(float(a) for a in expected_angles_deg)
    M = int(metal)
    donor_list = sorted(int(d) for d in donors if 0 <= int(d) < n and int(d) != M)
    if len(donor_list) < 2:
        return 0.0, grad

    total = 0.0
    Rm = R[M]
    rad_per_deg = math.pi / 180.0
    deg_per_rad = 180.0 / math.pi
    for ii in range(len(donor_list)):
        i = donor_list[ii]
        d_i = R[i] - Rm
        norm_i = float(np.linalg.norm(d_i))
        if norm_i < _NORM_EPS:
            continue
        inv_ni = 1.0 / norm_i
        inv_ni2 = inv_ni * inv_ni
        for jj in range(ii + 1, len(donor_list)):
            j = donor_list[jj]
            d_j = R[j] - Rm
            norm_j = float(np.linalg.norm(d_j))
            if norm_j < _NORM_EPS:
                continue
            inv_nj = 1.0 / norm_j
            inv_nj2 = inv_nj * inv_nj
            denom = norm_i * norm_j
            c = float(np.dot(d_i, d_j)) / denom
            # Clamp to the valid arccos domain.
            if c > 1.0:
                c = 1.0
            elif c < -1.0:
                c = -1.0
            theta_rad = math.acos(c)
            theta_deg = theta_rad * deg_per_rad
            _, delta_deg = _nearest_expected(theta_deg, expected)
            abs_delta = abs(delta_deg)
            if abs_delta <= tau:
                continue
            gap_deg = abs_delta - tau
            total += float(weight) * gap_deg * gap_deg
            # ---- analytical gradient ----
            # dL/d|delta|       = 2 w (|delta| - tau)
            # d|delta|/d_delta  = sign(delta)
            # so dL/d_theta_deg = 2 w gap_deg * sign(delta)
            dL_dtheta_deg = 2.0 * float(weight) * gap_deg * (
                1.0 if delta_deg >= 0.0 else -1.0
            )
            dL_dtheta_rad = dL_dtheta_deg * deg_per_rad  # dtheta_deg/dtheta_rad = 180/pi
            # dtheta/dc = -1 / sqrt(1 - c^2); guard against the saddle at
            # c = +-1 (perfectly linear / coincident donor pair).
            sin_theta_sq = max(1.0 - c * c, _ANGLE_EPS)
            dtheta_dc = -1.0 / math.sqrt(sin_theta_sq)
            dL_dc = dL_dtheta_rad * dtheta_dc
            # Chain rule for c = (d_i . d_j) / (|d_i| |d_j|):
            #   dc/d(d_i) = d_j / (|d_i||d_j|) - c * d_i / |d_i|^2
            #   dc/d(d_j) = d_i / (|d_i||d_j|) - c * d_j / |d_j|^2
            # Atoms i and j live in R-space; the metal sits at R[M]; the
            # diff vectors are d_i = R[i] - R[M] etc., so
            #   d(d_i)/dR[i]  = +I, d(d_i)/dR[M] = -I.  Same for j.
            inv_denom = 1.0 / denom
            dc_d_di = d_j * inv_denom - c * d_i * inv_ni2
            dc_d_dj = d_i * inv_denom - c * d_j * inv_nj2
            gi = dL_dc * dc_d_di
            gj = dL_dc * dc_d_dj
            grad[i] += gi
            grad[j] += gj
            grad[M] -= (gi + gj)

    if not np.isfinite(total):
        # Defensive: a divide-by-zero downstream should never escape, but
        # if it does we silently no-op so the polish does not crash.
        return 0.0, np.zeros_like(R)
    return float(total), grad


__all__ = [
    "DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT",
    "DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG",
    "donor_angle_polyhedron_active",
    "_resolve_donor_angle_polyhedron_weight",
    "_resolve_donor_angle_polyhedron_threshold_deg",
    "expected_angles_for_geometry",
    "get_expected_angles",
    "donor_angle_polyhedron_value_and_grad",
]
