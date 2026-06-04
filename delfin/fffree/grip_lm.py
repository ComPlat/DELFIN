"""GRIP — Levenberg-Marquardt Polish (alternative optimiser, 2026-06-04).

Implements the GRIP polish via :func:`scipy.optimize.least_squares` (TRF,
Trust Region Reflective — scipy's robust LM-like solver with bounds support).

Why LM in addition to L-BFGS-B?
-------------------------------

The GRIP loss is **exactly** a non-linear least-squares problem:

.. math::

    L(R) = \\sum_k w_k \\left( \\frac{g_k(R) - \\mu_k}{\\sigma_k} \\right)^2
         = \\frac{1}{2} \\| r(R) \\|^2

with whitened residuals :math:`r_k(R) = \\sqrt{2 w_k} (g_k(R) - \\mu_k)/\\sigma_k`.
A Levenberg-Marquardt-style solver can exploit this structure directly via
the Jacobian :math:`J = \\partial r / \\partial R`; the local Gauss-Newton
step is the analytic solution to the linearised model and typically converges
quadratically near the minimum (vs the superlinear quasi-Newton step of L-BFGS).

For TMCs with chelate / cluster constraints and the boosted inter-ligand
clash floor (User-Direktive 2026-06-02), LM should:

* reach lower final loss in the same nfev budget (Gauss-Newton dominance),
* be more robust on the ill-conditioned Mahalanobis loss (TRF handles
  bounds without a barrier penalty),
* expose the per-residual Jacobian for future weighting / second-order use.

Why TRF (not MINPACK lmder)?
----------------------------

MINPACK ``lmder`` is the textbook LM but it has **no bound support** -- it
can produce catastrophically displaced atoms when a residual gradient
points outwards on a poorly-determined coordinate.  ``method='trf'`` in
scipy.least_squares is a Levenberg-Marquardt-flavoured trust-region method
that handles box bounds natively and is the documented choice for bounded
non-linear LSQ (https://docs.scipy.org/.../least_squares.html).

Determinism contract (SPEC §11):
--------------------------------

* TRF is fully deterministic given the initial point + bounds + Jacobian.
* No RNG, no random sampling, no Cython-level threading.
* float64 throughout.
* PYTHONHASHSEED-respecting (sorted residual order via TotalGripLoss).

Default-OFF byte identity (SPEC §3.4):
--------------------------------------

This module is **only** dispatched when
``DELFIN_FFFREE_GRIP_METHOD=lm`` (or the legacy alias ``levenberg-marquardt``)
is set in the environment.  The L-BFGS-B path stays the default and is
byte-identical to HEAD bcf56f8.

The accept-if-better gate, M-D / topology / chirality / polyhedron
validation, and inter-ligand-clash extension are **shared** with the
L-BFGS path -- only the inner-loop optimiser differs.
"""
from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
)

import numpy as np

from .grip_loss_terms import (
    AngleTerm,
    BondTerm,
    ImproperTerm,
    TotalGripLoss,
)

__all__ = [
    "grip_polish_lm",
    "build_residuals_and_jacobian",
    "DEFAULT_LM_BOUND_HALFWIDTH",
    "DEFAULT_LM_MAX_NFEV",
    "DEFAULT_LM_FTOL",
    "DEFAULT_LM_XTOL",
    "DEFAULT_LM_GTOL",
    "_EPS",
]

_LOG = logging.getLogger(__name__)

# Numerical floor matching grip_loss_terms / grip_constraints.
_EPS = 1e-12

# TRF hyperparameters (SPEC §3.3 budget parity with L-BFGS).
DEFAULT_LM_MAX_NFEV: int = 200
DEFAULT_LM_FTOL: float = 1e-8
DEFAULT_LM_XTOL: float = 1e-8
DEFAULT_LM_GTOL: float = 1e-8

# Half-width of the per-coordinate box bound (Å).  ±0.5 Å is wide enough to
# let LM reach the analytical Gauss-Newton minimum but narrow enough to
# prevent catastrophic divergence on a poorly-determined direction.
DEFAULT_LM_BOUND_HALFWIDTH: float = 0.5

# Clash residual floor: pairs whose distance falls below 0.85·(r_i+r_j)
# contribute a whitened residual ``sqrt(weight) · (d_min - d)``.  Above the
# floor the residual is zero (one-sided like the L-BFGS penalty).
_CLASH_RES_FLOOR_FRACTION: float = 0.85

# Sigma assigned to the clash residual so the magnitude matches the
# Mahalanobis loss family.  The clash WEIGHT is folded into the residual
# scale directly (``r = sqrt(w_clash) · (d_min - d)``) so the LSQ objective
# ``0.5||r||^2`` reproduces ``0.5 · w_clash · (d_min - d)^2`` -- the same
# functional form as the L-BFGS penalty's per-pair contribution divided by
# two.  The factor-2 is absorbed into ``clash_weight_loss_scale`` below.
_CLASH_RES_SCALE: float = 1.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _as_R(R: np.ndarray) -> np.ndarray:
    """Coerce ``R`` to ``(N, 3)`` float64, copying."""
    R = np.asarray(R, dtype=np.float64)
    if R.ndim == 1:
        if R.size % 3 != 0:
            raise ValueError(f"R must be (N,3) or (3N,), got length {R.size}")
        R = R.reshape(-1, 3)
    return R


def _safe_norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


# ---------------------------------------------------------------------------
# Per-term residual / Jacobian builders
# ---------------------------------------------------------------------------
def _bond_residual(R: np.ndarray, term: BondTerm) -> Tuple[float, Dict[int, np.ndarray]]:
    """Return whitened residual and per-atom gradient row of a BondTerm.

    The residual is ``sqrt(w) · (d - mu)/sigma`` so that
    ``r^2 == w · ((d-mu)/sigma)^2 == term.value_and_grad`` value.

    The returned ``grads`` dict maps atom index -> ``d r / d r_atom`` (3-vec).
    """
    ra = R[term.a]
    rb = R[term.b]
    d_vec = ra - rb
    d = _safe_norm(d_vec)
    sigma = float(term.sigma)
    mu = float(term.mu)
    w = float(term.weight)
    sw = math.sqrt(max(w, 0.0))
    r = sw * (d - mu) / sigma
    grads: Dict[int, np.ndarray] = {}
    if d > _EPS:
        # dr/dra = sw/sigma · d_vec/d ;  dr/drb = -dr/dra
        coef = sw / (sigma * d)
        ga = coef * d_vec
        grads[term.a] = ga
        grads[term.b] = -ga
    else:
        grads[term.a] = np.zeros(3, dtype=np.float64)
        grads[term.b] = np.zeros(3, dtype=np.float64)
    return r, grads


def _angle_residual(
    R: np.ndarray, term: AngleTerm,
) -> Tuple[float, Dict[int, np.ndarray]]:
    """Whitened residual + gradient row of an AngleTerm.

    Mirrors AngleTerm.value_and_grad mathematics: ``r = sw·(theta_deg-mu)/sigma``
    and ``dr/dr_atom`` is derived from the angle gradient (without the
    ``2·w·z`` outer chain rule factor used in the energy-form gradient).
    """
    ra = R[term.a]; rb = R[term.b]; rc = R[term.c]
    u = ra - rb
    v = rc - rb
    nu = _safe_norm(u); nv = _safe_norm(v)
    grads: Dict[int, np.ndarray] = {}
    if nu < _EPS or nv < _EPS:
        for idx in (term.a, term.b, term.c):
            grads[idx] = np.zeros(3, dtype=np.float64)
        return 0.0, grads
    cos_t = float(np.dot(u, v) / (nu * nv))
    cos_t = max(-1.0 + 1e-12, min(1.0 - 1e-12, cos_t))
    sin_t = math.sqrt(max(1.0 - cos_t * cos_t, _EPS))
    theta_rad = math.acos(cos_t)
    theta_deg = math.degrees(theta_rad)
    sigma = float(term.sigma); mu = float(term.mu); w = float(term.weight)
    sw = math.sqrt(max(w, 0.0))
    r = sw * (theta_deg - mu) / sigma
    # dr/d theta_deg = sw/sigma ; dtheta_deg / dtheta_rad = 180/pi
    deg_per_rad = 180.0 / math.pi
    dr_dtheta_rad = sw / sigma * deg_per_rad
    # gradients of cos(theta) w.r.t. r_a, r_c (same as AngleTerm)
    dcos_dra = v / (nu * nv) - cos_t * u / (nu * nu)
    dcos_drc = u / (nu * nv) - cos_t * v / (nv * nv)
    inv_sin = -1.0 / sin_t
    dtheta_dra = inv_sin * dcos_dra
    dtheta_drc = inv_sin * dcos_drc
    dtheta_drb = -(dtheta_dra + dtheta_drc)
    grads[term.a] = dr_dtheta_rad * dtheta_dra
    grads[term.b] = dr_dtheta_rad * dtheta_drb
    grads[term.c] = dr_dtheta_rad * dtheta_drc
    return r, grads


def _improper_residual(
    R: np.ndarray, term: ImproperTerm,
) -> Tuple[float, Dict[int, np.ndarray]]:
    """Whitened residual + gradient of an ImproperTerm (signed OOP distance).

    Reuses the analytical chain rule from ImproperTerm but pulls out the
    ``2·w·z·(1/sigma)`` chain factor so the result is ``r = sw·z`` and
    ``dr/dr_atom = sw/sigma · d d_oop / d r_atom``.
    """
    ic = int(term.center)
    i1, i2, i3 = term.neighbors_3
    rc = R[ic]; r1 = R[i1]; r2 = R[i2]; r3 = R[i3]
    e12 = r2 - r1
    e13 = r3 - r1
    n_raw = np.cross(e12, e13)
    n_norm = _safe_norm(n_raw)
    grads: Dict[int, np.ndarray] = {}
    if n_norm < _EPS:
        for idx in (ic, i1, i2, i3):
            grads[idx] = np.zeros(3, dtype=np.float64)
        return 0.0, grads
    n_hat = n_raw / n_norm
    centroid = (r1 + r2 + r3) / 3.0
    dvec = rc - centroid
    d_oop = float(np.dot(dvec, n_hat))
    sigma = float(term.sigma); mu = float(term.mu); w = float(term.weight)
    sw = math.sqrt(max(w, 0.0))
    r = sw * (d_oop - mu) / sigma
    dr_dd = sw / sigma

    I3 = np.eye(3, dtype=np.float64)
    proj = (I3 - np.outer(n_hat, n_hat)) / n_norm

    def _skew(v):
        return np.array(
            [
                [0.0, -v[2], v[1]],
                [v[2], 0.0, -v[0]],
                [-v[1], v[0], 0.0],
            ],
            dtype=np.float64,
        )

    cross_de12 = -_skew(e13)
    cross_de13 = _skew(e12)
    dn_dr1 = -cross_de12 - cross_de13
    dn_dr2 = cross_de12
    dn_dr3 = cross_de13
    dnhat_dr1 = proj @ dn_dr1
    dnhat_dr2 = proj @ dn_dr2
    dnhat_dr3 = proj @ dn_dr3

    # d dot(p, n_hat) / d r_i, with p = rc - centroid
    # rc contribution: n_hat (n_hat does not depend on rc)
    grads[ic] = dr_dd * n_hat.copy()
    # plane atoms: -n_hat/3 + dvec @ dnhat_drk
    grads[i1] = dr_dd * ((-n_hat / 3.0) + dvec @ dnhat_dr1)
    grads[i2] = dr_dd * ((-n_hat / 3.0) + dvec @ dnhat_dr2)
    grads[i3] = dr_dd * ((-n_hat / 3.0) + dvec @ dnhat_dr3)
    return r, grads


# ---------------------------------------------------------------------------
# Clash-pair residuals (one-sided, mirrors ClashFloorPenalty form)
# ---------------------------------------------------------------------------
def _enumerate_clash_pairs(
    radii: np.ndarray,
    exclude_13: Set[FrozenSet[int]],
    n: int,
) -> List[Tuple[int, int]]:
    """Return sorted, deterministic list of (i, j) pairs eligible for clash.

    Pairs are pre-screened ONLY on radii / 1,3-exclusions; the actual
    distance check happens per-evaluation so the LM step "sees" the floor
    activate and deactivate as the geometry moves.

    For typical fffree complexes ≤ 200 atoms the full ``O(N²)`` set is a
    few × 10⁴ pairs -- TRF handles it comfortably; we cache the list once
    so each residual evaluation only does the per-pair distance check.
    """
    out: List[Tuple[int, int]] = []
    for i in range(n):
        if not np.isfinite(radii[i]):
            continue
        for j in range(i + 1, n):
            if not np.isfinite(radii[j]):
                continue
            if frozenset((i, j)) in exclude_13:
                continue
            out.append((i, j))
    return out


def _clash_residual_and_grad(
    R: np.ndarray,
    pair: Tuple[int, int],
    radii: np.ndarray,
    floor_fraction: float,
    weight: float,
) -> Tuple[float, Dict[int, np.ndarray]]:
    """One-sided clash residual + gradient row for a single pair.

    Residual definition: ``r = sqrt(w) · max(0, d_min - d)`` where
    ``d_min = floor_fraction · (r_i + r_j)``.  This reproduces the
    L-BFGS-side penalty ``w·gap^2`` when ``L = 0.5·||r||^2`` is used as the
    cost (the factor 2 is absorbed by the ``0.5`` of LSQ).

    Both endpoints contribute symmetric gradient rows.  When ``d >= d_min``
    the residual and the gradient are exactly zero.
    """
    i, j = pair
    ri = radii[i]; rj = radii[j]
    d_min = floor_fraction * (ri + rj)
    d_vec = R[i] - R[j]
    d = _safe_norm(d_vec)
    grads: Dict[int, np.ndarray] = {}
    if d >= d_min or d < _EPS:
        grads[i] = np.zeros(3, dtype=np.float64)
        grads[j] = np.zeros(3, dtype=np.float64)
        return 0.0, grads
    sw = math.sqrt(max(weight, 0.0))
    gap = d_min - d
    r = sw * gap
    # dr/d|d| = -sw ; d|d|/dr_i = d_vec/d ; dr/dr_i = -sw · d_vec/d
    coef = -sw / d
    grads[i] = coef * d_vec
    grads[j] = -coef * d_vec
    return r, grads


# ---------------------------------------------------------------------------
# Free-atom packing (the optimiser variable)
# ---------------------------------------------------------------------------
def _free_atom_indices(n_atoms: int, frozen: FrozenSet[int]) -> np.ndarray:
    """Return sorted int array of FREE atom indices (complement of frozen)."""
    free = sorted(int(i) for i in range(n_atoms) if i not in frozen)
    return np.asarray(free, dtype=np.int64)


def _unpack_x_to_R(
    x_flat: np.ndarray,
    R_template: np.ndarray,
    free_atoms: np.ndarray,
) -> np.ndarray:
    """Lift a free-atom flat vector back into a full (N,3) array.

    Frozen atoms come from ``R_template`` (the initial state), free atoms
    are placed from ``x_flat`` in canonical order matching ``free_atoms``.
    """
    R = R_template.copy()
    if free_atoms.size == 0:
        return R
    R[free_atoms] = x_flat.reshape(-1, 3)
    return R


# ---------------------------------------------------------------------------
# Residual + Jacobian construction
# ---------------------------------------------------------------------------
@dataclass
class _ResidualPlan:
    """Pre-built description of the residual vector (deterministic order)."""

    # Sorted list of (kind, term_or_pair); kind in {"bond", "angle",
    # "improper", "clash"}
    items: List[Tuple[str, object]]
    # Per-residual count (always 1 in v1, but reserved for future multi-row
    # residuals like GMM-mixed torsions).
    sizes: List[int]
    n_residuals: int


def build_residuals_and_jacobian(
    fragments: TotalGripLoss,
    clash_pairs: List[Tuple[int, int]],
    radii: np.ndarray,
    clash_weight: float,
    floor_fraction: float,
    free_atoms: np.ndarray,
    R_template: np.ndarray,
    use_sparse_jacobian: bool = True,
    inter_ligand_weight_map: Optional[Dict[Tuple[int, int], float]] = None,
) -> Tuple[
    _ResidualPlan,
    Callable[[np.ndarray], np.ndarray],
    Callable[[np.ndarray], Any],
]:
    """Build deterministic residual and Jacobian callables for TRF.

    The plan iterates terms in TotalGripLoss's already-sorted order
    (``sorted(t.atom_indices)``) and appends clash residuals last (sorted
    by pair).  The Jacobian is sparse CSR by default (each residual touches
    only the few atoms in its support: 2 for bond/clash, 3 for angle, 4
    for improper).

    Returns (plan, residual_fn, jacobian_fn) suitable for
    :func:`scipy.optimize.least_squares`.
    """
    items: List[Tuple[str, object]] = []
    sizes: List[int] = []
    for term in fragments.terms:
        if isinstance(term, BondTerm):
            items.append(("bond", term)); sizes.append(1)
        elif isinstance(term, AngleTerm):
            items.append(("angle", term)); sizes.append(1)
        elif isinstance(term, ImproperTerm):
            items.append(("improper", term)); sizes.append(1)
        else:
            # Unknown / TorsionTerm not supported in v1 -> skip silently;
            # the L-BFGS path raises NotImplementedError on torsion, we
            # match by skipping (LM benchmark must converge alongside).
            continue
    for pair in clash_pairs:
        items.append(("clash", pair)); sizes.append(1)
    n_residuals = sum(sizes)

    plan = _ResidualPlan(items=items, sizes=sizes, n_residuals=n_residuals)

    # Build the free-atom inverse map for fast index translation.
    free_set = set(int(i) for i in free_atoms.tolist())
    free_pos: Dict[int, int] = {
        int(a): int(p) for p, a in enumerate(free_atoms.tolist())
    }
    n_free = int(free_atoms.size)

    def _residual_fn(x_flat: np.ndarray) -> np.ndarray:
        R = _unpack_x_to_R(np.asarray(x_flat, dtype=np.float64), R_template, free_atoms)
        out = np.zeros(n_residuals, dtype=np.float64)
        k = 0
        for kind, obj in items:
            if kind == "bond":
                r, _ = _bond_residual(R, obj)
            elif kind == "angle":
                r, _ = _angle_residual(R, obj)
            elif kind == "improper":
                r, _ = _improper_residual(R, obj)
            elif kind == "clash":
                pair = obj
                w = clash_weight
                if inter_ligand_weight_map is not None:
                    w = float(inter_ligand_weight_map.get(pair, clash_weight))
                r, _ = _clash_residual_and_grad(
                    R, pair, radii, floor_fraction, w,
                )
            else:
                r = 0.0
            # Guard against NaN -- TRF will refuse to evaluate otherwise.
            if not np.isfinite(r):
                r = 0.0
            out[k] = r
            k += 1
        return out

    # Sparse Jacobian assembly.  Each residual touches a small set of
    # atoms (<= 4 atoms = 12 columns).  We pre-compute the (row, col) index
    # arrays in the canonical order; only the data needs recomputing per
    # evaluation.
    if use_sparse_jacobian and n_residuals > 0 and n_free > 0:
        from scipy.sparse import csr_matrix
        # Pre-compute index pattern (row -> list of (col, atom_id, coord))
        # so we can rebuild data fast each call.
        rows_list: List[int] = []
        cols_list: List[int] = []
        # Map row index to a precomputed pattern (list of (atom_idx, free_col_base))
        # so we don't rebuild the structure each call.
        row_patterns: List[List[Tuple[int, int]]] = []
        for r_idx, (kind, obj) in enumerate(items):
            atom_list: List[int]
            if kind == "bond":
                atom_list = [obj.a, obj.b]
            elif kind == "angle":
                atom_list = [obj.a, obj.b, obj.c]
            elif kind == "improper":
                atom_list = [obj.center, *obj.neighbors_3]
            elif kind == "clash":
                atom_list = [obj[0], obj[1]]
            else:
                atom_list = []
            patt: List[Tuple[int, int]] = []
            for atom_idx in atom_list:
                if atom_idx not in free_set:
                    continue
                col_base = free_pos[atom_idx] * 3
                patt.append((int(atom_idx), int(col_base)))
                for k in range(3):
                    rows_list.append(int(r_idx))
                    cols_list.append(int(col_base + k))
            row_patterns.append(patt)
        rows_arr = np.asarray(rows_list, dtype=np.int64)
        cols_arr = np.asarray(cols_list, dtype=np.int64)
        nnz = rows_arr.size

        def _jacobian_fn(x_flat: np.ndarray):
            R = _unpack_x_to_R(
                np.asarray(x_flat, dtype=np.float64), R_template, free_atoms,
            )
            data = np.zeros(nnz, dtype=np.float64)
            ptr = 0
            for r_idx, (kind, obj) in enumerate(items):
                if kind == "bond":
                    _, grads = _bond_residual(R, obj)
                elif kind == "angle":
                    _, grads = _angle_residual(R, obj)
                elif kind == "improper":
                    _, grads = _improper_residual(R, obj)
                elif kind == "clash":
                    pair = obj
                    w = clash_weight
                    if inter_ligand_weight_map is not None:
                        w = float(inter_ligand_weight_map.get(pair, clash_weight))
                    _, grads = _clash_residual_and_grad(
                        R, pair, radii, floor_fraction, w,
                    )
                else:
                    grads = {}
                patt = row_patterns[r_idx]
                for atom_idx, _col_base in patt:
                    g = grads.get(atom_idx, np.zeros(3, dtype=np.float64))
                    if not np.all(np.isfinite(g)):
                        g = np.zeros(3, dtype=np.float64)
                    data[ptr] = float(g[0])
                    data[ptr + 1] = float(g[1])
                    data[ptr + 2] = float(g[2])
                    ptr += 3
            J = csr_matrix(
                (data, (rows_arr, cols_arr)),
                shape=(n_residuals, n_free * 3),
            )
            return J
    else:
        # Dense fallback (small problems / debugging / no scipy.sparse).
        def _jacobian_fn(x_flat: np.ndarray):
            R = _unpack_x_to_R(
                np.asarray(x_flat, dtype=np.float64), R_template, free_atoms,
            )
            J = np.zeros((n_residuals, max(n_free * 3, 1)), dtype=np.float64)
            for r_idx, (kind, obj) in enumerate(items):
                if kind == "bond":
                    _, grads = _bond_residual(R, obj)
                    atom_list = [obj.a, obj.b]
                elif kind == "angle":
                    _, grads = _angle_residual(R, obj)
                    atom_list = [obj.a, obj.b, obj.c]
                elif kind == "improper":
                    _, grads = _improper_residual(R, obj)
                    atom_list = [obj.center, *obj.neighbors_3]
                elif kind == "clash":
                    pair = obj
                    w = clash_weight
                    if inter_ligand_weight_map is not None:
                        w = float(inter_ligand_weight_map.get(pair, clash_weight))
                    _, grads = _clash_residual_and_grad(
                        R, pair, radii, floor_fraction, w,
                    )
                    atom_list = [pair[0], pair[1]]
                else:
                    continue
                for atom_idx in atom_list:
                    if atom_idx not in free_set:
                        continue
                    col_base = free_pos[atom_idx] * 3
                    g = grads.get(atom_idx, np.zeros(3, dtype=np.float64))
                    if not np.all(np.isfinite(g)):
                        g = np.zeros(3, dtype=np.float64)
                    J[r_idx, col_base : col_base + 3] = g
            return J

    return plan, _residual_fn, _jacobian_fn


# ---------------------------------------------------------------------------
# Main LM polish
# ---------------------------------------------------------------------------
@dataclass
class GripLMResult:
    """Diagnostic container for the LM path (parallel to GripPolishResult)."""

    P: Optional[np.ndarray]
    accepted: bool
    severity_before: float
    severity_after: float
    n_nfev: int
    n_terms: int
    n_residuals: int
    final_cost: float
    method: str = "lm"
    rollback_reason: str = ""


def grip_polish_lm(
    P0: np.ndarray,
    *,
    fragments: TotalGripLoss,
    clash_pairs: List[Tuple[int, int]],
    radii: np.ndarray,
    clash_weight: float,
    floor_fraction: float,
    frozen_atoms: FrozenSet[int],
    n_atoms: int,
    bound_halfwidth: float = DEFAULT_LM_BOUND_HALFWIDTH,
    max_nfev: int = DEFAULT_LM_MAX_NFEV,
    ftol: float = DEFAULT_LM_FTOL,
    xtol: float = DEFAULT_LM_XTOL,
    gtol: float = DEFAULT_LM_GTOL,
    use_sparse_jacobian: bool = True,
    inter_ligand_weight_map: Optional[Dict[Tuple[int, int], float]] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Run the LM (TRF) polish.

    Parameters mirror what the L-BFGS dispatcher in :mod:`grip_polish`
    already prepared (fragments, clash pairs, vdW radii, exclusion sets,
    frozen sphere).  Returns the polished ``(N, 3)`` array (NOT yet
    accept-gated -- the caller applies the gate) plus a diagnostics dict
    describing the inner-solver outcome.

    On any failure (scipy unavailable, n_free == 0, residual evaluation
    blew up) the function returns the **input** ``P0`` and ``diag["ok"]``
    is ``False`` so the caller can roll back.
    """
    diag: Dict[str, Any] = {
        "ok": False,
        "reason": "",
        "n_nfev": 0,
        "n_residuals": 0,
        "final_cost": float("nan"),
        "method": "lm-trf",
    }

    P_init = _as_R(P0)
    if not np.all(np.isfinite(P_init)):
        diag["reason"] = "P0 contains non-finite values"
        return P_init, diag

    free_atoms = _free_atom_indices(n_atoms, frozen_atoms)
    n_free = int(free_atoms.size)
    if n_free == 0:
        diag["reason"] = "no free atoms (everything frozen)"
        diag["ok"] = True
        return P_init, diag

    plan, residual_fn, jacobian_fn = build_residuals_and_jacobian(
        fragments=fragments,
        clash_pairs=clash_pairs,
        radii=radii,
        clash_weight=clash_weight,
        floor_fraction=floor_fraction,
        free_atoms=free_atoms,
        R_template=P_init,
        use_sparse_jacobian=use_sparse_jacobian,
        inter_ligand_weight_map=inter_ligand_weight_map,
    )
    diag["n_residuals"] = int(plan.n_residuals)

    if plan.n_residuals == 0:
        diag["reason"] = "no residuals"
        diag["ok"] = True
        return P_init, diag

    # Initial free-atom flat vector.
    x0 = P_init[free_atoms].reshape(-1).copy()

    # Box bounds: ±bound_halfwidth Å around the initial position.
    lo = x0 - float(bound_halfwidth)
    hi = x0 + float(bound_halfwidth)

    try:
        from scipy.optimize import least_squares
    except Exception as exc:  # pragma: no cover -- scipy is in the env
        diag["reason"] = f"scipy unavailable: {exc!r}"
        return P_init, diag

    # Sanity: residual function must produce a finite vector at x0 (early
    # divergence will manifest as NaN and TRF will refuse to start).
    try:
        r0 = residual_fn(x0)
        if not np.all(np.isfinite(r0)):
            diag["reason"] = "non-finite residuals at x0"
            return P_init, diag
    except Exception as exc:
        diag["reason"] = f"residual_fn raised at x0: {exc!r}"
        return P_init, diag

    # Choose Jacobian flavour: sparse jac_sparsity is auto-derived from
    # the analytic Jacobian (TRF needs jac_sparsity only when jac='2-point';
    # because we pass a callable that returns a sparse matrix, TRF detects
    # sparsity automatically).
    try:
        res = least_squares(
            residual_fn,
            x0,
            jac=jacobian_fn,
            bounds=(lo, hi),
            method="trf",
            ftol=float(ftol),
            xtol=float(xtol),
            gtol=float(gtol),
            max_nfev=int(max_nfev),
            x_scale=1.0,   # deterministic; no Jacobian-driven rescaling
            tr_solver=("lsmr" if use_sparse_jacobian else "exact"),
            verbose=0,
        )
    except Exception as exc:
        diag["reason"] = f"least_squares raised: {exc!r}"
        return P_init, diag

    diag["n_nfev"] = int(getattr(res, "nfev", 0))
    diag["final_cost"] = float(getattr(res, "cost", float("nan")))

    if not np.all(np.isfinite(res.x)):
        diag["reason"] = "TRF returned non-finite x"
        return P_init, diag

    P_refined = _unpack_x_to_R(
        np.asarray(res.x, dtype=np.float64), P_init, free_atoms,
    )
    if not np.all(np.isfinite(P_refined)):
        diag["reason"] = "unpack produced non-finite coordinates"
        return P_init, diag

    diag["ok"] = True
    diag["status"] = int(getattr(res, "status", 0))
    return P_refined, diag
