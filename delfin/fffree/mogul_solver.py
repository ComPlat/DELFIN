"""Mogul-DG Phase B — Projected L-BFGS distance-geometry solver.

This module provides :func:`solve_dg`, the projected L-BFGS optimiser that
takes an initial point cloud :math:`P_{\\mathrm{init}} \\in \\mathbb{R}^{n \\times 3}`,
a pair of hard bounds matrices :math:`(L, U)` from
:mod:`delfin.fffree.mogul_bounds`, and a caller-supplied empirical severity
function :math:`L_{\\mathrm{emp}}(P)`, then minimises the augmented loss

.. math::

    L_{\\mathrm{total}}(P) = L_{\\mathrm{emp}}(P) + w_{\\mathrm{bnd}}
        \\sum_{(i,j)} \\phi(d_{ij}; l_{ij}, u_{ij})

where :math:`\\phi` is a one-sided quadratic that is zero inside
:math:`[l_{ij}, u_{ij}]` and grows quadratically outside.  This is the
classical Bertsekas projected-gradient formulation: in the limit of
:math:`w_{\\mathrm{bnd}} \\to \\infty` the gradient component normal to the
feasibility boundary is exactly the projection of the unconstrained
gradient onto the feasible set.  In practice a finite but large weight
gives the same answer at a fraction of the numerical noise.

Algorithm overview (SPEC §2.3, §2.4, §6, §7)
--------------------------------------------

1. **Augmented loss**.  Wrap the caller's severity_fn into a closure that
   adds the bounds-violation penalty.  Gradient is analytical: each pair
   :math:`(i, j)` outside its window contributes a force along the pair
   axis that pushes the two atoms back into the window.

2. **Multi-restart**.  Run L-BFGS-B from :math:`P_{\\mathrm{init}}` once
   (restart 0), then up to ``n_restarts - 1`` additional times from
   deterministically perturbed copies of :math:`P_{\\mathrm{init}}`.
   Frozen atoms (metal + donors, optional) are not perturbed.  The
   restart with the lowest augmented loss + zero bounds violations wins.

3. **M-D drift guard**.  After each restart, the maximum change of any
   metal-donor distance between the initial and refined geometry is
   checked.  Restarts that drift any M-D distance by more than
   ``md_drift_tol`` (default 0.5 Å) are discarded.

4. **Feasibility fallback** (SPEC §2.4).  If *no* restart converges with
   zero bounds violations, the bounds matrix is relaxed hierarchically:

   * relaxation tier 1: multiply 1,3-angle bound widths by ``1.5``
   * relaxation tier 2: multiply 1,4-torsion / non-bonded vdW bound
     widths by ``2.0`` (lower-bound only)
   * relaxation tier 3: shrink vdW floor (multiply lower-bound by
     ``0.7`` → ``0.8``)
   * M-D and D-D bounds **stay hard** at every tier.

   The classification of a pair as "1,3-angle", "1,4-torsion", "vdW",
   "M-D", "D-D" is provided by the caller through the optional
   ``pair_classes`` mapping; when omitted, the relaxer treats every
   non-M-D pair as a generic non-bonded pair (relaxed at tier 2/3).

5. **Determinism**.  All randomness is seeded from ``seed +
   restart_idx``; the L-BFGS solver itself is deterministic for a fixed
   input order and identical numpy build.  Two runs at
   ``PYTHONHASHSEED=0`` with the same arguments must produce
   byte-identical output (enforced by the test suite).

Public API
----------
* :func:`solve_dg` — the public solver entry point.
* :class:`SolverInfo` — diagnostic dataclass (mirrors the ``info`` dict).
* :func:`bounds_violation_penalty` — augmented-loss helper (also useful
  for unit tests of the projection mechanic).
* :func:`compute_md_drift` — M-D drift guard helper.

This module is the FF-free, universal solver for whole-complex
distance-geometry.  No SMILES strings, no element-specific constants,
no force-field calls — everything operates on bounds + a caller-supplied
severity callable.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.optimize import minimize as _scipy_minimize  # type: ignore
    _HAS_SCIPY = True
except Exception:  # pragma: no cover - scipy is required everywhere DELFIN runs
    _HAS_SCIPY = False


__all__ = [
    "solve_dg",
    "SolverInfo",
    "bounds_violation_penalty",
    "compute_md_drift",
    "DEFAULT_BOUNDS_WEIGHT",
    "DEFAULT_MAX_ITER",
    "DEFAULT_TOL",
    "DEFAULT_N_RESTARTS",
    "DEFAULT_SEED",
    "DEFAULT_MD_DRIFT_TOL",
    "DEFAULT_PERTURB_AMPLITUDE",
]


# ---------------------------------------------------------------------------
# Defaults & env knobs.  All defaults preserve byte-identity with the SPEC
# Phase B contract (Section 4 Phase B + Section 7).
# ---------------------------------------------------------------------------
DEFAULT_BOUNDS_WEIGHT: float = 50.0     # quadratic-penalty weight for hard bounds
DEFAULT_MAX_ITER: int = 1000             # L-BFGS iteration cap (SPEC §3.3)
DEFAULT_TOL: float = 1e-5                # |grad| convergence threshold (SPEC §2.3)
DEFAULT_N_RESTARTS: int = 3              # SPEC §2.3 multi-restart count
DEFAULT_SEED: int = 42                   # SPEC §7 — deterministic restart seed
DEFAULT_MD_DRIFT_TOL: float = 0.5        # Å, SPEC §6 — M-D drift guard
DEFAULT_PERTURB_AMPLITUDE: float = 0.10  # Å, perturbation amplitude per restart

# Upper-bound sentinel used by mogul_bounds (`_UB_INF`).  Pairs with
# ``upper >= 0.5 * _UB_INF`` are treated as effectively unbounded above and
# do not contribute an *upper*-bound penalty (only the lower-bound side
# matters, e.g. for non-bonded vdW floors and counter-ion separation).
_UB_INF_THRESHOLD: float = 5e5


# ---------------------------------------------------------------------------
# Diagnostic info container
# ---------------------------------------------------------------------------
@dataclass
class SolverInfo:
    """Container mirroring the spec's ``info`` dict."""

    n_iter: int = 0
    final_loss: float = float("inf")
    final_emp_loss: float = float("inf")
    final_penalty: float = float("inf")
    max_md_drift: float = 0.0
    n_bounds_violated: int = 0
    bounds_relaxed: List[Tuple[int, int, str]] = field(default_factory=list)
    restart_used: int = -1
    converged: bool = False
    failed: bool = False
    reason: str = ""
    n_restarts_tried: int = 0
    relaxation_tier: int = 0
    grad_norm: float = float("inf")

    def as_dict(self) -> Dict[str, Any]:
        d = dict(self.__dict__)
        d["bounds_relaxed"] = list(self.bounds_relaxed)
        return d


# ---------------------------------------------------------------------------
# Bounds-violation penalty (the projection mechanic, implemented as a smooth
# augmented loss).  This is reusable in tests + external callers that want
# to inspect feasibility without running the full solver.
# ---------------------------------------------------------------------------
def bounds_violation_penalty(
    P: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    weight: float = DEFAULT_BOUNDS_WEIGHT,
    eps: float = 1e-9,
) -> Tuple[float, np.ndarray]:
    """Compute the bounds-violation penalty + analytical gradient.

    The penalty is the sum over all ordered pairs ``i < j`` of:

    .. math::

        \\phi(d_{ij}) = \\tfrac{1}{2} w \\big(\\max(0, l_{ij} - d_{ij})^2
                                          + \\max(0, d_{ij} - u_{ij})^2\\big)

    Pairs with ``upper >= _UB_INF_THRESHOLD`` skip the upper-bound term
    (effectively unbounded above — counter-ions, distant non-bonded
    contacts).  The gradient is analytical and propagates symmetrically
    to both atoms in each violating pair.

    Parameters
    ----------
    P : (n, 3) ndarray
        Current point cloud.
    lower, upper : (n, n) ndarray
        Symmetric hard-bound matrices from
        :func:`delfin.fffree.mogul_bounds.build_bounds_matrix`.
    weight : float
        Quadratic-penalty weight :math:`w`.  Default 50.0; large values
        approach an exact projection but increase L-BFGS condition
        number.
    eps : float
        Small floor to avoid division by zero when two atoms collapse.

    Returns
    -------
    loss : float
        Scalar penalty value.
    grad : (n, 3) ndarray
        Analytical gradient :math:`\\partial \\phi / \\partial p_i`.
    """
    n = int(P.shape[0])
    grad = np.zeros_like(P)
    if n < 2:
        return 0.0, grad

    # Vectorised pair-distance computation (memory ~O(n^2) — acceptable
    # for n < 200 which covers every reasonable TMC).
    diff = P[:, None, :] - P[None, :, :]               # (n, n, 3)
    dist = np.linalg.norm(diff, axis=-1)               # (n, n)
    dist_safe = np.maximum(dist, eps)

    # Lower-bound violations: where dist < lower.  Active mask + magnitudes.
    iu = np.triu_indices(n, k=1)
    d_pair = dist[iu]
    l_pair = lower[iu]
    u_pair = upper[iu]

    # Lower-bound violation
    viol_lo = np.maximum(0.0, l_pair - d_pair)
    # Upper-bound violation — only for pairs with a finite upper bound
    upper_finite = u_pair < _UB_INF_THRESHOLD
    viol_hi = np.where(upper_finite, np.maximum(0.0, d_pair - u_pair), 0.0)

    loss = 0.5 * float(weight) * float(np.sum(viol_lo ** 2 + viol_hi ** 2))

    # Gradient.  For each pair (i, j) with d = ||p_i - p_j||:
    #   if d < l_ij:  phi = 0.5 w (l - d)^2
    #       dphi/dp_i = -w (l - d) * (p_i - p_j) / d
    #   if d > u_ij:  phi = 0.5 w (d - u)^2
    #       dphi/dp_i =  w (d - u) * (p_i - p_j) / d
    # Sum contributions over all violating pairs symmetrically.
    if loss > 0.0:
        i_idx = iu[0]
        j_idx = iu[1]
        # Direction vector p_i - p_j, length d
        dir_vec = diff[i_idx, j_idx]                       # (Npairs, 3)
        d_safe = dist_safe[i_idx, j_idx][:, None]          # (Npairs, 1)
        # combined signed-magnitude coefficient
        coef = (-viol_lo + viol_hi)[:, None]               # (Npairs, 1)
        force = float(weight) * coef * dir_vec / d_safe    # (Npairs, 3)
        # Scatter-add to atom i (+force) and atom j (-force)
        np.add.at(grad, i_idx, force)
        np.add.at(grad, j_idx, -force)

    return loss, grad


# ---------------------------------------------------------------------------
# M-D drift guard (SPEC §6)
# ---------------------------------------------------------------------------
def compute_md_drift(
    P_init: np.ndarray,
    P_final: np.ndarray,
    md_pairs: Optional[Sequence[Tuple[int, int]]],
) -> float:
    """Return the maximum |Δd_MD| in Å across the supplied M-D pairs.

    Returns 0.0 when ``md_pairs`` is empty or None (no guard active).
    """
    if not md_pairs:
        return 0.0
    P0 = np.asarray(P_init, dtype=np.float64)
    P1 = np.asarray(P_final, dtype=np.float64)
    max_drift = 0.0
    for (i, j) in md_pairs:
        ii, jj = int(i), int(j)
        d0 = float(np.linalg.norm(P0[ii] - P0[jj]))
        d1 = float(np.linalg.norm(P1[ii] - P1[jj]))
        drift = abs(d1 - d0)
        if drift > max_drift:
            max_drift = drift
    return float(max_drift)


# ---------------------------------------------------------------------------
# Determinism helpers — restart perturbation seeded by (seed, restart_idx).
# Mirrors the pattern used in :mod:`delfin.fffree.grip_polish` for byte
# identity with the wider GRIP determinism contract.
# ---------------------------------------------------------------------------
def _perturb_for_restart(
    P_init: np.ndarray,
    restart_idx: int,
    seed: int,
    amplitude: float,
    frozen_indices: Optional[np.ndarray],
) -> np.ndarray:
    """Generate a deterministic perturbation of ``P_init`` for restart ``k``.

    Restart 0 returns ``P_init`` verbatim.  For ``k >= 1`` a Mersenne-Twister
    seeded by ``seed + k`` produces a uniform :math:`[-A, A]` displacement
    per coordinate; frozen atoms get zero displacement so coordination
    constraints stay exactly satisfied.

    Determinism contract: same ``(P_init.shape, restart_idx, seed,
    amplitude, frozen_indices)`` → byte-identical output across two runs.
    """
    if restart_idx == 0 or amplitude <= 0.0:
        return np.asarray(P_init, dtype=np.float64).copy()
    rs_seed = (int(seed) + int(restart_idx)) & 0x7FFFFFFF
    rng = np.random.RandomState(seed=rs_seed)
    delta = rng.uniform(
        -float(amplitude), float(amplitude), size=P_init.shape
    ).astype(np.float64)
    if frozen_indices is not None and frozen_indices.size:
        delta[frozen_indices] = 0.0
    return np.asarray(P_init, dtype=np.float64) + delta


# ---------------------------------------------------------------------------
# Bounds violation count (post-solve diagnostic, also used for feasibility
# fallback trigger).
# ---------------------------------------------------------------------------
def _count_bounds_violations(
    P: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    abs_tol: float = 1e-3,
) -> int:
    """Count pairs (i<j) where |Δ| > abs_tol on either side."""
    n = int(P.shape[0])
    if n < 2:
        return 0
    diff = P[:, None, :] - P[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    iu = np.triu_indices(n, k=1)
    d_pair = dist[iu]
    l_pair = lower[iu]
    u_pair = upper[iu]
    upper_finite = u_pair < _UB_INF_THRESHOLD
    viol_lo = (l_pair - d_pair) > abs_tol
    viol_hi = upper_finite & ((d_pair - u_pair) > abs_tol)
    return int(np.sum(viol_lo | viol_hi))


# ---------------------------------------------------------------------------
# Feasibility-fallback bounds relaxation (SPEC §2.4)
# ---------------------------------------------------------------------------
def _relax_bounds(
    lower: np.ndarray,
    upper: np.ndarray,
    pair_classes: Optional[Mapping[Tuple[int, int], str]],
    md_pairs: Optional[Sequence[Tuple[int, int]]],
    dd_pairs: Optional[Sequence[Tuple[int, int]]],
    tier: int,
    info: SolverInfo,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return relaxed bounds matrices for a given relaxation tier.

    Tier 1: angle (1,3) bound widths × 1.5.
    Tier 2: torsion / non-bonded vdW lower bounds × 0.7-0.8 and widths × 2.0.
    Tier 3: vdW lower bound multiplied by 0.7 (final permissive tier).

    M-D and D-D pairs (passed explicitly OR classified as 'md' / 'dd' in
    ``pair_classes``) are NEVER relaxed.  When ``pair_classes`` is None,
    every non-M-D, non-D-D pair is treated as 'nonbonded' for tier 2/3.
    """
    L = lower.copy()
    U = upper.copy()
    n = int(L.shape[0])

    md_set = set(
        (min(int(i), int(j)), max(int(i), int(j))) for (i, j) in (md_pairs or [])
    )
    dd_set = set(
        (min(int(i), int(j)), max(int(i), int(j))) for (i, j) in (dd_pairs or [])
    )

    def _class(i: int, j: int) -> str:
        key = (min(i, j), max(i, j))
        if key in md_set:
            return "md"
        if key in dd_set:
            return "dd"
        if pair_classes is not None and key in pair_classes:
            return str(pair_classes[key])
        return "nonbonded"

    multipliers: Dict[str, Tuple[float, float]] = {}
    if tier >= 1:
        # angle widths × 1.5 — symmetric expansion around the midpoint
        multipliers["angle"] = (1.5, 1.5)
    if tier >= 2:
        # 1,4 torsion + non-bonded widths × 2.0, lower bound × 0.8
        multipliers["torsion"] = (2.0, 0.8)
        multipliers["nonbonded"] = (2.0, 0.8)
    if tier >= 3:
        # vdW lower bound × 0.7 (most permissive)
        multipliers["nonbonded"] = (2.0, 0.7)
        multipliers["vdw"] = (2.0, 0.7)

    if not multipliers:
        return L, U

    for i in range(n):
        for j in range(i + 1, n):
            cls = _class(i, j)
            if cls in ("md", "dd", "bonded"):
                continue
            mult = multipliers.get(cls)
            if mult is None:
                continue
            width_mult, low_mult = mult
            lo = float(L[i, j])
            hi = float(U[i, j])
            mid = 0.5 * (lo + hi) if hi < _UB_INF_THRESHOLD else lo + 0.5
            half = (hi - lo) * 0.5 if hi < _UB_INF_THRESHOLD else 0.5
            new_half = half * float(width_mult)
            new_lo = max(0.1, mid - new_half) * float(low_mult)
            new_hi = mid + new_half if hi < _UB_INF_THRESHOLD else hi
            new_lo = min(new_lo, lo)  # only relax (widen) — never tighten
            new_hi = max(new_hi, hi)
            if new_lo != lo or new_hi != hi:
                L[i, j] = L[j, i] = new_lo
                U[i, j] = U[j, i] = new_hi
                info.bounds_relaxed.append((int(i), int(j), f"tier{tier}:{cls}"))

    return L, U


# ---------------------------------------------------------------------------
# Augmented loss closure (severity + bounds penalty)
# ---------------------------------------------------------------------------
def _make_augmented_loss(
    severity_fn: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    lower: np.ndarray,
    upper: np.ndarray,
    bounds_weight: float,
    n_atoms: int,
) -> Callable[[np.ndarray], Tuple[float, np.ndarray]]:
    """Build the scipy-compatible objective ``f(x) -> (loss, grad_flat)``."""
    def _objective(x: np.ndarray) -> Tuple[float, np.ndarray]:
        P = np.asarray(x, dtype=np.float64).reshape(n_atoms, 3)
        try:
            emp_loss, emp_grad = severity_fn(P)
        except Exception as exc:  # fail-open contract
            # Return a high but finite loss and zero gradient — the solver
            # will treat this restart as non-converging and move on.
            return float(1e12), np.zeros(n_atoms * 3, dtype=np.float64)
        emp_grad = np.asarray(emp_grad, dtype=np.float64).reshape(n_atoms, 3)
        pen_loss, pen_grad = bounds_violation_penalty(
            P, lower, upper, weight=bounds_weight,
        )
        total_loss = float(emp_loss) + float(pen_loss)
        total_grad = (emp_grad + pen_grad).reshape(-1)
        return float(total_loss), total_grad

    return _objective


# ---------------------------------------------------------------------------
# Single L-BFGS-B call from a given start.  Returns (P_refined, n_iter,
# final_loss, grad_norm).  Failures are swallowed: a NaN start or solver
# breakdown returns the input untouched with `failed=True` signalled via
# inf loss.
# ---------------------------------------------------------------------------
def _single_lbfgs_call(
    P_start: np.ndarray,
    objective: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    max_iter: int,
    tol: float,
) -> Tuple[np.ndarray, int, float, float]:
    """Run scipy.minimize(L-BFGS-B) and return diagnostics."""
    if not _HAS_SCIPY:
        # Fallback: pure gradient-descent (deterministic) when scipy is missing.
        return _gradient_descent(P_start, objective, max_iter, tol)
    n_atoms = int(P_start.shape[0])
    x0 = P_start.reshape(-1).astype(np.float64)
    try:
        res = _scipy_minimize(
            objective,
            x0,
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": int(max_iter),
                "gtol": float(tol),
                "ftol": float(tol),
            },
        )
    except Exception:
        return P_start.copy(), 0, float("inf"), float("inf")
    x_opt = np.asarray(res.x, dtype=np.float64).reshape(n_atoms, 3)
    n_iter = int(getattr(res, "nit", 0))
    final_loss = float(getattr(res, "fun", float("inf")))
    grad = np.asarray(getattr(res, "jac", np.zeros(x_opt.size)), dtype=np.float64)
    grad_norm = float(np.linalg.norm(grad))
    return x_opt, n_iter, final_loss, grad_norm


def _gradient_descent(
    P_start: np.ndarray,
    objective: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    max_iter: int,
    tol: float,
    lr: float = 1e-3,
) -> Tuple[np.ndarray, int, float, float]:
    """Deterministic gradient-descent fallback when scipy is unavailable."""
    P = P_start.copy()
    n_atoms = int(P.shape[0])
    last_loss = float("inf")
    grad_norm = float("inf")
    for it in range(int(max_iter)):
        loss, grad = objective(P.reshape(-1))
        grad = grad.reshape(n_atoms, 3)
        grad_norm = float(np.linalg.norm(grad))
        if grad_norm < tol:
            return P, it, loss, grad_norm
        P = P - lr * grad
        last_loss = loss
    return P, int(max_iter), float(last_loss), grad_norm


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def solve_dg(
    P_init: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    severity_fn: Callable[[np.ndarray], Tuple[float, np.ndarray]],
    *,
    max_iter: int = DEFAULT_MAX_ITER,
    tol: float = DEFAULT_TOL,
    n_restarts: int = DEFAULT_N_RESTARTS,
    seed: int = DEFAULT_SEED,
    bounds_weight: float = DEFAULT_BOUNDS_WEIGHT,
    md_pairs: Optional[Sequence[Tuple[int, int]]] = None,
    dd_pairs: Optional[Sequence[Tuple[int, int]]] = None,
    pair_classes: Optional[Mapping[Tuple[int, int], str]] = None,
    frozen_indices: Optional[Sequence[int]] = None,
    md_drift_tol: float = DEFAULT_MD_DRIFT_TOL,
    perturb_amplitude: float = DEFAULT_PERTURB_AMPLITUDE,
    allow_relaxation: bool = True,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Projected L-BFGS solver for the whole-complex distance-geometry problem.

    See module-level docstring for the algorithm overview.

    Parameters
    ----------
    P_init : (n, 3) ndarray
        Initial point cloud.  Typically comes from ``assemble_complex`` or
        a previous polish.  Must not contain NaNs.
    lower, upper : (n, n) ndarray
        Symmetric hard-bound matrices from
        :func:`delfin.fffree.mogul_bounds.build_bounds_matrix`.  Pairs
        with ``upper >= 5e5`` are treated as effectively unbounded above.
    severity_fn : callable
        ``severity_fn(P) -> (loss, grad)`` where ``P`` is the ``(n, 3)``
        position array and ``grad`` is ``(n, 3)``.  This is the empirical
        Mahalanobis objective from Phase C (or any other ``f, gf`` pair).
    max_iter : int
        L-BFGS iteration cap.  SPEC default 1000.
    tol : float
        Convergence on |grad|.  SPEC default 1e-5.
    n_restarts : int
        Number of L-BFGS calls; restart 0 starts from ``P_init`` exactly,
        restart k ≥ 1 starts from a deterministically perturbed copy.
        SPEC default 3.
    seed : int
        Master seed for the restart perturbation.  Restart k uses
        ``seed + k``.
    bounds_weight : float
        Quadratic-penalty weight :math:`w` for bounds violations.  Larger
        approaches an exact projection but slows L-BFGS conditioning.
    md_pairs : sequence of (i, j), optional
        Metal-donor pairs.  When provided, the post-solve M-D drift guard
        rejects restarts whose maximum ``|d_md_after - d_md_before|`` >
        ``md_drift_tol``.  M-D pairs are also held hard during the
        feasibility-fallback relaxation tiers.
    dd_pairs : sequence of (i, j), optional
        Donor-donor pairs.  Held hard during relaxation tiers.
    pair_classes : mapping (i, j) -> str, optional
        Per-pair classification (``"bonded"``, ``"angle"``, ``"torsion"``,
        ``"vdw"``, ``"nonbonded"``, ``"md"``, ``"dd"``).  Drives the
        relaxation-tier targeting.  When omitted, the relaxer uses
        conservative defaults (every non-M-D, non-D-D pair gets the
        ``nonbonded`` treatment).
    frozen_indices : sequence of int, optional
        Atom indices that must not be perturbed across restarts (typically
        the metal + donors so the coordination geometry is preserved
        across restarts).
    md_drift_tol : float
        Maximum allowed |Δd_MD| in Å.  Default 0.5 (SPEC §6).
    perturb_amplitude : float
        Uniform perturbation amplitude in Å for restart ≥ 1.
    allow_relaxation : bool
        If True (default), failed restarts trigger hierarchical bounds
        relaxation before giving up.  Set False to short-circuit relaxer
        in unit tests of the infeasibility-detection path.

    Returns
    -------
    P_solved : (n, 3) ndarray
        The best converged geometry, or ``P_init.copy()`` on total
        failure (fail-open contract).
    info : dict
        Diagnostic counters; see :class:`SolverInfo` for fields.  On
        failure ``info['failed'] = True`` and the caller should fall
        back to the rigid-fit path (SPEC §13 Q3).
    """
    P_init = np.asarray(P_init, dtype=np.float64)
    lower = np.asarray(lower, dtype=np.float64)
    upper = np.asarray(upper, dtype=np.float64)

    info = SolverInfo()

    # ------------------------------------------------------------------ sanity
    if P_init.ndim != 2 or P_init.shape[1] != 3:
        info.failed = True
        info.reason = "P_init must have shape (n, 3)"
        return P_init.copy(), info.as_dict()
    n_atoms = int(P_init.shape[0])
    if n_atoms < 2:
        info.failed = True
        info.reason = "fewer than 2 atoms"
        return P_init.copy(), info.as_dict()
    if lower.shape != (n_atoms, n_atoms) or upper.shape != (n_atoms, n_atoms):
        info.failed = True
        info.reason = "bounds shape mismatch"
        return P_init.copy(), info.as_dict()
    if not np.all(np.isfinite(P_init)):
        info.failed = True
        info.reason = "non-finite P_init"
        return P_init.copy(), info.as_dict()

    # Normalise optional arguments
    if frozen_indices is not None:
        frozen_arr = np.asarray(sorted(set(int(i) for i in frozen_indices)), dtype=np.int64)
    else:
        frozen_arr = np.array([], dtype=np.int64)

    md_pairs_norm: List[Tuple[int, int]] = []
    if md_pairs is not None:
        for (i, j) in md_pairs:
            a, b = int(i), int(j)
            if a == b:
                continue
            md_pairs_norm.append((min(a, b), max(a, b)))
        md_pairs_norm = sorted(set(md_pairs_norm))

    dd_pairs_norm: List[Tuple[int, int]] = []
    if dd_pairs is not None:
        for (i, j) in dd_pairs:
            a, b = int(i), int(j)
            if a == b:
                continue
            dd_pairs_norm.append((min(a, b), max(a, b)))
        dd_pairs_norm = sorted(set(dd_pairs_norm))

    n_restarts = max(1, int(n_restarts))

    # ------------------------------------------------------------------
    # Restart loop, with optional relaxation tier escalation
    # ------------------------------------------------------------------
    best_P: Optional[np.ndarray] = None
    best_loss: float = float("inf")
    best_restart: int = -1
    best_n_iter: int = 0
    best_grad_norm: float = float("inf")
    best_drift: float = 0.0
    best_violations: int = -1

    max_relax_tier = 3 if allow_relaxation else 0

    for tier in range(max_relax_tier + 1):
        if tier == 0:
            L_use, U_use = lower, upper
        else:
            L_use, U_use = _relax_bounds(
                lower, upper, pair_classes, md_pairs_norm, dd_pairs_norm,
                tier, info,
            )
            info.relaxation_tier = tier

        objective = _make_augmented_loss(
            severity_fn, L_use, U_use, bounds_weight, n_atoms,
        )

        for restart_idx in range(n_restarts):
            P_start = _perturb_for_restart(
                P_init, restart_idx, seed, perturb_amplitude, frozen_arr,
            )
            P_try, n_iter_try, loss_try, grad_norm_try = _single_lbfgs_call(
                P_start, objective, max_iter, tol,
            )
            info.n_restarts_tried += 1

            # Skip non-finite results outright
            if not np.all(np.isfinite(P_try)):
                continue

            # M-D drift guard
            drift = compute_md_drift(P_init, P_try, md_pairs_norm)
            if md_drift_tol > 0.0 and drift > md_drift_tol:
                continue

            # Count bounds violations (use ORIGINAL bounds — relaxed-tier
            # successes are still infeasible on the original spec).
            n_viol = _count_bounds_violations(P_try, lower, upper)

            # Prefer feasible-on-original then lowest loss.
            is_better = False
            if best_violations < 0:
                is_better = True
            elif n_viol == 0 and best_violations > 0:
                is_better = True
            elif n_viol == best_violations and loss_try < best_loss:
                is_better = True
            elif n_viol < best_violations:
                is_better = True

            if is_better:
                best_P = P_try
                best_loss = loss_try
                best_restart = restart_idx
                best_n_iter = n_iter_try
                best_grad_norm = grad_norm_try
                best_drift = drift
                best_violations = n_viol

        # If we have a strictly feasible solution on the ORIGINAL bounds,
        # stop escalating relaxation tiers.
        if best_P is not None and best_violations == 0:
            break

    # ------------------------------------------------------------------
    # Wrap up info dict
    # ------------------------------------------------------------------
    if best_P is None:
        info.failed = True
        info.reason = "all restarts produced non-finite output or exceeded M-D drift"
        return P_init.copy(), info.as_dict()

    # Split final loss into empirical vs penalty contributions for forensik
    try:
        emp_loss, _ = severity_fn(best_P)
    except Exception:
        emp_loss = float("nan")
    pen_loss, _ = bounds_violation_penalty(
        best_P, lower, upper, weight=bounds_weight,
    )

    info.n_iter = int(best_n_iter)
    info.final_loss = float(best_loss)
    info.final_emp_loss = float(emp_loss)
    info.final_penalty = float(pen_loss)
    info.max_md_drift = float(best_drift)
    info.n_bounds_violated = int(max(0, best_violations))
    info.restart_used = int(best_restart)
    info.grad_norm = float(best_grad_norm)
    # Converged when grad-norm is below tol AND zero violations on the
    # original bounds.  Otherwise we honour the call but flag.
    info.converged = bool(
        info.grad_norm <= max(tol * 10.0, 1e-3) and info.n_bounds_violated == 0
    )
    info.failed = bool(info.n_bounds_violated > 0)
    if info.failed:
        info.reason = (
            f"{info.n_bounds_violated} bounds violations remain after "
            f"{info.n_restarts_tried} restarts / tier {info.relaxation_tier}"
        )

    return best_P, info.as_dict()


# ---------------------------------------------------------------------------
# Module-level smoke test (run only when executed directly, never on import)
# ---------------------------------------------------------------------------
def _self_check() -> None:
    """Tiny end-to-end smoke for local debugging.

    Three atoms on a line; target distances 1.0 / 1.0 / 2.0 (= bonded sum).
    The solver should pull a random initial cloud onto the equilibrium.
    """
    rng = np.random.RandomState(0)
    P0 = rng.uniform(-1.0, 1.0, size=(3, 3))
    L = np.zeros((3, 3))
    U = np.full((3, 3), 1e6)
    targets = {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 2.0}
    for (i, j), d in targets.items():
        L[i, j] = L[j, i] = d - 0.02
        U[i, j] = U[j, i] = d + 0.02

    def sev_fn(P: np.ndarray) -> Tuple[float, np.ndarray]:
        # Zero empirical term — the bounds penalty alone drives the solve.
        return 0.0, np.zeros_like(P)

    P_out, info = solve_dg(P0, L, U, sev_fn, n_restarts=2)
    assert info["n_bounds_violated"] == 0, info


if __name__ == "__main__":  # pragma: no cover
    _self_check()
    print("mogul_solver self-check passed.")
