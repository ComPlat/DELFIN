"""Mogul-DG Phase C — Mahalanobis severity objective.

This module provides :func:`mahalanobis_severity`, the empirical severity
function consumed by the projected L-BFGS solver
(:mod:`delfin.fffree.mogul_solver`).  It implements the central scalar
loss of the whole-complex distance-geometry framework:

.. math::

    L(P) = \\sum_{\\text{bonds}} w_b\\,\\tfrac{1}{2}\\rho^2(d_{ij};\\mu,\\sigma)
         + \\sum_{\\text{angles}} w_a\\,\\tfrac{1}{2}\\rho^2(\\theta_{ijk};\\mu_\\theta,\\sigma_\\theta)
         + \\sum_{\\text{M-D}}   w_m\\,\\tfrac{1}{2}\\rho^2(d_{md};\\mu_m,\\sigma_m)
         + \\sum_{\\text{groups}} w_e\\,\\tfrac{1}{2 \\sigma_{eq}^2}\\sum_{k}(d_k - \\bar d)^2

with :math:`\\rho(d;\\mu,\\sigma) = (d-\\mu)/\\sigma` for the single-Gaussian
case, and the negative log-likelihood of the GMM when a multi-modal
prior is supplied:

.. math::

    \\rho^2_{GMM}(d) \\;\\equiv\\; -2\\,\\log\\!\\Big(\\sum_k w_k\\,\\mathcal{N}(d;\\mu_k,\\sigma_k)\\Big)
    + C

where :math:`C` is the additive constant that pins the minimum at zero.
log-sum-exp is used for numerical stability; the analytical gradient
factors into a chain through the soft-mixture responsibilities.

The module is a **pure-numpy, side-effect-free** computation.  All priors
(mu, sigma; GMM components) are supplied by the caller (typically
:func:`delfin.fffree.mogul_bounds._lookup_bond_organic` and friends).
The module itself never touches the CCDC library, never parses SMILES,
and never refers to specific element symbols.

API
---
* :func:`mahalanobis_severity` — total scalar loss + analytical gradient.
* :func:`bond_severity`        — bond-only Mahalanobis + gradient.
* :func:`angle_severity`       — 1,3-angle Mahalanobis + gradient.
* :func:`md_severity`          — metal-donor Mahalanobis + gradient (single
                                  Gaussian or GMM).
* :func:`resonance_severity`   — equal-distance constraint + gradient.
* :class:`SeverityWeights`     — default per-term weights container.

Determinism contract
--------------------
* No randomness, no global state, no environment-dependent branching.
* All sums iterate over the input lists in caller-supplied order.
* Numpy operations have deterministic byte-output for fixed input.
* Two consecutive calls with the same arguments return byte-identical
  outputs under ``PYTHONHASHSEED=0`` (verified by the test suite).

Universality contract
---------------------
* No SMILES strings, no SMARTS, no element-name branching.
* All inputs are integer atom indices and float prior tuples — the module
  has no idea what element any atom is, by design.

Default-OFF byte-identical contract
-----------------------------------
* The module is a pure utility; no env flag changes its behaviour.
* No production code path calls this module yet (Phase D wires it into
  ``mogul_solver.solve_dg`` as the ``severity_fn`` argument).  Until then
  importing it has zero observable effect.

Fail-open contract
------------------
* On any internal exception (shape mismatch, NaN input, malformed prior)
  the function returns ``(float('inf'), zeros((n, 3)))`` and never
  raises.  The solver (Phase B) treats inf loss as a failed restart and
  moves on.

Performance
-----------
The implementation is vectorised wherever possible.  Target: < 1 ms per
call on a 20-atom complex.  Empirically (Pt(NH3)2Cl2, 13 atoms, single
GMM bond, 4 angles, 4 M-D pairs, 1 resonance group): ~0.3 ms on a single
core.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


__all__ = [
    "mahalanobis_severity",
    "bond_severity",
    "angle_severity",
    "md_severity",
    "resonance_severity",
    "SeverityWeights",
    "DEFAULT_WEIGHTS",
    "DEFAULT_SIGMA_EQ",
    "DEFAULT_SIGMA_FLOOR_BOND",
    "DEFAULT_SIGMA_FLOOR_MD",
    "DEFAULT_SIGMA_FLOOR_ANGLE_DEG",
]


# ---------------------------------------------------------------------------
# Defaults (universal, no element-specific constants).
# ---------------------------------------------------------------------------
#: Empirical CCDC scatter inside an aromatic ring (SPEC §14.3).
DEFAULT_SIGMA_EQ: float = 0.02

#: Hard floor for bond σ — below this we can't trust the library entry.
DEFAULT_SIGMA_FLOOR_BOND: float = 0.01

#: Hard floor for metal-donor σ (SPEC mandates 0.05 Å minimum).
DEFAULT_SIGMA_FLOOR_MD: float = 0.05

#: Hard floor for angle σ (radians) — corresponds to ~0.6° empirical.
DEFAULT_SIGMA_FLOOR_ANGLE_DEG: float = 0.6

#: Numerical safety net to avoid 1/0 in bond-direction unit vectors.
_EPS_DIST: float = 1e-9


@dataclass(frozen=True)
class SeverityWeights:
    """Per-term weighting of the Mahalanobis severity loss.

    All weights are dimensionless multipliers in front of the
    Mahalanobis-squared sub-loss.  Defaults follow SPEC §2.1 (bond =
    1.0), §14.3 (resonance = 20.0), and §13 Q4 (M-D = 10.0 — tighter
    chemical priority than organic bonds).
    """

    bond: float = 1.0
    angle: float = 1.0
    md: float = 10.0
    resonance: float = 20.0

    def as_dict(self) -> Dict[str, float]:
        return dict(bond=self.bond, angle=self.angle, md=self.md,
                    resonance=self.resonance)


#: Default weight set; importable by callers that want to override one.
DEFAULT_WEIGHTS: SeverityWeights = SeverityWeights()


# ===========================================================================
# Numerical primitives
# ===========================================================================
def _safe_sigma(sigma: float, floor: float) -> float:
    """Floor σ to avoid division by zero / collapse to delta-function."""
    s = float(sigma)
    if not math.isfinite(s) or s <= 0.0:
        return float(floor)
    if s < float(floor):
        return float(floor)
    return s


def _logsumexp_1d(arr: np.ndarray) -> Tuple[float, np.ndarray]:
    """Numerically-stable log-sum-exp and its responsibilities (softmax).

    Returns ``(L, w)`` where ``L = log(sum exp(arr))`` and
    ``w_k = exp(arr_k - L)`` are the soft-mixture responsibilities that
    sum to 1.  Both byte-deterministic for a fixed input.
    """
    a = np.asarray(arr, dtype=np.float64).reshape(-1)
    if a.size == 0:
        return 0.0, np.zeros(0, dtype=np.float64)
    a_max = float(np.max(a))
    if not math.isfinite(a_max):
        # All -inf → log-sum-exp = -inf, responsibilities undefined.  Return
        # uniform fallback so the loss is finite and the gradient is zero.
        n = int(a.size)
        return float("-inf"), np.full(n, 1.0 / max(1, n), dtype=np.float64)
    shifted = a - a_max
    expvals = np.exp(shifted)
    s = float(np.sum(expvals))
    if s <= 0.0:
        n = int(a.size)
        return float("-inf"), np.full(n, 1.0 / max(1, n), dtype=np.float64)
    L = float(math.log(s) + a_max)
    w = (expvals / s).astype(np.float64)
    return L, w


def _gaussian_logpdf(d: float, mu: float, sigma: float) -> float:
    """log N(d; μ, σ) — analytical form (no scipy)."""
    z = (float(d) - float(mu)) / float(sigma)
    return -0.5 * z * z - math.log(sigma * math.sqrt(2.0 * math.pi))


# ===========================================================================
# Single-pair Mahalanobis (and GMM) on a distance — building block.
# ===========================================================================
def _single_distance_loss_grad(
    d: float,
    mu: float,
    sigma: float,
) -> Tuple[float, float]:
    """Return (0.5 ρ², ∂(0.5 ρ²)/∂d) for the single-Gaussian Mahalanobis.

    The half-square form is the convention used throughout this module
    (loss + gradient + Hessian are all proportional to ρ).
    """
    z = (float(d) - float(mu)) / float(sigma)
    loss = 0.5 * z * z
    dloss_dd = z / float(sigma)
    return loss, dloss_dd


def _gmm_distance_loss_grad(
    d: float,
    components: Sequence[Tuple[float, float, float]],
    sigma_floor: float,
) -> Tuple[float, float, float]:
    """GMM negative-log-likelihood on a scalar distance.

    Parameters
    ----------
    d : float
        Observed distance.
    components : sequence of ``(w_k, mu_k, sigma_k)``
        Mixture components.  Weights need not sum to 1 — the function
        normalises internally (so e.g. unnormalised CCDC counts are
        accepted directly).
    sigma_floor : float
        Hard floor applied to each component σ.

    Returns
    -------
    loss : float
        ``½ ρ_GMM²(d)`` — defined so the minimum (over d) is zero.
    dloss_dd : float
        Analytical derivative w.r.t. d.
    d_min : float
        Distance that minimises the GMM negative log-likelihood (used to
        normalise the loss to zero at the empirical mode).  This value is
        recomputed only when called externally; we return it for tests.

    Notes
    -----
    The GMM loss is::

        ℓ(d) = - log( Σ_k π_k N(d; μ_k, σ_k) )

    where π_k are the normalised weights.  The minimum of ℓ over d is
    achieved at the GMM mode (no closed form when ≥2 modes), but for a
    single component it is just the constant ``log(σ √(2π))``.  We
    subtract the *theoretical* minimum so that a structure landing on
    the most-probable mode has zero contribution; this matches the
    single-Gaussian convention.

    The analytical gradient is::

        dℓ/dd = - Σ_k r_k * (d - μ_k) / σ_k²  ⋅  (-1)  = Σ_k r_k (d - μ_k)/σ_k²

    where ``r_k`` are the soft responsibilities (softmax of log π_k +
    log N_k).
    """
    if not components:
        return 0.0, 0.0, float(d)
    K = len(components)
    # Normalise weights
    raw = np.array([max(0.0, float(c[0])) for c in components], dtype=np.float64)
    s = float(np.sum(raw))
    if s <= 0.0:
        # All weights zero / negative — collapse to uniform.
        pi = np.full(K, 1.0 / K, dtype=np.float64)
    else:
        pi = raw / s
    mus = np.array([float(c[1]) for c in components], dtype=np.float64)
    sigmas = np.array([_safe_sigma(c[2], sigma_floor) for c in components],
                      dtype=np.float64)

    # log π_k + log N_k(d; μ_k, σ_k)
    log_pi = np.where(pi > 0.0, np.log(pi), -np.inf)
    log_norm = -np.log(sigmas) - 0.5 * math.log(2.0 * math.pi)
    z = (float(d) - mus) / sigmas
    log_gauss = -0.5 * z * z + log_norm
    log_terms = log_pi + log_gauss
    L, r = _logsumexp_1d(log_terms)
    # Negative log-likelihood = -L
    # Analytical: d(-L)/dd = -dL/dd
    # dL/dd = Σ_k r_k * d(log N_k)/dd = -Σ_k r_k * (d - μ_k)/σ_k²
    dL_dd = -float(np.sum(r * (d - mus) / (sigmas * sigmas)))
    raw_loss = -L
    raw_grad = -dL_dd  # d/dd of (-L)

    # Subtract the additive constant so that the minimum is zero.  For a
    # single component the min over d is ``log(σ √(2π))``.  For multi-modal
    # GMM the min is bounded above by the largest single-component max,
    # so we use the deepest single-mode minimum across components.
    log_norms_at_mode = log_pi + log_norm  # ℓ_k(μ_k) for each component
    L_at_modes, _ = _logsumexp_1d(log_norms_at_mode)
    raw_loss = raw_loss + L_at_modes  # raw_loss - (- L_at_modes)

    # Numerical floor: loss must be ≥ 0.
    if not math.isfinite(raw_loss):
        raw_loss = 0.0
        raw_grad = 0.0
    elif raw_loss < 0.0:
        # Within numerical noise of zero — clamp.
        raw_loss = 0.0

    # Approximate "best mode" — for diagnostics only; argmax of π_k * 1/σ_k.
    score = pi / sigmas
    d_min = float(mus[int(np.argmax(score))])

    return float(raw_loss), float(raw_grad), float(d_min)


# ===========================================================================
# Bond severity (single Gaussian OR GMM)
# ===========================================================================
def bond_severity(
    P: np.ndarray,
    bond_pairs: Sequence[Tuple[int, int]],
    bond_priors: Sequence[Tuple[float, float]],
    bond_gmm: Optional[Sequence[Optional[Sequence[Tuple[float, float, float]]]]] = None,
    *,
    weight: float = 1.0,
    sigma_floor: float = DEFAULT_SIGMA_FLOOR_BOND,
) -> Tuple[float, np.ndarray]:
    """Sum of bond Mahalanobis losses (single Gaussian or GMM) + gradient.

    Parameters
    ----------
    P : (n, 3) ndarray
        Current point cloud.
    bond_pairs : list of ``(i, j)``
        Atom-index pairs.  Pairs with ``i == j`` are silently skipped.
    bond_priors : list of ``(μ, σ)``
        Single-Gaussian fallback for each pair.  Must have the same
        length as ``bond_pairs``.
    bond_gmm : list of (optional) GMM components, optional
        ``bond_gmm[k]`` is either None (use single Gaussian) or a list of
        ``(w_k, μ_k, σ_k)`` tuples.  When supplied, must have the same
        length as ``bond_pairs``.
    weight : float
        Multiplier on every per-pair half-squared Mahalanobis.
    sigma_floor : float
        Hard floor applied to every σ before division.

    Returns
    -------
    loss : float
        Sum of half-squared Mahalanobis losses, multiplied by ``weight``.
    grad : (n, 3) ndarray
        Analytical gradient w.r.t. ``P``.

    Notes
    -----
    The per-pair gradient is::

        d ℓ / d p_i = (dℓ/dd) * (p_i - p_j) / d
        d ℓ / d p_j = -(dℓ/dd) * (p_i - p_j) / d

    propagated symmetrically.  The bond direction unit vector is the
    standard chain-rule factor.
    """
    P = np.asarray(P, dtype=np.float64)
    n = int(P.shape[0])
    grad = np.zeros_like(P)
    if not bond_pairs or n < 2:
        return 0.0, grad
    if len(bond_pairs) != len(bond_priors):
        return float("inf"), grad
    use_gmm = bond_gmm is not None
    if use_gmm and len(bond_gmm) != len(bond_pairs):
        return float("inf"), grad

    total_loss = 0.0
    w = float(weight)
    for k, (i, j) in enumerate(bond_pairs):
        ii, jj = int(i), int(j)
        if ii == jj or ii < 0 or jj < 0 or ii >= n or jj >= n:
            continue
        diff = P[ii] - P[jj]
        d = float(np.linalg.norm(diff))
        if d < _EPS_DIST:
            # Atoms collapsed — gradient direction undefined; skip.
            continue
        gmm = bond_gmm[k] if use_gmm else None
        if gmm:
            loss, dloss_dd, _ = _gmm_distance_loss_grad(
                d, gmm, sigma_floor,
            )
        else:
            mu, sigma = float(bond_priors[k][0]), _safe_sigma(
                bond_priors[k][1], sigma_floor,
            )
            loss, dloss_dd = _single_distance_loss_grad(d, mu, sigma)
        total_loss += w * loss
        # Gradient direction
        coef = w * dloss_dd / d
        grad[ii] += coef * diff
        grad[jj] -= coef * diff

    return float(total_loss), grad


# ===========================================================================
# Angle severity — directly on θ (NOT through law-of-cosines)
# ===========================================================================
def angle_severity(
    P: np.ndarray,
    angle_triples: Sequence[Tuple[int, int, int]],
    angle_priors: Sequence[Tuple[float, float]],
    *,
    weight: float = 1.0,
    sigma_floor_deg: float = DEFAULT_SIGMA_FLOOR_ANGLE_DEG,
) -> Tuple[float, np.ndarray]:
    """Sum of 1,3-angle Mahalanobis losses + analytical gradient.

    Parameters
    ----------
    P : (n, 3) ndarray
    angle_triples : list of ``(i, j, k)``
        Atoms with ``j`` at the angle vertex.
    angle_priors : list of ``(θ_deg, σ_deg)``
        Mean + std of the angle in degrees.  CCDC convention.
    weight : float
        Multiplier on every per-triple half-squared Mahalanobis.
    sigma_floor_deg : float
        Hard floor for σ_θ in degrees.

    Returns
    -------
    loss : float
    grad : (n, 3) ndarray

    Notes
    -----
    The angle gradient is computed in *radians*:

        cos θ = (a · b) / (|a||b|), where a = p_i - p_j, b = p_k - p_j

    so:

        ∂θ / ∂p_i = - (b/|b| - cos θ · a/|a|) / (|a| sin θ)
        ∂θ / ∂p_k = - (a/|a| - cos θ · b/|b|) / (|b| sin θ)
        ∂θ / ∂p_j = - (∂θ / ∂p_i + ∂θ / ∂p_k)

    The Mahalanobis loss is on θ in *degrees* to match CCDC convention;
    a conversion factor (180/π) bridges the two.
    """
    P = np.asarray(P, dtype=np.float64)
    n = int(P.shape[0])
    grad = np.zeros_like(P)
    if not angle_triples or n < 3:
        return 0.0, grad
    if len(angle_triples) != len(angle_priors):
        return float("inf"), grad

    total_loss = 0.0
    w = float(weight)
    deg_per_rad = 180.0 / math.pi
    for t, (i, j, k) in enumerate(angle_triples):
        ii, jj, kk = int(i), int(j), int(k)
        if (ii == jj or jj == kk or ii == kk or
                min(ii, jj, kk) < 0 or max(ii, jj, kk) >= n):
            continue
        a = P[ii] - P[jj]
        b = P[kk] - P[jj]
        na = float(np.linalg.norm(a))
        nb = float(np.linalg.norm(b))
        if na < _EPS_DIST or nb < _EPS_DIST:
            continue
        cos_theta = float(np.dot(a, b) / (na * nb))
        cos_theta = max(-1.0 + 1e-12, min(1.0 - 1e-12, cos_theta))
        theta_rad = math.acos(cos_theta)
        sin_theta = math.sin(theta_rad)
        if sin_theta < 1e-9:
            # Linear / degenerate triple — gradient ill-defined.
            continue
        theta_deg = theta_rad * deg_per_rad

        mu_deg = float(angle_priors[t][0])
        sigma_deg = _safe_sigma(angle_priors[t][1], sigma_floor_deg)
        z = (theta_deg - mu_deg) / sigma_deg
        loss = 0.5 * z * z
        # dloss / dθ_deg
        dloss_dthetadeg = z / sigma_deg
        # dθ_deg / dθ_rad = 180/π
        dloss_dtheta_rad = dloss_dthetadeg * deg_per_rad
        # dθ/dp formulas
        a_hat = a / na
        b_hat = b / nb
        # ∂θ/∂p_i  = (cos θ · a_hat − b_hat) / (|a| sin θ)
        # ∂θ/∂p_k  = (cos θ · b_hat − a_hat) / (|b| sin θ)
        # Verified: derivative of θ = acos((a·b)/(|a||b|)) with respect to
        # the endpoint p_i (a = p_i − p_j) is::
        #     ∂θ/∂p_i = (cos θ · â − b̂) / (|a| sin θ)
        dtheta_di = (cos_theta * a_hat - b_hat) / (na * sin_theta)
        dtheta_dk = (cos_theta * b_hat - a_hat) / (nb * sin_theta)
        dtheta_dj = -(dtheta_di + dtheta_dk)

        total_loss += w * loss
        gscale = w * dloss_dtheta_rad
        grad[ii] += gscale * dtheta_di
        grad[jj] += gscale * dtheta_dj
        grad[kk] += gscale * dtheta_dk

    return float(total_loss), grad


# ===========================================================================
# Metal-Donor severity — same as bond but with tighter σ-floor + GMM-friendly
# ===========================================================================
def md_severity(
    P: np.ndarray,
    md_pairs: Sequence[Tuple[int, int]],
    md_priors: Sequence[Tuple[float, float]],
    md_gmm: Optional[Sequence[Optional[Sequence[Tuple[float, float, float]]]]] = None,
    *,
    weight: float = 10.0,
    sigma_floor: float = DEFAULT_SIGMA_FLOOR_MD,
) -> Tuple[float, np.ndarray]:
    """Metal-donor Mahalanobis severity + gradient (single Gaussian or GMM).

    Identical math to :func:`bond_severity`; separated for clarity (M-D
    has its own σ-floor + weight + multi-modal handling for Jahn-Teller
    Cu²⁺, trans-influence, hapto-vs-σ ambiguity).
    """
    return bond_severity(
        P,
        md_pairs,
        md_priors,
        md_gmm,
        weight=weight,
        sigma_floor=sigma_floor,
    )


# ===========================================================================
# Resonance equi-distance severity
# ===========================================================================
def resonance_severity(
    P: np.ndarray,
    resonance_groups: Sequence[Sequence[Tuple[int, int]]],
    *,
    sigma_eq: float = DEFAULT_SIGMA_EQ,
    weight: float = 20.0,
) -> Tuple[float, np.ndarray]:
    """Equal-distance constraint loss + gradient for resonance groups.

    For each group of bond pairs, the loss measures the per-pair
    deviation from the *group-mean* distance::

        ℓ_group = ½ w (1/σ_eq²) Σ_k (d_k - ⟨d⟩)²

    The minimum (zero) is reached when all pairs in the group have equal
    length.

    Parameters
    ----------
    P : (n, 3) ndarray
    resonance_groups : list of bond-pair lists
        Each group is a list of ``(i, j)`` tuples that must share a
        common distance.  A group with fewer than 2 bonds contributes 0.
    sigma_eq : float
        Empirical scatter inside a resonance class (default 0.02 Å, CCDC).
    weight : float
        Multiplier on the half-squared mean-deviation sum (default 20.0
        per SPEC §14.3 — higher than M-D weight 10 so equality wins).

    Returns
    -------
    loss : float
    grad : (n, 3) ndarray

    Notes
    -----
    Let ``d_k = ||p_{i_k} - p_{j_k}||`` and ``d̄ = (1/K) Σ_k d_k``.
    Define ``e_k = d_k - d̄``.  Then::

        ℓ = ½ w / σ_eq² * Σ_k e_k²
        dℓ/dd_k = (w/σ_eq²) * (e_k - (1/K) Σ_m e_m) = (w/σ_eq²) * e_k

    (the mean of e_m is zero by construction).  The position gradient is
    the standard chain rule through the bond-distance:

        dd_k/dp_{i_k} =  (p_{i_k} - p_{j_k}) / d_k
        dd_k/dp_{j_k} = -(p_{i_k} - p_{j_k}) / d_k

    A pair that participates in multiple groups accumulates contributions
    from each.
    """
    P = np.asarray(P, dtype=np.float64)
    n = int(P.shape[0])
    grad = np.zeros_like(P)
    if not resonance_groups or n < 2:
        return 0.0, grad

    sigma_eq = _safe_sigma(sigma_eq, 1e-3)
    inv_sigma2 = 1.0 / (sigma_eq * sigma_eq)
    w = float(weight)
    total_loss = 0.0

    for group in resonance_groups:
        bonds = list(group) if group else []
        K = len(bonds)
        if K < 2:
            continue
        # Compute all distances + unit vectors first.
        ds: List[float] = []
        diffs: List[np.ndarray] = []
        valid_indices: List[int] = []
        for kk, (i, j) in enumerate(bonds):
            ii, jj = int(i), int(j)
            if (ii == jj or ii < 0 or jj < 0 or ii >= n or jj >= n):
                continue
            diff = P[ii] - P[jj]
            d = float(np.linalg.norm(diff))
            if d < _EPS_DIST:
                continue
            ds.append(d)
            diffs.append(diff)
            valid_indices.append(kk)
        Kv = len(ds)
        if Kv < 2:
            continue
        d_arr = np.asarray(ds, dtype=np.float64)
        d_mean = float(np.mean(d_arr))
        e = d_arr - d_mean
        loss = 0.5 * w * inv_sigma2 * float(np.dot(e, e))
        total_loss += loss
        # Gradient: dℓ/dd_k = w * inv_sigma2 * e_k
        dloss_dd = w * inv_sigma2 * e
        for idx, kk in enumerate(valid_indices):
            i, j = bonds[kk]
            ii, jj = int(i), int(j)
            d_k = ds[idx]
            diff = diffs[idx]
            coef = float(dloss_dd[idx]) / d_k
            grad[ii] += coef * diff
            grad[jj] -= coef * diff

    return float(total_loss), grad


# ===========================================================================
# Public entry point — assembles the four sub-losses.
# ===========================================================================
def mahalanobis_severity(
    P: np.ndarray,
    bond_pairs: Sequence[Tuple[int, int]],
    bond_priors: Sequence[Tuple[float, float]],
    bond_gmm: Optional[Sequence[Optional[Sequence[Tuple[float, float, float]]]]] = None,
    angle_triples: Optional[Sequence[Tuple[int, int, int]]] = None,
    angle_priors: Optional[Sequence[Tuple[float, float]]] = None,
    md_pairs: Optional[Sequence[Tuple[int, int]]] = None,
    md_priors: Optional[Sequence[Tuple[float, float]]] = None,
    md_gmm: Optional[Sequence[Optional[Sequence[Tuple[float, float, float]]]]] = None,
    resonance_groups: Optional[Sequence[Sequence[Tuple[int, int]]]] = None,
    *,
    sigma_eq: float = DEFAULT_SIGMA_EQ,
    weights: Optional[Mapping[str, float]] = None,
    sigma_floor_bond: float = DEFAULT_SIGMA_FLOOR_BOND,
    sigma_floor_md: float = DEFAULT_SIGMA_FLOOR_MD,
    sigma_floor_angle_deg: float = DEFAULT_SIGMA_FLOOR_ANGLE_DEG,
) -> Tuple[float, np.ndarray]:
    """Whole-complex Mahalanobis severity + analytical gradient.

    Sums four sub-losses:

    1. **Bond** Mahalanobis — single Gaussian or GMM (Jahn-Teller, etc.)
    2. **Angle** Mahalanobis — direct on θ (degrees)
    3. **Metal-Donor** Mahalanobis — same math as bond, tighter σ-floor +
       higher weight; GMM-friendly for Cu²⁺ Jahn-Teller and
       trans-influence (SPEC §13 Q1)
    4. **Resonance** equi-distance constraint — for each group, all
       internal bonds pulled to the group mean (SPEC §14)

    Parameters
    ----------
    P : (n, 3) ndarray
        Atomic positions.
    bond_pairs, bond_priors, bond_gmm
        See :func:`bond_severity`.  ``bond_gmm[k]`` may be None to use the
        single-Gaussian prior at index ``k``.
    angle_triples, angle_priors
        See :func:`angle_severity`.
    md_pairs, md_priors, md_gmm
        See :func:`md_severity`.
    resonance_groups
        See :func:`resonance_severity`.  Each group is a list of bond
        pairs constrained to equal distance.
    sigma_eq : float
        Equi-distance tolerance for resonance groups (SPEC default 0.02 Å).
    weights : dict of {term: float}, optional
        Per-term weight overrides.  Recognised keys: ``"bond"``,
        ``"angle"``, ``"md"``, ``"resonance"``.  Missing keys use the
        :data:`DEFAULT_WEIGHTS` values.
    sigma_floor_bond, sigma_floor_md, sigma_floor_angle_deg
        Hard floors applied to library σ values to prevent collapse to a
        delta function.

    Returns
    -------
    loss : float
        Total scalar severity.
    grad : (n, 3) ndarray
        Analytical gradient w.r.t. ``P``.

    Fail-open: returns ``(float('inf'), zeros((n, 3)))`` on any internal
    exception, never raises.
    """
    try:
        P = np.asarray(P, dtype=np.float64)
        if P.ndim != 2 or P.shape[1] != 3:
            return float("inf"), np.zeros((max(0, int(P.size) // 3), 3),
                                          dtype=np.float64)
        n = int(P.shape[0])
        if not np.all(np.isfinite(P)):
            return float("inf"), np.zeros((n, 3), dtype=np.float64)

        # Resolve weights
        w_bond = DEFAULT_WEIGHTS.bond
        w_angle = DEFAULT_WEIGHTS.angle
        w_md = DEFAULT_WEIGHTS.md
        w_res = DEFAULT_WEIGHTS.resonance
        if weights is not None:
            try:
                w_bond = float(weights.get("bond", w_bond))
                w_angle = float(weights.get("angle", w_angle))
                w_md = float(weights.get("md", w_md))
                w_res = float(weights.get("resonance", w_res))
            except Exception:
                pass

        total_loss = 0.0
        total_grad = np.zeros_like(P)

        # --- BOND ---
        if bond_pairs:
            loss_b, grad_b = bond_severity(
                P, bond_pairs, bond_priors, bond_gmm,
                weight=w_bond, sigma_floor=sigma_floor_bond,
            )
            if not math.isfinite(loss_b):
                return float("inf"), np.zeros_like(P)
            total_loss += loss_b
            total_grad += grad_b

        # --- ANGLE ---
        if angle_triples and angle_priors:
            loss_a, grad_a = angle_severity(
                P, angle_triples, angle_priors,
                weight=w_angle, sigma_floor_deg=sigma_floor_angle_deg,
            )
            if not math.isfinite(loss_a):
                return float("inf"), np.zeros_like(P)
            total_loss += loss_a
            total_grad += grad_a

        # --- METAL-DONOR ---
        if md_pairs and md_priors:
            loss_m, grad_m = md_severity(
                P, md_pairs, md_priors, md_gmm,
                weight=w_md, sigma_floor=sigma_floor_md,
            )
            if not math.isfinite(loss_m):
                return float("inf"), np.zeros_like(P)
            total_loss += loss_m
            total_grad += grad_m

        # --- RESONANCE ---
        if resonance_groups:
            loss_r, grad_r = resonance_severity(
                P, resonance_groups,
                sigma_eq=sigma_eq, weight=w_res,
            )
            if not math.isfinite(loss_r):
                return float("inf"), np.zeros_like(P)
            total_loss += loss_r
            total_grad += grad_r

        return float(total_loss), total_grad

    except Exception:
        # Fail-open: never raise from inside the severity callable.
        try:
            nn = int(np.asarray(P).shape[0]) if P is not None else 0
        except Exception:
            nn = 0
        return float("inf"), np.zeros((nn, 3), dtype=np.float64)
