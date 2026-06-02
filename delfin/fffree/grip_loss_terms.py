"""GRIP — Loss + Gradient Terms (Phase 2, v1).

Implements the per-fragment terms of the Mahalanobis-style loss

.. math::

    L(R) = \\sum_{f} w_f \\left( \\frac{q_f(R) - \\mu_f}{\\sigma_f} \\right)^2

where :math:`q_f` is the geometric observable (bond length, valence angle,
out-of-plane displacement) and :math:`(\\mu_f, \\sigma_f)` come from the
CCDC-grounded library (:mod:`grip_mogul_lookup`).

Each term provides ``value_and_grad(R)`` returning the scalar loss and a
``(n_atoms, 3)`` gradient array.  All gradients are analytical (closed-form
chain rule) and validated against central finite differences in the tests.

Determinism contract (SPEC §3.3, §11):

* float64 everywhere
* no RNG, no nondeterministic library calls
* sorted iteration order over terms in :class:`TotalGripLoss`
* numerically stable: a small ``_EPS`` floor guards every divide-by-norm

Force-field-free contract (SPEC §2.2):

* no Lennard-Jones, no Coulomb, no parametric ``k_b``/``k_theta`` springs
* the ``1/sigma**2`` weighting is *data-derived* (empirical scatter), not
  parameterised — that is the philosophical line vs. a force field
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "BondTerm",
    "AngleTerm",
    "ImproperTerm",
    "TorsionTerm",
    "TotalGripLoss",
    "default_weights",
    "sparse_downweight",
]

# Numerical floor used in every gradient that involves dividing by a vector norm.
# Chosen so float64 catastrophic cancellation cannot produce a NaN gradient.
_EPS = 1e-12

# Default per-class weights (SPEC §2.4).
_W_BOND = 1.0
_W_ANGLE = 0.5
_W_IMPROPER = 2.0
_W_TORSION = 0.2  # for Phase 1b


def default_weights() -> dict:
    """Return the canonical per-class weight table (read-only reference)."""
    return {
        "bond": _W_BOND,
        "angle": _W_ANGLE,
        "improper": _W_IMPROPER,
        "torsion": _W_TORSION,
    }


def sparse_downweight(weight: float, n: int, n_floor: int = 5) -> float:
    """Apply the sparse-bin downweighting rule from SPEC §2.4.

    Fragments with ``n < n_floor`` CCDC observations get scaled down so they
    cannot overshoot well-sampled fragments.
    """
    if n_floor <= 0:
        return float(weight)
    return float(weight) * min(1.0, float(n) / float(n_floor))


# ---------------------------------------------------------------------------
# Bond term
# ---------------------------------------------------------------------------
@dataclass
class BondTerm:
    """Mahalanobis loss on a bond length.

    Parameters
    ----------
    a, b : int
        Atom indices defining the bond.
    mu : float
        Library median bond length (Å).
    sigma : float
        Library robust sigma of bond length (Å); must be > 0.
    weight : float, default 1.0
        Per-term weight (already includes any sparse downweight scaling).
    """

    a: int
    b: int
    mu: float
    sigma: float
    weight: float = _W_BOND

    # Indices the gradient touches — used by TotalGripLoss for fast scatter.
    @property
    def atom_indices(self) -> Tuple[int, int]:
        return (int(self.a), int(self.b))

    def value_and_grad(
        self,
        R: np.ndarray,
        grad_out: Optional[np.ndarray] = None,
    ) -> Tuple[float, np.ndarray]:
        """Return ``(loss, grad)`` for this term at coordinates ``R``.

        ``grad_out`` is filled in-place if supplied (shape ``R.shape``,
        accumulator semantics — caller is responsible for zero-init).  Otherwise
        a fresh zero-initialised gradient is allocated.
        """
        R = np.asarray(R, dtype=np.float64)
        if grad_out is None:
            grad = np.zeros_like(R)
            accumulate = False
        else:
            grad = grad_out
            accumulate = True

        ra = R[self.a]
        rb = R[self.b]
        d_vec = ra - rb
        d = float(np.linalg.norm(d_vec))
        sigma = float(self.sigma)
        mu = float(self.mu)
        weight = float(self.weight)

        # Residual in sigma units.
        z = (d - mu) / sigma
        loss = weight * (z * z)

        # dL/d(ra) = 2*w*z * (1/sigma) * d|d|/d(ra) = 2*w*z/sigma * (ra-rb)/|d|
        # dL/d(rb) = -dL/d(ra)
        if d > _EPS:
            coef = 2.0 * weight * z / (sigma * d)
            ga = coef * d_vec
            gb = -ga
        else:
            # Degenerate (coincident atoms) — leave gradient zero for stability.
            ga = np.zeros(3, dtype=np.float64)
            gb = np.zeros(3, dtype=np.float64)

        if accumulate:
            grad[self.a] += ga
            grad[self.b] += gb
        else:
            grad[self.a] = ga
            grad[self.b] = gb

        return loss, grad


# ---------------------------------------------------------------------------
# Angle term
# ---------------------------------------------------------------------------
@dataclass
class AngleTerm:
    """Mahalanobis loss on a valence angle a-b-c (centered at ``b``).

    The loss is on the angle measured in DEGREES (matching the library's
    storage unit).  Gradients are computed in radians internally then
    converted via the chain rule.

    Numerical-stability note (see test ``test_angle_term_gradient_finite_diff``):
    near angles of 0 or pi the derivative of ``acos`` is singular.  We clamp
    the cosine to ``[-1+eps, 1-eps]`` before differentiating and zero the
    gradient on truly degenerate configurations (any zero-length leg).
    """

    a: int
    b: int
    c: int
    mu: float        # degrees
    sigma: float     # degrees
    weight: float = _W_ANGLE

    @property
    def atom_indices(self) -> Tuple[int, int, int]:
        return (int(self.a), int(self.b), int(self.c))

    def value_and_grad(
        self,
        R: np.ndarray,
        grad_out: Optional[np.ndarray] = None,
    ) -> Tuple[float, np.ndarray]:
        R = np.asarray(R, dtype=np.float64)
        if grad_out is None:
            grad = np.zeros_like(R)
            accumulate = False
        else:
            grad = grad_out
            accumulate = True

        ra = R[self.a]
        rb = R[self.b]
        rc = R[self.c]

        u = ra - rb
        v = rc - rb
        nu = float(np.linalg.norm(u))
        nv = float(np.linalg.norm(v))

        if nu < _EPS or nv < _EPS:
            # Degenerate — angle undefined. Skip gradient.
            return 0.0, grad

        # Clamp cos for numerical stability near 0/pi.
        cos_theta = float(np.dot(u, v) / (nu * nv))
        cos_theta = max(-1.0 + 1e-12, min(1.0 - 1e-12, cos_theta))
        sin_theta = math.sqrt(max(1.0 - cos_theta * cos_theta, _EPS))
        theta_rad = math.acos(cos_theta)
        theta_deg = math.degrees(theta_rad)

        sigma = float(self.sigma)
        mu = float(self.mu)
        weight = float(self.weight)
        z = (theta_deg - mu) / sigma
        loss = weight * (z * z)

        # d theta_deg / d theta_rad = 180/pi.
        # dL/d theta_deg = 2*w*z/sigma -> dL/d theta_rad = 2*w*z/sigma * (180/pi).
        # Standard angle-gradient (Bondi/Wilson form):
        #   d cos(theta)/d ra = (v/(nu*nv) - cos(theta)*u/(nu**2)) — careful chain
        # Derive via theta = acos(cos_t), d theta/d cos_t = -1/sin(theta), so
        #   d theta/d ra = -1/sin(theta) * d cos_t/d ra.
        deg_per_rad = 180.0 / math.pi
        dL_dtheta_rad = 2.0 * weight * z / sigma * deg_per_rad

        # Gradients of cos(theta) w.r.t. ra, rc:
        #   d cos_t/d ra = v/(nu*nv) - cos_t * u/(nu*nu)
        #   d cos_t/d rc = u/(nu*nv) - cos_t * v/(nv*nv)
        dcos_dra = v / (nu * nv) - cos_theta * u / (nu * nu)
        dcos_drc = u / (nu * nv) - cos_theta * v / (nv * nv)

        # d theta/d r = -1/sin(theta) * d cos_t / d r
        inv_sin = -1.0 / sin_theta
        dtheta_dra = inv_sin * dcos_dra
        dtheta_drc = inv_sin * dcos_drc
        # Translation invariance: d/d rb = -(d/d ra + d/d rc)
        dtheta_drb = -(dtheta_dra + dtheta_drc)

        ga = dL_dtheta_rad * dtheta_dra
        gb = dL_dtheta_rad * dtheta_drb
        gc = dL_dtheta_rad * dtheta_drc

        if accumulate:
            grad[self.a] += ga
            grad[self.b] += gb
            grad[self.c] += gc
        else:
            grad[self.a] = ga
            grad[self.b] = gb
            grad[self.c] = gc

        return loss, grad


# ---------------------------------------------------------------------------
# Improper term
# ---------------------------------------------------------------------------
@dataclass
class ImproperTerm:
    """Mahalanobis loss on an out-of-plane displacement at a 3-coordinated center.

    Geometry::

        improper d_oop = dot( r_center - centroid(neighbors_3), unit_normal )

    where ``unit_normal`` is the unit normal of the best plane through the 3
    neighbours (their cross-product, normalised).  The sign is *signed* OOP
    distance which preserves chirality information.  The library stores the
    signed value (matching the Mogul "improper" convention used at build
    time), so the loss reproduces it.

    The library stores improper as an angle (degrees, see grip_build_mogul_lib).
    Therefore we convert d_oop to a pseudo-angle via the convention used by
    the build script: improper_angle = signed dihedral C-N1-N2-N3 etc.  For
    Phase 2 we use the SIGNED out-of-plane DISTANCE (Å) as the observable,
    which is mathematically equivalent (up to a known transformation) for
    fixed neighbour-triangle geometry and has a simple analytic gradient.

    The choice of observable (distance vs angle) is documented in design
    notes; the term is mathematically self-consistent (mu/sigma supplied by
    the caller in matching units).
    """

    center: int
    neighbors_3: Tuple[int, int, int]
    mu: float
    sigma: float
    weight: float = _W_IMPROPER

    def __post_init__(self):
        if len(self.neighbors_3) != 3:
            raise ValueError(
                f"ImproperTerm requires exactly 3 neighbours, got {len(self.neighbors_3)}"
            )
        self.neighbors_3 = tuple(int(i) for i in self.neighbors_3)

    @property
    def atom_indices(self) -> Tuple[int, int, int, int]:
        return (int(self.center), *self.neighbors_3)

    def value_and_grad(
        self,
        R: np.ndarray,
        grad_out: Optional[np.ndarray] = None,
    ) -> Tuple[float, np.ndarray]:
        R = np.asarray(R, dtype=np.float64)
        if grad_out is None:
            grad = np.zeros_like(R)
            accumulate = False
        else:
            grad = grad_out
            accumulate = True

        ic = int(self.center)
        i1, i2, i3 = self.neighbors_3
        rc = R[ic]
        r1 = R[i1]
        r2 = R[i2]
        r3 = R[i3]

        # Plane spanned by (r2 - r1) and (r3 - r1); normal n_raw = cross.
        e12 = r2 - r1
        e13 = r3 - r1
        n_raw = np.cross(e12, e13)
        n_norm = float(np.linalg.norm(n_raw))

        if n_norm < _EPS:
            # Collinear neighbours — plane ill-defined. Skip gradient.
            return 0.0, grad

        n_hat = n_raw / n_norm
        centroid = (r1 + r2 + r3) / 3.0
        dvec = rc - centroid
        d_oop = float(np.dot(dvec, n_hat))

        sigma = float(self.sigma)
        mu = float(self.mu)
        weight = float(self.weight)
        z = (d_oop - mu) / sigma
        loss = weight * (z * z)

        dL_dd = 2.0 * weight * z / sigma

        # Gradient of d_oop = dot(rc - centroid, n_hat) w.r.t. each atom.
        # Use d(n_hat)/d(r_i) = (I - n_hat n_hat^T) / |n_raw| * d(n_raw)/d(r_i)
        # where d(n_raw)/d(r_i) is the matrix from cross-product derivatives.
        #
        # For dot(p, n_hat) where p = rc - centroid:
        #   d/d r_i [ dot(p, n_hat) ] = dot(dp/d r_i, n_hat) + dot(p, dn_hat/d r_i)

        # dp/d rc = I, dp/d r1=dp/d r2=dp/d r3 = -I/3
        # d dot(p, n_hat)/d rc directly = n_hat   (n_hat doesn't depend on rc)
        d_dop_drc = n_hat.copy()

        # For the three plane atoms, both p (via centroid) and n_hat change.
        # n_hat = n_raw / |n_raw|; d n_hat / d r_i = (I - n_hat*n_hat^T) / |n_raw| · d n_raw / d r_i
        I3 = np.eye(3, dtype=np.float64)
        proj = (I3 - np.outer(n_hat, n_hat)) / n_norm

        # d n_raw / d r1 = d cross(e12, e13) / d r1
        # e12 = r2 - r1, e13 = r3 - r1
        # n_raw_k = eps_kab e12_a e13_b
        # d n_raw / d r1 = -(d_cross_e12_e13/d e12) - (d_cross_e12_e13/d e13)
        # where d cross(e12, e13)/d e12 represented as matrix [e13]_x acting on right,
        # but since cross(a,b) is linear in both:
        #   d (a x b)/d a = -[b]_x  (skew of b, see e.g. Petersen "Matrix Cookbook")
        #   d (a x b)/d b =  [a]_x
        def skew(v):
            return np.array(
                [
                    [0.0, -v[2], v[1]],
                    [v[2], 0.0, -v[0]],
                    [-v[1], v[0], 0.0],
                ],
                dtype=np.float64,
            )

        cross_de12 = -skew(e13)  # d(e12 x e13)/d e12
        cross_de13 = skew(e12)   # d(e12 x e13)/d e13

        # d n_raw / d r1 = -cross_de12 - cross_de13
        # d n_raw / d r2 =  cross_de12
        # d n_raw / d r3 =  cross_de13
        dn_dr1 = -cross_de12 - cross_de13
        dn_dr2 = cross_de12
        dn_dr3 = cross_de13

        dnhat_dr1 = proj @ dn_dr1
        dnhat_dr2 = proj @ dn_dr2
        dnhat_dr3 = proj @ dn_dr3

        # d p / d r1, r2, r3 = -I/3 each (centroid contribution)
        dp_dri = -I3 / 3.0

        # Compose d dot(p, n_hat) / d r_i = (dp/d r_i)^T n_hat + p^T (dn_hat/d r_i)
        # Here gradient w.r.t. r_i is a 3-vector.
        # (dp/d r_i)^T n_hat = -n_hat/3 for each plane atom.
        # p^T (dn_hat/d r_i) -- left-multiply: vector p @ matrix M gives row vector;
        # the gradient component is p @ dnhat_dri (1x3).
        g1 = (-n_hat / 3.0) + dvec @ dnhat_dr1
        g2 = (-n_hat / 3.0) + dvec @ dnhat_dr2
        g3 = (-n_hat / 3.0) + dvec @ dnhat_dr3

        gc = dL_dd * d_dop_drc
        g1f = dL_dd * g1
        g2f = dL_dd * g2
        g3f = dL_dd * g3

        if accumulate:
            grad[ic] += gc
            grad[i1] += g1f
            grad[i2] += g2f
            grad[i3] += g3f
        else:
            grad[ic] = gc
            grad[i1] = g1f
            grad[i2] = g2f
            grad[i3] = g3f

        return loss, grad


# ---------------------------------------------------------------------------
# Torsion term — Phase 1b stub
# ---------------------------------------------------------------------------
@dataclass
class TorsionTerm:
    """Multi-modal dihedral loss term — NOT IMPLEMENTED in v1.

    The Phase 1 library (``grip_lib_v1.npz``) does not carry torsion
    statistics; those require a CIF-loop extractor that is scheduled for
    Phase 1b (post-CCDC-API).  Construction is allowed (for API symmetry
    in builders), but :meth:`value_and_grad` raises immediately so a caller
    that wires a torsion accidentally fails loud rather than producing a
    silently wrong loss.
    """

    a: int
    b: int
    c: int
    d: int
    gmm_params: object = None   # placeholder for (pi_k, mu_k, sigma_k) tuples
    weight: float = _W_TORSION

    @property
    def atom_indices(self) -> Tuple[int, int, int, int]:
        return (int(self.a), int(self.b), int(self.c), int(self.d))

    def value_and_grad(self, R: np.ndarray, grad_out: Optional[np.ndarray] = None):
        raise NotImplementedError(
            "TorsionTerm requires the Phase 1b torsion library (grip_lib_v2). "
            "v1 has no torsion stats; use BondTerm/AngleTerm/ImproperTerm only."
        )


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------
@dataclass
class TotalGripLoss:
    """Aggregates a deterministic, sorted list of fragment terms.

    Terms are sorted by the lexicographic order of their ``atom_indices``
    tuple on construction so the evaluation order is reproducible across
    runs and platforms.  The aggregator's :meth:`evaluate` returns the
    summed scalar loss and the summed gradient (shape ``(n_atoms, 3)``).
    """

    terms: List[object] = field(default_factory=list)

    def __post_init__(self):
        # Sort for determinism. Sort key is the atom_indices tuple — this
        # gives a stable, content-derived order independent of insertion
        # order or Python's dict-hash randomisation.
        self.terms = sorted(self.terms, key=lambda t: tuple(t.atom_indices))

    def add(self, term) -> None:
        """Append a term and re-sort (rare, prefer constructing once)."""
        self.terms.append(term)
        self.terms.sort(key=lambda t: tuple(t.atom_indices))

    def __len__(self) -> int:
        return len(self.terms)

    def evaluate(self, R: np.ndarray) -> Tuple[float, np.ndarray]:
        """Compute the total loss and gradient at ``R``.

        ``R`` may be a flat ``(3N,)`` or shaped ``(N, 3)`` array; the
        returned gradient matches the input shape.
        """
        R_in = np.asarray(R, dtype=np.float64)
        flat = R_in.ndim == 1
        if flat:
            if R_in.size % 3 != 0:
                raise ValueError(
                    f"R must be (N,3) or (3N,), got length {R_in.size}"
                )
            R_arr = R_in.reshape(-1, 3)
        else:
            R_arr = R_in

        grad = np.zeros_like(R_arr)
        total = 0.0
        # Pre-sorted: walk in the canonical order.
        for term in self.terms:
            l, _ = term.value_and_grad(R_arr, grad_out=grad)
            total += float(l)

        if flat:
            return float(total), grad.reshape(-1)
        return float(total), grad
