"""GRIP — M-D unfreeze loss terms (Phase 2, 2026-06-06).

Phase 2 of the GRIP M+D unfreeze adds *coordination-sphere* loss terms so the
metal + donor atoms can move under L-BFGS while still being held to their
CCDC-empirical / polyhedron-template ideals.

Three families of terms ship here:

* :class:`MDBondTerm` — Mahalanobis Gaussian on the M-D bond length around
  the covalent-radii sum (or a CCDC-pooled median when a v5 library is
  available).  σ = 0.05 Å by default, weight = 10.0 — high enough to pull
  M-D back to the ideal but not so high that it blows up the loss when an
  initial guess is already 0.3 Å off.

* :class:`MDXAngleTerm` — Mahalanobis Gaussian on the M-D-X angle where X is
  a heavy neighbour of donor D in the molecular graph.  The target angle
  comes from the donor's coordination ideal (sp³ → 109.5°, sp² → 120°,
  sp linear → 180°), giving a deterministic, FF-free angular pressure on
  the M-D vector through the donor's neighbour set.  σ = 5°, weight = 3.0.

* :class:`DMDAngleTerm` — Mahalanobis Gaussian on the D-M-D inter-donor
  angle, computed from the polyhedron reference vectors of the geometry
  (`polyhedra.ref_vectors`).  σ = 5°, weight = 5.0.  This is the term that
  pulls the donor set onto the ideal octahedron / TBP / etc. — the
  polyhedron CShM is currently only used as a rollback validator; this
  term moves it into the gradient.

Hapto branch:
* :class:`HaptoCentroidTerm` — for a π-coordinated ring (3+ same-element
  heavy donors in the same ring within 3.0 Å of the metal), the loss is
  on the *centroid-to-metal* distance (not the individual M-C distances)
  with σ = 0.05 Å and the configured M-π distance.

Determinism contract (matches :mod:`grip_loss_terms`):

* float64 everywhere
* no RNG, no hash-dependent iteration
* terms expose ``atom_indices`` so :class:`TotalGripLoss` can sort them
* analytical gradients validated against central finite differences in the
  Phase 2 test suite

FF-FREE invariant: every weight + sigma here is a Gaussian *prior*
(empirical scatter from CCDC), not a parameterised spring constant.  No
Lennard-Jones, no Coulomb, no UFF / MMFF parameters are referenced.

Robustness: every public builder is wrapped in try / except — when a target
cannot be computed the term is simply not emitted (fail open).  Callers must
never rely on every term being present; they must always run with whatever
:func:`build_md_loss_terms` returns.
"""
from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

_LOG = logging.getLogger(__name__)

# Numerical floor — match grip_loss_terms.
_EPS = 1e-12

# Default Gaussian widths (sigma) -- the values used when no CCDC sample
# resolves to a finer width.  Tuned conservatively so the gradient is
# strong enough to pull but cannot overwhelm a tight polyhedron CShM.
DEFAULT_MD_SIGMA: float = 0.05    # Å — matches the M-D-invariant tol
DEFAULT_MDX_SIGMA: float = 5.0    # degrees
DEFAULT_DMD_SIGMA: float = 5.0    # degrees

# Default per-class weights.  Resolution: explicit kwarg > env > default.
DEFAULT_MD_WEIGHT: float = 10.0
DEFAULT_MDX_WEIGHT: float = 3.0
DEFAULT_DMD_WEIGHT: float = 5.0

_MD_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_MD_WEIGHT"
_MDX_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_MDX_WEIGHT"
_DMD_WEIGHT_ENV: str = "DELFIN_FFFREE_GRIP_DMD_WEIGHT"

# Hapto-cluster detection: a donor counts as hapto when the metal has
# 3+ same-element heavy neighbours of that element within HAPTO_RADIUS.
HAPTO_RADIUS: float = 3.0   # Å (matches the spec)
HAPTO_MIN_COUNT: int = 3


# ---------------------------------------------------------------------------
# Env-flag resolvers
# ---------------------------------------------------------------------------
_UNFREEZE_ENV: str = "DELFIN_FFFREE_GRIP_UNFREEZE_MD"


def unfreeze_md_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_GRIP_UNFREEZE_MD`` is on (default OFF).

    Phase 2 master switch.  When unset / 0 the polish remains byte-identical
    with the legacy frozen-sphere path.  Resolve once per call so multi-proc
    pools that pass env through fork() see the flag.
    """
    raw = os.environ.get(_UNFREEZE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _resolve_pos_float(env_name: str, default: float) -> float:
    """Read a positive finite float from env; fall back to ``default``."""
    raw = os.environ.get(env_name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if math.isfinite(v) and v > 0.0:
            return v
    except (TypeError, ValueError):
        pass
    return float(default)


def _resolve_md_weight() -> float:
    return _resolve_pos_float(_MD_WEIGHT_ENV, DEFAULT_MD_WEIGHT)


def _resolve_mdx_weight() -> float:
    return _resolve_pos_float(_MDX_WEIGHT_ENV, DEFAULT_MDX_WEIGHT)


def _resolve_dmd_weight() -> float:
    return _resolve_pos_float(_DMD_WEIGHT_ENV, DEFAULT_DMD_WEIGHT)


# ---------------------------------------------------------------------------
# Helpers: symbols + neighbours
# ---------------------------------------------------------------------------
def _symbols(mol) -> List[str]:
    """Return per-atom element symbols (deterministic, indexed)."""
    try:
        return [str(a.GetSymbol()) for a in mol.GetAtoms()]
    except Exception:
        return []


def _heavy_neighbours(mol, idx: int) -> List[int]:
    """Sorted heavy-neighbour indices of atom ``idx`` in ``mol``."""
    out: List[int] = []
    try:
        atom = mol.GetAtomWithIdx(int(idx))
    except Exception:
        return out
    try:
        for nb in atom.GetNeighbors():
            try:
                if nb.GetSymbol() != "H":
                    out.append(int(nb.GetIdx()))
            except Exception:
                continue
    except Exception:
        return out
    return sorted(set(out))


# Covalent-radii table (subset of :data:`polyhedra.COV` -- duplicated here
# so this module has no circular import on polyhedra at load time).
_COV_RADII: Dict[str, float] = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39, "Se": 1.20, "As": 1.19,
    "B": 0.84, "Si": 1.11, "Te": 1.38,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "La": 2.07, "Lu": 1.87,
    "Sn": 1.39, "Pb": 1.46, "Bi": 1.48, "Sb": 1.39,
}


def _cov_sum(z1: str, z2: str) -> Optional[float]:
    """Return covalent-radii sum for two elements; ``None`` if missing."""
    r1 = _COV_RADII.get(str(z1))
    r2 = _COV_RADII.get(str(z2))
    if r1 is None or r2 is None:
        return None
    return float(r1 + r2)


def _md_target_from_library(
    library,
    z_metal: str,
    z_donor: str,
) -> Optional[Tuple[float, float, int]]:
    """Best-effort CCDC-empirical M-D bond ``(mu, sigma, n)``.

    The library API (:class:`grip_mogul_lookup.GripLibrary`) takes
    hybridisation strings which we do not have at the M-D level, so we use
    wildcards.  Robust: any exception -> ``None`` (caller falls back to
    cov-sum + DEFAULT_MD_SIGMA).
    """
    if library is None:
        return None
    try:
        hit = library.lookup_bond(
            str(z_metal), "*", str(z_donor), "*",
            ring_size_min=-1, in_aromatic=False,
        )
    except Exception:
        return None
    if hit is None:
        return None
    try:
        mu, sigma, n = hit
        mu_f = float(mu)
        sg_f = float(sigma)
        n_i = int(n)
    except Exception:
        return None
    if not (math.isfinite(mu_f) and math.isfinite(sg_f)):
        return None
    if mu_f <= 0.0 or sg_f <= 0.0:
        return None
    return (mu_f, sg_f, n_i)


def _md_target(
    library,
    z_metal: str,
    z_donor: str,
    init_distance: float,
) -> Tuple[float, float]:
    """Resolve ``(mu, sigma)`` for an M-D Gaussian.

    Resolution order:
      1. CCDC library median if available (and finite + positive).
      2. Covalent-radii sum (polyhedra.COV-equivalent).
      3. Initial distance ``init_distance`` (last resort; sigma defaults).
    """
    hit = _md_target_from_library(library, z_metal, z_donor)
    if hit is not None:
        return (hit[0], max(hit[1], DEFAULT_MD_SIGMA * 0.5))
    cov = _cov_sum(z_metal, z_donor)
    if cov is not None and 0.5 < cov < 5.0:
        return (cov, DEFAULT_MD_SIGMA)
    if math.isfinite(init_distance) and init_distance > 0.0:
        return (float(init_distance), DEFAULT_MD_SIGMA)
    return (2.0, DEFAULT_MD_SIGMA)


# ---------------------------------------------------------------------------
# Hapto detection (per-mission spec: 3+ same-element heavy neighbours of M
# within HAPTO_RADIUS).  This is *position-independent* heuristic for the
# build step: we use the molecular graph + the initial M-D distances on
# ``P_init`` to decide whether a donor is part of a hapto cluster.
# ---------------------------------------------------------------------------
def detect_hapto_donor_clusters(
    mol,
    metal_idx: int,
    donors: Sequence[int],
    P_init: np.ndarray,
    radius: float = HAPTO_RADIUS,
    min_count: int = HAPTO_MIN_COUNT,
) -> List[Tuple[str, Tuple[int, ...]]]:
    """Return ``[(element, (donor_indices,)), ...]`` for each hapto cluster.

    A cluster is a set of donor atoms with the SAME element whose count is
    ``>= min_count`` AND whose initial distance to the metal is ``<= radius``.
    Pure-functional; sorted output for determinism.

    The cluster element is the shared element symbol (e.g. ``"C"`` for a Cp
    ring).  Returned tuples of donor indices are sorted.
    """
    out: List[Tuple[str, Tuple[int, ...]]] = []
    syms = _symbols(mol)
    if not syms:
        return out
    try:
        metal_pos = np.asarray(P_init[int(metal_idx)], dtype=np.float64)
    except Exception:
        return out

    # Group donors by element.
    by_elem: Dict[str, List[int]] = {}
    for d in donors:
        try:
            di = int(d)
        except Exception:
            continue
        if di < 0 or di >= len(syms):
            continue
        try:
            dist = float(np.linalg.norm(P_init[di] - metal_pos))
        except Exception:
            continue
        if not math.isfinite(dist) or dist > float(radius):
            continue
        elem = str(syms[di])
        if elem == "H":  # H atoms cannot hapto-coordinate
            continue
        by_elem.setdefault(elem, []).append(di)

    for elem in sorted(by_elem.keys()):
        group = sorted(set(by_elem[elem]))
        if len(group) >= int(min_count):
            out.append((elem, tuple(group)))
    return sorted(out)


# ---------------------------------------------------------------------------
# Loss terms
# ---------------------------------------------------------------------------
@dataclass
class MDBondTerm:
    """Mahalanobis Gaussian on the metal-donor distance.

    Parameters mirror :class:`grip_loss_terms.BondTerm` so the term can be
    added to the same :class:`TotalGripLoss` aggregator.  The gradient is
    analytical (same algebra as ``BondTerm``).
    """

    metal: int
    donor: int
    mu: float
    sigma: float
    weight: float = DEFAULT_MD_WEIGHT

    @property
    def atom_indices(self) -> Tuple[int, int]:
        # Sorted for deterministic agg-order — the BondTerm convention.
        a, b = int(self.metal), int(self.donor)
        return (min(a, b), max(a, b))

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
        mi, di = int(self.metal), int(self.donor)
        rm = R[mi]
        rd = R[di]
        d_vec = rd - rm
        d = float(np.linalg.norm(d_vec))
        sigma = float(self.sigma)
        if sigma <= 0.0 or not math.isfinite(d):
            return 0.0, grad
        z = (d - float(self.mu)) / sigma
        w = float(self.weight)
        loss = w * (z * z)
        if d > _EPS:
            coef = 2.0 * w * z / (sigma * d)
            gd = coef * d_vec
            gm = -gd
        else:
            gd = np.zeros(3, dtype=np.float64)
            gm = np.zeros(3, dtype=np.float64)
        if accumulate:
            grad[mi] += gm
            grad[di] += gd
        else:
            grad[mi] = gm
            grad[di] = gd
        return loss, grad


@dataclass
class MDXAngleTerm:
    """Gaussian on the M-D-X angle (centred at donor).

    Same algebra as :class:`grip_loss_terms.AngleTerm`; the only difference
    is intent — the target ``mu`` here comes from the donor's *coordination*
    ideal (sp³ → 109.5°, sp² → 120°, sp → 180°) rather than from a CCDC
    histogram of free-ligand angles.
    """

    metal: int
    donor: int
    x: int
    mu: float        # degrees
    sigma: float     # degrees
    weight: float = DEFAULT_MDX_WEIGHT

    @property
    def atom_indices(self) -> Tuple[int, int, int]:
        return (int(self.metal), int(self.donor), int(self.x))

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
        mi, di, xi = int(self.metal), int(self.donor), int(self.x)
        rm, rd, rx = R[mi], R[di], R[xi]
        u = rm - rd
        v = rx - rd
        un = float(np.linalg.norm(u))
        vn = float(np.linalg.norm(v))
        if un < _EPS or vn < _EPS:
            return 0.0, grad
        cos_t = float(np.dot(u, v)) / (un * vn)
        # Numerical clamp matches AngleTerm.
        cos_clamped = max(-1.0 + 1e-12, min(1.0 - 1e-12, cos_t))
        theta = math.acos(cos_clamped)
        theta_deg = math.degrees(theta)
        sigma = float(self.sigma)
        if sigma <= 0.0:
            return 0.0, grad
        z = (theta_deg - float(self.mu)) / sigma
        w = float(self.weight)
        loss = w * (z * z)
        # dL/dtheta_deg = 2 w z / sigma
        dL_dtheta_deg = 2.0 * w * z / sigma
        # dtheta_deg/dtheta_rad = 180/pi
        # dtheta_rad/dcos_t = -1 / sin(theta)
        sin_t = math.sqrt(max(1.0 - cos_clamped * cos_clamped, 1e-24))
        d_theta_d_cos = -1.0 / sin_t
        # dcos_t / dr_m, dr_d, dr_x.
        # cos_t = (u . v) / (|u| |v|)  with u = rm - rd, v = rx - rd.
        u_hat = u / un
        v_hat = v / vn
        # d cos_t / dr_m =  (v_hat - cos_t * u_hat) / |u|
        # d cos_t / dr_x =  (u_hat - cos_t * v_hat) / |v|
        # d cos_t / dr_d = -(d cos/dr_m + d cos/dr_x)
        dcos_drm = (v_hat - cos_clamped * u_hat) / un
        dcos_drx = (u_hat - cos_clamped * v_hat) / vn
        dcos_drd = -(dcos_drm + dcos_drx)
        chain = dL_dtheta_deg * (180.0 / math.pi) * d_theta_d_cos
        gm = chain * dcos_drm
        gd = chain * dcos_drd
        gx = chain * dcos_drx
        if accumulate:
            grad[mi] += gm
            grad[di] += gd
            grad[xi] += gx
        else:
            grad[mi] = gm
            grad[di] = gd
            grad[xi] = gx
        return loss, grad


@dataclass
class DMDAngleTerm:
    """Gaussian on the D-M-D inter-donor angle (centred at metal).

    Target ``mu`` comes from the polyhedron reference vectors
    (``polyhedra.ref_vectors(geom)``) for the assigned geometry: the angle
    between the i-th and j-th reference vectors in degrees.  When the
    polyhedron is unknown the term is not emitted by :func:`build_md_loss_terms`.
    """

    metal: int
    donor_i: int
    donor_j: int
    mu: float        # degrees
    sigma: float     # degrees
    weight: float = DEFAULT_DMD_WEIGHT

    @property
    def atom_indices(self) -> Tuple[int, int, int]:
        a, b = int(self.donor_i), int(self.donor_j)
        if a > b:
            a, b = b, a
        return (a, int(self.metal), b)

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
        mi = int(self.metal)
        di = int(self.donor_i)
        dj = int(self.donor_j)
        rm, ri, rj = R[mi], R[di], R[dj]
        u = ri - rm
        v = rj - rm
        un = float(np.linalg.norm(u))
        vn = float(np.linalg.norm(v))
        if un < _EPS or vn < _EPS:
            return 0.0, grad
        cos_t = float(np.dot(u, v)) / (un * vn)
        cos_clamped = max(-1.0 + 1e-12, min(1.0 - 1e-12, cos_t))
        theta = math.acos(cos_clamped)
        theta_deg = math.degrees(theta)
        sigma = float(self.sigma)
        if sigma <= 0.0:
            return 0.0, grad
        z = (theta_deg - float(self.mu)) / sigma
        w = float(self.weight)
        loss = w * (z * z)
        dL_dtheta_deg = 2.0 * w * z / sigma
        sin_t = math.sqrt(max(1.0 - cos_clamped * cos_clamped, 1e-24))
        d_theta_d_cos = -1.0 / sin_t
        u_hat = u / un
        v_hat = v / vn
        dcos_dri = (v_hat - cos_clamped * u_hat) / un
        dcos_drj = (u_hat - cos_clamped * v_hat) / vn
        dcos_drm = -(dcos_dri + dcos_drj)
        chain = dL_dtheta_deg * (180.0 / math.pi) * d_theta_d_cos
        gm = chain * dcos_drm
        gi = chain * dcos_dri
        gj = chain * dcos_drj
        if accumulate:
            grad[mi] += gm
            grad[di] += gi
            grad[dj] += gj
        else:
            grad[mi] = gm
            grad[di] = gi
            grad[dj] = gj
        return loss, grad


@dataclass
class HaptoCentroidTerm:
    """Gaussian on the centroid-to-metal distance of a hapto-π ring.

    For a hapto ring of donor atoms ``ring_atoms``, the centroid is the
    arithmetic mean of their positions; the loss is a Gaussian on
    ``|M - centroid|`` around the target distance ``mu``.  This is the
    correct shape for π-coordination — it does not pull each individual
    M-C distance to an idealised value (which would over-constrain the
    geometry), only the metal-to-ring-plane distance.

    Gradient is analytical: the centroid is linear in the ring atoms so
    the gradient w.r.t. each ring atom is ``-1/k`` of the metal gradient.
    """

    metal: int
    ring_atoms: Tuple[int, ...]
    mu: float
    sigma: float
    weight: float = DEFAULT_MD_WEIGHT

    @property
    def atom_indices(self) -> Tuple[int, ...]:
        return (int(self.metal),) + tuple(sorted(int(a) for a in self.ring_atoms))

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
        ra = [int(a) for a in self.ring_atoms]
        k = len(ra)
        if k == 0:
            return 0.0, grad
        mi = int(self.metal)
        centroid = np.mean(R[ra], axis=0)
        d_vec = centroid - R[mi]
        d = float(np.linalg.norm(d_vec))
        sigma = float(self.sigma)
        if sigma <= 0.0 or not math.isfinite(d):
            return 0.0, grad
        z = (d - float(self.mu)) / sigma
        w = float(self.weight)
        loss = w * (z * z)
        if d > _EPS:
            coef = 2.0 * w * z / (sigma * d)
            # dL/d centroid = coef * (centroid - r_m) = coef * d_vec
            g_centroid = coef * d_vec
            g_m = -g_centroid
            g_ring_each = g_centroid / float(k)
        else:
            g_m = np.zeros(3, dtype=np.float64)
            g_ring_each = np.zeros(3, dtype=np.float64)
        if accumulate:
            grad[mi] += g_m
            for a in ra:
                grad[a] += g_ring_each
        else:
            grad[mi] = g_m
            for a in ra:
                grad[a] = g_ring_each
        return loss, grad


# ---------------------------------------------------------------------------
# Donor coordination ideal: target X-D-M angle based on donor hybridisation
# ---------------------------------------------------------------------------
def _donor_coordination_angle_deg(mol, donor: int) -> Optional[float]:
    """Best-effort ideal M-D-X angle (degrees) from donor hybridisation.

    Resolution: try RDKit ``GetHybridization`` first; fall back to a
    heavy-degree heuristic (deg 3 -> 109.5°, deg 2 -> 120°, deg 1 -> 180°).
    Returns ``None`` if the donor cannot be resolved.
    """
    try:
        atom = mol.GetAtomWithIdx(int(donor))
    except Exception:
        return None
    try:
        h = atom.GetHybridization()
        h_str = str(h).upper()
    except Exception:
        h_str = ""
    if "SP3" in h_str:
        return 109.471  # Td
    if "SP2" in h_str:
        return 120.0
    if h_str == "SP" or "SP_" in h_str:
        return 180.0
    # Fallback: heavy-neighbour degree.
    nbrs = _heavy_neighbours(mol, donor)
    deg = len(nbrs)
    if deg >= 3:
        return 109.471
    if deg == 2:
        return 120.0
    if deg == 1:
        return 180.0
    return None


# ---------------------------------------------------------------------------
# Polyhedron D-M-D angle lookup
# ---------------------------------------------------------------------------
def _polyhedron_dmd_angles(
    geom: str,
    n_donors: int,
) -> Optional[np.ndarray]:
    """Return the ``(n, n)`` matrix of pairwise polyhedron-vertex angles
    (degrees) for ``geom``, or ``None`` when unavailable.

    The matrix is indexed by donor SLOT (the same order the caller used to
    place the donors on the polyhedron vertices).  ``mat[i,j]`` is the
    angle ``deg(v_i, v_j)``; the diagonal is zero.
    """
    try:
        from .polyhedra import ref_vectors
    except Exception:
        return None
    try:
        v = ref_vectors(str(geom))
    except Exception:
        return None
    if v is None:
        return None
    v = np.asarray(v, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 3:
        return None
    if v.shape[0] < n_donors:
        return None
    # Use only the first n_donors rows (caller assigns donors to slots
    # in that order).
    sub = v[:n_donors]
    norms = np.linalg.norm(sub, axis=1, keepdims=True)
    norms = np.where(norms < _EPS, 1.0, norms)
    unit = sub / norms
    cos_mat = unit @ unit.T
    cos_mat = np.clip(cos_mat, -1.0, 1.0)
    return np.degrees(np.arccos(cos_mat))


# ---------------------------------------------------------------------------
# Public builder
# ---------------------------------------------------------------------------
def build_md_loss_terms(
    mol,
    P_init: np.ndarray,
    metal: int,
    donors: Sequence[int],
    geom: str = "",
    *,
    library=None,
    md_weight: Optional[float] = None,
    mdx_weight: Optional[float] = None,
    dmd_weight: Optional[float] = None,
    md_sigma: float = DEFAULT_MD_SIGMA,
    mdx_sigma: float = DEFAULT_MDX_SIGMA,
    dmd_sigma: float = DEFAULT_DMD_SIGMA,
    include_hapto: bool = True,
) -> List[object]:
    """Build the list of Phase 2 M+D loss terms.

    Returns a (possibly empty) list of terms.  Each is added to the existing
    :class:`TotalGripLoss` aggregator — Phase 2 terms compose with the
    Phase 1 (bond / angle / improper) terms cleanly because they all share
    the same ``value_and_grad`` interface.

    The function is robust by construction: every per-donor / per-pair
    sub-step is wrapped in try / except.  When any single target cannot be
    computed the corresponding term is simply omitted (fail open).

    Parameters
    ----------
    mol : RDKit Mol-like
        Molecular graph.
    P_init : ndarray (N, 3)
        Initial coordinates; used to compute fall-back distances when the
        library / cov-sum lookups fail.
    metal : int
        Metal atom index.
    donors : sequence of int
        Donor atom indices (in the same slot order used by the polyhedron
        placement).
    geom : str
        Coordination geometry name (e.g. ``"OC-6"``).  Empty disables the
        D-M-D angle term.
    library : GripLibrary, optional
        CCDC-empirical library.  When provided, M-D targets prefer the
        library median over cov-sum.
    md_weight, mdx_weight, dmd_weight : float, optional
        Override the per-class weights.  ``None`` -> env -> default.
    md_sigma, mdx_sigma, dmd_sigma : float
        Gaussian widths.  Defaults match SPEC.
    include_hapto : bool
        If True, hapto donor clusters are routed to a single
        :class:`HaptoCentroidTerm` instead of per-donor :class:`MDBondTerm`.
    """
    terms: List[object] = []

    P_init = np.asarray(P_init, dtype=np.float64)
    if P_init.ndim != 2 or P_init.shape[1] != 3:
        return terms
    syms = _symbols(mol)
    if not syms:
        return terms
    n_atoms = len(syms)
    metal = int(metal)
    if metal < 0 or metal >= n_atoms:
        return terms
    donor_list = [int(d) for d in donors if 0 <= int(d) < n_atoms and int(d) != metal]
    if not donor_list:
        return terms

    z_metal = str(syms[metal])
    w_md = float(_resolve_md_weight() if md_weight is None else md_weight)
    w_mdx = float(_resolve_mdx_weight() if mdx_weight is None else mdx_weight)
    w_dmd = float(_resolve_dmd_weight() if dmd_weight is None else dmd_weight)

    # ------------------------------------------------------------------
    # Hapto-cluster detection: collapse hapto donor sets to a single
    # centroid term and exclude them from per-donor MDBondTerm.
    # ------------------------------------------------------------------
    hapto_clusters: List[Tuple[str, Tuple[int, ...]]] = []
    hapto_donors_excluded: Set[int] = set()
    if include_hapto:
        try:
            hapto_clusters = detect_hapto_donor_clusters(
                mol, metal, donor_list, P_init,
            )
        except Exception:
            hapto_clusters = []
        for _elem, ring in hapto_clusters:
            for a in ring:
                hapto_donors_excluded.add(int(a))

    for elem, ring in hapto_clusters:
        try:
            ring_pos = P_init[list(ring)]
            centroid = np.mean(ring_pos, axis=0)
            init_d = float(np.linalg.norm(centroid - P_init[metal]))
            if not math.isfinite(init_d) or init_d <= 0.0:
                continue
            # Target: cov-sum + a small radial offset (typical η-bond ~ 2.0 Å)
            cov = _cov_sum(z_metal, elem)
            mu_target = cov if cov is not None else init_d
            terms.append(HaptoCentroidTerm(
                metal=metal,
                ring_atoms=tuple(int(a) for a in ring),
                mu=float(mu_target),
                sigma=float(md_sigma),
                weight=w_md,
            ))
        except Exception:
            continue

    # ------------------------------------------------------------------
    # Per-donor M-D bond Gaussian.
    # ------------------------------------------------------------------
    for d in donor_list:
        if d in hapto_donors_excluded:
            continue
        try:
            z_d = str(syms[d])
            init_d = float(np.linalg.norm(P_init[d] - P_init[metal]))
            mu, sig_lib = _md_target(library, z_metal, z_d, init_d)
            sigma_eff = max(float(md_sigma), sig_lib)
            terms.append(MDBondTerm(
                metal=metal,
                donor=d,
                mu=float(mu),
                sigma=float(sigma_eff),
                weight=w_md,
            ))
        except Exception:
            continue

    # ------------------------------------------------------------------
    # M-D-X angle Gaussians (per donor, per heavy neighbour X).  Skip
    # hapto donors (their angular environment is described by the
    # centroid term, not by individual M-D-X angles).
    # ------------------------------------------------------------------
    for d in donor_list:
        if d in hapto_donors_excluded:
            continue
        try:
            ideal_deg = _donor_coordination_angle_deg(mol, d)
            if ideal_deg is None:
                continue
            xs = _heavy_neighbours(mol, d)
            # Exclude the metal itself from the neighbour list (the M-D-X
            # angle is for X != M).
            xs = [x for x in xs if x != metal]
            for x in xs:
                terms.append(MDXAngleTerm(
                    metal=metal,
                    donor=d,
                    x=int(x),
                    mu=float(ideal_deg),
                    sigma=float(mdx_sigma),
                    weight=w_mdx,
                ))
        except Exception:
            continue

    # ------------------------------------------------------------------
    # D-M-D angle Gaussians from polyhedron template.
    #
    # We match donors to polyhedron vertices in the SAME order the caller
    # supplied them — this is the convention assemble_complex.py uses when
    # placing donors on the polyhedron.  Hapto donors are excluded from
    # the pairing (a hapto ring has no D-M-D angle to enforce against
    # another donor — the ring centroid is the relevant axis).
    # ------------------------------------------------------------------
    if geom:
        try:
            # The polyhedron ref_vectors has slots equal to ``len(donor_list)``
            # in the standard case.  We pass the full donor list (in caller
            # order) so the slot mapping is the same one the placement uses.
            non_hapto_donors = [d for d in donor_list if d not in hapto_donors_excluded]
            if len(non_hapto_donors) >= 2:
                # Recompute the polyhedron matrix on the FULL polyhedron
                # then slice by the non-hapto subset's *index in donor_list*
                # (the slot order).
                slot_of: Dict[int, int] = {d: i for i, d in enumerate(donor_list)}
                mat = _polyhedron_dmd_angles(geom, len(donor_list))
                if mat is not None:
                    for ai in range(len(non_hapto_donors)):
                        for aj in range(ai + 1, len(non_hapto_donors)):
                            da = non_hapto_donors[ai]
                            db = non_hapto_donors[aj]
                            si = slot_of[da]
                            sj = slot_of[db]
                            mu = float(mat[si, sj])
                            if not math.isfinite(mu):
                                continue
                            terms.append(DMDAngleTerm(
                                metal=metal,
                                donor_i=da,
                                donor_j=db,
                                mu=mu,
                                sigma=float(dmd_sigma),
                                weight=w_dmd,
                            ))
        except Exception:
            pass

    return terms


# ---------------------------------------------------------------------------
# Diagnostic: max M-D drift after polish (used by the fail-open guard in
# grip_polish so a Phase 2 polish that moves M-D >0.5 Å reverts to frozen).
# ---------------------------------------------------------------------------
def max_md_drift(
    P_before: np.ndarray,
    P_after: np.ndarray,
    metal: int,
    donors: Sequence[int],
) -> float:
    """Return ``max | |M-D|_after - |M-D|_before |`` over all donors (Å).

    Used by the fail-open Phase 2 guard: when this value exceeds the safety
    threshold (0.5 Å) the polish is reverted to the frozen-sphere result.
    """
    P_before = np.asarray(P_before, dtype=np.float64)
    P_after = np.asarray(P_after, dtype=np.float64)
    try:
        mi = int(metal)
        worst = 0.0
        for d in donors:
            di = int(d)
            db = float(np.linalg.norm(P_before[di] - P_before[mi]))
            da = float(np.linalg.norm(P_after[di] - P_after[mi]))
            if not (math.isfinite(db) and math.isfinite(da)):
                return float("inf")
            diff = abs(da - db)
            if diff > worst:
                worst = diff
        return float(worst)
    except Exception:
        return float("inf")


__all__ = [
    "unfreeze_md_active",
    "DEFAULT_MD_SIGMA",
    "DEFAULT_MDX_SIGMA",
    "DEFAULT_DMD_SIGMA",
    "DEFAULT_MD_WEIGHT",
    "DEFAULT_MDX_WEIGHT",
    "DEFAULT_DMD_WEIGHT",
    "HAPTO_RADIUS",
    "HAPTO_MIN_COUNT",
    "MDBondTerm",
    "MDXAngleTerm",
    "DMDAngleTerm",
    "HaptoCentroidTerm",
    "build_md_loss_terms",
    "detect_hapto_donor_clusters",
    "max_md_drift",
]
