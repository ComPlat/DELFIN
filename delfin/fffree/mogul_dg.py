"""Mogul-DG Phase D — Public ``mogul_embed`` driver.

This module ties Phases A (bounds), B (solver), and C (severity) together
into the public API surface that the converter dispatch and the embed
fallback consume:

* :func:`mogul_embed` — accepts a SMILES (or pre-decomposed mol + metal +
  donors) and returns ``(syms, P, info)`` or ``None`` on infeasibility.
* :func:`mogul_embed_from_assembled` — takes a pre-built ``(syms, P_init,
  mol, metal_idx, donor_idxs)`` topology and refines it via the same
  whole-complex DG machinery.  Used by the converter path to replace
  rigid-fit with the Mogul-informed solver.
* :func:`is_active` — single source of truth for "should we activate this
  path?"  Default OFF byte-identical (see SPEC §3.3).
* :func:`should_replace_etkdg`, :func:`should_replace_rigid` — finer
  granularity flags for the two integration points.

Determinism contract
--------------------
* All randomness seeded by ``DELFIN_FFFREE_MOGUL_DG_SEED`` (default 42).
* ``PYTHONHASHSEED=0`` is honoured by Phase A/B/C.  No additional global
  state introduced here.
* Two runs with the same SMILES + same env → byte-identical output
  (verified by ``tests/test_mogul_dg_integration.py``).

Fail-open contract
------------------
* Every internal exception is swallowed and returned as ``None``.
* The caller (converter_backend or embed_fallback) treats ``None`` as
  "fall back to the legacy path" — production stays safe even when the
  bounds/solver/severity stacks are misconfigured.

Universality contract
---------------------
* No SMILES-specific code paths anywhere in this module.
* No element-name branches: every chemistry decision flows through the
  graph + CCDC library lookup downstream.

Public env flags (default unset = OFF, byte-identical to current pipeline)
-------------------------------------------------------------------------
``DELFIN_FFFREE_MOGUL_DG``                     master toggle
``DELFIN_FFFREE_MOGUL_DG_REPLACE_ETKDG``       use mogul_dg INSTEAD of ETKDG in F2
``DELFIN_FFFREE_MOGUL_DG_REPLACE_RIGID``       use mogul_dg INSTEAD of rigid-fit
``DELFIN_FFFREE_MOGUL_DG_BOUNDS_TOL``          σ-mult for bond bounds (default 3.0)
``DELFIN_FFFREE_MOGUL_DG_MD_BOUNDS_TOL``       σ-mult for M-D bounds (default 2.0)
``DELFIN_FFFREE_MOGUL_DG_MAX_ITER``            L-BFGS iteration cap (default 1000)
``DELFIN_FFFREE_MOGUL_DG_SEED``                determinism seed (default 42)
``DELFIN_FFFREE_MOGUL_DG_N_RESTARTS``          multi-restart count (default 3)
``DELFIN_FFFREE_MOGUL_DG_MD_DRIFT_TOL``        M-D drift guard Å (default 0.5)
``DELFIN_FFFREE_MOGUL_DG_DONOR_SHELL_VERIFY``  enable first-shell integrity verification
``DELFIN_FFFREE_MOGUL_DG_REJECT_H_FIRST_SHELL`` push spurious first-shell H atoms outward
``DELFIN_FFFREE_MOGUL_DG_DONOR_REPAIR_RETRIES`` repair retries for missing donors (default 3)
``DELFIN_FFFREE_MOGUL_DG_MD_TOL``              Å tolerance vs ideal M-D target (default 0.4)

All env flags are read on every call (no module-level caching) so the
tests can flip them dynamically without import-side side effects.
"""

from __future__ import annotations

import math
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from delfin._bond_decollapse import _is_metal as _bd_is_metal
from delfin.fffree.grip_mogul_lookup import GripLibrary
from delfin.fffree.mogul_bounds import (
    build_bounds_matrix,
    detect_resonance_groups,
    hyb_str,
)
from delfin.fffree.mogul_severity import (
    DEFAULT_SIGMA_EQ,
    mahalanobis_severity,
)
from delfin.fffree.mogul_solver import (
    DEFAULT_MAX_ITER,
    DEFAULT_MD_DRIFT_TOL,
    DEFAULT_N_RESTARTS,
    DEFAULT_SEED,
    DEFAULT_TOL,
    solve_dg,
)


__all__ = [
    "mogul_embed",
    "mogul_embed_from_assembled",
    "is_active",
    "should_replace_etkdg",
    "should_replace_rigid",
    "verify_donor_shell_integrity",
    "DEFAULT_MD_DRIFT_TOL",
]


# ---------------------------------------------------------------------------
# Env flag accessors — fresh-read on every call.
# ---------------------------------------------------------------------------
def _env_truthy(name: str) -> bool:
    raw = os.environ.get(name, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        v = int(raw)
        if v > 0:
            return v
    except (TypeError, ValueError):
        pass
    return int(default)


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if math.isfinite(v) and v > 0.0:
            return v
    except (TypeError, ValueError):
        pass
    return float(default)


def is_active() -> bool:
    """Return True iff any Mogul-DG env flag is set.

    Default OFF byte-identical — when no env flag is set, this function
    returns False and the caller stays on the legacy path.
    """
    return (
        _env_truthy("DELFIN_FFFREE_MOGUL_DG")
        or _env_truthy("DELFIN_FFFREE_MOGUL_DG_REPLACE_ETKDG")
        or _env_truthy("DELFIN_FFFREE_MOGUL_DG_REPLACE_RIGID")
    )


def should_replace_etkdg() -> bool:
    """True iff ETKDG should be replaced in the F2 embed fallback."""
    return _env_truthy("DELFIN_FFFREE_MOGUL_DG_REPLACE_ETKDG")


def should_replace_rigid() -> bool:
    """True iff rigid-fit should be replaced in the native converter path."""
    return _env_truthy("DELFIN_FFFREE_MOGUL_DG_REPLACE_RIGID")


# ---------------------------------------------------------------------------
# Helper: build the bond / angle / md / resonance priors from the CCDC library
# ---------------------------------------------------------------------------
def _bonds_sorted(mol) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for b in mol.GetBonds():
        a = int(b.GetBeginAtomIdx())
        c = int(b.GetEndAtomIdx())
        if a == c:
            continue
        out.append((min(a, c), max(a, c)))
    out.sort()
    return out


def _angle_triples_sorted(mol) -> List[Tuple[int, int, int]]:
    """All 1,3 triples (i, j, k) with j at the vertex; deterministic order."""
    out: List[Tuple[int, int, int]] = []
    seen = set()
    n = mol.GetNumAtoms()
    for j in range(n):
        try:
            atom = mol.GetAtomWithIdx(int(j))
        except Exception:
            continue
        nbrs = sorted(int(nb.GetIdx()) for nb in atom.GetNeighbors())
        for a_pos in range(len(nbrs)):
            for b_pos in range(a_pos + 1, len(nbrs)):
                i, k = nbrs[a_pos], nbrs[b_pos]
                lo, hi = (i, k) if i < k else (k, i)
                key = (j, lo, hi)
                if key in seen:
                    continue
                seen.add(key)
                out.append((lo, j, hi))
    out.sort(key=lambda t: (t[1], t[0], t[2]))
    return out


def _collect_priors(
    syms: Sequence[str],
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
    grip_lib: Optional[GripLibrary],
    cod_lib: Optional[GripLibrary],
    lower: np.ndarray,
    upper: np.ndarray,
) -> Dict[str, Any]:
    """Extract bond / angle / M-D / resonance priors for the severity function.

    Each prior is a (μ, σ) pulled from the MIDPOINT and half-width of the
    bounds matrix.  This is correct by construction: Phase A already
    encoded the empirical μ/σ via μ ± k·σ → midpoint = μ, half-width =
    k·σ.  We invert that here to recover (μ, σ) consistent with the
    bounds.

    Returns a dict ready to spread into :func:`mahalanobis_severity`.
    """
    bonds = _bonds_sorted(mol)
    triples = _angle_triples_sorted(mol)
    donor_set = set(int(d) for d in donor_idxs if int(d) != int(metal_idx))

    # σ-multipliers (recovered from env, same defaults as mogul_bounds)
    bond_mult = _env_float("DELFIN_FFFREE_MOGUL_DG_BOUNDS_TOL", 3.0)
    md_mult = _env_float("DELFIN_FFFREE_MOGUL_DG_MD_BOUNDS_TOL", 2.0)

    bond_pairs: List[Tuple[int, int]] = []
    bond_priors: List[Tuple[float, float]] = []
    md_pairs: List[Tuple[int, int]] = []
    md_priors: List[Tuple[float, float]] = []
    angle_triples_out: List[Tuple[int, int, int]] = []
    angle_priors: List[Tuple[float, float]] = []

    # --- bond / md split ---
    for (i, j) in bonds:
        lo = float(lower[i, j])
        hi = float(upper[i, j])
        if hi >= 5e5:
            # No upper bound on this pair (unusual for a bond) — skip.
            continue
        is_md = (
            (i == int(metal_idx) and j in donor_set)
            or (j == int(metal_idx) and i in donor_set)
        )
        mu = 0.5 * (lo + hi)
        if is_md:
            sigma = max(0.05, (hi - lo) / (2.0 * md_mult))
            md_pairs.append((i, j))
            md_priors.append((mu, sigma))
        else:
            sigma = max(0.01, (hi - lo) / (2.0 * bond_mult))
            bond_pairs.append((i, j))
            bond_priors.append((mu, sigma))

    # --- angle priors: read from the 1,3 distance bounds ---
    # The angle bound was encoded as a 1,3 distance window via law-of-cosines.
    # For the severity we want θ in degrees.  We compute the *current* θ
    # implied by the midpoint distance using the bonded means d12, d23, then
    # treat the half-width as σ_θ via the same propagation.  This is an
    # approximation good to first order; for the severity it is more important
    # that θ_mean be RIGHT than that σ_θ be perfectly recovered.
    for (i, j, k) in triples:
        if int(j) == int(metal_idx):
            continue  # D-M-D handled by donor-donor bounds
        lo = float(lower[i, k])
        hi = float(upper[i, k])
        if hi >= 5e5:
            continue  # no 1,3 bound on this pair
        d_mid = 0.5 * (lo + hi)
        d_band = max(1e-3, 0.5 * (hi - lo))
        # Recover d12 and d23 from the bonded means (use bond bounds we just
        # collected, or fall back to (lo+hi)/2 of the bonded entries).
        d12_lo = float(lower[i, j])
        d12_hi = float(upper[i, j])
        d23_lo = float(lower[j, k])
        d23_hi = float(upper[j, k])
        if d12_hi >= 5e5 or d23_hi >= 5e5:
            continue
        d12 = 0.5 * (d12_lo + d12_hi)
        d23 = 0.5 * (d23_lo + d23_hi)
        # law of cosines: cos θ = (d12² + d23² - d13²) / (2·d12·d23)
        try:
            cos_theta = (d12 * d12 + d23 * d23 - d_mid * d_mid) / (2.0 * d12 * d23)
            cos_theta = max(-1.0 + 1e-6, min(1.0 - 1e-6, cos_theta))
            theta_deg = math.degrees(math.acos(cos_theta))
        except Exception:
            continue
        # σ_θ via finite-difference around d_mid → d_mid ± d_band
        try:
            cos_lo = (d12 * d12 + d23 * d23 - (d_mid + d_band) ** 2) / (
                2.0 * d12 * d23
            )
            cos_hi = (d12 * d12 + d23 * d23 - (d_mid - d_band) ** 2) / (
                2.0 * d12 * d23
            )
            cos_lo = max(-1.0 + 1e-6, min(1.0 - 1e-6, cos_lo))
            cos_hi = max(-1.0 + 1e-6, min(1.0 - 1e-6, cos_hi))
            theta_lo_deg = math.degrees(math.acos(cos_hi))  # cos decreasing
            theta_hi_deg = math.degrees(math.acos(cos_lo))
            sigma_theta_deg = max(1.0, 0.5 * (theta_hi_deg - theta_lo_deg))
        except Exception:
            sigma_theta_deg = 5.0
        angle_triples_out.append((i, j, k))
        angle_priors.append((theta_deg, sigma_theta_deg))

    # --- resonance groups → list of bond-pair lists ---
    try:
        res_groups_atoms = detect_resonance_groups(
            mol, metal_idx=int(metal_idx), use_automorphism=True,
        )
    except Exception:
        res_groups_atoms = []
    resonance_groups: List[List[Tuple[int, int]]] = []
    bond_set = set(bonds)
    if res_groups_atoms:
        # Build an adjacency for parent-radial bond resolution.
        n = int(len(syms))
        adj: Dict[int, List[int]] = {i: [] for i in range(n)}
        for (a, b) in bonds:
            adj[a].append(b)
            adj[b].append(a)
        for group in res_groups_atoms:
            atoms = sorted(int(a) for a in group)
            equi: List[Tuple[int, int]] = []
            group_set = set(atoms)
            for ap in range(len(atoms)):
                for bp in range(ap + 1, len(atoms)):
                    a, b = atoms[ap], atoms[bp]
                    if (a, b) in bond_set:
                        equi.append((a, b))
            # parent-radial pairs (e.g. carboxylate C-O × 2)
            parent_candidates: set = set()
            for x in atoms:
                for nb in adj.get(x, []):
                    if nb not in group_set:
                        parent_candidates.add(int(nb))
            for parent in sorted(parent_candidates):
                if _bd_is_metal(str(syms[parent])):
                    continue
                kids = sorted(int(k) for k in adj.get(parent, []) if k in group_set)
                if len(kids) >= 2:
                    for k in kids:
                        equi.append((min(parent, k), max(parent, k)))
            equi = sorted(set(equi))
            if len(equi) >= 2:
                resonance_groups.append(equi)

    return {
        "bond_pairs": bond_pairs,
        "bond_priors": bond_priors,
        "angle_triples": angle_triples_out,
        "angle_priors": angle_priors,
        "md_pairs": md_pairs,
        "md_priors": md_priors,
        "resonance_groups": resonance_groups,
    }


def _build_severity_callable(priors: Dict[str, Any]):
    """Wrap :func:`mahalanobis_severity` into a closure ``f(P) -> (loss, grad)``."""
    bond_pairs = priors["bond_pairs"]
    bond_priors = priors["bond_priors"]
    angle_triples = priors["angle_triples"]
    angle_priors = priors["angle_priors"]
    md_pairs = priors["md_pairs"]
    md_priors = priors["md_priors"]
    resonance_groups = priors["resonance_groups"]

    def _severity(P: np.ndarray) -> Tuple[float, np.ndarray]:
        return mahalanobis_severity(
            P,
            bond_pairs=bond_pairs,
            bond_priors=bond_priors,
            angle_triples=angle_triples,
            angle_priors=angle_priors,
            md_pairs=md_pairs,
            md_priors=md_priors,
            resonance_groups=resonance_groups,
            sigma_eq=DEFAULT_SIGMA_EQ,
        )

    return _severity


# ---------------------------------------------------------------------------
# Deterministic initial-cloud generator (used when no P_init is supplied)
# ---------------------------------------------------------------------------
def _random_initial_cloud(n_atoms: int, seed: int) -> np.ndarray:
    """Deterministic random scatter inside a sphere of radius ~ n_atoms^(1/3) Å.

    This is a *last-resort* P_init when no rigid-fit / ETKDG fallback is
    available.  The solver will pull the cloud onto the bounds-feasible
    set; the seed makes the result byte-deterministic.
    """
    rng = np.random.RandomState(seed=int(seed) & 0x7FFFFFFF)
    r = max(1.5, float(n_atoms) ** (1.0 / 3.0))
    P = rng.uniform(-r, r, size=(int(n_atoms), 3)).astype(np.float64)
    # Centre on origin so the metal lives near (0,0,0) — easier for the
    # solver to converge.
    P -= np.mean(P, axis=0, keepdims=True)
    return P


def _md_drift_acceptable(
    P_init: np.ndarray,
    P_final: np.ndarray,
    md_pairs: Sequence[Tuple[int, int]],
    tol: float,
) -> bool:
    """Check post-solve M-D drift against the SPEC §6 0.5 Å guard."""
    if not md_pairs or tol <= 0.0:
        return True
    for (i, j) in md_pairs:
        ii, jj = int(i), int(j)
        d0 = float(np.linalg.norm(P_init[ii] - P_init[jj]))
        d1 = float(np.linalg.norm(P_final[ii] - P_final[jj]))
        if abs(d1 - d0) > tol:
            return False
    return True


# ---------------------------------------------------------------------------
# Donor-shell integrity helpers (Task: TUMQAT Cl-drift / spurious-H first-shell)
# ---------------------------------------------------------------------------
# First-shell physical window (Å).  These are universal lower/upper guards
# only — the per-donor μ from the bounds matrix supplies the precise target.
_MD_SHELL_LOWER_A = 1.6
_MD_SHELL_UPPER_A = 2.6


def _md_target_from_bounds(
    lower: np.ndarray,
    upper: np.ndarray,
    metal_idx: int,
    donor_idx: int,
) -> Optional[float]:
    """Return the empirical M-D midpoint encoded by the bounds matrix.

    Returns ``None`` if the bounds entry is unbounded (sentinel ≥ 5e5).
    """
    lo = float(lower[int(metal_idx), int(donor_idx)])
    hi = float(upper[int(metal_idx), int(donor_idx)])
    if hi >= 5e5:
        return None
    return 0.5 * (lo + hi)


def verify_donor_shell_integrity(
    P: np.ndarray,
    syms: Sequence[str],
    metal_idx: int,
    declared_donors: Sequence[int],
    md_targets: Optional[Dict[int, float]] = None,
    md_tol: float = 0.4,
    shell_lower: float = _MD_SHELL_LOWER_A,
    shell_upper: float = _MD_SHELL_UPPER_A,
) -> Tuple[bool, List[str]]:
    """Verify post-embed first-shell integrity for declared donors.

    Universal (geometry + graph only — no SMILES patterns).  Checks:

    1.  Every declared donor sits in the physical M-D window
        ``[shell_lower, shell_upper]`` Å.
    2.  For each declared donor with an empirical μ in ``md_targets``, the
        actual M-D distance is within ``md_tol`` of μ.
    3.  No spurious H atom (not declared as a donor) is in the first shell
        — defined as closer to the metal than the worst declared heavy
        donor (or ``shell_upper`` Å when no declared donor exists).
    4.  Realised CN equals the declared donor count when only the
        physical shell window is considered.

    Returns ``(is_intact, violations)`` where ``violations`` is a list of
    human-readable strings (empty when intact).
    """
    violations: List[str] = []
    P_arr = np.asarray(P, dtype=np.float64)
    m = int(metal_idx)
    declared = sorted(set(int(d) for d in declared_donors if int(d) != m))
    if not declared:
        return True, []

    # 1 + 2: declared donor windows
    worst_md = float(shell_upper)
    for d in declared:
        dist = float(np.linalg.norm(P_arr[m] - P_arr[d]))
        if not (shell_lower <= dist <= shell_upper):
            violations.append(
                f"declared donor idx={d} sym={syms[d]} out of shell window: "
                f"d(M-D)={dist:.3f} Å not in [{shell_lower:.2f}, {shell_upper:.2f}]"
            )
        worst_md = max(worst_md, dist)
        if md_targets and d in md_targets and md_targets[d] is not None:
            mu = float(md_targets[d])
            if abs(dist - mu) > md_tol:
                violations.append(
                    f"declared donor idx={d} sym={syms[d]} far from target: "
                    f"d(M-D)={dist:.3f} Å, μ={mu:.3f} Å, |Δ|>{md_tol:.2f}"
                )

    # 3: spurious H atoms in first shell
    declared_set = set(declared)
    # Define "first shell" geometrically: closer than max(worst_md, shell_upper).
    shell_cut = max(worst_md, float(shell_upper))
    for i, s in enumerate(syms):
        if i == m or i in declared_set:
            continue
        if str(s).strip() != "H":
            continue
        dist = float(np.linalg.norm(P_arr[m] - P_arr[i]))
        if dist <= shell_cut:
            violations.append(
                f"spurious H idx={i} in first shell at d={dist:.3f} Å "
                f"(shell_cut={shell_cut:.3f} Å)"
            )

    # 4: realised CN matches declared count
    heavy_in_shell = sum(
        1
        for i, s in enumerate(syms)
        if i != m
        and str(s).strip() != "H"
        and float(np.linalg.norm(P_arr[m] - P_arr[i])) <= shell_upper
    )
    if heavy_in_shell != len(declared):
        violations.append(
            f"realised heavy CN={heavy_in_shell} != declared donor count="
            f"{len(declared)}"
        )

    return (len(violations) == 0), violations


def _bonded_heavy_parent(
    mol,
    h_idx: int,
    metal_idx: int,
) -> Optional[int]:
    """Return the bonded heavy parent atom index of an H, or None.

    Universal — uses the molecular graph only.  Excludes the metal even
    if the H is bonded to it (which would be unusual).
    """
    try:
        atom = mol.GetAtomWithIdx(int(h_idx))
    except Exception:
        return None
    for nb in atom.GetNeighbors():
        nbi = int(nb.GetIdx())
        if nbi == int(metal_idx):
            continue
        if str(nb.GetSymbol()).strip() == "H":
            continue
        return nbi
    return None


def _push_h_out_of_first_shell(
    P: np.ndarray,
    syms: Sequence[str],
    mol,
    metal_idx: int,
    declared_donors: Sequence[int],
    shell_upper: float = _MD_SHELL_UPPER_A,
    ch_length: float = 1.09,
) -> Tuple[np.ndarray, int]:
    """Repair spurious first-shell H atoms by moving them to their parent's C-H sphere.

    For every H atom that (a) is not a declared donor and (b) sits inside
    the first shell, relocate it to ``P[parent] + ch_length * dir(parent → H_original)``,
    a deterministic graph-only operation.  H atoms with no heavy parent
    are pushed radially outward to ``2 * shell_upper`` Å along the
    current ``metal → H`` direction.

    Returns ``(P_repaired, n_pushed)``.
    """
    P_new = np.asarray(P, dtype=np.float64).copy()
    m = int(metal_idx)
    declared_set = set(int(d) for d in declared_donors if int(d) != m)
    metal_pos = P_new[m]
    n_pushed = 0
    # Deterministic iteration order (sorted indices).
    for i in sorted(range(len(syms))):
        if i == m or i in declared_set:
            continue
        if str(syms[i]).strip() != "H":
            continue
        dist = float(np.linalg.norm(P_new[i] - metal_pos))
        if dist > float(shell_upper):
            continue
        parent = _bonded_heavy_parent(mol, i, m) if mol is not None else None
        if parent is not None:
            parent_pos = P_new[parent]
            # Direction priority:
            #   (1) current parent → H vector if it already points AWAY from
            #       the metal (cos(parent→H, parent→metal) < 0).
            #   (2) otherwise, the parent → away-from-metal direction
            #       (reflect H to the far side of the parent).
            #   (3) fallback unit vector when the parent sits on the metal.
            ph_vec = P_new[i] - parent_pos
            ph_norm = float(np.linalg.norm(ph_vec))
            pm_vec = metal_pos - parent_pos
            pm_norm = float(np.linalg.norm(pm_vec))
            direction: Optional[np.ndarray] = None
            if ph_norm > 1e-6 and pm_norm > 1e-6:
                cos_pm = float(np.dot(ph_vec, pm_vec) / (ph_norm * pm_norm))
                if cos_pm < 0.0:
                    direction = ph_vec / ph_norm
            if direction is None:
                if pm_norm > 1e-6:
                    # away-from-metal direction
                    direction = -pm_vec / pm_norm
                elif ph_norm > 1e-6:
                    direction = ph_vec / ph_norm
                else:
                    direction = np.array([1.0, 0.0, 0.0], dtype=np.float64)
            P_new[i] = parent_pos + float(ch_length) * direction
        else:
            # No parent — radial push outward.
            d_vec = P_new[i] - metal_pos
            d_norm = float(np.linalg.norm(d_vec))
            if d_norm > 1e-6:
                direction = d_vec / d_norm
            else:
                direction = np.array([1.0, 0.0, 0.0], dtype=np.float64)
            P_new[i] = metal_pos + 2.0 * float(shell_upper) * direction
        n_pushed += 1
    return P_new, n_pushed


def _tighten_bounds_for_missing_donor(
    lower: np.ndarray,
    upper: np.ndarray,
    metal_idx: int,
    donor_idx: int,
    half_width_factor: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Halve the M-D bounds half-width around its current midpoint.

    Tightens the M-D window for one donor that drifted out of the first
    shell.  Returns NEW arrays (does not mutate inputs).  When the
    current bound is unbounded (sentinel ≥ 5e5) the bounds are returned
    unchanged.
    """
    lo = float(lower[int(metal_idx), int(donor_idx)])
    hi = float(upper[int(metal_idx), int(donor_idx)])
    if hi >= 5e5 or hi <= lo:
        return lower, upper
    mu = 0.5 * (lo + hi)
    half = max(0.01, 0.5 * (hi - lo) * float(half_width_factor))
    new_lo = mu - half
    new_hi = mu + half
    L = lower.copy()
    U = upper.copy()
    L[int(metal_idx), int(donor_idx)] = new_lo
    L[int(donor_idx), int(metal_idx)] = new_lo
    U[int(metal_idx), int(donor_idx)] = new_hi
    U[int(donor_idx), int(metal_idx)] = new_hi
    return L, U


# ---------------------------------------------------------------------------
# Public driver — from already-assembled topology
# ---------------------------------------------------------------------------
def mogul_embed_from_assembled(
    syms: Sequence[str],
    P_init: np.ndarray,
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
    *,
    grip_lib: Optional[GripLibrary] = None,
    cod_lib: Optional[GripLibrary] = None,
    geometry: Optional[str] = None,
    max_iter: Optional[int] = None,
    seed: Optional[int] = None,
    tol: float = DEFAULT_TOL,
    n_restarts: Optional[int] = None,
    md_drift_tol: Optional[float] = None,
) -> Optional[Tuple[List[str], np.ndarray, Dict[str, Any]]]:
    """Refine an already-assembled (syms, P_init, mol, metal, donors) tuple.

    This is the integration entry point for the native converter path:
    when ``DELFIN_FFFREE_MOGUL_DG_REPLACE_RIGID=1`` is set, callers pass
    the result of ``assemble_from_config`` (rigid-fit) as ``P_init`` and
    receive a Mogul-refined geometry back.  If anything fails — bounds
    construction, solver, M-D drift guard — we return ``None`` so the
    caller silently keeps the rigid-fit result (fail-open contract).

    Parameters
    ----------
    syms : sequence of str
        Element symbols.
    P_init : (n_atoms, 3) ndarray
        Initial coordinates (typically from rigid-fit).
    mol : RDKit Mol
        Molecule with bonds + hybridisation populated.
    metal_idx : int
        Metal atom index in ``syms`` / ``mol``.
    donor_idxs : sequence of int
        Donor atom indices.
    grip_lib, cod_lib : GripLibrary, optional
        Pre-loaded libraries.  Default: load the pinned defaults.
    geometry : str, optional
        Polyhedron name (e.g. ``"OC-6 octahedron"``).
    max_iter, seed, n_restarts, md_drift_tol : optional
        Solver knobs.  Defaults come from env or SPEC defaults.

    Returns
    -------
    (syms, P, info) tuple, or ``None`` on failure.
    """
    try:
        P_init_arr = np.asarray(P_init, dtype=np.float64)
        if P_init_arr.ndim != 2 or P_init_arr.shape[1] != 3:
            return None
        n_atoms = int(P_init_arr.shape[0])
        if n_atoms != len(syms):
            return None
        if not np.all(np.isfinite(P_init_arr)):
            return None

        # Build bounds matrix from the graph + library
        lower, upper, bounds_info = build_bounds_matrix(
            list(syms),
            mol,
            int(metal_idx),
            list(donor_idxs),
            grip_lib=grip_lib,
            cod_lib=cod_lib,
            geometry=geometry,
        )

        # Extract priors from bounds — single source of truth for (μ, σ)
        priors = _collect_priors(
            list(syms), mol, int(metal_idx), list(donor_idxs),
            grip_lib, cod_lib, lower, upper,
        )
        severity_fn = _build_severity_callable(priors)

        # M-D pairs for the drift guard
        md_pairs_list = [
            (min(int(metal_idx), int(d)), max(int(metal_idx), int(d)))
            for d in donor_idxs if int(d) != int(metal_idx)
        ]

        # Frozen indices: metal stays put across restarts so the
        # coordination-sphere anchor doesn't drift between restarts.
        frozen_indices = [int(metal_idx)]

        # Solver knobs (env > default)
        _max_iter = max_iter if max_iter is not None else _env_int(
            "DELFIN_FFFREE_MOGUL_DG_MAX_ITER", DEFAULT_MAX_ITER,
        )
        _seed = seed if seed is not None else _env_int(
            "DELFIN_FFFREE_MOGUL_DG_SEED", DEFAULT_SEED,
        )
        _n_restarts = n_restarts if n_restarts is not None else _env_int(
            "DELFIN_FFFREE_MOGUL_DG_N_RESTARTS", DEFAULT_N_RESTARTS,
        )
        _md_drift_tol = md_drift_tol if md_drift_tol is not None else _env_float(
            "DELFIN_FFFREE_MOGUL_DG_MD_DRIFT_TOL", DEFAULT_MD_DRIFT_TOL,
        )

        # Donor-shell integrity / repair env flags (default-OFF byte-identical).
        _verify_on = _env_truthy("DELFIN_FFFREE_MOGUL_DG_DONOR_SHELL_VERIFY")
        _reject_h_on = _env_truthy("DELFIN_FFFREE_MOGUL_DG_REJECT_H_FIRST_SHELL")
        _repair_retries = _env_int(
            "DELFIN_FFFREE_MOGUL_DG_DONOR_REPAIR_RETRIES", 3,
        )
        _md_tol_a = _env_float("DELFIN_FFFREE_MOGUL_DG_MD_TOL", 0.4)
        # Pre-compute per-donor empirical μ targets from the bounds matrix
        # (single source of truth — Phase A already encoded the (μ, σ) here).
        declared_donor_list = [
            int(d) for d in donor_idxs if int(d) != int(metal_idx)
        ]
        md_target_map: Dict[int, float] = {}
        for d in declared_donor_list:
            mu = _md_target_from_bounds(lower, upper, int(metal_idx), d)
            if mu is not None:
                md_target_map[d] = float(mu)

        # Repair loop — solver + verify + tighten bounds.  When verification
        # is disabled the loop runs exactly once (preserving byte-identical
        # default-OFF behaviour).
        max_attempts = max(1, int(_repair_retries) + 1) if _verify_on else 1
        L_cur = lower
        U_cur = upper
        P_solved: Optional[np.ndarray] = None
        solver_info: Dict[str, Any] = {}
        repair_history: List[Dict[str, Any]] = []
        for attempt in range(max_attempts):
            # Re-derive the severity callable when bounds were tightened so
            # the priors stay in sync.  First attempt uses the originals.
            if attempt == 0:
                _sev_fn = severity_fn
            else:
                _priors_iter = _collect_priors(
                    list(syms), mol, int(metal_idx), list(donor_idxs),
                    grip_lib, cod_lib, L_cur, U_cur,
                )
                _sev_fn = _build_severity_callable(_priors_iter)

            P_solved, solver_info = solve_dg(
                P_init_arr,
                L_cur,
                U_cur,
                _sev_fn,
                max_iter=int(_max_iter),
                tol=float(tol),
                n_restarts=int(_n_restarts),
                seed=int(_seed),
                md_pairs=md_pairs_list,
                frozen_indices=frozen_indices,
                md_drift_tol=float(_md_drift_tol),
            )

            # Validate basic solver outcome
            if solver_info.get("failed", False):
                return None
            if not np.all(np.isfinite(P_solved)):
                return None
            if P_solved.shape != P_init_arr.shape:
                return None

            if not _verify_on:
                break  # legacy single-pass behaviour

            ok, violations = verify_donor_shell_integrity(
                P_solved, list(syms), int(metal_idx),
                declared_donor_list,
                md_targets=md_target_map,
                md_tol=float(_md_tol_a),
            )
            repair_history.append(
                {"attempt": attempt, "ok": ok, "violations": list(violations)}
            )
            if ok:
                break

            # Identify donors that drifted out of the physical window (or far
            # from μ) and tighten their bounds for the next attempt.
            tightened = False
            for d in declared_donor_list:
                dist = float(
                    np.linalg.norm(P_solved[int(metal_idx)] - P_solved[d])
                )
                out_of_shell = not (_MD_SHELL_LOWER_A <= dist <= _MD_SHELL_UPPER_A)
                mu = md_target_map.get(d)
                far_from_mu = mu is not None and abs(dist - mu) > float(_md_tol_a)
                if out_of_shell or far_from_mu:
                    L_cur, U_cur = _tighten_bounds_for_missing_donor(
                        L_cur, U_cur, int(metal_idx), d,
                        half_width_factor=0.5,
                    )
                    tightened = True
            if not tightened:
                # Bounds already exhausted — fall through.
                break

        # Final integrity verification + spurious-H rejection (optional).
        n_h_pushed = 0
        final_violations: List[str] = []
        final_ok = True
        if _verify_on:
            final_ok, final_violations = verify_donor_shell_integrity(
                P_solved, list(syms), int(metal_idx),
                declared_donor_list,
                md_targets=md_target_map,
                md_tol=float(_md_tol_a),
            )
        if _reject_h_on:
            P_solved, n_h_pushed = _push_h_out_of_first_shell(
                P_solved, list(syms), mol, int(metal_idx),
                declared_donor_list,
                shell_upper=_MD_SHELL_UPPER_A,
            )
            # Re-verify after H rejection (informational only).
            if _verify_on:
                final_ok, final_violations = verify_donor_shell_integrity(
                    P_solved, list(syms), int(metal_idx),
                    declared_donor_list,
                    md_targets=md_target_map,
                    md_tol=float(_md_tol_a),
                )

        # When verification is ON and STILL fails after all retries, treat
        # this as infeasibility and fall back to the caller's legacy path.
        # When OFF, downstream behaviour is unchanged.
        if _verify_on and not final_ok:
            # Repair retries exhausted — signal infeasibility.
            return None

        # Final M-D drift guard against the ORIGINAL P_init (the solver
        # may already enforce this; we double-check at the API boundary).
        if md_pairs_list and not _md_drift_acceptable(
            P_init_arr, P_solved, md_pairs_list, float(_md_drift_tol),
        ):
            return None

        info: Dict[str, Any] = {
            "n_atoms": n_atoms,
            "bounds": bounds_info,
            "solver": solver_info,
            "n_md_pairs": len(md_pairs_list),
            "n_bond_pairs": len(priors["bond_pairs"]),
            "n_angle_triples": len(priors["angle_triples"]),
            "n_resonance_groups": len(priors["resonance_groups"]),
            "max_md_drift": solver_info.get("max_md_drift", 0.0),
            "final_loss": solver_info.get("final_loss", float("inf")),
            "final_emp_loss": solver_info.get("final_emp_loss", float("inf")),
            "converged": solver_info.get("converged", False),
            "source": "mogul_dg",
        }
        if _verify_on or _reject_h_on:
            info["donor_shell_verify"] = bool(_verify_on)
            info["donor_shell_reject_h"] = bool(_reject_h_on)
            info["donor_shell_ok"] = bool(final_ok)
            info["donor_shell_violations"] = list(final_violations)
            info["donor_shell_repair_attempts"] = len(repair_history)
            info["donor_shell_repair_history"] = repair_history
            info["donor_shell_h_pushed"] = int(n_h_pushed)
        return list(syms), np.asarray(P_solved, dtype=np.float64), info
    except Exception:
        # Fail-open
        return None


# ---------------------------------------------------------------------------
# Public driver — from SMILES
# ---------------------------------------------------------------------------
def _decompose_smiles_to_topology(smiles: str):
    """Decompose a SMILES into (syms, mol, metal_idx, donor_idxs, geometry).

    Uses the project's ``_prepare_mol_for_embedding`` and the standard
    metal-detection routine.  Returns None on any parse failure.
    """
    try:
        from delfin.smiles_converter import _prepare_mol_for_embedding
    except Exception:
        return None
    try:
        mol = _prepare_mol_for_embedding(smiles, hapto_approx=False)
    except Exception:
        mol = None
    if mol is None:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol = Chem.AddHs(mol)
        except Exception:
            return None
    if mol is None:
        return None
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    metal_indices = [i for i, s in enumerate(syms) if _bd_is_metal(s)]
    if not metal_indices:
        # No metal — mogul_dg is for TMCs; signal back so caller falls back.
        return None
    metal_idx = int(metal_indices[0])
    donors = sorted(
        int(nb.GetIdx())
        for nb in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
        if not _bd_is_metal(nb.GetSymbol())
    )
    if not donors:
        return None
    # Geometry hint via decompose (best-effort; None is acceptable)
    geometry: Optional[str] = None
    try:
        from delfin.fffree import decompose as DEC
        d = DEC.decompose(smiles)
        if d is not None:
            geometry = d.get("geometry")
    except Exception:
        pass
    return syms, mol, metal_idx, donors, geometry


def mogul_embed(
    smiles: str,
    *,
    grip_lib: Optional[GripLibrary] = None,
    cod_lib: Optional[GripLibrary] = None,
    max_iter: Optional[int] = None,
    seed: Optional[int] = None,
    tol: float = DEFAULT_TOL,
    n_restarts: Optional[int] = None,
    md_drift_tol: Optional[float] = None,
    P_init: Optional[np.ndarray] = None,
) -> Optional[Tuple[List[str], np.ndarray, Dict[str, Any]]]:
    """Mogul-informed whole-complex distance-geometry embed (public API).

    Parameters
    ----------
    smiles : str
        Input SMILES.
    grip_lib, cod_lib : GripLibrary, optional
        Pre-loaded empirical libraries.  Default loads the pinned files.
    max_iter, seed, n_restarts, md_drift_tol : optional
        Solver knobs (env > default).
    P_init : (n_atoms, 3) ndarray, optional
        Pre-built initial cloud (e.g. from a rigid-fit).  When omitted, a
        deterministic random cloud is used as the starting point.

    Returns
    -------
    ``(syms, P, info)`` on success, ``None`` on infeasibility.  Never
    raises — every internal exception is caught and converted to None.
    """
    try:
        topo = _decompose_smiles_to_topology(smiles)
        if topo is None:
            return None
        syms, mol, metal_idx, donors, geometry = topo

        # If no P_init supplied, use a deterministic random cloud and disable
        # the M-D drift guard (there's no "true" M-D distance to drift away
        # from; the solver is meant to PULL the random cloud onto the
        # bounds-feasible set, which by definition moves M-D distances).
        _drift_for_random = 0.0  # 0 = guard disabled
        if P_init is None:
            _seed = seed if seed is not None else _env_int(
                "DELFIN_FFFREE_MOGUL_DG_SEED", DEFAULT_SEED,
            )
            P_init_arr = _random_initial_cloud(len(syms), int(_seed))
            effective_drift_tol = (
                md_drift_tol if md_drift_tol is not None else _drift_for_random
            )
        else:
            P_init_arr = np.asarray(P_init, dtype=np.float64)
            effective_drift_tol = md_drift_tol

        return mogul_embed_from_assembled(
            syms,
            P_init_arr,
            mol,
            metal_idx,
            donors,
            grip_lib=grip_lib,
            cod_lib=cod_lib,
            geometry=geometry,
            max_iter=max_iter,
            seed=seed,
            tol=tol,
            n_restarts=n_restarts,
            md_drift_tol=effective_drift_tol,
        )
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Convenience: XYZ formatter (matches converter_backend._xyz exactly)
# ---------------------------------------------------------------------------
def xyz_block(syms: Sequence[str], P: np.ndarray) -> str:
    """Canonical headerless XYZ atom block — matches converter_backend._xyz()."""
    return "\n".join(
        f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}"
        for s, (x, y, z) in zip(syms, P)
    )


# ---------------------------------------------------------------------------
# Module-level smoke (run only when executed directly)
# ---------------------------------------------------------------------------
def _self_check() -> None:  # pragma: no cover
    """Tiny end-to-end smoke for local debugging."""
    res = mogul_embed("N[Pt](N)(Cl)Cl")
    if res is None:
        print("mogul_dg self-check: None (library may be missing)")
        return
    syms, P, info = res
    print(f"mogul_dg self-check: {len(syms)} atoms, "
          f"final_loss={info['final_loss']:.3f}, "
          f"converged={info['converged']}")


if __name__ == "__main__":  # pragma: no cover
    _self_check()
