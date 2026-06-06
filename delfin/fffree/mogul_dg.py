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

        # Run the projected L-BFGS solver
        P_solved, solver_info = solve_dg(
            P_init_arr,
            lower,
            upper,
            severity_fn,
            max_iter=int(_max_iter),
            tol=float(tol),
            n_restarts=int(_n_restarts),
            seed=int(_seed),
            md_pairs=md_pairs_list,
            frozen_indices=frozen_indices,
            md_drift_tol=float(_md_drift_tol),
        )

        # Validate
        if solver_info.get("failed", False):
            return None
        if not np.all(np.isfinite(P_solved)):
            return None
        if P_solved.shape != P_init_arr.shape:
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
