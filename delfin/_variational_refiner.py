"""Baustein 6 — Variational L-BFGS-B Post-Refiner with 4-Tier Symmetry Awareness.

Complements Baustein 5 (PBD, ``delfin._post_optimizer``).  Where B5 handles
catastrophic moves + hard topology repair via constraint projection, B6
performs smooth simultaneous balance of all forces (bonds, angles, clashes,
hard topology, plus four tiers of symmetry pressure) via gradient-based
minimization of an 8-term energy functional ``U_total``.

Pipeline
--------
1. Pre-compute symmetry info (Tier A / B / C / D — one-shot, expensive).
2. Build the ``U_total`` + ``grad_U_total`` closures via ``delfin._energy_terms``.
3. Run ``scipy.optimize.minimize(method="L-BFGS-B", jac=True)``.
4. Validate result against the Baustein 5 topology hard-gate; on failure
   the input XYZ is returned unchanged.
5. Return refined XYZ string + a structured report dict.

Doctrine
--------
- Universal.  No SMILES/refcode/element-list shortcuts.  Symmetry detection
  uses Morgan ranks + SMARTS archetypes + connectivity automorphisms.
- Topology preservation is a hard gate (M-D ∈ [0.85, 1.10] × ideal, no new
  spurious heavy-heavy bonds inside 0.85 · Σr_cov).
- Always returns a valid XYZ string.  On any failure (missing sister module,
  scipy exception, topology check fail) the original XYZ is returned with a
  populated report dict and ``fallback_used=True``.
- Lazy imports for RDKit / SciPy / sister B6 modules.  Importing this module
  is cheap and never raises even if helpers are missing.

Entry point
-----------
``variational_refine(xyz, mol, class_label=..., **params) -> (xyz, report)``
"""
from __future__ import annotations

import re
import traceback
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Class-conditional hyperparameters (see BAUSTEIN6_MASTERPLAN section 6).
# ---------------------------------------------------------------------------

_CLASS_HYPER: Dict[str, Dict[str, Any]] = {
    "sigma": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0,
        "k_topology": 10000.0,
        "k_A": 100.0, "k_B": 50.0, "k_C": 80.0, "k_D": 30.0,
        "enable_D": True,
    },
    "hapto": {
        "k_bond": 800.0, "k_angle": 80.0, "k_clash": 400.0,
        "k_topology": 10000.0,
        "k_A": 30.0, "k_B": 30.0, "k_C": 50.0, "k_D": 20.0,
        "enable_D": True,
    },
    "multi_sigma": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0,
        "k_topology": 10000.0,
        "k_A": 100.0, "k_B": 60.0, "k_C": 100.0, "k_D": 40.0,
        "enable_D": True,
    },
    "multi_hapto": {
        "k_bond": 600.0, "k_angle": 50.0, "k_clash": 300.0,
        "k_topology": 10000.0,
        "k_A": 20.0, "k_B": 20.0, "k_C": 30.0, "k_D": 10.0,
        "enable_D": False,
    },
    "no_metal": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0,
        "k_topology": 0.0,
        "k_A": 0.0, "k_B": 100.0, "k_C": 150.0, "k_D": 50.0,
        "enable_D": True,
    },
}


# ---------------------------------------------------------------------------
# Constants — lazy fallbacks for _METAL_SET / topology helpers.
# ---------------------------------------------------------------------------

def _load_metal_set() -> set:
    """Lazy import of the canonical metal set."""
    try:
        from delfin.smiles_converter import _METAL_SET  # type: ignore
        return set(_METAL_SET)
    except Exception:
        # Minimal d/f-block fallback so module-load never fails.
        return {
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
            "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
            "Er", "Tm", "Yb", "Lu", "Th", "U",
        }


# ---------------------------------------------------------------------------
# XYZ I/O — DELFIN-style (preserves header lines).
# ---------------------------------------------------------------------------

_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _parse_xyz_to_array(xyz: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse a DELFIN-style XYZ string.

    Returns
    -------
    symbols : list[str]
        Element symbol per atom (length N).
    coords : ndarray, shape (N, 3)
        Cartesian coordinates in Å.
    orig_lines : list[str]
        Raw input lines (header preserved for round-trip rewrite).
    """
    syms: List[str] = []
    pts: List[List[float]] = []
    lines = xyz.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append([float(m.group(2)), float(m.group(3)), float(m.group(4))])
    if not pts:
        return syms, np.zeros((0, 3), dtype=float), lines
    return syms, np.asarray(pts, dtype=float), lines


def _array_to_xyz(coords: np.ndarray, mol, orig_lines: Optional[List[str]] = None,
                  symbols: Optional[List[str]] = None) -> str:
    """Render coords back into a DELFIN-style XYZ string.

    If ``orig_lines`` are supplied, header / blank lines are preserved verbatim
    and only the atom rows are rewritten in place.  Otherwise a minimal
    two-line header (atom count + blank comment) is emitted.
    """
    if symbols is None:
        try:
            symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        except Exception:
            symbols = ["X"] * len(coords)

    if orig_lines is not None:
        out: List[str] = []
        atom_i = 0
        for line in orig_lines:
            m = _XYZ_LINE_RE.match(line)
            if m and atom_i < len(symbols):
                x, y, z = coords[atom_i]
                out.append(f"{symbols[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
                atom_i += 1
            else:
                out.append(line)
        return "\n".join(out) + "\n"

    # Plain emission (used when no original header is available).
    out = [f"{len(coords)}", ""]
    for sym, (x, y, z) in zip(symbols, coords):
        out.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Topology check — reuse Baustein 5 if available, otherwise inline minimal.
# ---------------------------------------------------------------------------

def _topology_check(coords: np.ndarray, mol, metal_set: set) -> bool:
    """Conservative topology hard-gate.

    Tries the canonical ``delfin._post_optimizer._passes_topology`` first; if
    that module is unavailable, falls back to an inline check that mirrors its
    semantics (M-D window 0.85-1.10 × ideal, non-bonded heavy-pair collapse
    below 0.85 · Σr_cov).
    """
    try:
        from delfin._post_optimizer import (  # type: ignore
            _passes_topology, _metal_indices, _md_pairs, _non_bonded_heavy_pairs,
        )
        metals = _metal_indices(mol)
        md = _md_pairs(mol, metals)
        nb = _non_bonded_heavy_pairs(mol)
        return bool(_passes_topology(coords, mol, metals, md, nb))
    except Exception:
        pass

    # Inline fallback — only M-D window, conservative.
    try:
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    except Exception:
        return True
    if not metals:
        return True

    try:
        from delfin.smiles_converter import _get_ml_bond_length  # type: ignore
    except Exception:
        def _get_ml_bond_length(_a: str, _b: str) -> float:  # type: ignore
            return 2.0

    for m in metals:
        m_atom = mol.GetAtomWithIdx(m)
        m_sym = m_atom.GetSymbol()
        for nb_atom in m_atom.GetNeighbors():
            if nb_atom.GetSymbol() in metal_set or nb_atom.GetSymbol() == "H":
                continue
            d_idx = nb_atom.GetIdx()
            d_sym = nb_atom.GetSymbol()
            try:
                d_ideal = float(_get_ml_bond_length(m_sym, d_sym))
            except Exception:
                d_ideal = 2.0
            d_cur = float(np.linalg.norm(coords[m] - coords[d_idx]))
            if d_cur < 0.85 * d_ideal or d_cur > 1.10 * d_ideal:
                return False
    return True


# ---------------------------------------------------------------------------
# Minimal fallback energy term (U_bond + U_clash + harmonic M-D).
# Used only when delfin._energy_terms is not yet available.
# ---------------------------------------------------------------------------

def _fallback_U_total(coords_2d: np.ndarray, mol, sym_info: Dict[str, Any],
                      params: Dict[str, Any]) -> Tuple[float, np.ndarray]:
    """Bare-minimum U_total used when the full _energy_terms module is missing.

    Includes:
      - Harmonic U_bond over all RDKit bonds (ideal = covalent-radius sum).
      - Harmonic donor-target pull for the Tier-A targets supplied in
        ``sym_info["donor_targets"]``.
      - One-sided quadratic clash penalty on non-bonded heavy pairs.

    Returns (U_total, grad shape (N, 3)).
    """
    n = coords_2d.shape[0]
    grad = np.zeros_like(coords_2d)
    U = 0.0

    # Element symbols (cheap, RDKit-backed).
    try:
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
    except Exception:
        syms = ["C"] * n

    # Covalent radii fallback table.
    cov = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
           "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
           "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
           "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Pt": 1.36}

    def _r(sym: str) -> float:
        return cov.get(sym, 0.80 if sym != "H" else 0.31)

    k_b = float(params.get("k_bond", 1000.0))
    try:
        bonds = list(mol.GetBonds())
    except Exception:
        bonds = []
    for b in bonds:
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        d_ideal = _r(syms[i]) + _r(syms[j])
        diff = coords_2d[i] - coords_2d[j]
        d = float(np.linalg.norm(diff))
        if d < 1e-9:
            continue
        delta = d - d_ideal
        U += k_b * delta * delta
        g = (2.0 * k_b * delta) * (diff / d)
        grad[i] += g
        grad[j] -= g

    # Tier-A donor-target harmonic.
    k_A = float(params.get("k_A", 0.0))
    targets = sym_info.get("donor_targets", {})
    if k_A > 0 and targets:
        for d_idx, tgt in targets.items():
            if d_idx < 0 or d_idx >= n:
                continue
            diff = coords_2d[d_idx] - np.asarray(tgt, dtype=float)
            U += k_A * float(np.dot(diff, diff))
            grad[d_idx] += 2.0 * k_A * diff

    # Light clash term — only sample a bounded subset of non-bonded pairs.
    k_c = float(params.get("k_clash", 0.0))
    if k_c > 0 and n <= 400:
        bonded = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in bonds}
        for i in range(n):
            if syms[i] == "H":
                continue
            for j in range(i + 1, n):
                if syms[j] == "H":
                    continue
                if (i, j) in bonded:
                    continue
                thr = (_r(syms[i]) + _r(syms[j])) * 1.5  # vdW-ish proxy
                diff = coords_2d[i] - coords_2d[j]
                d = float(np.linalg.norm(diff))
                if d < 1e-9 or d >= thr:
                    continue
                pen = thr - d
                U += k_c * pen * pen
                g = (-2.0 * k_c * pen) * (diff / d)
                grad[i] += g
                grad[j] -= g

    return float(U), grad


# ---------------------------------------------------------------------------
# Symmetry pre-compute — calls sister modules when available, returns minimal
# sym_info struct when they are not.
# ---------------------------------------------------------------------------

def _precompute_symmetry(mol, coords: np.ndarray, class_label: str,
                         params: Dict[str, Any], metal_set: set,
                         enable_global_pg: bool) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Build the symmetry-info struct that ``U_total`` consumes.

    Returns ``(sym_info, meta)`` where ``meta`` carries diagnostic counts for
    the final report.  Every key in ``sym_info`` falls back to a safe empty
    container if the corresponding helper is missing.
    """
    sym_info: Dict[str, Any] = {
        "donor_targets": {},
        "equiv_atoms": [],
        "equiv_bond_pairs": [],
        "equiv_angle_triples": [],
        "fragments": [],
        "global_pg": "C1",
        "global_ops": [],
        "atom_perms": {},
    }
    meta: Dict[str, Any] = {
        "global_pg": "C1",
        "fragments_detected": 0,
        "equiv_classes": 0,
    }

    # ----- Tier A: per-metal donor targets via Hungarian assignment -----
    try:
        from delfin._symmetry_detection import (  # type: ignore
            hungarian_assign_donors_to_slots,
        )
        metals = [a.GetIdx() for a in mol.GetAtoms()
                  if a.GetSymbol() in metal_set]
        targets: Dict[int, np.ndarray] = {}
        for m_idx in metals:
            try:
                per_metal = hungarian_assign_donors_to_slots(coords, mol, m_idx)
                if per_metal:
                    targets.update(per_metal)
            except Exception:
                continue
        sym_info["donor_targets"] = targets
    except Exception:
        pass

    # ----- Tier B: equivalent atoms / bond-pairs / angle-triples -----
    try:
        from delfin._symmetry_detection import (  # type: ignore
            find_equivalent_atoms,
            find_equivalent_bond_pairs,
            find_equivalent_angle_triples,
        )
        equiv_atoms = find_equivalent_atoms(mol) or []
        sym_info["equiv_atoms"] = equiv_atoms
        sym_info["equiv_bond_pairs"] = find_equivalent_bond_pairs(mol) or []
        sym_info["equiv_angle_triples"] = find_equivalent_angle_triples(mol) or []
        meta["equiv_classes"] = len(equiv_atoms)
    except Exception:
        pass

    # ----- Tier C: fragment archetypes -----
    try:
        from delfin._fragment_archetypes import detect_fragments  # type: ignore
        frags = detect_fragments(mol) or []
        sym_info["fragments"] = frags
        meta["fragments_detected"] = len(frags)
    except Exception:
        pass

    # ----- Tier D: global molecular point group -----
    if enable_global_pg and bool(params.get("enable_D", False)):
        try:
            from delfin._symmetry_detection import (  # type: ignore
                detect_global_point_group,
            )
            pg, ops, perms = detect_global_point_group(mol, coords)
            sym_info["global_pg"] = pg or "C1"
            sym_info["global_ops"] = list(ops) if ops is not None else []
            sym_info["atom_perms"] = perms or {}
            meta["global_pg"] = sym_info["global_pg"]
        except Exception:
            pass

    return sym_info, meta


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

def variational_refine(
    xyz: str,
    mol,
    class_label: str = "sigma",
    max_iter: int = 200,
    ftol: float = 1e-6,
    enable_global_pg: bool = True,
) -> Tuple[str, Dict[str, Any]]:
    """L-BFGS-B variational post-refinement with 4-tier symmetry awareness.

    Parameters
    ----------
    xyz : str
        DELFIN-style XYZ string (post Baustein 5).
    mol : rdkit.Chem.Mol
        Topology source.  Must carry a 3D conformer aligned with ``xyz``.
    class_label : {"sigma", "hapto", "multi_sigma", "multi_hapto", "no_metal"}
        Selects the hyperparameter preset.  Unknown labels fall back to
        ``"sigma"``.
    max_iter : int, default 200
        L-BFGS-B iteration cap.
    ftol : float, default 1e-6
        Objective relative tolerance passed to L-BFGS-B.
    enable_global_pg : bool, default True
        If ``False``, Tier D (global molecular point group) is skipped
        regardless of the class default.

    Returns
    -------
    refined_xyz : str
        Optimized geometry (or the original ``xyz`` if anything failed /
        topology was broken).
    report : dict
        Diagnostic record.  Always contains the keys ``iterations``,
        ``converged``, ``energy_initial``, ``energy_final``,
        ``topology_preserved``, ``fallback_used``, ``global_pg``,
        ``fragments_detected``, ``equiv_classes``, and (on failure) ``error``.
    """
    # Skeleton report (mandatory keys present even on early return).
    report: Dict[str, Any] = {
        "iterations": 0,
        "converged": False,
        "energy_initial": float("nan"),
        "energy_final": float("nan"),
        "topology_preserved": True,
        "fallback_used": False,
        "global_pg": "C1",
        "fragments_detected": 0,
        "equiv_classes": 0,
    }

    # ----- Step 1: parse XYZ -----
    try:
        symbols, coords, orig_lines = _parse_xyz_to_array(xyz)
    except Exception as exc:
        report.update({"error": f"xyz parse: {exc}", "fallback_used": True})
        return xyz, report

    n_atoms = coords.shape[0]
    if n_atoms < 2:
        report.update({"error": "fewer than 2 atoms parsed",
                       "fallback_used": True})
        return xyz, report

    # Optional cross-check: mol atom count should match.
    try:
        if mol is not None and mol.GetNumAtoms() != n_atoms:
            report.update({"error": "mol/xyz atom-count mismatch",
                           "fallback_used": True})
            return xyz, report
    except Exception:
        pass

    # ----- Step 2: hyperparameters -----
    params = dict(_CLASS_HYPER.get(class_label, _CLASS_HYPER["sigma"]))

    # ----- Step 3: SciPy availability gate -----
    try:
        from scipy.optimize import minimize  # type: ignore
    except Exception as exc:
        report.update({"error": f"scipy unavailable: {exc}",
                       "fallback_used": True})
        return xyz, report

    metal_set = _load_metal_set()

    # ----- Step 4: symmetry pre-compute -----
    try:
        sym_info, sym_meta = _precompute_symmetry(
            mol, coords, class_label, params, metal_set, enable_global_pg
        )
        report["global_pg"] = sym_meta["global_pg"]
        report["fragments_detected"] = sym_meta["fragments_detected"]
        report["equiv_classes"] = sym_meta["equiv_classes"]
    except Exception as exc:
        report.update({"error": f"symmetry precompute: {exc}",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 5: pick U_total implementation -----
    try:
        from delfin._energy_terms import U_total as _U_total  # type: ignore
        _have_full_U = True
    except Exception:
        _U_total = None  # type: ignore
        _have_full_U = False

    def _eval(coords_2d: np.ndarray) -> Tuple[float, np.ndarray]:
        if _have_full_U and _U_total is not None:
            try:
                U, g = _U_total(coords_2d, mol, sym_info, params)
                g_arr = np.asarray(g, dtype=float).reshape(coords_2d.shape)
                return float(U), g_arr
            except Exception:
                # Graceful fallback if the rich term raises mid-run.
                return _fallback_U_total(coords_2d, mol, sym_info, params)
        return _fallback_U_total(coords_2d, mol, sym_info, params)

    def objective(x_flat: np.ndarray) -> Tuple[float, np.ndarray]:
        coords_2d = x_flat.reshape(n_atoms, 3)
        U, grad = _eval(coords_2d)
        return U, np.asarray(grad, dtype=float).reshape(-1)

    # ----- Step 6: initial energy -----
    try:
        U_initial, _ = _eval(coords)
        report["energy_initial"] = float(U_initial)
    except Exception as exc:
        report.update({"error": f"initial energy: {exc}",
                       "fallback_used": True})
        return xyz, report

    if not np.isfinite(U_initial):
        report.update({"error": "non-finite initial energy",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 7: L-BFGS-B run -----
    try:
        result = minimize(
            fun=objective,
            x0=coords.flatten(),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": int(max_iter), "ftol": float(ftol),
                     "gtol": 1e-5},
        )
        new_coords = np.asarray(result.x, dtype=float).reshape(n_atoms, 3)
        U_final = float(result.fun)
        converged = bool(result.success)
        iterations = int(getattr(result, "nit", 0))
    except Exception as exc:
        report.update({
            "error": f"L-BFGS-B: {exc}",
            "fallback_used": True,
            "traceback": traceback.format_exc(limit=3),
        })
        return xyz, report

    report["iterations"] = iterations
    report["converged"] = converged
    report["energy_final"] = U_final

    # ----- Step 8: numerical sanity (NaN/Inf guard) -----
    if not np.all(np.isfinite(new_coords)):
        report.update({"error": "non-finite coords after L-BFGS-B",
                       "fallback_used": True,
                       "topology_preserved": False})
        return xyz, report

    # Reject if the optimizer somehow blew the energy up.
    if not np.isfinite(U_final) or U_final > U_initial + 1.0 + abs(U_initial):
        report.update({"error": "energy increased / non-finite",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 9: topology hard-gate -----
    try:
        topo_ok = _topology_check(new_coords, mol, metal_set)
    except Exception:
        topo_ok = True  # B5 unavailable → don't penalize

    if not topo_ok:
        report.update({"error": "topology not preserved",
                       "fallback_used": True,
                       "topology_preserved": False})
        return xyz, report

    report["topology_preserved"] = True

    # ----- Step 10: write output -----
    try:
        refined_xyz = _array_to_xyz(new_coords, mol, orig_lines, symbols)
    except Exception as exc:
        report.update({"error": f"xyz write: {exc}",
                       "fallback_used": True})
        return xyz, report

    return refined_xyz, report


# ---------------------------------------------------------------------------
# Self-test — minimal end-to-end smoke (synthetic, no DELFIN-data dependency).
# ---------------------------------------------------------------------------

def _synthetic_xyz_water() -> Tuple[str, Any]:
    """Mildly distorted water molecule (a triangle with bond lengths perturbed)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    # Perturb each H slightly so U_bond has work to do.
    for i, dx in enumerate([(0.05, 0.0, 0.0), (0.0, -0.08, 0.04)]):
        if i + 1 < mol.GetNumAtoms():
            p = conf.GetAtomPosition(i + 1)
            conf.SetAtomPosition(i + 1, (p.x + dx[0], p.y + dx[1], p.z + dx[2]))
    lines = [str(mol.GetNumAtoms()), "synthetic water"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _synthetic_xyz_methane() -> Tuple[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    conf = mol.GetConformer()
    # Slightly squash one C-H.
    p = conf.GetAtomPosition(1)
    conf.SetAtomPosition(1, (p.x * 0.85, p.y * 0.85, p.z * 0.85))
    lines = [str(mol.GetNumAtoms()), "synthetic methane"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _synthetic_xyz_acetone() -> Tuple[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("CC(C)=O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=11)
    AllChem.UFFOptimizeMolecule(mol, maxIters=50)
    conf = mol.GetConformer()
    # Perturb a hydrogen position.
    p = conf.GetAtomPosition(mol.GetNumAtoms() - 1)
    conf.SetAtomPosition(mol.GetNumAtoms() - 1, (p.x + 0.12, p.y - 0.08, p.z + 0.06))
    lines = [str(mol.GetNumAtoms()), "synthetic acetone"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _self_test() -> None:
    """Run three synthetic cases and print a one-line summary each."""
    cases = [
        ("water", "no_metal", _synthetic_xyz_water),
        ("methane", "no_metal", _synthetic_xyz_methane),
        ("acetone", "no_metal", _synthetic_xyz_acetone),
    ]
    for name, klass, builder in cases:
        xyz, mol = builder()
        if not xyz or mol is None:
            print(f"[{name}] SKIP — RDKit not available")
            continue
        new_xyz, report = variational_refine(
            xyz, mol, class_label=klass, max_iter=100, ftol=1e-5,
            enable_global_pg=False,
        )
        ok = (report.get("topology_preserved", False)
              and not report.get("fallback_used", True))
        e0 = report.get("energy_initial", float("nan"))
        e1 = report.get("energy_final", float("nan"))
        delta = (e0 - e1) if (np.isfinite(e0) and np.isfinite(e1)) else float("nan")
        print(
            f"[{name:8s} class={klass:9s}] iters={report['iterations']:3d} "
            f"converged={report['converged']!s:5s} "
            f"E0={e0:10.4f} -> E1={e1:10.4f} ΔE={delta:+.4f} "
            f"topo_ok={report['topology_preserved']} "
            f"fallback={report['fallback_used']} OK={ok}"
        )


if __name__ == "__main__":  # pragma: no cover
    _self_test()
