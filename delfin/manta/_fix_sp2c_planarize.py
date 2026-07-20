"""SP2C-PLANARIZE — universal sp2-CARBON planarisation fixer.

Sibling of ``_fix_sp2n_planarize`` (which handles sp2 NITROGEN) for the exact
gap it left open: an sp2 CARBON — an acyclic azomethine / imine / vinyl / enone
``=CH-`` or ``=CR-`` — that the force-field-free ETKDG embed has left
**pyramidalised** (the carbon pushed out of the plane of its three substituents;
angle-sum drifts toward the sp3 fallback ~328° instead of the planar-sp2 360°).

WHY this is needed: the MANTA manifold is force-field-free by default, so the
planarity of a NON-ring conjugated sp2 carbon rests only on ETKDGv3's *soft*
basic-knowledge terms (a different per-fold seed => a different residual
out-of-plane), and the one relaxation that would harden it (MMFF's out-of-plane
term) is disabled on the chelate path.  The ring aromatic passes
(``_arom_planarize``) only touch ring atoms; ``_fix_sp2n_planarize`` is
N-only (``if atom.GetSymbol() != "N": continue``).  So a backbone ``N=CH-C``
carbon is planarised by NOTHING — observed as NAYKOQ C36 pyramidal (Walsh ~19°,
angle-sum 324°) in 17 of 30 folds while 13 land flat.

Touch-rules (chemical safety, mirrors the sp2n / F25 doctrine, UNIVERSAL — pure
graph + element rules, never SMILES-specific):
  - Only a 3-coordinate carbon that RDKit perceives sp2 AND carries a double bond
    to N/C/O (the robust pi marker) — a genuine planar centre.
  - Skip ring C (aromatic/ring planarity is owned by the ring passes).
  - Skip a metal-bonded C (a carbene / sigma-alkyl donor obeys coordination).
  - Geometry-only: project the carbon onto the best-fit plane of its three
    neighbours (bond lengths shift only by the small residual; no retarget).
  - Per-atom rollback if the projection introduces a NEW non-bonded clash OR does
    not actually improve planarity OR yields a non-finite coordinate (never-worse).
  - Deterministic: closed-form projection, no random state, no FF call.

Reuses the sp2n module's vendored helpers (parse/format/oop/clash/metal) so the
two fixers stay bit-for-bit consistent.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

# Share the sp2n module's helpers verbatim (DRY; keeps the two fixers consistent).
from delfin.manta._fix_sp2n_planarize import (
    _parse_xyz, _format_xyz, _oop_distance, _introduces_new_clash, _is_metal_sym,
)


def detect_planar_sp2c_groups(mol) -> List[Dict]:
    """Every acyclic, non-metal-bonded, 3-coordinate sp2 CARBON that must be planar.

    Each entry: {"c_idx", "nbr_idxs": [i, j, k]}.  Requires RDKit hybridization SP2
    AND a double bond from the carbon to N/C/O (imine/azomethine/vinyl/enone) — the
    robust pi marker.  Ring carbons are skipped (owned by the aromatic ring passes),
    as are metal-bonded carbons (carbene / sigma-donor obey coordination geometry)."""
    try:
        from rdkit import Chem
    except Exception:
        return []
    out: List[Dict] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        if atom.IsInRing():
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        if any(_is_metal_sym(nb.GetSymbol()) for nb in nbrs):
            continue
        try:
            if atom.GetHybridization() != Chem.HybridizationType.SP2:
                continue
        except Exception:
            continue
        # pi evidence: a genuine double bond from this carbon to N/C/O.
        has_double = False
        for b in atom.GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE and \
                    b.GetOtherAtom(atom).GetSymbol() in ("N", "C", "O"):
                has_double = True
                break
        if not has_double:
            continue
        out.append({"c_idx": atom.GetIdx(),
                    "nbr_idxs": [nb.GetIdx() for nb in nbrs]})
    return out


def planarize_sp2_carbon(xyz: str, mol,
                         oop_threshold_A: float = 0.20,
                         ) -> Tuple[str, Dict]:
    """Flatten pyramidalised acyclic sp2 CARBON.

    For each detected sp2 carbon whose distance from the plane of its three
    neighbours exceeds ``oop_threshold_A``, project the carbon into that plane
    (geometry-only; bond lengths shift only by the small residual).  Per-atom
    rollback on new clash / no improvement / non-finite output (never-worse).

    Returns ``(new_xyz, report)``.  report keys: ``n_candidates``, ``n_fixed``,
    ``max_oop_before_A``, ``per_c``."""
    syms, pts, orig_lines = _parse_xyz(xyz)
    report: Dict = {"n_candidates": 0, "n_fixed": 0,
                    "max_oop_before_A": 0.0, "per_c": []}
    if pts.shape[0] == 0 or mol is None:
        return xyz, report
    if mol.GetNumAtoms() != len(syms):
        return xyz, report

    try:
        groups = detect_planar_sp2c_groups(mol)
    except Exception:
        return xyz, report
    report["n_candidates"] = len(groups)
    if not groups:
        return xyz, report

    new_pts = pts.copy()

    for g in groups:
        c_idx = g["c_idx"]
        nbr = g["nbr_idxs"]
        try:
            if len(nbr) != 3:
                continue
            p_c = new_pts[c_idx]
            p_a, p_b, p_d = (new_pts[nbr[0]], new_pts[nbr[1]], new_pts[nbr[2]])
            oop = _oop_distance(p_c, p_a, p_b, p_d)
            if oop is None:
                continue
            report["max_oop_before_A"] = max(report["max_oop_before_A"], float(oop))
            if oop <= oop_threshold_A:
                report["per_c"].append({"c_idx": c_idx, "fixed": False,
                                        "reason": "within_tolerance",
                                        "oop_A": round(float(oop), 3)})
                continue
            nrm = np.cross(p_b - p_a, p_d - p_a)
            nn = float(np.linalg.norm(nrm))
            if nn < 1e-6:
                continue
            normal = nrm / nn
            # Drop the carbon's out-of-plane component (project onto the plane).
            disp = float(np.dot(p_c - p_a, normal))
            new_c = p_c - disp * normal
            cand = new_pts.copy()
            cand[c_idx] = new_c
            if not np.all(np.isfinite(cand[c_idx])):
                continue
            if _introduces_new_clash(new_pts, cand, syms, [c_idx]):
                report["per_c"].append({"c_idx": c_idx, "fixed": False,
                                        "reason": "rollback_clash"})
                continue
            new_oop = _oop_distance(cand[c_idx], cand[nbr[0]], cand[nbr[1]],
                                    cand[nbr[2]])
            if new_oop is not None and new_oop > float(oop) + 1e-6:
                report["per_c"].append({"c_idx": c_idx, "fixed": False,
                                        "reason": "no_improvement"})
                continue
            new_pts = cand
            report["n_fixed"] += 1
            report["per_c"].append({
                "c_idx": c_idx, "fixed": True,
                "oop_before_A": round(float(oop), 3),
                "oop_after_A": round(float(new_oop), 3)
                if new_oop is not None else None})
        except Exception:
            continue

    if report["n_fixed"] == 0:
        return xyz, report
    if not np.all(np.isfinite(new_pts)):
        return xyz, report
    return _format_xyz(orig_lines, syms, new_pts), report
