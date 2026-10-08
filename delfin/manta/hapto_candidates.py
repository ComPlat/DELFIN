"""Hapto candidate topology check, quality scores, RDKit UFF refinement and best-candidate selection of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _mol_to_xyz,
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    _hapto_candidate_collapsed_bonds,
    _hapto_seat_rigid_enabled,
)
from delfin.manta.geometry_quality import (
    _enforce_metal_topology,
    _geometry_quality_score,
    _has_bad_geometry,
    _has_ligand_intertwining,
)
from delfin.manta.hybrid_assembly import (
    _enforce_donor_pi_coplanarity,
    _final_clash_resolution,
    _fix_secondary_metal_distances,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    RDKIT_AVAILABLE,
    _METAL_SET,
    _target_mc_dist,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.topology_checks import (
    _fragment_topology_ok,
    _fragment_topology_relaxed_fallback_ok,
    _has_atom_clash,
    _has_severe_covalent_distortion,
    _has_unphysical_metal_nonbonded_contact,
    _heavy_graph_exact_match_ok,
    _heavy_local_signature_match_ok,
    _metal_aware_coordination_ok,
    _no_spurious_bonds,
    _roundtrip_ring_count_ok,
)

logger = get_logger("delfin.smiles_converter")


def _hapto_candidate_topology_ok(
    xyz_delfin: str,
    original_smiles: str,
    mol=None,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> bool:
    """Return True when a hapto candidate preserves the input topology."""
    try:
        if mol is not None and not _metal_aware_coordination_ok(mol, conf_id, hapto_groups):
            return False
        if not _heavy_graph_exact_match_ok(xyz_delfin, original_smiles):
            return False
        if not _heavy_local_signature_match_ok(xyz_delfin, original_smiles):
            return False
        if not _roundtrip_ring_count_ok(xyz_delfin, original_smiles):
            return False
        if not _no_spurious_bonds(xyz_delfin, original_smiles):
            return False
        if not _fragment_topology_ok(xyz_delfin, original_smiles):
            if not _fragment_topology_relaxed_fallback_ok(xyz_delfin, original_smiles):
                return False
        return True
    except Exception:
        # If topology perception is unavailable (for example without OB),
        # keep the candidate and let graph/geometry checks decide.
        return True


def _hapto_geometry_quality_score(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> float:
    """Return a hapto-specific geometry penalty (lower = better)."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return 0.0
    try:
        import numpy as np
        conf = mol.GetConformer(conf_id)
    except Exception:
        return 0.0

    def _metal_atomic_number(sym: str) -> int:
        try:
            if RDKIT_AVAILABLE:
                return int(Chem.GetPeriodicTable().GetAtomicNumber(sym))
        except Exception:
            pass
        return 0

    def _is_f_block_like(sym: str) -> bool:
        z = _metal_atomic_number(sym)
        return (
            57 <= z <= 71
            or 89 <= z <= 103
            or sym in {'Y', 'Sc'}
        )

    def _hapto_weight_profile(metal_sym: str, eta: int, n_groups_for_metal: int) -> Dict[str, float]:
        profile = {
            'mc': 180.0,
            'radial': 30.0,
            'plane': 120.0,
            'axis': 40.0,
            'lateral': 35.0,
            'min_pair_angle': 90.0,
            'pair_sep': 0.45,
            'equal_eta_mc_std': 40.0,
        }

        if eta >= 5:
            profile.update({
                'mc': 220.0,
                'radial': 36.0,
                'plane': 165.0,
                'axis': 72.0,
                'lateral': 55.0,
                'min_pair_angle': 118.0,
                'pair_sep': 0.80,
                'equal_eta_mc_std': 55.0,
            })
        elif eta == 4:
            profile.update({
                'mc': 195.0,
                'radial': 32.0,
                'plane': 135.0,
                'axis': 52.0,
                'lateral': 42.0,
                'min_pair_angle': 102.0,
                'pair_sep': 0.60,
                'equal_eta_mc_std': 45.0,
            })
        elif eta == 3:
            profile.update({
                'mc': 150.0,
                'radial': 24.0,
                'plane': 82.0,
                'axis': 24.0,
                'lateral': 20.0,
                'min_pair_angle': 74.0,
                'pair_sep': 0.28,
                'equal_eta_mc_std': 28.0,
            })
        else:
            profile.update({
                'mc': 135.0,
                'radial': 18.0,
                'plane': 55.0,
                'axis': 14.0,
                'lateral': 12.0,
                'min_pair_angle': 60.0,
                'pair_sep': 0.18,
                'equal_eta_mc_std': 18.0,
            })

        if _is_f_block_like(metal_sym):
            profile['mc'] *= 0.80
            profile['radial'] *= 0.75
            profile['plane'] *= 0.45
            profile['axis'] *= 0.35
            profile['lateral'] *= 0.35
            profile['min_pair_angle'] = min(profile['min_pair_angle'], 72.0 if eta >= 5 else 58.0)
            profile['pair_sep'] *= 0.40
            profile['equal_eta_mc_std'] *= 0.70
        elif metal_sym in {'Fe', 'Co', 'Ni', 'Ru', 'Rh', 'Ir', 'Os', 'Mo', 'W', 'Re'} and eta >= 5:
            profile['plane'] *= 1.10
            profile['axis'] *= 1.20
            profile['lateral'] *= 1.15

        if n_groups_for_metal >= 2 and eta >= 4 and not _is_f_block_like(metal_sym):
            profile['axis'] *= 1.10
            profile['lateral'] *= 1.10
            profile['pair_sep'] *= 1.15

        return profile

    penalty = 0.0
    groups_by_metal: Dict[int, List[Tuple[int, List[int]]]] = {}
    for metal_idx, group_atoms in hapto_groups:
        groups_by_metal.setdefault(int(metal_idx), []).append((int(metal_idx), list(group_atoms)))

    centroid_records: Dict[int, List[Tuple[int, np.ndarray, np.ndarray, float, Dict[str, float]]]] = {}
    for metal_idx, group_atoms in hapto_groups:
        if not group_atoms:
            continue
        try:
            metal_idx = int(metal_idx)
            metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            profile = _hapto_weight_profile(metal_sym, len(group_atoms), len(groups_by_metal.get(metal_idx, [])))
            metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
            pts = np.array([conf.GetAtomPosition(int(atom_idx)) for atom_idx in group_atoms], dtype=float)
        except Exception:
            continue
        if pts.ndim != 2 or pts.shape[0] < 3:
            continue

        centroid = pts.mean(axis=0)
        mc_dist = float(np.linalg.norm(metal_pos - centroid))
        target_mc = float(_target_mc_dist(metal_sym, len(group_atoms)))
        penalty += profile['mc'] * (mc_dist - target_mc) ** 2

        mc_atom_dists = np.linalg.norm(pts - metal_pos, axis=1)
        if mc_atom_dists.size:
            penalty += profile['radial'] * float(np.std(mc_atom_dists) ** 2)

        q = pts - centroid
        try:
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
        except Exception:
            continue
        normal = vh[-1]
        n_norm = float(np.linalg.norm(normal))
        if n_norm < 1e-12:
            continue
        normal = normal / n_norm
        plane_dev = np.abs(q @ normal)
        plane_rms = float(np.sqrt(np.mean(plane_dev * plane_dev)))
        penalty += profile['plane'] * plane_rms * plane_rms

        axis = metal_pos - centroid
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis = axis / axis_norm
            cosang = abs(float(np.dot(normal, axis)))
            penalty += profile['axis'] * (1.0 - cosang) ** 2
            lateral_offset = float(np.linalg.norm((metal_pos - centroid) - np.dot(metal_pos - centroid, normal) * normal))
            penalty += profile['lateral'] * lateral_offset * lateral_offset
            centroid_records.setdefault(metal_idx, []).append(
                (len(group_atoms), centroid, axis, mc_dist, profile)
            )

    for metal_idx, records in centroid_records.items():
        if len(records) < 2:
            continue
        for i in range(len(records)):
            eta_i, centroid_i, axis_i, mc_i, profile_i = records[i]
            for j in range(i + 1, len(records)):
                eta_j, centroid_j, axis_j, mc_j, profile_j = records[j]
                cosang = max(-1.0, min(1.0, float(np.dot(axis_i, axis_j))))
                angle = math.degrees(math.acos(cosang))
                min_pair_angle = min(profile_i['min_pair_angle'], profile_j['min_pair_angle'])
                if angle < min_pair_angle:
                    gap = (min_pair_angle - angle) / max(min_pair_angle, 1.0)
                    penalty += 220.0 * min(profile_i['pair_sep'], profile_j['pair_sep']) * gap * gap
                if eta_i == eta_j:
                    penalty += min(profile_i['equal_eta_mc_std'], profile_j['equal_eta_mc_std']) * (mc_i - mc_j) ** 2
    return penalty


def _hapto_candidate_quality_score(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> float:
    """Return a penalty score for a hapto candidate (lower = better)."""
    if not RDKIT_AVAILABLE or mol is None:
        return float("inf")

    score = 0.0
    try:
        score += float(_geometry_quality_score(mol, conf_id))
    except Exception:
        score += 1.0e6

    try:
        score += float(_hapto_geometry_quality_score(mol, conf_id, hapto_groups))
    except Exception:
        score += 1.0e5

    try:
        if _has_severe_covalent_distortion(mol, conf_id):
            score += 1.0e6
    except Exception:
        score += 1.0e5

    try:
        if _has_bad_geometry(mol, conf_id):
            score += 2.5e5
    except Exception:
        score += 1.0e5

    try:
        if _has_atom_clash(mol, conf_id, min_dist=0.80):
            score += 2.5e5
    except Exception:
        pass

    try:
        if _has_unphysical_metal_nonbonded_contact(mol, conf_id):
            score += 2.0e5
    except Exception:
        pass

    try:
        if _has_ligand_intertwining(mol, conf_id):
            score += 8.0e4
    except Exception:
        pass

    return score


def _hapto_mol_from_xyz_template(mol_template, xyz_delfin: str):
    """Map DELFIN XYZ coordinates back onto a template molecule."""
    if not RDKIT_AVAILABLE or mol_template is None:
        return None
    try:
        mol_tmp = Chem.Mol(mol_template)
        mol_tmp.RemoveAllConformers()
        conf = _xyz_to_rdkit_conformer(mol_tmp, xyz_delfin)
        if conf is None:
            return None
        mol_tmp.AddConformer(conf, assignId=True)
        return mol_tmp
    except Exception:
        return None


def _refine_hapto_candidate_with_rdkit_uff(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
):
    """Run local RDKit UFF while freezing the hapto scaffold."""
    if not RDKIT_AVAILABLE or mol is None:
        return None
    try:
        work = Chem.Mol(mol)
        if work.GetNumConformers() == 0:
            return None
        ff = AllChem.UFFGetMoleculeForceField(work, confId=0)
        if ff is None:
            return None

        fixed: set = set()
        for atom in work.GetAtoms():
            if atom.GetSymbol() in _METAL_SET:
                fixed.add(atom.GetIdx())
        for metal_idx, catoms in hapto_groups:
            fixed.add(metal_idx)
            fixed.update(catoms)
        for atom_idx in sorted(fixed):
            ff.AddFixedPoint(int(atom_idx))

        ff.Minimize(maxIts=600)
        try:
            _enforce_donor_pi_coplanarity(work, 0, hapto_groups)
        except Exception:
            pass
        try:
            _fix_secondary_metal_distances(work, 0, hapto_groups)
        except Exception:
            pass
        try:
            _final_clash_resolution(work, 0, hapto_groups)
        except Exception:
            pass
        return work
    except Exception:
        return None


def _select_best_hapto_candidate(
    smiles: str,
    hapto_groups: List[Tuple[int, List[int]]],
    candidates: List[Tuple[str, object]],
    *,
    apply_uff: bool,
):
    """Select the best topology-preserving hapto candidate."""
    accepted: List[Tuple[float, object, str]] = []
    relaxed: List[Tuple[float, object, str]] = []

    # Iter-26b (Task #14, hapto): the OpenBabel-UFF refinement leg is harmful
    # for hapto/multi_hapto on unparametrized 4d/5d metals — it cannot type the
    # metal cation (logged "Unrecognized atom type"), folds ring substituents /
    # donors together (controlled forensic: clean builder geometry -> 6 H-H
    # superpositions at 0.0-0.02 A after OB-UFF; Fe-Cp-imine even SIGSEGVs).
    # When DELFIN_HAPTO_NO_OB_UFF=1, suppress the OB-UFF leg and keep ONLY the
    # metal+eta-frozen RDKit-UFF refinement.  This function only runs on the
    # hapto path, so a plain env flag is class-correct.  Default OFF until the
    # CCDC-calibrated metric + smoke validate it.
    _suppress_ob_uff = os.environ.get("DELFIN_HAPTO_NO_OB_UFF", "0") == "1"

    for label, candidate_mol in candidates:
        if candidate_mol is None or candidate_mol.GetNumConformers() == 0:
            continue

        try:
            xyz_candidate = _mol_to_xyz(candidate_mol)
        except Exception:
            continue

        topo_ok = _hapto_candidate_topology_ok(
            xyz_candidate,
            smiles,
            mol=candidate_mol,
            conf_id=0,
            hapto_groups=hapto_groups,
        )
        base_score = _hapto_candidate_quality_score(candidate_mol, 0, hapto_groups)
        target_bucket = accepted if topo_ok else relaxed
        target_bucket.append((base_score, candidate_mol, label))

        if not apply_uff:
            continue

        refined_trials: List[Tuple[str, object]] = []
        rdkit_refined = _refine_hapto_candidate_with_rdkit_uff(candidate_mol, hapto_groups)
        if rdkit_refined is not None:
            refined_trials.append((f"{label}+rdkit-uff", rdkit_refined))

        if not _suppress_ob_uff:
            try:
                xyz_refined = _optimize_xyz_openbabel_safe(
                    xyz_candidate,
                    mol_template=candidate_mol,
                    smiles=smiles,
                    steps=750,
                    apply_template_constraints=True,
                )
            except Exception:
                xyz_refined = xyz_candidate
            if xyz_refined and xyz_refined != xyz_candidate:
                ob_refined = _hapto_mol_from_xyz_template(candidate_mol, xyz_refined)
                if ob_refined is not None:
                    refined_trials.append((f"{label}+ob-uff", ob_refined))

        for refined_label, refined_mol in refined_trials:
            try:
                xyz_refined = _mol_to_xyz(refined_mol)
            except Exception:
                continue
            if not _hapto_candidate_topology_ok(
                xyz_refined,
                smiles,
                mol=refined_mol,
                conf_id=0,
                hapto_groups=hapto_groups,
            ):
                continue
            refined_score = _hapto_candidate_quality_score(refined_mol, 0, hapto_groups)
            if refined_score + 1e-6 < base_score:
                accepted.append((refined_score, refined_mol, refined_label))

    pool = accepted if accepted else relaxed
    if not pool:
        return None
    # Iter-26b (Task #14, hapto): prefer the clean RIGID scaffold when its
    # topology is OK.  The analytical _build_hapto_scaffold places the eta-ring
    # + substituents + donors geometrically (forensic: ~0 superpositions), but
    # it loses the score-based selection to ETKDG/UFF candidates that fold
    # substituents together.  When DELFIN_HAPTO_PREFER_SCAFFOLD=1, return the RAW
    # scaffold (NOT a UFF-refined variant) if it is in the topology-OK pool;
    # otherwise fall through to normal scoring (safe — never picks a
    # topology-broken scaffold).  Default OFF until CCDC-metric + smoke validate.
    if os.environ.get("DELFIN_HAPTO_PREFER_SCAFFOLD", "0") == "1" and accepted:
        scaffold_entries = [it for it in accepted if it[2] == "scaffold"]
        if scaffold_entries:
            best_score, best_mol, best_label = scaffold_entries[0]
            logger.info(
                "Hapto: preferring raw scaffold candidate (score=%.2f, pool=%d)",
                best_score, len(pool),
            )
            _enforce_metal_topology(best_mol, conf_id=0)
            return best_mol
    pool.sort(key=lambda item: (item[0], item[2]))
    best_score, best_mol, best_label = pool[0]

    # ---- THE RIGID SEATING MAY ARRIVE (DELFIN_FFFREE_HAPTO_SEAT_RIGID) ----------
    # Default 0 -> this block is a single `if` and changes no byte.
    #
    # The score knows no collapse term (rationale and numbers are at
    # `_hapto_seat_rigid_enabled`).  If the ATOM-WISE build (`scaffold`) wins,
    # although a RIGID build (`hybrid*`) stands in the SAME pool and carries STRICTLY
    # fewer collapsed bonds, then take the rigid one.
    #
    # Three self-restrictions, so this does not become a blank cheque:
    #   1. only WITHIN the same pool -- a topologically broken candidate can
    #      thus never displace a topologically sound one.
    #   2. only with STRICTLY less collapse.  BAKLAB (rigid 54 : atom-wise 51) thereby
    #      stays with the atom-wise build on its own -- the rule limits itself.
    #   3. `None` (not measurable) is NO verdict and leaves the selection untouched.
    if _hapto_seat_rigid_enabled() and str(best_label).startswith("scaffold"):
        _n_atomwise = _hapto_candidate_collapsed_bonds(best_mol)
        if _n_atomwise is not None:
            _bester_starr = None
            for _sc, _mo, _lb in pool:
                if not str(_lb).startswith("hybrid"):
                    continue
                _n_starr = _hapto_candidate_collapsed_bonds(_mo)
                if _n_starr is None or _n_starr >= _n_atomwise:
                    continue
                if _bester_starr is None or _n_starr < _bester_starr[0]:
                    _bester_starr = (_n_starr, _sc, _mo, _lb)
            if _bester_starr is not None:
                _n_starr, best_score, best_mol, best_label = _bester_starr
                logger.info(
                    "HAPTO_SEAT_RIGID: starrer Sitz %s statt atomweisem %s "
                    "(Kollaps %d statt %d, Score %.2f statt %.2f)",
                    best_label, pool[0][2], _n_starr, _n_atomwise,
                    best_score, pool[0][0],
                )

    logger.info(
        "Selected hapto candidate %s (score=%.2f, strict_topology=%s, pool=%d)",
        best_label,
        best_score,
        bool(accepted),
        len(pool),
    )
    _enforce_metal_topology(best_mol, conf_id=0)
    return best_mol
