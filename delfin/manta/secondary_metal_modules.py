"""Secondary-metal coordination modules of multi-metal hapto complexes: geometry fits, pose optimisation and module assembly in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _mol_to_xyz,
)
from delfin.manta.converter_flags import (
    _all_polyhedra_codes,
)
from delfin.manta.hapto_scaffold import (
    _build_hapto_scaffold,
    _correct_hapto_geometry,
    _propagate_non_hapto_atoms,
)
from delfin.manta.hybrid_fragments import (
    _HybridHaptoDecomposition,
    _HybridHaptoFragment,
    _PrimaryOrganometalModule,
    _embed_hybrid_fragment,
    _hapto_primary_donor_indices,
    _primary_organometal_module_quality_ok,
    _rotation_matrix_from_vectors,
)
from delfin.manta.isomer_labels import (
    _TOPO_GEOMETRY_VECTORS,
)
from delfin.manta.ml_tables import (
    Chem,
    Point3D,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_SET,
    _PREFERRED_CN4_GEOMETRY,
    _get_ml_bond_length,
    _secondary_donor_fit_weight,
    _secondary_donor_target_length,
)

logger = get_logger("delfin.smiles_converter")


def _align_hybrid_fragment_to_targets(
    scaffold_mol,
    fragment: _HybridHaptoFragment,
    embedded_fragment,
    target_atom_indices: List[int],
    target_positions: Dict[int, object],
    reference_center=None,
    metal_idx: Optional[int] = None,
    exact_target_count: int = 0,
) -> bool:
    """Rigidly align a fragment to explicit target coordinates."""
    if not RDKIT_AVAILABLE or scaffold_mol is None or embedded_fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        scaffold_conf = scaffold_mol.GetConformer()
        fragment_conf = embedded_fragment.GetConformer()
    except Exception:
        return False

    anchor_pairs: List[Tuple[int, int]] = []
    seen_original: set = set()
    for original_idx in target_atom_indices:
        if original_idx in seen_original:
            continue
        seen_original.add(original_idx)
        fragment_idx = fragment.original_to_fragment.get(original_idx)
        if fragment_idx is None:
            continue
        target_value = target_positions.get(original_idx)
        if target_value is None:
            continue
        target_vec = np.asarray(target_value, dtype=float)
        if target_vec.shape != (3,) or not np.all(np.isfinite(target_vec)):
            continue
        anchor_pairs.append((fragment_idx, original_idx))

    if not anchor_pairs:
        return False

    src = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(fragment_idx).x,
                fragment_conf.GetAtomPosition(fragment_idx).y,
                fragment_conf.GetAtomPosition(fragment_idx).z,
            ]
            for fragment_idx, _original_idx in anchor_pairs
        ],
        dtype=float,
    )
    tgt = np.asarray(
        [np.asarray(target_positions[original_idx], dtype=float) for _fragment_idx, original_idx in anchor_pairs],
        dtype=float,
    )
    all_coords = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(i).x,
                fragment_conf.GetAtomPosition(i).y,
                fragment_conf.GetAtomPosition(i).z,
            ]
            for i in range(embedded_fragment.GetNumAtoms())
        ],
        dtype=float,
    )
    src_coord_map = {
        original_idx: all_coords[fragment_idx]
        for fragment_idx, original_idx in fragment.fragment_to_original.items()
    }
    src_planar_frame = _fragment_planar_donor_frame(
        scaffold_mol,
        fragment.atom_indices,
        target_atom_indices,
        src_coord_map,
    )

    if len(anchor_pairs) == 1:
        if reference_center is not None:
            ref = np.asarray(reference_center, dtype=float)
            frag_center = all_coords.mean(axis=0)
            if src_planar_frame is not None:
                src_dir = np.asarray(src_planar_frame["outward"], dtype=float)
                tgt_dir = ref - tgt[0]
            else:
                src_dir = frag_center - src[0]
                tgt_dir = tgt[0] - ref
            rot = _rotation_matrix_from_vectors(src_dir, tgt_dir)
            aligned = (all_coords - src[0]) @ rot.T + tgt[0]
        else:
            aligned = all_coords + (tgt[0] - src[0])
    elif len(anchor_pairs) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        src_mid = 0.5 * (src[0] + src[1])
        tgt_mid = 0.5 * (tgt[0] + tgt[1])
        aligned = (all_coords - src_mid) @ rot.T + tgt_mid
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_coords - src_centroid) @ rot.T + tgt_centroid

    if reference_center is not None and metal_idx is not None and len(anchor_pairs) in (1, 2):
        ref = np.asarray(reference_center, dtype=float)

        def _alignment_score(candidate_coords) -> float:
            donor_err = 0.0
            coord_map = {}
            donor_original_indices = []
            for pair_idx, (fragment_idx, original_idx) in enumerate(anchor_pairs):
                donor_original_indices.append(original_idx)
                donor_err += float(np.sum((candidate_coords[fragment_idx] - tgt[pair_idx]) ** 2))
            for fragment_idx, original_idx in fragment.fragment_to_original.items():
                coord_map[original_idx] = np.asarray(candidate_coords[fragment_idx], dtype=float)
            return (
                25.0 * donor_err
                + _secondary_non_donor_contact_penalty(
                    scaffold_mol,
                    metal_idx,
                    ref,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                )
                + _fragment_donor_selectivity_penalty(
                    scaffold_mol,
                    metal_idx,
                    ref,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                )
                + _fragment_environment_clash_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    coord_map,
                    ignore_atom_indices={metal_idx},
                    primary_metal_idx=metal_idx,
                )
                + _planar_fragment_metal_coplanarity_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    ref,
                    coord_map,
                )
                + _planar_fragment_donor_approach_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    ref,
                    coord_map,
                )
                + _fragment_bite_direction_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                    ref,
                )
            )

        aligned_coord_map = {
            original_idx: np.asarray(aligned[fragment_idx], dtype=float)
            for fragment_idx, original_idx in fragment.fragment_to_original.items()
        }
        aligned_planar_frame = _fragment_planar_donor_frame(
            scaffold_mol,
            fragment.atom_indices,
            [original_idx for _fragment_idx, original_idx in anchor_pairs],
            aligned_coord_map,
        )

        if len(anchor_pairs) == 1:
            pivot = tgt[0]
            axis = ref - tgt[0]
        else:
            pivot = tgt[0]
            axis = tgt[1] - tgt[0]
        axis_candidates = []
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis_candidates.append(axis / axis_norm)
        if aligned_planar_frame is not None:
            normal = np.asarray(aligned_planar_frame["normal"], dtype=float)
            normal_norm = float(np.linalg.norm(normal))
            if normal_norm > 1e-12:
                axis_candidates.append(normal / normal_norm)
            outward = np.asarray(aligned_planar_frame["outward"], dtype=float)
            outward_norm = float(np.linalg.norm(outward))
            if outward_norm > 1e-12:
                axis_candidates.append(outward / outward_norm)
            if axis_candidates:
                cross_axis = np.cross(axis_candidates[0], normal if normal_norm > 1e-12 else outward)
                cross_norm = float(np.linalg.norm(cross_axis))
                if cross_norm > 1e-12:
                    axis_candidates.append(cross_axis / cross_norm)

        if axis_candidates:
            best_aligned = aligned.copy()
            best_score = _alignment_score(best_aligned)
            step = 5 if aligned_planar_frame is not None else 15
            for axis in axis_candidates:
                for angle_deg in range(-180, 181, step):
                    if angle_deg == 0:
                        continue
                    rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
                    trial = (aligned - pivot) @ rot.T + pivot
                    trial_score = _alignment_score(trial)
                    if trial_score + 1e-9 < best_score:
                        best_score = trial_score
                        best_aligned = trial
            aligned = best_aligned

    for pair_idx, (fragment_idx, _original_idx) in enumerate(anchor_pairs[:max(0, exact_target_count)]):
        aligned[fragment_idx] = tgt[pair_idx]

    for fragment_idx, original_idx in fragment.fragment_to_original.items():
        scaffold_conf.SetAtomPosition(
            original_idx,
            Point3D(
                float(aligned[fragment_idx, 0]),
                float(aligned[fragment_idx, 1]),
                float(aligned[fragment_idx, 2]),
            ),
        )
    return True


def _secondary_metal_geometry_codes(n_coord: int, metal_symbol: str) -> List[str]:
    """Return geometry candidates for explicit secondary-metal placement."""
    if n_coord <= 1:
        return []
    # Iter-2 Subagent A: legacy list (HEAD-bit-exact when DELFIN_ALL_POLYHEDRA=0)
    # passed as the seed `legacy_codes` so env=0 path is unmodified.  When
    # env=1, helper appends missing polyhedra (SS for CN4, TPR for CN6, etc.).
    if n_coord == 2:
        legacy = ['LIN']
    elif n_coord == 3:
        legacy = ['TP', 'TS']
    elif n_coord == 4:
        pref = _PREFERRED_CN4_GEOMETRY.get(metal_symbol, 'SQ')
        other = 'TH' if pref == 'SQ' else 'SQ'
        legacy = [pref, other]
    elif n_coord == 5:
        legacy = ['TBP', 'SP']
    elif n_coord == 6:
        legacy = ['OH']
    elif n_coord == 7:
        legacy = ['PBP']
    elif n_coord == 8:
        legacy = ['SAP', 'DD']
    else:
        return []
    return _all_polyhedra_codes(n_coord, metal_symbol, legacy)


def _fit_secondary_metal_model_to_targets(
    model_positions,
    target_positions,
    fit_indices: List[int],
):
    """Fit a full coordination model onto explicit donor targets."""
    import numpy as np

    model = np.asarray(model_positions, dtype=float)
    targets = np.asarray(target_positions, dtype=float)
    if model.ndim != 2 or targets.ndim != 2 or model.shape != targets.shape:
        return None
    if not fit_indices:
        return None

    model_fit = model[fit_indices]
    target_fit = targets[fit_indices]

    if len(fit_indices) == 1:
        delta = target_fit[0] - model_fit[0]
        transformed = model + delta
        metal_pos = delta
        return transformed, metal_pos

    if len(fit_indices) == 2:
        rot = _rotation_matrix_from_vectors(
            model_fit[1] - model_fit[0],
            target_fit[1] - target_fit[0],
        )
        model_mid = 0.5 * (model_fit[0] + model_fit[1])
        target_mid = 0.5 * (target_fit[0] + target_fit[1])
        transformed = (model - model_mid) @ rot.T + target_mid
        metal_pos = (-model_mid) @ rot.T + target_mid
        return transformed, metal_pos

    model_centroid = model_fit.mean(axis=0)
    target_centroid = target_fit.mean(axis=0)
    covariance = (model_fit - model_centroid).T @ (target_fit - target_centroid)
    try:
        u, _s, vt = np.linalg.svd(covariance)
    except Exception:
        return None
    rot = vt.T @ u.T
    if float(np.linalg.det(rot)) < 0.0:
        vt[-1, :] *= -1.0
        rot = vt.T @ u.T
    transformed = (model - model_centroid) @ rot.T + target_centroid
    metal_pos = (-model_centroid) @ rot.T + target_centroid
    return transformed, metal_pos


def _axis_angle_rotation_matrix(axis, angle_rad: float):
    """Return the rotation matrix for a rotation around ``axis``."""
    import numpy as np

    axis_arr = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(axis_arr))
    if norm < 1e-12:
        return np.eye(3)
    axis_arr /= norm
    ux, uy, uz = axis_arr
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array([
        [c + ux * ux * one_c, ux * uy * one_c - uz * s, ux * uz * one_c + uy * s],
        [uy * ux * one_c + uz * s, c + uy * uy * one_c, uy * uz * one_c - ux * s],
        [uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c],
    ], dtype=float)


def _project_displacements_to_rigid_body(
    positions: Dict[int, object],
    displacement_map: Dict[int, object],
    atom_indices: List[int],
) -> Dict[int, object]:
    """Approximate a set of atomic displacements by one rigid-body motion."""
    import numpy as np

    ordered = [atom_idx for atom_idx in atom_indices if atom_idx in positions]
    if len(ordered) < 2:
        return {
            atom_idx: np.asarray(
                displacement_map.get(atom_idx, np.zeros(3, dtype=float)),
                dtype=float,
            )
            for atom_idx in ordered
        }

    pts = np.asarray([positions[atom_idx] for atom_idx in ordered], dtype=float)
    disp = np.asarray(
        [displacement_map.get(atom_idx, np.zeros(3, dtype=float)) for atom_idx in ordered],
        dtype=float,
    )
    centroid = pts.mean(axis=0)

    a = np.zeros((3 * len(ordered), 6), dtype=float)
    b = disp.reshape(-1)
    ident = np.eye(3, dtype=float)
    for row_idx, atom_idx in enumerate(ordered):
        rel = np.asarray(positions[atom_idx], dtype=float) - centroid
        cross_block = np.array(
            [
                [0.0, rel[2], -rel[1]],
                [-rel[2], 0.0, rel[0]],
                [rel[1], -rel[0], 0.0],
            ],
            dtype=float,
        )
        start = 3 * row_idx
        a[start:start + 3, 0:3] = ident
        a[start:start + 3, 3:6] = cross_block

    try:
        sol, *_rest = np.linalg.lstsq(a, b, rcond=None)
    except Exception:
        return {
            atom_idx: np.asarray(
                displacement_map.get(atom_idx, np.zeros(3, dtype=float)),
                dtype=float,
            )
            for atom_idx in ordered
        }

    trans = np.asarray(sol[:3], dtype=float)
    omega = np.asarray(sol[3:6], dtype=float)
    projected: Dict[int, object] = {}
    for atom_idx in ordered:
        rel = np.asarray(positions[atom_idx], dtype=float) - centroid
        projected[atom_idx] = trans + np.cross(omega, rel)
    return projected


def _secondary_non_donor_contact_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
) -> float:
    """Penalty for non-donor atoms collapsing onto a secondary metal center."""
    import numpy as np

    if mol is None or not atom_indices:
        return 0.0

    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    donor_set = set(donor_indices)
    atom_set = set(atom_indices)
    alpha_atoms: set = set()
    beta_atoms: set = set()

    for donor_idx in donor_indices:
        atom = mol.GetAtomWithIdx(donor_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in atom_set and nbr_idx not in donor_set and nbr.GetAtomicNum() > 1:
                alpha_atoms.add(nbr_idx)
    for atom_idx in alpha_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in atom_set and nbr_idx not in donor_set and nbr_idx not in alpha_atoms and nbr.GetAtomicNum() > 1:
                beta_atoms.add(nbr_idx)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _min_allowed(atom_idx: int, atom_sym: str) -> float:
        ml_ref = float(_get_ml_bond_length(metal_sym, atom_sym))
        min_allowed = max(1.70, 0.92 * ml_ref)
        if atom_idx in alpha_atoms:
            min_allowed = max(min_allowed, 2.05)
        elif atom_idx in beta_atoms:
            min_allowed = max(min_allowed, 1.85)
        atom = mol.GetAtomWithIdx(atom_idx)
        if _is_planar_like(atom):
            if atom_idx in alpha_atoms:
                min_allowed = max(min_allowed, 2.18)
            elif atom_idx in beta_atoms:
                min_allowed = max(min_allowed, 1.98)
            else:
                min_allowed = max(min_allowed, 1.82)
        return min_allowed

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    penalty = 0.0
    for atom_idx in atom_indices:
        if atom_idx in donor_set:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        coord = coord_map.get(atom_idx)
        if coord is None:
            continue
        atom_pos = np.asarray(coord, dtype=float)
        dist = float(np.linalg.norm(atom_pos - metal_pos))
        min_allowed = _min_allowed(atom_idx, atom.GetSymbol())
        if dist < min_allowed:
            gap = min_allowed - dist
            if _is_planar_like(atom):
                weight = 40.0 if atom_idx in alpha_atoms else 26.0
            else:
                weight = 28.0 if atom_idx in alpha_atoms else 18.0
            penalty += weight * gap * gap
        if dist < 1.55:
            penalty += 120.0 * (1.55 - dist) ** 2
    return penalty


def _secondary_fragment_donor_selectivity_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    donor_indices: List[int],
    donor_to_fragment: Dict[int, _HybridHaptoFragment],
    coord_map: Dict[int, object],
) -> float:
    """Penalize non-donor atoms approaching a secondary metal as closely as the true donor."""
    import numpy as np

    if mol is None or not donor_indices or not coord_map:
        return 0.0

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    donor_set = set(donor_indices)
    penalty = 0.0

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _is_carbonyl_like_carbon(atom_idx: int) -> bool:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            return False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() < 1.5:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetSymbol() == 'O':
                return True
        return False

    for donor_idx in donor_indices:
        donor_pos = coord_map.get(donor_idx)
        fragment = donor_to_fragment.get(donor_idx)
        if donor_pos is None or fragment is None:
            continue
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = donor_atom.GetSymbol()
        donor_dist = float(np.linalg.norm(np.asarray(donor_pos, dtype=float) - metal_pos))

        for atom_idx in fragment.atom_indices:
            if atom_idx in donor_set or atom_idx == metal_idx:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = coord_map.get(atom_idx)
            if atom_pos is None:
                continue
            atom_dist = float(np.linalg.norm(np.asarray(atom_pos, dtype=float) - metal_pos))
            try:
                topo_sep = len(Chem.GetShortestPath(mol, donor_idx, atom_idx)) - 1
            except Exception:
                topo_sep = 99
            if topo_sep < 1 or topo_sep > 3:
                continue

            planar_like = _is_planar_like(atom)
            delta = 0.14
            weight = 10.0
            if topo_sep == 1:
                delta = 0.28 if donor_sym == 'O' else 0.32
                weight = 28.0 if donor_sym == 'O' else 34.0
            elif topo_sep == 2:
                delta = 0.20 if donor_sym == 'O' else 0.24
                weight = 18.0 if donor_sym == 'O' else 24.0
            elif topo_sep == 3:
                delta = 0.10 if donor_sym == 'O' else 0.12
                weight = 10.0 if donor_sym == 'O' else 12.0

            if planar_like:
                delta += 0.08
                weight += 10.0
            if donor_sym == 'N' and atom.GetSymbol() == 'C':
                delta += 0.06
                weight += 8.0
            if donor_sym == 'O' and atom.GetSymbol() == 'C':
                delta += 0.05
                weight += 6.0
                if _is_carbonyl_like_carbon(atom_idx):
                    delta += 0.16
                    weight += 18.0

            preferred_min = donor_dist + delta
            if atom_dist < preferred_min:
                gap = preferred_min - atom_dist
                penalty += weight * gap * gap
            if atom_dist < donor_dist + 0.02:
                gap = donor_dist + 0.02 - atom_dist
                penalty += (34.0 + weight) * gap * gap

    return penalty


def _fragment_donor_selectivity_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
) -> float:
    """Penalize non-donor atoms in one donor fragment approaching as closely as the true donor."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    donor_set = set(donor_indices)
    fragment_set = set(fragment_atom_indices)
    penalty = 0.0

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _is_carbonyl_like_carbon(atom_idx: int) -> bool:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            return False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() < 1.5:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetSymbol() == 'O':
                return True
        return False

    for donor_idx in donor_indices:
        donor_pos = coord_map.get(donor_idx)
        if donor_pos is None:
            continue
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = donor_atom.GetSymbol()
        donor_dist = float(np.linalg.norm(np.asarray(donor_pos, dtype=float) - metal_pos))
        for atom_idx in fragment_atom_indices:
            if atom_idx not in fragment_set or atom_idx in donor_set or atom_idx == metal_idx:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = coord_map.get(atom_idx)
            if atom_pos is None:
                continue
            atom_dist = float(np.linalg.norm(np.asarray(atom_pos, dtype=float) - metal_pos))
            try:
                topo_sep = len(Chem.GetShortestPath(mol, donor_idx, atom_idx)) - 1
            except Exception:
                topo_sep = 99
            if topo_sep < 1 or topo_sep > 3:
                continue

            planar_like = _is_planar_like(atom)
            delta = 0.14
            weight = 10.0
            if topo_sep == 1:
                delta = 0.28 if donor_sym == 'O' else 0.32
                weight = 28.0 if donor_sym == 'O' else 34.0
            elif topo_sep == 2:
                delta = 0.20 if donor_sym == 'O' else 0.24
                weight = 18.0 if donor_sym == 'O' else 24.0
            elif topo_sep == 3:
                delta = 0.10 if donor_sym == 'O' else 0.12
                weight = 10.0 if donor_sym == 'O' else 12.0

            if planar_like:
                delta += 0.08
                weight += 10.0
            if donor_sym == 'N' and atom.GetSymbol() == 'C':
                delta += 0.06
                weight += 8.0
            if donor_sym == 'O' and atom.GetSymbol() == 'C':
                delta += 0.05
                weight += 6.0
                if _is_carbonyl_like_carbon(atom_idx):
                    delta += 0.16
                    weight += 18.0

            preferred_min = donor_dist + delta
            if atom_dist < preferred_min:
                gap = preferred_min - atom_dist
                penalty += weight * gap * gap
            if atom_dist < donor_dist + 0.02:
                gap = donor_dist + 0.02 - atom_dist
                penalty += (34.0 + weight) * gap * gap

    return penalty


def _planar_fragment_metal_coplanarity_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    metal_pos,
    coord_map: Dict[int, object],
) -> float:
    """Penalty when a planar donor fragment is not approximately coplanar with the metal."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    donor_set = set(donor_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = tuple(
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    )
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in set(planar_atoms)
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return 0.0

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return 0.0
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return 0.0
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return 0.0
    normal = normal / normal_norm

    weight = 42.0 + 12.0 * len(donor_overlap)
    if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'N' for donor_idx in donor_overlap):
        weight += 14.0
    if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O' for donor_idx in donor_overlap):
        weight += 14.0
    if any(mol.GetAtomWithIdx(donor_idx).GetIsAromatic() for donor_idx in donor_overlap):
        weight += 12.0

    plane_offset = float(np.dot(np.asarray(metal_pos, dtype=float) - centroid, normal))
    penalty = weight * plane_offset * plane_offset
    if abs(plane_offset) > 0.28:
        penalty += 75.0 * (abs(plane_offset) - 0.28) ** 2
    return penalty


def _fragment_planar_donor_frame(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
):
    """Return a simple planar donor frame for rigid early docking."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return None

    fragment_set = set(fragment_atom_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = [
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    ]
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in fragment_set and donor_idx in planar_atoms
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return None

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return None
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return None
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return None
    normal = normal / normal_norm

    donor_pts = np.asarray([coord_map[atom_idx] for atom_idx in donor_overlap], dtype=float)
    donor_centroid = donor_pts.mean(axis=0)

    body_pts = []
    for atom_idx in planar_atoms:
        if atom_idx in donor_overlap:
            continue
        body_pts.append(np.asarray(coord_map[atom_idx], dtype=float))
    if body_pts:
        body_centroid = np.asarray(body_pts, dtype=float).mean(axis=0)
    else:
        body_centroid = centroid

    outward = donor_centroid - body_centroid
    outward = outward - float(np.dot(outward, normal)) * normal
    outward_norm = float(np.linalg.norm(outward))
    if outward_norm < 1e-12:
        outward = donor_centroid - centroid
        outward = outward - float(np.dot(outward, normal)) * normal
        outward_norm = float(np.linalg.norm(outward))
    if outward_norm < 1e-12:
        return None
    outward = outward / outward_norm

    return {
        "planar_atoms": tuple(planar_atoms),
        "donor_overlap": tuple(donor_overlap),
        "centroid": centroid,
        "donor_centroid": donor_centroid,
        "normal": normal,
        "outward": outward,
    }


def _planar_fragment_donor_approach_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    metal_pos,
    coord_map: Dict[int, object],
) -> float:
    """Penalty when a planar donor fragment presents adjacent pi atoms instead of the donor side."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    fragment_set = set(fragment_atom_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = tuple(
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    )
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in fragment_set and donor_idx in planar_atoms
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return 0.0

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return 0.0
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return 0.0
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return 0.0
    normal = normal / normal_norm

    penalty = 0.0
    metal_vec_ref = np.asarray(metal_pos, dtype=float)
    for donor_idx in donor_overlap:
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_pos = coord_map.get(donor_idx)
        if donor_pos is None:
            continue
        donor_pos = np.asarray(donor_pos, dtype=float)

        outward_terms = []
        for nbr in donor_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in fragment_set:
                continue
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            if not _is_planar_like(nbr_atom):
                continue
            nbr_pos = coord_map.get(nbr_idx)
            if nbr_pos is None:
                continue
            vec = donor_pos - np.asarray(nbr_pos, dtype=float)
            vec = vec - float(np.dot(vec, normal)) * normal
            vec_norm = float(np.linalg.norm(vec))
            if vec_norm < 1e-8:
                continue
            outward_terms.append(vec / vec_norm)

        if outward_terms:
            outward_vec = np.sum(np.asarray(outward_terms, dtype=float), axis=0)
        else:
            outward_vec = donor_pos - centroid
            outward_vec = outward_vec - float(np.dot(outward_vec, normal)) * normal
        outward_norm = float(np.linalg.norm(outward_vec))
        if outward_norm < 1e-8:
            continue
        outward_unit = outward_vec / outward_norm

        metal_vec = metal_vec_ref - donor_pos
        metal_in_plane = metal_vec - float(np.dot(metal_vec, normal)) * normal
        metal_in_plane_norm = float(np.linalg.norm(metal_in_plane))
        if metal_in_plane_norm < 1e-8:
            continue
        metal_unit = metal_in_plane / metal_in_plane_norm

        cosang = float(np.dot(outward_unit, metal_unit))
        donor_sym = donor_atom.GetSymbol()
        min_cos = 0.18
        weight = 10.0
        if donor_sym == 'N':
            min_cos = 0.34 if donor_atom.GetIsAromatic() else 0.26
            weight = 18.0 if donor_atom.GetIsAromatic() else 14.0
        elif donor_sym == 'O':
            min_cos = 0.28
            weight = 16.0
            if any(bond.GetIsConjugated() for bond in donor_atom.GetBonds()):
                min_cos += 0.06
                weight += 6.0
        elif donor_sym == 'S':
            min_cos = 0.22
            weight = 12.0

        if cosang < min_cos:
            penalty += weight * (min_cos - cosang) ** 2
        if cosang < -0.05:
            penalty += (weight + 18.0) * (abs(cosang) + 0.05) ** 2

    return penalty


def _fragment_environment_clash_penalty(
    mol,
    fragment_atom_indices: List[int],
    coord_map: Dict[int, object],
    ignore_atom_indices: Optional[set] = None,
    primary_metal_idx: Optional[int] = None,
) -> float:
    """Penalty for a candidate rigid fragment pose colliding with the fixed environment."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not coord_map:
        return 0.0

    try:
        conf = mol.GetConformer()
    except Exception:
        return 0.0

    ignore = set(ignore_atom_indices or set())
    fragment_set = set(fragment_atom_indices)
    penalty = 0.0

    for atom_idx in fragment_atom_indices:
        if atom_idx in ignore:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        sym_i = atom.GetSymbol()
        if atom.GetAtomicNum() <= 1 or sym_i in _METAL_SET:
            continue
        pos_i = coord_map.get(atom_idx)
        if pos_i is None:
            continue
        pos_i = np.asarray(pos_i, dtype=float)
        if pos_i.shape != (3,) or not np.all(np.isfinite(pos_i)):
            continue
        for other_idx in range(mol.GetNumAtoms()):
            if other_idx in fragment_set or other_idx in ignore:
                continue
            other = mol.GetAtomWithIdx(other_idx)
            sym_j = other.GetSymbol()
            if other.GetAtomicNum() <= 1 or sym_j in _METAL_SET:
                continue
            if mol.GetBondBetweenAtoms(atom_idx, other_idx) is not None:
                continue
            other_pos = conf.GetAtomPosition(other_idx)
            pos_j = np.array([other_pos.x, other_pos.y, other_pos.z], dtype=float)
            dist = float(np.linalg.norm(pos_i - pos_j))
            other_metal_bound = any(
                nbr.GetSymbol() in _METAL_SET and nbr.GetIdx() != primary_metal_idx
                for nbr in other.GetNeighbors()
            )
            min_dist = max(
                1.76,
                _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.34,
            )
            if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                min_dist = max(min_dist, 1.90)
            if other_metal_bound:
                min_dist = max(min_dist, 2.12)
            if dist < min_dist:
                weight = 18.0 if other_metal_bound else 10.0
                penalty += weight * (min_dist - dist) ** 2
            if dist < 1.55:
                weight = 120.0 if other_metal_bound else 70.0
                penalty += weight * (1.55 - dist) ** 2

    return penalty


def _fragment_bite_direction_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
    metal_pos,
) -> float:
    """Penalty when a multidentate rigid fragment presents its back side to the metal."""
    import numpy as np

    if mol is None or len(donor_indices) < 2 or not coord_map:
        return 0.0

    donor_coords = []
    body_coords = []
    fallback_coords = []
    for atom_idx in fragment_atom_indices:
        coord = coord_map.get(atom_idx)
        if coord is None:
            continue
        arr = np.asarray(coord, dtype=float)
        if arr.shape != (3,) or not np.all(np.isfinite(arr)):
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        fallback_coords.append(arr)
        if atom_idx in donor_indices:
            donor_coords.append(arr)
        else:
            body_coords.append(arr)

    if len(donor_coords) < 2:
        return 0.0
    if not body_coords:
        body_coords = fallback_coords
    if not body_coords:
        return 0.0

    donor_centroid = np.mean(np.asarray(donor_coords, dtype=float), axis=0)
    body_centroid = np.mean(np.asarray(body_coords, dtype=float), axis=0)
    pocket_vec = donor_centroid - body_centroid
    metal_vec = np.asarray(metal_pos, dtype=float) - donor_centroid
    pocket_norm = float(np.linalg.norm(pocket_vec))
    metal_norm = float(np.linalg.norm(metal_vec))
    if pocket_norm < 1e-8 or metal_norm < 1e-8:
        return 0.0

    pocket_unit = pocket_vec / pocket_norm
    metal_unit = metal_vec / metal_norm
    cosang = float(np.dot(pocket_unit, metal_unit))
    if cosang >= 0.55:
        return 0.0
    gap = 0.55 - cosang
    return 6.0 * gap * gap


def _fit_secondary_metal_geometry_from_donors(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    constrained_indices: Optional[List[int]] = None,
    bite_distance_constraints: Optional[List[Tuple[int, int, float]]] = None,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
):
    """Fit an ideal coordination geometry while optionally honoring anchor donors."""
    fits = _enumerate_secondary_metal_geometry_fits(
        donor_positions,
        donor_symbols,
        metal_symbol,
        constrained_indices=constrained_indices,
        bite_distance_constraints=bite_distance_constraints,
        donor_target_lengths=donor_target_lengths,
        donor_fit_weights=donor_fit_weights,
        max_candidates=1,
    )
    if not fits:
        return None
    metal_pos, geom_code, transformed, _score = fits[0]
    return metal_pos, geom_code, transformed


def _enumerate_secondary_metal_geometry_fits(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    constrained_indices: Optional[List[int]] = None,
    bite_distance_constraints: Optional[List[Tuple[int, int, float]]] = None,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
    max_candidates: int = 4,
) -> List[Tuple[object, str, object, float]]:
    """Return ranked geometry fits for a secondary metal donor set."""
    try:
        import itertools
        import numpy as np
    except ImportError:
        return []

    donors = np.asarray(donor_positions, dtype=float)
    if donors.ndim != 2 or donors.shape[0] != len(donor_symbols) or donors.shape[0] < 2:
        return []

    n_coord = donors.shape[0]
    target_lengths = [
        float(donor_target_lengths[i]) if donor_target_lengths and i < len(donor_target_lengths)
        else float(_get_ml_bond_length(metal_symbol, donor_symbols[i]))
        for i in range(n_coord)
    ]
    fit_weights = np.asarray(
        [
            float(donor_fit_weights[i]) if donor_fit_weights and i < len(donor_fit_weights) else 1.0
            for i in range(n_coord)
        ],
        dtype=float,
    )
    geometry_codes = _secondary_metal_geometry_codes(n_coord, metal_symbol)
    if not geometry_codes:
        return []

    fit_indices = sorted({
        int(idx) for idx in (constrained_indices or [])
        if 0 <= int(idx) < n_coord
    })
    if len(fit_indices) < 2:
        fit_indices = list(range(n_coord))

    nonpreferred_cn4_penalty = 0.0
    if n_coord == 4:
        preferred_cn4 = _PREFERRED_CN4_GEOMETRY.get(metal_symbol)
        if preferred_cn4:
            if metal_symbol in {'Pt', 'Pd', 'Rh', 'Ir', 'Au'}:
                nonpreferred_cn4_penalty = 0.24
            elif metal_symbol == 'Ni':
                nonpreferred_cn4_penalty = 0.12
            else:
                nonpreferred_cn4_penalty = 0.16
        else:
            preferred_cn4 = None
    else:
        preferred_cn4 = None

    ranked_candidates: List[Tuple[float, float, float, float, object, str, object]] = []
    for geom_code in geometry_codes:
        vectors = _TOPO_GEOMETRY_VECTORS.get(geom_code)
        if not vectors or len(vectors) != n_coord:
            continue
        normed_vectors = []
        for vec in vectors:
            arr = np.asarray(vec, dtype=float)
            norm = float(np.linalg.norm(arr))
            if norm < 1e-12:
                break
            normed_vectors.append(arr / norm)
        if len(normed_vectors) != n_coord:
            continue

        if n_coord <= 6:
            permutations = itertools.permutations(range(n_coord))
        else:
            identity = tuple(range(n_coord))
            permutations = [identity, tuple(reversed(identity))]

        for perm in permutations:
            model = np.asarray(
                [
                    normed_vectors[perm[i]] * target_lengths[i]
                    for i in range(n_coord)
                ],
                dtype=float,
            )
            fitted = _fit_secondary_metal_model_to_targets(model, donors, fit_indices)
            if fitted is None:
                continue
            transformed, metal_pos = fitted
            fit_resid = np.sum((transformed[fit_indices] - donors[fit_indices]) ** 2, axis=1)
            fit_w = fit_weights[fit_indices]
            rmsd_fit = float(np.sqrt(np.sum(fit_w * fit_resid) / max(np.sum(fit_w), 1e-12)))
            all_resid = np.sum((transformed - donors) ** 2, axis=1)
            rmsd_all = float(np.sqrt(np.sum(fit_weights * all_resid) / max(np.sum(fit_weights), 1e-12)))
            bite_rmsd = 0.0
            if bite_distance_constraints:
                bite_errors = []
                for idx_i, idx_j, expected_dist in bite_distance_constraints:
                    if idx_i < 0 or idx_j < 0 or idx_i >= n_coord or idx_j >= n_coord:
                        continue
                    actual_dist = float(np.linalg.norm(transformed[idx_i] - transformed[idx_j]))
                    bite_errors.append((actual_dist - expected_dist) ** 2)
                if bite_errors:
                    bite_rmsd = float(np.sqrt(np.mean(bite_errors)))
                score = rmsd_fit + 0.10 * rmsd_all + 0.55 * bite_rmsd
                if preferred_cn4 is not None and geom_code != preferred_cn4:
                    score += nonpreferred_cn4_penalty
                ranked_candidates.append(
                    (
                        float(score),
                        float(rmsd_fit),
                        float(rmsd_all),
                        float(bite_rmsd),
                        np.asarray(metal_pos, dtype=float),
                        geom_code,
                        np.asarray(transformed, dtype=float),
                    )
                )

    if not ranked_candidates:
        return []

    ranked_candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    selected: List[Tuple[object, str, object, float]] = []

    def _is_distinct_geometry_candidate(
        transformed,
        geom_code: str,
        existing: List[Tuple[object, str, object, float]],
    ) -> bool:
        transformed_arr = np.asarray(transformed, dtype=float)
        for _pos, prev_geom, prev_transformed, _score in existing:
            if geom_code != prev_geom:
                continue
            prev_arr = np.asarray(prev_transformed, dtype=float)
            if prev_arr.shape != transformed_arr.shape:
                continue
            deltas = np.linalg.norm(prev_arr - transformed_arr, axis=1)
            if float(np.max(deltas)) < 0.22:
                return False
        return True

    first = ranked_candidates[0]
    selected.append((first[4], first[5], first[6], first[0]))
    seen_geometries = {first[5]}

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[5]
        if geom_code in seen_geometries:
            continue
        if not _is_distinct_geometry_candidate(candidate[6], geom_code, selected):
            continue
        selected.append((candidate[4], geom_code, candidate[6], candidate[0]))
        seen_geometries.add(geom_code)

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[5]
        if not _is_distinct_geometry_candidate(candidate[6], geom_code, selected):
            continue
        selected.append((candidate[4], geom_code, candidate[6], candidate[0]))

    return selected[:max(1, int(max_candidates))]


def _fit_secondary_metal_position_from_donors(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
):
    """Fit a secondary metal center to already-placed donor atoms."""
    fits = _enumerate_secondary_metal_position_fits(
        donor_positions,
        donor_symbols,
        metal_symbol,
        donor_target_lengths=donor_target_lengths,
        donor_fit_weights=donor_fit_weights,
        max_candidates=1,
    )
    if not fits:
        return None
    metal_pos, geom_code, _score = fits[0]
    return metal_pos, geom_code


def _enumerate_secondary_metal_position_fits(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
    max_candidates: int = 4,
) -> List[Tuple[object, str, float]]:
    """Return ranked metal-position fits for already placed donors."""
    try:
        import itertools
        import numpy as np
    except ImportError:
        return []

    donors = np.asarray(donor_positions, dtype=float)
    if donors.ndim != 2 or donors.shape[0] != len(donor_symbols) or donors.shape[0] < 2:
        return []

    n_coord = donors.shape[0]
    target_lengths = [
        float(donor_target_lengths[i]) if donor_target_lengths and i < len(donor_target_lengths)
        else float(_get_ml_bond_length(metal_symbol, donor_symbols[i]))
        for i in range(n_coord)
    ]
    fit_weights = np.asarray(
        [
            float(donor_fit_weights[i]) if donor_fit_weights and i < len(donor_fit_weights) else 1.0
            for i in range(n_coord)
        ],
        dtype=float,
    )
    geometry_codes = _secondary_metal_geometry_codes(n_coord, metal_symbol)
    if not geometry_codes:
        return []

    ranked_candidates: List[Tuple[float, float, object, str]] = []
    for geom_code in geometry_codes:
        vectors = _TOPO_GEOMETRY_VECTORS.get(geom_code)
        if not vectors or len(vectors) != n_coord:
            continue
        normed_vectors = []
        for vec in vectors:
            arr = np.asarray(vec, dtype=float)
            norm = float(np.linalg.norm(arr))
            if norm < 1e-12:
                break
            normed_vectors.append(arr / norm)
        if len(normed_vectors) != n_coord:
            continue

        if n_coord <= 6:
            permutations = itertools.permutations(range(n_coord))
        else:
            identity = tuple(range(n_coord))
            permutations = [identity, tuple(reversed(identity))]

        for perm in permutations:
            model = np.asarray(
                [
                    normed_vectors[perm[i]] * target_lengths[i]
                    for i in range(n_coord)
                ],
                dtype=float,
            )
            model_centroid = model.mean(axis=0)
            donor_centroid = donors.mean(axis=0)
            covariance = (model - model_centroid).T @ (donors - donor_centroid)
            try:
                u, _s, vt = np.linalg.svd(covariance)
            except Exception:
                continue
            rot = vt.T @ u.T
            if float(np.linalg.det(rot)) < 0.0:
                vt[-1, :] *= -1.0
                rot = vt.T @ u.T
            fitted = (model - model_centroid) @ rot.T + donor_centroid
            resid = np.sum((fitted - donors) ** 2, axis=1)
            rmsd = float(np.sqrt(np.sum(fit_weights * resid) / max(np.sum(fit_weights), 1e-12)))
            metal_pos = donor_centroid - model_centroid @ rot.T
            score = rmsd
            ranked_candidates.append(
                (float(score), float(rmsd), np.asarray(metal_pos, dtype=float), geom_code)
            )

    if not ranked_candidates:
        return []

    ranked_candidates.sort(key=lambda item: (item[0], item[1]))

    selected: List[Tuple[object, str, float]] = []

    def _is_distinct_position_candidate(
        metal_pos,
        geom_code: str,
        existing: List[Tuple[object, str, float]],
    ) -> bool:
        metal_arr = np.asarray(metal_pos, dtype=float)
        for prev_pos, prev_geom, _score in existing:
            if geom_code != prev_geom:
                continue
            prev_arr = np.asarray(prev_pos, dtype=float)
            if float(np.linalg.norm(prev_arr - metal_arr)) < 0.20:
                return False
        return True

    first = ranked_candidates[0]
    selected.append((first[2], first[3], first[0]))
    seen_geometries = {first[3]}

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[3]
        if geom_code in seen_geometries:
            continue
        if not _is_distinct_position_candidate(candidate[2], geom_code, selected):
            continue
        selected.append((candidate[2], geom_code, candidate[0]))
        seen_geometries.add(geom_code)

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[3]
        if not _is_distinct_position_candidate(candidate[2], geom_code, selected):
            continue
        selected.append((candidate[2], geom_code, candidate[0]))

    return selected[:max(1, int(max_candidates))]


def _secondary_module_fragment_is_movable(
    fragment: _HybridHaptoFragment,
    metal_idx: int,
) -> bool:
    """Return whether a fragment may be moved as part of one secondary-metal module."""
    return (
        metal_idx in fragment.metal_neighbor_indices
        and len(fragment.metal_neighbor_indices) == 1
        and not fragment.bridging_donor_indices
        and not fragment.hapto_group_ids
    )


def _find_scaffold_secondary_branch(
    mol,
    fragment: _HybridHaptoFragment,
    hapto_groups: List[Tuple[int, List[int]]],
    donor_idx: int,
) -> Optional[Tuple[int, set]]:
    """Return a rotatable donor branch off a fixed hapto/bridge scaffold, if present."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None:
        return None
    if donor_idx not in fragment.atom_indices or donor_idx in fragment.bridging_donor_indices:
        return None

    fragment_atoms = set(fragment.atom_indices)
    core_seed_atoms: set = set(fragment.bridging_donor_indices)
    for group_id in fragment.hapto_group_ids:
        if 0 <= int(group_id) < len(hapto_groups):
            _metal_idx, group_atoms = hapto_groups[int(group_id)]
            core_seed_atoms.update(atom_idx for atom_idx in group_atoms if atom_idx in fragment_atoms)
    for other_donor in fragment.donor_atom_indices:
        if other_donor != donor_idx and other_donor in fragment_atoms:
            core_seed_atoms.add(other_donor)
    if donor_idx in core_seed_atoms:
        return None

    movable_pool = fragment_atoms - core_seed_atoms
    if donor_idx not in movable_pool:
        return None

    branch_atoms: set = set()
    queue = [donor_idx]
    while queue:
        atom_idx = queue.pop()
        if atom_idx in branch_atoms or atom_idx not in movable_pool:
            continue
        branch_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in movable_pool and nbr_idx not in branch_atoms:
                queue.append(nbr_idx)

    if donor_idx not in branch_atoms or len(branch_atoms) < 2:
        return None

    attachment_atoms: set = set()
    for atom_idx in branch_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_seed_atoms:
                attachment_atoms.add(nbr_idx)
    interface_atoms = {
        atom_idx
        for atom_idx in branch_atoms
        if any(nbr.GetIdx() in core_seed_atoms for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors())
    }
    if len(interface_atoms) == 1:
        pivot_idx = next(iter(interface_atoms))
        branch_atoms = set(branch_atoms)
        branch_atoms.discard(pivot_idx)
    else:
        if len(attachment_atoms) != 1:
            return None
        pivot_idx = next(iter(attachment_atoms))

    if donor_idx == pivot_idx or donor_idx not in branch_atoms:
        return None
    return pivot_idx, branch_atoms


def _reorient_scaffold_secondary_branch(
    mol,
    fragment: _HybridHaptoFragment,
    hapto_groups: List[Tuple[int, List[int]]],
    metal_idx: int,
    donor_idx: int,
    donor_target,
) -> bool:
    """Rotate a branch off a fixed hapto/bridge core toward a secondary-metal target."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    branch_info = _find_scaffold_secondary_branch(
        mol,
        fragment,
        hapto_groups,
        donor_idx,
    )
    if branch_info is None:
        return False
    pivot_idx, branch_atoms = branch_info

    pivot_pos = np.array(conf.GetAtomPosition(pivot_idx), dtype=float)
    current_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in branch_atoms
    }
    donor_target_vec = np.asarray(donor_target, dtype=float)
    if donor_target_vec.shape != (3,) or not np.all(np.isfinite(donor_target_vec)):
        return False

    current_donor_vec = current_coords[donor_idx] - pivot_pos
    target_donor_vec = donor_target_vec - pivot_pos
    if float(np.linalg.norm(current_donor_vec)) < 1e-12 or float(np.linalg.norm(target_donor_vec)) < 1e-12:
        return False

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    base_rot = _rotation_matrix_from_vectors(current_donor_vec, target_donor_vec)
    rotated_coords = {
        atom_idx: (coord - pivot_pos) @ base_rot.T + pivot_pos
        for atom_idx, coord in current_coords.items()
    }
    branch_frame = _fragment_planar_donor_frame(
        mol,
        sorted(branch_atoms),
        [donor_idx],
        rotated_coords,
    )
    if branch_frame is not None:
        outward_src = np.asarray(branch_frame["outward"], dtype=float)
        metal_dir = metal_pos - donor_target_vec
        metal_dir = metal_dir - float(np.dot(metal_dir, branch_frame["normal"])) * np.asarray(branch_frame["normal"], dtype=float)
        if float(np.linalg.norm(metal_dir)) > 1e-8:
            outward_tgt = metal_dir / float(np.linalg.norm(metal_dir))
            donor_center = rotated_coords[donor_idx]
            branch_rot = _rotation_matrix_from_vectors(outward_src, outward_tgt)
            trial_coords = {
                atom_idx: (coord - donor_center) @ branch_rot.T + donor_center
                for atom_idx, coord in rotated_coords.items()
            }
            if donor_idx in trial_coords:
                donor_shift = donor_target_vec - trial_coords[donor_idx]
                trial_coords = {
                    atom_idx: coord + donor_shift
                    for atom_idx, coord in trial_coords.items()
                }
            rotated_coords = trial_coords

    def _score(coord_map: Dict[int, object]) -> Tuple[float, float]:
        donor_err = float(np.linalg.norm(coord_map[donor_idx] - donor_target_vec))
        penalty = _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(branch_atoms),
            [donor_idx],
            coord_map,
        )
        penalty += _fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(branch_atoms),
            [donor_idx],
            coord_map,
        )
        penalty += _planar_fragment_metal_coplanarity_penalty(
            mol,
            sorted(branch_atoms),
            [donor_idx],
            metal_pos,
            coord_map,
        )
        penalty += _planar_fragment_donor_approach_penalty(
            mol,
            sorted(branch_atoms),
            [donor_idx],
            metal_pos,
            coord_map,
        )
        return donor_err, 30.0 * donor_err * donor_err + penalty

    current_err, current_score = _score(current_coords)
    best_coords = rotated_coords
    best_err, best_score = _score(rotated_coords)

    axis = rotated_coords[donor_idx] - pivot_pos
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm > 1e-12:
        axis /= axis_norm
        step = 5 if branch_frame is not None else 15
        for angle_deg in range(-180, 181, step):
            if angle_deg == 0:
                continue
            rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
            trial_coords = {
                atom_idx: (coord - pivot_pos) @ rot.T + pivot_pos
                for atom_idx, coord in rotated_coords.items()
            }
            trial_err, trial_score = _score(trial_coords)
            if (
                trial_score + 1e-9 < best_score
                or (
                    abs(trial_score - best_score) < 1e-9
                    and trial_err + 1e-6 < best_err
                )
            ):
                best_coords = trial_coords
                best_err = trial_err
                best_score = trial_score

    if best_score + 1e-9 >= current_score and best_err + 1e-6 >= current_err:
        return False

    for atom_idx, coord in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
    return True


def _optimize_secondary_fragment_pose(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
) -> bool:
    """Improve a rigid donor fragment around one anchored donor without distorting it."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not donor_target_map:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    donor_indices = [
        donor_idx for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map
    ]
    if len(donor_indices) < 2:
        return False

    base_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in fragment.atom_indices
    }
    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    target_lengths = {
        donor_idx: float(_get_ml_bond_length(metal_sym, mol.GetAtomWithIdx(donor_idx).GetSymbol()))
        for donor_idx in donor_indices
    }

    donor_errors = {
        donor_idx: abs(float(np.linalg.norm(base_coords[donor_idx] - metal_pos)) - target_lengths[donor_idx])
        for donor_idx in donor_indices
    }
    pivot_idx = min(donor_indices, key=lambda donor_idx: donor_errors[donor_idx])
    pivot_pos = np.array(base_coords[pivot_idx], dtype=float)
    current_errors = donor_errors.copy()

    def _score(coord_map: Dict[int, object]) -> float:
        score = 0.0
        for donor_idx in donor_indices:
            donor_pos = coord_map[donor_idx]
            target_len = target_lengths[donor_idx]
            err = float(np.linalg.norm(donor_pos - metal_pos)) - target_len
            weight = 0.3 if donor_idx == pivot_idx else 3.0
            score += weight * err * err
        score += _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            fragment.atom_indices,
            donor_indices,
            coord_map,
        )
        score += _fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            fragment.atom_indices,
            donor_indices,
            coord_map,
        )
        score += _fragment_environment_clash_penalty(
            mol,
            fragment.atom_indices,
            coord_map,
            ignore_atom_indices={metal_idx},
            primary_metal_idx=metal_idx,
        )
        score += _planar_fragment_metal_coplanarity_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            metal_pos,
            coord_map,
        )
        score += _planar_fragment_donor_approach_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            metal_pos,
            coord_map,
        )
        score += _fragment_bite_direction_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            coord_map,
            metal_pos,
        )
        return score

    best_coords = {idx: vec.copy() for idx, vec in base_coords.items()}
    best_score = _score(best_coords)
    best_errors = current_errors.copy()
    donor_axis = metal_pos - pivot_pos
    donor_axis_norm = float(np.linalg.norm(donor_axis))
    if donor_axis_norm < 1e-12:
        return False
    donor_axis /= donor_axis_norm

    rng = np.random.default_rng(
        int((metal_idx + 1) * 1000 + min(fragment.atom_indices) + len(fragment.atom_indices))
    )
    axes = [donor_axis]
    for donor_idx in donor_indices:
        if donor_idx == pivot_idx:
            continue
        pair_axis = np.asarray(base_coords[donor_idx] - pivot_pos, dtype=float)
        pair_norm = float(np.linalg.norm(pair_axis))
        if pair_norm < 1e-12:
            continue
        axes.insert(0, pair_axis / pair_norm)
        break
    for _ in range(192):
        axis = rng.normal(size=3)
        norm = float(np.linalg.norm(axis))
        if norm < 1e-12:
            continue
        axes.append(axis / norm)

    improved = False
    for axis in axes:
        for angle_rad in (
            math.radians(v) for v in (-180, -150, -120, -90, -75, -60, -45, -30, -20, -10, 10, 20, 30, 45, 60, 75, 90, 120, 150, 180)
        ):
            rot = _axis_angle_rotation_matrix(axis, angle_rad)
            trial_coords = {}
            for atom_idx, vec in base_coords.items():
                if atom_idx == pivot_idx:
                    trial_coords[atom_idx] = vec.copy()
                    continue
                trial_coords[atom_idx] = (vec - pivot_pos) @ rot.T + pivot_pos
            trial_score = _score(trial_coords)
            trial_errors = {
                donor_idx: abs(float(np.linalg.norm(trial_coords[donor_idx] - metal_pos)) - target_lengths[donor_idx])
                for donor_idx in donor_indices
            }
            trial_max_nonpivot = max(
                (err for donor_idx, err in trial_errors.items() if donor_idx != pivot_idx),
                default=0.0,
            )
            best_max_nonpivot = max(
                (err for donor_idx, err in best_errors.items() if donor_idx != pivot_idx),
                default=0.0,
            )
            if (
                trial_score + 1e-9 < best_score
                and trial_max_nonpivot <= best_max_nonpivot + 1e-4
            ):
                best_score = trial_score
                best_coords = trial_coords
                best_errors = trial_errors
                improved = True

    if not improved:
        return False

    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(
            atom_idx,
            Point3D(float(vec[0]), float(vec[1]), float(vec[2])),
        )
    return True


def _secondary_fragment_prefers_rigid_pose(
    mol,
    fragment: _HybridHaptoFragment,
    targeted_donor_indices: List[int],
) -> bool:
    """Return whether a CN=4 secondary-metal fragment should keep its rigid aligned pose."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not targeted_donor_indices:
        return False

    targeted_donor_set = set(targeted_donor_indices)
    if len(targeted_donor_set) >= 2:
        return True

    planar_like_atoms = 0
    for atom_idx in fragment.atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        if atom.GetIsAromatic():
            planar_like_atoms += 1
            continue
        try:
            hyb = atom.GetHybridization()
        except Exception:
            continue
        if hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }:
            planar_like_atoms += 1

    return planar_like_atoms >= 4 and bool(targeted_donor_set & set(fragment.donor_atom_indices))


_HAPTO_SECONDARY_CN4_GEOMETRIES = {'SQ', 'TH', 'TET'}


_HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES = {'SQ', 'TET'}


def _relieve_secondary_oo_chelate_contacts(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
) -> bool:
    """Rotate a rigid O,O chelate around the O-O axis to relieve non-donor contacts."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not donor_target_map:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    donor_indices = [
        donor_idx for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map and mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O'
    ]
    if len(donor_indices) < 2:
        return False

    axis_start = np.array(conf.GetAtomPosition(donor_indices[0]), dtype=float)
    axis_end = np.array(conf.GetAtomPosition(donor_indices[1]), dtype=float)
    axis = axis_end - axis_start
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        return False
    axis = axis / axis_norm

    rotatable_atoms = [
        atom_idx for atom_idx in fragment.atom_indices
        if atom_idx not in donor_indices
        and mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
    ]
    if not rotatable_atoms:
        return False

    base_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in fragment.atom_indices
    }
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)

    def _score(coord_map: Dict[int, object]) -> float:
        return (
            _secondary_non_donor_contact_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_donor_selectivity_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _planar_fragment_metal_coplanarity_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _fragment_environment_clash_penalty(
                mol,
                fragment.atom_indices,
                coord_map,
                ignore_atom_indices={metal_idx},
                primary_metal_idx=metal_idx,
            )
            + _fragment_bite_direction_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                coord_map,
                metal_pos,
            )
            + _planar_fragment_metal_coplanarity_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
        )

    best_coords = {atom_idx: vec.copy() for atom_idx, vec in base_coords.items()}
    best_score = _score(best_coords)
    improved = False
    pivot = 0.5 * (axis_start + axis_end)

    for angle_deg in (-30, -25, -20, -15, -10, -5, 5, 10, 15, 20, 25, 30):
        rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
        trial_coords = {atom_idx: vec.copy() for atom_idx, vec in base_coords.items()}
        for atom_idx in rotatable_atoms:
            trial_coords[atom_idx] = (base_coords[atom_idx] - pivot) @ rot.T + pivot
        trial_score = _score(trial_coords)
        if trial_score + 1e-9 < best_score:
            best_score = trial_score
            best_coords = trial_coords
            improved = True

    if not improved:
        return False

    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))
    return True


def _restore_secondary_rigid_fragment_geometry(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
    reference_coords: Dict[int, object],
) -> bool:
    """Restore a rigid secondary fragment against its saved pre-optimizer geometry."""
    if (
        not RDKIT_AVAILABLE
        or mol is None
        or fragment is None
        or not donor_target_map
        or not reference_coords
    ):
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    donor_indices = [
        donor_idx
        for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map and donor_idx in reference_coords
    ]
    if len(donor_indices) < 2:
        return False

    atom_indices = [
        atom_idx for atom_idx in fragment.atom_indices
        if atom_idx in reference_coords
    ]
    if len(atom_indices) < 3:
        return False

    src = np.asarray([reference_coords[donor_idx] for donor_idx in donor_indices], dtype=float)
    tgt = np.asarray(
        [np.array(conf.GetAtomPosition(donor_idx), dtype=float) for donor_idx in donor_indices],
        dtype=float,
    )
    all_ref = np.asarray([reference_coords[atom_idx] for atom_idx in atom_indices], dtype=float)

    if len(donor_indices) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        src_mid = 0.5 * (src[0] + src[1])
        tgt_mid = 0.5 * (tgt[0] + tgt[1])
        aligned = (all_ref - src_mid) @ rot.T + tgt_mid
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_ref - src_centroid) @ rot.T + tgt_centroid

    donor_set = set(donor_indices)
    current_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in atom_indices
    }
    rigid_coords = {
        atom_idx: aligned[idx].copy()
        for idx, atom_idx in enumerate(atom_indices)
    }
    for row_idx, donor_idx in enumerate(donor_indices):
        rigid_coords[donor_idx] = tgt[row_idx].copy()

    heavy_atoms = [
        atom_idx for atom_idx in atom_indices
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
    ]
    if len(heavy_atoms) < 3:
        return False

    ref_pair_targets: List[Tuple[int, int, float]] = []
    for i in range(len(heavy_atoms)):
        for j in range(i + 1, len(heavy_atoms)):
            atom_i = heavy_atoms[i]
            atom_j = heavy_atoms[j]
            ref_pair_targets.append(
                (
                    atom_i,
                    atom_j,
                    float(
                        np.linalg.norm(
                            np.asarray(reference_coords[atom_i], dtype=float)
                            - np.asarray(reference_coords[atom_j], dtype=float)
                        )
                    ),
                )
            )

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    planar_atoms: List[int] = []
    for atom_idx in heavy_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            planar_atoms.append(atom_idx)
            continue
        try:
            hyb = atom.GetHybridization()
        except Exception:
            continue
        if hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }:
            planar_atoms.append(atom_idx)

    def _external_score(coord_map: Dict[int, object]) -> float:
        return (
            _secondary_non_donor_contact_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_donor_selectivity_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_environment_clash_penalty(
                mol,
                fragment.atom_indices,
                coord_map,
                ignore_atom_indices={metal_idx},
                primary_metal_idx=metal_idx,
            )
            + _fragment_bite_direction_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                coord_map,
                metal_pos,
            )
        )

    def _shape_penalty(coord_map: Dict[int, object]) -> float:
        penalty = 0.0
        for atom_i, atom_j, target_len in ref_pair_targets:
            vec_i = np.asarray(coord_map.get(atom_i), dtype=float)
            vec_j = np.asarray(coord_map.get(atom_j), dtype=float)
            if vec_i.shape != (3,) or vec_j.shape != (3,):
                continue
            penalty += (float(np.linalg.norm(vec_i - vec_j)) - target_len) ** 2
        return penalty

    def _planarity_penalty(coord_map: Dict[int, object]) -> float:
        if len(planar_atoms) < 4:
            return 0.0
        pts = np.asarray([coord_map[atom_idx] for atom_idx in planar_atoms], dtype=float)
        if pts.shape[0] < 4:
            return 0.0
        centroid = pts.mean(axis=0)
        try:
            _u, _s, vt = np.linalg.svd(pts - centroid)
        except Exception:
            return 0.0
        normal = vt[-1]
        normal_norm = float(np.linalg.norm(normal))
        if normal_norm < 1e-12:
            return 0.0
        normal = normal / normal_norm
        rms = float(np.sqrt(np.mean(((pts - centroid) @ normal) ** 2)))
        return rms * rms

    def _total_score(coord_map: Dict[int, object]) -> float:
        return (
            _external_score(coord_map)
            + 6.0 * _shape_penalty(coord_map)
            + 28.0 * _planarity_penalty(coord_map)
        )

    current_score = _total_score(current_coords)
    best_coords = {atom_idx: vec.copy() for atom_idx, vec in rigid_coords.items()}
    best_score = _total_score(best_coords)

    if len(donor_indices) == 2:
        axis = tgt[1] - tgt[0]
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis = axis / axis_norm
            pivot = 0.5 * (tgt[0] + tgt[1])
            rotatable_atoms = [atom_idx for atom_idx in atom_indices if atom_idx not in donor_set]
            for angle_deg in range(-180, 181, 15):
                if angle_deg == 0:
                    continue
                rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
                trial_coords = {
                    atom_idx: vec.copy() for atom_idx, vec in rigid_coords.items()
                }
                for atom_idx in rotatable_atoms:
                    trial_coords[atom_idx] = (rigid_coords[atom_idx] - pivot) @ rot.T + pivot
                trial_score = _total_score(trial_coords)
                if trial_score + 1e-9 < best_score:
                    best_score = trial_score
                    best_coords = trial_coords

    if best_score + 1e-9 >= current_score:
        return False

    for row_idx, donor_idx in enumerate(donor_indices):
        best_coords[donor_idx] = tgt[row_idx].copy()
    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))
    return True


def _optimize_secondary_metal_module_local(
    mol,
    fragments: List[_HybridHaptoFragment],
    hapto_groups: List[Tuple[int, List[int]]],
    metal_idx: int,
    donor_indices: List[int],
    donor_to_fragment: Dict[int, _HybridHaptoFragment],
    geometry_code: Optional[str] = None,
    geometry_metal_pos=None,
    geometry_target_map: Optional[Dict[int, object]] = None,
) -> set:
    """Locally relax one hapto-coupled secondary-metal module with the hapto core fixed."""
    if not RDKIT_AVAILABLE or mol is None or not donor_indices:
        return set()
    try:
        import numpy as np
    except ImportError:
        return set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set()

    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    donor_set = set(donor_indices)
    module_atoms: set = {metal_idx} | donor_set
    movable_atoms: set = {metal_idx}
    donor_pair_targets: List[Tuple[int, int, float]] = []
    rigid_pair_targets: List[Tuple[int, int, float, float]] = []
    seen_donor_pairs: set = set()
    rigid_units: List[Tuple[int, ...]] = []
    seen_units: set = set()
    metal_repulsion_targets: List[Tuple[int, float]] = []
    atom_to_fragment_idx: Dict[int, int] = {}

    hapto_metals = {idx for idx, _grp in hapto_groups}
    hapto_atoms = {atom_idx for _idx, grp in hapto_groups for atom_idx in grp}
    for frag_idx, fragment in enumerate(fragments):
        for atom_idx in fragment.atom_indices:
            atom_to_fragment_idx[atom_idx] = frag_idx
    for atom in mol.GetAtoms():
        other_idx = atom.GetIdx()
        if other_idx == metal_idx or atom.GetSymbol() not in _METAL_SET:
            continue
        min_dist = max(
            3.05,
            _COVALENT_RADII.get(metal_sym, 1.35) + _COVALENT_RADII.get(atom.GetSymbol(), 1.35) + 0.45,
        )
        metal_repulsion_targets.append((other_idx, float(min_dist)))

    for donor_idx in donor_indices:
        fragment = donor_to_fragment.get(donor_idx)
        if fragment is None:
            continue
        module_atoms.update(fragment.atom_indices)
        fragment_donors = [
            idx for idx in fragment.donor_atom_indices
            if idx in donor_set
        ]
        for i in range(len(fragment_donors)):
            for j in range(i + 1, len(fragment_donors)):
                pair = tuple(sorted((fragment_donors[i], fragment_donors[j])))
                if pair in seen_donor_pairs:
                    continue
                pos_i = np.array(conf.GetAtomPosition(pair[0]), dtype=float)
                pos_j = np.array(conf.GetAtomPosition(pair[1]), dtype=float)
                donor_pair_targets.append((pair[0], pair[1], float(np.linalg.norm(pos_i - pos_j))))
                seen_donor_pairs.add(pair)

        if _secondary_module_fragment_is_movable(fragment, metal_idx):
            movable_atoms.update(fragment.atom_indices)
            unit = tuple(
                sorted(
                    atom_idx for atom_idx in fragment.atom_indices
                    if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                    and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                )
            )
            if len(unit) >= 3 and unit not in seen_units:
                rigid_units.append(unit)
                seen_units.add(unit)
            continue
        if fragment.use_scaffold_only:
            branch_info = _find_scaffold_secondary_branch(
                mol,
                fragment,
                hapto_groups,
                donor_idx,
            )
            if branch_info is not None:
                pivot_idx, branch_atoms = branch_info
                module_atoms.add(pivot_idx)
                module_atoms.update(branch_atoms)
                movable_atoms.update(branch_atoms)
                unit = tuple(
                    sorted(
                        atom_idx for atom_idx in ({pivot_idx} | set(branch_atoms))
                        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                    )
                )
                if len(unit) >= 3 and unit not in seen_units:
                    rigid_units.append(unit)
                    seen_units.add(unit)

    movable_atoms -= hapto_atoms
    movable_atoms -= (hapto_metals - {metal_idx})
    movable_atoms.add(metal_idx)
    movable_atoms &= module_atoms
    if len(movable_atoms) <= 1:
        return set()

    bonded_pairs: set = set()
    bond_targets: List[Tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx not in module_atoms or end_idx not in module_atoms:
            continue
        bonded_pairs.add((min(begin_idx, end_idx), max(begin_idx, end_idx)))
        if begin_idx in movable_atoms or end_idx in movable_atoms:
            bond_targets.append((begin_idx, end_idx, _hybrid_bond_target_length(bond)))

    def _donor_target_profile(donor_idx: int) -> Tuple[float, float]:
        atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = atom.GetSymbol()
        target_len = float(_secondary_donor_target_length(mol, metal_sym, donor_idx))
        weight = 2.0 * float(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
        formal_charge = int(atom.GetFormalCharge())
        aromatic = bool(atom.GetIsAromatic())
        try:
            hyb = atom.GetHybridization()
        except Exception:
            hyb = None
        planar_like = aromatic or hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

        if donor_sym == 'N':
            if planar_like:
                weight += 0.35
            if formal_charge > 0:
                target_len += 0.05
                weight -= 0.45
            elif formal_charge < 0:
                target_len -= 0.04
                weight += 0.20
        elif donor_sym == 'O':
            if any(b.GetBondTypeAsDouble() >= 1.5 for b in atom.GetBonds()):
                weight += 0.20
            if any(b.GetIsConjugated() for b in atom.GetBonds()):
                weight += 0.15
            if formal_charge < 0:
                target_len -= 0.08
                weight += 0.70
            elif formal_charge > 0:
                target_len += 0.04
                weight -= 0.20
        elif donor_sym == 'C':
            if formal_charge < 0:
                target_len -= 0.06
                weight += 0.40
            elif formal_charge > 0:
                target_len += 0.05
                weight -= 0.30

        return target_len, max(weight, 0.85)

    ml_targets = [
        (
            donor_idx,
            *_donor_target_profile(donor_idx),
        )
        for donor_idx in donor_indices
    ]
    cn4_geometry_pair_targets: List[Tuple[int, int, float, float]] = []
    cn4_square_planar_normal = None
    if (
        geometry_code in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES
        and len(donor_indices) == 4
        and geometry_metal_pos is not None
        and geometry_target_map
    ):
        pair_weight = 0.55 if geometry_code == 'SQ' else 0.42
        for i in range(len(donor_indices)):
            pos_i = geometry_target_map.get(donor_indices[i])
            if pos_i is None:
                cn4_geometry_pair_targets = []
                break
            pos_i = np.asarray(pos_i, dtype=float)
            if pos_i.shape != (3,) or not np.all(np.isfinite(pos_i)):
                cn4_geometry_pair_targets = []
                break
            for j in range(i + 1, len(donor_indices)):
                pos_j = geometry_target_map.get(donor_indices[j])
                if pos_j is None:
                    cn4_geometry_pair_targets = []
                    break
                pos_j = np.asarray(pos_j, dtype=float)
                if pos_j.shape != (3,) or not np.all(np.isfinite(pos_j)):
                    cn4_geometry_pair_targets = []
                    break
                target_dist = float(np.linalg.norm(pos_i - pos_j))
                cn4_geometry_pair_targets.append(
                    (donor_indices[i], donor_indices[j], target_dist, pair_weight)
                )
            if not cn4_geometry_pair_targets and i > 0:
                break
        if geometry_code == 'SQ' and len(cn4_geometry_pair_targets) >= 2:
            target_vecs = []
            for donor_idx in donor_indices:
                target_pos = geometry_target_map.get(donor_idx)
                if target_pos is None:
                    target_vecs = []
                    break
                target_arr = np.asarray(target_pos, dtype=float) - np.asarray(geometry_metal_pos, dtype=float)
                if target_arr.shape != (3,) or not np.all(np.isfinite(target_arr)):
                    target_vecs = []
                    break
                target_vecs.append(target_arr)
            best_norm = 0.0
            for i in range(len(target_vecs)):
                for j in range(i + 1, len(target_vecs)):
                    normal = np.cross(target_vecs[i], target_vecs[j])
                    norm = float(np.linalg.norm(normal))
                    if norm > best_norm:
                        best_norm = norm
                        cn4_square_planar_normal = normal / norm

    initial_positions = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in movable_atoms
    }
    initial_positions.update(
        {
            atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
            for unit in rigid_units
            for atom_idx in unit
            if atom_idx not in initial_positions
        }
    )

    for unit in rigid_units:
        unit_atoms = list(unit)
        unit_planar = False
        planar_like_count = 0
        for atom_idx in unit_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            if atom.GetIsAromatic():
                planar_like_count += 1
                continue
            try:
                hyb = atom.GetHybridization()
            except Exception:
                continue
            if hyb in {
                Chem.rdchem.HybridizationType.SP,
                Chem.rdchem.HybridizationType.SP2,
            }:
                planar_like_count += 1
        unit_planar = planar_like_count >= 4

        if len(unit_atoms) <= 6 or (unit_planar and len(unit_atoms) <= 14):
            pair_iter = [
                (unit_atoms[i], unit_atoms[j])
                for i in range(len(unit_atoms))
                for j in range(i + 1, len(unit_atoms))
            ]
        else:
            pair_iter = []
            unit_set = set(unit_atoms)
            for atom_idx in unit_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in unit_set and nbr_idx > atom_idx:
                        pair_iter.append((atom_idx, nbr_idx))
                    if nbr_idx not in unit_set:
                        continue
                    for nbr2 in nbr.GetNeighbors():
                        nbr2_idx = nbr2.GetIdx()
                        if nbr2_idx in unit_set and nbr2_idx > atom_idx and nbr2_idx != atom_idx:
                            pair = (atom_idx, nbr2_idx)
                            if pair not in pair_iter:
                                pair_iter.append(pair)
        seen_pairs = set()
        for atom_i, atom_j in pair_iter:
            pair = (min(atom_i, atom_j), max(atom_i, atom_j))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            target = float(np.linalg.norm(initial_positions[pair[0]] - initial_positions[pair[1]]))
            bond = mol.GetBondBetweenAtoms(pair[0], pair[1])
            if bond is not None:
                weight = 1.60
            else:
                weight = 0.85 if unit_planar else 0.55
            rigid_pair_targets.append((pair[0], pair[1], target, weight))

    def _is_planar_like_atom(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_units: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
    for unit in rigid_units:
        unit_set = set(unit)
        movable_overlap = unit_set & movable_atoms
        if len(movable_overlap) < 2:
            continue
        donor_overlap = [
            donor_idx for donor_idx in donor_indices
            if donor_idx in unit_set and _is_planar_like_atom(mol.GetAtomWithIdx(donor_idx))
        ]
        if not donor_overlap:
            continue
        planar_atoms = tuple(
            atom_idx for atom_idx in unit
            if _is_planar_like_atom(mol.GetAtomWithIdx(atom_idx))
        )
        if len(planar_atoms) < 4:
            continue
        pts = np.asarray([initial_positions[atom_idx] for atom_idx in planar_atoms], dtype=float)
        centroid = pts.mean(axis=0)
        try:
            _u, _s, vt = np.linalg.svd(pts - centroid)
        except Exception:
            continue
        normal = vt[-1]
        rms_offset = float(np.sqrt(np.mean(((pts - centroid) @ normal) ** 2)))
        if rms_offset <= 0.18:
            planar_units.append((tuple(sorted(unit)), planar_atoms))

    metal_coplanar_units: List[Tuple[Tuple[int, ...], float]] = []
    seen_metal_coplanar_units: set = set()
    for unit in rigid_units:
        unit_set = set(unit)
        donor_overlap = [
            donor_idx for donor_idx in donor_indices
            if donor_idx in unit_set and _is_planar_like_atom(mol.GetAtomWithIdx(donor_idx))
        ]
        if not donor_overlap:
            continue
        planar_atoms = tuple(
            atom_idx for atom_idx in unit
            if _is_planar_like_atom(mol.GetAtomWithIdx(atom_idx))
        )
        if len(planar_atoms) < 4:
            continue
        planar_key = tuple(sorted(planar_atoms))
        if planar_key in seen_metal_coplanar_units:
            continue
        seen_metal_coplanar_units.add(planar_key)
        weight = 34.0 + 10.0 * len(donor_overlap)
        if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'N' for donor_idx in donor_overlap):
            weight += 12.0
        if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O' for donor_idx in donor_overlap):
            weight += 12.0
        if any(mol.GetAtomWithIdx(donor_idx).GetIsAromatic() for donor_idx in donor_overlap):
            weight += 10.0
        metal_coplanar_units.append((planar_key, weight))

    def _gp(atom_idx: int):
        pos = conf.GetAtomPosition(atom_idx)
        return np.array([pos.x, pos.y, pos.z], dtype=float)

    def _sp(atom_idx: int, vec):
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))

    def _score() -> float:
        metal_pos = _gp(metal_idx)
        coord_map = {
            atom_idx: _gp(atom_idx)
            for atom_idx in module_atoms
            if atom_idx != metal_idx
        }
        score = 0.0
        for begin_idx, end_idx, target_len in bond_targets:
            dist = float(np.linalg.norm(_gp(begin_idx) - _gp(end_idx)))
            score += 0.45 * (dist - target_len) ** 2
        for donor_idx, target_len, weight in ml_targets:
            dist = float(np.linalg.norm(_gp(donor_idx) - metal_pos))
            score += weight * (dist - target_len) ** 2
        for donor_i, donor_j, target_len in donor_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            score += 0.20 * (dist - target_len) ** 2
        for donor_i, donor_j, target_len, weight in cn4_geometry_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            score += weight * (dist - target_len) ** 2
        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            for donor_idx in donor_indices:
                offset = float(np.dot(_gp(donor_idx) - metal_pos, normal))
                score += 0.35 * offset * offset
        for planar_atoms, weight in metal_coplanar_units:
            pts = []
            for atom_idx in planar_atoms:
                arr = coord_map.get(atom_idx)
                if arr is None:
                    arr = _gp(atom_idx)
                pts.append(np.asarray(arr, dtype=float))
            if len(pts) < 4:
                continue
            pts_arr = np.asarray(pts, dtype=float)
            centroid = pts_arr.mean(axis=0)
            try:
                _u, _s, vt = np.linalg.svd(pts_arr - centroid)
            except Exception:
                continue
            normal = vt[-1]
            normal_norm = float(np.linalg.norm(normal))
            if normal_norm < 1e-12:
                continue
            normal = normal / normal_norm
            plane_offset = float(np.dot(metal_pos - centroid, normal))
            score += weight * plane_offset * plane_offset
            if abs(plane_offset) > 0.28:
                score += 75.0 * (abs(plane_offset) - 0.28) ** 2
        seen_fragment_planes: set = set()
        for donor_idx in donor_indices:
            fragment = donor_to_fragment.get(donor_idx)
            if fragment is None:
                continue
            fragment_key = tuple(sorted(fragment.atom_indices))
            if fragment_key in seen_fragment_planes:
                continue
            seen_fragment_planes.add(fragment_key)
            fragment_donors = [
                idx for idx in fragment.donor_atom_indices
                if idx in donor_set
            ]
            if not fragment_donors:
                continue
            score += _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                fragment_donors,
                metal_pos,
                coord_map,
            )
        for atom_i, atom_j, target_len, weight in rigid_pair_targets:
            dist = float(np.linalg.norm(_gp(atom_i) - _gp(atom_j)))
            score += weight * (dist - target_len) ** 2
        for other_metal_idx, min_dist in metal_repulsion_targets:
            dist = float(np.linalg.norm(_gp(other_metal_idx) - metal_pos))
            if dist < min_dist:
                score += 16.0 * (min_dist - dist) ** 2
        module_atom_list = sorted(module_atoms - {metal_idx})
        for idx_i, atom_i in enumerate(module_atom_list):
            frag_i = atom_to_fragment_idx.get(atom_i, -1)
            sym_i = mol.GetAtomWithIdx(atom_i).GetSymbol()
            if sym_i == 'H' or sym_i in _METAL_SET:
                continue
            pos_i = _gp(atom_i)
            for atom_j in module_atom_list[idx_i + 1:]:
                if (min(atom_i, atom_j), max(atom_i, atom_j)) in bonded_pairs:
                    continue
                frag_j = atom_to_fragment_idx.get(atom_j, -1)
                if frag_i == frag_j:
                    continue
                sym_j = mol.GetAtomWithIdx(atom_j).GetSymbol()
                if sym_j == 'H' or sym_j in _METAL_SET:
                    continue
                pos_j = _gp(atom_j)
                dist = float(np.linalg.norm(pos_i - pos_j))
                min_dist = max(
                    1.78,
                    _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.36,
                )
                if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                    min_dist = max(min_dist, 1.90)
                if dist < min_dist:
                    score += 12.0 * (min_dist - dist) ** 2
        score += _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(module_atoms - {metal_idx}),
            donor_indices,
            coord_map,
        )
        score += _secondary_fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            donor_indices,
            donor_to_fragment,
            coord_map,
        )
        return score

    def _cn4_geometry_deviation() -> float:
        if len(donor_indices) != 4 or geometry_code not in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES:
            return 0.0
        deviation = 0.0
        metal_pos = _gp(metal_idx)
        for donor_i, donor_j, target_len, _weight in cn4_geometry_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            deviation += (dist - target_len) ** 2
        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            for donor_idx in donor_indices:
                offset = float(np.dot(_gp(donor_idx) - metal_pos, normal))
                deviation += 0.5 * offset * offset
        return deviation

    best_positions = {atom_idx: vec.copy() for atom_idx, vec in initial_positions.items()}
    start_score = _score()
    best_score = start_score
    start_geometry_deviation = _cn4_geometry_deviation()
    best_geometry_deviation = start_geometry_deviation
    rng = np.random.default_rng(int(1009 * (metal_idx + 1) + len(module_atoms)))

    alpha_atoms: set = set()
    beta_atoms: set = set()
    for donor_idx in donor_indices:
        atom = mol.GetAtomWithIdx(donor_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in module_atoms and nbr_idx not in donor_set and nbr.GetAtomicNum() > 1:
                alpha_atoms.add(nbr_idx)
    for atom_idx in list(alpha_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if (
                nbr_idx in module_atoms
                and nbr_idx not in donor_set
                and nbr_idx not in alpha_atoms
                and nbr.GetAtomicNum() > 1
            ):
                beta_atoms.add(nbr_idx)

    for _pass in range(120):
        displacements: Dict[int, np.ndarray] = {
            atom_idx: np.zeros(3, dtype=float) for atom_idx in movable_atoms
        }
        active = False

        for begin_idx, end_idx, target_len in bond_targets:
            begin_pos = _gp(begin_idx)
            end_pos = _gp(end_idx)
            diff = end_pos - begin_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            strength = 0.42 if ((begin_idx in movable_atoms) ^ (end_idx in movable_atoms)) else 0.24
            move = strength * delta * unit
            if begin_idx in movable_atoms and end_idx in movable_atoms:
                displacements[begin_idx] = displacements[begin_idx] + 0.5 * move
                displacements[end_idx] = displacements[end_idx] - 0.5 * move
            elif begin_idx in movable_atoms:
                displacements[begin_idx] = displacements[begin_idx] + move
            elif end_idx in movable_atoms:
                displacements[end_idx] = displacements[end_idx] - move

        metal_pos = _gp(metal_idx)
        for donor_idx, target_len, _weight in ml_targets:
            donor_pos = _gp(donor_idx)
            diff = donor_pos - metal_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            strength = 0.40 if donor_idx not in movable_atoms else 0.22
            move = strength * delta * unit
            if metal_idx in movable_atoms and donor_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] + 0.55 * move
                displacements[donor_idx] = displacements[donor_idx] - 0.45 * move
            elif metal_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] + move
            elif donor_idx in movable_atoms:
                displacements[donor_idx] = displacements[donor_idx] - move

        for donor_i, donor_j, target_len in donor_pair_targets:
            pos_i = _gp(donor_i)
            pos_j = _gp(donor_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.015:
                continue
            active = True
            move = 0.12 * delta * unit
            if donor_i in movable_atoms and donor_j in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + 0.5 * move
                displacements[donor_j] = displacements[donor_j] - 0.5 * move
            elif donor_i in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + move
            elif donor_j in movable_atoms:
                displacements[donor_j] = displacements[donor_j] - move

        for donor_i, donor_j, target_len, weight in cn4_geometry_pair_targets:
            pos_i = _gp(donor_i)
            pos_j = _gp(donor_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            move = min(0.18, 0.08 + 0.10 * weight) * delta * unit
            if donor_i in movable_atoms and donor_j in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + 0.5 * move
                displacements[donor_j] = displacements[donor_j] - 0.5 * move
            elif donor_i in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + move
            elif donor_j in movable_atoms:
                displacements[donor_j] = displacements[donor_j] - move

        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            metal_pos = _gp(metal_idx)
            for donor_idx in donor_indices:
                donor_pos = _gp(donor_idx)
                offset = float(np.dot(donor_pos - metal_pos, normal))
                if abs(offset) < 0.01:
                    continue
                active = True
                move = 0.16 * offset * normal
                if donor_idx in movable_atoms and metal_idx in movable_atoms:
                    displacements[donor_idx] = displacements[donor_idx] - 0.55 * move
                    displacements[metal_idx] = displacements[metal_idx] + 0.45 * move
                elif donor_idx in movable_atoms:
                    displacements[donor_idx] = displacements[donor_idx] - move
                elif metal_idx in movable_atoms:
                    displacements[metal_idx] = displacements[metal_idx] + move

        for atom_i, atom_j, target_len, weight in rigid_pair_targets:
            pos_i = _gp(atom_i)
            pos_j = _gp(atom_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.008:
                continue
            active = True
            move = min(0.18, 0.10 + 0.04 * weight) * delta * unit
            if atom_i in movable_atoms and atom_j in movable_atoms:
                displacements[atom_i] = displacements[atom_i] + 0.5 * move
                displacements[atom_j] = displacements[atom_j] - 0.5 * move
            elif atom_i in movable_atoms:
                displacements[atom_i] = displacements[atom_i] + move
            elif atom_j in movable_atoms:
                displacements[atom_j] = displacements[atom_j] - move

        metal_pos = _gp(metal_idx)
        for other_metal_idx, min_dist in metal_repulsion_targets:
            other_pos = _gp(other_metal_idx)
            diff = metal_pos - other_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            if dist >= min_dist:
                continue
            active = True
            unit = diff / max(dist, 1e-12)
            gap = min_dist - dist
            displacements[metal_idx] = displacements[metal_idx] + 0.70 * gap * unit

        for atom_idx in sorted(module_atoms - donor_set - {metal_idx}):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = _gp(atom_idx)
            diff = atom_pos - metal_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            min_allowed = max(1.70, 0.92 * float(_get_ml_bond_length(metal_sym, atom.GetSymbol())))
            if atom_idx in alpha_atoms:
                if _is_planar_like_atom(atom):
                    min_allowed = max(min_allowed, 2.18)
                else:
                    min_allowed = max(min_allowed, 2.05)
            elif atom_idx in beta_atoms:
                if _is_planar_like_atom(atom):
                    min_allowed = max(min_allowed, 1.96)
                else:
                    min_allowed = max(min_allowed, 1.85)
            if dist >= min_allowed:
                continue
            active = True
            gap = min_allowed - dist
            move = 0.60 * gap * unit
            if metal_idx in movable_atoms and atom_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] - 0.65 * move
                displacements[atom_idx] = displacements[atom_idx] + 0.35 * move
            elif metal_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] - move
            elif atom_idx in movable_atoms:
                displacements[atom_idx] = displacements[atom_idx] + move

        for atom_i in sorted(module_atoms):
            for atom_j in sorted(module_atoms):
                if atom_j <= atom_i:
                    continue
                if (min(atom_i, atom_j), max(atom_i, atom_j)) in bonded_pairs:
                    continue
                if atom_i not in movable_atoms and atom_j not in movable_atoms:
                    continue
                pos_i = _gp(atom_i)
                pos_j = _gp(atom_j)
                diff = pos_j - pos_i
                dist = float(np.linalg.norm(diff))
                if dist < 1e-8:
                    diff = rng.standard_normal(3)
                    dist = float(np.linalg.norm(diff))
                unit = diff / max(dist, 1e-12)
                sym_i = mol.GetAtomWithIdx(atom_i).GetSymbol()
                sym_j = mol.GetAtomWithIdx(atom_j).GetSymbol()
                if sym_i in _METAL_SET or sym_j in _METAL_SET:
                    min_dist = 2.0
                elif sym_i == 'H' and sym_j == 'H':
                    min_dist = 1.5
                elif sym_i == 'H' or sym_j == 'H':
                    min_dist = 1.0
                else:
                    min_dist = 1.18
                    frag_i = atom_to_fragment_idx.get(atom_i, -1)
                    frag_j = atom_to_fragment_idx.get(atom_j, -1)
                    if frag_i != frag_j:
                        min_dist = max(
                            min_dist,
                            _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.36,
                            1.78,
                        )
                        if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                            min_dist = max(min_dist, 1.90)
                if dist >= min_dist:
                    continue
                active = True
                gap = min_dist - dist
                move_scale = 0.38 if atom_to_fragment_idx.get(atom_i, -1) != atom_to_fragment_idx.get(atom_j, -1) else 0.18
                move = move_scale * gap * unit
                if atom_i in movable_atoms and atom_j in movable_atoms:
                    displacements[atom_i] = displacements[atom_i] - 0.5 * move
                    displacements[atom_j] = displacements[atom_j] + 0.5 * move
                elif atom_i in movable_atoms:
                    displacements[atom_i] = displacements[atom_i] - move
                elif atom_j in movable_atoms:
                    displacements[atom_j] = displacements[atom_j] + move

        if not active:
            break

        capped_displacements: Dict[int, np.ndarray] = {}
        for atom_idx, disp in displacements.items():
            norm = float(np.linalg.norm(disp))
            if norm < 1e-8:
                continue
            if norm > 0.22:
                disp = disp * (0.22 / norm)
            capped_displacements[atom_idx] = disp

        for unit, planar_atoms in planar_units:
            if not any(atom_idx in capped_displacements for atom_idx in unit):
                continue
            unit_positions = {
                atom_idx: _gp(atom_idx)
                for atom_idx in unit
            }
            unit_displacements = {
                atom_idx: capped_displacements.get(atom_idx, np.zeros(3, dtype=float))
                for atom_idx in unit
            }
            rigid_projected = _project_displacements_to_rigid_body(
                unit_positions,
                unit_displacements,
                list(unit),
            )
            for atom_idx in unit:
                if atom_idx in capped_displacements and atom_idx in rigid_projected:
                    capped_displacements[atom_idx] = np.asarray(rigid_projected[atom_idx], dtype=float)

        for atom_idx, disp in capped_displacements.items():
            _sp(atom_idx, _gp(atom_idx) + disp)

        trial_score = _score()
        if trial_score + 1e-9 < best_score:
            best_score = trial_score
            best_geometry_deviation = _cn4_geometry_deviation()
            best_positions = {
                atom_idx: _gp(atom_idx).copy()
                for atom_idx in movable_atoms
            }

    if best_score + 1e-6 >= start_score:
        for atom_idx, vec in initial_positions.items():
            _sp(atom_idx, vec)
        return set()

    if len(donor_indices) == 4 and geometry_code in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES:
        allowed_geometry_deviation = max(
            start_geometry_deviation + 0.03,
            1.35 * start_geometry_deviation + 1e-6,
        )
        if best_geometry_deviation > allowed_geometry_deviation:
            for atom_idx, vec in initial_positions.items():
                _sp(atom_idx, vec)
            logger.info(
                "Hybrid hapto rejected local secondary metal module optimization for %s%d: CN=4 %s geometry deviation %.3f -> %.3f",
                metal_sym,
                metal_idx,
                geometry_code,
                start_geometry_deviation,
                best_geometry_deviation,
            )
            return set()

    for atom_idx, vec in best_positions.items():
        _sp(atom_idx, vec)

    logger.info(
        "Hybrid hapto locally optimized secondary metal module %s%d: %.3f -> %.3f",
        metal_sym,
        metal_idx,
        start_score,
        best_score,
    )
    return set(movable_atoms)


def _assemble_secondary_metal_coordination_modules(
    mol,
    decomposition: _HybridHaptoDecomposition,
    hapto_groups: List[Tuple[int, List[int]]],
    fit_variant_plan: Optional[Dict[int, int]] = None,
) -> Tuple[set, set, set]:
    """Build secondary metal coordination spheres from donor fragments."""
    if not RDKIT_AVAILABLE or mol is None or decomposition is None or not hapto_groups:
        return set(), set(), set()
    try:
        import numpy as np
    except ImportError:
        return set(), set(), set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set(), set(), set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    placed_secondary: set = set()
    module_atoms: set = set()
    relaxable_atoms: set = set()

    donor_to_fragment: Dict[int, _HybridHaptoFragment] = {}
    embedded_fragment_cache: Dict[Tuple[int, ...], object] = {}
    for fragment in decomposition.fragments:
        for donor_idx in fragment.donor_atom_indices:
            donor_to_fragment[donor_idx] = fragment

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        metal_sym = atom.GetSymbol()
        if metal_sym not in _METAL_SET or metal_idx in hapto_metals:
            continue
        variant_rank = max(0, int((fit_variant_plan or {}).get(metal_idx, 0)))

        donor_indices = [
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() not in _METAL_SET and nbr.GetAtomicNum() > 1
        ]
        if len(donor_indices) == 0:
            continue
        if len(donor_indices) == 1:
            # Single-donor fallback: place metal along donor→COM direction
            donor_idx = donor_indices[0]
            donor_pos = conf.GetAtomPosition(donor_idx)
            dp = np.array([donor_pos.x, donor_pos.y, donor_pos.z])
            all_pos = []
            for a2 in mol.GetAtoms():
                if a2.GetSymbol() not in _METAL_SET:
                    p2 = conf.GetAtomPosition(a2.GetIdx())
                    all_pos.append(np.array([p2.x, p2.y, p2.z]))
            com = np.mean(all_pos, axis=0) if all_pos else dp
            direction = dp - com
            norm = np.linalg.norm(direction)
            if norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
                norm = 1.0
            direction = direction / norm
            try:
                bond_len = _secondary_donor_target_length(mol, metal_sym, donor_idx)
            except Exception:
                bond_len = 2.1
            metal_pos = dp + direction * bond_len
            conf.SetAtomPosition(metal_idx, Point3D(*metal_pos))
            placed_secondary.add(metal_idx)
            module_atoms.add(metal_idx)
            module_atoms.add(donor_idx)
            continue

        donor_positions = []
        donor_symbols = []
        donor_target_lengths = []
        donor_fit_weights = []
        heterodonor_count = 0
        carbon_donor_count = 0
        for donor_idx in donor_indices:
            pos = conf.GetAtomPosition(donor_idx)
            donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            donor_symbols.append(donor_sym)
            donor_target_lengths.append(_secondary_donor_target_length(mol, metal_sym, donor_idx))
            donor_fit_weights.append(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
            if donor_sym == 'C':
                carbon_donor_count += 1
            elif donor_sym in {'N', 'O', 'S', 'P'}:
                heterodonor_count += 1

        movable_fragments: List[_HybridHaptoFragment] = []
        seen_fragments: set = set()
        anchored_slots: List[int] = []
        supplemental_fit_slots: List[int] = []
        bite_distance_constraints: List[Tuple[int, int, float]] = []
        for slot_idx, donor_idx in enumerate(donor_indices):
            fragment = donor_to_fragment.get(donor_idx)
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            if fragment is None:
                anchored_slots.append(slot_idx)
                continue
            if heterodonor_count and carbon_donor_count and donor_sym == 'C':
                if fragment.bridging_donor_indices or fragment.use_scaffold_only:
                    donor_fit_weights[slot_idx] *= 0.72
                elif len(fragment.hapto_group_ids) > 0:
                    donor_fit_weights[slot_idx] *= 0.82
            if donor_sym in {'N', 'O', 'S', 'P'}:
                if len(fragment.donor_atom_indices) >= 2:
                    donor_fit_weights[slot_idx] *= 1.12
                if any(
                    mol.GetAtomWithIdx(atom_idx).GetIsAromatic()
                    or mol.GetAtomWithIdx(atom_idx).GetHybridization() in {
                        Chem.rdchem.HybridizationType.SP,
                        Chem.rdchem.HybridizationType.SP2,
                    }
                    for atom_idx in fragment.atom_indices
                    if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                    and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                ):
                    donor_fit_weights[slot_idx] *= 1.08
            if _secondary_module_fragment_is_movable(fragment, metal_idx):
                frag_key = tuple(fragment.atom_indices)
                if frag_key not in seen_fragments:
                    movable_fragments.append(fragment)
                    seen_fragments.add(frag_key)
                supplemental_fit_slots.append(slot_idx)
                continue
            if (
                fragment.use_scaffold_only
                and _find_scaffold_secondary_branch(
                    mol,
                    fragment,
                    hapto_groups,
                    donor_idx,
                ) is not None
            ):
                supplemental_fit_slots.append(slot_idx)
                continue
            anchored_slots.append(slot_idx)

        fit_slot_indices: Optional[List[int]]
        if anchored_slots:
            fit_slot_indices = list(anchored_slots)
            for slot_idx in supplemental_fit_slots:
                if slot_idx not in fit_slot_indices:
                    fit_slot_indices.append(slot_idx)
                if len(fit_slot_indices) >= 2:
                    break
        else:
            fit_slot_indices = None

        has_multidentate_movable_fragment = any(
            len([
                donor_idx
                for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_indices
            ]) >= 2
            for fragment in movable_fragments
        )
        has_secondary_bridge_branch = any(
            fragment.use_scaffold_only and metal_idx in fragment.metal_neighbor_indices
            for fragment in decomposition.fragments
        )
        if (
            len(donor_indices) == 4
            and has_multidentate_movable_fragment
            and has_secondary_bridge_branch
        ):
            fit_slot_indices = None

        for fragment in movable_fragments:
            frag_key = tuple(fragment.atom_indices)
            embedded_fragment = embedded_fragment_cache.get(frag_key)
            if embedded_fragment is None:
                embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
                if embedded_fragment is not None:
                    embedded_fragment_cache[frag_key] = embedded_fragment
            if embedded_fragment is None:
                continue
            try:
                frag_conf = embedded_fragment.GetConformer()
            except Exception:
                continue
            frag_donors = [
                donor_idx for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_indices
            ]
            if len(frag_donors) < 2:
                continue
            donor_slots = {
                donor_idx: donor_indices.index(donor_idx)
                for donor_idx in frag_donors
            }
            for i in range(len(frag_donors)):
                for j in range(i + 1, len(frag_donors)):
                    donor_i = frag_donors[i]
                    donor_j = frag_donors[j]
                    frag_i = fragment.original_to_fragment.get(donor_i)
                    frag_j = fragment.original_to_fragment.get(donor_j)
                    if frag_i is None or frag_j is None:
                        continue
                    pos_i = frag_conf.GetAtomPosition(frag_i)
                    pos_j = frag_conf.GetAtomPosition(frag_j)
                    expected_dist = math.sqrt(
                        (pos_i.x - pos_j.x) ** 2
                        + (pos_i.y - pos_j.y) ** 2
                        + (pos_i.z - pos_j.z) ** 2
                    )
                    bite_distance_constraints.append(
                        (donor_slots[donor_i], donor_slots[donor_j], float(expected_dist))
                    )

        fit_candidates = _enumerate_secondary_metal_geometry_fits(
            donor_positions,
            donor_symbols,
            metal_sym,
            constrained_indices=fit_slot_indices,
            bite_distance_constraints=bite_distance_constraints,
            donor_target_lengths=donor_target_lengths,
            donor_fit_weights=donor_fit_weights,
            max_candidates=max(variant_rank + 1, 4),
        )
        fit = fit_candidates[min(variant_rank, len(fit_candidates) - 1)] if fit_candidates else None
        if fit is None:
            fallback_candidates = _enumerate_secondary_metal_position_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            fallback = fallback_candidates[min(variant_rank, len(fallback_candidates) - 1)] if fallback_candidates else None
            if fallback is None:
                continue
            metal_pos, geom_code, _fit_score = fallback
            conf.SetAtomPosition(
                metal_idx,
                Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
            )
            placed_secondary.add(metal_idx)
            module_atoms.add(metal_idx)
            logger.info(
                "Hybrid hapto built secondary metal %s%d via fallback %s fit to %d donor(s)",
                metal_sym,
                metal_idx,
                geom_code,
                len(donor_indices),
            )
            continue

        metal_pos, geom_code, target_positions, _fit_score = fit
        donor_target_map: Dict[int, np.ndarray] = {
            donor_idx: np.asarray(target_positions[slot_idx], dtype=float)
            for slot_idx, donor_idx in enumerate(donor_indices)
        }
        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )

        reoriented_scaffold_branches = 0
        for fragment in decomposition.fragments:
            if not fragment.use_scaffold_only or metal_idx not in fragment.metal_neighbor_indices:
                continue
            for donor_idx in fragment.donor_atom_indices:
                if donor_idx not in donor_target_map:
                    continue
                branch_info = _find_scaffold_secondary_branch(
                    mol,
                    fragment,
                    hapto_groups,
                    donor_idx,
                )
                if branch_info is None:
                    continue
                try:
                    if _reorient_scaffold_secondary_branch(
                        mol,
                        fragment,
                        hapto_groups,
                        metal_idx,
                        donor_idx,
                        donor_target_map[donor_idx],
                    ):
                        reoriented_scaffold_branches += 1
                        pivot_idx, branch_atoms = branch_info
                        module_atoms.add(pivot_idx)
                        module_atoms.update(branch_atoms)
                        relaxable_atoms.update(branch_atoms)
                except Exception:
                    continue

        if reoriented_scaffold_branches:
            donor_positions = []
            for donor_idx in donor_indices:
                pos = conf.GetAtomPosition(donor_idx)
                donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            refit_candidates = _enumerate_secondary_metal_geometry_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                constrained_indices=fit_slot_indices,
                bite_distance_constraints=bite_distance_constraints,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            refit = refit_candidates[min(variant_rank, len(refit_candidates) - 1)] if refit_candidates else None
            if refit is not None:
                metal_pos, geom_code, target_positions, _fit_score = refit
                donor_target_map = {
                    donor_idx: np.asarray(target_positions[slot_idx], dtype=float)
                    for slot_idx, donor_idx in enumerate(donor_indices)
                }
                conf.SetAtomPosition(
                    metal_idx,
                    Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
                )

        moved_fragments = 0
        has_multidentate_movable_fragment = False
        rigid_reference_fragments: List[Tuple[_HybridHaptoFragment, Dict[int, np.ndarray]]] = []
        oo_relief_fragments: List[_HybridHaptoFragment] = []
        for fragment in movable_fragments:
            frag_key = tuple(fragment.atom_indices)
            embedded_fragment = embedded_fragment_cache.get(frag_key)
            if embedded_fragment is None:
                embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
                if embedded_fragment is not None:
                    embedded_fragment_cache[frag_key] = embedded_fragment
            if embedded_fragment is None:
                continue
            frag_target_indices = [
                donor_idx for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_target_map
            ]
            if not frag_target_indices:
                continue
            if len(frag_target_indices) >= 2:
                has_multidentate_movable_fragment = True
            if _align_hybrid_fragment_to_targets(
                mol,
                fragment,
                embedded_fragment,
                frag_target_indices,
                donor_target_map,
                reference_center=metal_pos,
                metal_idx=metal_idx,
                exact_target_count=0,
            ):
                moved_fragments += 1
                module_atoms.update(fragment.atom_indices)
                relaxable_atoms.update(fragment.atom_indices)
                keep_rigid_pose = (
                    len(donor_indices) == 4
                    and geom_code in _HAPTO_SECONDARY_CN4_GEOMETRIES
                    and has_secondary_bridge_branch
                    and _secondary_fragment_prefers_rigid_pose(
                        mol,
                        fragment,
                        frag_target_indices,
                    )
                )
                if not keep_rigid_pose:
                    try:
                        _optimize_secondary_fragment_pose(
                            mol,
                            fragment,
                            metal_idx,
                        donor_target_map,
                    )
                    except Exception:
                        pass
                else:
                    rigid_reference_fragments.append(
                        (
                            fragment,
                            {
                                atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
                                for atom_idx in fragment.atom_indices
                            },
                        )
                    )
                if keep_rigid_pose and sum(
                    1 for donor_idx in frag_target_indices
                    if mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O'
                ) >= 2:
                    oo_relief_fragments.append(fragment)

        allow_post_fragment_refit = not (
            len(donor_indices) == 4
            and geom_code in {'SQ', 'TET'}
            and has_multidentate_movable_fragment
            and has_secondary_bridge_branch
        )

        if moved_fragments and allow_post_fragment_refit:
            def _module_error(candidate_metal_pos) -> float:
                err = 0.0
                metal_arr = np.asarray(candidate_metal_pos, dtype=float)
                coord_map: Dict[int, np.ndarray] = {}
                for fragment in decomposition.fragments:
                    if metal_idx not in fragment.metal_neighbor_indices:
                        continue
                    for atom_idx in fragment.atom_indices:
                        if atom_idx in coord_map:
                            continue
                        pos = conf.GetAtomPosition(atom_idx)
                        coord_map[atom_idx] = np.array([pos.x, pos.y, pos.z], dtype=float)
                for donor_idx in donor_indices:
                    donor_pos = np.array(conf.GetAtomPosition(donor_idx), dtype=float)
                    target_len = float(_secondary_donor_target_length(mol, metal_sym, donor_idx))
                    weight = 2.0 * float(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
                    err += weight * (float(np.linalg.norm(donor_pos - metal_arr)) - target_len) ** 2
                seen_fragments: set = set()
                for donor_idx in donor_indices:
                    fragment = donor_to_fragment.get(donor_idx)
                    if fragment is None:
                        continue
                    frag_key = tuple(fragment.atom_indices)
                    if frag_key in seen_fragments:
                        continue
                    seen_fragments.add(frag_key)
                    fragment_donors = [
                        idx for idx in fragment.donor_atom_indices
                        if idx in set(donor_indices)
                    ]
                    if not fragment_donors:
                        continue
                    err += _secondary_non_donor_contact_penalty(
                        mol,
                        metal_idx,
                        metal_arr,
                        fragment.atom_indices,
                        fragment_donors,
                        coord_map,
                    )
                    err += _planar_fragment_metal_coplanarity_penalty(
                        mol,
                        fragment.atom_indices,
                        fragment_donors,
                        metal_arr,
                        coord_map,
                    )
                    err += _planar_fragment_donor_approach_penalty(
                        mol,
                        fragment.atom_indices,
                        fragment_donors,
                        metal_arr,
                        coord_map,
                    )
                    if len(fragment_donors) >= 2:
                        err += _fragment_bite_direction_penalty(
                            mol,
                            fragment.atom_indices,
                            fragment_donors,
                            coord_map,
                            metal_arr,
                        )
                return err

            current_metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
            best_metal_pos = current_metal_pos
            best_metal_err = _module_error(current_metal_pos)
            donor_positions = []
            for donor_idx in donor_indices:
                pos = conf.GetAtomPosition(donor_idx)
                donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            refit_candidates = _enumerate_secondary_metal_geometry_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                constrained_indices=fit_slot_indices,
                bite_distance_constraints=bite_distance_constraints,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            refit = refit_candidates[min(variant_rank, len(refit_candidates) - 1)] if refit_candidates else None
            if refit is not None:
                refit_metal_pos, refit_geom_code, _target_positions, _refit_score = refit
                refit_err = _module_error(refit_metal_pos)
                if refit_err + 1e-9 < best_metal_err:
                    best_metal_pos = np.asarray(refit_metal_pos, dtype=float)
                    best_metal_err = refit_err
                    geom_code = refit_geom_code
            metal_pos = best_metal_pos
        elif moved_fragments and not allow_post_fragment_refit:
            logger.info(
                "Hybrid hapto kept pre-fragment-fit metal position for %s%d: skipped post-fragment refit for bridged CN=4 module",
                metal_sym,
                metal_idx,
            )

        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )
        optimized_atoms: set = set()
        try:
            optimized_atoms = _optimize_secondary_metal_module_local(
                mol,
                decomposition.fragments,
                hapto_groups,
                metal_idx,
                donor_indices,
                donor_to_fragment,
                geometry_code=geom_code,
                geometry_metal_pos=metal_pos,
                geometry_target_map=donor_target_map,
            )
        except Exception as e:
            logger.debug("Hybrid hapto local secondary-metal optimization skipped: %s", e)
            optimized_atoms = set()
        if optimized_atoms:
            module_atoms.update(optimized_atoms)
            relaxable_atoms.update(optimized_atoms - {metal_idx})
            pos = conf.GetAtomPosition(metal_idx)
            metal_pos = np.array([pos.x, pos.y, pos.z], dtype=float)
        restored_rigid_count = 0
        for fragment, reference_coords in rigid_reference_fragments:
            try:
                if _restore_secondary_rigid_fragment_geometry(
                    mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                    reference_coords,
                ):
                    restored_rigid_count += 1
            except Exception:
                continue
        oo_relief_count = 0
        for fragment in oo_relief_fragments:
            try:
                if _relieve_secondary_oo_chelate_contacts(
                    mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                ):
                    oo_relief_count += 1
            except Exception:
                continue
        placed_secondary.add(metal_idx)
        module_atoms.add(metal_idx)
        module_atoms.update(donor_indices)
        logger.info(
            "Hybrid hapto built secondary metal module %s%d via %s fit to %d donor(s); moved %d fragment(s), reoriented %d scaffold branch(es), locally optimized %d atom(s), restored %d rigid fragment(s), relieved %d O,O chelate(s)",
            metal_sym,
            metal_idx,
            geom_code,
            len(donor_indices),
            moved_fragments,
            reoriented_scaffold_branches,
            len(optimized_atoms),
            restored_rigid_count,
            oo_relief_count,
        )

    return placed_secondary, module_atoms, relaxable_atoms


def _build_primary_organometal_hapto_module(
    mol,
    decomposition: Optional[_HybridHaptoDecomposition],
    hapto_groups: List[Tuple[int, List[int]]],
    module: Optional[_PrimaryOrganometalModule],
    preview_only: bool = False,
    preview_store: Optional[List[Tuple[str, str]]] = None,
):
    """Build a single-metal hapto organometal module before falling back to the legacy path."""
    if (
        not RDKIT_AVAILABLE
        or mol is None
        or decomposition is None
        or module is None
        or not hapto_groups
    ):
        return None
    try:
        import numpy as np
    except ImportError:
        return None

    scaffold_mol = Chem.Mol(mol)
    if not _build_hapto_scaffold(scaffold_mol, hapto_groups):
        return None

    try:
        _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    try:
        conf = scaffold_mol.GetConformer(0)
    except Exception:
        return None

    metal_idx = int(module.metal_idx)
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    donor_target_map: Dict[int, np.ndarray] = {}
    for donor_idx in module.donor_atom_indices:
        pos = conf.GetAtomPosition(donor_idx)
        donor_target_map[donor_idx] = np.array([pos.x, pos.y, pos.z], dtype=float)

    trusted_atoms: set = {
        atom_idx
        for group_idx in module.hapto_group_ids
        if 0 <= int(group_idx) < len(hapto_groups)
        for atom_idx in hapto_groups[int(group_idx)][1]
    }
    aligned_fragments = 0
    aligned_correlated = 0

    target_fragment_indices = set(module.correlated_fragment_indices) | set(module.terminal_fragment_indices)
    for frag_idx, fragment in enumerate(decomposition.fragments):
        if frag_idx not in target_fragment_indices:
            if any(group_id in module.hapto_group_ids for group_id in fragment.hapto_group_ids):
                trusted_atoms.update(fragment.atom_indices)
            continue

        embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded_fragment is None:
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        frag_target_indices = [
            donor_idx for donor_idx in fragment.donor_atom_indices
            if donor_idx in donor_target_map
        ]
        if not frag_target_indices:
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        if not _align_hybrid_fragment_to_targets(
            scaffold_mol,
            fragment,
            embedded_fragment,
            frag_target_indices,
            donor_target_map,
            reference_center=metal_pos,
            metal_idx=metal_idx,
            exact_target_count=0,
        ):
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        aligned_fragments += 1
        if frag_idx in module.correlated_fragment_indices:
            aligned_correlated += 1
        trusted_atoms.update(fragment.atom_indices)

        if (
            len(frag_target_indices) >= 2
            and not _secondary_fragment_prefers_rigid_pose(
                scaffold_mol,
                fragment,
                frag_target_indices,
            )
        ):
            try:
                _optimize_secondary_fragment_pose(
                    scaffold_mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                )
            except Exception:
                pass

    if aligned_correlated < len(module.correlated_fragment_indices):
        return None
    if aligned_fragments == 0:
        return None

    try:
        _propagate_non_hapto_atoms(
            scaffold_mol,
            0,
            hapto_groups,
            extra_fixed_indices=trusted_atoms | _hapto_primary_donor_indices(scaffold_mol, hapto_groups),
        )
    except Exception:
        pass

    if preview_store is not None:
        try:
            preview_xyz = _mol_to_xyz(scaffold_mol)
            preview_entry = (preview_xyz, 'quick-primary-preview')
            if preview_entry not in preview_store:
                preview_store.append(preview_entry)
        except Exception:
            pass

    if (
        not preview_only
        and not _primary_organometal_module_quality_ok(
            scaffold_mol,
            decomposition,
            module,
        )
    ):
        return None

    logger.info(
        "Hybrid hapto primary organometal module builder succeeded for %s%d with %d correlated fragment(s)%s",
        mol.GetAtomWithIdx(metal_idx).GetSymbol(),
        metal_idx,
        len(module.correlated_fragment_indices),
        " [preview]" if preview_only else "",
    )
    return scaffold_mol


def _place_secondary_metals_in_hapto_fragments(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> set:
    """Place non-hapto metals from their donor clouds in hapto-coupled systems."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return set()
    try:
        import numpy as np
    except ImportError:
        return set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    placed_secondary: set = set()

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        metal_sym = atom.GetSymbol()
        if metal_sym not in _METAL_SET or metal_idx in hapto_metals:
            continue

        donor_indices = [
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() not in _METAL_SET and nbr.GetAtomicNum() > 1
        ]
        if len(donor_indices) < 2:
            continue

        donor_positions = []
        donor_symbols = []
        donor_target_lengths = []
        donor_fit_weights = []
        skip = False
        for donor_idx in donor_indices:
            pos = conf.GetAtomPosition(donor_idx)
            arr = np.array([pos.x, pos.y, pos.z], dtype=float)
            if not np.all(np.isfinite(arr)):
                skip = True
                break
            donor_positions.append(arr)
            donor_symbols.append(mol.GetAtomWithIdx(donor_idx).GetSymbol())
            donor_target_lengths.append(_secondary_donor_target_length(mol, metal_sym, donor_idx))
            donor_fit_weights.append(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
        if skip:
            continue

        fit = _fit_secondary_metal_position_from_donors(
            donor_positions,
            donor_symbols,
            metal_sym,
            donor_target_lengths=donor_target_lengths,
            donor_fit_weights=donor_fit_weights,
        )
        if fit is None:
            continue
        metal_pos, geom_code = fit
        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )
        placed_secondary.add(metal_idx)
        logger.info(
            "Hybrid hapto placed secondary metal %s%d via %s fit to %d donor(s)",
            metal_sym,
            metal_idx,
            geom_code,
            len(donor_indices),
        )

    return placed_secondary


def _hybrid_bond_target_length(bond) -> float:
    """Approximate covalent bond target length for hybrid hapto relaxation."""
    begin_atom = bond.GetBeginAtom()
    end_atom = bond.GetEndAtom()
    begin_sym = begin_atom.GetSymbol()
    end_sym = end_atom.GetSymbol()

    if begin_sym in _METAL_SET or end_sym in _METAL_SET:
        metal_sym = begin_sym if begin_sym in _METAL_SET else end_sym
        donor_sym = end_sym if begin_sym in _METAL_SET else begin_sym
        return float(_get_ml_bond_length(metal_sym, donor_sym))

    pair = frozenset([begin_sym, end_sym])
    bond_type = bond.GetBondType()

    if 'H' in pair:
        other = end_sym if begin_sym == 'H' else begin_sym
        return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48, 'B': 1.19}.get(other, 1.08)

    if bond.GetIsAromatic() or bond_type == Chem.BondType.AROMATIC:
        aromatic_map = {
            frozenset(['C', 'C']): 1.39,
            frozenset(['C', 'N']): 1.35,
            frozenset(['C', 'O']): 1.36,
        }
        if pair in aromatic_map:
            return aromatic_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.12, 0.95)

    if bond_type == Chem.BondType.TRIPLE:
        triple_map = {
            frozenset(['C', 'C']): 1.20,
            frozenset(['C', 'N']): 1.16,
            frozenset(['N', 'N']): 1.10,
        }
        if pair in triple_map:
            return triple_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.22, 0.90)

    if bond_type == Chem.BondType.DOUBLE:
        double_map = {
            frozenset(['C', 'C']): 1.34,
            frozenset(['C', 'N']): 1.29,
            frozenset(['C', 'O']): 1.23,
            frozenset(['N', 'O']): 1.22,
            frozenset(['N', 'N']): 1.25,
        }
        if pair in double_map:
            return double_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.12, 0.95)

    single_map = {
        frozenset(['Si', 'C']): 1.87,
        frozenset(['C', 'C']): 1.50,
        frozenset(['C', 'N']): 1.47,
        frozenset(['C', 'O']): 1.43,
        frozenset(['C', 'Si']): 1.87,
        frozenset(['C', 'S']): 1.82,
        frozenset(['C', 'P']): 1.84,
        frozenset(['N', 'N']): 1.45,
        frozenset(['N', 'O']): 1.40,
        frozenset(['Si', 'Si']): 2.34,
    }
    if pair in single_map:
        return single_map[pair]
    return _COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76)
