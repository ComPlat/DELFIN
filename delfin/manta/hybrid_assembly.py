"""Assembly and refinement of hybrid hapto complexes, donor pi-coplanarity, sequential multi-metal hapto building, hapto previews and final clash resolution in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _mol_to_xyz,
)
from delfin.manta.converter_flags import (
    _class_conditional_flag,
)
from delfin.manta.hapto_scaffold import (
    _build_hapto_scaffold,
    _correct_hapto_geometry,
    _propagate_non_hapto_atoms,
)
from delfin.manta.hybrid_fragments import (
    _HybridHaptoDecomposition,
    _align_hybrid_fragment_onto_scaffold,
    _decompose_hapto_complex,
    _detect_primary_organometal_module,
    _embed_hybrid_fragment,
    _hapto_primary_donor_indices,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    Point3D,
    RDKIT_AVAILABLE,
    _METAL_SET,
    _get_ml_bond_length,
    _target_mc_dist,
)
from delfin.manta.secondary_metal_modules import (
    _align_hybrid_fragment_to_targets,
    _assemble_secondary_metal_coordination_modules,
    _build_primary_organometal_hapto_module,
    _enumerate_secondary_metal_position_fits,
    _hybrid_bond_target_length,
    _place_secondary_metals_in_hapto_fragments,
)

logger = get_logger("delfin.smiles_converter")


def _append_hapto_preview_xyz(
    mol,
    preview_store: Optional[List[Tuple[str, str]]],
    *,
    label: str,
) -> None:
    """Append the current XYZ of a preview candidate if available."""
    if preview_store is None or mol is None:
        return
    try:
        preview_xyz = _mol_to_xyz(mol)
    except Exception:
        return
    preview_entry = (preview_xyz, label)
    if preview_entry not in preview_store:
        preview_store.append(preview_entry)


def _store_hapto_preview_candidate(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    preview_store: Optional[List[Tuple[str, str]]],
    *,
    label: str,
) -> None:
    """Store a preview XYZ for a hapto scaffold candidate if available."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups or preview_store is None:
        return
    try:
        preview_mol = Chem.Mol(mol)
        if not _build_hapto_scaffold(preview_mol, hapto_groups):
            return
        try:
            _correct_hapto_geometry(preview_mol, 0, hapto_groups)
        except Exception:
            pass
        try:
            _propagate_non_hapto_atoms(
                preview_mol,
                0,
                hapto_groups,
                extra_fixed_indices=_hapto_primary_donor_indices(preview_mol, hapto_groups),
            )
        except Exception:
            pass
        _append_hapto_preview_xyz(
            preview_mol,
            preview_store,
            label=label,
        )
    except Exception:
        pass


def _refine_hybrid_hapto_complex(
    mol,
    decomposition: _HybridHaptoDecomposition,
    hapto_groups: List[Tuple[int, List[int]]],
    trusted_atom_indices: Optional[set] = None,
    secondary_metal_indices: Optional[set] = None,
    apply_final_hapto_correction: bool = True,
) -> bool:
    """Local hapto-only relaxation around a fixed coordination scaffold."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metals = set(decomposition.metal_indices)
    if secondary_metal_indices:
        metals = metals | set(secondary_metal_indices)
    hapto_atoms = {idx for _metal_idx, grp in hapto_groups for idx in grp}
    primary_hapto_donors = _hapto_primary_donor_indices(mol, hapto_groups)
    fixed = metals | hapto_atoms | primary_hapto_donors
    if trusted_atom_indices:
        fixed = fixed | set(trusted_atom_indices)
    movable = {idx for idx in range(mol.GetNumAtoms()) if idx not in fixed}
    if not movable:
        return True

    def _gp(idx):
        pos = conf.GetAtomPosition(idx)
        return np.array([pos.x, pos.y, pos.z], dtype=float)

    def _sp(idx, arr):
        conf.SetAtomPosition(idx, Point3D(float(arr[0]), float(arr[1]), float(arr[2])))

    bond_targets = []
    bonded_pairs: set = set()
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bonded_pairs.add((begin_idx, end_idx))
        bonded_pairs.add((end_idx, begin_idx))
        if begin_idx not in movable and end_idx not in movable:
            continue
        bond_targets.append((begin_idx, end_idx, _hybrid_bond_target_length(bond)))

    rng = np.random.default_rng(42)
    for _pass in range(80):
        displacements: Dict[int, np.ndarray] = {idx: np.zeros(3, dtype=float) for idx in movable}
        max_err = 0.0
        for begin_idx, end_idx, target in bond_targets:
            begin_pos = _gp(begin_idx)
            end_pos = _gp(end_idx)
            diff = end_pos - begin_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target
            if abs(delta) < 0.01:
                continue
            max_err = max(max_err, abs(delta))
            strength = 0.35 if ((begin_idx in fixed) ^ (end_idx in fixed)) else 0.22
            move = strength * delta * unit
            if begin_idx in movable and end_idx in movable:
                displacements[begin_idx] = displacements[begin_idx] + 0.5 * move
                displacements[end_idx] = displacements[end_idx] - 0.5 * move
            elif begin_idx in movable:
                displacements[begin_idx] = displacements[begin_idx] + move
            elif end_idx in movable:
                displacements[end_idx] = displacements[end_idx] - move

        if max_err < 0.04:
            break
        for atom_idx, disp in displacements.items():
            if float(np.linalg.norm(disp)) < 1e-6:
                continue
            _sp(atom_idx, _gp(atom_idx) + disp)

    for _pass in range(24):
        any_push = False
        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                if (i, j) in bonded_pairs:
                    continue
                if i not in movable and j not in movable:
                    continue
                pi = _gp(i)
                pj = _gp(j)
                dist = float(np.linalg.norm(pi - pj))
                sym_i = mol.GetAtomWithIdx(i).GetSymbol()
                sym_j = mol.GetAtomWithIdx(j).GetSymbol()

                if sym_i in _METAL_SET or sym_j in _METAL_SET:
                    min_dist = 2.0
                elif sym_i == 'H' and sym_j == 'H':
                    min_dist = 1.5
                elif sym_i == 'H' or sym_j == 'H':
                    min_dist = 1.0
                else:
                    min_dist = 1.18

                if dist >= min_dist:
                    continue
                if dist < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pj - pi) / dist

                gap = min_dist - dist
                if i in movable and j in movable:
                    _sp(i, pi - 0.5 * gap * direction)
                    _sp(j, pj + 0.5 * gap * direction)
                elif i in movable:
                    _sp(i, pi - gap * direction)
                elif j in movable:
                    _sp(j, pj + gap * direction)
                any_push = True
        if not any_push:
            break

    try:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        if ff is not None:
            for atom_idx in sorted(fixed):
                ff.AddFixedPoint(atom_idx)
            ff.Minimize(maxIts=250)
    except Exception:
        pass

    if apply_final_hapto_correction:
        try:
            _correct_hapto_geometry(mol, 0, hapto_groups)
        except Exception:
            pass
    return True


def _enforce_donor_pi_coplanarity(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
    max_rotation_deg: float = 75.0,
) -> bool:
    """Rigidly rotate planar donor fragments so the metal lies in the pi-plane.

    For each metal, identifies planar (aromatic / sp2-conjugated) fragments
    that bear donor atoms, then rotates each fragment as a rigid body so that
    the metal atom ends up in (or very near) the fragment's best-fit plane.

    * Bidentate+ donors  -> rotation around the donor-donor axis.
    * Monodentate donor   -> rotation around the single donor atom.

    Donor-metal bond lengths are preserved exactly because the rotation axis
    passes through the donor atom(s).
    """
    if not RDKIT_AVAILABLE or mol is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    hapto_atom_set: set = set()
    hapto_metal_set: set = set()
    if hapto_groups:
        for _mi, grp in hapto_groups:
            hapto_atom_set.update(grp)
            hapto_metal_set.add(_mi)

    max_rot = np.radians(max_rotation_deg)

    def _gp(idx):
        p = conf.GetAtomPosition(idx)
        return np.array([p.x, p.y, p.z], dtype=float)

    def _sp(idx, arr):
        conf.SetAtomPosition(
            idx, Point3D(float(arr[0]), float(arr[1]), float(arr[2])))

    def _rod(v, k, angle):
        c, s = np.cos(angle), np.sin(angle)
        return v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1.0 - c)

    def _is_planar_atom(a):
        if a.GetIsAromatic():
            return True
        try:
            hyb = a.GetHybridization()
        except Exception:
            return False
        if hyb in {Chem.rdchem.HybridizationType.SP,
                    Chem.rdchem.HybridizationType.SP2}:
            return True
        for bond in a.GetBonds():
            if bond.GetBondTypeAsDouble() >= 1.5 or bond.GetIsConjugated():
                return True
        return False

    changed = False

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        if atom.GetSymbol() not in _METAL_SET:
            continue
        # Skip hapto metals: piano-stool geometry should not be coplanarised
        if metal_idx in hapto_metal_set:
            continue

        # Non-hapto, non-metal, heavy-atom neighbours = candidate donors
        donors = []
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in hapto_atom_set or nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                continue
            donors.append(ni)
        if not donors:
            continue

        m_pos = _gp(metal_idx)
        processed: set = set()

        for start_donor in donors:
            if start_donor in processed:
                continue

            if not _is_planar_atom(mol.GetAtomWithIdx(start_donor)):
                processed.add(start_donor)
                continue

            # --- BFS: find connected planar core ---
            planar_core: set = {start_donor}
            bfs_q = [start_donor]
            while bfs_q:
                cur = bfs_q.pop(0)
                for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in planar_core or ni in hapto_atom_set:
                        continue
                    if nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                        continue
                    if _is_planar_atom(nbr):
                        planar_core.add(ni)
                        bfs_q.append(ni)

            if len(planar_core) < 3:
                processed.add(start_donor)
                continue

            frag_donors = [d for d in donors if d in planar_core]
            processed.update(frag_donors)
            frag_donor_set = set(frag_donors)

            # --- Best-fit plane (SVD) ---
            core_list = sorted(planar_core)
            pts = np.array([_gp(i) for i in core_list])
            centroid = pts.mean(axis=0)
            q = pts - centroid
            try:
                _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1].copy()
            n_len = float(np.linalg.norm(normal))
            if n_len < 1e-12:
                continue
            normal /= n_len

            # --- Rigid body: planar core + direct substituents + their H ---
            rigid_body: set = set(planar_core)
            for ci in list(planar_core):
                for nbr in mol.GetAtomWithIdx(ci).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in rigid_body or ni in hapto_atom_set or nbr.GetSymbol() in _METAL_SET:
                        continue
                    rigid_body.add(ni)
                    if nbr.GetAtomicNum() > 1:
                        for nbr2 in nbr.GetNeighbors():
                            ni2 = nbr2.GetIdx()
                            if ni2 not in rigid_body and nbr2.GetAtomicNum() == 1:
                                rigid_body.add(ni2)

            # ==========================================================
            # Bidentate+: rotate around the donor-donor axis
            # ==========================================================
            if len(frag_donors) >= 2:
                d1, d2 = frag_donors[0], frag_donors[1]
                d1_pos, d2_pos = _gp(d1), _gp(d2)
                axis = d2_pos - d1_pos
                axis_len = float(np.linalg.norm(axis))
                if axis_len < 1e-8:
                    continue
                u_ax = axis / axis_len

                # Gram-Schmidt: make normal exactly perpendicular to axis
                normal = normal - float(np.dot(normal, u_ax)) * u_ax
                n_len = float(np.linalg.norm(normal))
                if n_len < 1e-8:
                    continue
                normal /= n_len

                # a*cos(t) + b*sin(t) = 0  =>  t = atan2(-a, b)
                v_M = m_pos - d1_pos
                a_coeff = float(np.dot(v_M, normal))
                b_dir = np.cross(u_ax, normal)
                b_coeff = float(np.dot(v_M, b_dir))

                if abs(a_coeff) < 0.05:
                    continue

                theta = float(np.arctan2(-a_coeff, b_coeff))

                # Outward-side check: metal should face donors, not ring body
                body_atoms = [i for i in planar_core if i not in frag_donor_set]
                if body_atoms:
                    body_ctr = np.mean([_gp(i) for i in body_atoms], axis=0)
                    body_rot = d1_pos + _rod(body_ctr - d1_pos, u_ax, theta)
                    d_ctr = (d1_pos + d2_pos) / 2.0
                    outward = d_ctr - body_rot
                    if float(np.dot(m_pos - d_ctr, outward)) < 0:
                        theta += np.pi

                if abs(theta) > max_rot:
                    theta = float(np.sign(theta)) * max_rot

                for ai in rigid_body:
                    p = _gp(ai) - d1_pos
                    _sp(ai, d1_pos + _rod(p, u_ax, theta))
                changed = True

            # ==========================================================
            # Monodentate: rotate around the single donor
            # ==========================================================
            else:
                d_idx = frag_donors[0]
                d_pos = _gp(d_idx)

                v = m_pos - d_pos
                v_len = float(np.linalg.norm(v))
                if v_len < 1e-8:
                    continue

                h = float(np.dot(v, normal))
                if abs(h) < 0.05:
                    continue

                # Desired normal: component of normal perpendicular to v
                n_proj = normal - (float(np.dot(normal, v)) / (v_len ** 2)) * v
                n_proj_len = float(np.linalg.norm(n_proj))
                if n_proj_len < 1e-8:
                    continue
                desired = n_proj / n_proj_len

                rot_axis = np.cross(normal, desired)
                ra_len = float(np.linalg.norm(rot_axis))
                if ra_len < 1e-10:
                    continue
                rot_axis /= ra_len
                cos_a = float(np.clip(np.dot(normal, desired), -1.0, 1.0))
                theta = float(np.arccos(cos_a))

                if theta > max_rot:
                    theta = max_rot

                for ai in rigid_body:
                    p = _gp(ai) - d_pos
                    _sp(ai, d_pos + _rod(p, rot_axis, theta))
                changed = True

    return changed


def _build_multimetal_hapto_sequential(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    decomposition: '_HybridHaptoDecomposition',
    secondary_variant_plan: Optional[Dict[int, int]] = None,
):
    """Sequential multi-metal builder with rigid-body clash-free placement.

    Strategy:
      1. Build hapto scaffold (ferrocene / Cp rings + hapto metal)
      2. Embed & align bridge fragments (Procrustes + UFF relax)
      3. Place secondary metal(s) from donor positions
      4. Place remaining fragments as rigid bodies, rotating around
         metal-donor axis to avoid clashes with already-placed atoms
      5. Whole-fragment M-L correction (rigid shift, no bond breaking)
      6. Final topology-preserving clash resolution
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups or decomposition is None:
        return None
    try:
        import numpy as np
        from rdkit.Geometry import Point3D
    except ImportError:
        return None

    # ---- helpers ----
    def _gp(conf, idx):
        p = conf.GetAtomPosition(idx)
        return np.array([p.x, p.y, p.z], dtype=float)

    def _sp(conf, idx, pos):
        conf.SetAtomPosition(idx, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    def _check_clash(conf, new_positions, placed_set, bonded_pairs, threshold=1.5):
        """Check if new_positions clash with placed atoms.
        Returns worst clash distance (0 = no clash, smaller = worse)."""
        worst = 999.0
        for new_idx, new_pos in new_positions.items():
            for placed_idx in placed_set:
                if placed_idx == new_idx:
                    continue
                pair = (min(new_idx, placed_idx), max(new_idx, placed_idx))
                if pair in bonded_pairs:
                    continue
                placed_pos = _gp(conf, placed_idx)
                d = float(np.linalg.norm(new_pos - placed_pos))
                if d < worst:
                    worst = d
        return worst

    def _rodrigues_rotate(points, axis, angle):
        """Rotate points around axis by angle (radians) using Rodrigues."""
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        return (points * cos_a
                + np.cross(axis, points) * sin_a
                + axis[np.newaxis, :] * (points @ axis[:, np.newaxis]) * (1 - cos_a))

    # ---- classify atoms ----
    hapto_metal_set = {mi for mi, _g in hapto_groups}
    all_hapto_carbons: set = set()
    hapto_carbons_by_metal: Dict[int, set] = {}
    for mi, grp in hapto_groups:
        all_hapto_carbons.update(grp)
        hapto_carbons_by_metal.setdefault(mi, set()).update(grp)

    all_metals: set = set()
    non_hapto_metals: set = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METAL_SET:
            all_metals.add(atom.GetIdx())
            if atom.GetIdx() not in hapto_metal_set:
                non_hapto_metals.add(atom.GetIdx())

    if not non_hapto_metals and len(hapto_metal_set) < 2:
        return None  # not a multi-metal system

    # ==== STEP 1: Build hapto scaffold (Fe + Cp rings) ====
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

    # Only trust hapto metals initially — secondary metals are NOT placed yet
    trusted: set = set(hapto_metal_set) | set(all_hapto_carbons)
    logger.info(
        "Sequential multi-metal: scaffold built, %d hapto atoms + %d hapto metals trusted",
        len(all_hapto_carbons), len(hapto_metal_set),
    )

    # Helper: Rodrigues rotation (used in Steps 2 and 4)
    def _rotate_fragment_around_axis(orig_positions, pivot, axis, angle):
        """Rotate a dict {idx: pos} around pivot+axis by angle radians."""
        ax = axis / max(float(np.linalg.norm(axis)), 1e-12)
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        result = {}
        for idx, pos in orig_positions.items():
            v = pos - pivot
            v_rot = (v * cos_a
                     + np.cross(ax, v) * sin_a
                     + ax * float(np.dot(ax, v)) * (1 - cos_a))
            result[idx] = pivot + v_rot
        return result

    # ==== STEP 2: Embed & align all ligand fragments ====
    # Strategy: ETKDG embed (gives correct ring geometry) → Procrustes align
    # to scaffold anchors (with anchor snap) → UFF relaxation with anchors
    # fixed to repair any bond stretching caused by the snap.
    for fragment in decomposition.fragments:
        # Pure hapto-ring fragments: trust scaffold geometry already
        if all(ai in all_hapto_carbons for ai in fragment.atom_indices):
            trusted.update(fragment.atom_indices)
            continue

        # Check if fragment has anchors in trusted set
        has_anchor = any(ai in trusted for ai in fragment.atom_indices)
        if not has_anchor:
            continue  # defer to Step 4

        # Embed fragment with ETKDG (correct ring geometry)
        embedded = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded is None:
            continue

        # Procrustes align to scaffold (with anchor snap for junction correctness)
        if not _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
            continue

        # UFF relaxation: fix anchor atoms, let non-anchors relax to fix
        # any bond stretching caused by the anchor snap.
        try:
            # Build a temporary molecule from the fragment with aligned coords
            frag_mol = Chem.Mol(fragment.fragment_mol)
            frag_mol.RemoveAllConformers()
            frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
            for frag_idx, orig_idx in fragment.fragment_to_original.items():
                p = conf.GetAtomPosition(orig_idx)
                frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
            frag_mol.AddConformer(frag_conf_obj, assignId=True)

            ff = AllChem.UFFGetMoleculeForceField(frag_mol)
            if ff is not None:
                # Fix anchor atoms (Cp carbons that must stay at scaffold positions)
                n_fixed = 0
                for orig_idx in fragment.anchor_atom_indices:
                    frag_idx = fragment.original_to_fragment.get(orig_idx)
                    if frag_idx is not None and orig_idx in trusted:
                        ff.AddFixedPoint(frag_idx)
                        n_fixed += 1
                e_before = ff.CalcEnergy()
                ff.Minimize(maxIts=500)
                e_after = ff.CalcEnergy()
                logger.debug(
                    "UFF relax: %d fixed, E %.1f -> %.1f",
                    n_fixed, e_before, e_after,
                )

                # Copy ALL relaxed coordinates back to scaffold
                # (anchors were fixed so their positions didn't change)
                relaxed_conf = frag_mol.GetConformer()
                for frag_idx, orig_idx in fragment.fragment_to_original.items():
                    p = relaxed_conf.GetAtomPosition(frag_idx)
                    conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
        except Exception as e:
            logger.debug("UFF relaxation failed for fragment: %s", e)

        # ---- Clash-aware torsion for bridge fragments ----
        # If this fragment contains hapto carbons, the non-hapto substituents
        # (e.g., pyridine ring) may clash with the OTHER Cp ring of the
        # metallocene.  Fix by rotating the non-hapto substructure around
        # the bond connecting it to the hapto ring.
        frag_hapto = [ai for ai in fragment.atom_indices if ai in all_hapto_carbons]
        frag_non_hapto_heavy = [ai for ai in fragment.atom_indices
                                if ai not in all_hapto_carbons
                                and ai not in all_metals
                                and mol.GetAtomWithIdx(ai).GetAtomicNum() > 1]
        if frag_hapto and frag_non_hapto_heavy:
            # Find the OTHER Cp ring's atoms (for clash checking)
            other_cp_atoms: set = set()
            for hm in hapto_metal_set:
                for nbr in mol.GetAtomWithIdx(hm).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in all_hapto_carbons and ni not in frag_hapto:
                        other_cp_atoms.add(ni)

            # Find the bond connecting hapto → non-hapto HEAVY atom
            bridge_bond = None
            for hi in frag_hapto:
                for nbr in mol.GetAtomWithIdx(hi).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in frag_non_hapto_heavy:
                        bridge_bond = (hi, ni)
                        break
                if bridge_bond is not None:
                    break

            if bridge_bond is not None and other_cp_atoms:
                pivot_idx, next_idx = bridge_bond
                pivot_pos = _gp(conf, pivot_idx)
                next_pos = _gp(conf, next_idx)
                rot_axis = next_pos - pivot_pos
                rot_axis_len = float(np.linalg.norm(rot_axis))

                if rot_axis_len > 1e-6:
                    # Check all clash-relevant atoms (other Cp + all trusted
                    # not in this fragment)
                    clash_targets = (other_cp_atoms | trusted) - set(fragment.atom_indices)
                    # Collect atoms to rotate: all non-hapto atoms (heavy + H)
                    frag_non_hapto_all = [ai for ai in fragment.atom_indices
                                          if ai not in all_hapto_carbons
                                          and ai not in all_metals]
                    rot_atoms = list(frag_non_hapto_all)
                    orig_positions = {ai: _gp(conf, ai) for ai in rot_atoms}

                    def _score_rotation(angle):
                        if abs(angle) < 1e-9:
                            positions = orig_positions
                        else:
                            positions = _rotate_fragment_around_axis(
                                orig_positions, pivot_pos, rot_axis, angle)
                        min_d = 999.0
                        for ai, pos in positions.items():
                            for ti in clash_targets:
                                d = float(np.linalg.norm(pos - _gp(conf, ti)))
                                if d < min_d:
                                    min_d = d
                        return min_d

                    best_angle = 0.0
                    best_min_d = _score_rotation(0.0)
                    for step in range(1, 36):
                        angle = step * (2 * np.pi / 36)
                        min_d = _score_rotation(angle)
                        if min_d > best_min_d:
                            best_min_d = min_d
                            best_angle = angle

                    if best_angle != 0.0:
                        rotated = _rotate_fragment_around_axis(
                            orig_positions, pivot_pos, rot_axis, best_angle)
                        for ai, pos in rotated.items():
                            _sp(conf, ai, pos)
                        logger.info(
                            "Step 2 torsion: rotated %d non-hapto atoms by %.0f° "
                            "around %s%d-%s%d, min_d=%.2f",
                            len(rot_atoms), np.degrees(best_angle),
                            mol.GetAtomWithIdx(pivot_idx).GetSymbol(), pivot_idx,
                            mol.GetAtomWithIdx(next_idx).GetSymbol(), next_idx,
                            best_min_d,
                        )

        trusted.update(fragment.atom_indices)
        logger.info(
            "Sequential: embed+align+relax fragment (%d atoms, bridge=%s)",
            len(fragment.atom_indices), fragment.use_scaffold_only,
        )

    # Mark all atoms reachable from trusted set (excluding secondary metals)
    from collections import deque
    bfs_visited: set = set(trusted)
    bfs_q: deque = deque(trusted)
    while bfs_q:
        curr = bfs_q.popleft()
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in bfs_visited:
                continue
            if ni in non_hapto_metals:
                continue  # don't cross into secondary metals yet
            bfs_visited.add(ni)
            bfs_q.append(ni)
    trusted = bfs_visited
    logger.info("Sequential Step 2: %d atoms trusted after fragment alignment", len(trusted))

    # ==== STEP 3: Place secondary (non-hapto) metals ====
    for metal_idx in sorted(non_hapto_metals):
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        variant_rank = max(0, int((secondary_variant_plan or {}).get(metal_idx, 0)))
        # Collect donors that are already placed (in trusted set)
        donor_indices = []
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in all_metals:
                continue
            if nbr.GetAtomicNum() <= 1:
                continue
            donor_indices.append(ni)

        placed_donors = [d for d in donor_indices if d in trusted]
        unplaced_donors = [d for d in donor_indices if d not in trusted]

        if len(placed_donors) < 1:
            logger.debug(
                "Sequential: skipping %s%d (no placed donors yet)", metal_sym, metal_idx,
            )
            continue

        # Get positions and symbols of placed donors
        donor_positions = []
        donor_symbols = []
        target_lengths = []
        for d in placed_donors:
            donor_positions.append(_gp(conf, d))
            d_sym = mol.GetAtomWithIdx(d).GetSymbol()
            donor_symbols.append(d_sym)
            target_lengths.append(_get_ml_bond_length(metal_sym, d_sym))

        if len(placed_donors) >= 2:
            # For bimetallic complexes, the LIN/centroid fit from 2 donors
            # often places the metal BETWEEN close donors (inside the
            # ferrocene cage). Strategy:
            #   1. Try SVD fit first
            #   2. Check inter-metal distance; if too close to hapto metals,
            #      use outward placement from hapto center through donor midpoint
            hapto_center = np.mean([_gp(conf, mi) for mi in hapto_metal_set], axis=0)
            avg_len = np.mean(target_lengths)

            fit_candidates = _enumerate_secondary_metal_position_fits(
                donor_positions, donor_symbols, metal_sym,
                donor_target_lengths=target_lengths,
                max_candidates=max(variant_rank + 1, 4),
            )
            fit_result = (
                fit_candidates[min(variant_rank, len(fit_candidates) - 1)]
                if fit_candidates else None
            )
            geom_code = "SVD"
            if fit_result is not None:
                metal_pos, geom_code, _fit_score = fit_result
            else:
                metal_pos = np.mean(donor_positions, axis=0)

            # Check if fitted position is too close to any hapto metal
            too_close = False
            for hm in hapto_metal_set:
                hm_pos = _gp(conf, hm)
                if float(np.linalg.norm(metal_pos - hm_pos)) < 2.5:
                    too_close = True
                    break

            if too_close:
                # Intersection-circle search: for 2+ donors, the ideal
                # metal position lies on the intersection circle of the
                # donor spheres. Search on this circle first, then fall
                # back to multi-sphere search.
                best_score = -1e9
                best_cand = metal_pos
                geom_code = "SEARCH"

                # Build candidate list from intersection circle + sphere
                candidates: list = []

                if len(donor_positions) >= 2:
                    # Compute intersection circle of first two donors
                    c1 = np.array(donor_positions[0])
                    c2 = np.array(donor_positions[1])
                    r1 = target_lengths[0]
                    r2 = target_lengths[1]
                    axis_vec = c2 - c1
                    d_donors = float(np.linalg.norm(axis_vec))
                    if d_donors > 1e-6 and d_donors < r1 + r2:
                        axis_unit = axis_vec / d_donors
                        # Distance from c1 to intersection plane
                        h = (r1 ** 2 - r2 ** 2 + d_donors ** 2) / (2 * d_donors)
                        circle_r_sq = r1 ** 2 - h ** 2
                        if circle_r_sq > 0:
                            circle_r = np.sqrt(circle_r_sq)
                            circle_center = c1 + h * axis_unit
                            # Build orthonormal basis in the plane
                            arb = np.array([1.0, 0.0, 0.0])
                            if abs(float(np.dot(axis_unit, arb))) > 0.9:
                                arb = np.array([0.0, 1.0, 0.0])
                            u = np.cross(axis_unit, arb)
                            u = u / float(np.linalg.norm(u))
                            v = np.cross(axis_unit, u)
                            # Sample circle at 5° resolution
                            for deg in range(0, 360, 5):
                                angle = np.radians(deg)
                                cand = (circle_center
                                        + circle_r * (np.cos(angle) * u + np.sin(angle) * v))
                                candidates.append(cand)

                # Also sample spheres around each donor (10° resolution)
                for center_pos, center_len in zip(donor_positions, target_lengths):
                    center = np.array(center_pos)
                    for theta_deg in range(0, 180, 15):
                        for phi_deg in range(0, 360, 15):
                            theta = np.radians(theta_deg)
                            phi = np.radians(phi_deg)
                            d_vec = np.array([
                                np.sin(theta) * np.cos(phi),
                                np.sin(theta) * np.sin(phi),
                                np.cos(theta),
                            ])
                            candidates.append(center + d_vec * center_len)

                scored_candidates: List[Tuple[float, float, object]] = []
                if fit_result is not None:
                    min_clash_fit = min(
                        (float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - _gp(conf, ti)))
                         for ti in trusted if ti not in all_metals),
                        default=999.0,
                    )
                    min_hm_fit = min(
                        float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - _gp(conf, hm)))
                        for hm in hapto_metal_set
                    )
                    ml_err_fit = sum(
                        (float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - np.array(dp2))) - tl2) ** 2
                        for dp2, tl2 in zip(donor_positions, target_lengths)
                    )
                    scored_candidates.append(
                        (
                            -ml_err_fit * 2.0 + min(min_clash_fit, 3.0) * 0.3 + min(min_hm_fit, 4.0) * 0.12,
                            min_clash_fit,
                            np.asarray(metal_pos, dtype=float),
                        )
                    )

                for clash_thresh in (1.8, 1.0):
                    for cand in candidates:
                        min_hm = min(
                            float(np.linalg.norm(cand - _gp(conf, hm)))
                            for hm in hapto_metal_set
                        )
                        if min_hm < 2.0:
                            continue
                        ml_err = sum(
                            (float(np.linalg.norm(cand - np.array(dp2))) - tl2) ** 2
                            for dp2, tl2 in zip(donor_positions, target_lengths)
                        )
                        min_clash = min(
                            (float(np.linalg.norm(cand - _gp(conf, ti)))
                             for ti in trusted if ti not in all_metals),
                            default=999.0,
                        )
                        if min_clash < clash_thresh:
                            continue
                        score = -ml_err * 2.0 + min(min_clash, 3.0) * 0.3
                        score += min(min_hm, 4.0) * 0.12
                        scored_candidates.append((score, min_clash, np.asarray(cand, dtype=float)))
                    if scored_candidates:
                        break

                if scored_candidates:
                    scored_candidates.sort(key=lambda item: (-item[0], -item[1]))
                    unique_candidates: List[Tuple[float, float, object]] = []
                    for score, min_clash, cand in scored_candidates:
                        if any(
                            float(np.linalg.norm(np.asarray(cand, dtype=float) - np.asarray(prev_cand, dtype=float))) < 0.28
                            for _score_prev, _clash_prev, prev_cand in unique_candidates
                        ):
                            continue
                        unique_candidates.append((score, min_clash, cand))
                    chosen_idx = min(variant_rank, len(unique_candidates) - 1)
                    best_score, best_min_clash, best_cand = unique_candidates[chosen_idx]
                    metal_pos = np.asarray(best_cand, dtype=float)
                    logger.info(
                        "SEARCH: selected ranked position %d/%d, score=%.2f, "
                        "min_clash=%.2f, ml_dists=%s",
                        chosen_idx + 1,
                        len(unique_candidates),
                        best_score,
                        best_min_clash,
                        [f"{float(np.linalg.norm(best_cand - np.array(dp))):.2f}"
                         for dp in donor_positions],
                    )

            _sp(conf, metal_idx, metal_pos)
            trusted.add(metal_idx)
            logger.info(
                "Sequential: placed %s%d via %s fit from %d donors, "
                "dist to hapto=%.2f",
                metal_sym, metal_idx, geom_code, len(placed_donors),
                float(np.linalg.norm(metal_pos - hapto_center)),
            )
        elif len(placed_donors) == 1:
            # Single placed donor: place metal along donor→outward direction
            d_pos = np.array(donor_positions[0])
            # Direction: away from hapto metal(s)
            hapto_center = np.mean([_gp(conf, mi) for mi in hapto_metal_set], axis=0)
            direction = d_pos - hapto_center
            norm = np.linalg.norm(direction)
            if norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
            else:
                direction = direction / norm
            _sp(conf, metal_idx, d_pos + direction * target_lengths[0])
            trusted.add(metal_idx)
            logger.info(
                "Sequential: placed %s%d from single donor %s%d at %.2f A",
                metal_sym, metal_idx,
                donor_symbols[0], placed_donors[0], target_lengths[0],
            )

        # Now place unplaced donors around the newly placed metal
        if metal_idx in trusted and unplaced_donors:
            metal_pos = _gp(conf, metal_idx)
            # Compute used directions (from metal to placed donors)
            used_dirs = []
            for d in placed_donors:
                vec = _gp(conf, d) - metal_pos
                n = np.linalg.norm(vec)
                if n > 1e-6:
                    used_dirs.append(vec / n)

            for ud in unplaced_donors:
                ud_sym = mol.GetAtomWithIdx(ud).GetSymbol()
                bl = _get_ml_bond_length(metal_sym, ud_sym)
                # Find direction maximally separated from ALL used directions.
                # For each candidate direction on a sphere, compute the
                # minimum angle to any used direction. Pick the candidate
                # with the largest minimum angle.
                if not used_dirs:
                    new_dir = np.array([0.0, 0.0, 1.0])
                else:
                    best_dir = np.array([0.0, 0.0, 1.0])
                    best_min_angle = -1.0
                    for t_deg in range(0, 180, 15):
                        for p_deg in range(0, 360, 15):
                            t = np.radians(t_deg)
                            p = np.radians(p_deg)
                            cand = np.array([
                                np.sin(t) * np.cos(p),
                                np.sin(t) * np.sin(p),
                                np.cos(t),
                            ])
                            # Minimum angle to any used direction
                            min_ang = min(
                                np.arccos(np.clip(float(np.dot(cand, ud)), -1, 1))
                                for ud in used_dirs
                            )
                            if min_ang > best_min_angle:
                                best_min_angle = min_ang
                                best_dir = cand
                    new_dir = best_dir

                _sp(conf, ud, metal_pos + new_dir * bl)
                used_dirs.append(new_dir)
                trusted.add(ud)

    # ==== Build bonded-pairs set for clash detection ====
    bonded_pairs: set = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((min(a, b), max(a, b)))
    # Also treat hapto-metal ↔ hapto-carbon as bonded
    for mi, grp in hapto_groups:
        for ci in grp:
            bonded_pairs.add((min(mi, ci), max(mi, ci)))
    # Also 1-3 pairs (atoms 2 bonds apart) — these are NOT clashes
    pairs_13: set = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for nbr_a in mol.GetAtomWithIdx(a).GetNeighbors():
            ni = nbr_a.GetIdx()
            if ni != b:
                pairs_13.add((min(ni, b), max(ni, b)))
        for nbr_b in mol.GetAtomWithIdx(b).GetNeighbors():
            ni = nbr_b.GetIdx()
            if ni != a:
                pairs_13.add((min(ni, a), max(ni, a)))
    bonded_pairs |= pairs_13

    def _min_dist_to_placed(positions_dict, placed_set, skip_set=frozenset()):
        """Return minimum non-bonded distance between new atoms and placed atoms."""
        worst = 999.0
        for new_idx, new_pos in positions_dict.items():
            for placed_idx in placed_set:
                if placed_idx in skip_set:
                    continue
                if placed_idx == new_idx:
                    continue
                pair = (min(new_idx, placed_idx), max(new_idx, placed_idx))
                if pair in bonded_pairs:
                    continue
                placed_pos = _gp(conf, placed_idx)
                d = float(np.linalg.norm(new_pos - placed_pos))
                if d < worst:
                    worst = d
        return worst

    # ==== STEP 4: Embed & place remaining fragments WITH CLASH-AWARE ROTATION ====
    # For fragments connected to secondary metals whose scaffold positions are
    # NOT yet set, we cannot use _align_hybrid_fragment_onto_scaffold (it reads
    # uninitialised scaffold coordinates). Instead:
    #   1. Embed with ETKDG (correct internal geometry)
    #   2. Find donor atoms connecting to the secondary metal
    #   3. Translate fragment so donor centroid sits at correct M-L distance
    #   4. Rotate 360° around metal→donor axis, pick clash-free orientation
    for fragment in decomposition.fragments:
        if all(ai in trusted for ai in fragment.atom_indices):
            continue
        # Must have at least one anchor already placed
        has_anchor = any(ai in trusted for ai in fragment.atom_indices)
        if not has_anchor:
            # Check if this fragment has a metal neighbor that's trusted
            for ai in fragment.atom_indices:
                for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                    if nbr.GetIdx() in trusted and nbr.GetIdx() in all_metals:
                        has_anchor = True
                        break
                if has_anchor:
                    break
        if not has_anchor:
            continue

        embedded = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded is None:
            continue

        # Find donor atoms in this fragment that connect to a trusted metal
        donor_to_metal: Dict[int, int] = {}  # orig_donor_idx → metal_idx
        for ai in fragment.atom_indices:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in trusted and ni in all_metals and ni not in fragment.atom_indices:
                    donor_to_metal[ai] = ni

        # Also check anchor atoms already in trusted set
        anchors_in_trusted = [ai for ai in fragment.anchor_atom_indices if ai in trusted]

        # If we have anchor atoms in the scaffold (from Step 2), use standard alignment
        if anchors_in_trusted and not donor_to_metal:
            if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
                # UFF relaxation with fixed anchors
                try:
                    frag_mol = Chem.Mol(fragment.fragment_mol)
                    frag_mol.RemoveAllConformers()
                    frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
                    for frag_idx, orig_idx in fragment.fragment_to_original.items():
                        p = conf.GetAtomPosition(orig_idx)
                        frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
                    frag_mol.AddConformer(frag_conf_obj, assignId=True)
                    ff = AllChem.UFFGetMoleculeForceField(frag_mol)
                    if ff is not None:
                        for orig_idx in fragment.anchor_atom_indices:
                            frag_idx = fragment.original_to_fragment.get(orig_idx)
                            if frag_idx is not None and orig_idx in trusted:
                                ff.AddFixedPoint(frag_idx)
                        ff.Minimize(maxIts=500)
                        relaxed_conf = frag_mol.GetConformer()
                        for frag_idx, orig_idx in fragment.fragment_to_original.items():
                            if orig_idx in trusted and orig_idx in all_hapto_carbons:
                                continue
                            p = relaxed_conf.GetAtomPosition(frag_idx)
                            conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
                except Exception:
                    pass
                trusted.update(fragment.atom_indices)
                logger.info(
                    "Sequential Step 4: scaffold-aligned fragment (%d atoms)",
                    len(fragment.atom_indices),
                )
                continue

        # ---- Placement for metal-connected fragments ----
        # The donor atoms (e.g. O11, O15) were placed at M-L distance from the
        # metal in Step 3. Use those positions as alignment targets: Procrustes
        # align the ETKDG-embedded fragment so its donors match the Step 3
        # positions. This preserves internal fragment geometry while placing
        # donors at correct M-L distances.
        if not donor_to_metal:
            # Fallback: try standard alignment anyway
            if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
                trusted.update(fragment.atom_indices)
            continue

        donor_indices_list = list(donor_to_metal.keys())
        metal_idx_for_frag = list(donor_to_metal.values())[0]
        metal_pos = _gp(conf, metal_idx_for_frag)
        metal_sym = mol.GetAtomWithIdx(metal_idx_for_frag).GetSymbol()

        # Target positions for donor atoms = their current scaffold positions
        # (set in Step 3's unplaced-donor VSEPR placement)
        target_positions = {d: _gp(conf, d) for d in donor_indices_list}
        for d in donor_indices_list:
            tp = target_positions[d]
            td = float(np.linalg.norm(tp - metal_pos))
            logger.info("Step 4 target: %s%d at %.2f from %s%d, pos=%s",
                        mol.GetAtomWithIdx(d).GetSymbol(), d, td,
                        metal_sym, metal_idx_for_frag, tp)

        # Procrustes alignment of ETKDG embedding to donor targets
        aligned = _align_hybrid_fragment_to_targets(
            scaffold_mol, fragment, embedded,
            target_atom_indices=donor_indices_list,
            target_positions=target_positions,
            reference_center=metal_pos,
            metal_idx=metal_idx_for_frag,
        )
        if not aligned:
            # Fallback: try standard alignment
            logger.info("Step 4: Procrustes alignment FAILED, using scaffold alignment")
            _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded)

        # ---- Snap donors to their target positions ----
        # Procrustes alignment may leave donors at wrong M-L distances
        # because the fragment's donor-donor distance differs from the
        # chelate bite distance. Snap each donor atom to its target
        # position, then UFF-relax the body with donors fixed.
        for d_idx in donor_indices_list:
            tp = target_positions[d_idx]
            _sp(conf, d_idx, np.asarray(tp, dtype=float))
            d_dist = float(np.linalg.norm(np.asarray(tp) - metal_pos))
            logger.info("Step 4 snap: %s%d → %.2f from %s%d",
                        mol.GetAtomWithIdx(d_idx).GetSymbol(), d_idx,
                        d_dist, metal_sym, metal_idx_for_frag)

        # UFF relaxation: fix donor atoms at correct positions, relax
        # body to accommodate the snap (repairs stretched bonds)
        try:
            frag_mol = Chem.Mol(fragment.fragment_mol)
            frag_mol.RemoveAllConformers()
            frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
            for frag_idx, orig_idx in fragment.fragment_to_original.items():
                p = conf.GetAtomPosition(orig_idx)
                frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
            frag_mol.AddConformer(frag_conf_obj, assignId=True)
            ff = AllChem.UFFGetMoleculeForceField(frag_mol)
            if ff is not None:
                for d_idx in donor_indices_list:
                    frag_idx = fragment.original_to_fragment.get(d_idx)
                    if frag_idx is not None:
                        ff.AddFixedPoint(frag_idx)
                ff.Minimize(maxIts=500)
                relaxed_conf = frag_mol.GetConformer()
                for frag_idx, orig_idx in fragment.fragment_to_original.items():
                    p = relaxed_conf.GetAtomPosition(frag_idx)
                    conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
        except Exception:
            pass

        # ---- Clash-aware rotation around metal→donor axis ----
        # Rotate the ENTIRE fragment as a rigid body (including donors)
        # to avoid breaking donor-to-body bonds.
        all_frag_atoms = list(fragment.atom_indices)
        non_frag_trusted = trusted - set(all_frag_atoms)
        if all_frag_atoms and non_frag_trusted:
            # Rotation axis: metal → donor centroid
            donor_centroid = np.mean([_gp(conf, d) for d in donor_indices_list], axis=0)
            axis = donor_centroid - metal_pos
            if float(np.linalg.norm(axis)) < 1e-6:
                axis = np.array([0.0, 0.0, 1.0])
            pivot = metal_pos
            orig_positions = {ai: _gp(conf, ai) for ai in all_frag_atoms}

            best_angle = 0.0
            best_min_d = _min_dist_to_placed(orig_positions, non_frag_trusted)

            for step in range(1, 36):
                angle = step * (2 * np.pi / 36)
                rotated = _rotate_fragment_around_axis(orig_positions, pivot, axis, angle)
                min_d = _min_dist_to_placed(rotated, non_frag_trusted)
                if min_d > best_min_d:
                    best_min_d = min_d
                    best_angle = angle

            if best_angle != 0.0:
                rotated = _rotate_fragment_around_axis(orig_positions, pivot, axis, best_angle)
                for ai, pos in rotated.items():
                    _sp(conf, ai, pos)

            logger.info(
                "Sequential Step 4: placed fragment (%d atoms) on %s%d, "
                "rot=%.0f°, min_clash=%.2f",
                len(fragment.atom_indices), metal_sym, metal_idx_for_frag,
                np.degrees(best_angle), best_min_d,
            )

        trusted.update(fragment.atom_indices)

    # BFS-propagate any truly remaining atoms (H atoms, small substituents only)
    remaining = set(range(mol.GetNumAtoms())) - trusted
    if remaining:
        try:
            _propagate_non_hapto_atoms(
                scaffold_mol, 0, hapto_groups,
                extra_fixed_indices=trusted,
            )
        except Exception as e:
            logger.debug("Sequential Step 4 BFS fallback: %s", e)
    trusted = set(range(mol.GetNumAtoms()))

    # ==== STEP 5: Build rigid body cache (skip M-L re-fit to preserve layout) ====
    # The metal was placed in Step 3, and fragments in Step 4 — both with
    # clash-aware placement. Re-fitting the metal position now would risk
    # moving it into fragment atoms. M-L distances may not be perfect but
    # topology is preserved.
    stop_set = all_metals | all_hapto_carbons
    frag_body_cache: Dict[int, Set[int]] = {}

    for metal_idx in sorted(non_hapto_metals):
        donor_indices = [
            nbr.GetIdx() for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
            if nbr.GetIdx() not in all_metals and nbr.GetAtomicNum() > 1
        ]
        for d_idx in donor_indices:
            if d_idx in all_hapto_carbons or d_idx in frag_body_cache:
                continue
            body = {d_idx}
            bfs_q_rb = [d_idx]
            while bfs_q_rb:
                curr = bfs_q_rb.pop()
                for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in body or ni in stop_set:
                        continue
                    body.add(ni)
                    bfs_q_rb.append(ni)
            frag_body_cache[d_idx] = body

        logger.info(
            "Sequential Step 5: %s%d kept, dists=%s",
            mol.GetAtomWithIdx(metal_idx).GetSymbol(), metal_idx,
            [f"{float(np.linalg.norm(_gp(conf, d) - _gp(conf, metal_idx))):.2f}" for d in donor_indices],
        )

    # ==== STEP 6: Final fixes ====
    conf = scaffold_mol.GetConformer(0)
    n_atoms = scaffold_mol.GetNumAtoms()

    # 6a: Clash-resolution rotations DISABLED — they compound and create
    # new clashes worse than the original. Step 4 clash-aware rotation
    # handles initial placement; post-hoc rotation just moves problems around.

    # 6b: Fix H atoms with bad parent distances or clashes
    for _h_pass in range(3):
        any_fixed = False
        for i in range(n_atoms):
            if scaffold_mol.GetAtomWithIdx(i).GetAtomicNum() != 1:
                continue
            nbrs = [n.GetIdx() for n in scaffold_mol.GetAtomWithIdx(i).GetNeighbors()]
            if not nbrs:
                continue
            parent_idx = nbrs[0]
            parent_pos = _gp(conf, parent_idx)
            h_pos = _gp(conf, i)
            parent_dist = float(np.linalg.norm(h_pos - parent_pos))

            needs_fix = parent_dist < 0.8 or parent_dist > 1.5
            if not needs_fix:
                for j in range(n_atoms):
                    if j == i or j == parent_idx:
                        continue
                    d = float(np.linalg.norm(h_pos - _gp(conf, j)))
                    if d < (1.0 if scaffold_mol.GetAtomWithIdx(j).GetAtomicNum() > 1 else 1.2):
                        needs_fix = True
                        break
            if not needs_fix:
                continue
            any_fixed = True

            # Place H opposite to parent's other neighbors
            parent_nbrs = [n.GetIdx() for n in scaffold_mol.GetAtomWithIdx(parent_idx).GetNeighbors()
                           if n.GetIdx() != i]
            if parent_nbrs:
                avg_dir = np.zeros(3)
                for pn in parent_nbrs:
                    v = _gp(conf, pn) - parent_pos
                    vn = float(np.linalg.norm(v))
                    if vn > 1e-6:
                        avg_dir += v / vn
                avg_len = float(np.linalg.norm(avg_dir))
                h_dir = -avg_dir / avg_len if avg_len > 1e-6 else np.array([0.0, 0.0, 1.0])
            else:
                h_dir = np.array([0.0, 0.0, 1.0])

            # Try multiple rotations, pick best
            best_pos = parent_pos + h_dir * 1.09
            best_min_d = 0.0
            for j2 in range(n_atoms):
                if j2 == i or j2 == parent_idx:
                    continue
                d2 = float(np.linalg.norm(best_pos - _gp(conf, j2)))
                if best_min_d == 0.0 or d2 < best_min_d:
                    best_min_d = d2

            if best_min_d < 1.0:
                perp = np.cross(h_dir, np.array([1.0, 0.0, 0.0]))
                if float(np.linalg.norm(perp)) < 1e-6:
                    perp = np.cross(h_dir, np.array([0.0, 1.0, 0.0]))
                perp = perp / float(np.linalg.norm(perp))
                for angle in [1.047, 2.094, 3.142, 4.189, 5.236]:
                    cos_a, sin_a = np.cos(angle), np.sin(angle)
                    rot_vec = (h_dir * cos_a
                               + np.cross(perp, h_dir) * sin_a
                               + perp * float(np.dot(perp, h_dir)) * (1 - cos_a))
                    rot_vec = rot_vec / max(float(np.linalg.norm(rot_vec)), 1e-12)
                    cand = parent_pos + rot_vec * 1.09
                    min_d_cand = min(
                        (float(np.linalg.norm(cand - _gp(conf, j3)))
                         for j3 in range(n_atoms)
                         if j3 != i and j3 != parent_idx),
                        default=999.0,
                    )
                    if min_d_cand > best_min_d:
                        best_min_d = min_d_cand
                        best_pos = cand
                        if best_min_d >= 1.0:
                            break

            _sp(conf, i, best_pos)

        if not any_fixed:
            break

    # 6c: Skip pi-coplanarity for multi-metal — it moves atoms and
    # can re-introduce clashes that Step 6a resolved.

    # Phase 6A (2026-05-12): conditional multi-hapto centroid re-correction.
    # Per Wave-3B forensics (/tmp/wave3b_multihapto_v2.md): secondary-metal
    # placement drags shared hapto atoms (C ∈ Cp(M1) AND σ-donor(M2)) toward
    # M2, collapsing M1-centroid distance.  Fe-Cp-Ni-bimetallic observed:
    # Fe-Cp2 1.65 → 1.012 Å (-39% off target).  This is the CORRECT location
    # for the Wave-2C fix (Phase 5B was at line 9965 in hybrid fallback path,
    # which multi-hapto SMILES never reach — explaining why Phase 5B failed).
    # Conditional gate: only fires when an actual centroid is >20% off target
    # AND env-flag set, so multi-sigma scaffolds (no hapto_groups) and
    # already-correct geometries are untouched.
    if (hapto_groups
            and _class_conditional_flag(
                "DELFIN_MULTI_HAPTO_RECORRECT", mol, default=0,
                default_classes=("multi_hapto",),
            )):
        try:
            conf_check = scaffold_mol.GetConformer(0)
            needs_correction = False
            for mi, catoms in hapto_groups:
                if not catoms:
                    continue
                mpos = _gp(conf_check, mi)
                cx = sum(_gp(conf_check, c)[0] for c in catoms) / len(catoms)
                cy = sum(_gp(conf_check, c)[1] for c in catoms) / len(catoms)
                cz = sum(_gp(conf_check, c)[2] for c in catoms) / len(catoms)
                dx = mpos[0] - cx; dy = mpos[1] - cy; dz = mpos[2] - cz
                d = (dx*dx + dy*dy + dz*dz) ** 0.5
                metal_sym = scaffold_mol.GetAtomWithIdx(mi).GetSymbol()
                target = _target_mc_dist(metal_sym, len(catoms))
                if target > 0 and abs(d - target) / target > 0.20:
                    needs_correction = True
                    break
            if needs_correction:
                _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
                logger.debug(
                    "multi-hapto centroid re-correction applied "
                    "(env DELFIN_MULTI_HAPTO_RECORRECT=1)"
                )
        except Exception as exc:
            logger.debug("multi-hapto re-correction skipped: %s", exc)

    logger.info(
        "Sequential multi-metal: %d/%d atoms placed",
        len(trusted), mol.GetNumAtoms(),
    )
    return scaffold_mol


def _secondary_metal_variant_plans(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    max_alternates_per_metal: int = 2,
) -> List[Tuple[str, Dict[int, int]]]:
    """Return ranked secondary-metal variant plans for hybrid hapto building."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return [("hybrid", {})]

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    plans: List[Tuple[str, Dict[int, int]]] = [("hybrid", {})]
    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        if atom.GetSymbol() not in _METAL_SET or metal_idx in hapto_metals:
            continue
        donor_count = sum(
            1
            for nbr in atom.GetNeighbors()
            if nbr.GetAtomicNum() > 1 and nbr.GetSymbol() not in _METAL_SET
        )
        if donor_count < 2:
            continue
        for rank in range(1, max(1, int(max_alternates_per_metal)) + 1):
            plans.append(
                (f"hybrid-secondary-{atom.GetSymbol()}{metal_idx}-alt{rank}", {metal_idx: rank})
            )
    return plans


def _build_hybrid_hapto_complex(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    preview_store: Optional[List[Tuple[str, str]]] = None,
    secondary_variant_plan: Optional[Dict[int, int]] = None,
):
    """Assemble a hapto complex from a scaffold plus rigid ligand fragments."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return None

    decomposition = _decompose_hapto_complex(mol, hapto_groups)
    if decomposition is None:
        return None

    # Multi-metal: use sequential builder
    n_metals = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET)
    if n_metals >= 2 and decomposition is not None:
        try:
            seq_result = _build_multimetal_hapto_sequential(
                mol,
                hapto_groups,
                decomposition,
                secondary_variant_plan=secondary_variant_plan,
            )
        except Exception as e:
            logger.debug("Sequential multi-metal builder failed: %s", e)
            seq_result = None
        if seq_result is not None:
            return seq_result
        logger.info("Sequential multi-metal builder fell back to standard path")

    primary_module = _detect_primary_organometal_module(decomposition, hapto_groups)
    if primary_module is not None:
        try:
            primary_candidate = _build_primary_organometal_hapto_module(
                mol,
                decomposition,
                hapto_groups,
                primary_module,
                preview_store=preview_store,
            )
        except Exception as e:
            logger.debug("Hybrid hapto primary organometal builder skipped: %s", e)
            primary_candidate = None
        if primary_candidate is not None:
            return primary_candidate
        logger.info(
            "Hybrid hapto primary organometal module builder fell back to standard scaffold assembly for %s%d",
            mol.GetAtomWithIdx(primary_module.metal_idx).GetSymbol(),
            primary_module.metal_idx,
        )
        if preview_store is not None and not preview_store:
            _store_hapto_preview_candidate(
                mol,
                hapto_groups,
                preview_store,
                label='quick-hapto-preview',
            )
    elif preview_store is not None:
        _store_hapto_preview_candidate(
            mol,
            hapto_groups,
            preview_store,
            label='quick-hapto-preview',
        )

    scaffold_mol = Chem.Mol(mol)
    if not _build_hapto_scaffold(scaffold_mol, hapto_groups):
        return None

    aligned_fragments = 0
    scaffold_only_fragments = 0
    deferred_secondary_fragments = 0
    trusted_atoms: set = set()
    hapto_metal_indices = {metal_idx for metal_idx, _grp in hapto_groups}
    for fragment in decomposition.fragments:
        if fragment.use_scaffold_only:
            scaffold_only_fragments += 1
            trusted_atoms.update(fragment.atom_indices)
            logger.debug(
                "Hybrid hapto: keeping scaffold geometry for fragment "
                "(atoms=%d, metals=%d, bridging_donors=%d)",
                len(fragment.atom_indices),
                len(fragment.metal_neighbor_indices),
                len(fragment.bridging_donor_indices),
            )
            continue
        if any(metal_idx not in hapto_metal_indices for metal_idx in fragment.metal_neighbor_indices):
            deferred_secondary_fragments += 1
            logger.debug(
                "Hybrid hapto: deferring fragment with secondary-metal coordination "
                "(atoms=%d, metals=%s) to secondary module assembly",
                len(fragment.atom_indices),
                sorted(fragment.metal_neighbor_indices),
            )
            continue
        embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded_fragment is None:
            continue
        if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded_fragment):
            aligned_fragments += 1
            trusted_atoms.update(fragment.atom_indices)

    if aligned_fragments == 0 and scaffold_only_fragments == 0:
        return None

    try:
        _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    secondary_metals, secondary_module_atoms, relaxable_secondary_atoms = _assemble_secondary_metal_coordination_modules(
        scaffold_mol,
        decomposition,
        hapto_groups,
        fit_variant_plan=secondary_variant_plan,
    )
    trusted_atoms |= secondary_module_atoms

    if not secondary_metals:
        secondary_metals = _place_secondary_metals_in_hapto_fragments(
            scaffold_mol,
            hapto_groups,
        )

    if secondary_metals:
        logger.debug(
            "Hybrid hapto built %d explicit secondary metal module(s)",
            len(secondary_metals),
        )
    else:
        logger.debug("Hybrid hapto used donor-cloud fallback for secondary metals")

    if aligned_fragments + scaffold_only_fragments < len(decomposition.fragments):
        try:
            _propagate_non_hapto_atoms(
                scaffold_mol,
                0,
                hapto_groups,
                extra_fixed_indices=trusted_atoms | _hapto_primary_donor_indices(scaffold_mol, hapto_groups) | secondary_metals,
            )
        except Exception as e:
            logger.debug("Hybrid hapto non-hapto propagation skipped: %s", e)
        if not secondary_metals:
            secondary_metals |= _place_secondary_metals_in_hapto_fragments(
                scaffold_mol,
                hapto_groups,
            )

    try:
        _enforce_donor_pi_coplanarity(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    try:
        _refine_hybrid_hapto_complex(
            scaffold_mol,
            decomposition,
            hapto_groups,
            trusted_atom_indices=trusted_atoms,
            secondary_metal_indices=secondary_metals,
            apply_final_hapto_correction=not bool(secondary_metals),
        )
    except Exception as e:
        logger.debug("Hybrid hapto local refinement skipped: %s", e)

    try:
        _enforce_donor_pi_coplanarity(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    # Phase 5B REVERTED (2026-05-12): unconditional _correct_hapto_geometry
    # caused multi-sigma regression 42.9% → 0% (5 failed instead of 1)
    # because the function disturbs multi-sigma scaffolds that have no
    # hapto groups.  Per Wave-2C forensics the gate should be
    # CONDITIONAL on Cp-centroid > 20% off target, not unconditional.
    # Reverting to original gate.  Surgical conditional fix is future
    # work; multi-hapto class stays at 0% topology pass for now.
    if not secondary_metals:
        try:
            _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
        except Exception:
            pass

    # Final clash resolution: _enforce_donor_pi_coplanarity may have moved
    # atoms (especially H) into metal clash zones.
    try:
        _final_clash_resolution(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    logger.info(
        "Hybrid hapto assembly aligned %d fragment(s), scaffold-kept %d, deferred-secondary %d/%d",
        aligned_fragments,
        scaffold_only_fragments,
        deferred_secondary_fragments,
        len(decomposition.fragments),
    )

    # Bond-length repair for severely distorted non-metal bonds, then quality gate
    try:
        import numpy as np
        conf = scaffold_mol.GetConformer(0)
        _expected = {
            frozenset(['C', 'C']): 1.45, frozenset(['C', 'N']): 1.40,
            frozenset(['C', 'O']): 1.35, frozenset(['C', 'S']): 1.78,
            frozenset(['N', 'N']): 1.40, frozenset(['N', 'O']): 1.35,
            frozenset(['O', 'O']): 1.45, frozenset(['C', 'Si']): 1.87,
        }
        hapto_frozen = set()
        for _mi, catoms in hapto_groups:
            hapto_frozen.add(_mi)
            hapto_frozen.update(catoms)
        for _ai in range(scaffold_mol.GetNumAtoms()):
            if scaffold_mol.GetAtomWithIdx(_ai).GetSymbol() in _METAL_SET:
                hapto_frozen.add(_ai)

        def _gp2(i):
            p = conf.GetAtomPosition(i)
            return np.array([p.x, p.y, p.z])
        def _sp2(i, pos):
            conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

        # Build H groups
        h_of2: Dict[int, List[int]] = {}
        for _ai in range(scaffold_mol.GetNumAtoms()):
            a = scaffold_mol.GetAtomWithIdx(_ai)
            if a.GetAtomicNum() == 1:
                for nbr in a.GetNeighbors():
                    if nbr.GetAtomicNum() > 1:
                        h_of2.setdefault(nbr.GetIdx(), []).append(_ai)
                        break

        # Iterative bond-length repair (60 iterations)
        for _rep_pass in range(60):
            max_err = 0.0
            for bond in scaffold_mol.GetBonds():
                bi = bond.GetBeginAtomIdx()
                bj = bond.GetEndAtomIdx()
                si = scaffold_mol.GetAtomWithIdx(bi).GetSymbol()
                sj = scaffold_mol.GetAtomWithIdx(bj).GetSymbol()
                if si in _METAL_SET or sj in _METAL_SET:
                    continue
                pi, pj = _gp2(bi), _gp2(bj)
                d = float(np.linalg.norm(pj - pi))
                if si == 'H' or sj == 'H':
                    expected = 1.08
                else:
                    expected = _expected.get(frozenset([si, sj]), 1.50)
                err = abs(d - expected)
                if err < 0.05:
                    continue
                if err > max_err:
                    max_err = err
                if d < 1e-8:
                    continue
                direction = (pj - pi) / d
                correction = 0.3 * (d - expected)
                bi_frozen = bi in hapto_frozen
                bj_frozen = bj in hapto_frozen
                if bi_frozen and bj_frozen:
                    continue
                if bi_frozen:
                    delta = -correction * direction
                    _sp2(bj, pj + delta)
                    for hi in h_of2.get(bj, []):
                        if hi not in hapto_frozen:
                            _sp2(hi, _gp2(hi) + delta)
                elif bj_frozen:
                    delta = correction * direction
                    _sp2(bi, pi + delta)
                    for hi in h_of2.get(bi, []):
                        if hi not in hapto_frozen:
                            _sp2(hi, _gp2(hi) + delta)
                else:
                    _sp2(bi, pi + 0.5 * correction * direction)
                    _sp2(bj, pj - 0.5 * correction * direction)
            if max_err < 0.1:
                break

        # Re-center secondary metals after bond repair shifted donors
        _fix_secondary_metal_distances(scaffold_mol, 0, hapto_groups)

        # Re-run clash resolution after bond repair + metal recentering
        _final_clash_resolution(scaffold_mol, 0, hapto_groups)

        # Quality gate after repair
        for bond in scaffold_mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            si = scaffold_mol.GetAtomWithIdx(bi).GetSymbol()
            sj = scaffold_mol.GetAtomWithIdx(bj).GetSymbol()
            if si in _METAL_SET or sj in _METAL_SET:
                continue
            pi, pj = _gp2(bi), _gp2(bj)
            d = float(np.linalg.norm(pj - pi))
            if si == 'H' or sj == 'H':
                expected = 1.08
            else:
                expected = _expected.get(frozenset([si, sj]), 1.50)
            ratio = d / expected
            if ratio < 0.50 or ratio > 3.00:
                logger.info(
                    "Hybrid hapto rejected: bond %s[%d]-%s[%d] = %.3f A (expected ~%.2f)",
                    si, bi, sj, bj, d, expected,
                )
                return None
    except Exception:
        pass

    return scaffold_mol


def _final_clash_resolution(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
) -> None:
    """Push apart non-bonded atoms that are too close after all geometry steps.

    Metals and hapto ring atoms are frozen; everything else can move.
    When a heavy atom is pushed, its bonded H atoms move rigidly with it.
    Runs after _enforce_donor_pi_coplanarity which may introduce clashes.
    """
    import numpy as np

    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    n_atoms = mol.GetNumAtoms()

    # PERF (BYTE-IDENTICAL): mirror the conformer into a local float64 array.
    # _gp returned ``np.array([p.x, p.y, p.z])`` (float64) and _sp wrote
    # ``Point3D(float(x), float(y), float(z))``; the conformer stores doubles,
    # so the read→array→write round-trip is bit-identical (verified
    # array_equal).  The clash passes mutate positions IN PLACE mid-pass (a
    # push to one pair is visible to later pairs in the same pass), so the loop
    # is inherently sequential and is kept scalar in the exact original order —
    # only the per-(get/set) RDKit boundary crossing (the dominant Python cost,
    # ~13× a local-array access) is removed by reading/writing the local array
    # and flushing to the conformer once at the end.  The per-pair distance is
    # ``sqrt(diff.dot(diff))`` which is what ``np.linalg.norm`` reduces to on a
    # 1-D vector (same BLAS ddot kernel → bit-identical, verified array_equal)
    # and is the fastest scalar form (matmul is bit-identical too but slower for
    # a single 3-vector).
    coords = np.empty((n_atoms, 3), dtype=float)
    for i in range(n_atoms):
        p = conf.GetAtomPosition(i)
        coords[i, 0] = p.x
        coords[i, 1] = p.y
        coords[i, 2] = p.z

    def _gp(i):
        # fresh copy so caller-side ``+`` does not alias the stored row
        return coords[i].copy()

    def _sp(i, pos):
        coords[i, 0] = float(pos[0])
        coords[i, 1] = float(pos[1])
        coords[i, 2] = float(pos[2])

    def _move_with_h(ai, displacement):
        """Move atom ai and its bonded H atoms by displacement."""
        _sp(ai, _gp(ai) + displacement)
        for hi in h_of.get(ai, []):
            if hi not in frozen:
                _sp(hi, _gp(hi) + displacement)

    frozen = set()
    for mi, catoms in hapto_groups:
        frozen.add(mi)
        frozen.update(catoms)
    for ai in range(n_atoms):
        if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
            frozen.add(ai)

    bonded_pairs = set()
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((bi, bj))
        bonded_pairs.add((bj, bi))

    # Build H→parent and parent→H maps
    h_of: Dict[int, List[int]] = {}
    h_parent: Dict[int, int] = {}
    for ai in range(n_atoms):
        a = mol.GetAtomWithIdx(ai)
        if a.GetAtomicNum() == 1:
            for nbr in a.GetNeighbors():
                if nbr.GetAtomicNum() > 1:
                    h_parent[ai] = nbr.GetIdx()
                    h_of.setdefault(nbr.GetIdx(), []).append(ai)
                    break

    # Only check heavy atoms (H moves with parent)
    heavy = [i for i in range(n_atoms) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    rng = np.random.default_rng(99)

    # Hoisted per-call invariants (recomputed identically every pass before):
    # per-atom symbol + metal flag, and the H-clash candidate list.
    _sym = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
    _is_metal_atom = [s in _METAL_SET for s in _sym]
    metal_frozen = [j for j in frozen if _is_metal_atom[j]]
    h_candidates = []
    for i in range(n_atoms):
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 1 or i in frozen:
            continue
        parent = h_parent.get(i)
        if parent is not None and parent in frozen:
            continue  # parent is frozen, can't fix
        h_candidates.append((i, parent))

    for _pass in range(30):
        any_push = False
        # Heavy-heavy clashes
        for ii, i in enumerate(heavy):
            i_frozen = i in frozen
            i_metal = _is_metal_atom[i]
            for j in heavy[ii + 1:]:
                if (i, j) in bonded_pairs:
                    continue
                if i_frozen and j in frozen:
                    continue
                pi, pj = _gp(i), _gp(j)
                diff = pi - pj
                d = float(np.sqrt(diff.dot(diff)))
                if i_metal or _is_metal_atom[j]:
                    min_d = 2.0
                else:
                    min_d = 1.2
                if d >= min_d:
                    continue
                if d < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pj - pi) / d
                gap = min_d - d
                if i_frozen:
                    _move_with_h(j, gap * direction)
                elif j in frozen:
                    _move_with_h(i, -gap * direction)
                else:
                    _move_with_h(i, -0.5 * gap * direction)
                    _move_with_h(j, 0.5 * gap * direction)
                any_push = True

        # H-metal clashes: move the parent heavy atom (+ all its H) so
        # the offending H clears the metal without breaking C-H bonds.
        for i, parent in h_candidates:
            for j in metal_frozen:
                if (i, j) in bonded_pairs:
                    continue
                pi, pj = _gp(i), _gp(j)
                diff = pi - pj
                d = float(np.sqrt(diff.dot(diff)))
                if d >= 2.0:
                    continue
                if d < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pi - pj) / d
                displacement = (2.0 - d) * direction
                if parent is not None and parent not in frozen:
                    _move_with_h(parent, displacement)
                else:
                    _sp(i, _gp(i) + displacement)
                any_push = True

        if not any_push:
            break

    # Flush the local array back to the conformer (bit-identical write-back).
    for i in range(n_atoms):
        conf.SetAtomPosition(i, Point3D(
            float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])))


def _fix_secondary_metal_distances(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
) -> None:
    """Fix non-hapto metal positions by moving them into their donor cloud.

    ETKDG doesn't understand metal-ligand bonds, so metal positions are random.
    For each non-hapto metal, compute the optimal position given its donors'
    current positions and target M-L distances, then move the metal there.
    This avoids moving donor atoms (which may share fragments).
    """
    import numpy as np

    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, pos):
        conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    hapto_metals = {mi for mi, _ in hapto_groups}

    # Collect ALL hapto-bound atoms across every hapto group.  These belong
    # to another metal's η-coordination sphere and MUST NOT be dragged by
    # the donor-pull below — doing so collapses the Cp/arene ring onto the
    # wrong metal (e.g. Fe5 η5-Cp ring 2 in Fe-Cp-Ni-bimetallic case
    # diagnosed 2026-05-12 — see /tmp/agent3_multi_hapto_recovery.md).
    all_hapto_atoms: set = set()
    for _mi_h, _catoms in hapto_groups:
        for _a in _catoms:
            all_hapto_atoms.add(_a)

    for atom in mol.GetAtoms():
        mi = atom.GetIdx()
        msym = atom.GetSymbol()
        if msym not in _METAL_SET or mi in hapto_metals:
            continue

        donors = []
        targets = []
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                continue
            dsym = nbr.GetSymbol()
            donors.append(ni)
            targets.append(float(_get_ml_bond_length(msym, dsym)))

        if len(donors) < 2:
            continue

        # Iteratively move metal toward the position that minimizes
        # sum of (actual_dist - target_dist)^2
        mpos = _gp(mi)
        for _it in range(80):
            grad = np.zeros(3)
            for di, target in zip(donors, targets):
                dpos = _gp(di)
                vec = mpos - dpos
                dist = float(np.linalg.norm(vec))
                if dist < 1e-8:
                    continue
                grad += 2.0 * (dist - target) * vec / dist
            step = -0.15 * grad
            if float(np.linalg.norm(step)) < 0.001:
                break
            mpos = mpos + step

        _sp(mi, mpos)

        # If any donor is still >30% off target, move it toward the metal.
        # Use a limited BFS (max 6 atoms) to avoid dragging shared fragments.
        for di, target in zip(donors, targets):
            dpos = _gp(di)
            dist = float(np.linalg.norm(dpos - mpos))
            if dist < 1e-8 or abs(dist - target) / target < 0.30:
                continue
            # Skip drag if donor is a hapto atom of ANOTHER metal — pulling
            # it would shear the Cp/arene ring off its η-metal.  Verified
            # on Fe-Cp-Ni-bimetallic: Fe5-centroid₂ collapsed 1.65 → 1.012 Å
            # before this guard.  See /tmp/agent3_multi_hapto_recovery.md.
            if di in all_hapto_atoms:
                continue
            direction = (dpos - mpos) / dist
            delta = direction * (target - dist)
            # Limited BFS: donor + up to 5 non-metal, non-hapto neighbors
            local_frag = {di}
            bfs_q = [di]
            while bfs_q and len(local_frag) < 7:
                cur = bfs_q.pop(0)
                for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                    ni2 = nbr.GetIdx()
                    if ni2 in local_frag or ni2 in hapto_metals:
                        continue
                    if nbr.GetSymbol() in _METAL_SET:
                        continue
                    # Don't cross into another donor's territory
                    if ni2 in set(donors) and ni2 != di:
                        continue
                    # Belt-and-braces: don't drag hapto atoms of another
                    # metal even if BFS reaches them via a sigma neighbor.
                    if ni2 in all_hapto_atoms:
                        continue
                    local_frag.add(ni2)
                    bfs_q.append(ni2)
            for ai in local_frag:
                _sp(ai, _gp(ai) + delta)
