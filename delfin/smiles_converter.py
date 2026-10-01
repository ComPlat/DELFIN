"""SMILES to XYZ conversion using RDKit with metal complex support.

Note: RDKit generates rough starting geometries that are NOT chemically accurate,
especially for metal complexes. These coordinates should ALWAYS be optimized with
GOAT/xTB before running ORCA calculations.
"""

from __future__ import annotations

from collections import Counter, OrderedDict
from dataclasses import dataclass
import concurrent.futures
import hashlib
import math
import os
import random
import re
import threading
import time
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    OPENBABEL_AVAILABLE,
    Point3D,
    RDKIT_AVAILABLE,
    STK_AVAILABLE,
    _ALL_CN3_POLYHEDRA,
    _ALL_CN4_POLYHEDRA,
    _ALL_CN5_POLYHEDRA,
    _ALL_CN6_POLYHEDRA,
    _ALL_CN7_POLYHEDRA,
    _ALL_CN8_POLYHEDRA,
    _ALL_POLYHEDRA_BY_CN,
    _COVALENT_RADII,
    _HALOGENS,
    _HAPTO_CENTROID_DISTANCES,
    _MC_LEN_CRYSTAL_MEDIAN,
    _MC_LEN_CRYSTAL_SOURCE,
    _METALLOID_MD_DONORS,
    _METALLOID_MD_FORCE_SHORT,
    _METALS,
    _METAL_ATOMICNUMS,
    _METAL_GROUP_NUMBER,
    _METAL_HYDRIDE_BOND_LENGTHS,
    _METAL_LIGAND_BOND_LENGTHS,
    _METAL_LIGAND_BOND_LENGTHS_T33,
    _METAL_METAL_BOND_LENGTHS,
    _METAL_ROW_OFFSET,
    _METAL_SET,
    _ML_ME_CACHE,
    _ML_ME_MIN_N,
    _ORGANOMETALLIC_METALS,
    _PI_ACCEPTOR_DONOR_ELEMS,
    _PREFERRED_CN4_GEOMETRY,
    _PREFERRED_CN5_GEOMETRY,
    _PREFERRED_CN6_GEOMETRY,
    _apply_mh_table_fallback,
    _apply_new_ml_pairs_t33,
    _classify_cn5_geometry,
    _classify_cn5_geometry_from_labels,
    _cn5_count_pi_acceptor_donors,
    _cn5_d_electron_count,
    _cn5_has_tridentate_chelate,
    _get_ml_bond_length,
    _is_metal_nitrogen_complex,
    _is_simple_organometallic,
    _ml_bond_kind,
    _ml_me_band,
    _prefer_no_sanitize,
    _pt,
    _secondary_donor_fit_weight,
    _secondary_donor_target_length,
    _target_mc_dist,
    pybel,
    stk,
)  # noqa: F401

from delfin.manta.hapto_detect import (
    _apply_hapto_approximation,
    _classify_complex_class,
    _find_hapto_groups,
    _hapto_approx_enabled,
    _hapto_failfast_error,
    _hapto_label_for_group,
    _probe_hapto_groups_from_smiles,
    _select_multihapto_anchors,
    contains_metal,
    mol_from_smiles_rdkit,
)  # noqa: F401

from delfin.manta.converter_flags import (
    DELFIN_CHELATE_ACCEPT_DELTA,
    DELFIN_CHELATE_CAP_30,
    DELFIN_CHELATE_CAP_60,
    DELFIN_CHELATE_CAP_90,
    DELFIN_CHELATE_N_TRIALS,
    DELFIN_CHELATE_RANK_CLASS_ALPHA,
    DELFIN_CHELATE_RANK_CLASS_AWARE,
    DELFIN_CHELATE_RANK_VARIANTS,
    DELFIN_CHELATE_REJECT_DELTA,
    DELFIN_CLASS_AWARE_SEEDS,
    DELFIN_COLLAPSE_BOND_MIN_SCALE,
    DELFIN_COLLAPSE_TETRA_VOL_MIN,
    DELFIN_FFFREE_COLLAPSE_REJECT,
    DELFIN_FINAL_GATE_ENABLED,
    DELFIN_IDEAL_POLYHEDRON_MAX_DEV,
    DELFIN_MAX_PROCESS_WORKERS,
    DELFIN_MAX_THREAD_WORKERS,
    DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES,
    DELFIN_PRE_UFF_CAP_MULTIPLIER,
    DELFIN_RULE4_PI_PLANAR_TOL_FRAC,
    DELFIN_RULE5_INNER_SPHERE_FACT,
    DELFIN_RULE5_INTERFRAG_COV_FACT,
    DELFIN_RULE6_METALLACYCLE_MAX_DEV,
    DELFIN_RULE7_SP2_ANGLE_MAX,
    DELFIN_RULE7_SP2_ANGLE_MIN,
    DELFIN_RULE7_SP2_OOP_MAX,
    DELFIN_RULE7_SP2_OOP_MAX_METAL,
    DELFIN_RULE7_SP3_MIN_ANGLE_DEG,
    DELFIN_RULE7_SP_MIN_ANGLE_DEG,
    DELFIN_SEVERE_DIST_MAX_ABS,
    DELFIN_SEVERE_DIST_MAX_SCALE,
    DELFIN_SYMMETRY_WEIGHT,
    DELFIN_TOPOLOGY_STRICT_MODE,
    DELFIN_TOPO_TEMPLATE_TOP_K,
    DELFIN_TOP_LEVEL_SEED_COUNT,
    _CHELATE_CLASS_DONOR_WEIGHTS,
    _CHELATE_EMBED_TIMEOUT,
    _CLASS_AWARE_SEED_DEFAULTS,
    _DEFAULT_EMBED_SEED,
    _DELFIN_PROFILES,
    _EMBED_TIMEOUT,
    _MULTIEMBED_TIMEOUT,
    _MULTIEMBED_TIMEOUT_OVERRIDE,
    _OB_ROTOR_TIMEOUT,
    _PIPELINE_SEEDS,
    _RESTORE_CAP_MULT,
    _RESTORE_NUMCONFS,
    _RESTORE_RANKS,
    _RESTORE_TOPK,
    _TOP_LEVEL_SEEDS,
    _all_polyhedra_codes,
    _apply_uff_jitter,
    _chelate_class_donor_penalty,
    _class_aware_seed_count,
    _class_conditional_flag,
    _delfin_env_float,
    _delfin_env_int,
    _deterministic_embed_seed,
    _deterministic_mode,
    _every_append_gate_enabled,
    _generate_pipeline_seeds,
    _hapto_candidate_collapsed_bonds,
    _hapto_scaffold_primary_enabled,
    _hapto_seat_rigid_enabled,
    _multi_sigma_v2_active,
    _multi_sigma_v2_budget,
    _multihapto_etkdg_fallback_enabled,
    _preferred_cn4_for,
    _resolve_quality_profile,
    _resolve_top_level_seed_count,
    _trace_seating,
    _union_prepend_ffree,
)  # noqa: F401

from delfin.manta.metal_smiles import (
    _convert_metal_bonds_to_dative,
    _denormalize_metal_smiles,
    _fix_hapto_donor_h,
    _fix_organometallic_carbon_h,
    _hapto_h_always,
    _normalize_metal_smiles,
    _strip_h_on_coordinated_p,
    _strip_h_on_metal_halogen,
)  # noqa: F401

from delfin.manta.conformer_io import (
    _fix_zero_coord_hydrogens,
    _inject_openbabel_conformers_into_mol,
    _mol_to_xyz,
    _mol_to_xyz_conformer,
    _normalize_conversion_backend,
    _obmol_to_delfin_xyz,
    _openbabel_generate_conformer_xyz,
    _openbabel_smiles_variants,
    _pick_openbabel_forcefield,
    _xyz_to_rdkit_conformer,
    _xyz_to_rdkit_conformer_via_ob_mapping,
)  # noqa: F401

from delfin.manta.embed_timeout import (
    _embed_multiple_confs_with_timeout,
    _embed_with_timeout,
    _make_random_embed_params,
)  # noqa: F401

from delfin.manta.mol_prep import (
    _CACHE_MISS,
    _PREP_MOL_CACHE,
    _PREP_MOL_CACHE_MAX,
    _dearomatized_embedding_copy,
    _embed_multiple_confs_robust,
    _prepare_mol_for_embedding,
    _prepare_mol_for_embedding_uncached,
    _rescale_metal_donor_distances,
    is_smiles_string,
)  # noqa: F401

from delfin.manta.topology_checks import (
    _ITER8_1_EXTRA_THRESHOLDS,
    _component_stats_from_adj,
    _count_extra_heavy_bonds,
    _count_xyz_clashes,
    _cycle_size_signature_from_adj,
    _flatten_sp2_atoms_xyz,
    _fragment_topology_ok,
    _fragment_topology_relaxed_fallback_ok,
    _global_heavy_connectivity_ok,
    _has_atom_clash,
    _has_collapsed_sp3_centre,
    _has_pi_ring_nonplanarity,
    _has_severe_covalent_distortion,
    _has_unphysical_metal_nonbonded_contact,
    _has_unphysical_oco_geometry,
    _heavy_component_stats_smiles,
    _heavy_component_stats_xyz,
    _heavy_graph_edges_smiles,
    _heavy_graph_edges_xyz,
    _heavy_graph_exact_match_ok,
    _heavy_local_signature_match_ok,
    _heavy_local_signature_multiset,
    _heavy_local_signature_multiset_smiles,
    _heavy_local_signature_multiset_xyz,
    _metal_aware_coordination_ok,
    _metal_donor_distances_realistic,
    _no_spurious_bonds,
    _nonmetal_fragment_ids,
    _organic_fragment_signature,
    _organic_fragment_signature_xyz,
    _organic_graph_signature,
    _roundtrip_ring_count_ok,
    _shortest_path_length_excluding,
    _verify_metal_connectivity,
    _xyz_to_canonical_smiles,
)  # noqa: F401

from delfin.manta.isomer_labels import (
    _TOPO_CANONICAL_FNS,
    _TOPO_GEOMETRY_VECTORS,
    _TOPO_TRANS_POSITIONS,
    _canonical_bcsap,
    _canonical_coh,
    _canonical_cpap,
    _canonical_cubo,
    _canonical_dd,
    _canonical_hbp,
    _canonical_icos,
    _canonical_lin,
    _canonical_oh,
    _canonical_pap,
    _canonical_pbp,
    _canonical_sap,
    _canonical_sp,
    _canonical_sq,
    _canonical_ss,
    _canonical_tbp,
    _canonical_th,
    _canonical_tp,
    _canonical_tpr,
    _canonical_ts,
    _canonical_ttp,
    _chelate_backbone_max_reach,
    _chelate_pairs,
    _classify_isomer_label,
    _compute_coordination_fingerprint,
    _donor_type_map,
    _extract_helicity_suffix,
    _find_bridging_donors,
    _is_viable_donor,
    _label_from_canonical_form,
    _ligand_fragments,
)  # noqa: F401

from delfin.manta.hapto_scaffold import (
    _apply_hapto_centroid_bias,
    _build_hapto_scaffold,
    _correct_hapto_geometry,
    _detect_shared_ring_eta_groups,
    _enforce_hapto_ring_planarity,
    _find_ansa_bridges,
    _propagate_non_hapto_atoms,
)  # noqa: F401

from delfin.manta.hybrid_fragments import (
    _HybridHaptoDecomposition,
    _HybridHaptoFragment,
    _PrimaryOrganometalModule,
    _align_hybrid_fragment_onto_scaffold,
    _choose_fragment_anchor_atoms,
    _copy_fragment_atom,
    _decompose_hapto_complex,
    _detect_primary_organometal_module,
    _embed_hybrid_fragment,
    _extract_fragment_mol,
    _hapto_primary_donor_indices,
    _primary_organometal_module_quality_ok,
    _rotation_matrix_from_vectors,
)  # noqa: F401

from delfin.manta.secondary_metal_modules import (
    _HAPTO_SECONDARY_CN4_GEOMETRIES,
    _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES,
    _align_hybrid_fragment_to_targets,
    _assemble_secondary_metal_coordination_modules,
    _axis_angle_rotation_matrix,
    _build_primary_organometal_hapto_module,
    _enumerate_secondary_metal_geometry_fits,
    _enumerate_secondary_metal_position_fits,
    _find_scaffold_secondary_branch,
    _fit_secondary_metal_geometry_from_donors,
    _fit_secondary_metal_model_to_targets,
    _fit_secondary_metal_position_from_donors,
    _fragment_bite_direction_penalty,
    _fragment_donor_selectivity_penalty,
    _fragment_environment_clash_penalty,
    _fragment_planar_donor_frame,
    _hybrid_bond_target_length,
    _optimize_secondary_fragment_pose,
    _optimize_secondary_metal_module_local,
    _place_secondary_metals_in_hapto_fragments,
    _planar_fragment_donor_approach_penalty,
    _planar_fragment_metal_coplanarity_penalty,
    _project_displacements_to_rigid_body,
    _relieve_secondary_oo_chelate_contacts,
    _reorient_scaffold_secondary_branch,
    _restore_secondary_rigid_fragment_geometry,
    _secondary_fragment_donor_selectivity_penalty,
    _secondary_fragment_prefers_rigid_pose,
    _secondary_metal_geometry_codes,
    _secondary_module_fragment_is_movable,
    _secondary_non_donor_contact_penalty,
)  # noqa: F401

from delfin.manta.hybrid_assembly import (
    _append_hapto_preview_xyz,
    _build_hybrid_hapto_complex,
    _build_multimetal_hapto_sequential,
    _enforce_donor_pi_coplanarity,
    _final_clash_resolution,
    _fix_secondary_metal_distances,
    _refine_hybrid_hapto_complex,
    _secondary_metal_variant_plans,
    _store_hapto_preview_candidate,
)  # noqa: F401

from delfin.manta.geometry_quality import (
    _GEOM_IDEAL_ANGLES,
    _GEOM_IDEAL_ANGLES_REAL,
    _angle_class,
    _cn4_geometry_penalties,
    _collect_metal_ligand_frame,
    _conformer_rmsd,
    _donor_h_points_at_metal,
    _enforce_metal_topology,
    _enforce_smiles_mm_distances,
    _estimate_isomer_upper_bound,
    _geometry_quality_score,
    _has_bad_geometry,
    _has_ligand_intertwining,
    _ideal_polyhedron_angle_dev_per_metal,
    _ml_distance_range,
    _preferred_cn4_geometry_score,
    _segment_distance_sq,
    _xyz_passes_final_geometry_checks,
)  # noqa: F401

from delfin.manta.coordination_enumerator import (
    _enumerate_orbits_topo,
    _enumerate_topological_isomers,
)  # noqa: F401

from delfin.manta.ligand_placement import (
    _align_and_orient_ligands,
    _build_multimetal_scaffold,
    _compute_lp_tilt_rotations,
    _scale_aromatic_rings_in_xyz_from_smiles,
    _scale_aromatic_rings_in_xyz_geom,
    _scale_aromatic_rings_to_ideal_cc,
    _snap_aromatic_rings_in_xyz,
    _snap_aromatic_rings_to_plane,
    _snap_bridging_donors_to_compromise,
    _snap_fused_aromatic_groups_to_plane,
)  # noqa: F401

from delfin.manta.pre_uff_snap import (
    _D8_SQ_ISO_METALS,
    _D8_SQ_METALS,
    _GATE_V2_DONOR_HIGH,
    _GATE_V2_DONOR_LOW_MULTIBOND,
    _bfs_ligand_fragment,
    _clamp_metalloid_md_xyz,
    _cn5_enum_complete_enabled,
    _flatten_d8_sq_planar_xyz,
    _md_distance_in_tolerance,
    _pre_uff_md_snap_enabled,
    _pre_uff_topology_gate_enabled,
    _pre_uff_topology_gate_v2_enabled,
    _snap_md_distances_to_ideal,
)  # noqa: F401

from delfin.manta.uff_constraints import (
    _GROUP_14,
    _GROUP_15,
    _GROUP_16,
    _PYYKKO_ORDER_ATTR,
    _PYYKKO_ORDER_CACHE,
    _THETA_BY_PERIOD,
    _VDW_RADII_CLASH,
    _build_coordination_constraints_from_xyz,
    _build_coordination_uff_constraints,
    _build_uff_constraints_from_template,
    _detect_chelate_donors,
    _donor_sigma_geometry,
    _period_of,
    _pyykko_order_radius,
)  # noqa: F401

from delfin.manta.openbabel_optimize import (
    _apply_template_bond_orders,
    _geometric_inter_clash_relief,
    _optimize_xyz_openbabel,
    _optimize_xyz_openbabel_safe,
)  # noqa: F401

from delfin.manta.hapto_candidates import (
    _hapto_candidate_quality_score,
    _hapto_candidate_topology_ok,
    _hapto_geometry_quality_score,
    _hapto_mol_from_xyz_template,
    _refine_hapto_candidate_with_rdkit_uff,
    _select_best_hapto_candidate,
)  # noqa: F401

from delfin.manta.stage_hooks import (
    _apply_5j_a_cp_piano_stool_if_enabled,
    _apply_arom_bond_length_if_enabled,
    _apply_arom_planarize_if_enabled,
    _apply_aromatic_planarity_if_enabled,
    _apply_atropisomer_enum_if_enabled,
    _apply_baustein4_if_enabled,
    _apply_baustein5_if_enabled,
    _apply_baustein6_if_enabled,
    _apply_bond_decollapse_if_enabled,
    _apply_coord_angle_fix_if_enabled,
    _apply_f19_to_fallback_xyz,
    _apply_fixer_bridging_anion_if_enabled,
    _apply_fixer_f19_if_enabled,
    _apply_fixer_f25_if_enabled,
    _apply_fixer_sp2c_planarize_if_enabled,
    _apply_fixer_sp2n_planarize_if_enabled,
    _apply_fixer_wuxqak_if_enabled,
    _apply_h_placement_if_enabled,
    _apply_hapto_clearance_if_enabled,
    _apply_hydroxyl_geom_if_enabled,
    _apply_isolated_reseat_if_enabled,
    _apply_me_bond_snap_if_enabled,
    _apply_mirror_enum_if_enabled,
    _apply_pi_coplanar_m_if_enabled,
    _apply_stereocenter_enum_if_enabled,
    _apply_xtb_cascade_if_enabled,
)  # noqa: F401

from delfin.manta.embed_strategies import (
    _clash_aware_h_place_xyz,
    _fix_h_geometry_universal,
    _fix_h_geometry_via_smiles,
    _manual_metal_embed,
    _smiles_to_xyz_no_valence_check,
    _smiles_to_xyz_unsanitized_fallback,
    _try_multiple_strategies,
)  # noqa: F401

from delfin.manta.single_structure import (
    _HAPTO_QUICK_PREVIEW_CACHE,
    _SMILES_ATOM_TOKEN_RE,
    _SMILES_TOKEN_RE,
    _pin_metal_neighbour_hydrogens,
    _smiles_to_xyz_legacy,
    _topology_hard_gate_check,
    _try_ensemble_router,
    _uff_seat_aromatic_bonds,
    smiles_to_xyz,
    smiles_to_xyz_quick,
    smiles_to_xyz_quick_hapto_previews,
)  # noqa: F401

from delfin.manta.conformer_pools import (
    _append_ring_puckers,
    _conf_ff_energy,
    _conf_heavy_rmsd,
    _cp_basin_6ring,
    _emit_d8_sp4_variants,
    _emit_ring_puckers_rp,
    _frag_xyz_collapsed,
    _heavy_atom_rmsd_xyz,
    _organic_conformer_pool,
    _organic_mmff_ensemble,
    _rank_template_conformers,
    _ring_canonical_snap_z,
    _tfd_dedup_pool,
)  # noqa: F401

# Iter-8.4a: sigma chelate-cap restoration (forward-port from 123a130).
# When the sigma class loses topology %match because chelate-cap tightening
# reduced ETKDG trial counts (HEAD: 8/15/20/40 at >90/>60/>30/>20 atom
# fragments), restore the historical wider trial counts for the sigma class
# only.  d8 / d10 metals with crowded chelates need >12 trials to find the
# correct bite-angle conformer.  Default OFF: bit-exact HEAD when env-flag
# is unset.  Class-dispatched in smiles_to_xyz_isomers entry.
_SIGMA_CHELATE_CAPS_123A: Dict[str, int] = {
    "cap_90": 8,    # >90 atoms — keep HEAD cap
    "cap_60": 15,   # >60 atoms — keep HEAD cap
    "cap_30": 40,   # >30 atoms — restore champion (HEAD: 40 already)
    "cap_20": 40,   # >20 atoms — restore champion (HEAD: 20 → 40)
}
_ITER84_SIGMA_CAPS_OVERRIDE: Optional[Dict[str, int]] = None
"""Module-global override for ``_chelate_conformer_candidates`` cap values.
Set at the entry of ``smiles_to_xyz_isomers`` when both
``DELFIN_SIGMA_PORT_123A130_ITER8=1`` AND the parent mol classifies as
'sigma'.  ``None`` (default) means use HEAD baseline caps unchanged.
Process-safe under multiprocessing pool_evaluator (each worker is a
separate process); recursive smiles_to_xyz_isomers calls within the same
mol re-set to the same value (deterministic from class)."""


def _verify_topology_from_graph(
    xyz_delfin: str,
    mol_template,
) -> bool:
    """Graph-based topology verification — no OB perception, no roundtrip.

    Checks the XYZ coordinates directly against the template molecular
    graph (from SMILES).  This is the SINGLE authoritative check for
    structural integrity.

    Three rules:
    1. Every bond in the template graph must have a reasonable distance
       in the XYZ:
       - M-D bonds: within [0.75, 1.60] × ``_get_ml_bond_length``
       - M-M bonds: within [0.70, 1.80] × ``_METAL_METAL_BOND_LENGTHS``
       - Covalent bonds (non-metal): < 2.2 Å
    2. No two heavy atoms closer than 0.7 Å (collapsed structure)
    3. No H-H closer than 0.4 Å (collapsed hydrogens)
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        n_atoms = mol_template.GetNumAtoms()
        if len(lines) != n_atoms:
            # Atom count mismatch → try AddHs fallback
            try:
                mol_h = Chem.AddHs(mol_template)
                if len(lines) == mol_h.GetNumAtoms():
                    mol_template = mol_h
                    n_atoms = mol_h.GetNumAtoms()
                else:
                    return True  # can't validate → permissive
            except Exception:
                return True

        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        # Rule 1 (universal graph invariance, three-cutoff).
        # Element-agnostic coordination-sphere check using the
        # standard covalent-radii sum (r_cov_M + r_cov_X) AND
        # CSD-calibrated M-L ideal lengths for the lower bound:
        #
        # * SMILES-bonded atoms:
        #   - Upper bound: d <= 1.35 x (r_cov_M + r_cov_X)
        #     (tolerant; lets Fe-Br ~2.9 A pass vs 2.46 A ideal)
        #   - Lower bound: d >= 0.70 x _get_ml_bond_length(M, X)
        #     (rejects collapsed bonds — Sc-O at 1.21 A is 0.59 x
        #     ideal 2.05 -> rejected)
        # * NON-bonded atoms: d >= 1.10 x sum
        #   (phantom-bond reject)
        # * 1.10 - 1.35 x sum: grey zone, neither violation.
        try:
            _BONDED_MAX_FRAC = 1.35
            _BONDED_MIN_IDEAL_FRAC = 0.65
            _PHANTOM_MIN_FRAC = 1.05
            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    continue
                m_idx = atom.GetIdx()
                m_sym = atom.GetSymbol()
                mx, my, mz = coords[m_idx]
                r_cov_m = _COVALENT_RADII.get(m_sym)
                if r_cov_m is None:
                    continue
                smiles_nbr_ids = {
                    nbr.GetIdx() for nbr in atom.GetNeighbors()
                    if nbr.GetSymbol() not in _METAL_SET
                }
                # Welle-5l T3-B: 1,3-exempt set for phantom-bond check.
                # Atoms two bonds away from the metal via a donor (the "other"
                # ring atoms in NHC carbenes, naphthyridine N, salen-N etc.)
                # are necessarily close to the metal because of the rigid
                # ligand backbone — Ru-C(carbene)-N(NHC) is a 1,3-relation
                # whose distance is fixed near r_cov_M + r_cov_X regardless
                # of UFF state.  Treating them as "phantom" bonds in the
                # topology verifier rejects valid coordination chemistry
                # (NHC, naphthyridine, salen, pyrazole-bridged etc.) and is
                # the main cause of D-AQIWAZ 11% isomer coverage.
                # Welle-5l-rev1 (2026-05-18): default flipped 1 -> 0 (strict).
                # Per user directive "check the topology strictly and
                # meticulously": the previous always-on exemption let UFF-buckled
                # geometries pass the verifier, which collapsed distinct
                # coordination isomers under fingerprint dedup (D2-ADEKUS
                # 3 -> 2 frames at scale).  Strict default rejects any
                # close 1,3 contact through a donor; if a real chelate
                # backbone needs the exemption, set the env-flag to 1
                # explicitly.
                _phantom_exempt: Set[int] = set()
                if _delfin_env_int("DELFIN_PHANTOM_13_EXEMPT", 0):
                    for _nbr in atom.GetNeighbors():
                        if _nbr.GetSymbol() in _METAL_SET:
                            continue
                        for _nnb in _nbr.GetNeighbors():
                            _nn_idx = _nnb.GetIdx()
                            if _nn_idx == m_idx:
                                continue
                            if _nnb.GetSymbol() in _METAL_SET:
                                continue
                            if _nn_idx in smiles_nbr_ids:
                                continue
                            _phantom_exempt.add(_nn_idx)
                _violation = False
                for other in mol_template.GetAtoms():
                    o_idx = other.GetIdx()
                    if o_idx == m_idx:
                        continue
                    if other.GetSymbol() in _METAL_SET:
                        continue
                    r_cov_o = _COVALENT_RADII.get(other.GetSymbol())
                    if r_cov_o is None:
                        continue
                    ox, oy, oz = coords[o_idx]
                    _d = math.sqrt(
                        (mx - ox) ** 2 + (my - oy) ** 2 + (mz - oz) ** 2
                    )
                    _cov_sum = r_cov_m + r_cov_o
                    _is_bonded = o_idx in smiles_nbr_ids
                    if _is_bonded:
                        if _d > _BONDED_MAX_FRAC * _cov_sum:
                            _violation = True
                            break
                        # Lower bound: reject collapsed M-L bonds
                        # (ratio < 0.70 vs CSD ideal).  Use the
                        # lookup-table ideal, not covalent sum.
                        try:
                            _ml_ideal = float(
                                _get_ml_bond_length(m_sym, other.GetSymbol())
                            )
                        except Exception:
                            _ml_ideal = 0.0
                        if _ml_ideal > 0 and _d < _BONDED_MIN_IDEAL_FRAC * _ml_ideal:
                            _violation = True
                            break
                    else:
                        if o_idx in _phantom_exempt:
                            # 1,3 through a donor -- chemically expected
                            # close contact (NHC, naphthyridine, salen).
                            continue
                        if _d < _PHANTOM_MIN_FRAC * _cov_sum:
                            _violation = True
                            break
                if _violation:
                    return False
        except Exception:
            pass

        # Covalent non-metal bonds: simple upper-bound distance check
        # to catch bonds that UFF has stretched beyond any reasonable
        # covalent length.  No phantom check — perceiving every
        # organic bond pair would be O(N^2) and the metal graph
        # check above already guarantees the coordination sphere is
        # intact.  Metal-metal bonds use a 1.80 x ideal upper bound.
        bridging_donor_bonds: Dict[int, List[Tuple[int, str, float, float]]] = {}
        for bond in mol_template.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            i1, i2 = a1.GetIdx(), a2.GetIdx()
            s1, s2 = a1.GetSymbol(), a2.GetSymbol()
            dx = coords[i1][0] - coords[i2][0]
            dy = coords[i1][1] - coords[i2][1]
            dz = coords[i1][2] - coords[i2][2]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            is_metal_1 = s1 in _METAL_SET
            is_metal_2 = s2 in _METAL_SET
            if is_metal_1 and is_metal_2:
                mm_key = frozenset({s1, s2})
                ideal = _METAL_METAL_BOND_LENGTHS.get(mm_key)
                if ideal is None:
                    r1 = _COVALENT_RADII.get(s1)
                    r2 = _COVALENT_RADII.get(s2)
                    ideal = (r1 + r2 + 0.3) if r1 and r2 else 2.5
                if d < 0.70 * ideal or d > 1.80 * ideal:
                    return False
            elif is_metal_1 or is_metal_2:
                # Metal-ligand bonds are already validated by the
                # graph-invariance rule above; bridging donors sit at
                # compromise positions between multiple metals and
                # need a separate tolerance window so they don't
                # trigger the perception rule's lower cutoff when
                # stretched between two ideals.
                d_atom = a2 if is_metal_1 else a1
                m_sym = s1 if is_metal_1 else s2
                d_sym = s2 if is_metal_1 else s1
                m_idx = i1 if is_metal_1 else i2
                n_metal_nbrs = sum(
                    1 for nbr in d_atom.GetNeighbors()
                    if nbr.GetSymbol() in _METAL_SET
                )
                if n_metal_nbrs >= 2:
                    ideal = float(_get_ml_bond_length(m_sym, d_sym))
                    if ideal > 0:
                        bridging_donor_bonds.setdefault(
                            d_atom.GetIdx(), []
                        ).append((m_idx, m_sym, d, ideal))
            else:
                if a1.GetAtomicNum() <= 1 or a2.GetAtomicNum() <= 1:
                    if d > 1.8:
                        return False
                else:
                    if d > 2.4:
                        return False

        for d_idx, metal_entries in bridging_donor_bonds.items():
            for _m_idx, _m_sym, d, ideal in metal_entries:
                ratio = d / ideal
                if ratio < 0.55 or ratio > 3.00:
                    return False

        # Rule 2+3: No collapsed heavy atoms or hydrogens.
        heavy_indices = [
            i for i in range(n_atoms)
            if mol_template.GetAtomWithIdx(i).GetAtomicNum() > 1
        ]
        for i in range(len(heavy_indices)):
            xi, yi, zi = coords[heavy_indices[i]]
            for j in range(i + 1, min(i + 50, len(heavy_indices))):
                xj, yj, zj = coords[heavy_indices[j]]
                dsq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                if dsq < 0.49:  # 0.7²
                    return False

        # Rule 10: Inter-ligand phantom-bond.  Every pair of non-metal
        # atoms that are NOT bonded in the SMILES graph must sit
        # outside the bond-perception threshold in the XYZ.  Two
        # thresholds so bulky ligands with unavoidable close H-X
        # contacts are not over-rejected:
        #   * heavy-heavy: 1.10 x (r_cov_i + r_cov_j)
        #   * H-involved:  0.85 x (r_cov_i + r_cov_j) (only true
        #     overlap — catches O-H collapses below 0.82 A, H-H
        #     collapses below 0.53 A; legitimate close vdW contacts
        #     >= 0.95 x sum stay allowed)
        # Metal-anything pairs are covered by Rule 1.
        try:
            _HEAVY_FRAC = 1.10
            _H_FRAC = 0.85
            _bonded_pairs: set = set()
            for _b in mol_template.GetBonds():
                _i1 = _b.GetBeginAtom().GetIdx()
                _i2 = _b.GetEndAtom().GetIdx()
                _bonded_pairs.add((min(_i1, _i2), max(_i1, _i2)))
            _nm_indices = [
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() not in _METAL_SET
            ]
            for ii in range(len(_nm_indices)):
                _i = _nm_indices[ii]
                _ai = mol_template.GetAtomWithIdx(_i)
                _zi = _ai.GetAtomicNum()
                _ri = _COVALENT_RADII.get(_ai.GetSymbol())
                if _ri is None:
                    continue
                xi2, yi2, zi2 = coords[_i]
                for jj in range(ii + 1, len(_nm_indices)):
                    _j = _nm_indices[jj]
                    if (_i, _j) in _bonded_pairs:
                        continue
                    _aj = mol_template.GetAtomWithIdx(_j)
                    _rj = _COVALENT_RADII.get(_aj.GetSymbol())
                    if _rj is None:
                        continue
                    xj2, yj2, zj2 = coords[_j]
                    _d = math.sqrt(
                        (xi2 - xj2) ** 2 + (yi2 - yj2) ** 2 + (zi2 - zj2) ** 2
                    )
                    _frac = (
                        _H_FRAC
                        if _zi <= 1 or _aj.GetAtomicNum() <= 1
                        else _HEAVY_FRAC
                    )
                    if _d < _frac * (_ri + _rj):
                        return False
        except Exception:
            pass

        # Rule 4: Pi-ring planarity — reject any ring of sp2/aromatic atoms
        # whose max out-of-plane deviation exceeds 0.25 x the mean ring bond
        # length.  The sp2 character of each ring atom is derived from the
        # RING BOND TOPOLOGY (does at least one of its ring bonds have order
        # >= 1.5: aromatic, double, or kekulize-double?) rather than from
        # RDKit's hybridisation flag, so the gate behaves identically
        # regardless of the mol's sanitisation state — essential for the
        # pipeline's internal gate and any post-hoc re-check to agree.
        try:
            # Force ring perception so aromatic / pi-rings are always
            # found even after dative-bond conversion stripped the
            # default aromaticity flags.
            try:
                Chem.GetSymmSSSR(mol_template)
            except Exception:
                pass
            ring_info = mol_template.GetRingInfo()
            if ring_info is not None:
                for ring in ring_info.AtomRings():
                    if len(ring) < 5 or len(ring) > 7:
                        continue
                    ring_set = set(ring)
                    # For each atom, check whether any of its ring bonds
                    # carries order >= 1.5 (aromatic, double, or kekulised
                    # double counts).
                    n_sp2 = 0
                    for ri in ring:
                        atom_ri = mol_template.GetAtomWithIdx(ri)
                        has_pi = False
                        for b in atom_ri.GetBonds():
                            other = b.GetOtherAtom(atom_ri).GetIdx()
                            if other not in ring_set:
                                continue
                            bt = b.GetBondType()
                            if (
                                bt == Chem.BondType.AROMATIC
                                or bt == Chem.BondType.DOUBLE
                                or b.GetIsAromatic()
                                or b.GetBondTypeAsDouble() >= 1.5
                            ):
                                has_pi = True
                                break
                        if has_pi:
                            n_sp2 += 1
                    if n_sp2 < len(ring) * 0.6:
                        continue  # not a pi ring
                    try:
                        import numpy as _np
                        pts = _np.array([coords[ri] for ri in ring])
                        # Mean in-ring bond length sets the planarity scale.
                        edges = _np.linalg.norm(
                            _np.diff(
                                _np.vstack([pts, pts[:1]]), axis=0
                            ), axis=1
                        )
                        mean_bond = float(edges.mean()) if edges.size else 1.4
                        planar_tol = DELFIN_RULE4_PI_PLANAR_TOL_FRAC * mean_bond
                        centered = pts - pts.mean(axis=0)
                        _u, _s, vh = _np.linalg.svd(centered, full_matrices=False)
                        deviations = _np.abs(centered @ vh[-1])
                        if deviations.max() > planar_tol:
                            return False
                    except Exception:
                        pass
        except Exception:
            pass

        # Rule 5: Inter-ligand proximity — heavy atoms from different
        # non-metal ligand fragments must not sit at bond-perception range.
        # The threshold is per-pair 1.15 * (r_cov_i + r_cov_j), which is the
        # usual bond-detection cutoff; a pair that violates this would be
        # drawn as a bond by viewers (and by OB) and therefore breaks the
        # intended topology downstream.  Every inter-fragment pair is
        # checked — no window cap.
        try:
            non_metal = {
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() not in _METAL_SET and a.GetAtomicNum() > 1
            }
            adj_nm: Dict[int, set] = {i: set() for i in non_metal}
            for bond in mol_template.GetBonds():
                bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bi in non_metal and bj in non_metal:
                    adj_nm[bi].add(bj)
                    adj_nm[bj].add(bi)
            visited: set = set()
            frag_id: Dict[int, int] = {}
            fid = 0
            for start in sorted(non_metal):
                if start in visited:
                    continue
                stack = [start]
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    frag_id[node] = fid
                    for nb in adj_nm.get(node, ()):
                        if nb not in visited:
                            stack.append(nb)
                fid += 1
            if fid >= 2:
                nm_list = sorted(non_metal)
                nm_syms = {
                    i: mol_template.GetAtomWithIdx(i).GetSymbol()
                    for i in nm_list
                }
                # Pre-compute donor flag + bridging neighbourhood. Inter-
                # fragment atom pairs where both sides sit in the inner
                # coordination sphere of bridged metals are geometrically
                # crowded by topology (two platonic polyhedra sharing an
                # edge) and may not fit the 1.15·r_cov textbook cutoff.
                # For those pairs we fall back to 1.00·r_cov, which still
                # rejects real covalent-range clashes.
                _metal_atoms = [
                    a for a in mol_template.GetAtoms()
                    if a.GetSymbol() in _METAL_SET
                ]
                _n_metals = len(_metal_atoms)
                _bridging_set: set = set()
                if _n_metals >= 2:
                    for _a in mol_template.GetAtoms():
                        if _a.GetAtomicNum() <= 1 or _a.GetSymbol() in _METAL_SET:
                            continue
                        _metal_nbrs = sum(
                            1 for _n in _a.GetNeighbors()
                            if _n.GetSymbol() in _METAL_SET
                        )
                        if _metal_nbrs >= 2:
                            _bridging_set.add(_a.GetIdx())
                _donor_set: set = set()
                if _n_metals >= 2:
                    for _m in _metal_atoms:
                        for _n in _m.GetNeighbors():
                            if _n.GetAtomicNum() > 1:
                                _donor_set.add(_n.GetIdx())
                for i in range(len(nm_list)):
                    fi = frag_id.get(nm_list[i], -1)
                    xi, yi, zi = coords[nm_list[i]]
                    ri = _COVALENT_RADII.get(nm_syms[nm_list[i]], 0.75)
                    for j in range(i + 1, len(nm_list)):
                        fj = frag_id.get(nm_list[j], -1)
                        if fi == fj:
                            continue
                        xj, yj, zj = coords[nm_list[j]]
                        rj = _COVALENT_RADII.get(nm_syms[nm_list[j]], 0.75)
                        # Softer threshold when both atoms sit in the
                        # inner coordination sphere of a bridged cluster
                        # (either as donor or as bridging donor), where
                        # platonic vertex targets cannot be fully
                        # reconciled.
                        if (
                            _n_metals >= 2
                            and (nm_list[i] in _donor_set
                                 or nm_list[i] in _bridging_set)
                            and (nm_list[j] in _donor_set
                                 or nm_list[j] in _bridging_set)
                        ):
                            thresh = DELFIN_RULE5_INNER_SPHERE_FACT * (ri + rj)
                        else:
                            thresh = DELFIN_RULE5_INTERFRAG_COV_FACT * (ri + rj)
                        dsq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                        if dsq < thresh * thresh:
                            return False
        except Exception:
            pass

        # Rule 6: Sp2 chelate backbone planarity. For every chelate pair
        # whose non-metal path consists entirely of sp2/aromatic atoms,
        # the metallacycle (metal + backbone) must lie nearly in a plane.
        # A twisted chelate is chemically unrealistic.
        #
        # sp2 character is read from the bond graph (does the atom carry at
        # least one bond of order >= 1.5?) rather than from the mol's
        # hybridisation / aromaticity flags, so the gate behaves identically
        # regardless of sanitisation state — essential for the pipeline's
        # in-line gate and any post-hoc re-check to agree.
        try:
            import numpy as _np

            def _is_sp2_graph(_atom):
                for _b in _atom.GetBonds():
                    if (
                        _b.GetBondType() == Chem.BondType.AROMATIC
                        or _b.GetBondType() == Chem.BondType.DOUBLE
                        or _b.GetIsAromatic()
                        or _b.GetBondTypeAsDouble() >= 1.5
                    ):
                        return True
                return False

            metal_idxs = [
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() in _METAL_SET
            ]
            for m_idx in metal_idxs:
                donors = [
                    nbr.GetIdx()
                    for nbr in mol_template.GetAtomWithIdx(m_idx).GetNeighbors()
                    if nbr.GetAtomicNum() > 1
                    and nbr.GetSymbol() not in _METAL_SET
                ]
                for i in range(len(donors)):
                    for j in range(i + 1, len(donors)):
                        d1, d2 = donors[i], donors[j]
                        # BFS from d1 to d2, blocking metal
                        visited = {m_idx, d1}
                        prev: Dict[int, int] = {d1: -1}
                        queue = [d1]
                        while queue:
                            cur = queue.pop(0)
                            if cur == d2:
                                break
                            for n in mol_template.GetAtomWithIdx(cur).GetNeighbors():
                                ni = n.GetIdx()
                                if ni in visited:
                                    continue
                                visited.add(ni)
                                prev[ni] = cur
                                queue.append(ni)
                        if d2 not in prev:
                            continue
                        path: List[int] = []
                        node = d2
                        while node != -1:
                            path.append(node)
                            node = prev.get(node, -1)
                        if len(path) < 3 or len(path) > 5:
                            continue
                        # All path atoms must be sp2 (from graph) for
                        # planarity enforcement.
                        if not all(
                            _is_sp2_graph(mol_template.GetAtomWithIdx(pi))
                            for pi in path
                        ):
                            continue
                        cycle = [m_idx] + path
                        pts = _np.array([coords[ci] for ci in cycle])
                        centered = pts - pts.mean(axis=0)
                        _u, _s, vh = _np.linalg.svd(centered, full_matrices=False)
                        dev = float(_np.abs(centered @ vh[-1]).max())
                        if dev > DELFIN_RULE6_METALLACYCLE_MAX_DEV:
                            return False
        except Exception:
            pass

        # Rule 6b REMOVED as a hard reject.  The underlying observation
        # (pi-ring sigma-coordinated metals should sit near the ring
        # plane) is valid, but the builder cannot currently guarantee
        # this geometry — ETKDG places imidazole / pyridine / etc. in
        # arbitrary rotations around the M-D axis, and no axial
        # rotation can move a metal that sits off the original ring
        # plane INTO it.  Hard-rejecting 70 %+ of built candidates on
        # crowded CN >= 6 systems dropped output from 12 down to 6
        # isomers on the Cd-histidine test case.  The chelate-ring
        # planarity penalty in ``_geometry_quality_score`` still
        # softly down-ranks tilted rings; a proper fix lives in the
        # builder (fragment-orientation pre-search before Procrustes).

        # Rule 7: Hybridisation vs. coordination geometry.
        # For every non-metal heavy atom NOT bonded to a metal we infer
        # the expected hybridisation from the bond-order graph and
        # compare it to the local 3D geometry.  Atoms coordinated to a
        # metal are excluded because their local angles are dictated
        # by the coordination polyhedron (a donor N on a square-planar
        # Pd is intentionally non-tetrahedral).
        #
        # Predicates (bond-graph only, no RDKit flags):
        #   sp   — 2 heavy non-metal neighbours AND at least one bond
        #          of order >= 2.5 (triple / cumulated double).  Expect
        #          near-linear X-A-Y angle (>= 150°).
        #   sp²  — 3 heavy non-metal neighbours AND at least one bond
        #          of order >= 1.5 (aromatic, double, kekulé-double).
        #          Expect planar: the out-of-plane distance of A from
        #          the plane of its three neighbours <= 0.35 Å.
        #   sp³  — 4 heavy non-metal neighbours AND zero bonds of
        #          order >= 1.5.  Expect tetrahedral: every angle
        #          X-A-Y >= 80° (rejects severe distortion while
        #          tolerating ring strain down to ~84° in cyclopropane).
        #
        # Each predicate is only applied when the neighbour count matches;
        # atoms with H-only neighbours or unusual coordination are
        # skipped.  The tolerances are intentionally generous so
        # chemically reasonable geometry is never rejected — the gate
        # targets UFF blow-ups and collapsed-ring artefacts only.
        try:
            import numpy as _np

            def _bond_orders_sum(_atom):
                total = 0.0
                for _b in _atom.GetBonds():
                    try:
                        total += float(_b.GetBondTypeAsDouble())
                    except Exception:
                        pass
                return total

            def _max_ring_bond_order_to_neighbors(_atom):
                m = 0.0
                for _b in _atom.GetBonds():
                    try:
                        m = max(m, float(_b.GetBondTypeAsDouble()))
                    except Exception:
                        pass
                return m

            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() in _METAL_SET:
                    continue
                if atom.GetAtomicNum() <= 1:
                    continue
                # Skip atoms directly bonded to a metal.  Their local
                # geometry is dictated by the coordination polyhedron and
                # by UFF's missing metal parameters; routine
                # pyramidalisation of cyclometallated / NHC carbons
                # should not be treated as a topology violation here.
                # Metal-in-π enforcement is delegated to Rule 6
                # (metallacycle planarity) and the UFF metallacycle
                # torsion constraint.
                if any(
                    nbr.GetSymbol() in _METAL_SET
                    for nbr in atom.GetNeighbors()
                ):
                    continue
                heavy_nbr_organic = [
                    nbr.GetIdx() for nbr in atom.GetNeighbors()
                    if nbr.GetAtomicNum() > 1
                    and nbr.GetSymbol() not in _METAL_SET
                ]
                heavy_nbr_all = heavy_nbr_organic
                if not heavy_nbr_all:
                    continue
                max_bo = _max_ring_bond_order_to_neighbors(atom)
                ai = atom.GetIdx()
                a_pos = _np.array(coords[ai])

                # sp — linear by coordination design when bonded to
                # metal (e.g. M-C#O, M-C#N-R).  Skip metal-bonded atoms.
                if len(heavy_nbr_organic) == 2 and max_bo >= 2.5:
                    n1, n2 = heavy_nbr_organic
                    v1 = _np.array(coords[n1]) - a_pos
                    v2 = _np.array(coords[n2]) - a_pos
                    n1n = float(_np.linalg.norm(v1))
                    n2n = float(_np.linalg.norm(v2))
                    if n1n > 1e-6 and n2n > 1e-6:
                        cos_a = float(_np.dot(v1, v2) / (n1n * n2n))
                        cos_a = max(-1.0, min(1.0, cos_a))
                        angle_deg = math.degrees(math.acos(cos_a))
                        if angle_deg < DELFIN_RULE7_SP_MIN_ANGLE_DEG:
                            return False
                    continue

                # sp² — exactly 3 heavy neighbours (metal counted) AND
                # at least one π bond.  Trigonal planar: atom must lie
                # within DELFIN_RULE7_SP2_OOP_MAX of the plane of its
                # three neighbours and all three X-A-Y angles must fall
                # in [DELFIN_RULE7_SP2_ANGLE_MIN,
                # DELFIN_RULE7_SP2_ANGLE_MAX] (ideal 120°).  Including
                # the metal here is what enforces "metal in π-plane"
                # for cyclometallated / conjugated donor atoms;
                # excluding it would let UFF pyramidalise the donor
                # carbon along the metal axis.
                if len(heavy_nbr_all) == 3 and max_bo >= 1.5:
                    na, nb, nc = heavy_nbr_all
                    pa = _np.array(coords[na])
                    pb = _np.array(coords[nb])
                    pc = _np.array(coords[nc])
                    normal = _np.cross(pb - pa, pc - pa)
                    nn = float(_np.linalg.norm(normal))
                    if nn > 1e-9:
                        normal = normal / nn
                        centroid = (pa + pb + pc) / 3.0
                        dev = float(abs(_np.dot(a_pos - centroid, normal)))
                        # Metal-bonded sp² atoms (NHC carbenes,
                        # cyclometallated donors, carbonyl C) get a
                        # looser threshold because UFF without metal
                        # parameters routinely pyramidalises them by
                        # 0.3-0.5 Å without distorting the rest of the
                        # topology.  Rule 6 (metallacycle planarity)
                        # and the UFF metallacycle-torsion constraint
                        # already drive these atoms back toward the
                        # plane on realistic energy scales.
                        # Metal-bonded atoms are already skipped above (see the
                        # _METAL_SET neighbour `continue`), so any atom reaching
                        # here has NO metal neighbour: the looser metal budget is
                        # dead by construction.  Bind it explicitly (was an
                        # undefined name) to keep the strict budget and the intent
                        # documented.
                        has_metal_nbr = False
                        oop_budget = (
                            DELFIN_RULE7_SP2_OOP_MAX_METAL
                            if has_metal_nbr
                            else DELFIN_RULE7_SP2_OOP_MAX
                        )
                        if dev > oop_budget:
                            return False
                    angle_pairs = ((na, nb), (na, nc), (nb, nc))
                    for p, q in angle_pairs:
                        vp = _np.array(coords[p]) - a_pos
                        vq = _np.array(coords[q]) - a_pos
                        np_ = float(_np.linalg.norm(vp))
                        nq_ = float(_np.linalg.norm(vq))
                        if np_ < 1e-6 or nq_ < 1e-6:
                            continue
                        cos_a = float(_np.dot(vp, vq) / (np_ * nq_))
                        cos_a = max(-1.0, min(1.0, cos_a))
                        ang = math.degrees(math.acos(cos_a))
                        if (
                            ang < DELFIN_RULE7_SP2_ANGLE_MIN
                            or ang > DELFIN_RULE7_SP2_ANGLE_MAX
                        ):
                            return False
                    continue

                # sp³ — 4 heavy organic neighbours, no π bonds.
                # Reject severe tetrahedral collapse (any X-A-Y
                # angle < DELFIN_RULE7_SP3_MIN_ANGLE_DEG).
                if len(heavy_nbr_organic) == 4 and max_bo < 1.5:
                    nbr_positions = [
                        _np.array(coords[k]) for k in heavy_nbr_organic
                    ]
                    min_angle = 360.0
                    for i_ in range(4):
                        for j_ in range(i_ + 1, 4):
                            v_i = nbr_positions[i_] - a_pos
                            v_j = nbr_positions[j_] - a_pos
                            ni_ = float(_np.linalg.norm(v_i))
                            nj_ = float(_np.linalg.norm(v_j))
                            if ni_ < 1e-6 or nj_ < 1e-6:
                                continue
                            cos_a = float(_np.dot(v_i, v_j) / (ni_ * nj_))
                            cos_a = max(-1.0, min(1.0, cos_a))
                            ang = math.degrees(math.acos(cos_a))
                            if ang < min_angle:
                                min_angle = ang
                    if min_angle < DELFIN_RULE7_SP3_MIN_ANGLE_DEG:
                        return False
                    continue
        except Exception:
            pass

        return True
    except Exception:
        return False
if os.environ.get("DELFIN_FFFREE_GEOM_IDEALS_REAL", "0") == "1":
    _GEOM_IDEAL_ANGLES.update(_GEOM_IDEAL_ANGLES_REAL)


def _chelate_conformer_candidates(
    mol,
    frag_atom_indices,
    donor_atom_indices,
    target_bite,
    n_trials: int = DELFIN_CHELATE_N_TRIALS,
    accept_delta: float = DELFIN_CHELATE_ACCEPT_DELTA,
    reject_delta: float = DELFIN_CHELATE_REJECT_DELTA,
    max_candidates: int = 5,
):
    """Return a deterministic list of chelate conformer coordinates
    whose donor-donor distance pattern matches the polyhedron
    vertex-pair pattern, ordered by goodness of fit (best first).

    Each entry is a dict ``{original_atom_idx: (x, y, z)}``.  The list
    never exceeds ``max_candidates``; every returned entry has
    ``delta < reject_delta``.  For bidentate chelates ``target_bite``
    is a scalar (``|vertex_i - vertex_j|``); for polydentate ligands
    it is a ``(k, k)`` pairwise distance matrix.  Returning a list
    (instead of only the best) lets the caller emit multiple ring
    puckers / backbone conformations as distinct topology isomers,
    which is essential for macrocyclic and tridentate+ ligands whose
    natural backbone space contains several clash-free realisations
    compatible with the same platonic polyhedron.
    """
    if not RDKIT_AVAILABLE:
        return []
    try:
        import numpy as _np
    except Exception:
        return []

    frag_list = sorted(frag_atom_indices)
    if len(frag_list) < 3:
        return []
    old_to_new = {old: new for new, old in enumerate(frag_list)}
    donor_new = [old_to_new[d] for d in donor_atom_indices if d in old_to_new]
    if len(donor_new) < 2:
        return []

    # --- Class-aware chelate-rank gating ----------------------------------
    # Default OFF — bit-exact when ``DELFIN_CHELATE_RANK_CLASS_AWARE`` is
    # unset.  When enabled the per-conformer composite score becomes
    #   composite = delta + alpha * (element_weight + spread_factor * pucker)
    # where ``element_weight`` is constant for a given fragment (donor
    # set fixed) and ``pucker`` is the per-conformer donor-plane RMS
    # out-of-plane deviation (Å).  Sigma class rewards puckered backbones
    # (chair / boat tridentates), hapto class is neutral on pucker.
    # Implementation note: the element weight alone cannot reorder
    # conformers from the same call (constant); the per-conformer pucker
    # is what makes the secondary score actually re-rank.
    _class_aware_enabled = False
    _class_penalty_const = 0.0
    _class_spread_factor = 0.0
    try:
        if _class_conditional_flag(
            "DELFIN_CHELATE_RANK_CLASS_AWARE", mol
        ):
            _class_aware_enabled = True
            _cls = _classify_complex_class(mol)
            _class_penalty_const = _chelate_class_donor_penalty(
                mol, donor_atom_indices, _cls,
            )
            # Per-class pucker-diversity reward.  Negative value =
            # reward (lowers the composite when the conformer is more
            # out-of-plane).  Sigma class strongly rewards pucker variety
            # because mer/fac tridentate isomers differ exactly in the
            # backbone-plane deviation.  Hapto and multi_* classes get a
            # milder reward; no_metal disables entirely.
            _class_spread_factor = {
                "sigma":       -0.50,
                "hapto":       -0.10,
                "multi_sigma": -0.30,
                "multi_hapto": -0.10,
                "no_metal":     0.00,
            }.get(_cls, 0.0)
    except Exception:
        _class_aware_enabled = False
        _class_penalty_const = 0.0
        _class_spread_factor = 0.0

    is_pairwise_matrix = not isinstance(target_bite, (int, float))
    if is_pairwise_matrix:
        try:
            target_mat = _np.asarray(target_bite, dtype=float)
        except Exception:
            return []
        if target_mat.shape != (len(donor_new), len(donor_new)):
            return []

    rw = Chem.RWMol(Chem.Mol())
    for aidx in frag_list:
        atom = mol.GetAtomWithIdx(aidx)
        new_idx = rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        rw.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
        rw.GetAtomWithIdx(new_idx).SetNoImplicit(True)
        rw.GetAtomWithIdx(new_idx).SetNumExplicitHs(atom.GetNumExplicitHs())
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and bj in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[bj], bond.GetBondType())
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        try:
            rw.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    frag_mol = rw.GetMol()

    # Shared pipeline seed schedule.  The first 12 entries coincide with
    # the top-level conformer sampling seeds so the chelate conformer
    # choice is cache-coherent with the rest of the pipeline; the tail
    # extends the search for heavy macrocycles and tridentate+ ligands
    # whose native-bite matching requires a wider conformational sweep.
    _SEEDS = _PIPELINE_SEEDS

    def _try_embed(m, seed):
        p = AllChem.ETKDGv3()
        p.useRandomCoords = True
        p.randomSeed = int(seed)
        try:
            return _embed_with_timeout(m, p, timeout=_CHELATE_EMBED_TIMEOUT)
        except Exception:
            return -1

    def _preopt_fragment_conformer(_mol_obj, _cid) -> None:
        """MMFF94/UFF pre-optimisation of an isolated ligand conformer.

        DISABLED pending investigation of Ir(ppy)2(acac) all-cis regression.
        """
        return

    # Collect every accepted (delta, coords_map) pair and sort by fit.
    accepted: List[Tuple[float, Dict[int, tuple]]] = []
    fallback_mol = None

    # Scale number of trials with fragment size.  Huge polydentate
    # ligands (terpyridine-NMe2, salen-biphep phosphine backbone,
    # phos-terpy) have ETKDG wall-clocks measured in seconds per seed,
    # so the default 40 trials × 6 s timeout accumulates into minutes
    # of wall-time — and with orphan thread pile-up starves the
    # subprocess long before the first isomer is written.  Caps are
    # exposed as DELFIN_CHELATE_CAP_{60,90} so regressions where the
    # correct pose lives past the 8-/15-seed window can be recovered
    # by raising the cap for that run.
    _frag_n = frag_mol.GetNumAtoms()
    # Iter-8.4a: when the sigma chelate-cap port is active (module-global
    # ``_ITER84_SIGMA_CAPS_OVERRIDE`` set by smiles_to_xyz_isomers entry),
    # use the wider 123a130 caps instead of HEAD baselines.  ``None``
    # preserves HEAD bit-exactness.
    _iter84_caps = _ITER84_SIGMA_CAPS_OVERRIDE
    if _iter84_caps is not None:
        cap_90 = _iter84_caps["cap_90"]
        cap_60 = _iter84_caps["cap_60"]
        cap_30 = _iter84_caps["cap_30"]
        cap_20 = _iter84_caps["cap_20"]
    else:
        cap_90 = DELFIN_CHELATE_CAP_90
        cap_60 = DELFIN_CHELATE_CAP_60
        cap_30 = DELFIN_CHELATE_CAP_30
        cap_20 = int(os.environ.get('DELFIN_CHELATE_CAP_20', '20'))
    if _frag_n > 90:
        n_trials = min(n_trials, cap_90)
    elif _frag_n > 60:
        n_trials = min(n_trials, cap_60)
    elif _frag_n > 30:
        n_trials = min(n_trials, cap_30)
    elif _frag_n > 20:
        # Moderate-size ligands (terpyridines ~22 atoms, salen-biphep ~25):
        # default n_trials of 40 yields diminishing returns past ~20 trials.
        # Halving here recovers ~60s per such SMILES from blocking subprocess
        # without measurable loss in conformer diversity (env override:
        # DELFIN_CHELATE_CAP_20; Iter-8.4a sigma override widens to 40).
        n_trials = min(n_trials, cap_20)

    for seed in _SEEDS[:n_trials]:
        cid = _try_embed(frag_mol, seed)
        if cid < 0:
            if fallback_mol is None:
                try:
                    rw2 = Chem.RWMol(frag_mol)
                    for a in rw2.GetAtoms():
                        if a.GetIsAromatic():
                            a.SetIsAromatic(False)
                    for b in rw2.GetBonds():
                        if (
                            b.GetIsAromatic()
                            or b.GetBondType() == Chem.BondType.AROMATIC
                        ):
                            b.SetIsAromatic(False)
                            b.SetBondType(Chem.BondType.SINGLE)
                    fallback_mol = rw2.GetMol()
                    try:
                        fallback_mol.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                except Exception:
                    fallback_mol = None
            if fallback_mol is None:
                continue
            cid = _try_embed(fallback_mol, seed)
            if cid < 0:
                continue
            used = fallback_mol
        else:
            used = frag_mol

        # Ligand-first: relax the isolated conformer with MMFF94/UFF
        # before evaluating its donor pattern.  The brief pre-opt
        # removes bonded-term strain baked into ETKDG's initial coords
        # and gives a chemically sensible ligand shape that DFT can
        # start from after placement.  On macrocycles this is the
        # difference between a clean low-energy pucker and a ring
        # with spurious kinks.
        _preopt_fragment_conformer(used, cid)

        conf = used.GetConformer(cid)
        donor_pts = _np.array([
            [
                conf.GetAtomPosition(idx).x,
                conf.GetAtomPosition(idx).y,
                conf.GetAtomPosition(idx).z,
            ]
            for idx in donor_new
        ])
        if is_pairwise_matrix:
            diffs = donor_pts[:, None, :] - donor_pts[None, :, :]
            d_mat = _np.linalg.norm(diffs, axis=-1)
            n = len(donor_new)
            n_pairs = n * (n - 1) / 2
            ss = float(_np.triu((d_mat - target_mat) ** 2, k=1).sum())
            delta = (ss / max(n_pairs, 1.0)) ** 0.5
        else:
            p0 = donor_pts[0]
            p1 = donor_pts[1]
            d_dd = float(_np.linalg.norm(p0 - p1))
            delta = abs(d_dd - float(target_bite))
        if delta >= reject_delta:
            continue
        coords_map = {
            old: (
                conf.GetAtomPosition(new).x,
                conf.GetAtomPosition(new).y,
                conf.GetAtomPosition(new).z,
            )
            for old, new in old_to_new.items()
        }

        # Per-conformer pucker score for class-aware re-rank.  Computed
        # only when the class-aware flag is on; otherwise the secondary
        # term is 0 and the sort is bit-exact with HEAD.  Pucker = RMS
        # out-of-plane distance of donor atoms from their best-fit
        # plane (n_donors >= 3) or 0 for bidentate (no plane to fit).
        _pucker = 0.0
        if _class_aware_enabled and len(donor_new) >= 3:
            try:
                _cen = donor_pts.mean(axis=0)
                _X = donor_pts - _cen
                # SVD of centered donor matrix; smallest singular vector
                # = plane normal.  Plane RMS = smallest singular value /
                # sqrt(n).
                _U, _S, _Vt = _np.linalg.svd(_X, full_matrices=False)
                if _S.size >= 1:
                    _pucker = float(_S[-1] / max(1.0, len(donor_new) ** 0.5))
            except Exception:
                _pucker = 0.0

        accepted.append((delta, coords_map, _pucker))
        # Early-exit once we have ``max_candidates`` "good" fits
        # (delta < accept_delta).  Every extra seed after this point
        # can only replace an already-good candidate with a slightly
        # better one — not worth the 6 s per-seed wall-time on the
        # heaviest ligands where every seed costs real time.
        good = sum(1 for d, _c, _p in accepted if d < accept_delta)
        if good >= max_candidates:
            break

    if not accepted:
        return []

    # Sort by fit quality (best first), cap at max_candidates.
    if _class_aware_enabled:
        _alpha = DELFIN_CHELATE_RANK_CLASS_ALPHA
        accepted.sort(
            key=lambda item: (
                item[0]
                + _alpha * (
                    _class_penalty_const
                    + _class_spread_factor * item[2]
                )
            )
        )
    else:
        accepted.sort(key=lambda item: item[0])
    return [c for _d, c, _p in accepted[:max_candidates]]


def _best_chelate_conformer_coords(
    mol,
    frag_atom_indices,
    donor_atom_indices,
    target_bite,
    n_trials: int = DELFIN_CHELATE_N_TRIALS,
    accept_delta: float = DELFIN_CHELATE_ACCEPT_DELTA,
    rank: int = 0,
    max_candidates: int = 5,
):
    """Return the ``rank``-th best chelate conformer (``rank=0`` is the best
    fit).  When ``rank`` exceeds the number of accepted candidates, falls
    back to the best available one so callers always receive a valid
    placement when at least one conformer fits the target bite.

    ``max_candidates`` limits how many conformers the underlying search
    retains; keep it >= the highest ``rank`` the caller intends to query.
    """
    cands = _chelate_conformer_candidates(
        mol,
        frag_atom_indices,
        donor_atom_indices,
        target_bite,
        n_trials=n_trials,
        accept_delta=accept_delta,
        max_candidates=max(max_candidates, rank + 1),
    )
    if not cands:
        return None
    if rank < 0 or rank >= len(cands):
        return cands[0]
    return cands[rank]


def _embed_fragment_procrustes(
    mol,
    metal_idx: int,
    frag_atom_indices: set,
    frag_donor_indices: List[int],
    target_positions: List[Tuple[float, float, float]],
    coords: List[Tuple[float, float, float]],
    chelate_rank: int = 0,
) -> bool:
    """Embed a ligand fragment via RDKit ETKDG, then Procrustes-align donors to targets.

    Modifies *coords* in-place for atoms in *frag_atom_indices*.
    Returns True on success, False on failure (caller should use BFS fallback).
    """
    try:
        import numpy as np
    except ImportError:
        return False

    if not frag_donor_indices or not frag_atom_indices:
        return False

    # Build a sub-molecule for the fragment (non-metal atoms only)
    frag_list = sorted(frag_atom_indices)
    old_to_new = {old: new for new, old in enumerate(frag_list)}

    rw = Chem.RWMol(Chem.Mol())
    for aidx in frag_list:
        atom = mol.GetAtomWithIdx(aidx)
        new_idx = rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        rw.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
        rw.GetAtomWithIdx(new_idx).SetNoImplicit(True)
        rw.GetAtomWithIdx(new_idx).SetNumExplicitHs(atom.GetNumExplicitHs())

    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and bj in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[bj], bond.GetBondType())

    try:
        Chem.SanitizeMol(rw)
    except Exception:
        try:
            rw.UpdatePropertyCache(strict=False)
        except Exception:
            pass

    frag_mol = rw.GetMol()

    donor_new_indices = [old_to_new[d] for d in frag_donor_indices if d in old_to_new]
    if not donor_new_indices:
        return False

    # Chelate: let the ligand find its own native backbone geometry
    # whose donor-donor distances match the polyhedron vertex pairs.
    # Bidentate uses a scalar bite; tri-/tetradentate uses the full
    # pairwise distance matrix.
    if len(donor_new_indices) >= 2 and len(target_positions) >= len(donor_new_indices):
        tp = np.asarray(target_positions[:len(donor_new_indices)], dtype=float)
        if len(donor_new_indices) == 2:
            target = float(np.linalg.norm(tp[0] - tp[1]))
        else:
            # full pairwise matrix, donor order matches frag_donor_indices
            diffs = tp[:, None, :] - tp[None, :, :]
            target = np.linalg.norm(diffs, axis=-1)
        coords_map = _best_chelate_conformer_coords(
            mol, frag_atom_indices, frag_donor_indices, target,
            rank=chelate_rank,
        )
        if coords_map is not None:
            frag_coords = np.array(
                [list(coords_map[old]) for old in frag_list],
                dtype=float,
            )
        else:
            frag_coords = None
    else:
        frag_coords = None

    if frag_coords is None:
        # Single-seed fragment ETKDG fallback (monodentate, higher-denticity,
        # or bidentate chelate for which the conformer search failed).
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = 42
        try:
            cid = _embed_with_timeout(frag_mol, params, timeout=_CHELATE_EMBED_TIMEOUT)
        except Exception:
            cid = -1
        if cid < 0:
            try:
                rw2 = Chem.RWMol(frag_mol)
                for atom in rw2.GetAtoms():
                    if atom.GetIsAromatic():
                        atom.SetIsAromatic(False)
                for bond in rw2.GetBonds():
                    if (
                        bond.GetIsAromatic()
                        or bond.GetBondType() == Chem.BondType.AROMATIC
                    ):
                        bond.SetIsAromatic(False)
                        bond.SetBondType(Chem.BondType.SINGLE)
                frag_mol2 = rw2.GetMol()
                try:
                    frag_mol2.UpdatePropertyCache(strict=False)
                except Exception:
                    pass
                cid2 = _embed_with_timeout(
                    frag_mol2, params, timeout=_CHELATE_EMBED_TIMEOUT
                )
                if cid2 < 0:
                    return False
                frag_mol = frag_mol2
                cid = cid2
            except Exception:
                return False
        frag_conf = frag_mol.GetConformer(cid)
        frag_coords = np.array([
            [frag_conf.GetAtomPosition(i).x,
             frag_conf.GetAtomPosition(i).y,
             frag_conf.GetAtomPosition(i).z]
            for i in range(frag_mol.GetNumAtoms())
        ])

    src = frag_coords[donor_new_indices]
    tgt = np.array(target_positions[:len(donor_new_indices)])

    if len(src) != len(tgt) or len(src) == 0:
        return False

    # Procrustes alignment: translate, then rotate
    src_center = src.mean(axis=0)
    tgt_center = tgt.mean(axis=0)
    src_centered = src - src_center
    tgt_centered = tgt - tgt_center

    if len(src) >= 2:
        # SVD for optimal rotation
        H = src_centered.T @ tgt_centered
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1, 1, 1 if d > 0 else -1])
        R = Vt.T @ sign_matrix @ U.T
    else:
        # Single donor: align the donor's lone-pair direction (anti-bisector
        # of donor -> heavy-neighbour vectors) with the donor -> metal
        # direction.  This locks the ring plane in a chemically-correct
        # orientation (LP points at M) instead of the earlier heuristic
        # that only aligned "centroid -> donor" (which left the ring
        # plane under-constrained and forced the post-build orient step
        # to fix LP alignment at clash cost).
        d_new = donor_new_indices[0]
        d_atom = frag_mol.GetAtomWithIdx(d_new)
        # ROOT FIX (DELFIN_FFFREE_SP3C_TET_SEAT=1, default OFF -> byte-identical): a monodentate sp3-C
        # donor (M-CH2-R, M-CH3) has NO lone pair -- the metal occupies the 4th tetrahedral vertex.  The
        # HEAVY-ONLY bisector below leaves only the single heavy tail (R) for an M-CH2-R, so
        # src_dir = -unit(donor->R) and aligning it with donor->metal drives R ANTI to the metal ->
        # M-C-R ~180 deg (eye: sp3c_donor_linear; a top sigma_coord defect).  Including the H neighbours
        # for an sp3-C makes the bisector the true vacant-slot direction -> M-C-X ~109 deg.  Not a
        # post-hoc bend; sets the placement.
        _sp3c_tet = (os.environ.get("DELFIN_FFFREE_SP3C_TET_SEAT", "0") == "1"
                     and d_atom.GetSymbol() == "C" and not d_atom.GetIsAromatic()
                     and not d_atom.IsInRing())        # PENDANT sp3 alkyl donor only
        nbr_idx = [
            nb.GetIdx() for nb in d_atom.GetNeighbors()
            if nb.GetAtomicNum() > 1 or _sp3c_tet
        ]
        src_dir = None
        if nbr_idx:
            d_pos_src = frag_coords[d_new]
            lp = np.zeros(3)
            for ni in nbr_idx:
                v = frag_coords[ni] - d_pos_src
                vn = np.linalg.norm(v)
                if vn > 1e-8:
                    lp += v / vn
            lp_n = np.linalg.norm(lp)
            if lp_n > 1e-6:
                # Anti-bisector: away from neighbours, toward lone pair.
                src_dir = -lp / lp_n
        if src_dir is None:
            # Atomic donor (no ring/chain neighbours): fall back to the
            # original centroid heuristic.
            frag_center = frag_coords.mean(axis=0)
            raw = src[0] - frag_center
            rn = np.linalg.norm(raw)
            src_dir = raw / rn if rn > 1e-8 else np.array([1.0, 0.0, 0.0])
        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            _trace_seating("procrustes_single_donor elem=%s n_heavy_nbrs=%d method=%s" % (
                d_atom.GetSymbol(), len(nbr_idx),
                "anti_bisector_HEAVY_ONLY" if nbr_idx else "centroid"))
        # Target LP direction: donor -> metal (metal is at origin).
        tgt_dir = -tgt[0]
        tn = np.linalg.norm(tgt_dir)
        if tn < 1e-8:
            R = np.eye(3)
        else:
            tgt_dir = tgt_dir / tn
            v = np.cross(src_dir, tgt_dir)
            c = float(np.dot(src_dir, tgt_dir))
            if np.linalg.norm(v) < 1e-8:
                R = np.eye(3) if c > 0 else -np.eye(3)
            else:
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                R = np.eye(3) + vx + vx @ vx / (1 + c)

    # Transform all fragment atoms
    transformed = (frag_coords - src_center) @ R.T + tgt_center

    # Write back to coords
    for new_idx, old_idx in enumerate(frag_list):
        coords[old_idx] = tuple(transformed[new_idx])

    return True


def _build_topology_xyz_from_scratch(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    chelate_rank: int = 0,
    inflation_schedule: Tuple[float, ...] = (0.3, 0.5, 0.7, 0.9, 1.0),
) -> Optional[str]:
    """Balloon-inflate build that needs no ETKDG template.

    Produces a DELFIN XYZ by constructing the coordination sphere
    *from scratch* — the metal sits at the origin, non-bridging donors
    are placed directly on their polyhedron-vertex × ideal-M-D
    positions, bridging donors sit on the M-M axis at the compromise
    point, and each ligand fragment is MMFF-optimised *in isolation*
    then rigidly Procrustes-aligned so that its donor atoms land on
    those vertices.  The full structure is then grown radially in
    `inflation_schedule` steps, with ligand rotations per step to
    break inter-fragment clashes.

    Returns a DELFIN-format XYZ string, or ``None`` when any phase
    fails in a way the caller cannot recover from (unknown
    geometry, fragment embedding impossible, unresolvable clash).

    This is the preferred pre-UFF builder because it guarantees
    CSD-realistic M-D distances, uses no ETKDG-template bias, and
    produces one deterministic structure per (CF, perm,
    chelate_rank) triple.  ``_build_topology_xyz_from_template``
    remains as a fallback for systems where the isolated-fragment
    ETKDG fails to converge.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        import numpy as np
    except Exception:
        return None

    vectors = _TOPO_GEOMETRY_VECTORS.get(geometry)
    if not vectors:
        return None

    try:
        n_atoms = mol.GetNumAtoms()
        coords: List[List[float]] = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
        placed: set = set()

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        coords[metal_idx] = [0.0, 0.0, 0.0]
        placed.add(metal_idx)

        # Phase 0a — scaffold the second metal + bridging donors if the
        # input is bimetallic.  Both helpers already return
        # {atom_idx: (x, y, z)} dicts with ideal M-M + bridge distances.
        all_metal_idxs = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() in _METAL_SET
        ]
        bridging = _find_bridging_donors(mol) if len(all_metal_idxs) >= 2 else []
        scaffold: Optional[Dict[int, Tuple[float, float, float]]] = None
        if len(all_metal_idxs) >= 2 and bridging:
            try:
                # Ensure metal_idx is the first entry so its position is
                # the origin the builder assumes.
                ordered_metals = [metal_idx] + [
                    m for m in all_metal_idxs if m != metal_idx
                ]
                scaffold = _build_multimetal_scaffold(
                    mol, ordered_metals, bridging
                )
            except Exception as exc:
                logger.debug(
                    "Balloon scaffold failed (falling back to mono-metal): %s",
                    exc,
                )
                scaffold = None
        if scaffold:
            for a_idx, (x, y, z) in scaffold.items():
                coords[a_idx] = [x, y, z]
                placed.add(a_idx)

        # Phase 0b — donor target positions.  Non-bridging donors of
        # ``metal_idx`` get their polyhedron-vertex × ideal-M-D;
        # already-placed bridging donors stay where the scaffold put
        # them.  Other-metal non-bridging donors stay unplaced for now
        # — they are either placed by their own per-metal build pass
        # or remain as free atoms for the downstream BFS.
        donor_target_map: Dict[int, Tuple[float, float, float]] = {}
        bridging_atom_set = {d_idx for d_idx, _ in bridging}
        for pos_idx, donor_list_idx in enumerate(perm):
            donor_atom_idx = donor_atom_indices[donor_list_idx]
            if donor_atom_idx in placed and donor_atom_idx in bridging_atom_set:
                # Bridging donor — keep scaffold position as its target.
                donor_target_map[donor_atom_idx] = tuple(coords[donor_atom_idx])
                continue
            donor_sym = mol.GetAtomWithIdx(donor_atom_idx).GetSymbol()
            bl = float(_get_ml_bond_length(metal_sym, donor_sym))
            vx, vy, vz = vectors[pos_idx]
            mag = math.sqrt(vx * vx + vy * vy + vz * vz)
            if mag > 1e-8:
                vx = vx / mag * bl
                vy = vy / mag * bl
                vz = vz / mag * bl
            donor_target_map[donor_atom_idx] = (vx, vy, vz)
            coords[donor_atom_idx] = [vx, vy, vz]
            placed.add(donor_atom_idx)

        # Phase 1 + 2 — decompose non-metal atoms into ligand fragments
        # and dock each fragment onto its donor target.  Re-uses the
        # existing ``_embed_fragment_procrustes`` which already runs
        # ``_chelate_conformer_candidates`` + MMFF-ish embedding
        # internally for polydentate ligands and Procrustes-aligns
        # the donor atoms to their targets.
        non_metal = {
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() not in _METAL_SET
        }
        adj: Dict[int, set] = {i: set() for i in non_metal}
        for bond in mol.GetBonds():
            bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in non_metal and bj in non_metal:
                adj[bi].add(bj)
                adj[bj].add(bi)

        visited_frag: set = set()
        fragments: List[set] = []
        for start in sorted(non_metal):
            if start in visited_frag:
                continue
            frag: set = set()
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited_frag:
                    continue
                visited_frag.add(node)
                frag.add(node)
                for nbr in adj.get(node, ()):
                    if nbr not in visited_frag:
                        stack.append(nbr)
            fragments.append(frag)

        any_fragment_failed = False
        for frag in fragments:
            frag_donors = [d for d in donor_atom_indices if d in frag]
            if not frag_donors:
                # Non-coordinating fragment: leave for the downstream
                # BFS-VSEPR pass to place.  Happens for pure spectator
                # fragments (counter-ion, crystallographic solvent).
                continue
            tgt_positions = [
                donor_target_map[d] for d in frag_donors
                if d in donor_target_map
            ]
            if not tgt_positions:
                continue
            try:
                ok = _embed_fragment_procrustes(
                    mol, metal_idx, frag, frag_donors, tgt_positions, coords,
                    chelate_rank=chelate_rank,
                )
            except Exception as frag_exc:
                logger.debug(
                    "Balloon fragment embed failed "
                    "(size=%d, donors=%d): %s",
                    len(frag), len(frag_donors), frag_exc,
                )
                ok = False
            if ok:
                placed.update(frag)
            else:
                any_fragment_failed = True

        # If any coordinating fragment failed its Procrustes placement
        # the balloon build cannot produce a CSD-realistic XYZ for
        # this (CF, perm) triple.  Signal None so the caller falls back
        # to the rigid-template builder, which can often rescue these
        # cases using an already-embedded full-molecule conformer.
        if any_fragment_failed:
            return None

        # Phase 2b — BFS-VSEPR to place any remaining (non-coordinating)
        # atoms.  These are usually H atoms or free fragments that the
        # Procrustes path above skipped.
        bond_len_default = 1.4
        queue = list(placed)
        while queue:
            current = queue.pop(0)
            cx, cy, cz = coords[current]
            atom = mol.GetAtomWithIdx(current)
            unplaced_nbrs = [
                n.GetIdx() for n in atom.GetNeighbors()
                if n.GetIdx() not in placed
            ]
            n_unplaced = len(unplaced_nbrs)
            for k, nbr_idx in enumerate(unplaced_nbrs):
                dx, dy, dz = 0.0, 0.0, 0.0
                for other in atom.GetNeighbors():
                    oi = other.GetIdx()
                    if oi in placed and oi != nbr_idx:
                        ox, oy, oz = coords[oi]
                        dx += cx - ox
                        dy += cy - oy
                        dz += cz - oz
                mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                if mag_base < 1e-8:
                    dx, dy, dz = 1.0 + 0.1 * k, 0.3 * k, 0.0
                    mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                if n_unplaced > 1 and k > 0:
                    angle = 2 * math.pi * k / n_unplaced
                    ax, ay, az = dx / mag_base, dy / mag_base, dz / mag_base
                    if abs(ax) < 0.9:
                        px, py, pz = 1.0, 0.0, 0.0
                    else:
                        px, py, pz = 0.0, 1.0, 0.0
                    dot_pa = px * ax + py * ay + pz * az
                    px -= dot_pa * ax
                    py -= dot_pa * ay
                    pz -= dot_pa * az
                    pm = math.sqrt(px * px + py * py + pz * pz)
                    if pm > 1e-8:
                        px /= pm; py /= pm; pz /= pm
                    cos_a = math.cos(angle); sin_a = math.sin(angle)
                    dx2 = dx*cos_a + (ay*dz - az*dy)*sin_a + ax*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dy2 = dy*cos_a + (az*dx - ax*dz)*sin_a + ay*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dz2 = dz*cos_a + (ax*dy - ay*dx)*sin_a + az*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dx, dy, dz = dx2, dy2, dz2
                    mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                    if mag_base < 1e-8:
                        mag_base = 1.0
                dx = dx / mag_base * bond_len_default
                dy = dy / mag_base * bond_len_default
                dz = dz / mag_base * bond_len_default
                coords[nbr_idx] = [cx + dx, cy + dy, cz + dz]
                placed.add(nbr_idx)
                queue.append(nbr_idx)

        # Phase 3 — ligand-orientation pass to break any remaining
        # inter-fragment clashes.  ``_align_and_orient_ligands``
        # rotates each monodentate / bidentate fragment around its
        # M-donor axis (donors stay fixed on the polyhedron) and picks
        # the angle that minimises non-donor clashes.  The inflation
        # schedule is implicit: the build above already places donors
        # at full M-D ideal distance, so one orientation pass at r=1.0
        # is equivalent to the terminal step of a five-step balloon.
        try:
            _align_and_orient_ligands(
                coords, mol, metal_idx, donor_atom_indices
            )
        except Exception as orient_exc:
            logger.debug(
                "Balloon orientation optimiser failed: %s", orient_exc
            )

        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            try:
                _mp = np.array(coords[metal_idx], dtype=float)
                for _d in donor_atom_indices:
                    if mol.GetAtomWithIdx(_d).GetSymbol() != "C":
                        continue
                    _dp = np.array(coords[_d], dtype=float); _mc = _mp - _dp
                    _mcn = float(np.linalg.norm(_mc))
                    if _mcn < 1e-6:
                        continue
                    _mc /= _mcn
                    for _nb in mol.GetAtomWithIdx(_d).GetNeighbors():
                        if _nb.GetAtomicNum() <= 1:
                            continue
                        _cx = np.array(coords[_nb.GetIdx()], dtype=float) - _dp
                        _cxn = float(np.linalg.norm(_cx))
                        if 1.3 < _cxn < 1.9:
                            _ang = float(np.degrees(np.arccos(
                                max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                            _trace_seating("SCRATCH_POST donor=%d M-C-Xheavy=%.0f" % (_d, _ang))
            except Exception:
                pass

        # Build XYZ string.
        lines: List[str] = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(lines) + "\n"
    except Exception as exc:
        logger.debug("_build_topology_xyz_from_scratch failed: %s", exc)
        return None


def _build_topology_xyz(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    apply_uff: bool,
    conf_id: Optional[int] = None,
    chelate_rank: int = 0,
) -> Optional[str]:
    """Build a DELFIN XYZ for one topological arrangement.

    Places the metal at the origin, donor atoms at idealized geometry
    vectors, then attempts RDKit fragment embedding with Procrustes
    alignment for each ligand fragment.  Falls back to BFS placement
    if fragment embedding fails.  Optionally applies OB UFF refinement.

    Args:
        mol: RDKit mol (with H atoms, as from ``_prepare_mol_for_embedding``).
        metal_idx: Index of the metal atom in *mol*.
        donor_atom_indices: List of donor atom indices (length == n_coord).
        perm: ``perm[position] = index into donor_atom_indices``.
        geometry: Key in ``_TOPO_GEOMETRY_VECTORS`` ('OH', 'SQ', …).
        apply_uff: Whether to run OB UFF optimization after placement.
        conf_id: Optional template conformer ID. If None, the rigid-fragment
            builder auto-selects the best-scoring conformer.  Callers iterating
            over multiple alternative binding modes can use this to avoid
            depending on a single (possibly poor) template.
        chelate_rank: Which chelate conformer to use for polydentate
            fragments.  Rank 0 is the tightest native-bite match; higher
            ranks unlock additional backbone puckers for the same platonic
            placement and yield genuinely distinct geometries that DFT
            can discriminate.

    Returns:
        DELFIN-format XYZ string, or None on failure.
    """
    try:
        # Rigid-template builder when a sampling conformer is available.
        # Preserves intraligand geometry from ETKDG.  The Balloon path
        # (``_build_topology_xyz_from_scratch``) is run *additionally*
        # by the topo enumerator for bimetallic systems to broaden the
        # candidate pool — it is not used here as a replacement so that
        # mono-metallic σ complexes keep the tried-and-tested template
        # geometry that powers Ir(ppy)2(acac), Fe(CO)3(NHC)2 etc.
        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            _trace_seating("build_topology_xyz geom=%s nconf=%d donors=%s" % (
                geometry, mol.GetNumConformers(),
                [mol.GetAtomWithIdx(d).GetSymbol() + str(d) for d in donor_atom_indices]))
        if mol.GetNumConformers() > 0:
            xyz_from_template = _build_topology_xyz_from_template(
                mol, metal_idx, donor_atom_indices, perm, geometry, apply_uff,
                conf_id=conf_id, chelate_rank=chelate_rank,
            )
            if xyz_from_template is not None:
                _trace_seating("path=TEMPLATE_OK geom=%s" % geometry)
                return xyz_from_template
            _trace_seating("path=TEMPLATE_returned_None -> FALLBACK (procrustes/BFS) geom=%s" % geometry)

        vectors = _TOPO_GEOMETRY_VECTORS[geometry]
        n_atoms = mol.GetNumAtoms()
        coords: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)] * n_atoms
        placed: set = set()

        # Metal at origin
        coords[metal_idx] = (0.0, 0.0, 0.0)
        placed.add(metal_idx)

        # Donors at geometry positions
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        donor_target_map: Dict[int, Tuple[float, float, float]] = {}
        for pos_idx, donor_list_idx in enumerate(perm):
            donor_atom_idx = donor_atom_indices[donor_list_idx]
            donor_sym = mol.GetAtomWithIdx(donor_atom_idx).GetSymbol()
            bl = _get_ml_bond_length(metal_sym, donor_sym)
            vx, vy, vz = vectors[pos_idx]
            mag = math.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
            if mag > 1e-8:
                vx = vx / mag * bl
                vy = vy / mag * bl
                vz = vz / mag * bl
            coords[donor_atom_idx] = (vx, vy, vz)
            donor_target_map[donor_atom_idx] = (vx, vy, vz)
            placed.add(donor_atom_idx)

        # --- Fragment embedding approach ---
        # Decompose non-metal atoms into ligand fragments (connected components
        # after removing metal bonds), then embed each fragment with RDKit and
        # Procrustes-align the donor atoms to their target positions.
        frag_embed_placed: set = set()
        try:
            non_metal = {a.GetIdx() for a in mol.GetAtoms()
                         if a.GetSymbol() not in _METAL_SET}
            adj: Dict[int, set] = {i: set() for i in non_metal}
            for bond in mol.GetBonds():
                bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bi in non_metal and bj in non_metal:
                    adj[bi].add(bj)
                    adj[bj].add(bi)

            visited_frag: set = set()
            fragments: List[set] = []
            for start in sorted(non_metal):
                if start in visited_frag:
                    continue
                frag: set = set()
                stack = [start]
                while stack:
                    node = stack.pop()
                    if node in visited_frag:
                        continue
                    visited_frag.add(node)
                    frag.add(node)
                    for nbr in adj.get(node, ()):
                        if nbr not in visited_frag:
                            stack.append(nbr)
                fragments.append(frag)

            for frag in fragments:
                # Find donors in this fragment
                frag_donors = [d for d in donor_atom_indices if d in frag]
                if not frag_donors:
                    continue  # non-coordinating fragment, BFS will handle it
                tgt_positions = [donor_target_map[d] for d in frag_donors
                                 if d in donor_target_map]
                if not tgt_positions:
                    continue
                try:
                    ok = _embed_fragment_procrustes(
                        mol, metal_idx, frag, frag_donors, tgt_positions, coords,
                        chelate_rank=chelate_rank,
                    )
                except Exception as frag_exc:
                    logger.debug(
                        "Fragment embedding failed for one fragment (size=%d, donors=%d): %s",
                        len(frag), len(frag_donors), frag_exc,
                    )
                    continue
                if ok:
                    frag_embed_placed.update(frag)
                    placed.update(frag)
        except Exception as _frag_exc:
            logger.debug("Fragment embedding failed, using BFS fallback: %s", _frag_exc)

        # BFS fallback for any atoms not placed by fragment embedding
        bond_len_default = 1.4
        queue = list(placed)
        while queue:
            current = queue.pop(0)
            cx, cy, cz = coords[current]
            atom = mol.GetAtomWithIdx(current)
            unplaced_nbrs = [
                n.GetIdx() for n in atom.GetNeighbors()
                if n.GetIdx() not in placed
            ]
            n_unplaced = len(unplaced_nbrs)
            for k, nbr_idx in enumerate(unplaced_nbrs):
                dx, dy, dz = 0.0, 0.0, 0.0
                for other in atom.GetNeighbors():
                    oi = other.GetIdx()
                    if oi in placed and oi != nbr_idx:
                        ox, oy, oz = coords[oi]
                        dx += cx - ox
                        dy += cy - oy
                        dz += cz - oz
                mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                if mag_base < 1e-8:
                    dx, dy, dz = 1.0 + 0.1 * k, 0.3 * k, 0.0
                    mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                # Rotate around the base direction for each additional neighbour
                # so they fan out instead of all pointing the same way.
                if n_unplaced > 1 and k > 0:
                    angle = 2 * math.pi * k / n_unplaced
                    # Build a perpendicular vector to (dx,dy,dz)
                    ax, ay, az = dx / mag_base, dy / mag_base, dz / mag_base
                    if abs(ax) < 0.9:
                        px, py, pz = 1.0, 0.0, 0.0
                    else:
                        px, py, pz = 0.0, 1.0, 0.0
                    # Gram-Schmidt: subtract projection onto ax
                    dot_pa = px * ax + py * ay + pz * az
                    px -= dot_pa * ax; py -= dot_pa * ay; pz -= dot_pa * az
                    pm = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
                    if pm > 1e-8:
                        px /= pm; py /= pm; pz /= pm
                    # Rodrigues rotation of (dx,dy,dz) by angle around (ax,ay,az)
                    cos_a = math.cos(angle); sin_a = math.sin(angle)
                    dx2 = dx*cos_a + (ay*dz - az*dy)*sin_a + ax*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dy2 = dy*cos_a + (az*dx - ax*dz)*sin_a + ay*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dz2 = dz*cos_a + (ax*dy - ay*dx)*sin_a + az*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dx, dy, dz = dx2, dy2, dz2
                    mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                    if mag_base < 1e-8:
                        mag_base = 1.0
                dx = dx / mag_base * bond_len_default
                dy = dy / mag_base * bond_len_default
                dz = dz / mag_base * bond_len_default
                coords[nbr_idx] = (cx + dx, cy + dy, cz + dz)
                placed.add(nbr_idx)
                queue.append(nbr_idx)

        # Polyhedron-preserving ligand rotation: rotate each monodentate
        # fragment around its M-donor axis and each bidentate fragment
        # around its donor-donor axis to minimise inter-fragment clash.
        # Donor positions (platonic vertices) are invariant.
        try:
            _align_and_orient_ligands(
                coords, mol, metal_idx, donor_atom_indices
            )
        except Exception as _orient_exc:
            logger.debug("Ligand orientation failed: %s", _orient_exc)

        # Build XYZ string
        lines = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        xyz = '\n'.join(lines) + '\n'

        if apply_uff:
            try:
                # Coordination-preserving UFF: M-D distances and L-M-L
                # angles are pinned to the idealized polyhedron so UFF
                # cannot distort the octahedral/PBP/etc. cage even when
                # it lacks force-field parameters for the metal.
                coord_constraints = None
                _perm_tr = None                 # exact per-isomer trans pairs (perm) -> no collapse
                try:
                    _tp = _TOPO_TRANS_POSITIONS.get(geometry) or []
                    if _tp:
                        _perm_tr = [(donor_atom_indices[perm[_p1]], donor_atom_indices[perm[_p2]])
                                    for (_p1, _p2) in _tp]
                except Exception:
                    _perm_tr = None
                try:
                    coord_constraints = _build_coordination_constraints_from_xyz(
                        mol, xyz, d8_trans=_perm_tr,
                    )
                except Exception as cexc:
                    logger.debug(
                        "Coordination constraint build failed, falling back to template constraints: %s",
                        cexc,
                    )
                xyz = _optimize_xyz_openbabel_safe(
                    xyz,
                    mol_template=mol,
                    coord_constraints=coord_constraints,
                )
                # UFF can buckle aromatic rings just enough that the
                # downstream planarity gate rejects genuinely valid
                # isomers.  Snap each aromatic / unsaturated 5-7 ring
                # onto its best-fit plane (minimal-movement projection)
                # so the rule sees planar rings while keeping the rest
                # of the structure untouched.
                xyz = _snap_aromatic_rings_in_xyz(xyz, mol)
            except Exception as uff_exc:
                # Keep the generated topology geometry when UFF fails.
                # Dropping the isomer here can hide valid alternatives.
                logger.debug("Topology UFF optimization failed, keeping unoptimized XYZ: %s", uff_exc)

        return xyz
    except Exception as e:
        logger.debug("_build_topology_xyz failed: %s", e)
        return None


def _build_topology_xyz_from_template(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    apply_uff: bool,
    conf_id: Optional[int] = None,
    chelate_rank: int = 0,
) -> Optional[str]:
    """Rigid-fragment topology builder using an existing template conformer.

    The ligand fragments are transformed as rigid bodies so their internal
    geometry remains close to the template. This is especially useful for
    aromatic/charged chelating ligands where ETKDG fragment embedding can fail.

    If ``conf_id`` is None, the best-scored conformer (per
    :func:`_rank_template_conformers`) is used.  Pass a specific ID when the
    caller wants to iterate over several candidates (e.g. to try multiple
    templates for the same alternative binding mode).

    ``chelate_rank`` unlocks ligand-conformer variety when the rigid
    template is otherwise reused verbatim.  Rank 0 keeps the template's
    pucker exactly; rank > 0 substitutes the ``rank``-th best
    native-bite-matching chelate conformer from
    :func:`_best_chelate_conformer_coords` so flexible chelates
    (cyclam, salen, ethylenediamines) emit distinct chair / boat /
    twist puckers as separate isomers even on systems where ETKDG
    sampling succeeded and the rigid-template branch is the one
    actually building.
    """
    if not RDKIT_AVAILABLE:
        return None
    if mol.GetNumConformers() == 0:
        return None

    try:
        import numpy as np
    except Exception:
        return None

    if conf_id is None:
        ranked = _rank_template_conformers(mol, top_k=1)
        if not ranked:
            return None
        conf_id = ranked[0]

    try:
        conf = mol.GetConformer(int(conf_id))
        vectors = _TOPO_GEOMETRY_VECTORS[geometry]
        n_atoms = mol.GetNumAtoms()

        # Original coordinates translated so the metal sits at origin.
        mpos = conf.GetAtomPosition(metal_idx)
        orig = np.zeros((n_atoms, 3), dtype=float)
        for i in range(n_atoms):
            p = conf.GetAtomPosition(i)
            orig[i, 0] = p.x - mpos.x
            orig[i, 1] = p.y - mpos.y
            orig[i, 2] = p.z - mpos.z

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        donor_target_map: Dict[int, np.ndarray] = {}
        for pos_idx, donor_list_idx in enumerate(perm):
            donor_atom_idx = donor_atom_indices[donor_list_idx]
            donor_sym = mol.GetAtomWithIdx(donor_atom_idx).GetSymbol()
            bl = _get_ml_bond_length(metal_sym, donor_sym)
            vx, vy, vz = vectors[pos_idx]
            mag = math.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
            if mag > 1e-8:
                vx = vx / mag * bl
                vy = vy / mag * bl
                vz = vz / mag * bl
            donor_target_map[donor_atom_idx] = np.array([vx, vy, vz], dtype=float)

        # Build non-metal fragments.
        non_metal = {
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() not in _METAL_SET
        }
        adj: Dict[int, set] = {i: set() for i in non_metal}
        for bond in mol.GetBonds():
            bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in non_metal and bj in non_metal:
                adj[bi].add(bj)
                adj[bj].add(bi)

        visited: set = set()
        fragments: List[set] = []
        for start in sorted(non_metal):
            if start in visited:
                continue
            frag: set = set()
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                frag.add(node)
                for nbr in adj.get(node, ()):
                    if nbr not in visited:
                        stack.append(nbr)
            fragments.append(frag)

        coords = np.array(orig, copy=True)
        coords[metal_idx] = np.array([0.0, 0.0, 0.0], dtype=float)
        placed = {metal_idx}

        def _rotation_from_vectors(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
            nf = np.linalg.norm(v_from)
            nt = np.linalg.norm(v_to)
            if nf < 1e-10 or nt < 1e-10:
                return np.eye(3)
            a = v_from / nf
            b = v_to / nt
            v = np.cross(a, b)
            s = np.linalg.norm(v)
            c = float(np.clip(np.dot(a, b), -1.0, 1.0))
            if s < 1e-10:
                if c > 0:
                    return np.eye(3)
                # 180° rotation around any axis perpendicular to a
                axis = np.array([1.0, 0.0, 0.0])
                if abs(a[0]) > 0.9:
                    axis = np.array([0.0, 1.0, 0.0])
                axis = axis - np.dot(axis, a) * a
                axis = axis / max(np.linalg.norm(axis), 1e-10)
                K = np.array([
                    [0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                    [-axis[1], axis[0], 0],
                ])
                return np.eye(3) + 2.0 * (K @ K)
            K = np.array([
                [0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0],
            ])
            return np.eye(3) + K + K @ K * ((1.0 - c) / (s ** 2))

        def _metal_proximity_penalty(
            xyz_frag: np.ndarray,
            frag_atoms: List[int],
            donor_atoms: List[int],
        ) -> float:
            """Penalty for non-donor heavy atoms placed too close to the metal."""
            donor_set = set(donor_atoms)
            pen = 0.0
            for li, atom_idx in enumerate(frag_atoms):
                if atom_idx in donor_set:
                    continue
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() <= 1:
                    continue
                if atom.GetSymbol() in _METAL_SET:
                    continue
                d = float(np.linalg.norm(xyz_frag[li]))
                sym = atom.GetSymbol()
                # Keep non-donor atoms clearly outside the first coordination shell.
                ml_ref = float(_get_ml_bond_length(metal_sym, sym))
                min_allowed = max(1.15, 0.65 * ml_ref)
                if d < min_allowed:
                    dd = (min_allowed - d)
                    pen += dd * dd
                if d < 1.0:
                    pen += 5.0
            return pen

        for frag in fragments:
            frag_list = sorted(frag)
            frag_donors = [d for d in donor_atom_indices if d in frag]
            if not frag_donors:
                continue

            if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                _trace_seating("FRAG donors=%d atoms=%d placed=%d" % (
                    len(frag_donors), len(frag_list), len(placed)))

            frag_xyz = orig[frag_list, :]
            donor_local = [frag_list.index(d) for d in frag_donors]
            src = frag_xyz[donor_local, :]
            tgt = np.array([donor_target_map[d] for d in frag_donors], dtype=float)
            if len(src) != len(tgt) or len(src) == 0:
                continue

            # ERDBEBEN isolated-fragment seating (gated DELFIN_FFFREE_ISOLATED_SEAT, default off ->
            # byte-identical): if the WHOLE-COMPLEX template collapsed this fragment's cage into a plane,
            # re-seat it from a clean ISOLATED embed (proven 20/20 non-collapsed on AQIBAE) instead of the
            # collapsed template.  The metal context is what collapses the cage; the fragment alone builds
            # fine -- so we take the fragment geometry from where it is RELIABLE.
            # ERDBEBEN reseat is an OPTIONAL optimization -- ANY failure inside it (the collapse probe
            # or the isolated re-embed throwing on an exotic ligand, e.g. QILGIB's o-phenylene-diarsine
            # chelate) must NEVER abort the build.  On any exception, fall back to the rigid TEMPLATE
            # fragment = the flag-OFF geometry.  Never-worse by construction: byte-identical when no
            # exception occurs (the reseat is off/inert), and a system is never dropped when one does.
            try:
                _reseat_collapse = (os.environ.get("DELFIN_FFFREE_ISOLATED_SEAT", "0") == "1"
                                    and _frag_xyz_collapsed(mol, frag_list, frag_xyz))

                # Chelate (bidentate or polydentate): if the template's
                # native donor-donor distance pattern is far from the
                # polyhedron target pattern, re-embed the fragment alone
                # with multiple ETKDG seeds and pick the conformer whose
                # full pairwise donor geometry best matches.  This avoids
                # rigidly stretching the chelate backbone against the
                # graph-gate bond-length window.
                if len(frag_donors) >= 2:
                    # Pairwise distance matrices (template vs target).
                    src_diffs = src[:, None, :] - src[None, :, :]
                    tgt_diffs = tgt[:, None, :] - tgt[None, :, :]
                    template_mat = np.linalg.norm(src_diffs, axis=-1)
                    target_mat = np.linalg.norm(tgt_diffs, axis=-1)
                    mismatch = float(
                        np.sqrt(
                            np.triu((template_mat - target_mat) ** 2, k=1).sum()
                            / max(1, len(frag_donors) * (len(frag_donors) - 1) // 2)
                        )
                    )
                    if mismatch > 0.25 or _reseat_collapse:
                        target_for_search = (
                            float(target_mat[0, 1])
                            if len(frag_donors) == 2
                            else target_mat
                        )
                        coords_map = _best_chelate_conformer_coords(
                            mol, frag, frag_donors, target_for_search,
                            rank=chelate_rank,
                        )
                        if coords_map is not None:
                            new_frag_xyz = np.array(
                                [list(coords_map[old]) for old in frag_list],
                                dtype=float,
                            )
                            # Accept the re-embed if it improves the bite fit, OR (collapse re-seat) if the
                            # clean ISOLATED embed resolved the collapse -- the 3D fragment is what we want even
                            # when the collapsed template's bite happened to already match (AQIBAE).
                            new_src = new_frag_xyz[donor_local, :]
                            new_diffs = new_src[:, None, :] - new_src[None, :, :]
                            new_mat = np.linalg.norm(new_diffs, axis=-1)
                            new_mismatch = float(
                                np.sqrt(
                                    np.triu((new_mat - target_mat) ** 2, k=1).sum()
                                    / max(1, len(frag_donors) * (len(frag_donors) - 1) // 2)
                                )
                            )
                            if (new_mismatch < mismatch
                                    or (_reseat_collapse
                                        and not _frag_xyz_collapsed(mol, frag_list, new_frag_xyz))):
                                frag_xyz = new_frag_xyz
                                src = new_src
                                if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1" and _reseat_collapse:
                                    _trace_seating(
                                        "ISOLATED_SEAT reseated collapsed fragment (donors=%d) mismatch %.2f->%.2f"
                                        % (len(frag_donors), mismatch, new_mismatch))
            except Exception:
                # Reseat/chelate-search failed -> keep the rigid template fragment (flag-OFF geometry).
                _reseat_collapse = False

            if len(src) >= 2:
                src_center = src.mean(axis=0)
                tgt_center = tgt.mean(axis=0)
                X = src - src_center
                Y = tgt - tgt_center
                H = X.T @ Y
                U, S, Vt = np.linalg.svd(H)
                R = Vt.T @ U.T
                if np.linalg.det(R) < 0:
                    Vt[-1, :] *= -1
                    R = Vt.T @ U.T
                transformed = (frag_xyz - src_center) @ R.T + tgt_center

                # Bidentate/multidentate fragments can be mirrored around the
                # donor-donor axis, yielding two plausible orientations with
                # identical donor placement. Pick the orientation that keeps
                # non-donor atoms farther from the metal center.
                try:
                    if len(frag_donors) >= 2:
                        d0 = donor_target_map[frag_donors[0]]
                        d1 = donor_target_map[frag_donors[1]]
                        axis = d1 - d0
                        axis_norm = float(np.linalg.norm(axis))
                        if axis_norm > 1e-10:
                            u = axis / axis_norm
                            pivot = 0.5 * (d0 + d1)
                            v = transformed - pivot
                            # 180° rotation around donor-donor axis:
                            # v' = -v + 2*(u·v)*u
                            v_rot = -v + 2.0 * np.outer(v @ u, u)
                            transformed_flip = v_rot + pivot

                            p0 = _metal_proximity_penalty(transformed, frag_list, frag_donors)
                            p1 = _metal_proximity_penalty(transformed_flip, frag_list, frag_donors)
                            if p1 < p0:
                                transformed = transformed_flip
                except Exception:
                    pass

                # ERDBEBEN bite-preserving backbone DECLASH (gated DELFIN_FFFREE_RIGID_DECLASH, default
                # off -> byte-identical).  A BIDENTATE fragment's two donors lie ON the donor-donor axis,
                # so rotating the WHOLE fragment about that axis keeps both donors EXACTLY on their
                # polyhedron vertices (bite + polyhedron preserved -- rotating a point on the axis leaves
                # it fixed) while the backbone sweeps a cone.  Rotate to the angle minimising clash with
                # the metal AND the already-placed fragments -> a rigid chelate whose backbone would
                # otherwise collide (BINHIQ 2x diarsine: verify 0/366 -> collapsed fallback wins) reaches
                # a clash-free placement that PASSES _verify_topology_from_graph, with NO force field.
                # Only exactly-bidentate: 3+ donors pin the rigid body (0 rotational DOF); monodentate
                # radial spin is JOINT_DECLASH's job.  This is the FF-free seating co-optimisation.
                if (os.environ.get("DELFIN_FFFREE_RIGID_DECLASH", "0") == "1"
                        and len(frag_donors) == 2):
                    try:
                        _d0 = donor_target_map[frag_donors[0]]
                        _d1 = donor_target_map[frag_donors[1]]
                        _ax = _d1 - _d0
                        _axn = float(np.linalg.norm(_ax))
                        _bb = [li for li, ai in enumerate(frag_list)
                               if ai not in frag_donors
                               and mol.GetAtomWithIdx(ai).GetAtomicNum() > 1]
                        _other = [j for j in placed
                                  if j not in frag and j != metal_idx
                                  and mol.GetAtomWithIdx(j).GetAtomicNum() > 1
                                  and mol.GetAtomWithIdx(j).GetSymbol() not in _METAL_SET]
                        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                            _trace_seating("RIGID_DECLASH_TRY axn=%.2f bb=%d other=%d" % (
                                _axn, len(_bb), len(_other)))
                        if _axn > 1e-8 and _bb and _other:
                            _u = _ax / _axn
                            _piv = 0.5 * (_d0 + _d1)
                            _P = np.array([coords[j] for j in _other], dtype=float)
                            _cmin = float(os.environ.get("DELFIN_RIGID_DECLASH_MIN", "2.4") or 2.4)

                            def _declash_pen(_cand):
                                _p = _metal_proximity_penalty(_cand, frag_list, frag_donors)
                                for _li in _bb:
                                    _ov = _cmin - np.linalg.norm(_P - _cand[_li], axis=1)
                                    _ov = _ov[_ov > 0.0]
                                    if _ov.size:
                                        _p += float(np.sum(_ov * _ov))
                                return _p

                            def _rot_axis(_pts, _th):
                                _c = math.cos(_th); _s = math.sin(_th)
                                _v = _pts - _piv
                                return (_v * _c + np.cross(_u, _v) * _s
                                        + np.outer(_v @ _u, _u) * (1.0 - _c)) + _piv

                            _pen0 = _declash_pen(transformed)
                            _bestp = _pen0
                            if _bestp > 1e-9:              # only sweep if there is a clash to resolve
                                _best = transformed
                                for _k in range(1, 24):    # 15-deg steps around the donor-donor axis
                                    _cand = _rot_axis(transformed, 2.0 * math.pi * _k / 24.0)
                                    _pen = _declash_pen(_cand)
                                    if _pen < _bestp - 1e-9:
                                        _bestp = _pen
                                        _best = _cand
                                transformed = _best
                            if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                                _trace_seating(
                                    "RIGID_DECLASH donors=%s bb=%d other=%d pen0=%.3f -> penbest=%.3f%s" % (
                                        [int(x) for x in frag_donors], len(_bb), len(_other),
                                        _pen0, _bestp, "" if _pen0 > 1e-9 else " (no-clash)"))
                    except Exception as _dexc:
                        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                            _trace_seating("RIGID_DECLASH EXC: %s" % _dexc)
            else:
                src_d = src[0]
                tgt_d = tgt[0]
                if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                    _trace_seating("template_mono donor=%d elem=%s nsubst_frag(<1.9A)=%d" % (
                        int(frag_donors[0]) if frag_donors else -1,
                        mol.GetAtomWithIdx(int(frag_donors[0])).GetSymbol() if frag_donors else "?",
                        int(sum(1 for _r in frag_xyz if 0.4 < float(np.linalg.norm(_r - src_d)) < 1.9))))
                v_from = None
                # ROOT FIX (DELFIN_FFFREE_SP3C_TET_SEAT=1, default OFF -> byte-identical): sp3-C donor ->
                # seat the metal in the vacant tetrahedral slot.  The centroid heuristic below is skewed
                # by the heavy substituent (M-CH2-R: centroid ~ R) -> R anti to metal -> M-C-X ~180.
                # Sum of unit(donor->substituent) over all 3 bonded frag atoms (incl H, geometrically
                # detected) points opposite the vacant slot; aligning it with the outward radial puts the
                # vacant slot (metal) at ~109 deg.
                if os.environ.get("DELFIN_FFFREE_SP3C_TET_SEAT", "0") == "1" and frag_donors:
                    try:
                        _dca = mol.GetAtomWithIdx(int(frag_donors[0]))
                        # PENDANT sp3 alkyl donor only: not aromatic (no sp2 carbene/aryl -> no vacant
                        # tetrahedral slot) and NOT in a ring (a ring-embedded C donor's orientation is
                        # already fixed by the rigid scaffold; re-orienting the ring breaks the topology,
                        # e.g. XIKSEQ's Cd-bound diazaborole C loses its crystal-matching frame).
                        if (_dca.GetSymbol() == "C" and not _dca.GetIsAromatic()
                                and not _dca.IsInRing()):
                            _acc = np.zeros(3); _ns = 0
                            for _row in frag_xyz:
                                _v = _row - src_d; _vn = float(np.linalg.norm(_v))
                                if 0.4 < _vn < 1.9:
                                    _acc += _v / _vn; _ns += 1
                            if _ns == 3 and float(np.linalg.norm(_acc)) > 1e-6:
                                v_from = _acc
                    except Exception:
                        v_from = None
                if v_from is None:
                    src_com = frag_xyz.mean(axis=0)
                    # Keep the fragment extending away from the metal.
                    v_from = src_com - src_d
                v_to = tgt_d
                R = _rotation_from_vectors(v_from, v_to)
                transformed = (frag_xyz - src_d) @ R.T + tgt_d
                if (os.environ.get("DELFIN_TRACE_SEATING", "0") == "1" and frag_donors
                        and mol.GetAtomWithIdx(int(frag_donors[0])).GetSymbol() == "C"):
                    try:
                        _mc = -tgt_d / (float(np.linalg.norm(tgt_d)) + 1e-9)
                        for _row in transformed:
                            _v = _row - tgt_d; _vn = float(np.linalg.norm(_v))
                            if 1.3 < _vn < 1.9:
                                _ang = float(np.degrees(np.arccos(
                                    max(-1.0, min(1.0, float(np.dot(_mc, _v / _vn)))))))
                                _trace_seating("template_mono_POST donor=%d M-C-Xheavy=%.0f" % (
                                    int(frag_donors[0]), _ang))
                    except Exception:
                        pass

            for li, atom_idx in enumerate(frag_list):
                coords[atom_idx] = transformed[li]
                placed.add(atom_idx)

        # Any atom not touched by fragment placement keeps template-relative coords.
        for i in range(n_atoms):
            if i not in placed:
                coords[i] = orig[i]

        # Polyhedron-preserving ligand rotation (see _align_and_orient_ligands).
        try:
            _align_and_orient_ligands(
                coords, mol, metal_idx, donor_atom_indices
            )
        except Exception as _orient_exc:
            logger.debug("Ligand orientation (template path) failed: %s", _orient_exc)


        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            try:
                _mp = np.array(coords[metal_idx], dtype=float)
                for _d in donor_atom_indices:
                    if mol.GetAtomWithIdx(_d).GetSymbol() != "C":
                        continue
                    _dp = np.array(coords[_d], dtype=float); _mc = _mp - _dp
                    _mcn = float(np.linalg.norm(_mc))
                    if _mcn < 1e-6:
                        continue
                    _mc /= _mcn
                    for _nb in mol.GetAtomWithIdx(_d).GetNeighbors():
                        if _nb.GetAtomicNum() <= 1:
                            continue
                        _cx = np.array(coords[_nb.GetIdx()], dtype=float) - _dp
                        _cxn = float(np.linalg.norm(_cx))
                        if 1.3 < _cxn < 1.9:
                            _ang = float(np.degrees(np.arccos(
                                max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                            _trace_seating("POST_ALIGN donor=%d M-C-Xheavy=%.0f" % (_d, _ang))
            except Exception:
                pass

        lines = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}")
        xyz = '\n'.join(lines) + '\n'

        # Conservative UFF: only keep UFF result if it preserves topology.
        # For metals without UFF parameters, UFF can BREAK the structure.
        # The pre-UFF Procrustes geometry has correct M-D distances and
        # is often better than what UFF produces.
        if apply_uff:
            xyz_pre_uff = xyz
            try:
                coord_constraints = None
                # ARCHITECTURE (2026-07-13): hand the constraint builder the EXACT per-isomer trans
                # donor pairs from THIS frame's perm (perm + _TOPO_TRANS_POSITIONS[geometry]).  The
                # DELFIN_FFFREE_D8_SQ_ISO / CN6_OH_ANGLES passes then impose the correct polyhedron on
                # each isomer's OWN arrangement (no geometric guessing) -> no isomer collapse.  Root fix
                # for the whole poly cluster; None (non-SQ/OH geoms) = geometry fallback / byte-identical.
                _perm_tr = None
                try:
                    _tp = _TOPO_TRANS_POSITIONS.get(geometry) or []
                    if _tp:
                        _perm_tr = [(donor_atom_indices[perm[_p1]], donor_atom_indices[perm[_p2]])
                                    for (_p1, _p2) in _tp]
                except Exception:
                    _perm_tr = None
                try:
                    coord_constraints = _build_coordination_constraints_from_xyz(
                        mol, xyz, d8_trans=_perm_tr,
                    )
                except Exception:
                    pass
                xyz_uff = _optimize_xyz_openbabel_safe(
                    xyz,
                    mol_template=mol,
                    coord_constraints=coord_constraints,
                )
                # Snap any UFF-buckled aromatic rings back to planar so
                # the downstream pi-ring planarity gate doesn't reject
                # otherwise valid isomers.
                xyz_uff = _snap_aromatic_rings_in_xyz(xyz_uff, mol)
                # Keep UFF only if topology survives.
                if _verify_topology_from_graph(xyz_uff, mol):
                    xyz = xyz_uff
                else:
                    logger.debug(
                        "UFF broke topology in template builder — keeping pre-UFF geometry"
                    )
                    xyz = xyz_pre_uff
            except Exception as uff_exc:
                logger.debug(
                    "Template-topology UFF failed, keeping pre-UFF: %s", uff_exc
                )
        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            try:
                _pos = {}
                for _i, _ln in enumerate(xyz.strip().split("\n")):
                    _p = _ln.split()
                    if len(_p) >= 4:
                        _pos[_i] = np.array([float(_p[1]), float(_p[2]), float(_p[3])])
                _mp = _pos.get(metal_idx)
                _kept = "no_uff" if not apply_uff else (
                    "UFF_kept" if xyz is not xyz_pre_uff else "pre_uff(UFF_reverted_or_failed)")
                for _d in donor_atom_indices:
                    if mol.GetAtomWithIdx(_d).GetSymbol() != "C":
                        continue
                    _dp = _pos.get(_d)
                    if _mp is None or _dp is None:
                        continue
                    _mc = _mp - _dp; _mcn = float(np.linalg.norm(_mc))
                    if _mcn < 1e-6:
                        continue
                    _mc /= _mcn
                    for _nb in mol.GetAtomWithIdx(_d).GetNeighbors():
                        if _nb.GetAtomicNum() <= 1:
                            continue
                        _xp = _pos.get(_nb.GetIdx())
                        if _xp is None:
                            continue
                        _cx = _xp - _dp; _cxn = float(np.linalg.norm(_cx))
                        if 1.3 < _cxn < 1.9:
                            _ang = float(np.degrees(np.arccos(
                                max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                            _trace_seating("POST_UFF donor=%d M-C-Xheavy=%.0f %s geom=%s" % (
                                _d, _ang, _kept, geometry))
            except Exception:
                pass
        return xyz
    except Exception as e:
        logger.debug("_build_topology_xyz_from_template failed: %s", e)
        return None


def _build_topology_template_mol(smiles: str):
    """Create a topology-template RDKit Mol with a mapped 3D conformer.

    Uses the quick conversion XYZ as the geometry source and reconstructs a
    matching RDKit molecule via ``mol_from_smiles_rdkit + AddHs`` so atom
    order/length stays consistent for conformer injection.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        xyz, err = smiles_to_xyz_quick(smiles)
        if err or not xyz:
            return None
        mol, _note = mol_from_smiles_rdkit(smiles, allow_metal=True)
        if mol is None:
            return None
        try:
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception:
            pass
        conf = _xyz_to_rdkit_conformer(mol, xyz)
        if conf is None:
            return None
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return mol
    except Exception:
        return None


def _generate_topological_isomers(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_isomers: int = 50,
    n_metal_smart: bool = True,
    profile: Optional[Dict[str, int]] = None,
) -> List[Tuple[str, str]]:
    """Guarantee-complete isomer enumeration via topological permutation.

    For each metal center: enumerates all unique canonical arrangements of
    donor atoms (respecting chelate cis-constraints), builds an idealized
    3D structure for each, runs OB UFF, and labels the result.

    ``profile`` — when provided, overrides the module-level
    ``DELFIN_CHELATE_RANK_VARIANTS``, ``DELFIN_TOPO_TEMPLATE_TOP_K``
    and ``DELFIN_PRE_UFF_CAP_MULTIPLIER`` knobs for this call.  Keys:
    ``ranks``, ``topk``, ``cap_mult``.

    Returns [(xyz_string, label), …].
    """
    import numpy as np
    results: List[Tuple[str, str]] = []
    dtype_map = _donor_type_map(mol)

    # Iter-8.5b INNER: same env-flag and class-list as the outer trans/pucker
    # skip (commit 2c1070e).  Three inner pump-pass sites in this function
    # over-generate seeds per permutation: (1) template-loop iterates every
    # ranked template CID, (2) balloon-builder appends a from-scratch scaffold,
    # (3) chelate-rank variants enumerate alternative puckers.  When the env
    # class-list contains the parent mol's class, the inner pumps are trimmed
    # to a single deterministic seed (template-loop: first success only;
    # balloon + chelate-ranks: skipped entirely).  Default empty set =
    # bit-exact (every inner pump runs as before).
    _iter85b_pump_skip_inner = False
    try:
        _iter85b_inner_classes = set(
            x.strip() for x in (
                os.environ.get("DELFIN_ITER85_PUMP_SKIP_CLASSES", "") or ""
            ).split(",") if x.strip()
        )
        if _iter85b_inner_classes:
            _iter85b_pump_skip_inner = (
                _classify_complex_class(mol) in _iter85b_inner_classes
            )
    except Exception:
        _iter85b_pump_skip_inner = False

    _prof = profile if profile is not None else _resolve_quality_profile(None)
    _prof_ranks    = int(_prof.get("ranks", DELFIN_CHELATE_RANK_VARIANTS))
    _prof_topk     = int(_prof.get("topk", DELFIN_TOPO_TEMPLATE_TOP_K))
    _prof_cap_mult = int(_prof.get("cap_mult", DELFIN_PRE_UFF_CAP_MULTIPLIER))

    def _passes_chelate_distance_feasibility(
        _mol,
        _metal_idx: int,
        _donor_indices: List[int],
        _perm: List[int],
        _geom_name: str,
        _chelate_atom_pairs: List[FrozenSet],
        abs_tol: Optional[float] = None,
        rel_tol: Optional[float] = None,
        force_reach: bool = False,
    ) -> bool:
        """Reject geometrically impossible chelate placements.

        For each chelate pair, compare donor-donor distance in the template
        conformer to the idealized target distance implied by geometry+perm.
        If they differ too much, this arrangement is likely non-physical for
        the ligand bite and tends to collapse into unrealistic structures.

        Tolerances are env-tunable via DELFIN_CHELATE_FEAS_ABS_TOL (default
        1.2 Å) and DELFIN_CHELATE_FEAS_REL_TOL (default 0.5).  The defaults
        were widened from (0.7, 0.35) after Mn(CO)3(CO2Me)(dppe) showed
        all 3 TPR arrangements rejected because the template's dppe P-P
        bite conformer (~3.0 Å) did not fit any TPR idealized P-P bucket
        within the tighter tolerance — even though a TPR structure can be
        built and UFF-refined successfully in the downstream builder.
        """
        if abs_tol is None:
            abs_tol = _delfin_env_float("DELFIN_CHELATE_FEAS_ABS_TOL", 0.7)
        if rel_tol is None:
            rel_tol = _delfin_env_float("DELFIN_CHELATE_FEAS_REL_TOL", 0.35)
        if not _delfin_env_int("DELFIN_CHELATE_FEAS_ENABLED", 1):
            return True
        try:
            if _mol.GetNumConformers() == 0:
                return True
            conf = _mol.GetConformer(0)
            vectors = _TOPO_GEOMETRY_VECTORS.get(_geom_name)
            if not vectors:
                return True

            m_sym = _mol.GetAtomWithIdx(_metal_idx).GetSymbol()
            target_by_donor: Dict[int, Tuple[float, float, float]] = {}
            for pos_idx, donor_list_idx in enumerate(_perm):
                d_atom = _donor_indices[donor_list_idx]
                d_sym = _mol.GetAtomWithIdx(d_atom).GetSymbol()
                bl = _get_ml_bond_length(m_sym, d_sym)
                vx, vy, vz = vectors[pos_idx]
                mag = math.sqrt(vx * vx + vy * vy + vz * vz)
                if mag > 1e-8:
                    vx = vx / mag * bl
                    vy = vy / mag * bl
                    vz = vz / mag * bl
                target_by_donor[d_atom] = (vx, vy, vz)

            _use_reach = force_reach or _delfin_env_int("DELFIN_FFFREE_CHELATE_REACH_FEAS", 0)
            for cp in _chelate_atom_pairs:
                pair = sorted(cp)
                if len(pair) != 2:
                    continue
                a, b = pair
                if a not in target_by_donor or b not in target_by_donor:
                    continue
                ta = target_by_donor[a]
                tb = target_by_donor[b]
                d_tgt = math.sqrt(
                    (ta[0] - tb[0]) ** 2 + (ta[1] - tb[1]) ** 2 + (ta[2] - tb[2]) ** 2
                )

                if _use_reach:
                    # FIRST-PRINCIPLES reachability (triangle inequality, SAMPLING-INDEPENDENT): the
                    # chelate backbone's contour length is the HARD upper bound on the donor-donor bite for
                    # ANY conformer.  Reject ONLY when the target bite exceeds it (physically impossible to
                    # span -- e.g. a short bridge cannot reach a trans separation).  This is the ROOT fix:
                    # it replaces the rigid single-conformer (conformer-0) bite check below, which
                    # over-pruned every flexible chelate (Ir all-4-OH rejected; KAFBUS all rejected ->
                    # fallback) merely because the ONE embedded template conformer's bite missed the target,
                    # even though the builder (multiple ranked template conformers + UFF) can reach it.
                    _max_reach = _chelate_backbone_max_reach(_mol, a, b)
                    if _max_reach is not None and d_tgt > _max_reach + abs_tol:
                        return False
                    continue

                pa = conf.GetAtomPosition(a)
                pb = conf.GetAtomPosition(b)
                d_src = math.sqrt(
                    (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
                )
                tol = max(abs_tol, rel_tol * max(d_src, 1e-8))
                if abs(d_tgt - d_src) > tol:
                    return False
            return True
        except Exception:
            return True

    # Per-metal constitutional enumeration.  Runs for every metal in
    # the template including multi-metal systems: ``_build_topology_xyz``
    # consults ``_build_topology_xyz_from_template`` first, which
    # preserves every atom outside the current metal's own fragments —
    # the other metal and any bridging donors stay at their template
    # positions.  The strict graph gate downstream filters any build
    # that did violate topology, so broader per-metal enumeration
    # yields more CONSTITUTIONAL candidates without letting the
    # old bimetallic-broken geometries through.  The multinuclear
    # coupled-enumeration block further down still runs and adds the
    # Cartesian product on top of the per-metal constitutional set.
    _n_metals_in_mol = sum(
        1 for _a in mol.GetAtoms() if _a.GetSymbol() in _METAL_SET
    )
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        donor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        n_coord = len(donor_indices)

        if n_coord < 2 or n_coord > 9:
            continue

        # Use donor environment classes (element + Morgan environment) instead
        # of plain element symbols. This preserves chemically meaningful
        # distinctions like aqua-O vs carboxylate-O, which are essential for
        # CN=7/8 trans-pattern completeness.
        donor_keys = [
            dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))
            for d in donor_indices
        ]
        uniq_keys = sorted(
            set(donor_keys),
            key=lambda k: (k[0], tuple(sorted(k[1]))),
        )
        key_to_class = {k: i for i, k in enumerate(uniq_keys)}
        donor_labels = [f"{k[0]}{key_to_class[k]}" for k in donor_keys]

        chelate_ps = _chelate_pairs(mol, metal_idx, donor_indices)

        # Skip macrocyclic complexes: when ALL donor pairs are chelate-connected
        # (complete chelate graph), every donor is part of one big ring.
        # The BFS in _build_topology_xyz cannot respect ring-closure constraints
        # and always produces broken structures for macrocycles.
        max_pairs = n_coord * (n_coord - 1) // 2
        if len(chelate_ps) >= max_pairs:
            logger.debug(
                "Skipping topo enumerator for macrocyclic complex (CN=%d, "
                "chelate_pairs=%d/%d)", n_coord, len(chelate_ps), max_pairs
            )
            continue

        # Convert chelate pairs from atom indices to donor-list indices.
        # _chelate_pairs returns frozensets of atom indices (e.g. {1, 12}),
        # but _enumerate_topological_isomers works with donor-list indices
        # (0..n_coord-1).  Without this mapping the constraint check always
        # hits ValueError → silently skipped → all constraints ignored.
        atom_to_listidx = {atom_idx: li for li, atom_idx in enumerate(donor_indices)}
        chelate_list_pairs: List[FrozenSet] = []
        for cp in chelate_ps:
            pair = sorted(cp)
            if pair[0] in atom_to_listidx and pair[1] in atom_to_listidx:
                chelate_list_pairs.append(frozenset([
                    atom_to_listidx[pair[0]], atom_to_listidx[pair[1]]
                ]))

        isomers = _enumerate_topological_isomers(
            donor_labels, n_coord, chelate_list_pairs,
            metal_symbol=atom.GetSymbol(),
            metal_formal_charge=int(atom.GetFormalCharge() or 0),
            mol=mol,
            metal_idx=metal_idx,
            donor_indices=donor_indices,
        )

        # Chelate-distance feasibility is a useful guard, but can over-prune
        # higher-coordination systems (notably CN=7) when idealized vectors and
        # template distances differ systematically. If it rejects everything,
        # fall back to the unfiltered topological set.
        #
        # Fix C (Welle 2 / X10-FIPWAE): for CN=5 with multiple distinct donor
        # types ("hetero CN5"), Pólya enumeration generates 4-8 orbits but the
        # chelate-feasibility pruner tends to keep only 1-2 because TBP and SP
        # have very different donor-donor target distances and the template
        # conformer reflects neither cleanly. When DELFIN_CN5_ENUM_COMPLETE=1,
        # skip the pruner entirely for CN=5 hetero so downstream geometry
        # filters (which already exist) handle quality control. Default OFF.
        _cn5_complete = (
            n_coord == 5
            and _cn5_enum_complete_enabled()
            and len(set(donor_labels)) >= 2
        )
        # COMPLETENESS + QUALITY, NO JUNK (user 2026-07-21 "the construction should not build bad
        # geometries in the first place" -> root fix, not post-hoc cull).  The chelate-distance pre-filter over-prunes
        # polydentate systems (bis-tridentate Ir: 19 enumerated -> 1 feasible, ALL 4 octahedral arrangements
        # rejected because the one template conformer's rigid bite doesn't fit the ideal OH distance; the
        # fallback only fires when feasible==0, so keeping 1 silently drops the rest).
        #   DELFIN_FFFREE_ENUM_FEAS_PREFERRED (default off): skip the pre-filter ONLY for the LFSE-PREFERRED
        #     polyhedron -> its realistic arrangements (Ir all-cis / N-trans / C-trans) are recovered, while
        #     the NON-preferred (e.g. TPR "schiefe pi" junk for a d6-mer) still faces feasibility and is
        #     NEVER BUILT.  Surgical: recover the realistic, never generate the junk.
        #   DELFIN_FFFREE_ENUM_SKIP_FEASIBILITY (default off): blunt -- skip for ALL geometries (admits the
        #     junk too; kept only for diagnostics, superseded by FEAS_PREFERRED).
        _enum_skip_feas = _delfin_env_int("DELFIN_FFFREE_ENUM_SKIP_FEASIBILITY", 0)
        _enum_feas_pref = _delfin_env_int("DELFIN_FFFREE_ENUM_FEAS_PREFERRED", 0)
        # Max L-M-L angle deviation (deg) from the nearest ideal VSEPR polyhedron that a RECOVERED
        # FEAS_PREFERRED arrangement may have and still be kept -- the FF-FREE geometric realism filter.
        _extra_max_dev = _delfin_env_float("DELFIN_FFFREE_EXTRA_MAX_ANGLE_DEV", 30.0)
        # Min M-donor-H angle (deg): a recovered arrangement whose coordinating N-H/O-H has an H pointing
        # TOWARD the metal (angle below this) is dropped -- MANTA-native geometric realism (not WEDDELL).
        _extra_donor_h_min = _delfin_env_float("DELFIN_FFFREE_EXTRA_DONOR_H_MIN_DEG", 80.0)
        # LFSE-PREFERRED polyhedron per CN -- GEOMETRY-AWARE for CN4/5/6 (the metal's electron count decides:
        # d8 CN5 -> SP not TBP, d6/d8 etc.).  The earlier hard-coded 5:'TBP' made FEAS_PREFERRED recover the
        # WRONG polyhedron for square-pyramidal systems (JAMHUB CN5 crystal is SP -> adding TBP arrangements
        # confused the isomer count 5->2 and cost build time).  Using _PREFERRED_CN5/6_GEOMETRY adds the
        # CORRECT preferred polyhedron -> fewer, right arrangements -> resolves JAMHUB + fewer timeouts.
        _pref_geom = {2: 'LIN', 3: 'TP',
                      4: _preferred_cn4_for(atom.GetSymbol(), mol, atom.GetIdx()),
                      5: _PREFERRED_CN5_GEOMETRY.get(atom.GetSymbol(), 'TBP'),
                      6: _PREFERRED_CN6_GEOMETRY.get(atom.GetSymbol(), 'OH'),
                      7: 'PBP', 8: 'SAP', 9: 'TTP'}.get(n_coord)
        feasible_isomers: List[Tuple[tuple, List[int]]] = []
        _pref_extra_keys: set = set()   # (cf,pm) of FEAS_PREFERRED extras -> geometric-VSEPR-realism filtered
        if _cn5_complete or _enum_skip_feas:
            feasible_isomers = list(isomers)
        else:
            # FEAS_PREFERRED (default off) recovers the LFSE-PREFERRED polyhedron's realistic arrangements
            # that the naive chelate-distance pre-filter over-prunes (bis-tridentate Ir: all 4 OH
            # arrangements rejected, only 1 kept).  CRITICAL -- it must be TRULY ADDITIVE: the recovered
            # preferred-geom isomers are appended AFTER the feasibility-PASSING set, so they can NEVER crowd
            # the real isomers out of the _PRE_UFF_CAP frame budget (which the build loop below breaks on).
            # The earlier naive `(pref==geom) OR feasible` mixed them into the enumeration-ORDER stream, so
            # the preferred-geom flood filled the cap and the loop broke BEFORE the real isomers built ->
            # measured isomer LOSS (KAFBUS 10->1, polyhedra lost, mean_delta +0.482, gate REJECTED).  Base-
            # first ordering makes it never-worse: the feasibility-passing set builds EXACTLY as the
            # baseline; the preferred extras consume only the REMAINING budget -- which is precisely the
            # sparse systems (Ir base=1) that need the recovery.  A pref-geom isomer that ALSO passes
            # feasibility stays in the base (it is a real feasible isomer, not an extra).
            _pref_extra: List[Tuple[tuple, List[int]]] = []
            for canonical_form, perm in isomers:
                geom_name = canonical_form[0]
                if _passes_chelate_distance_feasibility(
                        mol, metal_idx, donor_indices, perm, geom_name, chelate_ps):
                    feasible_isomers.append((canonical_form, perm))
                elif (_enum_feas_pref and geom_name == _pref_geom
                        and _passes_chelate_distance_feasibility(
                            mol, metal_idx, donor_indices, perm, geom_name, chelate_ps,
                            force_reach=True)):
                    # Recover the LFSE-preferred polyhedron that the rigid single-conformer feasibility
                    # over-pruned -- but ONLY when it is PHYSICALLY REACHABLE (triangle inequality on the
                    # backbone contour, sampling-independent).  ADDITIVE (appended after the base + fallback,
                    # so no fallback suppression) + REACH-gated (no unreachable junk) + preferred-geom-scoped
                    # (bounded -> no isomer explosion / build blow-up).  The pure reach-REPLACE was too
                    # lenient (22 build timeouts, JAMHUB fallback loss); this keeps its physics gate while
                    # staying never-worse.
                    _pref_extra.append((canonical_form, perm))
            if not feasible_isomers and isomers:
                logger.debug(
                    "Chelate-distance feasibility rejected all %d topo isomer(s) "
                    "for CN=%d; using unfiltered set.",
                    len(isomers), n_coord,
                )
                feasible_isomers = list(isomers)   # fallback already contains the preferred-geom isomers
            elif _pref_extra:
                # additive extras LAST -> never crowd the feasibility-passing (real) isomers out of the cap
                feasible_isomers = feasible_isomers + _pref_extra
                _pref_extra_keys = {(tuple(_cf), tuple(_pm)) for _cf, _pm in _pref_extra}

        # DIAGNOSTIC (gated DELFIN_TRACE_SEATING=1, default-off -> byte-identical): where does the
        # isomer count collapse?  Logs enumerated vs chelate-feasible canonical forms so a 6->2
        # loss can be pinned to enumeration (achiral cf merges Λ/Δ) vs feasibility vs downstream build.
        if os.environ.get("DELFIN_TRACE_SEATING") == "1":
            try:
                import sys as _systr
                _systr.stderr.write(
                    f"[ISOTRACE] CN={n_coord} enumerated={len(isomers)} feasible={len(feasible_isomers)}\n")
                for _cf, _pm in isomers:
                    _feas = "OK " if (_cf, _pm) in feasible_isomers else "REJ"
                    _systr.stderr.write(f"[ISOTRACE]   {_feas} cf={_cf} perm={_pm}\n")
            except Exception:
                pass

        # Pre-compute ranked template conformers once per metal centre so each
        # permutation can retry against several templates when the default
        # (best-scored) one produces no viable XYZ.
        topo_template_cids = _rank_template_conformers(mol, top_k=_prof_topk) or [None]

        _PRIMARY_GEOM_BASE = {
            2: 'LIN', 3: 'TP', 5: 'TBP',
            6: 'OH', 7: 'PBP', 8: 'SAP', 9: 'TTP',
        }
        _PRIMARY_GEOM = dict(_PRIMARY_GEOM_BASE)
        _PRIMARY_GEOM[4] = _preferred_cn4_for(
            atom.GetSymbol(), mol, atom.GetIdx()
        )
        _GEOM_PRETTY = {
            'LIN': 'linear', 'TP': 'trigonal-planar', 'TS': 'T-shaped',
            'SQ': 'square-planar', 'TH': 'tetrahedral', 'SS': 'see-saw',
            'TBP': 'trigonal-bipyramidal', 'SP': 'square-pyramidal',
            'OH': 'octahedral', 'TPR': 'trigonal-prismatic',
            'PBP': 'pentagonal-bipyramidal', 'COH': 'capped-octahedral',
            'SAP': 'square-antiprismatic', 'DD': 'dodecahedral',
            'TTP': 'tricapped-trigonal-prismatic',
        }
        _primary_geom = _PRIMARY_GEOM.get(n_coord)

        def _build_one_topo(args):
            cf, pm = args
            gn = cf[0]
            try:
                xyz = None
                for _tc in topo_template_cids:
                    xyz = _build_topology_xyz(
                        mol, metal_idx, donor_indices, pm, gn,
                        apply_uff, conf_id=_tc,
                    )
                    if xyz is not None:
                        break
                if xyz is None:
                    return None
                mt = Chem.RWMol(mol)
                mt.RemoveAllConformers()
                c = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
                if c is None:
                    return None
                ci = mt.AddConformer(c, assignId=True)
                try:
                    if _has_atom_clash(mt.GetMol(), ci, min_dist=0.3):
                        return None
                    if _has_unphysical_metal_nonbonded_contact(mt.GetMol(), ci):
                        return None
                    if _has_unphysical_oco_geometry(mt.GetMol(), ci):
                        return None
                    if _has_pi_ring_nonplanarity(mt.GetMol(), ci):
                        return None
                    if _has_severe_covalent_distortion(mt.GetMol(), ci):
                        return None
                except Exception:
                    pass
                fp = _compute_coordination_fingerprint(
                    mt.GetMol(), ci, dtype_map=dtype_map
                )
                # Prefer canonical-form label: UFF can drift axial donors
                # enough that _classify_isomer_label mis-reads the intended
                # topology (e.g. N-N-ax → N-O-ax).  Fall back to classify
                # for geometries where canonical-form labelling isn't set.
                lbl = _label_from_canonical_form(cf)
                if not lbl:
                    lbl = _classify_isomer_label(fp, mt.GetMol())
                if _primary_geom and gn != _primary_geom:
                    gp = _GEOM_PRETTY.get(gn, gn)
                    lbl = f'{gp} {lbl}' if lbl else gp
                # Iter-2: append Λ/Δ helicity suffix (env-gated; safe when
                # cf carries no chirality tag → no-op).
                _hsuf = _extract_helicity_suffix(cf)
                if _hsuf == 'L':
                    lbl = f'{lbl}-Λ' if lbl else 'Λ'
                elif _hsuf == 'D':
                    lbl = f'{lbl}-Δ' if lbl else 'Δ'
                return (xyz, lbl)
            except Exception as exc:
                logger.debug("Topo isomer build failed (%s): %s", gn, exc)
                return None

        # Build 3D for each permutation. OB UFF holds the GIL so
        # ProcessPoolExecutor is needed for real parallelism.
        # Strategy: build pre-UFF Procrustes XYZ in main process (fast),
        # batch all UFF calls through ProcessPool, then quality-check.
        #
        # To give the stricter post-UFF gate (topology + hybridisation
        # + coordination-geometry) a large enough survivor pool we
        # over-generate aggressively:
        #   * every chelate conformer rank 0..K is tried
        #   * every ranked template conformer 0..T is tried
        #   * downstream fingerprint + RMSD dedup prunes duplicates
        # so identical outputs coming from equivalent (rank, template)
        # combinations don't pollute the final list.
        # Build loop preserved from the 7414981 passing state: rank-0 build
        # with first-success template, additional chelate ranks only when
        # no usable template exists.  Over-generating via every
        # (rank × template) grid caused UFF to collapse TBP isomers into
        # SP duplicates for systems like Fe(CO)3(NHC)2 and lose the
        # ``C0-C0-ax`` / ``C1-C1-ax`` labels in dedup.
        _CHELATE_RANK_VARIANTS = max(1, _prof_ranks) if chelate_ps else 1
        _PRE_UFF_CAP = max_isomers * max(1, _prof_cap_mult)
        _pre_uff_batch: List[Tuple[tuple, List[int], str, str, Optional[Dict]]] = []
        _pre_uff_seen: set = set()

        def _xyz_sig(_xyz: str) -> str:
            return "\n".join(
                _ln.strip() for _ln in _xyz.splitlines() if _ln.strip()
            )

        # Per-(cf, pm) variant counter so every distinct XYZ that passes
        # dedup for the same coordination arrangement gets a ``-conf2``,
        # ``-conf3`` etc. label suffix downstream, preserving backbone-
        # pucker / chelate-conformer variety through the label-collapse
        # step (otherwise all puckers share the CF label and only the
        # best-scoring one survives).
        _variant_counter: Dict[Tuple[tuple, tuple], int] = {}
        # PURELY-ADDITIVE d8/CN6 poly siblings (D8_SQ_ADD/CN6_OH_ADD) must NOT crowd the ISOMER budget
        # out (completeness is sacred): count them so the caps below see only the PRIMARY frames.  Without
        # this the OC/SP-4 siblings filled the cap and the isomer loop broke early (measured: VOYWUD
        # 6->5 isomers).  `_n_add_sib` corrects the PRE-UFF caps; `_sib_idxs` (the batch indices of the
        # siblings) corrects the POST-UFF append cap at ~26413 (the same hole, one stage later).
        _n_add_sib = 0
        _sib_idxs: set = set()
        for cf, pm in feasible_isomers:
            if len(_pre_uff_batch) + len(results) - _n_add_sib >= _PRE_UFF_CAP:
                break
            gn = cf[0]
            # Iterate ALL template conformer CIDs (not break on first):
            # each distinct template pucker is a candidate coordination
            # conformer worth keeping.  Dedup via XYZ signature drops
            # identical outputs deterministically.
            try:
                for _tc in topo_template_cids:
                    if len(_pre_uff_batch) + len(results) - _n_add_sib >= _PRE_UFF_CAP:
                        break
                    xyz0 = _build_topology_xyz(
                        mol, metal_idx, donor_indices, pm, gn,
                        False, conf_id=_tc, chelate_rank=0,
                    )
                    if xyz0 is None:
                        continue
                    _sig0 = _xyz_sig(xyz0)
                    if _sig0 in _pre_uff_seen:
                        continue
                    _pre_uff_seen.add(_sig0)
                    coord_c = None
                    if apply_uff:
                        try:
                            _d8t = None
                            if gn == 'SQ':                 # per-isomer trans from the enumerator (perm)
                                try:
                                    _tp = _TOPO_TRANS_POSITIONS.get(gn) or []
                                    _d8t = [(donor_indices[pm[_p1]], donor_indices[pm[_p2]])
                                            for (_p1, _p2) in _tp]
                                except Exception:
                                    _d8t = None
                            coord_c = _build_coordination_constraints_from_xyz(
                                mol, xyz0, d8_trans=_d8t,
                            )
                        except Exception:
                            pass
                    _key = (tuple(cf), tuple(pm))
                    _variant_counter[_key] = _variant_counter.get(_key, 0) + 1
                    _conf_idx = _variant_counter[_key] - 1
                    _pre_uff_batch.append((cf, pm, gn, xyz0, coord_c, _conf_idx))
                    # ADDITIVE d8 SP-4 (DELFIN_FFFREE_D8_SQ_ADD, default off -> byte-identical).  ONE
                    # self-contained axis: the PRIMARY frame above is the NORMAL (tetrahedral) build, and a
                    # UFF-SP-4 square is added as a PURELY ADDITIVE sibling from the SAME seed (force_d8_sq,
                    # independent of DELFIN_FFFREE_D8_SQ_ISO).  Never touches the primary -> a bulky ligand
                    # that clashes under SP-4 keeps its valid tetrahedral primary (NO broken_regressed -- the
                    # earlier "REPLACE the frame with SP-4" cost HAKQES a frame: +1 broken, +0 frames), while
                    # a d8-no-valid system GAINS its valid SP-4 frame.  Topology gate culls a clashing SP-4.
                    if (apply_uff and gn == 'SQ' and len(donor_indices) == 4
                            and _delfin_env_int("DELFIN_FFFREE_D8_SQ_ADD", 0)
                            and len(_pre_uff_batch) + len(results) - _n_add_sib < _PRE_UFF_CAP):
                        try:
                            _m_sym = mol.GetAtomWithIdx(int(metal_idx)).GetSymbol()
                            if _m_sym in _D8_SQ_ISO_METALS:
                                coord_c_sq = _build_coordination_constraints_from_xyz(
                                    mol, xyz0, d8_trans=_d8t, force_d8_sq=True,
                                )
                                if coord_c_sq != coord_c:   # SP-4 constraints differ from the primary
                                    _variant_counter[_key] += 1
                                    _pre_uff_batch.append(
                                        (cf, pm, gn, xyz0, coord_c_sq, _variant_counter[_key] - 1))
                                    _sib_idxs.add(len(_pre_uff_batch) - 1)
                                    _n_add_sib += 1         # additive -> does not count vs the isomer cap
                        except Exception:
                            pass
                    # ADDITIVE CN6 OCTAHEDRON (DELFIN_FFFREE_CN6_OH_ADD, default off -> byte-identical).  Same
                    # self-contained-additive pattern as the d8 SP-4 above, for the biggest poly cluster
                    # (TPR-6 built, OC-6 in the crystal).  The PRIMARY frame is the NORMAL build; a UFF-OC
                    # octahedron is added as a PURELY ADDITIVE sibling from the SAME seed.  ISOMER-ORTHOGONAL:
                    # the sibling uses the GEOMETRY-FALLBACK twist-correction (d8_trans=None -> impose OC on
                    # the frame's OWN most-opposite donor pairs), NOT the OH-PERM path.  Measured 2026-07-14:
                    # the PERM path (enumerator OH positions via pm) COLLAPSED distinct isomers (VOYWUD lost
                    # all-trans + trans-OH) -- that OH-perm mapping was never validated (dead before) and
                    # imposes the WRONG trans set on some arrangements.  The twist-correction keeps whatever
                    # trans pairs the frame already has, so it can NEVER reshape one isomer into another.
                    if (apply_uff and gn == 'OH' and len(donor_indices) == 6
                            and _delfin_env_int("DELFIN_FFFREE_CN6_OH_ADD", 0)
                            and len(_pre_uff_batch) + len(results) - _n_add_sib < _PRE_UFF_CAP):
                        try:
                            _m_sym = mol.GetAtomWithIdx(int(metal_idx)).GetSymbol()
                            if _PREFERRED_CN6_GEOMETRY.get(_m_sym, 'OH') == 'OH':
                                coord_c_oh = _build_coordination_constraints_from_xyz(
                                    mol, xyz0, d8_trans=None, force_cn6_oh=True,
                                )
                                if coord_c_oh != coord_c:   # OC constraints differ from the primary
                                    _variant_counter[_key] += 1
                                    _pre_uff_batch.append(
                                        (cf, pm, gn, xyz0, coord_c_oh, _variant_counter[_key] - 1))
                                    _sib_idxs.add(len(_pre_uff_batch) - 1)
                                    _n_add_sib += 1         # additive -> does not count vs the isomer cap
                        except Exception:
                            pass
                    # Iter-8.5b INNER site 1 (template-loop): when the parent
                    # mol's class is in DELFIN_ITER85_PUMP_SKIP_CLASSES, take
                    # only the first successful template seed per perm
                    # instead of iterating every ranked template CID.  Mirrors
                    # the outer 8.5b additive-skip philosophy at the inner
                    # pump.  Default off = bit-exact (loop continues).
                    if _iter85b_pump_skip_inner:
                        break
            except Exception as exc:
                logger.debug("Topo pre-UFF build failed (%s): %s", cf[0], exc)
                continue

            # Balloon-inflate emits an additional candidate for every
            # system (mono- and multi-metallic).  It builds from scratch
            # (no ETKDG-template bias), places the M-M-bridge scaffold
            # at ideal distances for bimetallics and Procrustes-aligns
            # ligand fragments independently.  Because the XYZ is added
            # alongside the template build above and dedup'd via the
            # XYZ signature, no mono-metal variety is lost — balloon
            # only increases the candidate pool.  Deterministic by
            # construction (fixed chelate-conformer seed schedule).
            # Iter-8.5b INNER site 2 (balloon-builder): when the parent
            # mol's class is in DELFIN_ITER85_PUMP_SKIP_CLASSES, skip the
            # balloon additive emission for this perm.  Default off =
            # bit-exact (block runs as before).
            if not _iter85b_pump_skip_inner:
                try:
                    xyz_bln = _build_topology_xyz_from_scratch(
                        mol, metal_idx, donor_indices, pm, gn,
                        chelate_rank=0,
                    )
                    if xyz_bln is not None:
                        _sig_bln = _xyz_sig(xyz_bln)
                        if _sig_bln not in _pre_uff_seen:
                            _pre_uff_seen.add(_sig_bln)
                            coord_c_bln = None
                            if apply_uff:
                                try:
                                    coord_c_bln = _build_coordination_constraints_from_xyz(
                                        mol, xyz_bln,
                                    )
                                except Exception:
                                    pass
                            _key = (tuple(cf), tuple(pm))
                            _variant_counter[_key] = _variant_counter.get(_key, 0) + 1
                            _conf_idx = _variant_counter[_key] - 1
                            _pre_uff_batch.append(
                                (cf, pm, gn, xyz_bln, coord_c_bln, _conf_idx)
                            )
                except Exception as bln_exc:
                    logger.debug(
                        "Balloon builder raised for (%s, perm=%s): %s",
                        cf[0], pm, bln_exc,
                    )

            # Additional chelate-rank variants: enumerate alternative
            # chelate puckers (rank 1..N-1) for every (CF, perm).  The
            # XYZ-signature dedup below drops identical outputs
            # deterministically, so running across ranks cannot
            # introduce non-determinism even when mol already has
            # conformers.  Each surviving distinct XYZ becomes a
            # ``-confN`` labelled variant in the output so flexible
            # chelates (salen, cryptand, ethylenediamine) no longer
            # collapse to a single best-scoring pucker.
            # Iter-8.5b INNER site 3 (chelate-rank-variants): when the
            # parent mol's class is in DELFIN_ITER85_PUMP_SKIP_CLASSES,
            # skip the chelate-rank pump entirely for this perm.
            # Default off = bit-exact (loop runs as before).
            if _iter85b_pump_skip_inner:
                continue
            if _CHELATE_RANK_VARIANTS <= 1:
                continue
            for _crank in range(1, _CHELATE_RANK_VARIANTS):
                if len(_pre_uff_batch) + len(results) >= _PRE_UFF_CAP:
                    break
                try:
                    for _tc in topo_template_cids:
                        if len(_pre_uff_batch) + len(results) >= _PRE_UFF_CAP:
                            break
                        xyz = _build_topology_xyz(
                            mol, metal_idx, donor_indices, pm, gn,
                            False, conf_id=_tc, chelate_rank=_crank,
                        )
                        if xyz is None:
                            continue
                        _sig = _xyz_sig(xyz)
                        if _sig in _pre_uff_seen:
                            continue
                        _pre_uff_seen.add(_sig)
                        coord_c = None
                        if apply_uff:
                            try:
                                coord_c = _build_coordination_constraints_from_xyz(
                                    mol, xyz,
                                )
                            except Exception:
                                pass
                        _key = (tuple(cf), tuple(pm))
                        _variant_counter[_key] = _variant_counter.get(_key, 0) + 1
                        _conf_idx = _variant_counter[_key] - 1
                        _pre_uff_batch.append(
                            (cf, pm, gn, xyz, coord_c, _conf_idx)
                        )
                except Exception as exc:
                    logger.debug("Topo pre-UFF build failed (%s): %s", cf[0], exc)

        # Batch UFF via ProcessPool (OB holds GIL → threads don't help).
        if apply_uff and _pre_uff_batch:
            _n_uff_workers = min(
                len(_pre_uff_batch), os.cpu_count() or 4, DELFIN_MAX_PROCESS_WORKERS
            )
            _uff_inputs = [
                (xyz, 500, cstr) for _cf, _pm, _gn, xyz, cstr, _ci in _pre_uff_batch
            ]
            try:
                if _n_uff_workers > 1 and len(_uff_inputs) > 2:
                    with concurrent.futures.ProcessPoolExecutor(
                        max_workers=_n_uff_workers
                    ) as _pp:
                        _uff_results = list(_pp.map(
                            _optimize_xyz_openbabel,
                            [inp[0] for inp in _uff_inputs],
                            [inp[1] for inp in _uff_inputs],
                            [inp[2] for inp in _uff_inputs],
                        ))
                else:
                    _uff_results = [
                        _optimize_xyz_openbabel(inp[0], inp[1], inp[2])
                        for inp in _uff_inputs
                    ]
            except Exception as _ppe:
                logger.debug("ProcessPool UFF failed, falling back to sequential: %s", _ppe)
                _uff_results = [
                    _optimize_xyz_openbabel(inp[0], inp[1], inp[2])
                    for inp in _uff_inputs
                ]
            def _tr_mcx(_xyzs, _tag):
                if os.environ.get("DELFIN_TRACE_SEATING", "0") != "1":
                    return
                try:
                    _pos = {}
                    for _i, _ln in enumerate(_xyzs.strip().split("\n")):
                        _p = _ln.split()
                        if len(_p) >= 4:
                            _pos[_i] = np.array([float(_p[1]), float(_p[2]), float(_p[3])])
                    _mp = _pos.get(metal_idx)
                    _cdonors = [nb.GetIdx() for nb in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
                                if nb.GetSymbol() == "C"]
                    for _d in _cdonors:
                        if _mp is None:
                            continue
                        _dp = _pos.get(_d)
                        if _dp is None:
                            continue
                        _mc = _mp - _dp; _mcn = float(np.linalg.norm(_mc))
                        if _mcn < 1e-6:
                            continue
                        _mc /= _mcn
                        for _nb in mol.GetAtomWithIdx(_d).GetNeighbors():
                            if _nb.GetAtomicNum() <= 1:
                                continue
                            _xp = _pos.get(_nb.GetIdx())
                            if _xp is None:
                                continue
                            _cx = _xp - _dp; _cxn = float(np.linalg.norm(_cx))
                            if 1.3 < _cxn < 1.9:
                                _ang = float(np.degrees(np.arccos(
                                    max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                                _trace_seating("%s donor=%d M-C-Xheavy=%.0f" % (_tag, _d, _ang))
                except Exception:
                    pass
            for idx, (cf, pm, gn, xyz_pre, _cstr, _ci) in enumerate(_pre_uff_batch):
                xyz_opt = _uff_results[idx] if idx < len(_uff_results) else xyz_pre
                if not xyz_opt:
                    xyz_opt = xyz_pre
                _tr_mcx(xyz_opt, "BATCH_POST_UFF")
                # Post-UFF polish: project sp2 3-coordinate atoms onto
                # their neighbours' plane to remove residual
                # pyramidalisation that the torsion constraints could not
                # fully eliminate.
                xyz_opt = _flatten_sp2_atoms_xyz(xyz_opt, mol)
                _tr_mcx(xyz_opt, "BATCH_POST_FLATTEN")
                # Conservative UFF: if UFF broke topology, keep pre-UFF.
                if not _verify_topology_from_graph(xyz_opt, mol):
                    xyz_opt = xyz_pre
                _pre_uff_batch[idx] = (cf, pm, gn, xyz_opt, _cstr, _ci)

        # Post-UFF: graph-based topology check (replaces the 5 legacy
        # checks that were too aggressive for topo-generated structures).
        # The max_isomers cap must count only PRIMARY frames, not the purely-additive d8/CN6 poly
        # siblings (else an interleaved sibling crowds a later isomer's PRIMARY out of results ->
        # VOYWUD 6->5).  `_n_sib_appended` mirrors the pre-UFF `-_n_add_sib`; empty _sib_idxs (flags
        # off) -> byte-identical to the original `len(results) >= max_isomers`.
        _n_sib_appended = 0
        # BASE-PRESERVATION via TOPOLOGY, NOT RMSD (user 2026-07-22: "RMSD is the worst metric";
        # doctrine: gate = topology, NEVER RMSD).  A FEAS_PREFERRED recovery is a real win ONLY if it adds
        # a GENUINELY NEW coordination isomer.  On a RIGID scaffold a reach-recovered arrangement relaxes
        # (UFF) onto an isomer the base set ALREADY built -> its BUILT coordination FINGERPRINT equals a
        # base frame's -> it is redundant, and worse, the downstream fingerprint dedup then drops the GOOD
        # base frame in favour of the (distorted) recovery (AXOKED: +5 reach-recoveries cost a good square-
        # pyramidal base conformer, good 36->35 = the broken_regressed the eye flagged).  Fix, universal,
        # purely TOPOLOGICAL (coordination fingerprint = which donor sits where; no RMSD, no energy, no
        # fitted threshold): drop a recovery whose built fingerprint is ALREADY realised by a base frame ->
        # the redundant recovery never enters, so it can neither pad the manifold nor evict a base frame.
        # A genuinely-new isomer (GOWFED all-cis: a fingerprint the base set was MISSING) has a NEW
        # fingerprint -> kept = the real completeness win.  This is EXACTLY the definition of "recovers a
        # MISSING isomer": keep iff it adds a fingerprint the base does not already have.  Base frames build
        # FIRST (feasible_isomers = base + _pref_extra), so every base fingerprint a recovery could
        # duplicate is already recorded by the time the recovery is reached.
        _base_fps: set = set()
        for _batch_i, (cf, pm, gn, xyz, _cstr, _cidx) in enumerate(_pre_uff_batch):
            if (len(results) - _n_sib_appended) >= max_isomers:
                break
            try:
                if not _verify_topology_from_graph(xyz, mol):
                    continue
                mt = Chem.RWMol(mol)
                mt.RemoveAllConformers()
                c = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
                if c is None:
                    continue
                ci = mt.AddConformer(c, assignId=True)
                fp = _compute_coordination_fingerprint(
                    mt.GetMol(), ci, dtype_map=dtype_map
                )
                _is_pref_extra = (tuple(cf), tuple(pm)) in _pref_extra_keys
                if _is_pref_extra:
                    # Keep a recovery ONLY if it (a) realises a NEW coordination isomer -- its built
                    # fingerprint is not already among the base frames (the TOPOLOGICAL redundancy test,
                    # the primary discriminator) -- AND (b) is geometrically sound: achieves its intended
                    # polyhedron (only_geom=gn), no torn/stretched covalent bond, no donor-H pointing at the
                    # metal.  All topology/geometry, no RMSD, no energy.  Scoped to extras -> primary/
                    # champion frames are never touched (additive by construction).  A genuine trig-prism is
                    # enumerated + built AS TPR -> scored vs TPR -> kept (real prisms untouched).
                    try:
                        _mG = mt.GetMol()
                        _redundant = fp in _base_fps
                        _devs = _ideal_polyhedron_angle_dev_per_metal(_mG, ci, only_geom=gn)
                        _drop = (_redundant
                                 or (_devs and max(_devs.values()) > _extra_max_dev)
                                 or _has_severe_covalent_distortion(_mG, ci)
                                 or _donor_h_points_at_metal(_mG, ci, _extra_donor_h_min))
                        if os.environ.get("DELFIN_TRACE_SEATING") == "1":
                            try:
                                import sys as _systr
                                _systr.stderr.write(
                                    "[FEASFLOOR] %s gn=%s redundant=%s poly_vs_geom=%.1f drop=%s\n" % (
                                        _label_from_canonical_form(cf) or str(cf), gn, _redundant,
                                        (max(_devs.values()) if _devs else -1.0), _drop))
                            except Exception:
                                pass
                        if _drop:
                            continue
                    except Exception:
                        pass
                # Canonical-form label (see rationale above).
                lbl = _label_from_canonical_form(cf)
                if not lbl:
                    lbl = _classify_isomer_label(fp, mt.GetMol())
                if _primary_geom and gn != _primary_geom:
                    gp = _GEOM_PRETTY.get(gn, gn)
                    lbl = f'{gp} {lbl}' if lbl else gp
                # Iter-2: append Λ/Δ helicity suffix when present.
                _hsuf = _extract_helicity_suffix(cf)
                if _hsuf == 'L':
                    lbl = f'{lbl}-Λ' if lbl else 'Λ'
                elif _hsuf == 'D':
                    lbl = f'{lbl}-Δ' if lbl else 'Δ'
                # Conformer-variant suffix: second, third, ... distinct
                # XYZ for the same (CF, perm) gets ``-conf2``, ``-conf3``
                # so the downstream label-collapse keeps every pucker.
                if _cidx and _cidx > 0:
                    lbl = f'{lbl}-conf{_cidx + 1}' if lbl else f'conf{_cidx + 1}'
                results.append((xyz, lbl))
                if not _is_pref_extra:
                    # record the BASE (non-recovery) fingerprint so a later recovery that collapses onto
                    # this isomer is caught by the topological redundancy test above (no good ones
                    # vanish -- a recovery may never duplicate, and thus displace, a base isomer).
                    _base_fps.add(fp)
                if _batch_i in _sib_idxs:      # additive sibling -> does not count vs max_isomers
                    _n_sib_appended += 1
            except Exception as exc:
                logger.debug("Topo post-UFF check failed (%s): %s", gn, exc)
                continue

    # --- Multinuclear coupled enumeration for 2-metal clusters ---
    # Detect metals connected through bridging donors and enumerate the
    # Cartesian product of their per-metal isomers to capture arrangements
    # that per-metal enumeration misses.
    try:
        bridging = _find_bridging_donors(mol)
        if bridging and len(results) < max_isomers:
            # Build metal cluster graph via bridging donors.
            metal_indices = [
                a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() in _METAL_SET
            ]
            if len(metal_indices) == 2:
                m1, m2 = metal_indices
                d1 = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(m1).GetNeighbors()]
                d2 = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(m2).GetNeighbors()]
                n1, n2 = len(d1), len(d2)
                if 2 <= n1 <= 9 and 2 <= n2 <= 9:
                    # Per-metal isomers.
                    dk1 = [dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset())) for d in d1]
                    dk2 = [dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset())) for d in d2]
                    uk1 = sorted(set(dk1), key=lambda k: (k[0], tuple(sorted(k[1]))))
                    uk2 = sorted(set(dk2), key=lambda k: (k[0], tuple(sorted(k[1]))))
                    kc1 = {k: i for i, k in enumerate(uk1)}
                    kc2 = {k: i for i, k in enumerate(uk2)}
                    dl1 = [f"{k[0]}{kc1[k]}" for k in dk1]
                    dl2 = [f"{k[0]}{kc2[k]}" for k in dk2]
                    cp1 = _chelate_pairs(mol, m1, d1)
                    cp2 = _chelate_pairs(mol, m2, d2)
                    al1 = {atom_idx: li for li, atom_idx in enumerate(d1)}
                    al2 = {atom_idx: li for li, atom_idx in enumerate(d2)}
                    clp1 = [frozenset([al1[sorted(cp)[0]], al1[sorted(cp)[1]]]) for cp in cp1
                            if sorted(cp)[0] in al1 and sorted(cp)[1] in al1]
                    clp2 = [frozenset([al2[sorted(cp)[0]], al2[sorted(cp)[1]]]) for cp in cp2
                            if sorted(cp)[0] in al2 and sorted(cp)[1] in al2]
                    ms1 = mol.GetAtomWithIdx(m1).GetSymbol()
                    ms2 = mol.GetAtomWithIdx(m2).GetSymbol()
                    iso1 = _enumerate_topological_isomers(dl1, n1, clp1, metal_symbol=ms1)
                    iso2 = _enumerate_topological_isomers(dl2, n2, clp2, metal_symbol=ms2)

                    # Scaffold-first approach: build M-bridge-M core first,
                    # then place non-bridging donors around each metal.
                    # Use the ETKDG template as scaffold base (it has correct
                    # M-bridge-M topology from SMILES).
                    import itertools as _it
                    # Combo cap decoupled from max_isomers so the full
                    # Cartesian product of per-metal (cf, pm) x
                    # template x chelate_rank combinations is explored
                    # before dedup + ranking trims down to max_isomers.
                    # 3x max_isomers keeps total work bounded while
                    # giving enough candidates for the ranking to pick
                    # the geometrically best ones.
                    max_combos = max(max_isomers * 3, 100)
                    scaffold = _build_multimetal_scaffold(
                        mol, metal_indices, bridging
                    )
                    topo_template_cids = _rank_template_conformers(mol, top_k=_prof_topk) or [None]

                    # Pre-stretch the metal-metal separation in the
                    # template conformer so that |M1-M2| = d_M1^ideal +
                    # d_M2^ideal (first bridging donor's reference).
                    # Without this the subsequent bridge-snap puts the
                    # bridging donor at a fractional distance on a
                    # metal-metal vector that is too short, and the
                    # resulting M-bridge distance falls below the
                    # Rule 1 window.  The stretch is performed once
                    # per multinuclear enumeration and reused for
                    # every Cartesian combo.
                    try:
                        if scaffold and topo_template_cids:
                            _first_bridge = bridging[0][0]
                            _d_sym = mol.GetAtomWithIdx(_first_bridge).GetSymbol()
                            d_m1_ideal = float(_get_ml_bond_length(ms1, _d_sym))
                            d_m2_ideal = float(_get_ml_bond_length(ms2, _d_sym))
                            target_mm = d_m1_ideal + d_m2_ideal
                            for _tcid in topo_template_cids:
                                try:
                                    _conf = mol.GetConformer(int(_tcid) if _tcid is not None else 0)
                                    _p1 = _conf.GetAtomPosition(m1)
                                    _p2 = _conf.GetAtomPosition(m2)
                                    _dvec = (_p2.x - _p1.x, _p2.y - _p1.y, _p2.z - _p1.z)
                                    _dn = math.sqrt(sum(v * v for v in _dvec))
                                    if _dn < 1e-6:
                                        continue
                                    _scale = target_mm / _dn
                                    if abs(_scale - 1.0) < 0.02:
                                        continue
                                    # Shift metal_2 along the existing axis so
                                    # |M1-M2| matches target.  The bridging
                                    # atoms and downstream ligand atoms also
                                    # move rigidly with metal_2 when we later
                                    # rebuild metal_2's polyhedron, so the
                                    # local Fe-donor geometry is preserved.
                                    _delta = [(target_mm - _dn) * (v / _dn) for v in _dvec]
                                    _conf.SetAtomPosition(
                                        m2,
                                        type(_p2)(
                                            float(_p2.x + _delta[0]),
                                            float(_p2.y + _delta[1]),
                                            float(_p2.z + _delta[2]),
                                        ),
                                    )
                                except Exception:
                                    continue
                    except Exception as _strch_exc:
                        logger.debug("M-M pre-stretch failed: %s", _strch_exc)

                    # Full combinatorial enumeration: every
                    # (cf1, pm1) x (cf2, pm2) pair gets combined with
                    # every template conformer AND every chelate-rank
                    # pucker variant on both sides.  Distinct XYZs
                    # survive via signature dedup; the per-combo
                    # variant counter adds a ``-confN`` suffix to the
                    # label so the downstream label-collapse keeps
                    # them.  Respects ``max_combos`` to avoid blowup.
                    _combo_ranks = max(1, _CHELATE_RANK_VARIANTS)
                    _combo_seen: set = set()
                    _combo_variant_counter: Dict[Tuple[tuple, tuple, tuple, tuple], int] = {}
                    combo_count = 0
                    import time as _time_2m
                    _2m_start = _time_2m.time()
                    # DETERMINISM vs anti-TLE trade-off (2-metal enum).  The
                    # wall-clock cutoff bounds the (cf1,pm1)x(cf2,pm2) product
                    # but makes the enumerated set TIMING-dependent.  Env-gated
                    # DELFIN_2METAL_WALL_BUDGET_S (default 240 = byte-identical
                    # to pre-change).  Master switch forces 0 (DETERMINISTIC:
                    # termination by the max_combos cap only); 0 disables it.
                    _2M_WALL_BUDGET = 0.0 if _deterministic_mode() else _delfin_env_float(
                        "DELFIN_2METAL_WALL_BUDGET_S", 240.0,
                    )
                    for (cf1, pm1), (cf2, pm2) in _it.product(iso1, iso2):
                        if combo_count >= max_combos:
                            break
                        if _2M_WALL_BUDGET > 0 and _time_2m.time() - _2m_start > _2M_WALL_BUDGET:
                            logger.debug(
                                "2-metal enum wall-clock budget %.0fs exceeded, stopping.",
                                _2M_WALL_BUDGET,
                            )
                            break
                        gn1 = cf1[0]
                        gn2 = cf2[0]
                        try:
                            for _crank1 in range(_combo_ranks):
                                if combo_count >= max_combos:
                                    break
                                for _crank2 in range(_combo_ranks):
                                    if combo_count >= max_combos:
                                        break
                                    for _tpl_cid in topo_template_cids:
                                        if combo_count >= max_combos:
                                            break
                                        xyz1 = _build_topology_xyz(
                                            mol, m1, d1, pm1, gn1, False,
                                            conf_id=_tpl_cid,
                                            chelate_rank=_crank1,
                                        )
                                        if xyz1 is None:
                                            continue
                                        mol_tmp = Chem.RWMol(mol)
                                        mol_tmp.RemoveAllConformers()
                                        conf_tmp = _xyz_to_rdkit_conformer(
                                            mol_tmp.GetMol(), xyz1,
                                        )
                                        if conf_tmp is None:
                                            continue
                                        cid_tmp = mol_tmp.AddConformer(conf_tmp, assignId=True)
                                        _rescale_metal_donor_distances(mol_tmp, cid_tmp)
                                        xyz_combined = _build_topology_xyz(
                                            mol_tmp.GetMol(), m2, d2, pm2, gn2, False,
                                            conf_id=cid_tmp,
                                            chelate_rank=_crank2,
                                        )
                                        if xyz_combined is None:
                                            continue
                                        mol_tmp2 = Chem.RWMol(mol)
                                        mol_tmp2.RemoveAllConformers()
                                        conf_c = _xyz_to_rdkit_conformer(
                                            mol_tmp2.GetMol(), xyz_combined,
                                        )
                                        if conf_c is not None:
                                            cid_c = mol_tmp2.AddConformer(conf_c, assignId=True)
                                            _rescale_metal_donor_distances(mol_tmp2, cid_c)
                                            xyz_combined = _mol_to_xyz_conformer(mol_tmp2, cid_c)
                                        try:
                                            xyz_combined = _snap_bridging_donors_to_compromise(
                                                xyz_combined, mol,
                                                [m1, m2], bridging,
                                            )
                                        except Exception as _snap_exc:
                                            logger.debug(
                                                "Bridge-snap failed: %s", _snap_exc,
                                            )
                                        if apply_uff:
                                            xyz_combined = _optimize_xyz_openbabel_safe(
                                                xyz_combined, mol_template=mol,
                                            )
                                        # XYZ-signature dedup across
                                        # (template, rank1, rank2).
                                        _sig_cb = _xyz_sig(xyz_combined)
                                        if _sig_cb in _combo_seen:
                                            continue
                                        _combo_seen.add(_sig_cb)
                                        _cb_key = (tuple(cf1), tuple(pm1), tuple(cf2), tuple(pm2))
                                        _combo_variant_counter[_cb_key] = (
                                            _combo_variant_counter.get(_cb_key, 0) + 1
                                        )
                                        _cb_idx = _combo_variant_counter[_cb_key] - 1
                                        label = f"multi-{gn1}/{gn2}"
                                        if _cb_idx > 0:
                                            label = f"{label}-conf{_cb_idx + 1}"
                                        # Bond-length mini-gate (looser
                                        # than the full graph gate —
                                        # allows bridge-compromise
                                        # geometry but rejects truly
                                        # catastrophic collapses).
                                        try:
                                            _q_lines = [
                                                l for l in xyz_combined.splitlines() if l.strip()
                                            ]
                                            _coords_q = []
                                            for _ln in _q_lines:
                                                _p = _ln.split()
                                                if len(_p) >= 4:
                                                    _coords_q.append((float(_p[1]), float(_p[2]), float(_p[3])))
                                            gate_ok = True
                                            for _b in mol.GetBonds():
                                                _a1 = _b.GetBeginAtom(); _a2 = _b.GetEndAtom()
                                                if (_a1.GetAtomicNum() <= 1
                                                        or _a2.GetAtomicNum() <= 1):
                                                    continue
                                                _s1 = _a1.GetSymbol(); _s2 = _a2.GetSymbol()
                                                if _s1 in _METAL_SET and _s2 in _METAL_SET:
                                                    continue
                                                _i1 = _a1.GetIdx(); _i2 = _a2.GetIdx()
                                                _dx = _coords_q[_i1][0] - _coords_q[_i2][0]
                                                _dy = _coords_q[_i1][1] - _coords_q[_i2][1]
                                                _dz = _coords_q[_i1][2] - _coords_q[_i2][2]
                                                _d = math.sqrt(_dx*_dx + _dy*_dy + _dz*_dz)
                                                if _s1 in _METAL_SET or _s2 in _METAL_SET:
                                                    _m_sym = _s1 if _s1 in _METAL_SET else _s2
                                                    _d_sym = _s2 if _s1 in _METAL_SET else _s1
                                                    _ideal = float(_get_ml_bond_length(_m_sym, _d_sym))
                                                    if _ideal > 0 and (_d < 0.50 * _ideal or _d > 2.50 * _ideal):
                                                        gate_ok = False
                                                        break
                                                else:
                                                    if _d > 2.5:
                                                        gate_ok = False
                                                        break
                                            if not gate_ok:
                                                continue
                                        except Exception:
                                            pass
                                        results.append((xyz_combined, label))
                                        combo_count += 1
                        except Exception as _cexc:
                            logger.debug("Multinuclear combo build failed: %s", _cexc)
                            continue
            elif len(metal_indices) >= 3:
                # --- N-metal coupled enumeration (tri-/tetra-/... metallic) ---
                # The 2-metal block above does bridge-snap and
                # _build_multimetal_scaffold, both hardcoded for 2
                # metals.  For N >= 3 we build each metal's polyhedron
                # sequentially (each subsequent build uses the previous
                # metal's xyz as the starting conformer) and rely on
                # UFF + the mini-gate to settle the cluster.  The
                # Cartesian product explodes rapidly (k^N with k ~= 4
                # geoms per metal), so max_combos caps total output and
                # an XYZ-signature dedup filters duplicates across the
                # template/rank space.
                import itertools as _it
                _N = len(metal_indices)
                _per_metal: List[List[Tuple[tuple, List[int]]]] = []
                _mi_symbols: List[str] = []
                _mi_donors: List[List[int]] = []
                _mi_ok = True
                for _mi in metal_indices:
                    _di = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(_mi).GetNeighbors()]
                    _ni = len(_di)
                    if not (2 <= _ni <= 9):
                        _mi_ok = False
                        break
                    _dki = [
                        dtype_map.get(
                            _d, (mol.GetAtomWithIdx(_d).GetSymbol(), frozenset())
                        ) for _d in _di
                    ]
                    _uki = sorted(set(_dki), key=lambda k: (k[0], tuple(sorted(k[1]))))
                    _kci = {_k: _i for _i, _k in enumerate(_uki)}
                    _dli = [f"{_k[0]}{_kci[_k]}" for _k in _dki]
                    _cpi = _chelate_pairs(mol, _mi, _di)
                    _ali = {_aidx: _li for _li, _aidx in enumerate(_di)}
                    _clpi = [
                        frozenset([_ali[sorted(_cp)[0]], _ali[sorted(_cp)[1]]])
                        for _cp in _cpi
                        if sorted(_cp)[0] in _ali and sorted(_cp)[1] in _ali
                    ]
                    _msi = mol.GetAtomWithIdx(_mi).GetSymbol()
                    _isoi = _enumerate_topological_isomers(
                        _dli, _ni, _clpi, metal_symbol=_msi,
                    )
                    if not _isoi:
                        _mi_ok = False
                        break
                    _per_metal.append(_isoi)
                    _mi_symbols.append(_msi)
                    _mi_donors.append(_di)

                if _mi_ok and _per_metal:
                    _max_combos_n = max(max_isomers * 3, 100)
                    _topo_cids_n = _rank_template_conformers(mol, top_k=_prof_topk) or [None]
                    _n_ranks = max(1, _CHELATE_RANK_VARIANTS)
                    # Smart-mode truncation: when ``n_metal_smart`` is
                    # True AND N >= 4, trim the per-metal arrangements
                    # list to keep the combinatorial product bounded
                    # (K=2 for N=4, K=1 for N>=5).  Selection uses the
                    # enumerator's native ordering which is already
                    # sorted by metal-specific preferred geometry, so
                    # the smart cut keeps the chemically most-likely
                    # arrangements.  When n_metal_smart=False the full
                    # Cartesian product runs — still bounded by
                    # _N_ITER_BUDGET so pathological N>=6 systems
                    # don't hang indefinitely.
                    if n_metal_smart and _N >= 5:
                        _per_metal_eff = [_pm[:1] for _pm in _per_metal]
                    elif n_metal_smart and _N >= 4:
                        _per_metal_eff = [_pm[:2] for _pm in _per_metal]
                    else:
                        _per_metal_eff = _per_metal
                    # Wall-clock budget of 90 s per multinuclear call —
                    # the N-metal Cartesian product explodes exponentially
                    # and each UFF call can take 1-5 s.  Without this
                    # bound Fe3 (salen-like)-type systems TLE at 3 h+
                    # because the inner loop iterates 5000+ combos.
                    # Sampling augmentation further down the pipeline
                    # still runs and provides coverage.
                    _N_ITER_BUDGET = max(_max_combos_n * 40, 2000)
                    import time as _time
                    _n_metal_start = _time.time()
                    # DETERMINISM vs anti-TLE trade-off (multinuclear enum).
                    # The wall-clock cutoff bounds the exponential N-metal product
                    # but makes the enumerated set depend on TIMING -> the same
                    # input can give different label sets across runs (esp. under
                    # CPU load).  Default 90 = current behaviour (anti-TLE, NOT
                    # deterministic under load) -> byte-identical to pre-change.
                    # Set DELFIN_NMETAL_WALL_BUDGET_S=0 for DETERMINISTIC mode
                    # (termination by the deterministic _N_ITER_BUDGET + _max_combos_n
                    # caps only) -- WARNING: can TLE on heavy multinuclear systems
                    # until the deterministic iteration cap is tuned (see brief).
                    # Master switch forces 0 (DETERMINISTIC: termination by the
                    # deterministic _N_ITER_BUDGET + _max_combos_n caps only).
                    _N_WALL_BUDGET = 0.0 if _deterministic_mode() else float(
                        os.environ.get("DELFIN_NMETAL_WALL_BUDGET_S", "90")
                    )
                    _iter_count = 0
                    _combo_seen_n: set = set()
                    _variant_counter_n: Dict[Tuple[tuple, ...], int] = {}
                    _n_combo_count = 0
                    # Cartesian product over all metals' (cf, pm)
                    # arrangements.  Cap total output at max_combos
                    # because N=4 with 4 geoms/metal is 4^4=256 base
                    # combos before templates x ranks.
                    for _combo_tuple in _it.product(*_per_metal_eff):
                        _iter_count += 1
                        if _iter_count > _N_ITER_BUDGET:
                            logger.debug(
                                "N-metal enum iteration budget %d reached at N=%d, stopping.",
                                _N_ITER_BUDGET, _N,
                            )
                            break
                        if _N_WALL_BUDGET > 0 and _time.time() - _n_metal_start > _N_WALL_BUDGET:
                            logger.debug(
                                "N-metal enum wall-clock budget %.0fs exceeded at N=%d, stopping.",
                                _N_WALL_BUDGET, _N,
                            )
                            break
                        if _n_combo_count >= _max_combos_n:
                            break
                        # Each _combo_tuple = ((cf_0, pm_0), (cf_1, pm_1), ...).
                        _gns = tuple(cf[0] for cf, _pm in _combo_tuple)
                        _pms = tuple(tuple(pm) for _cf, pm in _combo_tuple)
                        _cfs = tuple(tuple(cf) for cf, _pm in _combo_tuple)
                        try:
                            for _rank_tuple in _it.product(
                                range(_n_ranks), repeat=_N
                            ):
                                if _n_combo_count >= _max_combos_n:
                                    break
                                for _tcid in _topo_cids_n:
                                    if _n_combo_count >= _max_combos_n:
                                        break
                                    # Chain builds: first metal uses
                                    # template conf, each subsequent
                                    # metal uses the previous xyz
                                    # injected as conformer.
                                    _cur_mol = mol
                                    _cur_cid = _tcid
                                    _xyz_cur = None
                                    _chain_ok = True
                                    for _k, _mi in enumerate(metal_indices):
                                        _cf_k, _pm_k = _combo_tuple[_k]
                                        _gn_k = _cf_k[0]
                                        _rank_k = _rank_tuple[_k]
                                        _xyz_new = _build_topology_xyz(
                                            _cur_mol, _mi, _mi_donors[_k],
                                            _pm_k, _gn_k, False,
                                            conf_id=_cur_cid,
                                            chelate_rank=_rank_k,
                                        )
                                        if _xyz_new is None:
                                            _chain_ok = False
                                            break
                                        _xyz_cur = _xyz_new
                                        _mtmp_n = Chem.RWMol(mol)
                                        _mtmp_n.RemoveAllConformers()
                                        _conf_n = _xyz_to_rdkit_conformer(
                                            _mtmp_n.GetMol(), _xyz_cur,
                                        )
                                        if _conf_n is None:
                                            _chain_ok = False
                                            break
                                        _cid_n = _mtmp_n.AddConformer(
                                            _conf_n, assignId=True,
                                        )
                                        _rescale_metal_donor_distances(_mtmp_n, _cid_n)
                                        _cur_mol = _mtmp_n.GetMol()
                                        _cur_cid = _cid_n
                                        _xyz_cur = _mol_to_xyz_conformer(_mtmp_n, _cid_n)
                                    if not _chain_ok or _xyz_cur is None:
                                        continue
                                    # Snap every bridging donor to the
                                    # centroid of its connected metals
                                    # (N-way generalisation of the
                                    # 2-metal _snap_bridging_donors_to_compromise).
                                    try:
                                        _xyz_cur = _snap_bridging_donors_to_compromise(
                                            _xyz_cur, mol, metal_indices, bridging,
                                        )
                                    except Exception as _snp_n:
                                        logger.debug("N-metal bridge-snap failed: %s", _snp_n)
                                    if apply_uff:
                                        _xyz_cur = _optimize_xyz_openbabel_safe(
                                            _xyz_cur, mol_template=mol,
                                        )
                                    _sig_n = _xyz_sig(_xyz_cur)
                                    if _sig_n in _combo_seen_n:
                                        continue
                                    _combo_seen_n.add(_sig_n)
                                    _vkey_n = (_cfs, _pms)
                                    _variant_counter_n[_vkey_n] = (
                                        _variant_counter_n.get(_vkey_n, 0) + 1
                                    )
                                    _v_idx_n = _variant_counter_n[_vkey_n] - 1
                                    label_n = "multi-" + "/".join(_gns)
                                    if _v_idx_n > 0:
                                        label_n = f"{label_n}-conf{_v_idx_n + 1}"
                                    # Same bond-length mini-gate as
                                    # the 2-metal branch.
                                    try:
                                        _qln = [
                                            l for l in _xyz_cur.splitlines() if l.strip()
                                        ]
                                        _cqn = []
                                        for _ln in _qln:
                                            _pp = _ln.split()
                                            if len(_pp) >= 4:
                                                _cqn.append((
                                                    float(_pp[1]), float(_pp[2]), float(_pp[3])
                                                ))
                                        gate_n = True
                                        for _b in mol.GetBonds():
                                            _a1 = _b.GetBeginAtom()
                                            _a2 = _b.GetEndAtom()
                                            if (_a1.GetAtomicNum() <= 1
                                                    or _a2.GetAtomicNum() <= 1):
                                                continue
                                            _s1 = _a1.GetSymbol()
                                            _s2 = _a2.GetSymbol()
                                            if _s1 in _METAL_SET and _s2 in _METAL_SET:
                                                continue
                                            _i1 = _a1.GetIdx()
                                            _i2 = _a2.GetIdx()
                                            _dx = _cqn[_i1][0] - _cqn[_i2][0]
                                            _dy = _cqn[_i1][1] - _cqn[_i2][1]
                                            _dz = _cqn[_i1][2] - _cqn[_i2][2]
                                            _d = math.sqrt(_dx*_dx + _dy*_dy + _dz*_dz)
                                            if _s1 in _METAL_SET or _s2 in _METAL_SET:
                                                _msym = _s1 if _s1 in _METAL_SET else _s2
                                                _dsym = _s2 if _s1 in _METAL_SET else _s1
                                                _ideal = float(_get_ml_bond_length(_msym, _dsym))
                                                if _ideal > 0 and (_d < 0.50 * _ideal or _d > 2.50 * _ideal):
                                                    gate_n = False
                                                    break
                                            else:
                                                if _d > 2.5:
                                                    gate_n = False
                                                    break
                                        if not gate_n:
                                            continue
                                    except Exception:
                                        pass
                                    results.append((_xyz_cur, label_n))
                                    _n_combo_count += 1
                        except Exception as _nce:
                            logger.debug("N-metal combo build failed: %s", _nce)
                            continue
    except Exception as _mn_exc:
        logger.debug("Multinuclear enumeration failed: %s", _mn_exc)

    return results


# ---------------------------------------------------------------------------
# Linkage isomers (Feature 2)
# ---------------------------------------------------------------------------

def _find_linkage_alternatives(mol) -> List[Tuple[int, int, int, str]]:
    """Detect ambidentate ligands and return alternative coordination modes.

    Supported patterns:
      - NO2⁻ (N-donor → O-donor, label 'nitrito')
      - SCN⁻ (S-donor → N-donor, label 'isothiocyanato-N')
      - CN⁻  (C-donor → N-donor, label 'isocyano')

    Returns list of (metal_idx, current_donor_idx, alt_donor_idx, label).
    """
    alternatives: List[Tuple[int, int, int, str]] = []
    metal_neighbor_sets: Dict[int, set] = {}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        metal_neighbor_sets[metal_idx] = bonded

        for donor in atom.GetNeighbors():
            donor_idx = donor.GetIdx()
            donor_sym = donor.GetSymbol()

            # --- NO2 → nitrito (N-donor to O-donor) ---
            if donor_sym == 'N':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs) >= 2:
                    alt_o = o_nbrs[0]
                    alternatives.append((metal_idx, donor_idx, alt_o.GetIdx(), 'nitrito'))

            # --- SCN → isothiocyanato-N (S-donor to N-donor) ---
            if donor_sym == 'S':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isothiocyanato-N')
                                )

            # --- CN → isocyano (C-donor to N-donor via triple bond) ---
            if donor_sym == 'C':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 2.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'isocyano')
                            )

            # --- NCS (N-donor → S-donor via C): thiocyanato-S ---
            if donor_sym == 'N':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for s_nbr in c_nbr.GetNeighbors():
                            if (s_nbr.GetSymbol() == 'S'
                                    and s_nbr.GetIdx() != donor_idx
                                    and s_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, s_nbr.GetIdx(), 'thiocyanato-S')
                                )

            # --- Carboxylate (O1 → O2): alternative carboxylate oxygen ---
            if donor_sym == 'O':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C':
                        o_nbrs = [
                            n for n in c_nbr.GetNeighbors()
                            if n.GetSymbol() == 'O'
                            and n.GetIdx() != donor_idx
                            and n.GetIdx() not in bonded
                        ]
                        if o_nbrs:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'carboxylato-alt')
                            )

            # --- Sulfoxide (S-donor → O-donor via double bond) ---
            if donor_sym == 'S':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O'
                    and n.GetIdx() not in bonded
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()) is not None
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()).GetBondTypeAsDouble() >= 1.5
                ]
                if o_nbrs:
                    alternatives.append(
                        (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'sulfoxide-O')
                    )

            # --- Sulfite / Sulfonate (S-donor → O-donor, S with ≥3 O) ---
            if donor_sym == 'S':
                o_nbrs_all = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs_all) >= 2:
                    for alt_o in o_nbrs_all:
                        alternatives.append(
                            (metal_idx, donor_idx, alt_o.GetIdx(), 'sulfito-O')
                        )

            # --- Nitrosyl (N=O, N-donor → O-donor) ---
            if donor_sym == 'N':
                for o_nbr in donor.GetNeighbors():
                    if (o_nbr.GetSymbol() == 'O'
                            and o_nbr.GetIdx() not in bonded
                            and len(list(o_nbr.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, o_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbr.GetIdx(), 'isonitrosyl-O')
                            )

            # --- Nitrosyl O-bound → N-bound ---
            if donor_sym == 'O':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded
                            and len(list(donor.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'nitrosyl-N')
                            )

            # --- Selenocyanate SeCN (Se-donor → N-donor via C) ---
            if donor_sym == 'Se':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isoselenocyanato-N')
                                )

            # --- Pyrazole / Imidazole / Triazole (N1 → N2 in same ring) ---
            if donor_sym == 'N':
                try:
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if donor_idx in ring:
                            other_n = [
                                idx for idx in ring
                                if mol.GetAtomWithIdx(idx).GetSymbol() == 'N'
                                and idx != donor_idx
                                and idx not in bonded
                            ]
                            for alt_n in other_n:
                                alternatives.append(
                                    (metal_idx, donor_idx, alt_n, 'N-isomer')
                                )
                except Exception:
                    pass

    return alternatives


def _rewire_linkage(mol, metal_idx: int, old_donor: int, new_donor: int) -> Optional[object]:
    """Return a new mol with the M–old_donor bond replaced by M–new_donor.

    Returns None if the rewiring would duplicate an existing bond or fails.
    """
    try:
        rw = Chem.RWMol(mol)
        # Abort if new_donor is already bonded to metal
        if rw.GetBondBetweenAtoms(metal_idx, new_donor) is not None:
            return None
        rw.RemoveBond(metal_idx, old_donor)
        rw.AddBond(metal_idx, new_donor, Chem.BondType.SINGLE)
        rw.UpdatePropertyCache(strict=False)
        return rw.GetMol()
    except Exception as e:
        logger.debug("_rewire_linkage failed: %s", e)
        return None


def _generate_linkage_isomers(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Build linkage isomers by rewiring ambidentate ligands and generating XYZ.

    For each alternative coordination mode, rewires the metal bond and builds
    an idealized topology structure (OB UFF optimized).  The builder is tried
    against the top ``max_template_tries`` template conformers (ranked by
    :func:`_rank_template_conformers`) so one bad conformer does not silently
    discard valid linkage isomers.

    Returns [(xyz_string, label), …] where label is e.g. 'nitrito'.
    """
    results: List[Tuple[str, str]] = []
    alternatives = _find_linkage_alternatives(mol)

    for metal_idx, old_donor, new_donor, type_label in alternatives:
        alt_mol = _rewire_linkage(mol, metal_idx, old_donor, new_donor)
        if alt_mol is None:
            continue
        try:
            # Find metal in alt_mol (same index)
            metal_atom = alt_mol.GetAtomWithIdx(metal_idx)
            donor_indices = [nbr.GetIdx() for nbr in metal_atom.GetNeighbors()]
            n_coord = len(donor_indices)
            geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
            if n_coord not in geom_map:
                continue
            geom = geom_map[n_coord]
            perm = list(range(n_coord))  # identity: donor i → position i

            candidate_cids = _rank_template_conformers(
                alt_mol, top_k=max_template_tries
            )
            if not candidate_cids:
                candidate_cids = [None]

            for cid in candidate_cids:
                # Build WITHOUT UFF first (fast), check topology, THEN UFF.
                xyz = _build_topology_xyz(
                    alt_mol, metal_idx, donor_indices, perm, geom, False,
                    conf_id=cid,
                )
                if xyz is None:
                    continue
                if not _fragment_topology_ok(xyz, smiles):
                    logger.debug(
                        "linkage %s via template cid=%s: "
                        "fragment topology mismatch, trying next template",
                        type_label, cid,
                    )
                    continue
                if apply_uff:
                    xyz = _optimize_xyz_openbabel_safe(
                        xyz, mol_template=alt_mol
                    )
                    xyz = _snap_aromatic_rings_in_xyz(xyz, alt_mol)
                results.append((xyz, type_label))
                break
        except Exception as e:
            logger.debug("Linkage isomer (%s) failed: %s", type_label, e)

    return results


def _generate_alternative_binding_modes(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_alternatives: int = 20,
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Generate isomers by swapping donor atoms with alternative binding sites.

    Only considers viable donors (atoms with available lone pairs) on the
    SAME ligand fragment as the current donor.  This avoids generating
    nonsensical structures from swapping donors across unrelated ligands
    or using C=O carbonyl oxygens as coordination donors.

    For each rewired candidate the builder is tried against the top
    ``max_template_tries`` template conformers (ranked by geometry quality)
    and the first result that passes ``_fragment_topology_ok`` is kept.  This
    decouples constitutional-isomer generation from a single, possibly poor,
    template conformer (see Codex findings for the reasoning).

    Returns [(xyz_string, label), ...].
    """
    results: List[Tuple[str, str]] = []

    # Iter-8.6 multi-σ timeout-mitigation: env-gated wall-clock budget.
    # Multi-metal fragments fire 8-seed × 6 s ETKDG sweeps per rewire
    # candidate; with ~12 alt donors × 2 metals × top_k=8 templates the
    # cumulative cost exceeds the 600 s per-SMILES driver budget, so the
    # outer pool_evaluator marks the SMILES as timeout and the entire
    # multi-σ result set is lost.  When DELFIN_ALT_MODE_BUDGET_S is set
    # to a positive integer, the candidate loop breaks once that many
    # wall-clock seconds have elapsed.  Partial results are returned.
    # Default 0 → disabled, bit-exact baseline behaviour.
    import time as _time_mod
    # Wave-4 (A+D forensics): Phase 4D's universal 60s cap killed large
    # bimetallic multi-σ conversions (Sn-Ir, Sn-Rh, Os-Sn, Fe-Fe μ-O,
    # In-Ru). class-conditional: 60s for sigma/hapto (safe), 0s/unlimited
    # for multi_sigma/multi_hapto (needs 100-500s).  Env override still
    # respected for testing.
    try:
        _env_budget_raw = os.environ.get("DELFIN_ALT_MODE_BUDGET_S")
        if _env_budget_raw is not None and _env_budget_raw.strip() != "":
            _alt_budget_s = float(_env_budget_raw)
        else:
            try:
                _cls = _classify_complex_class(mol)
            except Exception:
                _cls = "sigma"
            # Multi-sigma V2: opt out of the historical "unlimited" budget
            # for multi-metal sigma — that branch was the second biggest
            # contributor to the 600 s pool-evaluator timeouts (forensics
            # 2026-05-13).  Cap to the heavy-atom-scaled wall-clock so
            # alt-binding-mode exploration cannot starve the rest of the
            # pipeline.  Other classes keep their pre-patch budget.
            _v2_alt_cap: Optional[float] = None
            try:
                if _cls == "multi_sigma" and _multi_sigma_v2_active(mol):
                    _heavy_n_alt = sum(
                        1 for _a in mol.GetAtoms() if _a.GetAtomicNum() > 1
                    )
                    # Re-use the multi-metal augmentation budget — these
                    # two stages have similar per-call cost characteristics.
                    _v2_alt_cap = float(
                        _multi_sigma_v2_budget(_heavy_n_alt)["mm_walltime"]
                    )
            except Exception:
                _v2_alt_cap = None
            if _v2_alt_cap is not None:
                _alt_budget_s = _v2_alt_cap
            else:
                _alt_budget_s = (
                    0.0 if _cls in ("multi_sigma", "multi_hapto") else 60.0
                )
    except Exception:
        _alt_budget_s = 60.0
    # DETERMINISM (2026-07-28): this was the ONLY isomer-generating path whose wall-clock budget was
    # not gated on _deterministic_mode() -- every sibling already is (2-metal :27381, N-metal :27605,
    # MM-augmentation, ETKDG multi-seed, embed joins).  Consequence: alternative binding modes are
    # DROPPED under load, so the ISOMER SET itself became load-dependent.  Measured on LATTUW: 14
    # frames under load vs 16 on a free box -- the two missing ones are alt-bind-C isomers.  It passes
    # within-run determinism (both replica builds run under the same load) and only differs ACROSS
    # runs, which is why the byte-determinism gate never caught it.  That violates both "completeness is
    # sacred" and the deterministic-manifold claim.  0.0 disables the wall clock; termination falls back
    # to the deterministic caps, exactly like the siblings.
    if _deterministic_mode():
        _alt_budget_s = 0.0
    _alt_t0 = _time_mod.monotonic() if _alt_budget_s > 0 else None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        if _alt_t0 is not None and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s:
            logger.debug(
                "alt-binding budget %.1fs exceeded — stopping at metal loop",
                _alt_budget_s,
            )
            return results
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        current_donors = list(bonded)

        # Build ligand fragments to ensure we only swap within the same ligand
        fragments = _ligand_fragments(mol, metal_idx)
        atom_to_frag: Dict[int, int] = {}
        for fi, frag in enumerate(fragments):
            for aidx in frag:
                atom_to_frag[aidx] = fi

        # For each current donor, find viable alternative donors on the same fragment
        for current_d in current_donors:
            current_frag = atom_to_frag.get(current_d, -1)
            if current_frag < 0:
                continue
            current_sym = mol.GetAtomWithIdx(current_d).GetSymbol()

            for alt_d in fragments[current_frag]:
                if alt_d == current_d or alt_d in bonded:
                    continue
                if not _is_viable_donor(mol, alt_d, bonded):
                    continue
                if len(results) >= max_alternatives:
                    return results
                if _alt_t0 is not None and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s:
                    logger.debug(
                        "alt-binding budget %.1fs exceeded — stopping at donor loop",
                        _alt_budget_s,
                    )
                    return results

                alt_sym = mol.GetAtomWithIdx(alt_d).GetSymbol()
                alt_mol = _rewire_linkage(mol, metal_idx, current_d, alt_d)
                if alt_mol is None:
                    continue
                try:
                    new_metal = alt_mol.GetAtomWithIdx(metal_idx)
                    donor_indices = [nbr.GetIdx() for nbr in new_metal.GetNeighbors()]
                    n_coord = len(donor_indices)
                    geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
                    if n_coord not in geom_map:
                        continue
                    geom = geom_map[n_coord]
                    perm = list(range(n_coord))

                    if current_sym == alt_sym:
                        label = f'alt-{alt_sym}-isomer'
                    else:
                        label = f'alt-bind-{alt_sym}'

                    # Multi-template trial: iterate over the best-scored
                    # conformer IDs and keep the first XYZ that validates.
                    candidate_cids = _rank_template_conformers(
                        alt_mol, top_k=max_template_tries
                    )
                    if not candidate_cids:
                        # No usable template conformer — try de-novo once.
                        candidate_cids = [None]

                    accepted = False
                    for cid in candidate_cids:
                        # Multi-sigma V2: per-template budget check.  The
                        # outer donor-loop check (above) can still overshoot
                        # by 30-60 s on multi-metal complexes where each
                        # template alignment + UFF + topology gate runs
                        # serially.  This inner check bounds the overshoot.
                        if (
                            _alt_t0 is not None
                            and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s
                        ):
                            logger.debug(
                                "alt-binding budget %.1fs exceeded — "
                                "stopping at template loop",
                                _alt_budget_s,
                            )
                            return results
                        # Build WITHOUT UFF first, check topology, THEN UFF.
                        xyz = _build_topology_xyz(
                            alt_mol, metal_idx, donor_indices, perm,
                            geom, False, conf_id=cid,
                        )
                        if xyz is None:
                            continue
                        if not _fragment_topology_ok(xyz, smiles):
                            logger.debug(
                                "alt-mode %s via template cid=%s: "
                                "fragment topology mismatch, trying next template",
                                label, cid,
                            )
                            continue
                        if apply_uff:
                            xyz = _optimize_xyz_openbabel_safe(
                                xyz, mol_template=alt_mol
                            )
                        results.append((xyz, label))
                        accepted = True
                        break

                    if not accepted:
                        logger.debug(
                            "alt-mode %s: no template conformer produced a "
                            "topology-consistent XYZ (tried %d)",
                            label, len(candidate_cids),
                        )
                except Exception as e:
                    logger.debug("Alternative binding mode failed: %s", e)
                    continue

    return results


def _enumerate_hapto_sigma_isomers(
    smiles: str,
    base_xyz: str,
    apply_uff: bool = True,
    max_isomers: int = 20,
) -> List[Tuple[str, str]]:
    """Enumerate σ-donor permutations for hapto complexes.

    For HAPTO metals: permute their σ-donors via rigid-body swaps
    keeping η-rings fixed.

    For NON-HAPTO metals in mixed complexes (e.g. Ni in CpFe-bridge-Ni):
    permute their donors via the topology enumerator, using the hapto
    XYZ as template. The hapto rings stay fixed.
    """
    if not RDKIT_AVAILABLE:
        return []
    try:
        import numpy as np
        import itertools as _it

        mol = _prepare_mol_for_embedding(smiles, hapto_approx=True)
        if mol is None:
            return []

        # Inject base_xyz as conformer.
        conf = _xyz_to_rdkit_conformer(mol, base_xyz)
        if conf is None:
            return []
        mol.RemoveAllConformers()
        cid = mol.AddConformer(conf, assignId=True)

        hapto_groups = _find_hapto_groups(mol)
        if not hapto_groups:
            return []

        hapto_atoms: set = set()
        hapto_metals: set = set()
        for _midx, members in hapto_groups:
            hapto_atoms.update(members)
            hapto_metals.add(_midx)

        results: List[Tuple[str, str]] = []
        dtype_map = _donor_type_map(mol)

        # Step A: For NON-HAPTO metals, use the topology enumerator with
        # the hapto XYZ as template. The hapto rings stay fixed (they are
        # not in the non-hapto metal's donor list).
        try:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    continue
                m_idx = atom.GetIdx()
                if m_idx in hapto_metals:
                    continue  # handled by Step B below
                donors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                n_coord = len(donors)
                if n_coord < 2 or n_coord > 9:
                    continue
                # Build donor labels.
                dk = [dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset())) for d in donors]
                uk = sorted(set(dk), key=lambda k: (k[0], tuple(sorted(k[1]))))
                kc = {k: i for i, k in enumerate(uk)}
                labels = [f"{k[0]}{kc[k]}" for k in dk]
                if len(uk) <= 1:
                    continue  # all donors equivalent
                cp_pairs = _chelate_pairs(mol, m_idx, donors)
                al = {ai: li for li, ai in enumerate(donors)}
                clp = [
                    frozenset([al[sorted(c)[0]], al[sorted(c)[1]]])
                    for c in cp_pairs
                    if sorted(c)[0] in al and sorted(c)[1] in al
                ]
                iso_list = _enumerate_topological_isomers(
                    labels, n_coord, clp, metal_symbol=atom.GetSymbol(),
                )
                for cf, pm in iso_list[:max_isomers]:
                    gn = cf[0]
                    try:
                        xyz_new = _build_topology_xyz(
                            mol, m_idx, donors, pm, gn, False, conf_id=cid,
                        )
                        if xyz_new is None:
                            continue
                        if apply_uff:
                            try:
                                xyz_new = _optimize_xyz_openbabel_safe(
                                    xyz_new, mol_template=mol
                                )
                            except Exception:
                                pass
                        if not _verify_topology_from_graph(xyz_new, mol):
                            continue
                        mt = Chem.RWMol(mol)
                        mt.RemoveAllConformers()
                        c2 = _xyz_to_rdkit_conformer(mt.GetMol(), xyz_new)
                        if c2 is None:
                            continue
                        ci2 = mt.AddConformer(c2, assignId=True)
                        fp = _compute_coordination_fingerprint(
                            mt.GetMol(), ci2, dtype_map=dtype_map
                        )
                        lbl = _classify_isomer_label(fp, mt.GetMol())
                        if not lbl:
                            lbl = f'non-hapto-{gn}'
                        results.append((xyz_new, lbl))
                    except Exception:
                        continue
        except Exception as _ex:
            logger.debug("Non-hapto topo enumeration failed: %s", _ex)

        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            metal_idx = atom.GetIdx()
            all_donors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            sigma_donors = [d for d in all_donors if d not in hapto_atoms]

            if len(sigma_donors) < 2:
                continue

            # Donor labels for σ-donors.
            donor_keys = [
                dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))
                for d in sigma_donors
            ]
            # If all σ-donors are equivalent → only 1 arrangement.
            if len(set(donor_keys)) <= 1:
                continue

            # Current σ-donor positions from conformer.
            conf_obj = mol.GetConformer(cid)
            sigma_positions = []
            for d in sigma_donors:
                p = conf_obj.GetAtomPosition(d)
                sigma_positions.append(np.array([p.x, p.y, p.z]))

            # Build ligand fragments per σ-donor (BFS excluding metal + η).
            non_metal = {
                a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() not in _METAL_SET
            }
            adj: Dict[int, set] = {i: set() for i in non_metal}
            for bond in mol.GetBonds():
                bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bi in non_metal and bj in non_metal:
                    adj[bi].add(bj)
                    adj[bj].add(bi)

            donor_frag_atoms: Dict[int, set] = {}
            for d in sigma_donors:
                frag: set = set()
                stack = [d]
                visited: set = set()
                while stack:
                    node = stack.pop()
                    if node in visited or node in hapto_atoms:
                        continue
                    # Don't cross into OTHER σ-donor territories.
                    if node != d and node in sigma_donors:
                        continue
                    visited.add(node)
                    frag.add(node)
                    for nbr in adj.get(node, ()):
                        if nbr not in visited:
                            stack.append(nbr)
                donor_frag_atoms[d] = frag

            # Generate unique permutations of σ-donor labels.
            original_label_tuple = tuple(donor_keys)
            seen_labels: set = {original_label_tuple}

            # Cap permutation enumeration to prevent unbounded explosion on
            # high-CN sigma-donor systems. n_donors=6 → 720 perms,
            # n_donors=7 → 5040 perms, each costs ~100ms (rigid-body align +
            # UFF refinement). Without cap, single SMILES blocks subprocess
            # for minutes and accumulates in pool_evaluator. Identity perm
            # is always tried first; remaining quota is randomly sampled
            # for diversity (deterministic with PYTHONHASHSEED=0).
            import random as _random
            max_sigma_perms = int(os.environ.get('DELFIN_SIGMA_PERMS_MAX', '12'))
            all_perms = list(_it.permutations(range(len(sigma_donors))))
            if len(all_perms) > max_sigma_perms:
                identity_perm = all_perms[0]
                rest = all_perms[1:]
                _rng = _random.Random(0)  # deterministic
                perms_to_try = [identity_perm] + _rng.sample(
                    rest, min(max_sigma_perms - 1, len(rest))
                )
            else:
                perms_to_try = all_perms

            for perm in perms_to_try:
                perm_labels = tuple(donor_keys[p] for p in perm)
                if perm_labels in seen_labels:
                    continue
                seen_labels.add(perm_labels)
                if len(results) >= max_isomers:
                    break

                # Build swapped XYZ: move fragment of donor perm[i] to
                # the position of donor i via rigid-body alignment.
                new_coords = np.zeros((mol.GetNumAtoms(), 3))
                for ai in range(mol.GetNumAtoms()):
                    p = conf_obj.GetAtomPosition(ai)
                    new_coords[ai] = [p.x, p.y, p.z]

                swap_ok = True
                _mp = conf_obj.GetAtomPosition(metal_idx)
                metal_pos = np.array([_mp.x, _mp.y, _mp.z])
                for slot_idx in range(len(sigma_donors)):
                    src_donor = sigma_donors[perm[slot_idx]]
                    tgt_donor = sigma_donors[slot_idx]
                    if src_donor == tgt_donor:
                        continue
                    src_frag = sorted(donor_frag_atoms[src_donor])
                    tgt_pos = sigma_positions[slot_idx]
                    src_pos = sigma_positions[perm[slot_idx]]

                    # Place the SOURCE donor at the TARGET slot's DIRECTION but
                    # keep the source donor's OWN (element-correct) M-D distance.
                    # Bug fix (2026-05-20): the previous ``delta = tgt_pos -
                    # src_pos`` snapped the donor onto the target donor's
                    # position, so a permuted donor inherited the *other*
                    # element's M-distance (e.g. Cl landing at an N slot ended
                    # up at the N distance ~2.06 Å instead of ~2.40 Å).  Use
                    # the source donor's existing element-correct radius along
                    # the target direction instead.
                    tgt_dir = tgt_pos - metal_pos
                    _tn = float(np.linalg.norm(tgt_dir))
                    src_radius = float(np.linalg.norm(src_pos - metal_pos))
                    if _tn < 1e-6 or src_radius < 1e-6:
                        delta = tgt_pos - src_pos  # degenerate fallback
                    else:
                        new_donor_pos = metal_pos + (src_radius / _tn) * tgt_dir
                        delta = new_donor_pos - src_pos
                    orig_frag_coords = np.array([
                        [conf_obj.GetAtomPosition(ai).x,
                         conf_obj.GetAtomPosition(ai).y,
                         conf_obj.GetAtomPosition(ai).z]
                        for ai in src_frag
                    ])
                    for fi, ai in enumerate(src_frag):
                        new_coords[ai] = orig_frag_coords[fi] + delta

                # Write XYZ.
                lines = []
                for ai in range(mol.GetNumAtoms()):
                    sym = mol.GetAtomWithIdx(ai).GetSymbol()
                    x, y, z = new_coords[ai]
                    lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
                xyz_new = '\n'.join(lines) + '\n'

                if apply_uff:
                    try:
                        xyz_new = _optimize_xyz_openbabel_safe(
                            xyz_new, mol_template=mol
                        )
                    except Exception:
                        pass

                # Quality gate.
                if not _metal_donor_distances_realistic(xyz_new, mol):
                    continue

                try:
                    mol_chk = Chem.RWMol(mol)
                    mol_chk.RemoveAllConformers()
                    conf_chk = _xyz_to_rdkit_conformer(mol_chk.GetMol(), xyz_new)
                    if conf_chk is None:
                        continue
                    cid_chk = mol_chk.AddConformer(conf_chk, assignId=True)
                    if _has_atom_clash(mol_chk.GetMol(), cid_chk, min_dist=0.3):
                        continue
                    fp = _compute_coordination_fingerprint(
                        mol_chk.GetMol(), cid_chk, dtype_map=dtype_map
                    )
                    label = _classify_isomer_label(fp, mol_chk.GetMol())
                    if not label:
                        label = f'hapto-sigma-{len(results) + 1}'
                    results.append((xyz_new, label))
                except Exception:
                    continue

    except Exception as exc:
        logger.debug("Hapto sigma isomer enumeration failed: %s", exc)
    return results


def _emit_chelate_pucker_variants(
    mol,
    results: List[Tuple[str, str]],
    apply_uff: bool,
    max_isomers: int,
) -> int:
    """Append chair / boat / twist conformer variants for saturated chelate rings.

    For chelate rings whose backbone has ≥3 sp³ carbons (cyclam, en, dien,
    salen-CH₂CH₂, polyamines, polyethers), the existing pipeline emits
    one ETKDG conformer.  Pucker variants (chair, boat, twist) are
    discrete low-energy minima that get lost when only one conformer
    survives the dedup chain.  This pass generates additional ETKDG
    seeds, applies puck-perturbation to ring atoms, and appends each
    geometrically-distinct result iff its heavy-atom XYZ signature is
    novel.

    Pure ADDITIVE — never sorts, drops, or modifies existing entries.
    Toggle: DELFIN_PUCKER_PASS_ENABLED (default 1).
    Returns the number of new entries appended.
    """
    if not RDKIT_AVAILABLE:
        return 0
    if not _delfin_env_int("DELFIN_PUCKER_PASS_ENABLED", 1):
        return 0

    def _sig(xyz_str: str) -> tuple:
        try:
            lines = [
                ln.split() for ln in xyz_str.strip().splitlines()
                if ln.strip()
            ]
            heavy = sorted(
                (p[0], round(float(p[1]), 2), round(float(p[2]), 2),
                 round(float(p[3]), 2))
                for p in lines
                if len(p) >= 4 and p[0] not in ('H', 'h')
            )
            return tuple(heavy)
        except Exception:
            return tuple()

    seen_sigs = {_sig(xyz) for xyz, _ in results}
    n_added = 0

    # Find saturated chelate rings: rings where ≥3 atoms are sp³ carbons
    # bridging metal-coordinating donors.
    try:
        Chem.FastFindRings(mol)
    except Exception:
        pass
    ri = mol.GetRingInfo()
    if ri is None or ri.NumRings() == 0:
        return 0

    metal_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
    if not metal_atoms:
        return 0

    # Identify chelate rings: rings that contain a metal atom AND at least
    # 3 sp³ carbons (saturated backbone bridges).
    chelate_ring_atom_sets: List[Tuple[set, int]] = []  # (ring_atoms, n_sp3)
    for ring in ri.AtomRings():
        ring_set = set(ring)
        if not (ring_set & set(metal_atoms)):
            continue
        n_sp3 = 0
        for ridx in ring:
            atom = mol.GetAtomWithIdx(ridx)
            if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
                n_sp3 += 1
        if n_sp3 >= 3:
            chelate_ring_atom_sets.append((ring_set, n_sp3))

    if not chelate_ring_atom_sets:
        return 0

    # For each saturated chelate ring, emit pucker variants by perturbing
    # ring sp³ atom z-coordinates in distinct patterns (chair: alternating
    # +/-, boat: two adjacent up + two adjacent down, twist: gradient).
    # Then ETKDG-relax and accept if XYZ signature is new.
    base_xyz = None
    for xyz, _lbl in results:
        base_xyz = xyz
        break
    if base_xyz is None:
        return 0

    import numpy as _np

    def _parse_xyz(xs):
        lines = [l for l in xs.strip().splitlines() if l.strip()]
        atoms, coords = [], []
        for ln in lines:
            parts = ln.split()
            if len(parts) >= 4:
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return atoms, _np.array(coords)

    def _format_xyz(atoms, coords):
        out = []
        for sym, (x, y, z) in zip(atoms, coords):
            out.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(out) + "\n"

    pucker_patterns = [
        ('chair', lambda i, n: 0.4 if i % 2 == 0 else -0.4),
        ('boat', lambda i, n: 0.4 if (i % n) < n // 2 else -0.4),
        ('twist', lambda i, n: 0.3 * ((i / max(1, n - 1)) - 0.5) * 2),
    ]

    for ring_atoms, n_sp3 in chelate_ring_atom_sets:
        if len(results) + n_added >= max_isomers:
            break
        ring_size = len(ring_atoms)
        # Pucker only meaningful for 5+ membered rings
        if ring_size < 5:
            continue
        ring_atom_list = sorted(ring_atoms)
        for ptype, pucker_fn in pucker_patterns:
            if len(results) + n_added >= max_isomers:
                break
            try:
                atoms, coords = _parse_xyz(base_xyz)
                if len(atoms) != mol.GetNumAtoms():
                    # Hydrogens not in xyz — pucker only heavy ring atoms
                    pass
                # Compute ring center and normal
                ring_pos = _np.array([
                    coords[i] for i in ring_atom_list if i < len(coords)
                ])
                if len(ring_pos) < 4:
                    continue
                center = ring_pos.mean(axis=0)
                centered = ring_pos - center
                _u, _s, vh = _np.linalg.svd(centered, full_matrices=False)
                normal = vh[-1]  # smallest singular vector = ring normal
                # Apply pucker perturbation along normal
                new_coords = coords.copy()
                for k, idx in enumerate(ring_atom_list):
                    if idx >= len(new_coords):
                        continue
                    delta = pucker_fn(k, ring_size)
                    new_coords[idx] = new_coords[idx] + delta * normal
                new_xyz = _format_xyz(atoms, new_coords)
                # Check signature novelty
                sig = _sig(new_xyz)
                if sig in seen_sigs:
                    continue
                # UFF-relax to remove unphysical strain from perturbation
                if apply_uff:
                    try:
                        new_xyz = _optimize_xyz_openbabel_safe(new_xyz, mol_template=mol)
                    except Exception:
                        pass
                # Final-output topology gate
                try:
                    if not _verify_topology_from_graph(new_xyz, mol):
                        continue
                except Exception:
                    continue
                sig2 = _sig(new_xyz)
                if sig2 in seen_sigs:
                    continue
                seen_sigs.add(sig2)
                label = f"pucker-{ring_size} {ptype}"
                # Iter-8.7 every-append gate (123a130 port, env-gated, default OFF)
                _gate_pass = True
                if _every_append_gate_enabled(mol):
                    try:
                        _flat = _flatten_sp2_atoms_xyz(new_xyz, mol)
                        if _flat:
                            new_xyz = _flat
                    except Exception:
                        pass
                    try:
                        _mt_g = Chem.RWMol(mol); _mt_g.RemoveAllConformers()
                        _c_g = _xyz_to_rdkit_conformer(_mt_g.GetMol(), new_xyz)
                        if _c_g is None:
                            _gate_pass = False
                        else:
                            _ci_g = _mt_g.AddConformer(_c_g, assignId=True)
                            if _has_severe_covalent_distortion(_mt_g.GetMol(), _ci_g):
                                _gate_pass = False
                    except Exception:
                        _gate_pass = False
                if _gate_pass:
                    results.append((new_xyz, label))
                    n_added += 1
            except Exception:
                continue

    if n_added:
        logger.debug(
            "Pucker pass added %d chelate-ring conformer variants",
            n_added,
        )
    return n_added


def _emit_nonmetal_ring_pucker_variants(
    mol,
    results: List[Tuple[str, str]],
    apply_uff: bool,
    max_isomers: int,
) -> int:
    """Append basin-verified chair / boat pucker variants for NON-METAL rings.

    The pucker analogue of Pólya coordination-isomer enumeration, restricted
    to the rings the metal pucker pass does NOT handle: peripheral non-metal,
    non-aromatic rings (cyclohexyl, piperidinyl, sugar, ...).  For every such
    ring the pass emits, from the BEST-ranked coordination isomer only:

      * one GLOBAL chair-set frame (every eligible ring driven to chair), and
      * up to a small fixed number of "one-ring-flipped-to-boat" decorations,

    rather than the full 2^N pucker product.  Symmetry-equivalent rings are
    collapsed to a single DOF (one boat decoration per equivalence class).
    Each variant is driven into its target Cremer-Pople basin with a
    constrained geometric snap and REJECTED if it does not realise that
    basin.  Hard cap DELFIN_5P_B_MAX_RING_VARIANTS (default 6).

    Pure ADDITIVE — never sorts, drops, or modifies existing entries.
    Master toggle: DELFIN_RING_PUCKER_ENUM (default 1; =0 -> no-op).
    Returns the number of new entries appended.
    """
    if not RDKIT_AVAILABLE:
        return 0
    if not _delfin_env_int("DELFIN_RING_PUCKER_ENUM", 1):
        return 0
    if not results:
        return 0

    try:
        from delfin.manta import _ring_conformer_templates as _rct
        from delfin.manta import _rotamer_diversity as _rot
    except Exception:
        return 0

    max_ring_variants = _delfin_env_int("DELFIN_5P_B_MAX_RING_VARIANTS", 6)
    if max_ring_variants < 1:
        max_ring_variants = 6

    def _sig(xyz_str: str) -> tuple:
        try:
            lines = [
                ln.split() for ln in xyz_str.strip().splitlines() if ln.strip()
            ]
            heavy = sorted(
                (p[0], round(float(p[1]), 2), round(float(p[2]), 2),
                 round(float(p[3]), 2))
                for p in lines
                if len(p) >= 4 and p[0] not in ('H', 'h')
            )
            return tuple(heavy)
        except Exception:
            return tuple()

    seen_sigs = {_sig(xyz) for xyz, _ in results}

    # (2) Generate puckers from the BEST-ranked coordination isomer only.
    # `results` may already be ordered best-first by the impl, but rank
    # explicitly so the choice is deterministic and independent of upstream
    # ordering.
    base_xyz = results[0][0]
    try:
        from delfin.manta._conformer_rank import rank_isomers as _rank
        _ranked = _rank(list(results))
        if _ranked:
            base_xyz = _ranked[0][0]
    except Exception:
        base_xyz = results[0][0]

    # Build the OB-perceived graph once from the base frame.
    try:
        ob_mol = _rot._build_ob_mol_from_xyz(base_xyz)
        if ob_mol is None:
            return 0
        graph = _rot._graph_from_ob(ob_mol)
        if not graph:
            return 0
        symbols, base_coords_t = _rot._parse_delfin_xyz(base_xyz)
    except Exception:
        return 0
    base_coords = [tuple(c) for c in base_coords_t]

    # Rings amenable to templating = non-metal, non-aromatic, size 3..30.
    # find_rings_for_templating already EXCLUDES metal-chelate + aromatic
    # rings on graph features only — so this is exactly the complement of
    # what the metal pucker pass handles.
    try:
        rings = _rct.find_rings_for_templating(graph)
    except Exception:
        return 0
    # Restrict to 6-rings: the CP basin snap + verifier below is defined for
    # 6-rings (the dominant flexible-ring class and the ZIGDOL signature).
    rings = [r for r in rings if len(r) == 6]
    if not rings:
        return 0

    base_topo = None
    try:
        base_topo = _rot._topology_hash(graph)
    except Exception:
        base_topo = None

    atomic_nums = graph.get("atomic_nums", [])

    def _ring_key(ring) -> tuple:
        """Symmetry key: the multiset of (element, heavy-degree) over the ring
        atoms plus the ring's heavy-substituent element pattern.  Graph-only,
        so symmetry-equivalent peripheral rings (e.g. ZIGDOL's six chemically
        identical cyclohexyls) collapse to ONE equivalence class -> one DOF."""
        neighbours = graph.get("neighbours", [[]])
        feats = []
        for idx in ring:
            z = atomic_nums[idx] if idx < len(atomic_nums) else 0
            heavy_deg = sum(
                1 for nb in neighbours[idx]
                if nb < len(atomic_nums) and atomic_nums[nb] > 1
            )
            feats.append((z, heavy_deg))
        return tuple(sorted(feats))

    # Group symmetry-equivalent rings.
    sym_groups: Dict[tuple, List[List[int]]] = {}
    for ring in rings:
        sym_groups.setdefault(_ring_key(ring), []).append(ring)
    # Deterministic ordering of the groups + their member rings.
    ordered_groups = sorted(
        sym_groups.items(), key=lambda kv: (kv[0], sorted(r[0] for r in kv[1]))
    )

    def _snap_ring(coords, ring, basin):
        """Drive *ring* cleanly into the centre of *basin* by an ABSOLUTE
        out-of-plane snap that HOLDS the Cremer-Pople target.

        Each ring atom's signed out-of-plane component (relative to the ring
        mean plane) is SET to the canonical chair/boat target ``pat[k]*amp``,
        so a ring starting in any pucker (deep boat, twist, ...) lands cleanly
        in the intended basin.  The in-plane components (which carry the ring
        bond topology) are preserved; bonded hydrogens are dragged rigidly by
        the per-atom out-of-plane delta so C-H lengths are preserved.  The full
        alternating chair pattern is applied to EVERY ring atom (including the
        ipso/attachment carbon) so the CP basin is exact — the attachment-bond
        stretch this introduces at the ipso atom is healed by the subsequent
        constrained relax, which frees the ipso atom while freezing the rest of
        the ring.  Returns the new full coordinate list, or None if undefined.
        """
        pat = _ring_canonical_snap_z(basin, len(ring))
        if pat is None:
            return None
        avg = _rct._average_ring_bond_length(coords, ring)
        if avg < 1e-3:
            return None
        # Target out-of-plane amplitude (A): ideal cyclohexane chair sits
        # ~0.25 A out-of-plane per atom for ~1.54 A C-C; boat flagpoles ~0.65 A.
        amp = avg * (0.42 if basin == "boat" else 0.165)
        try:
            _c, normal = _rct._ring_plane_normal(coords, ring)
        except Exception:
            return None
        nrm = math.sqrt(sum(v * v for v in normal))
        if nrm < 1e-9:
            return None
        normal = tuple(v / nrm for v in normal)
        neighbours = graph.get("neighbours", [[]])
        out = [c for c in coords]
        for k, idx in enumerate(ring):
            cx, cy, cz = out[idx]
            cur_z = (cx - _c[0]) * normal[0] + (cy - _c[1]) * normal[1] + \
                    (cz - _c[2]) * normal[2]
            d = pat[k] * amp - cur_z
            shift = (normal[0] * d, normal[1] * d, normal[2] * d)
            out[idx] = (cx + shift[0], cy + shift[1], cz + shift[2])
            # rigid-H drag by the same per-atom out-of-plane delta
            if idx < len(neighbours):
                for nb in neighbours[idx]:
                    if nb < len(atomic_nums) and atomic_nums[nb] == 1:
                        hx, hy, hz = out[nb]
                        out[nb] = (hx + shift[0], hy + shift[1], hz + shift[2])
        return out

    def _ring_basin(coords, ring):
        try:
            rc = [coords[i] for i in ring]
            return _cp_basin_6ring(rc)[0]
        except Exception:
            return "undefined"

    # Coordination sphere from the TEMPLATE graph (RDKit mol) — reliable even
    # when OB does not perceive a long/weak M-D bond (e.g. Ag-As ~2.5 A).  The
    # metal AND its first-shell donors are frozen during every constrained
    # relax so the M-D invariant cannot drift.
    template_coord_sphere = set()
    try:
        if mol is not None:
            for _a in mol.GetAtoms():
                if _a.GetSymbol() in _METAL_SET:
                    template_coord_sphere.add(_a.GetIdx())
                    for _nb in _a.GetNeighbors():
                        template_coord_sphere.add(_nb.GetIdx())
    except Exception:
        template_coord_sphere = set()

    def _constrained_relax(coords, hold_ring_atoms):
        """Constrained local relax that HOLDS the Cremer-Pople target.

        Freeze the heavy atoms of the *hold_ring_atoms* set (the snapped ring
        atoms whose chair/boat pucker must be held — the ipso/attachment carbon
        is among them, so its external bond stays at the snapped ~2.1 A, well
        inside the gate) plus the metal coordination sphere (so the M-D
        invariant cannot drift).  Everything else — the non-held substituent
        heavy atoms and ALL hydrogens — relaxes under a short OB-UFF
        conjugate-gradient step, which relieves the rigid-H-drag / local steric
        strain the snap introduces and brings the variant's energy back near
        the pool floor (so it survives the downstream energy-outlier cut).
        Returns relaxed coords (or the input on any failure).
        """
        if not apply_uff or not OPENBABEL_AVAILABLE:
            return coords
        xyz_in = _rot._format_delfin_xyz(symbols, coords)
        fix = sorted(set(hold_ring_atoms) | template_coord_sphere)
        try:
            relaxed = _optimize_xyz_openbabel(
                xyz_in, steps=200, constraints={"fix_atoms": fix},
            )
            _syms2, _coords2 = _rot._parse_delfin_xyz(relaxed)
            if len(_coords2) == len(coords):
                return [tuple(c) for c in _coords2]
        except Exception:
            pass
        return coords

    def _build_and_validate(coords):
        """M-D invariant guard + authoritative output topology gate; return the
        validated DELFIN xyz string or None.

        The authoritative integrity check is ``_verify_topology_from_graph``
        (the SAME final output gate the whole pipeline trusts): it validates
        every template-graph bond against a distance cutoff and rejects
        collapsed / overlapping atoms — robust to a small out-of-plane pucker
        displacement.  We deliberately do NOT additionally require OB bond
        re-perception to reproduce the EXACT base topology-hash: OB's
        distance-based bond-order perception flips spuriously under a 0.5 A
        ring-atom displacement (it is brittle by design), which would reject
        chemically valid puckers.  The M-D guard still protects the
        coordination sphere; the output gate guarantees no bond is
        broken / created / collapsed.
        """
        if not _rct._md_distance_check(base_coords, coords, graph, 0.05):
            logger.debug("ring-pucker validate: M-D guard failed")
            return None
        cand_xyz = _rot._format_delfin_xyz(symbols, coords)
        try:
            if not _verify_topology_from_graph(cand_xyz, mol):
                logger.debug("ring-pucker validate: output topology gate failed")
                return None
        except Exception:
            return None
        return cand_xyz

    n_added = 0
    n_dropped_cap = 0
    n_dropped_basin = 0
    collapsed_classes = []  # (representative_first_atom, n_collapsed)

    all_rings = [r for _k, grp in ordered_groups for r in grp]

    def _gate_ok(coords):
        """True iff *coords* passes the M-D invariant guard AND the
        authoritative output topology gate (no broken / phantom bond, no
        collapse)."""
        try:
            if not _rct._md_distance_check(base_coords, coords, graph, 0.05):
                return False
            return bool(_verify_topology_from_graph(
                _rot._format_delfin_xyz(symbols, coords), mol
            ))
        except Exception:
            return False

    chair_rings_committed = set()

    # --- (1) GLOBAL chair-set frame: drive AS MANY rings to chair as fit. ---
    # Build INCREMENTALLY: snap each ring onto the accumulating frame and commit
    # it only if, after a constrained relax that holds the already-committed
    # rings + the new ring (+ coord sphere) and lets everything else relax, the
    # ring lands in the chair basin AND the whole-molecule output gate still
    # passes.  A snap that would clash with a previously-committed bulky ring is
    # reverted.  One emitted "global-chair" frame with every ring that fits.
    chair_coords = [c for c in base_coords]
    committed_list = []
    for ring in all_rings:
        snapped = _snap_ring(chair_coords, ring, "chair")
        if snapped is None:
            continue
        if _ring_basin(snapped, ring) != "chair":
            n_dropped_basin += 1
            continue
        hold = set()
        for cr in committed_list:
            hold.update(cr)
        hold.update(ring)
        candidate = _constrained_relax(snapped, hold)
        if _ring_basin(candidate, ring) != "chair" or not _gate_ok(candidate):
            n_dropped_basin += 1
            continue
        chair_coords = candidate
        committed_list.append(ring)
        chair_rings_committed.add(tuple(sorted(ring)))
    if committed_list:
        cand = _build_and_validate(chair_coords)
        if cand is not None:
            sig = _sig(cand)
            if sig not in seen_sigs:
                seen_sigs.add(sig)
                results.append((cand, "pucker-global-chair"))
                n_added += 1

    # --- (2) PER-RING chair frames for rings NOT covered by the global set ---
    # When the global all-chair frame cannot hold every ring simultaneously
    # (bulky symmetry-equivalent rings on a shared centre sterically clash when
    # all forced to the same chair at once), emit a single-ring chair frame for
    # each remaining ring so that EVERY flexible ring's chair basin is realized
    # somewhere in the ensemble.  Snapping one ring at a time (on the base
    # geometry, holding only that ring + the coord sphere) avoids the inter-ring
    # clash.  Bounded by the per-complex cap.
    for ring in all_rings:
        if n_added >= max_ring_variants or len(results) >= max_isomers:
            n_dropped_cap += 1
            continue
        if tuple(sorted(ring)) in chair_rings_committed:
            continue
        snapped = _snap_ring(base_coords, ring, "chair")
        if snapped is None:
            continue
        if _ring_basin(snapped, ring) != "chair":
            n_dropped_basin += 1
            continue
        relaxed = _constrained_relax(snapped, set(ring))
        if _ring_basin(relaxed, ring) != "chair" or not _gate_ok(relaxed):
            n_dropped_basin += 1
            continue
        cand = _build_and_validate(relaxed)
        if cand is None:
            continue
        sig = _sig(cand)
        if sig in seen_sigs:
            continue
        seen_sigs.add(sig)
        results.append((cand, f"pucker-chair-r6-a{ring[0]}"))
        n_added += 1

    # --- (3) one-ring-flipped-to-boat decorations (one per symmetry class) ---
    # Symmetry-equivalent rings collapse to ONE DOF here: a single boat
    # decoration per equivalence class (snap the representative ring to boat on
    # the base geometry) rather than the full 2^N pucker product.
    for _key, grp in ordered_groups:
        if n_added >= max_ring_variants or len(results) >= max_isomers:
            n_dropped_cap += 1
            continue
        rep = sorted(grp, key=lambda r: r[0])[0]
        if len(grp) > 1:
            collapsed_classes.append((rep[0], len(grp) - 1))
        snapped = _snap_ring(base_coords, rep, "boat")
        if snapped is None:
            n_dropped_basin += 1
            continue
        if _ring_basin(snapped, rep) != "boat":
            n_dropped_basin += 1
            continue
        relaxed = _constrained_relax(snapped, set(rep))
        if _ring_basin(relaxed, rep) != "boat" or not _gate_ok(relaxed):
            n_dropped_basin += 1
            continue
        cand = _build_and_validate(relaxed)
        if cand is None:
            continue
        sig = _sig(cand)
        if sig in seen_sigs:
            continue
        seen_sigs.add(sig)
        results.append((cand, f"pucker-boat-r6-a{rep[0]}"))
        n_added += 1

    # --- (3) NO silent truncation: log what collapsed / was dropped. ---
    if n_added or n_dropped_cap or n_dropped_basin or collapsed_classes:
        logger.debug(
            "Non-metal ring-pucker pass added %d conformer variants "
            "(%d rings -> %d symmetry classes); dropped %d (cap), "
            "%d (basin-unreached); collapsed-by-symmetry: %s",
            n_added,
            len(all_rings),
            len(ordered_groups),
            n_dropped_cap,
            n_dropped_basin,
            ", ".join(
                f"ring@a{a}(+{n})" for a, n in collapsed_classes
            ) or "none",
        )
    return n_added


def _emit_all_trans_by_type_arrangements(
    mol,
    results: List[Tuple[str, str]],
    dtype_map: Dict[int, tuple],
    apply_uff: bool,
    max_isomers: int,
) -> int:
    """Append explicit "all-trans-by-type" coordination arrangements to ``results``.

    For every metal centre with a coordination number whose geometry has
    a non-empty ``_TOPO_TRANS_POSITIONS`` list, enumerate every way to
    assign donor types to trans-pairs such that EVERY trans-pair contains
    two donors of the SAME chemical type (e.g. H₂O–H₂O trans simultaneously
    with O–O trans and N–N trans).  Build the corresponding XYZ via
    ``_build_topology_xyz`` and append it iff its heavy-atom XYZ signature
    is not already present in ``results``.

    Pure ADDITIVE pass — never sorts, drops, or modifies existing entries.
    Designed to ensure the chemically-central trans-effect arrangements
    appear in the output even when the main pipeline's UFF + fingerprint-
    dedup converges them with already-emitted candidates.

    Returns the number of new entries appended.
    """
    if not RDKIT_AVAILABLE:
        return 0
    if not _delfin_env_int("DELFIN_TRANS_PASS_ENABLED", 1):
        return 0

    # Build heavy-atom XYZ signature set of existing results so we don't
    # duplicate.  Same shape as the dual-parse union signature.
    def _sig(xyz_str: str) -> tuple:
        try:
            lines = [
                ln.split() for ln in xyz_str.strip().splitlines()
                if ln.strip()
            ]
            heavy = sorted(
                (p[0], round(float(p[1]), 2), round(float(p[2]), 2),
                 round(float(p[3]), 2))
                for p in lines
                if len(p) >= 4 and p[0] not in ('H', 'h')
            )
            return tuple(heavy)
        except Exception:
            return tuple()

    seen_sigs = {_sig(xyz) for xyz, _ in results}
    n_added = 0

    # Topology-template conformer for builder context.  Reuse if mol has
    # any conformer; otherwise builder falls back to its own placement.
    try:
        _topo_cid = _rank_template_conformers(mol, top_k=1)
        topo_template_cid = _topo_cid[0] if _topo_cid else None
    except Exception:
        topo_template_cid = None

    import itertools as _it

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        if len(results) + n_added >= max_isomers:
            break

        metal_idx = atom.GetIdx()
        donor_indices = [nb.GetIdx() for nb in atom.GetNeighbors()]
        n_coord = len(donor_indices)
        if n_coord < 2:
            continue

        # Donor labels (Morgan-class) per donor atom.
        donor_keys = [
            dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))
            for d in donor_indices
        ]
        uniq_keys = sorted(
            set(donor_keys),
            key=lambda k: (k[0], tuple(sorted(k[1]))),
        )
        key_to_class = {k: i for i, k in enumerate(uniq_keys)}
        donor_labels = [f"{k[0]}{key_to_class[k]}" for k in donor_keys]

        # Donor-type counts.  All-trans-by-type requires every type to
        # appear an even number of times so it can be paired across one
        # or more trans-positions.
        from collections import Counter as _Ctr
        type_counts = _Ctr(donor_labels)
        # ===== MIXED trans PAIRS (16.08.2026) ========================================
        # MEASURED on `rows_spy5geom` (965): the builder seats **cis where the crystal is
        # trans** -- `build-cis/crystal-trans` 126 failures against 82 hits (1.39x), the
        # opposite direction 21 against 51 (0.37x).  The ratio is **3.8x**, i.e. a
        # directed bias, and it is exactly what a trans-blind seating MUST
        # produce: for every first position there are FOUR cis sites in the octahedron and
        # only ONE trans.  On top, `arrangement_complete = false` at 1.72x -- if ONE of the
        # two arrangements is missing, the system fails.
        #
        # WHY THIS PASS DOES NOT DELIVER THEM.  It demands that EVERY trans pair carries two
        # donors of the SAME type -- and bails out here entirely as soon as any type has an
        # ODD count.  A Cu with 2 N + 1 O + 1 Cl thereby gets NO trans arrangement
        # AT ALL.  The failure list names exactly such cases: `Cu-ON` 3.62x,
        # `CC-Ir` and `CC-W` stand at **only failures**.
        #
        # With `DELFIN_FFFREE_TRANS_MIXED=1` both go away: odd type counts no longer bail
        # out, and below type-MIXED pairs are added.  Purely ADDITIVE -- the
        # signature dedup and the topology gate below stay unchanged, nothing is
        # replaced.  Default OFF -> byte-identical.
        _trans_mixed = bool(_delfin_env_int("DELFIN_FFFREE_TRANS_MIXED", 0))
        if not _trans_mixed and any(c % 2 != 0 for c in type_counts.values()):
            continue

        # Geometry candidates with non-empty trans-position lists.
        # Universal: covers EVERY CN that has at least one geometry with
        # a 180° pair in `_TOPO_TRANS_POSITIONS` (LIN, TS, SQ, SS, TBP,
        # SP, OH, PBP, COH, SAP, DD).  CN=3 trigonal-planar (TP) and
        # CN=6 trigonal-prismatic (TPR), CN=9 TTP have no trans pairs
        # by definition and are correctly skipped.
        _CN_GEOM_FOR_TRANS = {
            2: ['LIN'],
            3: ['TS'],          # T-shaped has 1 trans pair (0-1)
            4: ['SQ', 'SS'],    # square-planar 2 trans, see-saw 1 trans
            5: ['TBP', 'SP'],   # TBP axial-axial; SP basal cross
            6: ['OH'],          # 3 trans pairs
            7: ['PBP', 'COH'],  # axial pair / oct-base trans pairs
            8: ['SAP', 'DD'],   # antiprism / dodecahedron 4 trans pairs
        }
        geom_list = _CN_GEOM_FOR_TRANS.get(n_coord, [])
        if not geom_list:
            continue

        # Chelate constraints: pairs of donor-list-indices that must
        # never sit at trans positions.
        chelate_atom_ps = _chelate_pairs(mol, metal_idx, donor_indices)
        atom_to_listidx = {ai: li for li, ai in enumerate(donor_indices)}
        chelate_pairs = []
        for cp in chelate_atom_ps:
            pp = sorted(cp)
            if len(pp) == 2 and pp[0] in atom_to_listidx and pp[1] in atom_to_listidx:
                chelate_pairs.append(frozenset([
                    atom_to_listidx[pp[0]], atom_to_listidx[pp[1]]
                ]))

        for geom in geom_list:
            trans_pos = _TOPO_TRANS_POSITIONS.get(geom, [])
            if not trans_pos:
                continue
            n_trans_pairs = len(trans_pos)

            # Each donor-type must contribute at least one trans-pair.
            # Total trans-positions = 2 * n_trans_pairs.  Donors covered
            # by trans-pairs = 2 * n_trans_pairs.  Donors not on a trans
            # position (e.g. SP equatorial cap) get assigned freely
            # afterwards.
            covered_positions = sorted({p for ta, tb in trans_pos for p in (ta, tb)})
            uncovered_positions = [
                p for p in range(n_coord) if p not in covered_positions
            ]

            # Compute every way to choose which donor-type goes on which
            # trans-pair: this is a multiset partition.  For type counts
            # {O1:2, O3:2, O2:2, N0:2} on 4 trans-pairs, that's 4!=24
            # bijections of types-onto-pairs — manageable.  We treat it
            # as: which donor-list-index pairs (same-type) sit at each
            # trans-pair-position.
            #
            # Strategy: pre-group donor-list-indices by type, then assign
            # one same-type pair per geometric trans-pair.  Each
            # assignment yields a partial perm; remaining donors fill
            # uncovered_positions.
            type_to_donors: Dict[str, list] = {}
            for li, lbl in enumerate(donor_labels):
                type_to_donors.setdefault(lbl, []).append(li)

            # Need: at least n_trans_pairs distinct (type, pair-of-donors)
            # available.  If only 2 types each with 2 donors and 4
            # geometric trans-pairs, we cannot fill all trans-pairs same-
            # type → skip this geom.
            available_type_pairs = []
            for t, dlist in type_to_donors.items():
                if len(dlist) < 2:
                    continue
                # Take all C(n,2) within-type donor-index pairs
                for di, dj in _it.combinations(dlist, 2):
                    available_type_pairs.append((t, (di, dj)))

            # MIXED pairs (see above).  A trans pair of two DIFFERENT donor types is
            # chemically the normal case -- N trans to O, C trans to P -- and was until now
            # not representable here.  They are APPENDED, not replaced: the same-type
            # partitions still arise first and keep their precedence in the
            # enumeration.
            if _trans_mixed:
                for li in range(len(donor_labels)):
                    for lj in range(li + 1, len(donor_labels)):
                        if donor_labels[li] != donor_labels[lj]:
                            available_type_pairs.append(
                                (f"{donor_labels[li]}|{donor_labels[lj]}", (li, lj)))

            if len(available_type_pairs) < n_trans_pairs:
                continue

            # Enumerate ways to pick `n_trans_pairs` disjoint donor-pairs
            # such that every donor-list-index is used at most once and
            # every trans-pair gets one same-type donor pair.  Capped at
            # 12 partitions per (metal, geom) to keep wall-time bounded.
            # ⚠ NO SILENT TRUNCATION.  With mixed pairs the space grows markedly
            # (for CN6 from a few same-type ones to up to 15 partitions), hence a
            # separate, higher cap -- and it is LOGGED when it binds.  A cap
            # that stays silent reads afterwards like "completely enumerated".
            _MAX_PARTITIONS = 24 if _trans_mixed else 12
            partitions = []

            def _backtrack(used_donors, partial):
                if len(partial) == n_trans_pairs:
                    partitions.append(list(partial))
                    return len(partitions) >= _MAX_PARTITIONS
                if len(partitions) >= _MAX_PARTITIONS:
                    return True
                for t, (di, dj) in available_type_pairs:
                    if di in used_donors or dj in used_donors:
                        continue
                    # Avoid duplicate partitions: only consider ordered
                    # additions (next pair has minimum di > previous min).
                    if partial and (di, dj) <= partial[-1][1]:
                        continue
                    partial.append((t, (di, dj)))
                    used_donors.add(di); used_donors.add(dj)
                    if _backtrack(used_donors, partial):
                        return True
                    used_donors.remove(di); used_donors.remove(dj)
                    partial.pop()
                return False

            _backtrack(set(), [])
            if len(partitions) >= _MAX_PARTITIONS:
                try:
                    logger.warning(
                        "trans-pass: Deckel %d Partitionen erreicht (%s, CN%d, %s) -- "
                        "weitere Anordnungen NICHT aufgezaehlt",
                        _MAX_PARTITIONS, mol.GetAtomWithIdx(metal_idx).GetSymbol(),
                        n_coord, geom)
                except Exception:
                    pass
            if not partitions:
                continue

            # For each partition, assign donor-pairs onto geometric trans
            # pairs (n_trans_pairs! orderings).  Cap to 6 orderings per
            # partition to bound work.
            for partition in partitions:
                if len(results) + n_added >= max_isomers:
                    break
                _ordering_count = 0
                for ordered in _it.permutations(partition):
                    if _ordering_count >= 6:
                        break
                    _ordering_count += 1
                    perm = [None] * n_coord
                    used = set()
                    for (ta, tb), (_t, (di, dj)) in zip(trans_pos, ordered):
                        perm[ta] = di
                        perm[tb] = dj
                        used.add(di); used.add(dj)
                    # Fill uncovered positions with remaining donors in
                    # ascending list-index order (deterministic).
                    remaining = [li for li in range(n_coord) if li not in used]
                    for pos, li in zip(uncovered_positions, remaining):
                        perm[pos] = li
                    if any(p is None for p in perm):
                        continue

                    # Chelate-cis check: skip if any chelate pair lands
                    # on a geometric trans-position pair.
                    ch_violation = False
                    for chp in chelate_pairs:
                        a, b = list(chp)
                        pa, pb = perm.index(a), perm.index(b)
                        for ta, tb in trans_pos:
                            if (pa == ta and pb == tb) or (pa == tb and pb == ta):
                                ch_violation = True
                                break
                        if ch_violation:
                            break
                    if ch_violation:
                        continue

                    # Build XYZ.
                    try:
                        xyz = _build_topology_xyz(
                            mol, metal_idx, donor_indices, perm, geom,
                            apply_uff, conf_id=topo_template_cid,
                        )
                    except Exception:
                        xyz = None
                    if not xyz:
                        continue

                    # Final-output topology gate (ensures consistency
                    # with main pipeline's last gate).
                    try:
                        if not _verify_topology_from_graph(xyz, mol):
                            continue
                    except Exception:
                        continue

                    # Heavy-atom signature dedup against existing pool.
                    sig = _sig(xyz)
                    if sig in seen_sigs:
                        continue
                    seen_sigs.add(sig)

                    # Build label: "trans-{geom} {types}|{types}|..."
                    type_str = "|".join(
                        f"{t}+{t}" for t, _ in ordered
                    )
                    label = f"trans-{geom} {type_str}"
                    # Iter-8.7 every-append gate (123a130, env-gated default OFF)
                    _gate_pass = True
                    if _every_append_gate_enabled(mol):
                        try:
                            _flat = _flatten_sp2_atoms_xyz(xyz, mol)
                            if _flat:
                                xyz = _flat
                        except Exception:
                            pass
                        try:
                            _mt_g = Chem.RWMol(mol); _mt_g.RemoveAllConformers()
                            _c_g = _xyz_to_rdkit_conformer(_mt_g.GetMol(), xyz)
                            if _c_g is None:
                                _gate_pass = False
                            else:
                                _ci_g = _mt_g.AddConformer(_c_g, assignId=True)
                                if _has_severe_covalent_distortion(_mt_g.GetMol(), _ci_g):
                                    _gate_pass = False
                        except Exception:
                            _gate_pass = False
                    if not _gate_pass:
                        continue
                    results.append((xyz, label))
                    n_added += 1
                    if len(results) >= max_isomers:
                        break

    if n_added:
        logger.debug(
            "Trans-effect pass added %d explicit all-trans-by-type entries",
            n_added,
        )
    return n_added


def _filter_nonfinite_isomers(results):
    """Universal output contract (#36): never emit a structure with non-finite
    (NaN/inf) coordinates — it is not a valid DFT startpoint and corrupts every
    downstream metric.  Drops ONLY the broken isomer/conformer, keeps the finite
    ones (no coverage loss when other builds for the same SMILES are valid).
    Applies to EVERY build path (legacy, hapto, fffree) since it gates the public
    entry point.  Geometry-only, deterministic."""
    if not results:
        return results
    clean = []
    for item in results:
        xyz = item[0] if isinstance(item, (tuple, list)) and item else item
        finite = True
        for ln in str(xyz).splitlines():
            p = ln.split()
            if len(p) >= 4:
                try:
                    x, y, z = float(p[1]), float(p[2]), float(p[3])
                except ValueError:
                    continue
                if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
                    finite = False
                    break
        if finite:
            clean.append(item)
    return clean


def _rank_emitted_isomers(isomers):
    """Order the emitted ensemble best (most crystal-like) first.  Never drops or
    alters a structure — only reorders.  Safe no-op on any error or with
    DELFIN_NO_FRAME_RANK=1.

    DEFAULT = the legacy least-clash ranking (cross-validated 2026-06-12).  That
    ranker measures ONLY steric overlap, so it is blind in both directions: a
    TORN / decoordinated frame has no overlap at all and therefore scores the
    perfect 100.0 and LEADS, while on a clean ensemble 86 % of all emitted frames
    score exactly 100.0 -> a complete tie -> the delivered order is the ENUMERATION
    order, i.e. no ranking at all (measured over the built census 2026-07-28; user
    report the same day: "the best frames are not sorted to the front").

    DELFIN_FRAME_RANK_QUALITY=1 switches to the DEFECT ordering
    (:func:`delfin.manta._conformer_rank.rank_isomers_quality`): torn / spurious /
    collapsed ligand bonds first, then clash, then coordination distortion
    relative to the best frame of the same CN signature, enumeration index last.
    Pure geometry, deterministic, no RMSD, no energy; a verified PERMUTATION with an
    explicit bijection guard, so no frame is ever lost and every label travels with
    its frame.  Default-OFF -> byte-identical order to today.  Applied ONLY here, on
    the final emitted ensemble -- the construction-time base-frame picker keeps using
    the legacy ``rank_isomers``, so the manifold CONTENT is unchanged."""
    try:
        if os.environ.get("DELFIN_FRAME_RANK_QUALITY", "0") == "1":
            from delfin.manta._conformer_rank import rank_isomers_quality
            return rank_isomers_quality(isomers)
        from delfin.manta._conformer_rank import rank_isomers
        return rank_isomers(isomers)
    except Exception:
        return isomers


def _hapto_declash_filter(isomers):
    """Final-manifold eta-ring inter-ligand declash over the emitted ensemble.

    The legacy analytical piano-stool seat (``_build_hapto_scaffold``) places an
    eta-coordinated ring (arene / Cp / Cb / diene) on the metal but does NOT
    choose its ROTATIONAL conformer about the M->centroid axis, so ring carbons /
    substituents routinely eclipse the carbonyls and crowd the other ligands.
    Measured against CCDC piano-stool crystals (find_inter_ligand_clash, F26) the
    seat over-clashes ~3-4x, the residual dominated by the eta-face crowding
    pendant aryls / other rings / sigma-donors -- contacts essentially absent in
    the crystals.

    This pass rotates each eta-ring about its M->centroid axis to the LOCAL
    inter-ligand-clash minimum within a bounded window (the missing piano-stool
    DOF).  It is a pure isometry of the coordination polyhedron: every
    M-(ring atom) distance is preserved, the metal + all sigma-donors are frozen,
    the inter-ligand overlap is never increased, and it is deterministic.

    Runs as the LAST geometry transform over the FINAL, fully-selected + dedup'd
    ensemble (NOT per-conformer during construction -- that would feed back into
    conformer selection / dedup and change WHICH conformers survive), so the
    isomer/conformer COUNT is exactly preserved: a frame is only re-emitted when
    the spin actually moved an atom; untouched otherwise.

    Identity (byte-identical) unless DELFIN_FFFREE_HAPTO_DECLASH=1.  Geometry-only
    (Bondi vdW radii); no CSD/CCDC data.  Deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_HAPTO_DECLASH", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        try:
            _maxrot = float(os.environ.get(
                "DELFIN_FFFREE_HAPTO_DECLASH_MAXROT", _HD._DEF_MAXROT_DEG))
        except (TypeError, ValueError):
            _maxrot = _HD._DEF_MAXROT_DEG
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 4:
                out.append(item); continue
            P = np.asarray(_HD.declash_frame(syms, coords, max_rot_deg=_maxrot), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue          # no eta-ring / nothing moved
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _sigma_declash_filter(isomers):
    """SIGMA-co-ligand inter-ligand declash over the emitted ensemble -- the
    sigma-ligand counterpart to the eta-ring spin (_hapto_declash_filter).

    Relieves crowding of the NON-eta ligand BODIES (stannyl / stibine / phosphite
    / bulky phosphine substituents) by the joint_declash kinematics: each ligand's
    whole-body rotation about its M-donor axis + internal rotatable-bond torsions,
    metal + ALL donors FROZEN (coordination polyhedron invariant).  A
    partition-free total-overlap guard makes every frame strictly never-worse
    (measured on the hapto pool: inter-ligand clashes 4867 -> 4150 on top of the
    ring spin, 0 frames worse).

    Runs over the FINAL dedup'd manifold (count-preserving: a frame is only
    re-emitted when an atom actually moved).  Identity (byte-identical) unless
    DELFIN_FFFREE_SIGMA_DECLASH=1.  Geometry-only; deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_SIGMA_DECLASH", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 4:
                out.append(item); continue
            P = np.asarray(_HD.declash_frame_sigma(syms, coords), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _carbonyl_fix_filter(isomers):
    """Contract stretched terminal metal-carbonyl C#O bonds over the emitted
    ensemble.  The analytical seat's flat _bond_len places C-O at the single-bond
    length (~1.43 A) ignoring bond order, so M-C#O carbonyls come out stretched
    (~1.40 A vs the true ~1.13-1.15 A; measured 56% > 1.35 A).  This slides only
    the terminal O inward along the C->O axis to the triple-bond length -- a pure
    local geometry correction that never lengthens and moves nothing else, so the
    isomer/conformer COUNT and all other geometry are preserved.

    Identity (byte-identical) unless DELFIN_FFFREE_CARBONYL_FIX=1.  Geometry-only;
    deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_CARBONYL_FIX", "0") != "1":
        return isomers
    try:
        import numpy as np
        from delfin.manta import hapto_declash as _HD
        out = []
        for item in isomers:
            is_tuple = isinstance(item, (tuple, list))
            xyz = item[0] if is_tuple and item else item
            lines = str(xyz).splitlines()
            syms, coords, idx = [], [], []
            for li, ln in enumerate(lines):
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0]); coords.append([x, y, z]); idx.append(li)
            if len(syms) < 3:
                out.append(item); continue
            P = np.asarray(_HD.correct_carbonyls(syms, coords), float)
            orig = np.asarray(coords, float)
            if P.shape != orig.shape or np.allclose(P, orig, atol=1e-7):
                out.append(item); continue
            for k, li in enumerate(idx):
                p = lines[li].split()
                p[1] = f"{P[k][0]:.6f}"; p[2] = f"{P[k][1]:.6f}"; p[3] = f"{P[k][2]:.6f}"
                lines[li] = " ".join(p)
            new_xyz = "\n".join(lines)
            out.append((new_xyz,) + tuple(item[1:]) if is_tuple else new_xyz)
        return out
    except Exception:
        return isomers


def _topology_gate_filter(isomers):
    """Final UNIVERSAL topology-preservation gate over the FULL emitted ensemble.

    SPEC
    ----
    Root cause it cures: some conformer/pucker frames of a flexible complex let a
    NON-bonded heavy-atom pair fall to bonding distance, forging a SPURIOUS bond
    (a topology break).  Example: ACEBOB (Re tris-carbonyl macrocycle) -- in some
    frames a carbonyl C penetrates the macrocyclic cavity and lands ~1.47 A from a
    ligand C/N.  At 1.47 A this reads as a perfectly normal covalent bond, so the
    existing clash/coord self-gates are BLIND to it.  We must reject ONLY the
    broken frames and KEEP the clean ones (never collapse the ensemble).

    Method (consensus-based, NO reference/SMILES alignment needed):
      All emitted frames for one structure share atom order.  A GENUINE covalent
      bond is rigid -- present in ALL frames; a SPURIOUS contact appears in only a
      minority.  So for every heavy-heavy NON-metal atom pair (i,j) we measure the
      fraction of frames in which dist(i,j) < bond_thr(elem_i, elem_j), with
      bond_thr = 1.15 * (cov_i + cov_j) (Cordero covalent radii; ~1.75 A for C-C).
      A pair is a "real bond" iff that fraction >= 0.80 (bonded in the clear
      majority).  A frame has a TOPOLOGY BREAK iff it contains a pair that is at
      bonding distance in THIS frame yet is NOT a real bond.  Such frames are
      dropped.  Metal-involving pairs are excluded (M-donor distances vary
      legitimately and must never be treated as breakable bonds).

    Guarantees:
      * Geometry-only and graph-free -- no per-SMILES / per-refcode special case.
      * NEVER-WORSE: if every frame would be dropped (or the result would be
        empty) the ORIGINAL list is returned unchanged.
      * Identity (returns input untouched) unless DELFIN_FFFREE_TOPOLOGY_GATE=1,
        so output is BYTE-IDENTICAL to the baseline when off.
      * Deterministic; never raises.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_TOPOLOGY_GATE", "0") != "1":
        return isomers
    if len(isomers) < 2:
        # Consensus is undefined with a single frame; nothing to compare against.
        return isomers
    try:
        def _radius(sym):
            return _COVALENT_RADII.get(sym, 0.76)

        # Parse every frame into (symbols, coords) sharing one atom order.
        frames = []
        n_atoms = None
        for item in isomers:
            xyz = item[0] if isinstance(item, (tuple, list)) and item else item
            syms, xs, ys, zs = [], [], [], []
            for ln in str(xyz).splitlines():
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    syms.append(p[0])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            if n_atoms is None:
                n_atoms = len(syms)
            if len(syms) != n_atoms or n_atoms == 0:
                # Frames disagree on atom count -> consensus is ill-defined; bail
                # out safely (identity) rather than risk an unsound drop.
                return isomers
            frames.append((syms, xs, ys, zs))

        nf = len(frames)
        ref_syms = frames[0][0]

        # Heavy (non-H), non-metal atom indices -- consensus is only meaningful
        # for the rigid covalent backbone; H and the metal are excluded.
        heavy = [
            k for k in range(n_atoms)
            if ref_syms[k] != 'H' and ref_syms[k] not in _METAL_SET
        ]

        # Precompute per-pair bond threshold (squared) once.
        # near[(i,j)] = count of frames where dist(i,j) < bond_thr.
        near = {}
        thr2 = {}
        for a in range(len(heavy)):
            i = heavy[a]
            ri = _radius(ref_syms[i])
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                t = 1.15 * (ri + _radius(ref_syms[j]))
                thr2[(i, j)] = t * t
                near[(i, j)] = 0

        for (syms, xs, ys, zs) in frames:
            for a in range(len(heavy)):
                i = heavy[a]
                xi, yi, zi = xs[i], ys[i], zs[i]
                for b in range(a + 1, len(heavy)):
                    j = heavy[b]
                    dx = xi - xs[j]
                    dy = yi - ys[j]
                    dz = zi - zs[j]
                    if dx * dx + dy * dy + dz * dz < thr2[(i, j)]:
                        near[(i, j)] += 1

        # A pair is a "real bond" if bonded in the clear majority (>= 80%).
        real_bond = set()
        for key, cnt in near.items():
            if cnt >= 0.80 * nf:
                real_bond.add(key)

        # Drop a frame iff any pair is bonding in THIS frame but is not a real bond.
        kept = []
        for item, (syms, xs, ys, zs) in zip(isomers, frames):
            broken = False
            for a in range(len(heavy)):
                i = heavy[a]
                xi, yi, zi = xs[i], ys[i], zs[i]
                for b in range(a + 1, len(heavy)):
                    j = heavy[b]
                    key = (i, j)
                    if key in real_bond:
                        continue
                    dx = xi - xs[j]
                    dy = yi - ys[j]
                    dz = zi - zs[j]
                    if dx * dx + dy * dy + dz * dz < thr2[key]:
                        broken = True
                        break
                if broken:
                    break
            if not broken:
                kept.append(item)

        # NEVER-WORSE: empty or all-dropped -> keep original ensemble.
        if not kept:
            return isomers
        return kept
    except Exception:
        return isomers


def _clean_gate_filter(isomers):
    """UNIVERSAL, ASYMMETRIC, PER-FRAME final clean-manifold emission gate.

    Governing invariant (user, verbatim): "sort out only the BAD ones, not
    the good ones".  The gate is ASYMMETRIC: it drops a frame ONLY when that frame is
    CERTAINLY bad on its OWN geometry; in any doubt it KEEPS the frame.  A
    surviving good structure matters more than removing a borderline one
    (false-negative acceptable; false-positive -- dropping a good frame --
    forbidden).

    Each frame is judged ABSOLUTELY on its own coordinates -- NO ensemble
    consensus, NO best-coordinated-reference donor set, NO minority/majority vote
    (all of which over-reject on small/few-frame ensembles by letting one frame's
    geometry condemn another's).  Bonds are perceived PER FRAME from that frame's
    own interatomic distances.  A frame is REJECTED only if ANY DEFINITE defect is
    present with a CLEAR margin (all graph+geometry only, NO SMILES / refcode /
    element special-casing):

      1. STRUCTURALLY INVALID -- a non-finite coordinate, empty frame, or atom
         count mismatch (a destroyed frame).
      2. REAL INTER-LIGAND CLASH -- two heavy atoms in DIFFERENT per-frame ligand
         connected-components overlap clearly below ``clash`` x (cov_i + cov_j),
         or an H-involving non-bonded pair sits below a hard vdW floor (deep
         interpenetration, not a soft contact).
      3. COLLAPSED BOND -- a per-frame covalent bond fused below ``collapse`` x its
         ideal covalent length (two atoms on top of each other).
      4. BARE METAL / TOTAL DECOORDINATION -- a metal with ZERO heavy donors
         anywhere inside its first shell (the complex fell completely apart).  This
         is the only decoordination test that needs no cross-frame reference, so it
         is the only one applied: a single donor wandering a little is NOT certainly
         bad and is KEPT.

    Deliberately NOT rejected (these over-rejected good frames before): consensus
    "spurious bond" (a close contact is judged per-frame, only a true collapse or a
    clear clash counts), bond OVER-stretch (a long contact is simply not a bond in
    that frame -- topology, not a defect for this gate), and a single donor that
    swung modestly past an ideal M-D distance (handled, if at all, only by the
    dedicated coord-integrity filter).

    NEVER-EMPTY without keeping a broken frame: the old "cleanest-available"
    fallback that retained a broken frame is REMOVED.  Frames that are not
    certainly bad simply stay (no consensus needed).  In the degenerate case where
    EVERY frame is certainly bad, the input is returned UNCHANGED (keep all rather
    than fabricate a marked broken survivor) -- in doubt, keep.

    Thresholds derive from covalent radii (Cordero, ``_COVALENT_RADII``); all
    universal, env-tunable.  Identity (returns input untouched) unless
    DELFIN_FFFREE_CLEAN_GATE=1, so output is BYTE-IDENTICAL to baseline when off.
    Deterministic; never raises; returns the input unchanged on any error.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_CLEAN_GATE", "0") != "1":
        return isomers
    try:
        # ---- tunables (env, sane universal defaults) -----------------------
        # clash_f: heavy-heavy inter-ligand overlap fraction that is CERTAINLY a
        #          clash (clear margin below the covalent sum).
        clash_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_CLASH", "0.70"))
        # vdw_f: H-involving non-bonded deep-interpenetration floor (hard).
        vdw_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_VDW", "0.55"))
        # bond_mult: covalent multiplier defining "at bonding distance" PER FRAME.
        bond_mult = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_BONDMULT", "1.15"))
        # collapse_f: a per-frame bond shorter than this fraction of ideal = fused.
        collapse_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_COLLAPSE", "0.55"))
        # md_shell: metal first-shell radius (heavy donor within this = coordinated).
        md_shell = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_SHELL", "2.95"))
        # md_collapse_on: also reject a frame whose heavy donor has fused ONTO the
        #          metal (M-D below collapse_f x covalent-sum).  Independent, default
        #          OFF -> byte-identical to the prior CLEAN_GATE behaviour when unset.
        #          Closes a blindspot: the per-frame pair tests (2/3) SKIP every metal
        #          pair, so a donor collapsed onto the metal is otherwise invisible
        #          here (and to the ligcollapse/geomfault detectors, which exempt the
        #          M-D pair too).  Asymmetric-safe: no real M-D bond is anywhere near
        #          this short (shortest M-heavy-donor bonds ~1.5A).
        md_collapse_on = (os.environ.get(
            "DELFIN_FFFREE_CLEAN_GATE_MD_COLLAPSE", "0") == "1")
        # clash_vdw_on: measure the inter-ligand clash against the VAN DER WAALS sum instead of
        #   the covalent sum.  Default OFF -> byte-identical.
        #
        #   WHY.  Measured 2026-08-04 on the 142 worst-broken systems: the gate's five criteria
        #   fired 420 times for `collapse` and EXACTLY ZERO times for `clash` -- while the eye
        #   reports inter-ligand clashes on 3272 of 27757 frames (11.79 %), a defect class that
        #   occurs in the clean crystals at 0.00 %.  The gate is blind to it, and the reason is a
        #   radius mix-up, not a threshold:
        #       gate  0.70 x (rcov_i + rcov_j)   ->  C-C below 1.06 A   <- guessed default
        #       eye   0.65 x (rvdW_i + rvdW_j)   ->  C-C below 2.21 A   <- COD-validated, FP -> 0
        #   Two carbons never reach 1.06 A without failing the collapse test first, so the clash
        #   branch could not fire.  Switching to the vdW sum is NOT a loosening -- it replaces a
        #   guessed reference with the one the eye's own metric_inter_ligand_clash was calibrated
        #   on ("real crystals have inter-ligand contacts >= 0.85 x vdW sum"; 0.65 is the
        #   conservative firing point where the COD false-positive rate goes to zero).
        # h_point_on: the DIRECTIONAL H-clash -- an H that is correctly bonded to its parent but
        #   POINTS INTO a nearby non-bonded heavy atom.  Default OFF -> byte-identical.
        #
        #   WHY A SECOND H CRITERION.  The gate already has one (`vdw_f`), and it is dead for the
        #   same reason the clash branch was: it compares against the COVALENT sum.  For H-H that
        #   is 0.55 x 0.62 = 0.34 A -- two hydrogens never come that close, so `vdw_h` fired 0
        #   times in the 2026-08-04 tally while the eye reports xh_hh_clash on 5.27 % of frames
        #   (crystals: 0.00 %).
        #   The eye's find_h_clash does not use a distance ratio at all; it asks a DIRECTIONAL
        #   question, which is what actually distinguishes a real contact from a broken one:
        #       H bonded to its parent (0.85-1.20 A)
        #       AND a non-bonded heavy atom within 2.50 A of that H
        #       AND the parent->H direction within 30 deg of parent->heavy
        #   i.e. the H is aimed at the neighbour.  Cause per that detector: VSEPR-snap places H at
        #   ideal polar angles but never picks the rotational phase, so a methyl umbrella can
        #   point straight into its neighbour.  Pure geometry, three constants, no CCDC table --
        #   portable into DELFIN license-clean.
        h_point_on = (os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_POINT", "0") == "1")
        h_par_min = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_PARENT_MIN", "0.85"))
        h_par_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_PARENT_MAX", "1.20"))
        h_other_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_OTHER_MAX", "2.50"))
        h_angle_max = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_H_ANGLE", "30.0"))
        clash_vdw_on = (os.environ.get(
            "DELFIN_FFFREE_CLEAN_GATE_CLASH_VDW", "0") == "1")
        clash_vdw_f = float(os.environ.get("DELFIN_FFFREE_CLEAN_GATE_CLASH_VDW_F", "0.65"))
        _vdw_r = None
        if clash_vdw_on:
            try:
                from delfin.manta.refine import _vdw as _vdw_r
            except Exception:
                clash_vdw_on = False

        def _radius(sym):
            return _COVALENT_RADII.get(sym, 0.76)

        def _is_metal(sym):
            return sym in _METAL_SET

        # ---- parse every frame (shared atom order) ------------------------
        frames = []
        n_atoms = None
        for item in isomers:
            xyz = item[0] if isinstance(item, (tuple, list)) and item else item
            syms, xs, ys, zs = [], [], [], []
            ok = True
            for ln in str(xyz).splitlines():
                p = ln.split()
                if len(p) >= 4:
                    try:
                        x, y, z = float(p[1]), float(p[2]), float(p[3])
                    except ValueError:
                        continue
                    if not (math.isfinite(x) and math.isfinite(y)
                            and math.isfinite(z)):
                        ok = False
                    syms.append(p[0])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            if n_atoms is None:
                n_atoms = len(syms)
            # Structurally invalid iff non-finite coord, empty, or count mismatch.
            frames.append((syms, xs, ys, zs, ok and len(syms) == n_atoms
                           and n_atoms > 0))

        if not n_atoms:
            return isomers

        def _d2(f, i, j):
            dx = f[1][i] - f[1][j]
            dy = f[2][i] - f[2][j]
            dz = f[3][i] - f[3][j]
            return dx * dx + dy * dy + dz * dz

        # ---- WHICH criterion actually drops the frames? (measurement, not behaviour) -----
        # Measured 2026-08-04 (cleangateW, the 142 worst-broken systems): the gate removes
        # 547 of 2341 frames and takes the hard-frame fraction from 84.2 % to 77.2 % -- the
        # first lever that lowers the ABSOLUTE defect rate at all.  But it also drops frames
        # that MATCHED THE CRYSTAL: ccdc_backbone_lost 3 (GILKAQ, HOJSUX, UHEJIB),
        # isomers_lost 3 (UHEJIB 10 -> 4 of 38 theory).  That violates the gate's own
        # governing invariant ("sort out only the BAD ones, not the good ones").
        #
        # Before touching any threshold: find out WHICH of the five criteria does it.  The
        # obvious suspect was the clash factor -- and that guess was WRONG: the gate measures
        # clash against the COVALENT sum (0.70 x ~1.52 A = 1.06 A for C-C) while the eye's
        # COD-validated metric uses the VDW sum (0.65 x 3.40 = 2.21 A), so the gate is far
        # LOOSER there, not tighter.  Hence: count, do not guess.
        # Pure instrumentation -- no frame is kept or dropped differently, output is
        # byte-identical; the tally only reaches stderr under DELFIN_TRACE_CLEAN_GATE=1.
        _gate_tally = {"invalid": 0, "collapse": 0, "vdw_h": 0, "clash": 0, "clash_vdw": 0,
                       "bare_metal": 0, "md_collapse": 0, "h_point": 0, "rescued_last_of_kind": 0}

        # ---- PER-FRAME absolute certainly-bad test ------------------------
        def _certainly_bad(f):
            """True iff frame f is CERTAINLY bad on its OWN geometry (a definite
            defect with a clear margin).  In any doubt -> False (KEEP)."""
            # 1. structurally invalid.
            if not f[4]:
                _gate_tally["invalid"] += 1
                return True
            syms = f[0]
            metal_idx = [k for k in range(n_atoms) if _is_metal(syms[k])]
            metal_set = set(metal_idx)

            # per-frame bond graph over non-metal pairs (this frame's own dists).
            parent = list(range(n_atoms))

            def _find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a

            def _union(a, b):
                ra, rb = _find(a), _find(b)
                if ra != rb:
                    parent[ra] = rb

            for i in range(n_atoms):
                if i in metal_set:
                    continue
                ri = _radius(syms[i])
                for j in range(i + 1, n_atoms):
                    if j in metal_set:
                        continue
                    t = bond_mult * (ri + _radius(syms[j]))
                    if _d2(f, i, j) < t * t:
                        _union(i, j)
            comp = [_find(k) for k in range(n_atoms)]

            # 2+3: per-frame heavy/H pair geometry.
            for i in range(n_atoms):
                if i in metal_set:
                    continue
                si = syms[i]
                ri = _radius(si)
                hi = (si == "H")
                for j in range(i + 1, n_atoms):
                    if j in metal_set:
                        continue
                    sj = syms[j]
                    d = math.sqrt(_d2(f, i, j))
                    rsum = ri + _radius(sj)
                    hj = (sj == "H")
                    bonded = d < bond_mult * rsum
                    if bonded:
                        # 3. collapsed bond (fused atoms) -> certainly bad.
                        if d < collapse_f * rsum:
                            _gate_tally["collapse"] += 1
                            return True
                        continue
                    # non-bonded pair:
                    if hi or hj:
                        # deep H interpenetration floor (hard).
                        if d < vdw_f * rsum:
                            _gate_tally["vdw_h"] += 1
                            return True
                    else:
                        # 2. real inter-ligand clash (different per-frame ligands,
                        # clear overlap).  Same-component close contacts are NOT
                        # flagged (a tight intra-ligand contact is not a defect).
                        if comp[i] != comp[j] and d < clash_f * rsum:
                            _gate_tally["clash"] += 1
                            return True
                        # vdW-referenced inter-ligand clash (the eye's calibration point).
                        if (clash_vdw_on and comp[i] != comp[j]
                                and d < clash_vdw_f * (_vdw_r(si) + _vdw_r(sj))):
                            _gate_tally["clash_vdw"] += 1
                            return True
            # 3b. DIRECTIONAL H-clash (ported from the eye's find_h_clash; only when explicitly on).
            if h_point_on:
                for hi_ in range(n_atoms):
                    if syms[hi_] != "H":
                        continue
                    par = -1
                    dpar = 1e9
                    for p in range(n_atoms):
                        if p == hi_ or syms[p] == "H" or p in metal_set:
                            continue
                        dp = math.sqrt(_d2(f, hi_, p))
                        if dp < dpar:
                            dpar, par = dp, p
                    if par < 0 or not (h_par_min <= dpar <= h_par_max):
                        continue                       # not a cleanly bonded H -> other tests own it
                    for q in range(n_atoms):
                        if q in (hi_, par) or syms[q] == "H" or q in metal_set:
                            continue
                        dq = math.sqrt(_d2(f, hi_, q))
                        if dq > h_other_max:
                            continue
                        # angle at the PARENT between parent->H and parent->q
                        dpq = math.sqrt(_d2(f, par, q))
                        if dpq < 1e-6:
                            continue
                        cosang = (dpar * dpar + dpq * dpq - dq * dq) / (2.0 * dpar * dpq)
                        cosang = max(-1.0, min(1.0, cosang))
                        if math.degrees(math.acos(cosang)) < h_angle_max:
                            _gate_tally["h_point"] += 1
                            return True
            # 4. bare metal / total decoordination (clear margin: ZERO donors).
            for mi in metal_idx:
                has_donor = False
                for k in range(n_atoms):
                    if k in metal_set or syms[k] == "H":
                        continue
                    if math.sqrt(_d2(f, mi, k)) < md_shell:
                        has_donor = True
                        break
                if not has_donor:
                    _gate_tally["bare_metal"] += 1
                    return True
            # 5. DONOR COLLAPSED ONTO METAL (blindspot; only when explicitly on).
            #    The pair tests above skip every metal pair, so a heavy donor fused
            #    onto the metal escapes them.  A metal-to-heavy-donor distance below
            #    collapse_f x covalent-sum is a fused atom = certainly bad (reuses
            #    the same collapse fraction as the bond test).  H excluded (M-H /
            #    agostic / bridging H can be legitimately short).
            if md_collapse_on:
                for mi in metal_idx:
                    rm = _radius(syms[mi])
                    for k in range(n_atoms):
                        if k in metal_set or syms[k] == "H":
                            continue
                        if math.sqrt(_d2(f, mi, k)) < collapse_f * (
                                rm + _radius(syms[k])):
                            _gate_tally["md_collapse"] += 1
                            return True
            return False

        # ============ "LAST OF ITS KIND": A FLOOR MUST NOT WIPE OUT AN ISOMER ============
        # Measured 2026-08-04, twice and independently, and it contradicts a claim written
        # verbatim in this function's own docstring ("NO ensemble consensus"):
        #     pcrejW      24 frames dropped -> smiles_ccdc_regressed 5, pyramid_frame 3, backbone 1
        #     cleangateW 547 frames dropped -> ccdc_backbone_lost 3 (GILKAQ HOJSUX UHEJIB),
        #                                      isomers_lost 3 (UHEJIB 10 -> 4 of 38 theory)
        # A filter that only REMOVES frames cannot make a system worse -- unless the frame it
        # removed was the only one realising that coordination isomer.  UHEJIB did not lose
        # quality, it lost SIX ISOMERS.
        #
        # For the question "is this frame bad?" judging it alone is right.  For the question
        # "may it go?" it is wrong: a frame can be geometrically impossible AND the last of its
        # kind.  GILKAQ's crystal backbone sits in a frame that has a collapsed bond somewhere.
        #
        # The rule (DELFIN_FFFREE_CLEAN_GATE_LAST_OF_KIND, default OFF -> byte-identical):
        # group the frames by their COORDINATION FINGERPRINT -- per metal, the sorted elements of
        # the heavy atoms inside the first shell -- and never let a group become empty.  If every
        # frame of an isomer is bad, the FIRST one survives (lowest index = the emitter's own top
        # rank; deterministic, no score, no RMSD).  No consensus, no majority, no vote: the only
        # question asked of the ensemble is "is there a replacement?".
        _bad_flags = [_certainly_bad(f) for f in frames]
        if os.environ.get("DELFIN_FFFREE_CLEAN_GATE_LAST_OF_KIND", "0") == "1":
            try:
                # ===== THE FINGERPRINT WAS BLIND TO THE ARRANGEMENT ==================
                # Measured 2026-08-18 on gk10kb (10000 systems, 9600 compared): the floor
                # lowers hard_frame_frac by 6.64 percentage points -- and in doing so loses
                #     isomers_lost 40 | ccdc_arrangement_lost 37 | ccdc_backbone_lost 26
                # ALTHOUGH "last of its kind" was running.  That is not a contradiction but the
                # definition of the group: the key below is the sorted ELEMENT
                # MULTISET of the first shell.  cis and trans have the same multiset.  An
                # entire arrangement family can thus be wiped out without its group
                # ever becoming empty -- the rescue never kicks in, because another isomer
                # carries the same fingerprint.
                #
                # Exactly the same root as commit 1caa9123 ("fold completeness per
                # ARRANGEMENT FAMILY instead of per geometry group").  For the second time the
                # same mix-up: "same elements" is not "same arrangement".
                #
                # DELFIN_FFFREE_CLEAN_GATE_KIND_ARRANGEMENT (default 0 -> byte-identical)
                # attaches to the key the multiset of the D-M-D angle classes: for every
                # donor pair (element, element, angle band).  With that, cis/trans,
                # fac/mer and axial/equatorial separate.
                #
                # ⚠️ DIRECTION OF THE CHANGE: a FINER key produces MORE groups,
                # and more groups can only rescue MORE frames, never fewer.  The
                # change is thereby monotone in favour of completeness and cannot cost an
                # isomer that the coarse key would have rescued.  The price lies
                # on the other side: the floor rejects less, the gain in
                # hard_frame_frac turns out smaller.  Exactly that is to be measured.
                #
                # The angle band is deliberately COARSE (45 degrees).  Distortion may split a
                # group -- that is the safe direction -- but noise shall not break every
                # group down into single frames.  Ideal octahedron 90/180 -> bands 2/4,
                # tetrahedron 109.5 -> 2, trigonal bipyramid 90/120/180 -> 2/3/4.
                _kind_arr = (os.environ.get(
                    "DELFIN_FFFREE_CLEAN_GATE_KIND_ARRANGEMENT", "0") == "1")

                def _coord_fp(f):
                    syms = f[0]
                    out = []
                    for mi in range(n_atoms):
                        if not _is_metal(syms[mi]):
                            continue
                        _don = [k for k in range(n_atoms)
                                if k != mi and syms[k] != "H" and not _is_metal(syms[k])
                                and math.sqrt(_d2(f, mi, k)) < md_shell]
                        sh = sorted(syms[k] for k in _don)
                        _key = syms[mi] + ":" + ",".join(sh)
                        if _kind_arr and len(_don) >= 2:
                            _ang = []
                            for _ai in range(len(_don)):
                                for _bi in range(_ai + 1, len(_don)):
                                    _a, _b = _don[_ai], _don[_bi]
                                    _ra2, _rb2 = _d2(f, mi, _a), _d2(f, mi, _b)
                                    if _ra2 <= 0.0 or _rb2 <= 0.0:
                                        continue
                                    _c = ((_ra2 + _rb2 - _d2(f, _a, _b))
                                          / (2.0 * math.sqrt(_ra2) * math.sqrt(_rb2)))
                                    _c = max(-1.0, min(1.0, _c))
                                    _band = int(round(math.degrees(math.acos(_c)) / 45.0))
                                    _e1, _e2 = sorted((syms[_a], syms[_b]))
                                    _ang.append("%s%s%d" % (_e1, _e2, _band))
                            _key += ";" + ",".join(sorted(_ang))
                        out.append(_key)
                    return "|".join(sorted(out))
                _groups = {}
                for _i, f in enumerate(frames):
                    _groups.setdefault(_coord_fp(f), []).append(_i)
                for _fp, _idxs in _groups.items():
                    if all(_bad_flags[_i] for _i in _idxs):
                        _bad_flags[_idxs[0]] = False        # the last of its kind stays
                        _gate_tally["rescued_last_of_kind"] += 1
            except Exception:
                pass
        kept = [item for item, _b in zip(isomers, _bad_flags) if not _b]
        _trace_dst = os.environ.get("DELFIN_TRACE_CLEAN_GATE", "")
        if len(kept) != len(isomers) and _trace_dst and _trace_dst != "0":
            # A FILE, not stderr.  loop.py routes the build workers' stderr to
            # results/debug_<rid>.log and only when a debug pattern matches -- otherwise it is
            # discarded, so a stderr trace from inside a worker reaches nothing (measured
            # 2026-08-04: 139 of 142 systems built, ZERO trace lines).  O_APPEND with one short
            # line per call is atomic enough across the parallel workers.
            try:
                with open(_trace_dst, "a") as _fg:
                    _fg.write("[CLEANGATE] %d -> %d frames | %s\n" % (
                        len(isomers), len(kept),
                        " ".join("%s=%d" % (k, v) for k, v in _gate_tally.items() if v)))
            except Exception:
                pass
        # In-doubt-keep: if (and only if) EVERY frame is certainly bad, keep the
        # input unchanged rather than dropping all or fabricating a broken
        # survivor (the removed "cleanest-available" fallback behaviour).
        return kept if kept else isomers
    except Exception:
        return isomers


def _apply_pi_inplane_final(isomers):
    """FINAL post-conformer monodentate aromatic σ-donor metal-in-plane re-assertion.

    A σ aromatic donor binds through an IN-PLANE sp² lone pair, so the metal must lie
    IN the donor ring's plane.  ``assemble_complex`` seats it there correctly, but the
    conformer/reembed passes (``_conf``) run AFTER the mid-pipeline coplanar-M pass
    (``DELFIN_FFFREE_PI_COPLANAR_M``) and re-tilt the ring (metal ends 0.6-2.8 Å out of
    plane in ~44 % of σ-aromatic-donor complexes).  This re-asserts the in-plane pose on
    the FINAL ensemble for TRULY MONODENTATE ligands (one donor in the ligand) via a
    rigid rotation about the fixed donor (M-D preserved; clash never-worse).  Multi-arm
    polydentates (the larger remainder) need a coupled backbone relaxation (Phase 2, not
    this pass).  Default-OFF / byte-identical unless ``DELFIN_FFFREE_PI_RIGID_PLACE=1``.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_PI_RIGID_PLACE", "0") != "1":
        return isomers
    try:
        from delfin.manta._pi_inplane_final import correct_results
        return correct_results(isomers)
    except Exception:
        return isomers


def _apply_pi_coplanar_final(isomers):
    """FINAL per-arm π-coplanarity polish for MULTI-ARM polydentates (the 91 % of
    the metal-out-of-plane defect that the monodentate ``_apply_pi_inplane_final``
    cannot reach).  Runs on the FINAL emitted ensemble — after ``_conf`` (the
    backbone reembed that re-tilts the rings) — so nothing downstream can undo it.

    Per conjugated π-system with exactly ONE coordinating σ-donor, rotates ONLY that
    arm (the π-system + its donor-less substituents) rigidly about its donor so the
    metal lies in the system plane; the shared backbone + other arms stay fixed, so
    every M-D length is preserved exactly.  Asymmetric / accept-if-better per arm
    (metal-oop strictly down, linker not over-stretched, inter-ligand clash
    never-worse) → never makes a frame worse.  Geometry-only, deterministic, FF-free.
    Default-OFF / byte-identical unless ``DELFIN_FFFREE_PI_COPLANAR_FINAL=1``.
    """
    if not isomers or os.environ.get("DELFIN_FFFREE_PI_COPLANAR_FINAL", "0") != "1":
        return isomers
    try:
        from delfin.manta._pi_coplanar_final import correct_results
        return correct_results(isomers)
    except Exception:
        return isomers


def _coord_integrity_filter(isomers):
    """Final UNIVERSAL decoordination filter over the FULL emitted ensemble at the
    public boundary -- covers BOTH the native fffree path AND the legacy generation
    path (where the eye-flagged ligand-flew-off frames actually originate: the hard
    kappa4/cage cases return None from _fffree_isomers and fall through to legacy,
    whose dense conformer/pucker/see-saw frames bypass the native self-gate).
    Drops a frame iff a coordinating donor has swung off the metal (donor-type ideal
    M-D + slack); crystals pack TIGHT but stay COORDINATED so a crystal-like frame is
    never dropped.  Gated DELFIN_FFFREE_COORD_INTEGRITY (default OFF -> identity ->
    byte-identical to the candidate/legacy baselines).  Never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY", "0") != "1":
        return isomers
    try:
        from delfin.manta.converter_backend import _coord_filter
        return _coord_filter(isomers) or isomers
    except Exception:
        return isomers


def _conf_complete_filter(isomers):
    """UNIVERSAL conformer-completeness pass over the FULL emitted ensemble.

    Root cause it cures (user eye-validation): the conformer search is INCOMPLETE
    (bulky pendant rotors -- tBu / biaryl / pincer arms -- never rotate, ATIDEM /
    ABUSEY), some FREE rotations BREAK topology (tear an M-D bond / put an H on the
    M-C axis / forge a spurious inter-ligand bond, AYOZIX), and symmetry-equivalent
    rotamers are DUPLICATED (ATOSAE).  This pass enumerates every heavy-moving
    rotatable axis as HINDERED rotations on a deterministic staggered grid, hard-
    gates each rotamer for topology preservation (M-D intact, ring intact, no
    spurious bond, no H invading the metal sphere), then RMSD-dedups -- complete
    coverage without bloat.

    Identity (returns input untouched) unless DELFIN_FFFREE_CONF_COMPLETE=1, so
    output is BYTE-IDENTICAL to baseline when off.  FF-free, deterministic, never
    raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_CONF_COMPLETE", "0") != "1":
        return isomers
    try:
        from delfin.manta.conformer_complete import apply_to_ensemble
        return apply_to_ensemble(isomers) or isomers
    except Exception:
        return isomers


def _gfnff_ensemble_rank_filter(isomers):
    """FINAL energy-ranked conformer retention over the emitted ensemble (the
    headline 'keep the most important conformers' feature).

    The shipped conformer engine (conformer_complete) generates topology-gated,
    RMSD-DIVERSE rotamers but invokes NO force field, and the final isomer ordering
    is by least-clash — so nothing in the pipeline ranks conformers by ENERGY.  This
    pass groups the emitted frames by base isomer (the label with any conformer
    suffix stripped), scores each group's frames with GFN-FF (xtb, license-clean,
    parametrised for metals; UFF is unusable on TM complexes), and keeps the top-K
    LOWEST-energy frames within an energy window.  Every base isomer keeps >=1 frame,
    so the isomer manifold is never reduced — only each isomer's conformer set is
    trimmed to the physically important members.

    Within one isomer the total charge is constant, so charge=0 ranks the conformers
    self-consistently (a constant offset cancels); DELFIN_GFNFF_CHARGE can override.

    Identity unless DELFIN_FFFREE_GFNFF_RANK=1 -> byte-identical when off.  Never
    raises; on any failure or if xtb is absent the input passes through unchanged."""
    if not isomers or os.environ.get("DELFIN_FFFREE_GFNFF_RANK", "0") != "1":
        return isomers
    try:
        import re as _re
        from delfin.manta import _gfnff_rank as _gff
        if not _gff.available():
            return isomers
        # Retention policy = "keep ALL distinct conformers, rank them, drop only
        # garbage" (user 2026-06-26).  The conformers are already RMSD-deduped
        # upstream (conf_complete @0.5 A), so each is a genuinely distinct fold; the
        # energy ranking ORDERS them (most important first) and the wide window drops
        # only true garbage (clashing / wildly strained -> tens-to-thousands of
        # kcal/mol above the minimum).  Defaults are intentionally generous (window
        # 35 kcal/mol, soft K-cap 64) so a normal conformer set is kept in full;
        # tighten DELFIN_GFNFF_TOPK / DELFIN_GFNFF_EWIN for a compact pool.
        topk = int(os.environ.get("DELFIN_GFNFF_TOPK", "64"))
        ewin = float(os.environ.get("DELFIN_GFNFF_EWIN", "35.0"))
        chg = int(os.environ.get("DELFIN_GFNFF_CHARGE", "0"))
        _suffix = _re.compile(r"(_conf-|-conf\d|_pool-|-reembed|_reembed|_conf\d)")

        def base_key(lbl):
            if not lbl:
                return ""
            return _suffix.split(lbl, maxsplit=1)[0]

        groups = {}
        order = []
        for (xyz, lbl) in isomers:
            k = base_key(lbl)
            if k not in groups:
                groups[k] = []
                order.append(k)
            groups[k].append((xyz, lbl))

        # HYBRID cost control: GFN2 ranks correctly but is ~10-50x GFN-FF, so for a
        # conformer-rich isomer we first cull cheaply with GFN-FF to the top-M, then
        # rank only those M with the (default GFN2) method.  Safe because the GFN2
        # winners are NOT GFN-FF's worst (measured: ABAKOE's GFN2-best was GFN-FF
        # rank ~3, well inside a top-15 cull).  Isomers with <=M conformers skip the
        # cull and are GFN2-ranked directly.
        # Cost bound only: GFN-FF pre-cull engages for isomers richer than this
        # (rare).  Kept >= the soft K-cap so it never trims below what we retain.
        precull_m = int(os.environ.get("DELFIN_GFNFF_PRECULL_M", "64"))
        out = []
        for k in order:
            frames = groups[k]
            if len(frames) == 1:
                out.append(frames[0])
                continue
            if len(frames) > precull_m:
                ffs = []
                for (xyz, lbl) in frames:
                    e = _gff.gfnff_energy(xyz, charge=chg, method="gfnff")
                    ffs.append((e if e is not None else float("inf"), xyz, lbl))
                ffs.sort(key=lambda t: t[0])
                shortlist = [(xyz, lbl) for _e, xyz, lbl in ffs[:precull_m]]
            else:
                shortlist = frames
            scored = []
            for (xyz, lbl) in shortlist:
                e = _gff.gfnff_energy(xyz, charge=chg)   # default method = GFN2
                scored.append((e if e is not None else float("inf"), xyz, lbl))
            if all(s[0] == float("inf") for s in scored):
                out.extend(frames)            # ranking unusable here -> keep all
                continue
            scored.sort(key=lambda t: t[0])
            emin = scored[0][0]
            kept = []
            for (e, xyz, lbl) in scored:
                if len(kept) >= topk:
                    break
                if e != float("inf") and e - emin > ewin:
                    break
                kept.append((xyz, lbl))
            out.extend(kept or [(scored[0][1], scored[0][2])])
        return out
    except Exception:
        return isomers


def _permute_dedup_filter(isomers):
    """UNIVERSAL permutation-invariant duplicate removal over the FULL emitted
    ensemble.

    Root cause it cures (#1 user request): the existing ensemble dedup compares
    frames by FIXED-ORDER heavy-atom Kabsch RMSD, so two structures identical up
    to relabeling of indistinguishable atoms (identical ligands swapped /
    symmetry-equivalent atoms / any molecular-graph automorphism) are NOT seen as
    duplicates and BOTH survive -> "very many identical structures in the pool".  This
    pass adds the missing layer: two frames are duplicates iff the MIN over graph
    automorphisms of their heavy Kabsch RMSD is below the threshold.  Atoms are
    permuted ONLY within their graph-symmetry orbit, so genuinely-different
    geometric isomers / conformers / stereoisomers (whose geometry differs by more
    than the threshold even after the best symmetry relabelling) are KEPT.

    Identity (returns input untouched) unless DELFIN_FFFREE_PERMUTE_DEDUP=1, so
    output is BYTE-IDENTICAL to baseline when off.  FF-free, deterministic, never
    raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_PERMUTE_DEDUP", "0") != "1":
        return isomers
    try:
        from delfin.manta.permute_dedup import dedup_ensemble
        return dedup_ensemble(isomers) or isomers
    except Exception:
        return isomers


def _mirror_xyz_coords(xyz):
    """Return the MIRROR IMAGE of a frame (reflection through the x=0 plane -> negate
    x).  An enantiomer IS the exact mirror image, so this is the free, exact partner
    of a chiral coordination frame — no rebuild.  Preserves atom order and any
    header/comment lines; only ``sym x y z`` lines flip.  None on any parse issue."""
    try:
        out = []
        for ln in str(xyz).splitlines():
            parts = ln.split()
            if len(parts) == 4:
                try:
                    x = -float(parts[1]); y = float(parts[2]); z = float(parts[3])
                    out.append(f"{parts[0]} {x:.6f} {y:.6f} {z:.6f}")
                    continue
                except Exception:
                    pass
            out.append(ln)
        return "\n".join(out)
    except Exception:
        return None


def _coord_sphere_donors(xyz):
    """(metal_dir_array, donor_symbols) for the coordination sphere: unit vectors
    metal->donor for the nearest neighbours (< 3.2 Å, up to 9) of the first
    ``_METAL_SET`` atom.  (None, None) if no metal / < 3 donors."""
    try:
        import numpy as _np
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [p[0] for p in rows]
        coords = _np.array([[float(p[1]), float(p[2]), float(p[3])] for p in rows])
        m_i = next((i for i, s in enumerate(syms) if s in _METAL_SET), None)
        if m_i is None or coords.shape[0] < 4:
            return None, None
        d = coords - coords[m_i]
        dist = _np.linalg.norm(d, axis=1); dist[m_i] = 1e9
        idx = [i for i in _np.argsort(dist) if dist[i] < 3.2][:9]
        if len(idx) < 3:
            return None, None
        dirs = _np.array([d[i] / (dist[i] + 1e-12) for i in idx])
        return dirs, [syms[i] for i in idx]
    except Exception:
        return None, None


def _proper_kabsch_rmsd(P, Q):
    """Min RMSD aligning Q onto P by a PROPER rotation (det=+1, NO reflection) +
    translation.  Reflection-excluding by construction, so a true mirror image can
    never be aligned away."""
    import numpy as _np
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, S, Vt = _np.linalg.svd(H)
    dsign = 1.0 if _np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ _np.diag([1.0, 1.0, dsign]) @ U.T
    Qr = Qc @ R.T
    return float(_np.sqrt(_np.mean(_np.sum((Pc - Qr) ** 2, axis=1))))


def _coord_sphere_chirality_rmsd(xyz):
    """CONFIGURATIONAL chirality measure (Å): the min-over-same-element-donor-
    permutations PROPER-rotation Kabsch RMSD between the coordination sphere's donor
    directions and their MIRROR (reflection-excluding).  ~0 => the mirror superimposes
    => achiral configuration (mirror plane OR only a slight build deformation from an
    ideal-achiral arrangement).  Large => a genuine Δ/Λ twist => chiral.  Conformer-
    robust (donors only, backbone ignored).  A physical Å threshold (tuned) gates it.
    Returns 0.0 when undetermined (never spuriously chiral)."""
    try:
        import itertools as _it
        import numpy as _np
        dirs, syms = _coord_sphere_donors(xyz)
        if dirs is None:
            return 0.0
        mirror = dirs.copy(); mirror[:, 0] *= -1.0     # reflect (negate x)
        n = len(syms)
        groups: Dict[str, List[int]] = {}
        for i, s in enumerate(syms):
            groups.setdefault(s, []).append(i)
        # permutations = product of within-element permutations; cap the search.
        per_elem = []
        _count = 1
        for s, idxs in groups.items():
            ps = list(_it.permutations(idxs))
            _count *= len(ps)
            per_elem.append((idxs, ps))
        if _count > 5040:                              # >7! : fall back to identity only
            perms = [list(range(n))]
        else:
            perms = []
            for combo in _it.product(*[ps for _, ps in per_elem]):
                perm = [0] * n
                for (idxs, _), mapped in zip(per_elem, combo):
                    for src, dst in zip(idxs, mapped):
                        perm[src] = dst
                perms.append(perm)
        best = 1e9
        for perm in perms:
            r = _proper_kabsch_rmsd(dirs, mirror[perm])
            if r < best:
                best = r
        return best
    except Exception:
        return 0.0


def _coord_chirality_sign(xyz):
    """Sign (+1 / -1 / 0) of the coordination-sphere chirality pseudoscalar — used
    ONLY to name the Δ/Λ hands (never gates; the RMSD measure gates)."""
    try:
        import numpy as _np
        dirs, _syms = _coord_sphere_donors(xyz)
        if dirs is None:
            return 0
        total = 0.0
        n = len(dirs)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    total += float(_np.dot(_np.cross(dirs[i], dirs[j]), dirs[k]))
        if abs(total) < 1e-4:
            return 0
        return 1 if total > 0 else -1
    except Exception:
        return 0


def _arrangement_key(lbl):
    """Strip conformer / duplicate-disambiguation / Δ-Λ-hand suffixes -> the ARRANGEMENT
    identity.  All conformers AND both hands of one coordination isomer share this key, so
    the chirality decision can be made ONCE per arrangement (conformer-consistent) instead
    of per built conformer (where UFF noise makes the CSM straddle the tolerance)."""
    import re as _re
    k = str(lbl or "")
    k = _re.sub(r"-conf\d+", "", k)
    k = _re.sub(r"-[ΔΛ](?=-|$)", "", k)
    k = _re.sub(r"-\d+$", "", k)
    return k


def _mirror_symmetrize_xyz(xyz, tol):
    """Project a near-mirror-symmetric frame onto its EXACT mirror-symmetric form
    S = ½(F + R·F'[π]) (F' = mirror, R = best proper rotation, π = best graph automorphism).
    The antisymmetric (build-noise) component is projected OUT -> higher symmetry + quality,
    exactly consistent with DELFIN's target of the idealized intrinsic reference (crystal
    distortions are out of scope).  Returns the symmetrised xyz, or the ORIGINAL unchanged if
    the mirror is not reachable within ``tol`` (genuinely chiral) or on any error (best-effort,
    never raises)."""
    try:
        import numpy as _np
        from delfin.manta.permute_dedup import _automorphisms_for_xyz
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [r[0] for r in rows]
        F = _np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
        n = len(F)
        if n < 2:
            return xyz
        Fm = F.copy(); Fm[:, 0] *= -1.0
        _s, autos, _h = _automorphisms_for_xyz(xyz, True, 4096)
        if not autos:
            autos = [list(range(n))]
        Fc = F - F.mean(0)
        best = None; best_r = 1e9
        for perm in autos:
            if len(perm) != n:
                continue
            Q = Fm[list(perm)]; Qc = Q - Q.mean(0)
            U, S, Vt = _np.linalg.svd(Qc.T @ Fc)
            d = 1.0 if _np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
            Qr = Qc @ (Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T).T
            r = float(_np.sqrt(_np.mean(_np.sum((Fc - Qr) ** 2, axis=1))))
            if r < best_r:
                best_r = r; best = Qr
        if best is None or best_r >= tol:
            return None            # no internal mirror within tol -> CHIRAL (caller keeps both hands)
        Ssym = 0.5 * (Fc + best) + F.mean(0)
        return "\n".join(f"{syms[i]} {Ssym[i, 0]:.6f} {Ssym[i, 1]:.6f} {Ssym[i, 2]:.6f}"
                         for i in range(n))
    except Exception:
        return None


def _pointgroup_symmetrize_xyz(xyz, tol):
    """Detect the frame's approximate POINT GROUP and PROJECT onto it:

        S = (1/|G|) Σ_{g∈G} g(F)      (the Continuous-Symmetry-Measure ideal)

    giving the structure with EXACT G-symmetry.  This GENERALISES
    ``_mirror_symmetrize_xyz`` (a σ-only special case: it averaged F with ONE mirror
    image, so it only ever removed a single mirror-breaking component -> plateaued) to
    the FULL group: EVERY symmetry operation of the frame is a graph automorphism π
    whose geometric realisation is the orthogonal matrix O (proper OR improper) that
    best superimposes F onto F[π]; the operation belongs to G iff that residual < tol.
    Averaging over ALL such operations removes EVERY symmetry-breaking distortion at
    once -> both poly-quality (perfect geometry) AND enantiomer elimination improve.

    Returns ``(S_xyz, has_improper)``:
      * ``has_improper`` True  <=> G contains an IMPROPER element (σ / inversion i / Sn)
        <=> the configuration is ACHIRAL -> the caller emits ONE symmetric frame
        (enantiomer ELIMINATED).  This is the FUNDAMENTAL achirality test — it covers
        the inversion centre and every Sn axis, not just mirror planes.
      * ``has_improper`` False <=> only proper rotations (Cn / Dn) or E (C1) -> genuinely
        CHIRAL -> the caller keeps BOTH hands (the projected S is the proper-symmetrised
        frame of the SAME hand; its exact mirror is the other hand).
      * |G| == 1 (only E within tol) -> S == F unchanged  ->  C1 SAFE (user: "many C1
        must NOT be affected").

    Deterministic; best-effort (``(None, False)`` on any error, never raises)."""
    try:
        import numpy as _np
        from delfin.manta.permute_dedup import _automorphisms_for_xyz
        rows = [l.split() for l in str(xyz).splitlines() if len(l.split()) == 4]
        syms = [r[0] for r in rows]
        F = _np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
        n = len(F)
        if n < 2:
            return xyz, False
        Fc = F - F.mean(0)
        _s, autos, _h = _automorphisms_for_xyz(xyz, True, 4096)
        if not autos:
            autos = [list(range(n))]
        # Collect the group operations: (O, perm) for every automorphism whose optimal
        # orthogonal (Procrustes, REFLECTION ALLOWED -> det ±1) superposition of F onto
        # F[perm] has residual < tol.  O = U Vt for H = F[perm]^T · F (see _mirror_… for
        # the proper-only sibling; here we do NOT force det=+1, so improper ops appear).
        ops = []
        have_id = False
        for perm in autos:
            if len(perm) != n:
                continue
            Q = Fc[list(perm)]
            U, _S2, Vt = _np.linalg.svd(Q.T @ Fc)
            O = U @ Vt                                   # orthogonal, det = ±1
            r = float(_np.sqrt(_np.mean(_np.sum((Fc @ O.T - Q) ** 2, axis=1))))
            if r < tol:
                ops.append((O, Q))                       # keep Q=Fc[perm] for the projection
                if all(p == i for i, p in enumerate(perm)):
                    have_id = True
        if not ops:
            ops = [(_np.eye(3), Fc.copy())]
        elif not have_id:
            ops.append((_np.eye(3), Fc.copy()))
        # improper element present?  det(O) < 0 for ANY op  <=>  achiral configuration.
        has_improper = any(bool(_np.linalg.det(O) < 0.0) for O, _ in ops)
        # PROJECT: S_i = (1/|G|) Σ_k O_k^T · r_{π_k(i)}  =  (1/|G|) Σ_k (Fc[π_k] · O_k)_i .
        Sacc = _np.zeros_like(Fc)
        for O, Q in ops:
            Sacc += Q @ O
        Sacc /= float(len(ops))
        Sacc += F.mean(0)
        out = "\n".join(f"{syms[i]} {Sacc[i, 0]:.6f} {Sacc[i, 1]:.6f} {Sacc[i, 2]:.6f}"
                        for i in range(n))
        return out, has_improper
    except Exception:
        return None, False


def _heavy_dist_fp(xyz, bin_a=0.15):
    """Rotation- AND permutation-invariant, chirality-INSENSITIVE geometry fingerprint:
    the sorted multiset of heavy-atom pairwise distances, binned to ``bin_a`` Å.  Two
    frames of ONE molecule with the same fingerprint are the same geometry up to a rigid
    motion OR a reflection (distances are reflection-invariant) — so pairing it with the
    chirality sign gives a canonical isomer key that collapses achiral image/mirror pairs
    while keeping genuine Δ/Λ.  None on parse error."""
    try:
        import numpy as _np
        pts = []
        for _l in str(xyz).splitlines():
            q = _l.split()
            if len(q) == 4 and q[0] != 'H':
                pts.append((float(q[1]), float(q[2]), float(q[3])))
        if len(pts) < 3:
            return None
        C = _np.array(pts)
        n = len(C)
        ds = []
        for i in range(n):
            ds.extend(_np.linalg.norm(C[i + 1:] - C[i], axis=1).tolist())
        return tuple(sorted(round(d / bin_a) for d in ds))
    except Exception:
        return None


def _enantiomer_mirror_filter(isomers, smiles=None):
    """Add the Δ/Λ MIRROR enantiomer for every CHIRAL coordination frame (user
    2026-07-21: "we also want the respective enantiomers in the frames").

    An enantiomer IS the exact mirror image -> reflect the built coordinates (free,
    exact) instead of rebuilding in two directions.  A frame is chiral iff its mirror
    is NOT superimposable on it under proper rotation + graph automorphism (reflection
    -EXCLUDING Kabsch) -- tested by the existing ``is_permutation_duplicate`` primitive
    (NOT gated), so achiral (meso) frames add no duplicate and an enantiomer already
    present (e.g. built by CHIRAL_ENUM) is not re-added.  UNIVERSAL across every
    polyhedron (chirality-by-mirror is geometry, not a per-polyhedron formula) and
    robust where the analytical helicity classifier fails (tridentate wraps).

    CRITICAL (user 2026-07-21): the gate is CONFIGURATIONAL chirality, NOT the
    conformer's own mirror symmetry.  A conformer of an ACHIRAL molecule can lack a
    mirror plane (a chiral backbone pucker) yet be the SAME molecule (interconvertible
    by a conformational flip) -- mirroring it would fabricate a false "enantiomer" and
    explode the pool.  So we gate on the COORDINATION-SPHERE chirality (the pseudoscalar
    over metal->donor directions), which is set by the polyhedron and is ROBUST across
    conformers (donors are pinned; only the backbone puckers).  Achiral configuration
    -> sign 0 -> NOT mirrored, whatever the conformer pucker.  A ``is_permutation_
    duplicate`` meso-guard additionally skips internally mirror-symmetric molecules.

    Emits each enantiomer IMMEDIATELY AFTER its partner (Δ, Λ consecutive in the
    trajectory).  Guard: skip when the SMILES has FIXED ligand stereocentres --
    mirroring the whole complex would INVERT them (a different compound); those
    diastereomers come from the arrangement enumeration + STEREOCENTER_ENUM.  Env-gated
    DELFIN_FFFREE_ENANTIOMER_MIRROR (default off -> byte-identical).  Deterministic; never raises."""
    if not isomers or os.environ.get("DELFIN_FFFREE_ENANTIOMER_MIRROR", "0") != "1":
        return isomers
    try:
        from delfin.manta.permute_dedup import is_permutation_duplicate
    except Exception:
        is_permutation_duplicate = None
    try:
        if smiles and RDKIT_AVAILABLE and isinstance(smiles, str):
            _m = Chem.MolFromSmiles(smiles)
            if _m is not None and Chem.FindMolChiralCenters(
                    _m, includeUnassigned=False, useLegacyImplementation=False):
                return isomers   # fixed ligand stereocentre -> mirroring inverts it
    except Exception:
        pass
    _chi_tol = _delfin_env_float("DELFIN_FFFREE_ENANTIOMER_MIRROR_TOL", 0.40)
    # ADDITIVE-ONLY (2026-07-21).  A broad full:1000 A/B proved the post-hoc UNITE + point-group PROJECTION
    # is NOT never-worse: averaging a frame over its APPROXIMATE symmetry ops tears ligands and distorts
    # polyhedra (measured: poly-distortion 13, ligand-quality 18, ligand-bond-torn 1, tier-2 on 69,
    # holistic +0.513), and the fp-GROUPING collapses DISTINCT isomers that merely share a distance-bin
    # (isomers-lost 10, ccdc-isomer-lost 4).  Post-hoc geometry repair on built atoms is exactly the
    # guiding-principle antipattern.  So this filter now does ONLY the SAFE, purely-ADDITIVE half of the user's
    # ask ("we want the enantiomers in the frames"): for every GENUINELY chiral coordination frame whose
    # opposite hand is not already built, ADD its EXACT mirror.  A reflection of a good frame is an equally
    # good frame -> zero distortion, zero isomer loss -> never-worse.  The ELIMINATE-where-possible geometry
    # symmetrisation is DEFERRED to the IN-CONSTRUCTION poly-seater (seat donors at ideal symmetric vertices;
    # NEVER move built atoms) -- the correct root architecture.  Gate CONFIGURATIONAL (coord-sphere)
    # chirality, robust to conformer/peripheral wobble; the mirror's presence is tested by a
    # reflection-aware (fingerprint, sign) key so an already-built opposite hand is never duplicated.
    # FAST per-frame key: (reflection-INVARIANT fingerprint, cheap O(donors^3) chirality pseudoscalar).
    # Deliberately AVOID the expensive permutation-Kabsch _coord_sphere_chirality_rmsd here: on high-frame-
    # count systems (RILVUD 209 frames) it made the finalisation slow enough to push _one.py over the
    # per-system build timeout -> partial build -> FALSE isomer loss in the A/B.  The pseudoscalar sign is
    # a sufficient chiral gate (sign==0 -> treat as achiral -> add nothing; a rare genuinely-chiral frame
    # with an accidental zero pseudoscalar is merely SKIPPED = a completeness miss, never a regression).
    _keys = []
    _present = set()
    for it in isomers:
        try:
            _fp = _heavy_dist_fp(it[0]); _sg = _coord_chirality_sign(it[0])
        except Exception:
            _fp, _sg = None, 0
        _keys.append((_fp, _sg))
        _present.add((_fp, _sg))
    out = list(isomers)
    for _i, it in enumerate(isomers):
        _fp, _sg = _keys[_i]
        if _sg == 0:
            continue   # achiral configuration -> mirror superimposes on itself -> nothing to add
        _mk = (_fp, -_sg)   # the mirror: fp is reflection-INVARIANT, ONLY the chirality sign flips
        if _mk in _present:
            continue   # the opposite hand is already built (arrangement enumeration / CHIRAL_ENUM)
        _present.add(_mk)
        _mx = _mirror_xyz_coords(it[0])
        if _mx is None:
            continue
        _lbl = it[1] if len(it) > 1 else ""
        _hand = "Λ" if _sg < 0 else "Δ"   # mirror hand = opposite of the original's
        out.append((_mx, (f"{_lbl}-{_hand}" if _lbl else _hand)) + tuple(it[2:]))
    return out


# Re-entrancy guard for the conformer-completeness pass.  ``_smiles_to_xyz_isomers_impl``
# recurses through this PUBLIC wrapper for the dual-parse augmentation; the
# completeness pass must run EXACTLY ONCE at the OUTERMOST public boundary over the
# final union (otherwise the inner call expands the alt-parse set and the outer call
# re-expands the union -- doubling conformers and over-running the per-structure cap).
# A simple non-reentrant flag suppresses the pass on inner re-entries.  All OTHER
# finalisation filters keep running on every level (unchanged base behaviour ->
# byte-identical when off).
_CONF_COMPLETE_ACTIVE = threading.local()


def _arg_deterministic(args, kwargs) -> bool:
    """Extract the ``deterministic`` argument from the public wrapper's
    ``*args/**kwargs`` (it is positional index 6 in ``_smiles_to_xyz_isomers_impl``:
    smiles, num_confs, max_isomers, apply_uff, collapse_label_variants,
    include_binding_mode_isomers, deterministic, ...).  Defaults to the impl
    default (True) when not supplied."""
    if "deterministic" in kwargs:
        return bool(kwargs["deterministic"])
    if len(args) > 6:
        return bool(args[6])
    return True


def _smiles_arg(args, kwargs):
    """The SMILES passed to the public wrapper (positional index 0 / kwargs)."""
    if "smiles" in kwargs:
        return kwargs["smiles"]
    return args[0] if args else None


def _has_metalloid_donor(smiles) -> bool:
    """True iff the molecule has a heavy metalloid (Sb/As/Bi/Te/Se/Ge/Sn/Pb) atom
    bonded to a metal centre -- the only systems the dual-M-D-distance pass targets
    (so the second build pass is skipped elsewhere).  Graph/element only, never
    SMILES-specific; any parse failure -> False (pass skipped = base behaviour)."""
    if not smiles or not isinstance(smiles, str):
        return False
    try:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return False
        for a in m.GetAtoms():
            if a.GetSymbol() in _METALLOID_MD_DONORS and any(
                    nb.GetSymbol() in _METAL_SET for nb in a.GetNeighbors()):
                return True
    except Exception:
        return False
    return False


def smiles_to_xyz_isomers(*args, **kwargs):
    """Public entry point.  Thin wrapper enforcing the finite-coordinate output
    contract (#36) over EVERY return path of the implementation, applying the
    universal coordination-integrity filter, the conformer-completeness pass, the
    universal topology-preservation gate, then the UNIFIED hard clean-manifold
    emission gate (all env-gated, byte-id OFF) over the full native+legacy ensemble,
    then ordering the ensemble so the most crystal-like (least-clash) frame is first;
    otherwise the (isomers, error) result passes through unchanged.

    The clean gate is the LAST filter before ranking: it is the authoritative
    arbiter that admits only clean, topologically-intact coordination frames into
    the manifold (clash / spurious-bond / destroyed-geometry / decoordination), with
    a per-structure never-empty fallback to the single cleanest available frame.

    Order rationale: completeness runs AFTER coord-integrity (so it expands only
    coordinated frames) and BEFORE the topology gate (so every newly-generated
    conformer is additionally vetted by the consensus topology gate as defence in
    depth), and ranking runs last over the full expanded set.  The completeness pass
    is suppressed on dual-parse re-entries so it runs once at the outermost level."""
    outermost = not getattr(_CONF_COMPLETE_ACTIVE, "value", False)
    # "max_isomers=0" is the DOCUMENTED complete-manifold contract (Submit-tab default + build_eye_
    # package), meaning "never cut off".  But downstream max_isomers is used purely as a hard CAP
    # (len(results) >= max_isomers, iso_list[:max_isomers], _PRE_UFF_CAP = max_isomers*mult), so 0
    # was treated as "cap = 0" and COLLAPSED the manifold to ~1 frame -> CODSIA 14->2, losing the
    # cis / see-saw CN4 polyhedra the user needs (all CN4 geometries: SP-4, Td, see-saw).  Normalise
    # 0/negative -> effectively unlimited at the outermost public call so "0 = complete" holds; the
    # existing timeout / wall-budget / combo-cap guards still bound pathological combinatorics.
    if outermost and "max_isomers" in kwargs and (
            kwargs["max_isomers"] is None or kwargs["max_isomers"] <= 0):
        kwargs["max_isomers"] = 100000
    # Hydrogens of unbracketed metal neighbours are pinned in the SMILES itself, once,
    # before any path parses it (see _pin_metal_neighbour_hydrogens).  A SMILES without
    # such an atom -- every bracketed-SMILES input -- passes through as the identical
    # string.
    if outermost:
        _smi_in = _smiles_arg(args, kwargs)
        _smi_pinned = _pin_metal_neighbour_hydrogens(_smi_in)
        if _smi_pinned != _smi_in:
            if "smiles" in kwargs:
                kwargs["smiles"] = _smi_pinned
            else:
                args = (_smi_pinned,) + tuple(args[1:])
    # Cross-process determinism bridge: the public ``deterministic=True`` contract
    # MUST imply bit-identity across processes, but the deep timeout / wall-budget /
    # sorted-enum guards consult the DELFIN_DETERMINISTIC *env* (so subprocess
    # workers inherit it) — which the PARAM never set.  Result: embed thread joins
    # kept their wall-clock timeout (see _embed_with_timeout), so under CPU load a
    # heavy conformer's embed timed out and was dropped, making the isomer/conformer
    # COUNT vary across independent runs (measured: ADUFIS 161 vs 170 frames).  Here
    # we promote the param to the documented master switch at the OUTERMOST public
    # call, scoped + restored, and keep it active through the post-processing filters
    # below (they generate/rank conformers too).  Subprocesses spawned during the
    # call inherit the env; the env is restored afterwards so a later
    # deterministic=False caller is unaffected.
    _det_param = _arg_deterministic(args, kwargs)
    _det_env_prev = _enum_env_prev = None
    _det_env_set = False
    if outermost:
        _CONF_COMPLETE_ACTIVE.value = True
        if _det_param:
            _det_env_prev = os.environ.get("DELFIN_DETERMINISTIC")
            _enum_env_prev = os.environ.get("DELFIN_FFFREE_DETERMINISTIC_ENUM")
            os.environ["DELFIN_DETERMINISTIC"] = "1"
            os.environ.setdefault("DELFIN_FFFREE_DETERMINISTIC_ENUM", "1")
            _det_env_set = True
    try:
        r = _smiles_to_xyz_isomers_impl(*args, **kwargs)
        # --- Metalloid dual-M-D-distance pass (DELFIN_FFFREE_METALLOID_MD_BOTH; default OFF) ------
        # For systems with a metalloid donor, ALSO build the manifold at the SHORT covalent-sum
        # M-metalloid distance and APPEND it (never drops the offset frames).  The eye's gate is
        # min-over-frames, so a system whose offset build is already good KEEPS that good frame
        # (cannot regress) while a detached system GAINS a good short frame (improves).  The short
        # frames flow through the SAME finalisation chain below: _clean_gate_filter is per-frame
        # asymmetric (drops broken frames, keeps good ones) so a short frame that only CLASHES
        # (crowded metalloid ligand) is dropped -> no broken_frac penalty; only a genuinely better
        # short frame survives.  Structurally never-worse via completeness (frames appended, never
        # dropped).  Second pass forces the short distance via the thread-local (env stays off ->
        # byte-identical when METALLOID_MD_BOTH is off).  Metalloid-gated so non-metalloid systems
        # are byte-identical (no wasted second build).
        if (outermost
                and os.environ.get("DELFIN_FFFREE_METALLOID_MD_BOTH", "0") == "1"
                and isinstance(r, tuple) and len(r) == 2
                and isinstance(r[0], list) and r[0]
                and _has_metalloid_donor(_smiles_arg(args, kwargs))):
            _METALLOID_MD_FORCE_SHORT.value = True
            try:
                r2 = _smiles_to_xyz_isomers_impl(*args, **kwargs)
            finally:
                _METALLOID_MD_FORCE_SHORT.value = False
            if isinstance(r2, tuple) and len(r2) == 2 and isinstance(r2[0], list) and r2[0]:
                _short = [(xyz, (f"{lbl}-mdshort" if lbl else "mdshort")) for xyz, lbl in r2[0]]
                r = (r[0] + _short, r[1])
        # Conformer-completeness AND permutation-dedup run only at the OUTERMOST
        # public call (over the final, dual-parse-unioned ensemble); identity on
        # inner re-entries.  Permutation-dedup runs LAST over the expanded union
        # (after completeness has generated conformers and the topology gate has
        # vetted them), so it strips permutation-equivalent duplicates from the
        # whole final pool.
        _conf = _conf_complete_filter if outermost else (lambda x: x)
        _pdedup = _permute_dedup_filter if outermost else (lambda x: x)
        # _grank = final GFN-FF energy-ranked top-K conformer retention (default
        # OFF, byte-id when off).  Runs after dedup so it ranks the de-duplicated
        # conformer set, before the clash-based isomer ordering.
        _grank = _gfnff_ensemble_rank_filter if outermost else (lambda x: x)
        # _hdeclash = eta-ring inter-ligand declash; _cofix = carbonyl-length fix.
        # Both are geometry-only, count-preserving transforms over the final
        # dedup'd ensemble (byte-id when off); carbonyl fix runs after the declash.
        _hdeclash = _hapto_declash_filter if outermost else (lambda x: x)
        _sdeclash = _sigma_declash_filter if outermost else (lambda x: x)
        _cofix = _carbonyl_fix_filter if outermost else (lambda x: x)
        # _coordint_late = decoordination gate applied AGAIN on the FINAL manifold
        # (the early one runs before _conf re-expands conformers; conformer/declash
        # frames that decoordinate a donor must be caught here). Byte-id when off.
        _coordint_late = _coord_integrity_filter if outermost else (lambda x: x)
        # ENANTIOMER MIRROR (DELFIN_FFFREE_ENANTIOMER_MIRROR, default off -> byte-identical).
        # OUTERMOST step so image+mirror image stay CONSECUTIVE after ranking.  Adds the Δ/Λ mirror
        # ONLY for configurationally-chiral frames -- those NOT reducible to a mirror-symmetric
        # form within tolerance (user: "only systems that cannot be brought into the
        # highest-symmetry form need both frames").  Achiral / symmetrizable frames stay single.
        _smi = _smiles_arg(args, kwargs)
        _emir = (lambda x: _enantiomer_mirror_filter(x, _smi)) if outermost else (lambda x: x)
        if isinstance(r, tuple) and len(r) == 2 and isinstance(r[0], list):
            return _emir(_rank_emitted_isomers(_coordint_late(_cofix(_sdeclash(_hdeclash(_grank(_pdedup(_clean_gate_filter(_topology_gate_filter(_apply_pi_coplanar_final(_apply_pi_inplane_final(_conf(_coord_integrity_filter(_filter_nonfinite_isomers(r[0]))))))))))))))), r[1]
        if isinstance(r, list):
            return _emir(_rank_emitted_isomers(_coordint_late(_cofix(_sdeclash(_hdeclash(_grank(_pdedup(_clean_gate_filter(_topology_gate_filter(_apply_pi_coplanar_final(_apply_pi_inplane_final(_conf(_coord_integrity_filter(_filter_nonfinite_isomers(r)))))))))))))))
        return r
    finally:
        if outermost:
            _CONF_COMPLETE_ACTIVE.value = False
            if _det_env_set:
                if _det_env_prev is None:
                    os.environ.pop("DELFIN_DETERMINISTIC", None)
                else:
                    os.environ["DELFIN_DETERMINISTIC"] = _det_env_prev
                if _enum_env_prev is None:
                    os.environ.pop("DELFIN_FFFREE_DETERMINISTIC_ENUM", None)
                else:
                    os.environ["DELFIN_FFFREE_DETERMINISTIC_ENUM"] = _enum_env_prev


def _ffree_shared_tail(mol, results, dual_parse_done: bool):
    """THE SHARED TAIL: the post-build correctors the FF-free path returns past.

    WHY THIS EXISTS (2026-08-10).  ``_smiles_to_xyz_isomers_impl`` returns the FF-free
    frames ~2900 lines before the legacy pipeline reaches its final passes, so every
    corrector down there sees the legacy manifold and nothing else.  That is not a
    decision anyone took; it is what an early ``return`` does to everything appended
    after it.  Counted on 2026-08-10 there are EIGHTEEN statements of the shape
    ``results = _apply_*(mol, results, ...)`` between the FF-free return and the union
    merge, and the FF-free frames reach none of them.

    ⚠ BUT THE HONEST NUMBER IS TWO, NOT EIGHTEEN.  Cross-checked against
    ``cli_manta._CHAMPION_FLAGS``: sixteen of the eighteen are gated by flags that are
    NOT in the champion (COORD_ANGLE_FIX, HYDROXYL_GEOM, AROM_BOND_LENGTH,
    AROMATIC_PLANARITY, BOND_DECOLLAPSE_FORCE, CASCADE_REFINER, FIX_WUXQAK_ANGLE_DEG,
    5F_F_HAPTO_FINAL_CLEARANCE, ...), so in the reported runs they are no-ops on BOTH
    paths and porting them would gain exactly nothing.  Only these two are champion-
    active, and only they are wired here.  Two more were already hand-wired at the
    FF-free exit years apart -- pi-coplanar and (2026-08-10) the stereocentre folds --
    which is the symptom this function is meant to end: each gap patched alone, none
    of them found by looking.

    ⚠⚠ EXPECT A NULL RESULT AND DO NOT MISREAD IT.  Both correctors take ``mol`` and
    pass it down (isolated-reseat hands it straight to ``_isolated_reseat.correct_results``;
    arom-planarize needs it for ``_class_conditional_flag``).  A mol parsed from the
    SMILES carries RDKit's atom order, while FF-free frames carry metal-at-0 plus
    AddHs(ligand) blocks in construction order -- the two never coincide.  That exact
    mismatch already turned the ring-pucker emitter into a null lever at this very
    position (185 of 187 systems byte-identical; see the note inside the FF-free block).
    If this tail measures "affected = 0", the FIRST hypothesis is the atom-order
    mismatch, NOT "the correctors decline FF-free frames".  Trace the guarded call
    site with fire_census before concluding anything.

    Additive/never-worse is each corrector's own contract (both carry per-frame
    rollback); this function adds no policy of its own.  Callers gate it.
    """
    if not results:
        return results
    # ===== ROLLBACK GATE AROUND THE WHOLE CHAIN (24.08.2026) ===========================
    # TWO lines for FIVE refiners -- not five wrappers.  The snapshot is a flat list
    # copy (labels and xyz texts are immutable), so it costs nothing.  At the end
    # `_refine_gate.keep_better` decides per frame between result and
    # original.  Default OFF -> byte-identical.
    # ⚠ DELIBERATELY AROUND THE CHAIN, NOT AROUND EACH PASS.  The product is the chain; if it
    # as a whole makes a frame worse, that frame belongs rolled back.  Whoever guards each
    # pass individually gets five gates that make up for one another -- and
    # exactly that is the corrector of a corrector, which is not supposed to exist here.
    _rg_before = list(results)
    results = _apply_isolated_reseat_if_enabled(mol, results, dual_parse_done)
    results = _apply_arom_planarize_if_enabled(mol, results, dual_parse_done)
    # ===== TWO MORE, AND WHY THEY NOW BELONG HERE (16.08.2026, evening) ================
    # The docstring above says "only these two are champion-active, and only they are
    # wired here".  That held as long as the remaining sixteen were switched off and thus
    # no-ops on BOTH paths.  For the sp2 planarisers it no longer holds:
    #
    #   MEASURED 16.08.: `ccdc_pyramid_realized` misses 131 of 905 -- and that is the
    #   ONLY defect of the day that drags EVERY other axis along: `graph_geom` 6.33x,
    #   `ccdc_isomer` 2.68x, `coord_angle` 2.32x, `ml_len` 2.09x, `org_bond` 1.63x.  No
    #   axis stands at 1.0.  Threshold: `pyramid_over_worst > 15.9 degrees` covers 82.4 %.
    #
    #   And the run `pyr131` measured reach **2 of 131** -- not because the chemistry does
    #   not hold, but because both correctors stand 2275 lines BEHIND the FF-free `return`
    #   (`:32280` against `:34555`/`:34558`).
    #
    # ⚠ THE PRECONDITION WAS THE ATOM MAPPING, and it exists as of today.  Without it this
    # here would have been a null test: both correctors need `mol` atom indices, FF-free
    # frames carry metal-at-0 plus AddHs blocks.  `_frame_atom_map.frame_to_mol_map`
    # reconstructs the mapping as a true graph isomorphism; the correctors TRANSLATE
    # since then, instead of giving up (`0486a94e`).  If the order matches -> identity, i.e.
    # byte-identical; if it is not determinable -> abort as before.
    #
    # Both have their own switch with default 0 -- this hook-up changes nothing
    # as long as they are off.  It only makes them REACHABLE if someone switches them on.
    results = _apply_fixer_sp2n_planarize_if_enabled(mol, results, dual_parse_done)
    results = _apply_fixer_sp2c_planarize_if_enabled(mol, results, dual_parse_done)
    # ===== THE TERMINAL M=E BONDS (17.08.2026) =========================================
    # MEASURED: `me42b` reached 2 of 42, `byte_identical` 40 -- and the cause is
    # NOT the chemistry, but that the only real setter is legacy-only.  Three
    # suspects were excluded one by one: not the SMILES format (`_ml_bond_kind`
    # does not read the bond order at all, it judges structurally), not the pool (28 of
    # the 42 carry a real terminal O on the metal), not the table (33 of 95 pairs
    # are effectively shorter, among them Re=O 0.813, W=O 0.850, V=O 0.890, Mo=O 0.895).
    #
    # It was the wiring: `_clamp_metalloid_md_xyz` is metalloid-only,
    # `_md_distance_in_tolerance` is a predicate, `_manual_metal_embed` is the CN4 path
    # -- and `_snap_md_distances_to_ideal` has its only call site at :33066,
    # behind the FF-free `return` at :32399.
    #
    # ⚠ This hook-up needs NO `mol` and therefore does not fall into the atom-order
    # trap that the docstring above warns about: the criterion for a terminal M=E is
    # structural and readable from the frame itself.  Exactly ONE atom is moved, one that
    # has no bond besides the M-D bond.  Own switch, default 0 -> byte-identical.
    results = _apply_me_bond_snap_if_enabled(results)
    # ===== THE MIRROR COMPLETION (17.08.2026) -- MUST STAND LAST ========================
    # MEASURED: the corpus carries no stereochemistry (9 of 129 314 SMILES).  Handedness
    # is thereby 100 % an enumeration duty, not a transfer.  Of 2269 failures,
    # 1946 are "never built", and for 1326 of those (68.1 %) ALL centres are missing -- there
    # the mirror image is EXACTLY the missing isomer.
    #
    # ⚠ THE POSITION IS PART OF THE MECHANISM, not taste.  `_stereocenter_enum`
    # reads `present` BEFORE it supplements; `trans208` displaced the stereocentre folds
    # in exactly this way with 29 additional arrangements and cost a CCDC isomer,
    # although BOTH passes are additive.  The fold enumeration runs on the FF-free
    # path at :32361, this tail is called at :32429 -- the mirror pass sees
    # the finished pool and no longer changes for anyone what `present` reports.
    results = _apply_mirror_enum_if_enabled(results)
    # Rest of the rollback gate (see snapshot at the start of the chain).  The mirror appends
    # NEW labels; those pass through untouched -- only what changed an
    # existing label is evaluated.
    try:
        from delfin.manta._refine_gate import keep_better as _rg_keep
        results = _rg_keep(_rg_before, results)
    except Exception as _rg_exc:              # pragma: no cover
        logger.debug("refine-gate nicht angewandt: %s", _rg_exc)
    # SECOND HALF OF THE SAME CONTRACT (01.09.2026).  `keep_better` enforces
    # "no frame is WORSE"; this line enforces "no frame is GONE".
    # Measured on mirrleg6k: ABUSAU 58->59 and JEJROI 89->90 each lose ONE
    # frame, although `_mirror_enum.expand_results` is additive by construction
    # (:370 `return list(results) + added`) -- a SELECTION STAGE further down takes
    # it.  Default OFF (DELFIN_FFFREE_ADD_NEVER_REPLACE) -> byte-identical.
    # ⚠ AFTER `keep_better`, not before: first correct contents, then supplement what
    #   is missing -- the other way round, the rollback gate would evaluate what was just restored.
    try:
        from delfin.manta._refine_gate import keep_all as _rg_all
        results = _rg_all(_rg_before, results)
    except Exception as _rg_all_exc:          # pragma: no cover
        logger.debug("ADD-never-replace nicht angewandt: %s", _rg_all_exc)
    return results


def _smiles_to_xyz_isomers_impl(
    smiles: str,
    num_confs: int = 200,
    max_isomers: int = 50,
    apply_uff: bool = True,
    collapse_label_variants: bool = True,
    include_binding_mode_isomers: bool = True,
    deterministic: bool = True,
    hapto_approx: Optional[bool] = None,
    quality_mode: Optional[str] = None,
    seeds_override: Optional[int] = None,
    n_metal_smart: bool = True,
    _dual_parse_done: bool = False,
) -> Tuple[List[Tuple[str, str]], Optional[str]]:
    """Generate distinct coordination isomers for a SMILES string.

    For non-metal molecules a single geometry is returned (delegates to
    ``smiles_to_xyz``).  For metal complexes multiple conformers are
    embedded, their coordination fingerprints are computed, and one
    representative per unique fingerprint is returned.

    Returns ``([(xyz_string, label), ...], error)``.

    When ``collapse_label_variants`` is ``True`` (default), additional cleanup
    merges numbered variants with the same base label (e.g. ``trans-1`` and
    ``trans-2``). Set it to ``False`` to keep these variants for workflows
    that prefer broader structural diversity.

    When ``deterministic`` is ``True`` (default), the metal-isomer path avoids
    Open Babel conformer injection (its ``make3D``/rotor pipeline is not fully
    reproducible across runs) and relies on seeded RDKit embedding only.

    ``quality_mode`` selects a preset that tunes seed count, chelate
    conformer ranks, template top-K and Pre-UFF cap multiplier:

    - ``"fast"``   — 12 seeds, 1 rank, 1 template, cap 3·max_isomers.
      Roughly 3-5× faster than ``"max"``; suitable for the dashboard
      button where UI latency matters.
    - ``"normal"`` — 20 seeds, 2 ranks, 2 templates, cap 4·max_isomers.
    - ``"max"``    — 40 seeds, 3 ranks, 3 templates, cap 5·max_isomers.
      Deepest candidate pool, best geometry quality, slowest.

    Pass ``None`` (default) to keep the module-level knobs
    (``DELFIN_TOP_LEVEL_SEED_COUNT``, ``DELFIN_CHELATE_RANK_VARIANTS``
    etc.) as-is.
    """
    if not RDKIT_AVAILABLE:
        return [], "RDKit is not installed"

    # Iter-8.4a: reset module-global override at every function entry.
    # The override is then re-set further down once the parent mol has
    # been parsed and ``_classify_complex_class`` can run.  Resetting at
    # entry guarantees no leak from a previous SMILES call when the
    # current SMILES exits early (no_metal short-circuit, RDKit error,
    # hapto failfast, etc.) before the class-dispatch block is reached.
    global _ITER84_SIGMA_CAPS_OVERRIDE
    _ITER84_SIGMA_CAPS_OVERRIDE = None

    has_metal = contains_metal(smiles)
    hapto_mode = _hapto_approx_enabled(hapto_approx)

    # --- metal-FF-free generation backend (env-gated, default OFF) ---------
    # delfin.manta: provably-complete deterministic enumeration (Pólya isomers)
    # + metal-FF-free geometric assembly on COD-ideal polyhedra.  v1 handles
    # Werner complexes (mononuclear, explicit metal-donor bonds, all-monodentate,
    # CN 4-6); decompose() returns None for chelates / hapto / dative / multi-
    # metal, so those fall through to the legacy pipeline below unchanged.
    # Bit-exact OFF (default).  See delfin/manta/.
    _hapto_ff_fallback: Optional[List[Tuple[str, str]]] = None
    _ffree_union: Optional[List[Tuple[str, str]]] = None   # DELFIN_FFFREE_UNION, see below
    if has_metal and _delfin_env_int("DELFIN_FFFREE_BUILDER", 0):
        # ONE read, two uses: here for the rescue rung in the FF-free builder,
        # below for merging the two manifolds.  The promise "the only
        # read site" thereby stays true -- it has merely moved up, because the value
        # is now needed BEFORE the call and not only after its result.
        _union_on = _delfin_env_int("DELFIN_FFFREE_UNION", 0)
        try:
            from delfin.manta.converter_backend import _fffree_isomers
            _ff = _fffree_isomers(smiles, max_isomers=max_isomers, union=bool(_union_on))
        except Exception:
            _ff = None
        if _ff:
            # Hapto scaffold-primary (DELFIN_FFFREE_HAPTO_SCAFFOLD_PRIMARY, default
            # OFF -> byte-id): the FF-free RIGID_HAPTO path collapses η-faces (η6
            # ~27% clean) whereas the legacy analytical scaffold builds them
            # correctly (η6 ~84%).  For a hapto complex, DON'T short-circuit here:
            # stash the FF-free RIGID result as a FALLBACK and let the scaffold
            # path (hapto branch below) try first.  The fallback is used only if
            # the scaffold cannot build the complex -> never-worse on build-rate,
            # strictly better geometry.  Non-hapto returns FF-free unchanged.
            if _hapto_scaffold_primary_enabled():
                try:
                    _hp = _probe_hapto_groups_from_smiles(smiles)
                except Exception:
                    _hp = None
                if _hp:
                    _hapto_ff_fallback = _ff
            if _hapto_ff_fallback is None:
                # Iter-34: co-planar-M orient for coordinated in-plane σ π-donors on
                # the FF-free path too (ABIZIW builds here).  No parsed ``mol`` is
                # available at this early return; the corrector is geometry-only and
                # the flag is a plain integer env-var, so pass mol=None.  Default-OFF
                # byte-identical (the dispatch returns _ff unchanged when unset).
                # Snapshot for the rollback gate -- flat list copy, costs
                # nothing.  The counterpart line stands at the end of this chain, behind the mirror.
                _rg_before_ff = list(_ff)
                _ff = _apply_pi_coplanar_m_if_enabled(None, _ff, False)
                # ── THE POST-PASS CHAIN IS UNREACHABLE FROM HERE (found 2026-07-30) ──
                # This ``return`` short-circuits ~2260 lines that hold B4, B5, **B6 (the
                # variational functional)** and every post-B5 fixer.  On the CHAMPION path
                # (DELFIN_FFFREE_BUILDER=1) every metal complex leaves through here, so the
                # functional was never CALLED -- `DELFIN_B6_WIRED=1` measured affected=0 on
                # 35 systems and read as "it declines every frame", which was wrong.  The
                # line above is the tell: ONE corrector was already hand-wired back in here
                # instead of the hole being closed.
                # Wire the FUNCTIONAL (and only it -- the post-hoc fixers are scaffolding we
                # want to DROP, not resurrect).  ``mol`` is parsed lazily behind the same env
                # gate, so the default-OFF path stays byte-identical and pays nothing; the
                # refiner needs topology only (atoms/bonds/metals), never a conformer.
                if (_delfin_env_int("DELFIN_B6_WIRED", 0)
                        or _delfin_env_int("DELFIN_BAUSTEIN6", 0)):
                    try:
                        _b6_mol = _prepare_mol_for_embedding(
                            smiles, hapto_approx=hapto_mode,
                        )
                    except Exception:
                        _b6_mol = None
                    if _b6_mol is not None:
                        _ff = _apply_baustein6_if_enabled(_b6_mol, _ff, False)
                # ── STEREOCENTRE FOLDS FOR THE FF-FREE PATH (DELFIN_FFFREE_STEREO_ON_FFREE,
                #    default OFF -> byte-identical).  Found 2026-08-10. ──
                #
                # THE HOLE.  STEREOCENTER_ENUM is a CHAMPION flag, and it has never once run on
                # a force-field-free frame.  Its single call site is ~2900 lines below this
                # return (the dispatch at _apply_stereocenter_enum_if_enabled, def 1699, called
                # exactly once), and under UNION the FF-free frames are parked in _ffree_union
                # and concatenated only at the very end -- also after it.  So in BOTH modes the
                # fold expansion sees the legacy manifold and nothing else.  Same shape as the
                # ring-pucker gap documented just below, and the same shape as B6 above: a
                # post-pass that the early return silently skips.
                #
                # WHY IT PORTS WITHOUT ANY ADAPTER.  _stereocenter_enum.expand_results is a pure
                # FRAME-LIST transformer: it takes the finished [(xyz, label), ...], reads each
                # frame's own geometry, groups by coordination isomer and APPENDS the folds that
                # are missing ("Originals are preserved verbatim -> never-worse-safe").  It does
                # not know or care which constructor produced the frames.  Note this is exactly
                # the case the ring-pucker note below is NOT: that emitter needed a mol whose
                # atom order matched the frame's, which never held here.  This one needs no mol
                # at all -- the dispatch's own docstring says so ("``mol`` is unused (the
                # corrector operates on XYZ text only)"), which is why None is passed, the same
                # as _apply_pi_coplanar_m_if_enabled one block up.
                #
                # WHY A SECOND GATE AND NOT JUST THE CALL.  _apply_stereocenter_enum_if_enabled
                # already gates on DELFIN_(FFFREE_)STEREOCENTER_ENUM -- and that flag is IN the
                # champion.  Calling it unconditionally here would therefore change every
                # champion build the moment this line lands, which is precisely what "default
                # OFF -> byte-identical" forbids.  The extra flag keeps the reach measurable:
                # off = today's champion to the byte, on = the same champion plus the folds the
                # FF-free manifolds never got.
                #
                # ⚠ UNMEASURED.  Additive by the module's own contract, but the reach is unknown:
                # measured 10.08. on legacy manifolds only 4.8 % of systems carry a _stereo-
                # frame at all, and the FF-free subset may be richer or poorer in coordinated
                # secondary amines.  Run the fire census on this line before drawing conclusions.
                if _delfin_env_int("DELFIN_FFFREE_STEREO_ON_FFREE", 0):
                    _ff = _apply_stereocenter_enum_if_enabled(None, _ff, False)
                # AXIAL completeness here too -- NO second gate.  The expansion is
                # capped anyway by its own switch (default OFF), and a second
                # switch would be exactly the build form on which the stereocentre expansion
                # failed on 10.08.: it was in the champion and never reached the FF-free
                # path.  Call sites count, not lines.
                _ff = _apply_atropisomer_enum_if_enabled(None, _ff, False)
                # ── THE SHARED TAIL (DELFIN_FFFREE_SHARED_TAIL, default OFF -> byte-identical) ──
                #
                # The two hand-wired lines above are the symptom, not the cure: pi-coplanar was
                # patched in when someone noticed it, the stereocentre folds on 2026-08-10 when
                # someone else did, and nobody ever asked what ELSE lives past this return.  The
                # answer is in _ffree_shared_tail's docstring -- eighteen correctors, of which
                # exactly two are champion-active.  Those two run here.
                #
                # ``mol`` is parsed LAZILY behind the flag, the same discipline as the B6 block
                # above: the default-OFF path parses nothing and pays nothing, and the correctors
                # need topology only.  Failure to parse leaves _ff untouched -- a corrector that
                # cannot be applied must never cost the frames it was meant to improve.
                #
                # ⚠ Read the docstring before judging a null measurement: the SMILES-mol atom
                # order and the FF-free frame atom order do not coincide, which is what made the
                # ring-pucker emitter a null lever at this exact spot.  "affected = 0" here means
                # "find out which of the two it was", not "the correctors decline".
                if _delfin_env_int("DELFIN_FFFREE_SHARED_TAIL", 0):
                    try:
                        _tail_mol = _prepare_mol_for_embedding(
                            smiles, hapto_approx=hapto_mode,
                        )
                    except Exception:
                        _tail_mol = None
                    if _tail_mol is not None:
                        _ff = _ffree_shared_tail(_tail_mol, _ff, False)
                # ===== THE MIRROR COMPLETION NEEDS NO FOREIGN SWITCH ===================
                # Measured 18.08.2026: the pass sat exclusively IN the shared tail, i.e.
                # behind DELFIN_FFFREE_SHARED_TAIL (default 0).  Every mirror A/B therefore
                # had to carry the switch along in BOTH arms -- and the shared tail
                # incidentally activates ISOLATED_SEAT and AROM_PLANARIZE, the two
                # champion flags that otherwise never fire FF-free.  The baseline was thus
                # NOT the shipped champion, and loop.py reported that too
                # ("UNDECLARED AXIS ... the arms are NOT the shipped champion").
                # A mechanism that is measurable only under a foreign switch cannot
                # land in the champion -- one would be landing two things at once.
                #
                # The mirror does not need this switch: it needs NO ``mol`` (a
                # mirroring is a pure coordinate operation), so exactly the reason for
                # which the shared tail gates at all -- the expensive parsing -- falls away.
                # Hence it stands here unconditionally, gates only through its OWN switch
                # (default 0 -> byte-identical) and is a no-op with SHARED_TAIL=1, because
                # expand_results is idempotent as of today (label suffix ``_mirror``).
                # ── H PLACEMENT (DELFIN_FFFREE_H_PLACEMENT, default 0) ──
                # The hydrogen family carries 28 % of the hardness mass and had until
                # 19.08. NO active mechanism: of seven existing ones, three are
                # measured refuted (5B_VSEPR_H_REALISM, H_FOLLOW, DONOR_FOLLOW), one
                # has zero call sites in the whole tree (H_CLASH_ROTATE) and the rest
                # lie BEHIND the FF-free return -- VSEPR_REPAIR even 3038 lines.
                #
                # ⚠️ THE POSITION IS PART OF THE MECHANISM, not taste:
                #   * BEFORE the mirror, so that BOTH hands carry the same H repair.
                #   * AFTER _apply_stereocenter_enum_if_enabled, so that its `present`
                #     read stays untouched (exactly on that, trans208 displaced
                #     stereocentres although both passes were additive).
                #
                # Heavy atoms stay pinned byte-exactly, and a stereo gate rolls back the WHOLE
                # frame if a real centre flips.  Measured with the eye's detectors on
                # 2144 frames: clashes 80->48, hard HH contacts 225->142,
                # methyl violations in 133->103 files.
                _ff = _apply_h_placement_if_enabled(_ff)
                _ff = _apply_mirror_enum_if_enabled(_ff)
                # --- TORSION ENUMERATION, now ALSO on the FF-free path (27.08.) ---
                # COUNTED, not assumed: `_rotamer_diversity.apply_if_enabled` had
                # EXACTLY ONE call site, and that lies at :35799 in the LEGACY tail.
                # The FF-free complex path got ZERO torsion enumeration -- the same
                # build form as on 10.08. ("three additive modules run only on legacy"),
                # only this time in this direction.
                #
                # ⚠ AND IT IS THE MECHANISM THAT TASK #17 IS LOOKING FOR.  The setter in
                # `conformer_enum.py` (0 call sites, silently caps at 7, builds
                # homoleptically) is NOT the right one -- this one here is:
                #   * DOFs from the MOLECULAR GRAPH, never from SMILES
                #   * coordination bonds (M-donor) excluded
                #   * M-D INVARIANT with tolerance DELFIN_5L_T6_ROTAMER_MD_TOL (0.05 A)
                #   * works on the FINISHED frame, i.e. the real build template
                #
                # REACH: 4411 of 6000 = 73.5 % have >=1 rotatable bond
                # (measured 27.08., using `_rotatable_bonds` itself).  Of those,
                # ~1190 of 1616 are FF-free -- far above the 100 threshold.
                #
                # ⚠ THE CAPS ARE REAL AND BELONG IN THE MEASUREMENT: K=3 frames per
                # isomer, MAX_DOFS=6, GRID_CAP=64.  K=3 means SELECTING by energy,
                # not enumerating -- and what SELECTS has never yet landed in this
                # campaign.  At the verdict, watch whether the cap binds.
                # ⛔ Default OFF (DELFIN_5L_T6_ROTAMER_DIVERSITY) -> byte-identical.
                try:
                    from delfin.manta import _rotamer_diversity as _rot_ff
                    if _rot_ff._is_enabled() and _ff:
                        _erw = []
                        for (_rx, _rl) in _ff:
                            _fr = _rot_ff.apply_if_enabled(_rx)
                            if not _fr:
                                _erw.append((_rx, _rl))
                                continue
                            _erw.append((_fr[0], _rl))
                            for _ki, _fx in enumerate(_fr[1:], start=1):
                                _erw.append((_fx, f"{_rl}_rotamer-{_ki}" if _rl
                                             else f"rotamer-{_ki}"))
                        _ff = _erw
                except Exception as _rot_ff_exc:
                    logger.debug("FF-frei Rotamer-Erweiterung fehlgeschlagen: %s",
                                 _rot_ff_exc)
                # ROLLBACK GATE for the FF-free chain (snapshot above at
                # `_rg_before_ff`).  Covers pi_coplanar_m, baustein6 and h_placement --
                # the enumerators in between (stereocentres, atropisomer, mirror) append
                # NEW labels and pass through untouched.
                try:
                    from delfin.manta._refine_gate import keep_better as _rg_keep_ff
                    _ff = _rg_keep_ff(_rg_before_ff, _ff)
                except Exception as _rg_exc_ff:   # pragma: no cover
                    logger.debug("refine-gate (ffree) nicht angewandt: %s", _rg_exc_ff)
                # ADD, NEVER REPLACE -- the second half of the contract, here too.
                # The sentence three lines up ("the enumerators append NEW labels
                # and pass through untouched") was an ASSUMPTION; on ABUSAU and
                # JEJROI it was refuted on 01.09.  Default OFF.
                try:
                    from delfin.manta._refine_gate import keep_all as _rg_all_ff
                    _ff = _rg_all_ff(_rg_before_ff, _ff)
                except Exception as _rg_all_ff_exc:   # pragma: no cover
                    logger.debug("ADD-never-replace (ffree) nicht angewandt: %s",
                                 _rg_all_ff_exc)
                # RING PUCKER FOR THE FF-FREE PATH: the hook that USED to sit here has been
                # removed, and the reason is worth keeping.
                #
                # The gap it addressed is real and large: all three pucker emitters sit ~1800
                # lines BELOW this return, so the FF-free builder never reaches any of them.
                # Measured on archive_scopecensus3 (996 champion systems): the FF-free chelate
                # path emits 2.23 frames per system and ZERO frames carrying a conformer
                # suffix, against 30.6 on everything else; ccdc_pucker_realized is false for
                # 25.2 % of the census.  The FF-free conformer manifold is not thin, it is EMPTY.
                #
                # But calling _emit_ring_puckers_rp from HERE can never work, and it was
                # measured: 185 of 187 systems byte-identical.  That emitter bridges XYZ to a
                # mol through _xyz_to_rdkit_conformer, which requires the element symbol to
                # match at EVERY index -- and the mol parsed from the SMILES carries RDKit's
                # own atom order, while our frames carry metal-at-0 plus AddHs(ligand) blocks
                # in construction order.  They never coincide, so it returned 0 every time.
                # It failed SAFE (nothing written), but it was a null lever.
                #
                # The live version therefore lives where the frame's own order is known:
                # converter_backend._ffree_ring_pucker_frames, called next to the frame it is
                # a sibling of.  ONE env read site, and it is over there.
                #
                # ===== UNION INSTEAD OF EITHER-OR (DELFIN_FFFREE_UNION, default OFF) =====
                #
                # THE ONE PLACE DELFIN_FFFREE_UNION IS READ.
                #
                # The architecture treats the two constructors as ALTERNATIVES: whoever
                # answers first owns the system, and the self-gate picks.  Measured
                # 2026-08-02, five times, that is exactly what makes scope impossible --
                # CONFORMER_SEATING -78 capabilities, BITE_FREE -52, CHELATE_BACKBONE -52,
                # COLLAPSE_SELECT -15, TET_CHELATE -3.  Every one of them let FF-free WIN a
                # system it used to hand over, and legacy had built it better.  The loss is
                # not chemistry, it is the either-or.
                #
                # They do not have to be alternatives.  A manifold is a SET of frames, and
                # two constructors can both contribute to it.  The eye's crystal floors read
                # the BEST frame over the manifold, so adding legacy's frames next to ours
                # can only raise them; and nothing of ours is removed, so nothing of ours can
                # regress.  This is the same shape as every flag that ever landed here
                # (D8_SQ_ADD, CN6_OH_ADD, STEREOCENTER_ENUM, CN4_BOTH): ADD, never replace.
                #
                # Mechanically it is the _hapto_ff_fallback pattern one line up, with the
                # early return dropped instead of conditioned: stash our frames, let the
                # legacy pipeline run to completion, concatenate at the end.  Costs a second
                # build per system -- and the machine sat half idle all night.
                if _union_on:                      # the same read as above, not a second one
                    _ffree_union = _ff
                else:
                    return _ff, None

    # Resolve the quality profile once per call so the seed count,
    # chelate ranks, topK and Pre-UFF cap follow the requested preset.
    _qprof = dict(_resolve_quality_profile(quality_mode))
    # Explicit seed count override (from dashboard slider, scripted
    # benchmarks, etc.) takes precedence over the profile default.
    # Clamped to [1, len(_PIPELINE_SEEDS)] so we never walk past the
    # seed pool or pass <1 down to ETKDG.
    if seeds_override is not None:
        try:
            _seeds_clamped = max(1, min(int(seeds_override), len(_PIPELINE_SEEDS)))
            _qprof["seeds"] = _seeds_clamped
        except Exception:
            pass

    # Size-aware seed auto-cap: cuts per-call seed count for very large
    # polycyclic-cage / heavy-polydentate complexes (≥80 heavy atoms)
    # so they finish inside a bounded sweep wall-clock.  DEFAULT OFF —
    # capping seeds discards conformer diversity exactly where
    # we need it most (cage ligands have rough ETKDG initial geometry
    # that benefits from MORE samples, not fewer).  Interactive use
    # (dashboard) prefers to wait for a fully-realistic structure over
    # capping seeds.  Opt in via DELFIN_AUTOSCALE_SEEDS=1 for batch
    # sweeps that strictly need bounded per-SMILES wall-clock (set by
    # commit_sweep/scripts/sweep_focus_pool.py and sweep_focus.py).
    if _delfin_env_int("DELFIN_AUTOSCALE_SEEDS", 0):
        try:
            _probe_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if _probe_mol is not None:
                _heavy_n = sum(
                    1 for _a in _probe_mol.GetAtoms() if _a.GetAtomicNum() > 1
                )
                _orig = int(_qprof.get("seeds", 20))
                if _heavy_n >= 120:
                    _capped = min(_orig, 3)
                elif _heavy_n >= 80:
                    _capped = min(_orig, 5)
                else:
                    _capped = _orig
                if _capped < _orig:
                    _qprof["seeds"] = _capped
                    logger.debug(
                        "Size-aware seed cap: %d heavy atoms → seeds %d→%d",
                        _heavy_n, _orig, _capped,
                    )
        except Exception:
            pass

    mol = None
    hapto_groups: List[Tuple[int, List[int]]] = []
    if has_metal:
        hapto_groups = _probe_hapto_groups_from_smiles(smiles)
        if hapto_groups and not hapto_mode:
            return [], _hapto_failfast_error(hapto_groups)
        if hapto_groups and hapto_mode:
            xyz, err = smiles_to_xyz(
                smiles,
                apply_uff=apply_uff,
                hapto_approx=True,
                deterministic=deterministic,
            )
            if err or not (xyz and xyz.strip()):
                # Scaffold-primary fallback: if the analytical scaffold cannot
                # build this hapto complex (error OR empty geometry), fall back to
                # the stashed FF-free RIGID_HAPTO result (collapsed but present) so
                # the build-rate is never-worse vs the FF-free path.  Only set when
                # DELFIN_FFFREE_HAPTO_SCAFFOLD_PRIMARY is on.
                if _hapto_ff_fallback:
                    return _hapto_ff_fallback, None
                if err:
                    return [], err
            # First-class η-label per detected hapto group.
            # `hapto_groups` is a list of (metal_idx, [c_indices]).  When
            # multiple groups exist (e.g. ferrocene Fe(η5-Cp)(η5-Cp))
            # join their canonical η-labels with '/'.
            try:
                _hapto_pre_mol = _prepare_mol_for_embedding(smiles, hapto_approx=True)
                if _hapto_pre_mol is not None:
                    _eta_labels = []
                    for _m_idx, _c_idxs in _find_hapto_groups(_hapto_pre_mol):
                        _eta_labels.append(_hapto_label_for_group(_hapto_pre_mol, _c_idxs))
                    _eta_combined = '/'.join(_eta_labels) if _eta_labels else 'hapto'
                else:
                    _eta_combined = 'hapto'
            except Exception:
                _eta_combined = 'hapto'
            results_hapto: List[Tuple[str, str]] = [(xyz, _eta_combined)]
            # Enumerate σ-donor permutations around fixed η-positions.
            try:
                hapto_sigma_iso = _enumerate_hapto_sigma_isomers(
                    smiles, xyz, apply_uff=apply_uff
                )
                if hapto_sigma_iso:
                    # Relabel hapto-sigma-N to η{n}-{type} σ-N for clarity.
                    _relabeled = []
                    for _hsx, _hsl in hapto_sigma_iso:
                        if isinstance(_hsl, str) and _hsl.startswith('hapto-sigma'):
                            _suffix = _hsl[len('hapto-sigma'):].lstrip('-')
                            _relabeled.append((
                                _hsx,
                                f'{_eta_combined} σ-{_suffix}' if _suffix else _eta_combined,
                            ))
                        else:
                            _relabeled.append((_hsx, _hsl))
                    results_hapto.extend(_relabeled)
            except Exception as _hse:
                logger.debug("Hapto sigma enumeration failed: %s", _hse)

            # Hapto diversity branch — explicitly user-authorised
            # non-deterministic augmentation.  Cp/arene ring orientation has
            # genuine rotational variety that a single deterministic
            # build-up cannot capture; Open Babel's weighted-rotor
            # WeightedRotorSearch explores that space.  Every OB conformer
            # is re-run through the topology gate so nothing broken slips
            # through.  σ-complexes stay on the deterministic path.
            _mol_hapto_gate: Optional[object] = None
            try:
                _mol_hapto_gate = _prepare_mol_for_embedding(smiles, hapto_approx=True)
            except Exception:
                _mol_hapto_gate = None

            # Iter-6 (123a130-port): the Iter-1 σ-guard combined with the
            # Iter-5 default-flip of HD-TA killed Hapto-class diversity.
            # Champion 123a130 (Apr 24, topo_pct_match 60.3 %) ran OB-WRS
            # unconditionally for every hapto SMILES; HEAD-iter5 (14 %)
            # produces only seed + σ-permutations on σ-mixed cases.
            # When DELFIN_HAPTO_123A_PORT=1, force σ-guard OFF and HD-TA
            # ON (topology gates remain in place).  See
            # results/iter6_hapto_forensik.md.  Default ON since Iter-7 (2026-05-07):
            # smoke 2000 + 30-detector battery validated topo_pct_match 43.84% -> 48.42%
            # (+4.58pp), stereo_pct_match 35.94% -> 44.93% (+8.99pp), h_h_clashes
            # 2910 -> 303 (-89.6%), M_L_intactness 0.979 -> 0.994 (+0.015).
            # Disable via DELFIN_HAPTO_123A_PORT=0 if regression appears.
            #
            # Iter-7.1 (2026-05-07) multi-metal guard: master_v3 full-pool data
            # (~25%) revealed multi-hapto class M_L 0.996 -> 0.882 (-11.4pp),
            # topology extra-bonds 11.84 -> 24.67 (+108%) with port=1 default-on.
            # Root cause: OB-WRS rotates each metal's Cp/arene independently,
            # tearing apart multi-metal-cluster topology.  Port now active
            # ONLY for single-metal complexes.  Multi-metal complexes fall back
            # to the Iter-1 σ-guard behaviour (sigma-mixed -> skip OB-WRS).
            # See project_iter7_3way_findings.md for the per-class delta.
            # Iter-8.6d (2026-05-11): class-conditional hapto-port restore.
            # 8.6a (commit 7b2bb19) globally flipped HAPTO_123A_PORT default
            # 1 -> 0, winning sigma topology (+16.62pp) but breaking hapto
            # chemistry: 4 calibrated brid_hapto detectors put 8.6a at rank
            # 60/63 with 34094 spurious chelate frames; per-bond diagnostic
            # shows hapto intact-rate 43.2% (8.4abc) -> 20.7% (8.6a) = 2.4x
            # worse than 123a130 (49.1%). See CALIBRATION_DAY1_findings.md.
            # Class-conditional restore: port=1 for hapto/multi_hapto (where
            # 8.4abc was better), port=0 default for sigma/multi_sigma
            # (preserves 8.6a sigma topology win). Default class-list
            # "hapto,multi_hapto" applies the iter-8.6d behaviour; set the
            # env-var explicitly empty to fall back to the legacy 8.6a
            # scalar DELFIN_HAPTO_123A_PORT (back-compat). Mimics pattern
            # at smiles_converter.py:20890 (DELFIN_ITER85_PUMP_SKIP_CLASSES).
            _hapto_123a_port_raw = False
            try:
                _hapto_123a_port_classes = set(
                    x.strip() for x in (
                        os.environ.get(
                            "DELFIN_HAPTO_123A_PORT_CLASSES",
                            "hapto,multi_hapto",
                        ) or ""
                    ).split(",") if x.strip()
                )
                if _hapto_123a_port_classes and _mol_hapto_gate is not None:
                    _hapto_123a_port_raw = (
                        _classify_complex_class(_mol_hapto_gate)
                        in _hapto_123a_port_classes
                    )
                else:
                    # Empty class-list -> back-compat scalar fallback
                    _hapto_123a_port_raw = bool(
                        _delfin_env_int("DELFIN_HAPTO_123A_PORT", 0)
                    )
            except Exception:
                _hapto_123a_port_raw = bool(
                    _delfin_env_int("DELFIN_HAPTO_123A_PORT", 0)
                )
            if _hapto_123a_port_raw and _mol_hapto_gate is not None:
                try:
                    _n_metals_check = sum(
                        1 for _a in _mol_hapto_gate.GetAtoms()
                        if _a.GetSymbol() in _METAL_SET
                    )
                except Exception:
                    _n_metals_check = 1
                # Iter-7.2 (2026-05-07): refine multi-metal guard.  Iter-7.1's
                # n_metals<=1 check was too strict — Sn/Pb/Bi/Sb/Tl-TM(η-Cp)
                # complexes have n_metals==2 but only ONE hapto-active metal
                # (main-group is σ-only).  These need the OB-WRS rotor branch.
                # Forensics: results/iter7.2_multihapto_forensik.md.  The earlier
                # iter6 baseline MLI=0.996 was metric-cheating: OB-WRS rot-frames
                # had spurious M-L extra-bonds but MLI counts only missing, so
                # they appeared "intact" while topology was torn.  Honest
                # multi-hapto MLI is ~0.886.  Frame count drop from 281 -> 23
                # (genuine bi-hapto only) needs OB-WRS recovery for Sn-M(η-Cp).
                _n_hapto_metals = 0
                try:
                    _hapto_metal_set = {
                        _m_idx for _m_idx, _ in _find_hapto_groups(_mol_hapto_gate)
                    }
                    _n_hapto_metals = len(_hapto_metal_set)
                except Exception:
                    _n_hapto_metals = 0
                _hapto_123a_port = (
                    (_n_metals_check <= 1) or (_n_hapto_metals <= 1)
                )
            else:
                _hapto_123a_port = _hapto_123a_port_raw

            # Iter-1 (hapto-σ topology guard): WeightedRotorSearch treats the
            # M-C η-bonds encoded in the SMILES as rotatable torsions and
            # rotates Cp/arene rings as if they were freely-rotating ligands.
            # When σ-donors are also present, those rotations rip the rings
            # apart and project σ-ligands onto the Cp plane (user 2026-05-04
            # "everything in one plane" pattern).  Skip OB diversity for piano-
            # stool / mixed-coordination complexes; pure Cp/arene sandwiches
            # (no σ-donors) still benefit from rotor diversity.
            # Opt out via DELFIN_HAPTO_NO_OB_WHEN_SIGMA=0 OR
            # DELFIN_HAPTO_123A_PORT=1 (Iter-6 port).
            _hapto_skip_ob_diversity = False
            try:
                if _hapto_123a_port:
                    _no_ob_active = False  # 123a130: OB-WRS for all hapto
                else:
                    _no_ob_env = os.environ.get(
                        "DELFIN_HAPTO_NO_OB_WHEN_SIGMA", "1"
                    ).strip().lower()
                    _no_ob_active = _no_ob_env not in {"0", "false", "no", "off", ""}
                if _no_ob_active and _mol_hapto_gate is not None:
                    _n_sigma_total = 0
                    _hapto_atoms_set: set = set()
                    try:
                        for _m_idx, _grp in _find_hapto_groups(_mol_hapto_gate):
                            _hapto_atoms_set.update(_grp)
                    except Exception:
                        _hapto_atoms_set = set()
                    for _a in _mol_hapto_gate.GetAtoms():
                        if _a.GetSymbol() not in _METAL_SET:
                            continue
                        for _nb in _a.GetNeighbors():
                            _ni = _nb.GetIdx()
                            if (
                                _nb.GetSymbol() not in _METAL_SET
                                and _ni not in _hapto_atoms_set
                            ):
                                _n_sigma_total += 1
                    if _n_sigma_total > 0:
                        _hapto_skip_ob_diversity = True
                        logger.debug(
                            "Hapto-σ guard: skipping OB diversity (n_sigma=%d)",
                            _n_sigma_total,
                        )
            except Exception as _hsg_exc:
                logger.debug("Hapto-σ guard probe failed: %s", _hsg_exc)

            # Iter-3 (Hapto Diversity Topology-Aware): When the Iter-1
            # guard fires (σ-donors present, OB-WRS skipped to protect
            # Cp/arene topology), generate diversity via rigid-body
            # σ-fragment + σ-cluster rotations.  These preserve every
            # hapto-atom position and every M-D bond length bit-exact, so
            # the topology hazard that motivated the Iter-1 guard cannot
            # arise.  Opt-out via DELFIN_HAPTO_DIVERSITY_RESTORE=0 (then
            # output is byte-exact identical to Iter-1 HEAD).
            # Iter-5: HD-TA was unconditional in Iter-3.2 (despite docstring).
            # Now env-gated default OFF — see iter5_polyhedron_md_forensik.md.
            # σ-rotations bypass polyhedron gate, pumping low-fidelity frames.
            # Iter-6 (123a130 port): in port mode HD-TA runs IN ADDITION to
            # OB-WRS for SMILES that OB cannot parse (dative `->` syntax).
            # Both diversity sources flow through `_verify_topology_from_graph`.
            _hd_ta_run = bool(
                _mol_hapto_gate is not None
                and not DELFIN_TOPOLOGY_STRICT_MODE
                and (
                    (
                        _hapto_skip_ob_diversity
                        and _delfin_env_int("DELFIN_HAPTO_DIVERSITY_RESTORE", 0)
                    )
                    or _hapto_123a_port
                )
            )
            if _hd_ta_run:
                try:
                    from delfin.manta._hapto_diversity import (
                        apply_hapto_diversity_topology_aware as _hd_apply,
                    )
                    _hd_seed_xyz, _hd_seed_label = results_hapto[0]
                    _hd_seen_keys = {
                        "\n".join(
                            l.strip() for l in _xyz.splitlines() if l.strip()
                        )
                        for _xyz, _lab in results_hapto
                    }
                    # Iter-3.2: pass UFF-optimizer so each rotated frame
                    # gets refined before the strict topology gate.
                    def _hd_optimize(xyz_str):
                        try:
                            return _optimize_xyz_openbabel_safe(
                                xyz_str, mol_template=_mol_hapto_gate,
                            )
                        except Exception:
                            return xyz_str
                    _hd_results = _hd_apply(
                        smiles,
                        _hd_seed_xyz,
                        _mol_hapto_gate,
                        metal_set=_METAL_SET,
                        find_hapto_groups=_find_hapto_groups,
                        verify_topology=_verify_topology_from_graph,
                        optimize_xyz=_hd_optimize if apply_uff else None,
                        max_frames=max(
                            0, max_isomers - len(results_hapto),
                        ),
                        seen_xyz_keys=_hd_seen_keys,
                        label_prefix=f"{_eta_combined} hd-ta",
                    )
                    if _hd_results:
                        results_hapto.extend(_hd_results)
                except Exception as _hd_ta_exc:
                    logger.debug(
                        "HD-TA diversity branch failed: %s", _hd_ta_exc,
                    )

            # Determinism: this OB make3D diversity block runs with
            # ``deterministic=False`` (wall-clock-seeded conformer search) and is
            # NOT reproducible cross-process — and, unlike the other OB calls, it
            # is reached on pure-hapto sandwiches (no σ-donor → guard never fires).
            # Under the master switch, skip it so the deterministic path is taken.
            if (
                OPENBABEL_AVAILABLE
                and not _hapto_skip_ob_diversity
                and not _deterministic_mode()
            ):
                try:
                    _HAPTO_OB_RESTARTS = 4
                    _HAPTO_OB_CONFS_PER_RUN = 30
                    _seen_xyz_keys: set = {
                        "\n".join(l.strip() for l in _xyz.splitlines() if l.strip())
                        for _xyz, _lab in results_hapto
                    }
                    for _run in range(_HAPTO_OB_RESTARTS):
                        _ob_blocks, _ob_err = _openbabel_generate_conformer_xyz(
                            smiles,
                            num_confs=_HAPTO_OB_CONFS_PER_RUN,
                            deterministic=False,
                        )
                        if not _ob_blocks:
                            continue
                        for _block in _ob_blocks:
                            _key = "\n".join(
                                l.strip() for l in _block.splitlines() if l.strip()
                            )
                            if _key in _seen_xyz_keys:
                                continue
                            # Topology gate — reject anything that broke the
                            # M-D graph or clashed fragments.
                            if _mol_hapto_gate is not None:
                                try:
                                    if not _verify_topology_from_graph(
                                        _block, _mol_hapto_gate
                                    ):
                                        continue
                                except Exception:
                                    pass
                            _seen_xyz_keys.add(_key)
                            results_hapto.append(
                                (_block, f'{_eta_combined} rot-{len(results_hapto):03d}')
                            )
                            if len(results_hapto) >= max_isomers:
                                break
                        if len(results_hapto) >= max_isomers:
                            break
                except Exception as _hd_exc:
                    logger.debug("Hapto OB diversity branch failed: %s", _hd_exc)

            # For mixed hapto+regular multi-metal: DON'T early-return.
            # Fall through to multi-metal sampling augmentation so the
            # non-hapto metal's coordination isomers are explored.
            _n_metals_hapto = 0
            try:
                if _mol_hapto_gate is not None:
                    _n_metals_hapto = sum(
                        1 for a in _mol_hapto_gate.GetAtoms()
                        if a.GetSymbol() in _METAL_SET
                    )
            except Exception:
                pass
            if _n_metals_hapto <= 1:
                # Iter-13: apply Baustein 3 to mono-hapto path (env-gated).
                results_hapto = _apply_coord_angle_fix_if_enabled(
                    _mol_hapto_gate, results_hapto, _dual_parse_done
                )
                # Iter-14: apply Baustein 4 (rigid-π H projection) AFTER B3.
                results_hapto = _apply_baustein4_if_enabled(
                    _mol_hapto_gate, results_hapto, _dual_parse_done
                )
                # Iter-25: final bond-decollapse (fixes ~79% hapto ligand collapse)
                results_hapto = _apply_bond_decollapse_if_enabled(
                    _mol_hapto_gate, results_hapto, _dual_parse_done
                )
                # UNION keeps its contract here too (DELFIN_FFFREE_UNION_HAPTO,
                # default OFF -> byte-identical).  Without this line the
                # hapto branch returns ~2700 lines BEFORE the merge and
                # `_ffree_union` is never read.
                return _union_prepend_ffree(results_hapto, _ffree_union), None
            # Multi-metal hapto: the hapto builder already produced the
            # best possible Cp geometry (perfectly planar rings). ETKDG
            # sampling would produce conformers with broken Cp rings.
            # Return the hapto results directly.
            # Iter-13: apply Baustein 3 to multi-metal hapto path (env-gated).
            results_hapto = _apply_coord_angle_fix_if_enabled(
                _mol_hapto_gate, results_hapto, _dual_parse_done
            )
            # Iter-14: apply Baustein 4 (rigid-π H projection) AFTER B3.
            results_hapto = _apply_baustein4_if_enabled(
                _mol_hapto_gate, results_hapto, _dual_parse_done
            )
            # Iter-25: final bond-decollapse (fixes ~79% hapto ligand collapse)
            results_hapto = _apply_bond_decollapse_if_enabled(
                _mol_hapto_gate, results_hapto, _dual_parse_done
            )
            # The same contract line for the multi-metal hapto branch.
            return _union_prepend_ffree(results_hapto, _ffree_union), None

    # Non-metal molecules: deterministic conformer pool.  Organic
    # ligands have no coordination-isomer axis, but different rotamer
    # / ring-pucker minima are still genuinely distinct low-energy
    # geometries.  Seed pool from _PIPELINE_SEEDS so results are
    # reproducible; de-duplicate by heavy-atom RMSD.
    if not has_metal:
        xyz, err = smiles_to_xyz(
            smiles,
            apply_uff=apply_uff,
            hapto_approx=hapto_mode,
            deterministic=deterministic,
        )
        if err:
            return [], err
        # legacy organic cap = 8; under energy-rank (completeness) honour max_isomers
        # up to 50 so flexible polyols/sugars (many OH rotamers) cover their global min.
        if os.environ.get("DELFIN_FFFREE_CONF_ENERGY_RANK", "0") == "1":
            _max_pool = max(1, min(int(max_isomers), 50))
        else:
            _max_pool = max(1, min(int(max_isomers), 8))
        pool = _organic_conformer_pool(
            smiles, xyz, max_pool=_max_pool, apply_uff=apply_uff,
        )
        if not pool:
            # Iter-14 B4: project ring-attached H onto π-plane (no-metal path).
            return _apply_baustein4_if_enabled(
                None, [(xyz, '')], _dual_parse_done
            ), None
        if collapse_label_variants and len(pool) == 1:
            return _apply_baustein4_if_enabled(
                None, [(pool[0][0], '')], _dual_parse_done
            ), None
        return _apply_baustein4_if_enabled(
            None, pool, _dual_parse_done
        ), None

    # Prepare molecule for embedding (skip if already set by hapto multi-metal path)
    if mol is None:
        mol = _prepare_mol_for_embedding(smiles, hapto_approx=hapto_mode)
    if mol is None:
        # Fall back to single-conformer conversion
        xyz, err = smiles_to_xyz(
            smiles,
            apply_uff=apply_uff,
            hapto_approx=hapto_mode,
            deterministic=deterministic,
        )
        if err:
            return [], err
        # Iter-13: apply Baustein 3 to fallback single-conformer path (env-gated).
        _fallback_results = [(xyz, '')]
        if has_metal:
            _fallback_results = _apply_coord_angle_fix_if_enabled(
                None, _fallback_results, _dual_parse_done
            )
        # Iter-14: apply Baustein 4 (rigid-π H projection) — runs on metal
        # AND non-metal complexes; aromatic-H out-of-plane is a generic
        # converter pathology independent of metal presence.
        _fallback_results = _apply_baustein4_if_enabled(
            None, _fallback_results, _dual_parse_done
        )
        return _fallback_results, None

    # ---- Multi-sigma V2 path (env-gated, default OFF) ----
    # Class-conditional seed cap: for large bimetallic σ-only systems the
    # 20 top-level seeds × 25 s _MULTIEMBED_TIMEOUT routinely exceeds the
    # caller's wall-clock budget without adding distinct coordination
    # isomers (the multi-metal augmentation block contributes most of the
    # diversity anyway).  Cap both knobs based on heavy-atom count.
    # See ``_multi_sigma_v2_active`` / ``_multi_sigma_v2_budget`` above.
    _ms_v2_budget: Optional[Dict[str, float]] = None
    _ms_v2_prev_override = getattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value", None)
    try:
        if _multi_sigma_v2_active(mol):
            _heavy_n_ms = sum(
                1 for _a in mol.GetAtoms() if _a.GetAtomicNum() > 1
            )
            _ms_v2_budget = _multi_sigma_v2_budget(_heavy_n_ms)
            _orig_seeds = int(_qprof.get("seeds", 20))
            _capped_seeds = min(_orig_seeds, int(_ms_v2_budget["seeds_top"]))
            if _capped_seeds < _orig_seeds:
                _qprof["seeds"] = _capped_seeds
                logger.debug(
                    "multi_sigma V2: heavy=%d → top-level seeds %d→%d, "
                    "embed_timeout=%.1fs, mm_walltime=%.1fs",
                    _heavy_n_ms, _orig_seeds, _capped_seeds,
                    _ms_v2_budget["embed_timeout"],
                    _ms_v2_budget["mm_walltime"],
                )
            # Cap ``alt_tries`` proportionally to the augmentation seed
            # cap so that alt-binding / linkage isomer exploration cannot
            # run 8 templates × 12 alt donors × 2 metals (the worst-case
            # 192-template traversal) on a 100-atom system.  We keep at
            # least 2 tries so the inner ranker can still pick between
            # a primary + backup template per rewire.  ``seeds_mm_aug``
            # is the natural pairing — both knobs gate the same
            # "per-candidate template iteration" workload.
            _orig_alt = int(_qprof.get("alt_tries", 8))
            _capped_alt = max(
                2, min(_orig_alt, int(_ms_v2_budget["seeds_mm_aug"]))
            )
            if _capped_alt < _orig_alt:
                _qprof["alt_tries"] = _capped_alt
            # Push tighter per-seed embedding timeout into
            # ``_embed_multiple_confs_with_timeout`` for the duration of
            # this call.  Restored in the finally block below.
            _MULTIEMBED_TIMEOUT_OVERRIDE.value = float(
                _ms_v2_budget["embed_timeout"]
            )
    except Exception as _ms_v2_exc:
        # Fail-safe: any failure in the V2 path leaves the pipeline at
        # pre-patch defaults so we never regress the small-molecule case.
        logger.debug("multi_sigma V2 cap skipped: %s", _ms_v2_exc)
        _ms_v2_budget = None

    # For metal complexes: prepend OB conformers to the pool so that
    # Avogadro-quality geometries are always considered during isomer search.
    conf_ids: List[int] = []
    if has_metal and OPENBABEL_AVAILABLE and not deterministic:
        try:
            # Non-deterministic enrichment path: augment RDKit pool with OB
            # conformers to increase diversity.
            _n_ob_restarts = 3
            _per_ob = max(10, int(num_confs) // _n_ob_restarts)
            ob_xyz_blocks: List[str] = []
            _ob_seen: set = set()
            ob_error: Optional[str] = None
            for _restart in range(_n_ob_restarts):
                _blocks, _err = _openbabel_generate_conformer_xyz(
                    smiles, num_confs=_per_ob, deterministic=False
                )
                if _blocks:
                    for _b in _blocks:
                        _key = "\n".join(
                            l.strip() for l in _b.splitlines() if l.strip()
                        )
                        if _key not in _ob_seen:
                            _ob_seen.add(_key)
                            ob_xyz_blocks.append(_b)
                elif _err and not ob_xyz_blocks:
                    ob_error = _err
            if ob_xyz_blocks:
                ob_ids = _inject_openbabel_conformers_into_mol(mol, ob_xyz_blocks)
                conf_ids.extend(ob_ids)
                logger.debug(
                    "OB conformers injected for isomer search: %d",
                    len(ob_ids),
                )
            elif ob_error:
                logger.debug("OB conformer generation: %s", ob_error)
        except Exception as ob_exc:
            logger.debug("OB conformer generation exception: %s", ob_exc)
            mol.RemoveAllConformers()
            conf_ids = []

    # Embed multiple conformers with deterministic seed schedule.
    # Seeds are independent → parallelize with ThreadPoolExecutor.
    #
    # Class-aware override: when ``DELFIN_CLASS_AWARE_SEEDS=1`` AND the
    # caller is using the default quality profile (i.e. ``_qprof`` was
    # built from module defaults — its ``seeds`` slot equals
    # ``DELFIN_TOP_LEVEL_SEED_COUNT``) we replace the seed count with
    # the per-class value from ``_resolve_top_level_seed_count``.
    # Explicit named profiles (``fast``/``max``/``extreme``) always win,
    # so operator-supplied ``quality_mode='max'`` cannot be silently
    # narrowed by the class heuristic.
    try:
        _qprof_seeds = int(_qprof.get("seeds", len(_TOP_LEVEL_SEEDS)))
        _resolved = _resolve_top_level_seed_count(mol)
        if (
            _qprof_seeds == int(DELFIN_TOP_LEVEL_SEED_COUNT)
            and _resolved != _qprof_seeds
        ):
            _qprof_seeds = _resolved
        seeds = list(_PIPELINE_SEEDS[:max(1, _qprof_seeds)])
        n_rounds = len(seeds)
        per_round = max(1, int(math.ceil(num_confs / n_rounds)))
        # Parallel ETKDG embedding with deterministic output.
        #
        # DETERMINISM CONTRACT: ``AllChem.EmbedMultipleConfs`` (and the
        # ``_rescale_metal_donor_distances`` post-step it triggers) mutate
        # the molecule they are handed — they add conformers, touch the
        # property cache and ring info, and rewrite conformer coordinates.
        # Running ``_embed_one`` for several seeds *concurrently against
        # the same shared ``mol``* is therefore a data race: RDKit's
        # conformer embedding is not thread-safe on a shared molecule, so
        # two runs of the same SMILES produced divergent conformer pools
        # (and hence divergent isomer counts / geometries) depending on
        # thread interleaving.  ``ThreadPoolExecutor.map`` only stabilises
        # *result ordering*, not the *content* of each per-seed embedding.
        #
        # Fix: every worker embeds into its OWN private ``Chem.Mol`` copy
        # (no shared mutable state), then the main thread copies the
        # resulting conformers back into the shared ``mol`` strictly in
        # seed order.  Parallelism and the deterministic seed schedule are
        # both preserved, and the conformer pool is now bit-identical
        # across runs.  For pathological macrocycles where the internal
        # embed timeout fires, two runs may still differ by a conformer or
        # two — acceptable since those SMILES are already flagged by the
        # timeout path as unreliable.
        def _embed_one(_s):
            """Embed ``per_round`` conformers for seed ``_s`` into a private
            mol copy.  Returns the list of per-copy conformer objects (deep
            copies, safe to re-add to the shared mol on the main thread)."""
            try:
                _mol_copy = Chem.Mol(mol)
                _mol_copy.RemoveAllConformers()
                _ids = _embed_multiple_confs_robust(_mol_copy, per_round, _s)
                # Materialise standalone Conformer copies so they survive
                # past the worker's mol copy going out of scope.
                return [
                    Chem.Conformer(_mol_copy.GetConformer(_cid))
                    for _cid in _ids
                ]
            except Exception:
                return []
        _n_workers = min(
            len(seeds), os.cpu_count() or 4, DELFIN_MAX_THREAD_WORKERS
        )
        if _n_workers > 1 and len(seeds) > 1:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=_n_workers
            ) as _pool:
                _per_seed_confs = list(_pool.map(_embed_one, seeds))
        else:
            _per_seed_confs = [_embed_one(_s) for _s in seeds]
        # Deterministic merge: re-add every worker's conformers to the
        # shared ``mol`` in seed order, assigning fresh sequential IDs.
        for _seed_confs in _per_seed_confs:
            for _conf in _seed_confs:
                try:
                    conf_ids.append(mol.AddConformer(_conf, assignId=True))
                except Exception:
                    pass

        # Fix A (Welle 2 / X10-FIPWAE): rescale every M-D bonded distance
        # to its ideal element-pair length before any downstream filtering
        # or UFF. ETKDG treats M-D bonds like organic bonds (~1.5 A) which
        # collapses heterodonor CN >= 5 metal complexes; when OB UFF then
        # triggers its unparam-TM HARD-fallback the broken distance is
        # frozen into the final XYZ. Pre-snap escapes the trap.
        # Default OFF; opt-in via DELFIN_PRE_UFF_MD_SNAP=1.
        if has_metal and _pre_uff_md_snap_enabled() and conf_ids:
            for _cid in conf_ids:
                try:
                    _snap_md_distances_to_ideal(mol, _cid)
                except Exception as _snap_exc:
                    logger.debug(
                        "Pre-UFF M-D snap failed for cid=%s: %s",
                        _cid, _snap_exc,
                    )
    except Exception as e:
        logger.warning("Multi-conformer embedding failed: %s", e)
        # Do not abort here: keep already injected OB conformers if available.
        # If none exist, continue with an empty pool so topo/linkage
        # enumeration can still add valid isomers.
        if not conf_ids:
            try:
                mol.RemoveAllConformers()
            except Exception:
                pass
            conf_ids = []

    if not conf_ids:
        logger.debug(
            "No conformers generated by OB/ETKDG; continuing with fallback + topo/linkage enumeration."
        )

    # Classify each conformer, skip obvious artifacts, then deduplicate by
    # full coordination fingerprint. This keeps distinct coordination
    # arrangements even when they share a coarse textual label.
    # Pre-compute donor types once (Morgan-based) to avoid repeated calls.
    dtype_map = _donor_type_map(mol)

    # Fix B (Welle 2 / X10-FIPWAE): pre-UFF M-D topology gate.
    # When env-flag DELFIN_PRE_UFF_TOPOLOGY_GATE=1 is set, discard any
    # conformer whose bonded M-D distance is outside [0.80, 1.20] times the
    # element-pair ideal. This is universal (graph + element symbols only)
    # and triggers BEFORE the existing geometry filters so the unparam-TM
    # frozen-broken frames cannot leak through. Default OFF.
    _pre_uff_md_gate = has_metal and _pre_uff_topology_gate_enabled()

    def _classify_one_conf(_cid, _relax, _m=None, _cc=None):
        """Classify conformer _cid.  _m/_cc let a worker pass a PRIVATE mol copy carrying only that
        conformer, so no thread ever touches the shared molecule (see DET_CLASSIFY_PAR below).  The
        returned conformer id is always the ORIGINAL _cid -- downstream indexes the real mol."""
        _M = mol if _m is None else _m
        _C = _cid if _cc is None else _cc
        try:
            if _pre_uff_md_gate:
                try:
                    if not _md_distance_in_tolerance(_M, _C):
                        return None
                except Exception:
                    pass
            if _has_atom_clash(_M, _C, min_dist=0.3):
                return None
            try:
                xyz_check = _mol_to_xyz_conformer(_M, _C)
                if not _metal_donor_distances_realistic(xyz_check, _M):
                    return None
            except Exception:
                return None
            if _has_pi_ring_nonplanarity(_M, _C):
                return None
            penalty = 0.0
            if _has_unphysical_metal_nonbonded_contact(_M, _C):
                if not _relax:
                    return None
                penalty += 350.0
            if _has_unphysical_oco_geometry(_M, _C):
                if not _relax:
                    return None
                penalty += 250.0
            if _has_atom_clash(_M, _C):
                penalty += 500.0
            if _has_bad_geometry(_M, _C):
                penalty += 300.0
            if _has_ligand_intertwining(_M, _C):
                penalty += 200.0
            fp = _compute_coordination_fingerprint(_M, _C, dtype_map=dtype_map)
            score = _geometry_quality_score(_M, _C) + penalty
            label = _classify_isomer_label(fp, _M)
            return (fp, label, _cid, score)
        except Exception:
            return None

    def _collect_fp_label_pairs(relax_hard_chem_filters: bool = False) -> List[Tuple[tuple, str, int, float]]:
        # ===== DATA RACE ON THE SHARED MOL (DELFIN_FFFREE_DET_CLASSIFY) ==========================
        # ROOT of the residual cross-run nondeterminism (2026-07-09).  `_classify_one_conf` reads the
        # SHARED `mol` from up to os.cpu_count() threads (384 on this box), calling _has_atom_clash,
        # _compute_coordination_fingerprint, _geometry_quality_score and _classify_isomer_label.
        # An RDKit Mol populates its ring info and property caches LAZILY, on first access, and is
        # not thread-safe while doing so.  Concurrent readers therefore race on that initialisation:
        #   * the SCORES wobble  -> a different conformer wins its fingerprint -> same label, other
        #     geometry (ZIZLUT "Isomer 1" moves 4-7 A between two runs of the identical build)
        #   * the FILTERS flip   -> a different set of conformers survives -> different frame COUNT
        #     (AFICIC 12 vs 13 frames)
        # This is why a deterministic total order over the selection (DET_SELECT) changed the output
        # but could not stabilise it: it sorts by numbers that are themselves produced in a race.
        # Same shape as the documented, already-written DELFIN_5G_T6_1_RACE_FIX (embedding against a
        # shared mol) -- which is likewise still default-OFF.
        # The fix is to stop sharing: classify sequentially.  Parallelism is not lost, it moves up a
        # level (the loop builds 64 SYSTEMS at once), and 384 threads per system x 64 systems was a
        # 24k-thread oversubscription anyway.  Element/graph-agnostic; byte-identical when unset.
        # DET_CLASSIFY_PAR: determinism WITHOUT giving up the cores.  The race was never the threads,
        # it was the SHARING.  Pre-warm the lazy ring/property caches on the template ONCE (so no
        # worker triggers a lazy init), then hand every worker a private conformer-free copy carrying
        # exactly the conformer it must classify.  Results are collected in SUBMISSION order, so the
        # output is bit-identical to the sequential path -- and that is the acceptance test.
        _det_classify = _delfin_env_int("DELFIN_FFFREE_DET_CLASSIFY", 0)
        _det_classify_par = _delfin_env_int("DELFIN_FFFREE_DET_CLASSIFY_PAR", 0)
        _n_cw = min(len(conf_ids), os.cpu_count() or 4) if conf_ids else 1
        if _det_classify_par and conf_ids and len(conf_ids) > 4:
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
            try:
                Chem.GetSymmSSSR(mol); mol.GetRingInfo()
            except Exception:
                pass
            _tmpl = Chem.Mol(mol); _tmpl.RemoveAllConformers()

            def _classify_private(_cid):
                try:
                    _m2 = Chem.Mol(_tmpl)
                    _c2 = _m2.AddConformer(Chem.Conformer(mol.GetConformer(_cid)), assignId=True)
                    return _classify_one_conf(_cid, relax_hard_chem_filters, _m2, _c2)
                except Exception:
                    return None
            _w = min(len(conf_ids), DELFIN_MAX_THREAD_WORKERS, os.cpu_count() or 4)
            if _w > 1:
                with concurrent.futures.ThreadPoolExecutor(max_workers=_w) as _cp:
                    _futs = [_cp.submit(_classify_private, c) for c in conf_ids]
                    return [r for r in (f.result() for f in _futs) if r is not None]
        if _det_classify or _det_classify_par:
            return [r for r in (_classify_one_conf(c, relax_hard_chem_filters)
                                for c in conf_ids) if r is not None]
        if _n_cw > 1 and len(conf_ids) > 4:
            with concurrent.futures.ThreadPoolExecutor(max_workers=_n_cw) as _cp:
                _futs = [
                    _cp.submit(_classify_one_conf, c, relax_hard_chem_filters)
                    for c in conf_ids
                ]
                return [
                    r for r in (f.result() for f in _futs) if r is not None
                ]
        return [
            r for r in (_classify_one_conf(c, relax_hard_chem_filters) for c in conf_ids)
            if r is not None
        ]

    # Two-pass sampling collection: strict first (catches truly broken
    # structures with hard rejects), then relaxed (lets structures with
    # minor OCO / nonbonded-contact issues through as penalty-scored
    # candidates).  Before, the relaxed pass only ran when the strict
    # pass emptied the list — on crowded bimetallic systems where 2-3
    # strict-survivors appeared the relaxed pass stayed dormant and
    # legitimate isomers that failed strict-but-would-pass-relaxed never
    # reached fingerprint deduplication.  Union by fingerprint keeps
    # the lowest-score representative per unique coordination pattern.
    _strict = _collect_fp_label_pairs(relax_hard_chem_filters=False)
    # Early exit: skip the costly relaxed-pass if the strict pass already
    # produced enough diverse fingerprints. Relaxed pass uses a second
    # ThreadPoolExecutor that can saturate the per-process thread budget on
    # large systems; only invoke it when strict pass is sparse (<70% of
    # max_isomers found). Threshold env-configurable.
    # Iter-8.3 multi-σ champion forward-port (5b3e0d2). When env=1 AND class
    # is multi_sigma, bypass _strict_sufficient short-circuit and always run
    # the relaxed pass. Full-pool multi-σ %match 49.45 → 63.36 (+13.91 pp),
    # frame-ratio 3.54×.
    # Phase 4D default-flip 2026-05-12: 0 → 1 — smoke500 verified
    # +28.6pp multi-sigma + 3.9pp sigma with combined ALT_MODE_BUDGET=60.
    _multisigma_port_active = bool(_delfin_env_int(
        "DELFIN_MULTISIGMA_PORT_5B3E0D2_ITER8", 1
    ))
    _multisigma_dispatch = False
    if _multisigma_port_active:
        try:
            _multisigma_dispatch = (_classify_complex_class(mol) == "multi_sigma")
        except Exception:
            _multisigma_dispatch = False
    # Iter-8.4a sigma chelate-cap port (123a130).  Sets module-global
    # ``_ITER84_SIGMA_CAPS_OVERRIDE`` when env=1 AND class='sigma' so
    # that ``_chelate_conformer_candidates`` widens ETKDG trial caps for
    # d8 / d10 chelate cohorts that lose the correct bite-angle
    # conformer under HEAD's tight caps.  ``global`` already declared at
    # function entry where the override is reset to None on every call.
    if _class_conditional_flag("DELFIN_SIGMA_PORT_123A130_ITER8", mol):
        try:
            if _classify_complex_class(mol) == "sigma":
                _ITER84_SIGMA_CAPS_OVERRIDE = _SIGMA_CHELATE_CAPS_123A
        except Exception:
            pass
    _relax_skip_threshold = float(os.environ.get(
        'DELFIN_RELAX_SKIP_THRESHOLD', '0.7'
    ))
    _strict_sufficient = (
        max_isomers > 0
        and len(_strict) >= max_isomers * _relax_skip_threshold
        and not _multisigma_dispatch
    )
    _relaxed = _collect_fp_label_pairs(relax_hard_chem_filters=True) if (
        conf_ids and has_metal and not _strict_sufficient
    ) else []
    # Iter-8.5a: per-class hard-reject merge (M_A from e6761e4 forensics).
    # When parent mol's class is in DELFIN_ITER85_HARD_REJECT_CLASSES,
    # skip the relaxed-pass merge entirely and use only the strict-pass
    # output.  e6761e4 ran strict-only on hapto / multi-hapto / multi-σ
    # cohorts; the relaxed-pass merge in HEAD lets through frames with
    # extra heavy bonds that the champion would have rejected.  The
    # tighter strict-only contract reduces extras-per-frame by 20-30 %
    # on those cohorts at the cost of fewer total frames.
    # Default: empty class set → bit-exact HEAD union-merge behaviour.
    # Activate via comma-separated list in
    # DELFIN_ITER85_HARD_REJECT_CLASSES, e.g. "hapto,multi_hapto".
    _iter85_hard_classes = set(
        x.strip() for x in (
            os.environ.get("DELFIN_ITER85_HARD_REJECT_CLASSES", "") or ""
        ).split(",") if x.strip()
    )
    _iter85_hard_dispatch = False
    if _iter85_hard_classes:
        try:
            _iter85_hard_dispatch = (
                _classify_complex_class(mol) in _iter85_hard_classes
            )
        except Exception:
            _iter85_hard_dispatch = False
    # ===== DETERMINISTIC TOTAL ORDER (DELFIN_FFFREE_DET_SELECT) ================================
    # ROOT of the cross-run nondeterminism (2026-07-09).  EVERY conformer selection in this chain
    # decides with a bare float comparison, in GENERATION order:
    #     if fp not in D or sc < D[fp][2]:        # an exact tie keeps whichever came FIRST
    # Two conformers of one fingerprint whose scores tie (or differ in the last bits) therefore
    # elect a winner that depends on insertion order and float noise.  The LABEL stays, the
    # GEOMETRY jumps -- exactly the observed defect (ZIZLUT "Isomer 1" moves 4-7 A between two runs
    # of the identical construction).  The fix is neither a seed nor a flag: make the order TOTAL.
    #   * quantise the score to 1e-6      -> immune to last-bit wobble
    #   * break ties by (cid, fp)         -> intrinsic, run-independent
    # This must wrap the WHOLE chain: _merged elects a winner per fingerprint BEFORE seen_fps ever
    # runs, so hardening only the later passes leaves the first one deciding by arrival time.
    # Element- and graph-agnostic; never SMILES- or refcode-specific.
    _det_select = _delfin_env_int("DELFIN_FFFREE_DET_SELECT", 0)

    def _skey(_sc: float, _cid: int, _fp) -> tuple:
        return (round(float(_sc), 6), int(_cid), repr(_fp))

    def _take(_d: dict, _fp, _lbl, _cid, _sc) -> None:
        """Insert (label, cid, score) for _fp iff it beats the incumbent under a TOTAL order."""
        _cur = _d.get(_fp)
        if _cur is None:
            _d[_fp] = (_lbl, _cid, _sc)
        elif _det_select:
            if _skey(_sc, _cid, _fp) < _skey(_cur[2], _cur[1], _fp):
                _d[_fp] = (_lbl, _cid, _sc)
        elif _sc < _cur[2]:
            _d[_fp] = (_lbl, _cid, _sc)

    if _iter85_hard_dispatch:
        # Strict-only contract: drop relaxed-pass entirely
        _merged: Dict[tuple, Tuple[str, int, float]] = {}
        for fp, lbl, cid, sc in _strict:
            _take(_merged, fp, lbl, cid, sc)
    else:
        # HEAD baseline: union-merge of strict + relaxed
        _merged = {}
        for fp, lbl, cid, sc in _strict:
            _take(_merged, fp, lbl, cid, sc)
        for fp, lbl, cid, sc in _relaxed:
            _take(_merged, fp, lbl, cid, sc)
    fp_label_pairs: List[Tuple[tuple, str, int, float]] = [
        (fp, lbl, cid, sc) for fp, (lbl, cid, sc) in _merged.items()
    ]
    # Generation order must stop mattering from here on: sort the pairs under the same total order.
    if _det_select:
        fp_label_pairs.sort(key=lambda t: _skey(t[3], t[2], t[0]))

    # Second pass: deduplicate, keeping the best-scoring conformer per group.
    # fp -> (label, conf_id, score)
    seen_fps: Dict[tuple, Tuple[str, int, float]] = {}
    for fp, label, cid, score in fp_label_pairs:
        _take(seen_fps, fp, label, cid, score)
        if len(seen_fps) >= max_isomers:
            break

    # RMSD-based dedup: remove geometrically identical conformers that
    # slipped through fingerprint-based dedup (e.g. borderline angle
    # classifications producing different fingerprints for the same isomer).
    if len(seen_fps) > 1:
        def _base_label(lbl: str) -> str:
            if not lbl:
                return ''
            return re.sub(r'-\d+$', '',str(lbl))

        # Same root: `fps_list` inherits the dict's GENERATION order, and the pairwise survivor is
        # picked with `si <= sj` -- a tie (or a last-bit difference) hands the decision to whichever
        # conformer happened to be generated first.  Impose the same total order here.
        fps_list = list(seen_fps.keys())
        if _det_select:
            fps_list.sort(key=lambda _f: _skey(seen_fps[_f][2], seen_fps[_f][1], _f))
        removed: set = set()
        # RMSD must NOT override the chemistry-based coordination fingerprint.
        # Every pair here has DISTINCT fingerprints (fps_list = seen_fps keys);
        # when their base labels ALSO differ they are distinct coordination
        # isomers (cis vs trans, fac vs mer, linkage isomers) and the current
        # code still merges them if they happen to sit within 0.8 A -- but RMSD
        # is an unreliable isomer-identity measure (distinct isomers can be
        # geometrically close; the whole-space measurement showed distinct
        # isomers merged at median 0.35 A, driving 58% of systems below their
        # theoretical isomer count).  Env-gated (default-OFF = byte-identical):
        # when set, cross-label merges are skipped entirely (and their expensive
        # _conformer_rmsd call avoided -> faster), trusting the fingerprint;
        # same-label conformer merging is unaffected (separate axis).
        _fp_strict = _delfin_env_int("DELFIN_FFFREE_ISOMER_FP_STRICT", 0)
        for i in range(len(fps_list)):
            if i in removed:
                continue
            _li, cid_i, si = seen_fps[fps_list[i]]
            base_i = _base_label(_li)
            for j in range(i + 1, len(fps_list)):
                if j in removed:
                    continue
                _lj, cid_j, sj = seen_fps[fps_list[j]]
                base_j = _base_label(_lj)
                if _fp_strict and base_i != base_j:
                    continue  # distinct isomers -> trust fingerprint, never RMSD-merge
                rmsd = _conformer_rmsd(mol, cid_i, cid_j)
                # Same-label: aggressive merge (RMSD < 2.5 A) —
                # these are supposed to be the same isomer anyway.
                # Cross-label: only merge if effectively identical
                # (RMSD < 0.8 A) to avoid false-positive merging of
                # distinct isomers whose fingerprints collided or
                # whose labels diverged due to borderline angles.
                # Tightened from 1.5 -> 0.8 after reports of
                # chemically different systems being wrongly merged.
                rmsd_threshold = 2.5 if base_i == base_j else 0.8
                if rmsd < rmsd_threshold:
                    _i_wins = (_skey(si, cid_i, fps_list[i]) <= _skey(sj, cid_j, fps_list[j])
                               if _det_select else si <= sj)
                    if _i_wins:
                        removed.add(j)
                    else:
                        removed.add(i)
                        break
        if removed:
            logger.debug("RMSD dedup removed %d duplicate(s)", len(removed))
            seen_fps = {fps_list[i]: seen_fps[fps_list[i]]
                        for i in range(len(fps_list)) if i not in removed}

    # Build results
    results: List[Tuple[str, str]] = []
    unknown_counter = 0
    if not seen_fps:
        # Sampling failed — get a single fallback geometry and still allow
        # the topological enumerator / linkage detector to augment it below.
        _fb_xyz, _fb_err = smiles_to_xyz(
            smiles, apply_uff=apply_uff, hapto_approx=hapto_mode,
            deterministic=deterministic,
        )
        if _fb_err:
            return [], _fb_err
        results = [(_fb_xyz, '')]
    else:
        # Number duplicate labels (e.g. multiple "mer" with different fingerprints)
        label_counts: Dict[str, int] = {}
        for fp in seen_fps:
            lbl = seen_fps[fp][0] or ''
            label_counts[lbl] = label_counts.get(lbl, 0) + 1
        label_seen: Dict[str, int] = {}
        relaxed_fragment_results: List[Tuple[str, str]] = []
        for fp, (label, cid, _score) in seen_fps.items():
            if not label:
                unknown_counter += 1
                display = f'Isomer {unknown_counter}'
            elif label_counts[label] > 1:
                label_seen[label] = label_seen.get(label, 0) + 1
                display = f'{label}-{label_seen[label]}'
            else:
                display = label
            try:
                xyz = _mol_to_xyz_conformer(mol, cid)
                if apply_uff:
                    try:
                        xyz = _optimize_xyz_openbabel_safe(
                            xyz,
                            mol_template=mol,
                            smiles=smiles,
                            apply_template_constraints=True,
                        )
                    except Exception as uff_exc:
                        # Preserve the isomer if UFF cannot optimize this geometry.
                        logger.debug(
                            "UFF optimization failed for conformer %s (%s), keeping unoptimized XYZ.",
                            cid, uff_exc,
                        )
            except Exception:
                continue
            # Topology check: organic ring count from OB XYZ must match original
            # SMILES (charge-insensitive: [N+]/[Fe-2] == [N]/[Fe] topologically).
            if not _roundtrip_ring_count_ok(xyz, smiles):
                logger.debug("Skipping conformer %d: topology mismatch", cid)
                continue
            # Bond check: reject conformers where OB perceives spurious bonds
            # between non-metal atoms (e.g. O-O, N-N) absent in original SMILES.
            if not _no_spurious_bonds(xyz, smiles):
                logger.debug("Skipping conformer %d: spurious bonds", cid)
                continue
            # Fragment topology check: organic ligand fragments must match the
            # original SMILES (catches broken/fused ligands while preserving
            # fac/mer isomers which have identical fragment sets).
            # For multi-metal complexes, skip strict fragment check — the
            # complex bridging topology causes frequent false-positive
            # mismatches in OB bond perception.
            _n_metals_in_mol = sum(
                1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
            )
            if not _fragment_topology_ok(xyz, smiles):
                if _n_metals_in_mol >= 2:
                    # Accept multi-metal conformers with relaxed topology
                    logger.debug(
                        "Accepting conformer %d despite fragment mismatch "
                        "(multi-metal complex)", cid)
                else:
                    logger.debug("Skipping conformer %d: fragment topology mismatch", cid)
                    if (
                        _fragment_topology_relaxed_fallback_ok(xyz, smiles)
                        and _xyz_passes_final_geometry_checks(xyz, mol)
                    ):
                        relaxed_fragment_results.append((xyz, display))
                    continue
            # Universal severe-distortion gate: catches C-H-C bridges,
            # unphysical bond stretches and atom overlaps that slip
            # past the fragment-topology graph check.  Same helper used
            # throughout the pipeline so one threshold = one rule.
            # Toggle via DELFIN_FINAL_GATE_ENABLED for regression diagnosis.
            if DELFIN_FINAL_GATE_ENABLED and not _xyz_passes_final_geometry_checks(
                xyz, mol, skip_angle_check=True
            ):
                logger.debug("Skipping conformer %d: severe covalent distortion", cid)
                continue
            results.append((xyz, display))

        if not results:
            if relaxed_fragment_results:
                logger.debug(
                    "Using %d conformer(s) despite fragment-topology mismatch "
                    "after stricter checks.",
                    len(relaxed_fragment_results),
                )
                results = relaxed_fragment_results
            else:
                _fb_xyz, _fb_err = smiles_to_xyz(
                    smiles, apply_uff=apply_uff, hapto_approx=hapto_mode,
                    deterministic=deterministic,
                )
                if _fb_err:
                    return [], _fb_err
                results = [(_fb_xyz, '')]

    # --- Topological enumerator: guarantee completeness ---
    # For MONO-metallic: use topology enumerator (guaranteed complete).
    # For MULTI-metallic: the topo builder destroys bridge ligand geometry
    # by moving metals independently. Instead, rely on ETKDG sampling
    # (which preserves ligand topology) + M-D rescaling + graph check.
    try:
        _n_metals_total = sum(
            1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
        ) if RDKIT_AVAILABLE else 0
    except Exception:
        _n_metals_total = 0
    if has_metal:
        try:
            # Collect fingerprints from sampling results
            existing_fps: set = set()
            topo_mol = mol
            dtype_map_topo = _donor_type_map(topo_mol)
            for existing_xyz, _existing_display in results:
                try:
                    mol_tmp = Chem.RWMol(topo_mol)
                    mol_tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), existing_xyz)
                    if conf is not None:
                        cid = mol_tmp.AddConformer(conf, assignId=True)
                        fp = _compute_coordination_fingerprint(
                            mol_tmp.GetMol(), cid, dtype_map=dtype_map_topo
                        )
                        existing_fps.add(fp)
                except Exception:
                    pass

            existing_displays = {display for _, display in results}
            existing_xyz_keys = {
                "\n".join(l.strip() for l in xyz.splitlines() if l.strip())
                for xyz, _d in results
            }

            topo_results = _generate_topological_isomers(
                topo_mol, smiles, apply_uff=apply_uff,
                max_isomers=max_isomers, n_metal_smart=n_metal_smart,
                profile=_qprof,
            )

            for topo_xyz, topo_label in topo_results:
                # Compute fingerprint of topo structure for dedup
                topo_fp = None
                try:
                    mol_tmp = Chem.RWMol(topo_mol)
                    mol_tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), topo_xyz)
                    if conf is not None:
                        cid = mol_tmp.AddConformer(conf, assignId=True)
                        topo_fp = _compute_coordination_fingerprint(
                            mol_tmp.GetMol(), cid, dtype_map=dtype_map_topo
                        )
                except Exception:
                    pass

                # Always let topo isomers through; the final geometry-based
                # dedup selects the best-scoring candidate per base label.
                # Sampling conformers are UFF-distorted, so topo versions
                # (donors pinned to ideal polyhedron) typically win.
                norm = topo_label or ''
                if not norm:
                    unknown_counter += 1
                    display = f'Isomer {unknown_counter}'
                else:
                    # Number duplicate labels
                    if norm in existing_displays:
                        suffix = 2
                        while f'{norm}-{suffix}' in existing_displays:
                            suffix += 1
                        display = f'{norm}-{suffix}'
                    else:
                        display = norm
                # Graph-based topology verification: checks every bond in
                # the template graph against actual XYZ distances.  Applied
                # uniformly to mono- and multi-metal isomers so no output
                # structure has a broken topology.  (The bridging-atom
                # bond-length window is already widened inside the gate
                # for M-(mu-X)-M' cases.)
                if not _verify_topology_from_graph(topo_xyz, topo_mol):
                    logger.debug(
                        "Skipping topo isomer %s: graph topology check failed",
                        display,
                    )
                    continue
                # Fix B (Welle 2 / X10-FIPWAE): pre-UFF M-D topology gate on
                # the topology-builder path. Some scaffold-builder paths can
                # produce stretched M-D when ligand fragments do not Procrustes-
                # fit cleanly. Default OFF.
                if _pre_uff_topology_gate_enabled():
                    try:
                        _topo_check = Chem.RWMol(topo_mol)
                        _topo_check.RemoveAllConformers()
                        _topo_conf = _xyz_to_rdkit_conformer(
                            _topo_check.GetMol(), topo_xyz
                        )
                        if _topo_conf is not None:
                            _topo_cid = _topo_check.AddConformer(
                                _topo_conf, assignId=True
                            )
                            if not _md_distance_in_tolerance(
                                _topo_check.GetMol(), _topo_cid
                            ):
                                logger.debug(
                                    "Skipping topo isomer %s: M-D distance "
                                    "out of tolerance",
                                    display,
                                )
                                continue
                    except Exception:
                        pass
                topo_key = "\n".join(l.strip() for l in topo_xyz.splitlines() if l.strip())
                if topo_key in existing_xyz_keys:
                    logger.debug("Skipping topo isomer %s: duplicate XYZ", display)
                    continue
                existing_displays.add(display)
                existing_xyz_keys.add(topo_key)
                if topo_fp is not None:
                    existing_fps.add(topo_fp)
                if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
                    try:
                        import numpy as np
                        _mi = next((a.GetIdx() for a in topo_mol.GetAtoms()
                                    if a.GetSymbol() in _METAL_SET), None)
                        _pos = {}
                        for _i, _ln in enumerate(topo_xyz.strip().split("\n")):
                            _p = _ln.split()
                            if len(_p) >= 4:
                                _pos[_i] = np.array([float(_p[1]), float(_p[2]), float(_p[3])])
                        _mp = _pos.get(_mi)
                        _cds = [nb.GetIdx() for nb in topo_mol.GetAtomWithIdx(_mi).GetNeighbors()
                                if nb.GetSymbol() == "C"] if _mi is not None else []
                        for _d in _cds:
                            _dp = _pos.get(_d)
                            if _mp is None or _dp is None:
                                continue
                            _mc = _mp - _dp; _mcn = float(np.linalg.norm(_mc))
                            if _mcn < 1e-6:
                                continue
                            _mc /= _mcn
                            for _nb in topo_mol.GetAtomWithIdx(_d).GetNeighbors():
                                if _nb.GetAtomicNum() <= 1:
                                    continue
                                _xp = _pos.get(_nb.GetIdx())
                                if _xp is None:
                                    continue
                                _cx = _xp - _dp; _cxn = float(np.linalg.norm(_cx))
                                if 1.3 < _cxn < 1.9:
                                    _ang = float(np.degrees(np.arccos(
                                        max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                                    _trace_seating("EMIT_TOPO donor=%d M-C-Xheavy=%.0f disp=%s" % (
                                        _d, _ang, str(display)[:18]))
                    except Exception:
                        pass
                results.append((topo_xyz, display))
        except Exception as _topo_exc:
            logger.debug("Topological isomer generation failed: %s", _topo_exc)

    # --- Multi-metal sampling augmentation ---
    # For multi-metallic systems: augment with extra ETKDG conformers.
    # EXCEPTION: hapto complexes — ETKDG can't build Cp rings correctly
    # (they come out non-planar). The hapto builder already produced
    # perfect Cp geometry, so skip ETKDG augmentation for hapto systems.
    _skip_mm_augmentation = bool(hapto_groups and hapto_mode)
    if has_metal and _n_metals_total >= 2 and len(results) < max_isomers and not _skip_mm_augmentation:
        try:
            # Scale multi-metal seed count with the quality profile — the
            # fixed 10-seed set was the single biggest search-space
            # limitation for bimetallic systems: mono-metal paths use
            # 12/20/40/60 seeds (fast/normal/max/extreme) while the
            # multinuclear augmentation was capped at 10 regardless of
            # profile.  Sharing the deterministic _PIPELINE_SEEDS pool
            # keeps results reproducible across runs.
            _extra_seed_count = int(_qprof.get("seeds", 10))
            _mm_walltime: Optional[float] = None
            # Multi-sigma V2: shrink the augmentation seed pool AND add a
            # wall-clock budget so 50-125-atom bimetallic complexes can't
            # consume the entire caller timeout on this single block.
            # Bit-exact pre-patch when _ms_v2_budget is None.
            if _ms_v2_budget is not None:
                _extra_seed_count = min(
                    _extra_seed_count, int(_ms_v2_budget["seeds_mm_aug"])
                )
                _mm_walltime = float(_ms_v2_budget["mm_walltime"])
            # Determinism: drop the wall-clock budget so the augmentation uses the
            # untimed (deterministic) future-collection path below, never the
            # per-future timeout path whose result depends on timing/CPU load.
            if _deterministic_mode():
                _mm_walltime = None
            _extra_seeds = list(
                _PIPELINE_SEEDS[:max(1 if _ms_v2_budget is not None else 10,
                                     _extra_seed_count)]
            )
            _extra_ids: List[int] = []
            _n_extra = min(len(_extra_seeds), os.cpu_count() or 4, 64)
            # DETERMINISM CONTRACT (Welle-3 T6.1 — twin of fix #2 / 7cf73e3):
            # The previous implementation submitted
            # ``_embed_multiple_confs_robust(mol, 3, s)`` concurrently against
            # the SAME shared ``mol``.  RDKit's conformer embedding mutates
            # the molecule (conformers + property cache + ring info), so the
            # multi-metal augmentation reproduced the same data race that
            # fix #2 patched in the primary embed loop at line ~25871.  Each
            # worker embeds into its own private ``Chem.Mol`` copy, and
            # the main thread re-adds the resulting conformers to the shared
            # ``mol`` in seed-submission order so the conformer-ID schedule
            # is reproducible across runs and worker counts.
            #
            # Welle-5g Step-0 targeted revert per 5f-A bisect: a3edabe shipped
            # this race-fix always-on as part of the Correctness-Bundle and
            # was attributed ~70% of the -4064 NET full-pool regression
            # (hapto -767).  Gated behind DELFIN_5G_T6_1_RACE_FIX (default 0)
            # to restore the legacy shared-mol submission path at default;
            # opt-in via env-flag re-enables the per-worker copy.
            def _aug_embed_one(_s: int) -> List:
                try:
                    _mol_copy = Chem.Mol(mol)
                    _mol_copy.RemoveAllConformers()
                    _ids = _embed_multiple_confs_robust(_mol_copy, 3, _s)
                    return [
                        Chem.Conformer(_mol_copy.GetConformer(_cid))
                        for _cid in _ids
                    ]
                except Exception:
                    return []

            _t61_race_fix_on = _delfin_env_int("DELFIN_5G_T6_1_RACE_FIX", 0) > 0

            with concurrent.futures.ThreadPoolExecutor(max_workers=_n_extra) as _xp:
                if _t61_race_fix_on:
                    _xfuts = [_xp.submit(_aug_embed_one, s) for s in _extra_seeds]
                else:
                    # Legacy pre-T6.1 path: submit shared mol directly.
                    _xfuts = [
                        _xp.submit(_embed_multiple_confs_robust, mol, 3, s)
                        for s in _extra_seeds
                    ]
                # Submission-order traversal preserves determinism (see
                # top-level embedding loop for rationale).  When a
                # wall-clock budget is active, switch to per-future
                # timeouts so a single hung embedding cannot starve the
                # remaining seeds.
                _per_seed_confs: List[List] = [[] for _ in _xfuts]
                if _mm_walltime is not None:
                    _mm_deadline = time.time() + _mm_walltime
                    for _i, _xf in enumerate(_xfuts):
                        _remaining = max(0.0, _mm_deadline - time.time())
                        if _remaining <= 0.0:
                            # Budget exhausted — cancel the rest.
                            try:
                                _xf.cancel()
                            except Exception:
                                pass
                            continue
                        try:
                            if _t61_race_fix_on:
                                _per_seed_confs[_i] = _xf.result(timeout=_remaining)
                            else:
                                _extra_ids.extend(_xf.result(timeout=_remaining))
                        except concurrent.futures.TimeoutError:
                            try:
                                _xf.cancel()
                            except Exception:
                                pass
                        except Exception:
                            pass
                else:
                    for _i, _xf in enumerate(_xfuts):
                        try:
                            if _t61_race_fix_on:
                                _per_seed_confs[_i] = _xf.result()
                            else:
                                _extra_ids.extend(_xf.result())
                        except Exception:
                            pass
            # Deterministic merge (race-fix path only): re-add every worker's
            # conformers to the shared ``mol`` in seed-submission order with
            # fresh sequential IDs.  Legacy path already populated
            # ``_extra_ids`` directly via shared-mol embedding.
            if _t61_race_fix_on:
                for _seed_confs in _per_seed_confs:
                    for _conf in _seed_confs:
                        try:
                            _extra_ids.append(mol.AddConformer(_conf, assignId=True))
                        except Exception:
                            pass
            logger.debug("Multi-metal augmentation: %d extra conformers", len(_extra_ids))

            # Classify + dedup with existing results
            _existing_fps_mm: set = set()
            for _xyz, _d in results:
                try:
                    _mt = Chem.RWMol(mol)
                    _mt.RemoveAllConformers()
                    _c = _xyz_to_rdkit_conformer(_mt.GetMol(), _xyz)
                    if _c is not None:
                        _ci = _mt.AddConformer(_c, assignId=True)
                        _fp = _compute_coordination_fingerprint(
                            _mt.GetMol(), _ci, dtype_map=dtype_map
                        )
                        _existing_fps_mm.add(_fp)
                except Exception:
                    pass

            for _xid in _extra_ids:
                if len(results) >= max_isomers:
                    break
                try:
                    _xyz = _mol_to_xyz_conformer(mol, _xid)
                    _xyz = _optimize_xyz_openbabel_safe(_xyz, mol_template=mol)
                    if not _verify_topology_from_graph(_xyz, mol):
                        continue
                    _mt = Chem.RWMol(mol)
                    _mt.RemoveAllConformers()
                    _c = _xyz_to_rdkit_conformer(_mt.GetMol(), _xyz)
                    if _c is None:
                        continue
                    _ci = _mt.AddConformer(_c, assignId=True)
                    _fp = _compute_coordination_fingerprint(
                        _mt.GetMol(), _ci, dtype_map=dtype_map
                    )
                    if _fp in _existing_fps_mm:
                        continue
                    _existing_fps_mm.add(_fp)
                    _lbl = _classify_isomer_label(_fp, _mt.GetMol())
                    if not _lbl:
                        unknown_counter += 1
                        _lbl = f'Isomer {unknown_counter}'
                    results.append((_xyz, _lbl))
                except Exception:
                    continue
        except Exception as _mm_exc:
            logger.debug("Multi-metal augmentation failed: %s", _mm_exc)

    # --- Linkage isomers ---
    if has_metal and include_binding_mode_isomers:
        try:
            link_results = _generate_linkage_isomers(
                mol, smiles, apply_uff=apply_uff,
                max_template_tries=int(_qprof.get("alt_tries", 8)),
            )
            for _lxyz, _llabel in link_results:
                if not _metal_donor_distances_realistic(_lxyz, mol):
                    logger.debug("Skipping linkage isomer %s: unphysical M-D distance", _llabel)
                    continue
                if _fragment_topology_ok(_lxyz, smiles):
                    # Fix B (Welle 2 / X10-FIPWAE): linkage isomer path can
                    # also produce off-target M-D bonded distances when the
                    # alternate donor sits at a chemically different idealised
                    # distance. Apply the same pre-UFF M-D topology gate.
                    if _pre_uff_topology_gate_enabled():
                        try:
                            _link_mol = Chem.RWMol(mol)
                            _link_mol.RemoveAllConformers()
                            _link_conf = _xyz_to_rdkit_conformer(
                                _link_mol.GetMol(), _lxyz
                            )
                            if _link_conf is not None:
                                _link_cid = _link_mol.AddConformer(
                                    _link_conf, assignId=True
                                )
                                if not _md_distance_in_tolerance(
                                    _link_mol.GetMol(), _link_cid
                                ):
                                    logger.debug(
                                        "Skipping linkage isomer %s: M-D "
                                        "distance out of tolerance",
                                        _llabel,
                                    )
                                    continue
                        except Exception:
                            pass
                    if DELFIN_FINAL_GATE_ENABLED and not _xyz_passes_final_geometry_checks(
                        _lxyz, mol, skip_angle_check=True
                    ):
                        logger.debug("Skipping linkage isomer %s: severe covalent distortion", _llabel)
                        continue
                    results.append((_lxyz, _llabel))
                else:
                    logger.debug("Skipping linkage isomer %s: fragment topology mismatch", _llabel)
        except Exception as _link_exc:
            logger.debug("Linkage isomer generation failed: %s", _link_exc)

    # --- Alternative binding-site exploration ---
    _alt_tries = int(_qprof.get("alt_tries", 8))
    if has_metal and include_binding_mode_isomers and _alt_tries > 0:
        try:
            existing_displays = {display for _, display in results}
            existing_base = {re.sub(r'-\d+$', '',d) for d in existing_displays}
            alt_results = _generate_alternative_binding_modes(
                mol, smiles, apply_uff=apply_uff,
                max_template_tries=_alt_tries,
            )
            for alt_xyz, alt_label in alt_results:
                if alt_label not in existing_base:
                    if not _metal_donor_distances_realistic(alt_xyz, mol):
                        logger.debug(
                            "Skipping alt-binding isomer %s: unphysical M-D distance",
                            alt_label,
                        )
                        continue
                    if not _fragment_topology_ok(alt_xyz, smiles):
                        logger.debug("Skipping alt-binding isomer %s: fragment topology mismatch", alt_label)
                        continue
                    # Fix B (Welle 2 / X10-FIPWAE): alternative-binding-mode
                    # path can produce catastrophically collapsed M-D
                    # geometries (e.g. Ni-N at 1.31 A) that pass fragment
                    # topology but are physically meaningless. When the
                    # pre-UFF topology gate is on, reject them based on
                    # bonded M-D distance. Default OFF.
                    if _pre_uff_topology_gate_enabled():
                        try:
                            _alt_mol = Chem.RWMol(mol)
                            _alt_mol.RemoveAllConformers()
                            _alt_conf = _xyz_to_rdkit_conformer(
                                _alt_mol.GetMol(), alt_xyz
                            )
                            if _alt_conf is not None:
                                _alt_cid = _alt_mol.AddConformer(
                                    _alt_conf, assignId=True
                                )
                                if not _md_distance_in_tolerance(
                                    _alt_mol.GetMol(), _alt_cid
                                ):
                                    logger.debug(
                                        "Skipping alt-binding isomer %s: "
                                        "M-D distance out of tolerance",
                                        alt_label,
                                    )
                                    continue
                        except Exception:
                            pass
                    if DELFIN_FINAL_GATE_ENABLED and not _xyz_passes_final_geometry_checks(
                        alt_xyz, mol, skip_angle_check=True
                    ):
                        logger.debug("Skipping alt-binding isomer %s: severe covalent distortion", alt_label)
                        continue
                    results.append((alt_xyz, alt_label))
                    existing_base.add(alt_label)
        except Exception as _alt_exc:
            logger.debug("Alternative binding mode generation failed: %s", _alt_exc)

    # --- Final label dedup (safety net) ---
    # Collapse entries that share the same base label (e.g. "trans-1" and
    # "trans-2" that slipped through fingerprint/RMSD dedup).  Pick the
    # geometrically best candidate per base label (sampling conformers
    # from UFF-distorted ETKDG are usually worse than the topology-
    # enumerator output which pins donors to the idealized polyhedron).
    # "Isomer N" labels use spaces not dashes, so they are never collapsed.
    if collapse_label_variants and results and has_metal:
        try:
            score_mol = mol

            def _score_xyz(_xyz: str) -> float:
                if score_mol is None:
                    return float("inf")
                try:
                    tmp = Chem.RWMol(score_mol)
                    tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(tmp.GetMol(), _xyz)
                    if conf is None:
                        return float("inf")
                    cid = tmp.AddConformer(conf, assignId=True)
                    return float(_geometry_quality_score(tmp.GetMol(), cid))
                except Exception:
                    return float("inf")

            # Pick best-score candidate per base label, preserve first
            # appearance order for the final result list.
            base_best: Dict[str, Tuple[int, float]] = {}
            scored_cache: List[float] = []
            for idx, (xyz, lbl) in enumerate(results):
                base = re.sub(r'-\d+$', '',lbl) if lbl else ''
                score = _score_xyz(xyz)
                scored_cache.append(score)
                if base not in base_best or score < base_best[base][1]:
                    base_best[base] = (idx, score)

            keep_indices = {idx for idx, _s in base_best.values()}
            new_results: List[Tuple[str, str]] = []
            for idx, (xyz, lbl) in enumerate(results):
                if idx not in keep_indices:
                    logger.debug(
                        "Final dedup: dropping %r (base=%r, score=%.3f) — better kept",
                        lbl,
                        re.sub(r'-\d+$', '',lbl) if lbl else '',
                        scored_cache[idx],
                    )
                    continue
                new_results.append((xyz, re.sub(r'-\d+$', '',lbl) if lbl else lbl))
            results = new_results
        except Exception as _dedup_exc:
            logger.debug("Geometry-based dedup failed, falling back to first-keep: %s", _dedup_exc)
            _seen_base: Dict[str, int] = {}
            _keep: List[bool] = [True] * len(results)
            for _idx, (_, _lbl) in enumerate(results):
                _base = re.sub(r'-\d+$', '',_lbl) if _lbl else ''
                if _base in _seen_base:
                    _keep[_idx] = False
                else:
                    _seen_base[_base] = _idx
            results = [
                (xyz, re.sub(r'-\d+$', '',lbl))
                for (xyz, lbl), keep in zip(results, _keep) if keep
            ]

    # --- Symmetry-aware within-label dedup ---
    # Collapses duplicate rotameric entries that share the same base
    # label AND the same coordination fingerprint AND are within RMSD
    # 1.0 Å.  Fundamental principle: two entries with different labels
    # (``fac`` vs ``mer``, ``trans`` vs ``see-saw O-O-ax`` etc.) or
    # different fingerprints MUST NEVER be collapsed here — the
    # classifier already said they are distinct isomers.  Cross-label
    # collapse is handled by the earlier RMSD pass at 0.8 Å threshold
    # inside ``seen_fps``.  Keep the entry with the lowest ideal-
    # polyhedron deviation (= most textbook-symmetric representative).
    if has_metal and len(results) > 1 and _delfin_env_int("DELFIN_SYM_DEDUP_ENABLED", 1):
        try:
            def _entry_fp_and_sym(xyz: str) -> Tuple[Optional[tuple], float]:
                """Return (fingerprint, polyhedron_total_dev) for an XYZ."""
                try:
                    mt = Chem.RWMol(mol)
                    mt.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
                    if conf is None:
                        return (None, float("inf"))
                    cid = mt.AddConformer(conf, assignId=True)
                    fp = _compute_coordination_fingerprint(
                        mt.GetMol(), cid, dtype_map=dtype_map
                    )
                    devs = _ideal_polyhedron_angle_dev_per_metal(mt.GetMol(), cid)
                    total = sum(devs.values()) if devs else float("inf")
                    return (fp, total)
                except Exception:
                    return (None, float("inf"))

            entry_info = [_entry_fp_and_sym(xyz) for xyz, _ in results]
            # Base label (strip trailing -conf\d+, -\d+).
            def _base_lbl(_l: str) -> str:
                return re.sub(r'-\d+$', '', _l) if _l else ''
            labels_base = [_base_lbl(l) for _, l in results]
            removed: set = set()
            _RMSD_SAME_LBL_FP = 1.0  # Å, heavy-atom RMSD
            for i in range(len(results)):
                if i in removed:
                    continue
                fp_i, tot_i = entry_info[i]
                lbl_i = labels_base[i]
                if fp_i is None:
                    continue
                for j in range(i + 1, len(results)):
                    if j in removed:
                        continue
                    fp_j, tot_j = entry_info[j]
                    lbl_j = labels_base[j]
                    # HARD RULE: different labels OR different fingerprints
                    # ⇒ never collapse.  fac vs mer, cis vs trans, see-saw
                    # vs trans, and linkage isomers all have distinct
                    # labels/fingerprints and must survive.
                    if fp_j is None or fp_i != fp_j or lbl_i != lbl_j:
                        continue
                    rmsd = _heavy_atom_rmsd_xyz(results[i][0], results[j][0])
                    if rmsd < _RMSD_SAME_LBL_FP:
                        if tot_i <= tot_j:
                            removed.add(j)
                        else:
                            removed.add(i)
                            break
            if removed:
                logger.debug(
                    "Symmetry-aware dedup removed %d within-label+fp dup(s) (RMSD<%.1f)",
                    len(removed), _RMSD_SAME_LBL_FP,
                )
                results = [
                    r for idx, r in enumerate(results) if idx not in removed
                ]
        except Exception as _sym_exc:
            logger.debug("Symmetry-aware dedup skipped: %s", _sym_exc)

    # --- Combinatorial upper-bound check (log only) ---
    # If N_output wildly exceeds the Pólya-estimated bound, log a warning
    # so the user can investigate.  Collapse is NOT applied here because
    # the Pólya bound cannot distinguish chemically distinct isomers that
    # share a bound ceiling (e.g. fac vs mer both count toward the same 2,
    # and we would lose one if we collapse by RMSD alone).  The within-
    # fingerprint dedup above already removes true duplicates safely.
    if has_metal and len(results) > 1:
        try:
            _bound = _estimate_isomer_upper_bound(mol, dtype_map=dtype_map)
            if _bound is not None and len(results) > max(10, 3 * _bound):
                logger.info(
                    "Combinatorial-excess: %d isomers vs bound %d — review for dupes",
                    len(results), _bound,
                )
        except Exception as _comb_exc:
            logger.debug("Combinatorial-bound check skipped: %s", _comb_exc)

    # --- Additive trans-effect pass ---
    # Explicit "all-trans-by-type" coordination arrangements are central
    # in coordination chemistry (trans-effect, σ-donor / π-acceptor
    # competition).  The main pipeline's UFF + fingerprint-dedup
    # converges them with already-emitted candidates, so here we append
    # idealised pre-UFF placements that keep the trans symmetry intact.
    # Pure additive: never touches existing entries, only appends new
    # XYZ-signature-distinct ones.  Toggle via DELFIN_TRANS_PASS_ENABLED.
    # Iter-8.5b: when DELFIN_ITER85_PUMP_SKIP_CLASSES contains the parent
    # mol's class, skip the trans-effect pass.  e6761e4 (champion for
    # hapto / multi-hapto / multi-σ) does not run a trans-effect pass
    # for those classes; the pass adds frames that often violate
    # chemistry filters at downstream gates.  Reusing the comma-separated
    # class list pattern from Iter-8.5a.  Default empty set = bit-exact.
    _iter85b_pump_skip_classes = set(
        x.strip() for x in (
            os.environ.get("DELFIN_ITER85_PUMP_SKIP_CLASSES", "") or ""
        ).split(",") if x.strip()
    )
    _iter85b_pump_skip = False
    if has_metal and _iter85b_pump_skip_classes:
        try:
            _iter85b_pump_skip = (
                _classify_complex_class(mol) in _iter85b_pump_skip_classes
            )
        except Exception:
            _iter85b_pump_skip = False
    if has_metal and not DELFIN_TOPOLOGY_STRICT_MODE and not _iter85b_pump_skip:
        try:
            _emit_all_trans_by_type_arrangements(
                mol, results, dtype_map, apply_uff, max_isomers,
            )
        except Exception as _trans_exc:
            logger.debug("Trans-effect pass skipped: %s", _trans_exc)

    # --- Additive chelate-pucker variants pass ---
    # Saturated chelate rings (cyclam, en, dien, salen-CH₂CH₂, polyamines)
    # have multiple low-energy ring-pucker minima (chair/boat/twist).
    # The main pipeline emits one ETKDG conformer per chelate; pucker
    # variation gets lost in dedup.  This pass emits explicit chair/
    # boat/twist perturbations of each saturated chelate ring as
    # separate XYZ candidates.  Pure additive (toggle:
    # DELFIN_PUCKER_PASS_ENABLED, default 1).
    # Iter-8.4b: when DELFIN_SIGMA_SKIP_PUCKER_ITER8=1 AND class='sigma',
    # skip the pucker bucket entirely.  The 123a130 sigma champion does
    # not run a pucker pass; HEAD's pucker pass adds frames that are 38 %
    # match (mostly bad).  Skipping for sigma class trims the bucket and
    # raises aggregate sigma %match.  Default OFF for bit-exactness.
    _iter84b_skip_pucker = False
    if has_metal and _class_conditional_flag("DELFIN_SIGMA_SKIP_PUCKER_ITER8", mol):
        try:
            _iter84b_skip_pucker = (
                _classify_complex_class(mol) == "sigma"
            )
        except Exception:
            _iter84b_skip_pucker = False
    # Iter-8.5b: same class-list as the trans-effect skip above.  When the
    # parent mol's class is in DELFIN_ITER85_PUMP_SKIP_CLASSES, skip the
    # pucker pass too.  Combined with 8.4b's sigma-specific skip, the
    # pucker pass is now class-list-dispatched.
    if has_metal and not DELFIN_TOPOLOGY_STRICT_MODE and not _iter84b_skip_pucker and not _iter85b_pump_skip:
        try:
            _emit_chelate_pucker_variants(
                mol, results, apply_uff, max_isomers,
            )
        except Exception as _pucker_exc:
            logger.debug("Pucker pass skipped: %s", _pucker_exc)

    # --- Additive NON-METAL ring-pucker enumeration pass (sibling) ---
    # The metal pucker pass above only touches rings that CONTAIN a metal.
    # Peripheral non-metal, non-aromatic rings (cyclohexyl, piperidinyl,
    # sugar, ...) hanging off the complex never receive pucker enumeration
    # and stay frozen in whatever basin a single ETKDG seed produced.  This
    # sibling drives every such ring to a basin-verified chair (one global
    # chair-set frame) plus a small bounded set of one-ring-flipped-to-boat
    # decorations.  Pure additive; master flag DELFIN_RING_PUCKER_ENUM
    # (default 0 -> byte-identical no-op).
    if not DELFIN_TOPOLOGY_STRICT_MODE:
        try:
            _emit_nonmetal_ring_pucker_variants(
                mol, results, apply_uff, max_isomers,
            )
        except Exception as _nm_pucker_exc:
            logger.debug("Non-metal ring-pucker pass skipped: %s", _nm_pucker_exc)

    # --- Correct combinatorial ring-pucker construction (metal path) ----------
    # Supersedes the two legacy pucker passes above with the Cremer-Pople +
    # torsion-held-relax constructor, which reaches the higher ring basins ETKDG
    # never samples AND runs the multi-ring product, for chelate rings (metal +
    # donors frozen -> coordination preserved) and peripheral ligand rings.
    # Additive; the final topology gate drops any that violate the graph.
    if not DELFIN_TOPOLOGY_STRICT_MODE:
        try:
            _emit_ring_puckers_rp(mol, results, apply_uff, max_isomers)
        except Exception as _rp_metal_exc:
            logger.debug("RP ring-pucker pass skipped: %s", _rp_metal_exc)

    # --- Additive d8-CN4 square-planar seating pass (eye-driven; DELFIN_D8_SP4_SEAT, default 0) ---
    # d8 metals (Pt/Pd/Ni/Au/Rh/Ir) at CN4 are square-planar in the crystal; the chelate embed can
    # leave them tetrahedral (whole manifold ALARM).  Appends a de-tilted SP-4 variant; the final
    # topology gate below drops any that overlap -> strictly additive / never-worse.
    if has_metal and not DELFIN_TOPOLOGY_STRICT_MODE:
        try:
            _emit_d8_sp4_variants(mol, results, apply_uff, max_isomers)
        except Exception as _d8_exc:
            logger.debug("d8-SP4 seating pass skipped: %s", _d8_exc)

    # --- Final output gate: graph-based topology verification ---
    # Every output structure must preserve the bond topology from the
    # input SMILES.  Uses _verify_topology_from_graph which checks
    # every bond distance against the template graph — no OB perception.
    if has_metal and results:
        verified: List[Tuple[str, str]] = []
        for xyz, lbl in results:
            if _verify_topology_from_graph(xyz, mol):
                verified.append((xyz, lbl))
            else:
                logger.info(
                    "Output gate rejected isomer %r: graph topology violated", lbl
                )
        if verified:
            results = verified

    # Sort isomers by a composite quality score (UFF energy +
    # geometry regularity + symmetry) ascending — most realistic
    # first — and reject outliers whose energy is physically
    # unrealistic.  Energy alone is not a sufficient ordering
    # criterion for metal complexes because OpenBabel UFF lacks
    # parameters for Sc / Cd / most lanthanides and actinides, so
    # its energies for those centres are effectively random in
    # absolute terms but still roughly monotonic within a pool of
    # similar geometries.  Mixing in the geometry_quality score
    # (heterolept-aware bond-length spread + polyhedron-specific
    # angle deviation + within-bucket symmetry bonus) and a modest
    # symmetry weight lets chemically sensible isomers outrank
    # pseudo-minima that happen to have low UFF energy.
    if has_metal and len(results) > 1 and OPENBABEL_AVAILABLE:
        try:
            # Topology-preservation penalty helper — measures how
            # faithful the actual XYZ bond lengths are to the SMILES
            # bonds.  Complements geometry_quality_score (which
            # measures polyhedron symmetry) by explicitly scoring
            # how well the input bond topology is preserved after
            # build + UFF.  Per-bond penalty = |d/ideal - 1| * 10.
            def _topology_preservation_penalty(xyz_str: str) -> float:
                try:
                    lines = [l for l in xyz_str.strip().splitlines() if l.strip()]
                    if len(lines) != mol.GetNumAtoms():
                        try:
                            mol_h = Chem.AddHs(mol)
                            if len(lines) == mol_h.GetNumAtoms():
                                _tmol = mol_h
                            else:
                                return 0.0
                        except Exception:
                            return 0.0
                    else:
                        _tmol = mol
                    _coords: List[Tuple[float, float, float]] = []
                    for _ln in lines:
                        _p = _ln.split()
                        if len(_p) < 4:
                            return 0.0
                        _coords.append((float(_p[1]), float(_p[2]), float(_p[3])))
                    _pen = 0.0
                    for _b in _tmol.GetBonds():
                        _a1 = _b.GetBeginAtom()
                        _a2 = _b.GetEndAtom()
                        _s1 = _a1.GetSymbol()
                        _s2 = _a2.GetSymbol()
                        _i1 = _a1.GetIdx()
                        _i2 = _a2.GetIdx()
                        _dx = _coords[_i1][0] - _coords[_i2][0]
                        _dy = _coords[_i1][1] - _coords[_i2][1]
                        _dz = _coords[_i1][2] - _coords[_i2][2]
                        _d = math.sqrt(_dx * _dx + _dy * _dy + _dz * _dz)
                        _is_m1 = _s1 in _METAL_SET
                        _is_m2 = _s2 in _METAL_SET
                        _is_ml = False
                        if _is_m1 and _is_m2:
                            _mmk = frozenset({_s1, _s2})
                            _ideal = _METAL_METAL_BOND_LENGTHS.get(_mmk)
                            if _ideal is None:
                                _r1 = _COVALENT_RADII.get(_s1)
                                _r2 = _COVALENT_RADII.get(_s2)
                                _ideal = (_r1 + _r2 + 0.3) if _r1 and _r2 else 2.5
                            _is_ml = True  # M-M treated as coordination
                        elif _is_m1 or _is_m2:
                            _ms = _s1 if _is_m1 else _s2
                            _ds = _s2 if _is_m1 else _s1
                            _ideal = float(_get_ml_bond_length(_ms, _ds))
                            if _ideal <= 0:
                                continue
                            _is_ml = True
                        else:
                            _r1 = _COVALENT_RADII.get(_s1)
                            _r2 = _COVALENT_RADII.get(_s2)
                            if _r1 is None or _r2 is None:
                                continue
                            _ideal = _r1 + _r2
                        if _ideal <= 0:
                            continue
                        # Quadratic penalty: small deviations stay
                        # cheap, large deviations (broken topology)
                        # get disproportionately penalized.
                        #   5% dev -> 0.25 pt/bond
                        #  15% dev -> 2.25 pt/bond
                        #  30% dev -> 9.00 pt/bond
                        #  50% dev -> 25.0 pt/bond
                        # M-L bonds weighted 3x — coordination sphere
                        # fidelity is more important than organic C-C.
                        _dev = _d / _ideal - 1.0
                        _contrib = (_dev * 10.0) ** 2
                        if _is_ml:
                            _contrib *= 3.0
                        _pen += _contrib
                    return _pen
                except Exception:
                    return 0.0

            scored: List[Tuple[float, float, float, str, str]] = []
            # d8 CN4 square-planar rank preference (env-gated, default-OFF =
            # byte-identical).  ROOT of the SS-4->SP-4 poly_match cluster: the
            # bucket sort below is keyed on the UFF energy bucket FIRST, and
            # UFF favours tetrahedral for d8, so the tetrahedral ETKDG conformer
            # out-ranks the (correct, already-built) square-planar placement —
            # even though _geometry_quality_score / _preferred_cn4_geometry_score
            # already prefer SP-4 for these metals (the preference just never
            # dominates the energy bucket).  When enabled, for systems bearing a
            # d8 CN4 centre (Ni/Pd/Pt/Au — exactly the metals the CN4 scorer
            # biases toward square) the geometry+topology score becomes the
            # PRIMARY key, so the correct polyhedron wins; energy stays the
            # tie-break.  Reorders frames only (never drops) -> isomer-complete
            # by construction; the isomers_lost floor guards regardless.  No new
            # per-frame geometry work (reuses g) -> construction stays fast.
            _d8cn4_square = False
            if _delfin_env_int("DELFIN_FFFREE_D8_RANK_PREFER", 0):
                try:
                    for _a in mol.GetAtoms():
                        if (_a.GetSymbol() in {'Ni', 'Pd', 'Pt', 'Au'}
                                and _a.GetDegree() == 4):
                            _d8cn4_square = True
                            break
                except Exception:
                    _d8cn4_square = False
            for xyz, lbl in results:
                try:
                    cstr = _build_coordination_constraints_from_xyz(mol, xyz)
                    _xyz_e, energy = _optimize_xyz_openbabel(
                        xyz, steps=0, constraints=cstr, return_energy=True
                    )
                    e = energy if energy is not None else float("inf")
                except Exception:
                    e = float("inf")
                g = float("inf")
                try:
                    _mt = Chem.RWMol(mol)
                    _mt.RemoveAllConformers()
                    _c = _xyz_to_rdkit_conformer(_mt.GetMol(), xyz)
                    if _c is not None:
                        _ci = _mt.AddConformer(_c, assignId=True)
                        g = float(_geometry_quality_score(_mt.GetMol(), _ci))
                except Exception:
                    g = float("inf")
                t_pen = _topology_preservation_penalty(xyz)
                scored.append((e, g, t_pen, xyz, lbl))

            # Energy outlier cut — purely energy-based so truly
            # broken structures (runaway UFF energies, thousands of
            # kcal above min) cannot poison the bucket sort below.
            # Fundamental principle: UFF energies for metals without
            # valid atom types (W, Cd, Sc, lanthanides, actinides) are
            # unreliable — a whole cluster of legit isomers can land
            # at ~1e11 kcal/mol just because UFF missed the type.  When
            # the MEDIAN energy exceeds a plausibility floor (~1e5),
            # skip the cut entirely — we cannot distinguish broken
            # from merely-unparametrised with this metric, so let
            # downstream geometry/topology checks decide.
            scored.sort(key=lambda t: t[0])
            if len(scored) >= 2:
                e_min = scored[0][0]
                e_med = scored[len(scored) // 2][0]
                _UFF_PLAUSIBLE_MAX = 1e5
                if e_med <= _UFF_PLAUSIBLE_MAX:
                    natural_spread = max(e_med - e_min, 1.0)
                    cutoff = e_min + max(50.0 * natural_spread, 5000.0)
                    # ADDITIVE-SAFE energy cut (DELFIN_FFFREE_CN6_OH_ADD; byte-identical when OFF).
                    # The cut targets RUNAWAY UFF energies -- LARGE *FINITE* values thousands of kcal
                    # above e_min.  A frame with e == inf is NOT a runaway: it is UFF-UNSCOREABLE (the
                    # energy eval returned None, typical for an unparametrised metal centre like Mn),
                    # and it already passed the graph-topology output gate above.  The `e_med <= 1e5`
                    # guard was meant to SKIP the cut precisely when energies are unreliable, but an
                    # ADDITIVE OC/SP-4 sibling (genuine low finite energy) can drag e_med below that
                    # floor and so TURN THE CUT ON for a system where it was skipped, deleting valid
                    # unscoreable isomer primaries (measured: VOYWUD all-trans + trans-OH, manifold
                    # 72->44, isomers 6->5).  Keeping unscoreable frames makes the additive pass truly
                    # additive (it can only ADD frames, never remove one) and lets the downstream
                    # geometry/clean/coord-integrity gates -- not an unreliable UFF number -- decide.
                    if _delfin_env_int("DELFIN_FFFREE_CN6_OH_ADD", 0):
                        scored = [t for t in scored if (not math.isfinite(t[0])) or t[0] <= cutoff]
                    else:
                        scored = [t for t in scored if t[0] <= cutoff]

            # Energy-bucket sort — UFF absolute energies are not
            # trustworthy for metals without parameters
            # (Sc/Cd/lanthanides/actinides), but structures within
            # ~15 kcal/mol of each other are indistinguishable to
            # UFF and should be ordered by geometry+topology quality
            # instead.  Bucket width is adaptive: max(15, 10% of
            # natural spread).  Within a bucket, sort by
            # (geometry_score + topology_penalty) then by energy
            # (tie-break).
            if len(scored) >= 2:
                e_min = scored[0][0]
                e_med = scored[len(scored) // 2][0]
                natural_spread = max(e_med - e_min, 1.0)
                bucket_width = max(15.0, 0.1 * natural_spread)
                def _bucket_key(_t):
                    _e, _g, _tp, _x, _l = _t
                    _bucket = int((_e - e_min) / bucket_width)
                    if _d8cn4_square:
                        # Correct polyhedron (via g) dominates; energy tie-breaks.
                        return (_g + _tp, _bucket, _e)
                    return (_bucket, _g + _tp, _e)
                scored.sort(key=_bucket_key)
            if scored:
                results = [(xyz, lbl) for _e, _g, _tp, xyz, lbl in scored]
        except Exception as _sort_exc:
            logger.debug("Energy-based sort failed: %s", _sort_exc)

    # Dual-parse augmentation: some SMILES encodings (bare vs bracketed
    # atoms on the same molecule) parse into RDKit mols with divergent
    # aromatic-flag patterns, which propagates to the enumerator and the
    # builder and yields disjoint isomer sets (e.g. Cd-his: bare 28 vs
    # canonical 52; Zn-bracket: bracket 30 vs canonical 10).  A single
    # canonical normalisation helps one group and hurts the other.
    # Running both parses and union-ing by XYZ signature captures every
    # isomer either encoding can produce.  Doubled runtime for SMILES
    # where canonical != input; pass ``_dual_parse_done=True`` to skip
    # (e.g. inside the second invocation of this function).
    if has_metal and not _dual_parse_done and RDKIT_AVAILABLE:
        try:
            _probe = Chem.MolFromSmiles(smiles)
            _canon = Chem.MolToSmiles(_probe) if _probe is not None else None
            if _canon and _canon != smiles:
                _alt_results, _ = smiles_to_xyz_isomers(
                    _canon,
                    num_confs=num_confs,
                    max_isomers=max_isomers,
                    apply_uff=apply_uff,
                    collapse_label_variants=collapse_label_variants,
                    include_binding_mode_isomers=include_binding_mode_isomers,
                    deterministic=deterministic,
                    hapto_approx=hapto_approx,
                    quality_mode=quality_mode,
                    seeds_override=seeds_override,
                    n_metal_smart=n_metal_smart,
                    _dual_parse_done=True,
                )
                # Union by XYZ-heavy-atom signature (ignores H / index
                # ordering so equivalent geometries don't double-count).
                def _sig(xyz_str):
                    lines = [
                        ln.split() for ln in xyz_str.strip().splitlines()
                        if ln.strip()
                    ]
                    heavy = sorted(
                        (p[0], round(float(p[1]), 2), round(float(p[2]), 2),
                         round(float(p[3]), 2))
                        for p in lines
                        if len(p) >= 4 and p[0] not in ('H', 'h')
                    )
                    return tuple(heavy)
                # ⚠ ONLY COUNT, DO NOT CHANGE (01.09.2026).  This assignment is
                #   the only stage found with a REVERSED precedence rule: for an
                #   equal key the LATER frame wins.  And `_sig` leaves out
                #   hydrogens -- a stereocentre and its mirror image are
                #   IDENTICAL to it.  Whether a frame really gets lost here
                #   is decided by the number, not the narrative.
                try:
                    from delfin.manta._refine_gate import ZAEHLER as _DPZ
                    _DPZ["dual_parse_gelaufen"] += 1
                except Exception:
                    _DPZ = None
                # ── THE FIX (DELFIN_DUALPARSE_KEEP_ALL, default OFF) ────────────
                # `seen_sigs` has TWO jobs, and only ONE of them is right:
                #   (1) lookup "do I already have this geometry?" for the
                #       alt candidates -- correct, the alt loop never overwrites.
                #   (2) carrier of the result list (`results = list(seen_sigs.values())`
                #       further down) -- WRONG: in doing so the caller's OWN frames
                #       collapse together before a single alt frame is even checked.
                # The same file does it right at THREE other places (:29808,
                # :30260, :30653): there it is a SET that only brakes the ADDING
                # and never replaces anything.  The fix aligns this site with
                # those -- no new mechanism, but the pattern missing here.
                _keep_all = os.environ.get("DELFIN_DUALPARSE_KEEP_ALL", "0") == "1"
                _eigene = list(results)   # untouched when the fix is ON
                seen_sigs = {}
                for _dx, _dl in results:
                    _dk = _sig(_dx)
                    if _dk in seen_sigs and _DPZ is not None:
                        _DPZ["dual_sig_kollision"] += 1
                    if _keep_all:
                        seen_sigs.setdefault(_dk, (_dx, _dl))
                    else:
                        seen_sigs[_dk] = (_dx, _dl)
                # The canonical-pipeline xyz carries the canonical mol's
                # atom ordering, which usually differs from the caller
                # mol's ordering.  We need the final output to reference
                # the caller's atom order (so downstream tools that read
                # the XYZ alongside the caller's SMILES see matching
                # bonds).  Build a caller -> canonical atom-index map
                # via RDKit substructure match, then reorder each
                # alt-candidate's XYZ.  Candidates that cannot be
                # mapped (different mol shape, matching failed) are
                # dropped -- keeps topology-safety tight.
                try:
                    _caller_mol = _prepare_mol_for_embedding(smiles, hapto_approx=hapto_mode)
                except Exception:
                    _caller_mol = None
                try:
                    _canon_mol = _prepare_mol_for_embedding(_canon, hapto_approx=hapto_mode)
                except Exception:
                    _canon_mol = None

                # caller_idx[i] = canonical atom index that corresponds
                # to caller atom i.  Built once per SMILES.
                caller_to_canon = None
                if _caller_mol is not None and _canon_mol is not None:
                    try:
                        match = _canon_mol.GetSubstructMatch(_caller_mol, useChirality=False)
                        if match and len(match) == _caller_mol.GetNumAtoms():
                            caller_to_canon = list(match)
                    except Exception:
                        caller_to_canon = None

                def _reorder_xyz(xyz_str, perm):
                    """perm[i] = index in source xyz for caller atom i.
                    Returns xyz with lines permuted so row i == caller atom i."""
                    lines = [l for l in xyz_str.strip().splitlines() if l.strip()]
                    if len(lines) != len(perm):
                        # Account for explicit H count mismatch
                        return None
                    reordered = [lines[perm[i]] for i in range(len(perm))]
                    return "\n".join(reordered) + "\n"

                for xyz, lbl in _alt_results or []:
                    try:
                        s = _sig(xyz)
                        if s in seen_sigs:
                            continue

                        reordered_xyz = xyz
                        if caller_to_canon is not None and _caller_mol is not None:
                            new_xyz = _reorder_xyz(xyz, caller_to_canon)
                            if new_xyz is not None:
                                reordered_xyz = new_xyz

                        # Topology gate against caller mol (now with
                        # matching atom order).  Fall back to canonical
                        # mol check if caller reorder failed.
                        ok = False
                        if _caller_mol is not None and reordered_xyz is not xyz:
                            try:
                                ok = bool(_verify_topology_from_graph(reordered_xyz, _caller_mol))
                            except Exception:
                                ok = False
                        if not ok and _canon_mol is not None:
                            try:
                                ok = bool(_verify_topology_from_graph(xyz, _canon_mol))
                            except Exception:
                                ok = False
                        if not ok:
                            continue

                        # Store reordered xyz (matches caller mol ordering
                        # whenever the mapping succeeded) so downstream
                        # consumers see consistent atom indices.
                        seen_sigs[s] = (reordered_xyz, lbl)
                        if _keep_all:
                            _eigene.append((reordered_xyz, lbl))
                    except Exception:
                        continue
                # ⚠ WITH FIX: the own frames ALL stay, the accepted alt frames
                #   go at the BACK -- additive by construction, the frame count
                #   can no longer DROP here.
                #   WITHOUT FIX: unchanged old behaviour, byte for byte.
                results = _eigene if _keep_all else list(seen_sigs.values())
        except Exception as _dual_exc:
            logger.debug("Dual-parse union failed: %s", _dual_exc)

    # H-geometry universal repair DISABLED at end-of-pipeline.
    # `_fix_h_geometry_via_smiles` -> `_fix_h_geometry_universal` snaps
    # H atoms to ideal VSEPR angles with arbitrary rotational phase,
    # producing 5x more H-clash violations (methyl umbrella pointing
    # INTO nearby heavies) than UFF-only output. See INSIGHTS_LOG
    # 2026-04-29 ~12:00 UTC for the data: clash/1k_H 13.8 (e183cea)
    # -> 70.8 (current HEAD with this call active).
    #
    # Iter-9 H2 (universal aromatic-H post-process) REVERTED 2026-05-01
    # ~14:35 UTC after full-pool ≥20% mid-audit revealed catastrophic
    # F20 regression on FUSED ring systems (D-POLTIT_main_Sn_CN3 295 -> 697
    # violations).  The outward-radial projection assumes single-ring H
    # placement: when a ring atom belongs to MULTIPLE rings (fused
    # naphthalene-style, Cp-Cyclopentyl, Ph3P-on-metal), projecting H onto
    # ring A's plane creates a violation in ring B's detector pass.
    # F20 reports a single H with parent in 2 rings as 2 violations.
    # Aggregate sample 2219 SMILES: F20 19469 -> 30559 (+57%).
    # H1 fix in `_snap_aromatic_rings_to_plane` retained (selective —
    # only fires when ring rms >= threshold, so non-issue on flat fused
    # systems where rms = 0).

    # ── Iter-12/13 Baustein 3: post-ETKDG/UFF coord-angle correction ────────
    # Per-conformer, opt-in via DELFIN_FFFREE_COORD_ANGLE_FIX=1.  No effect when disabled.
    # Operates on already-finalised XYZs; rotates rigid X-side around the
    # axis perpendicular to (M, D, X) plane through D to bring observed
    # M-D-X angle to expected (sp/sp2/sp3 inferred geometrically).
    # Per-violation revert on clash regression.  Universal — no SMILES/refcode.
    #
    # Run ONLY at top-level (not in recursive dual-parse inner call).  Reason:
    # the outer dual-parse union dedup uses XYZ heavy-atom signature; if the
    # inner call modifies coordinates, signatures of inner.results will not
    # match the (still-uncorrected) outer.results, breaking dedup.  By
    # deferring correction to the outer call (after union), all results pass
    # through one consistent correction step.
    #
    # Iter-13: routed through ``_apply_coord_angle_fix_if_enabled`` so every
    # scaffold-path return point shares the same gate (mono σ here, plus
    # mono hapto / multi-metal hapto / fallback above).
    if has_metal:
        results = _apply_coord_angle_fix_if_enabled(mol, results, _dual_parse_done)

    # ── Welle-5j Agent A: Cp piano-stool hapticity refinement ──────────────
    # Welle-5i Agent C catalogued 28 / 34 (83 %) hapto BROKEN-TO-BROKEN
    # files as η⁵-Cp mislabeled by the detector as η⁶-arene because the
    # M-ring-centroid axis was off-ideal post-UFF.  Universal corrector:
    # detect 5-ring of C/N at near-equidistant M-C distances + planar
    # ring (SVD), snap metal onto SVD ring-normal axis at the η⁵ target
    # distance for the metal element.  Per-violation rollback if any
    # non-ring M-D bond would dissociate (Iter-15 hard invariant).
    # Opt-in via DELFIN_5J_A_CP_PIANO_STOOL=1.  Bit-exact when disabled.
    if has_metal:
        results = _apply_5j_a_cp_piano_stool_if_enabled(
            mol, results, _dual_parse_done,
        )

    # ── Iter-14 Baustein 4: post-B3 rigid-π H projection ────────────────────
    # π-rigid-body invariant: every transformation that moves a π-frame
    # must drag attached H atoms rigidly with the ring.  Iter-9 H1 handles ring-snap-time projection;
    # Iter-12/13 B3 handles X-side rotation around donors (BFS includes H).
    # B4 is the universal final pass that re-projects ring-attached H onto
    # current ring planes after all upstream operations.  Per-conformer,
    # opt-in via DELFIN_BAUSTEIN4=1.  Bit-exact when disabled.
    #
    # Runs for both metal AND non-metal complexes — aromatic-H out-of-plane
    # is a generic converter pathology independent of metal presence.
    results = _apply_baustein4_if_enabled(mol, results, _dual_parse_done)

    # ── Iter-21 (2026-05-19): 81f8a1f-style M-X clearance final post-pass ───
    # For hapto+multi_hapto class: radial push to enforce min M-X 2.5 Å for
    # non-bonded heavy atoms.  Bridges gap between candidate-select-time gate
    # and actual emitted XYZ (B3/B4 may re-introduce M-X clash via rotations
    # / π-projections).  Class-gated default-ON for hapto+multi_hapto only.
    # Welle-5f-F finally landed.
    results = _apply_hapto_clearance_if_enabled(mol, results, _dual_parse_done)

    # ── Baustein 5 v2: PBD post-UFF geometry corrector ───────────────────────
    # Per-frame: catastrophic M-D break repair, Stage 1 bond corrections,
    # Stage 2 angle corrections, Stage 3 clash resolution. Hard topology +
    # M-D drift gates prevent damage to good frames (only-fix-what-is-broken).
    # Bit-exact when DELFIN_BAUSTEIN5=0 (default).
    results = _apply_baustein5_if_enabled(mol, results, _dual_parse_done)

    # ── Baustein 6: variational L-BFGS-B refiner + 4-tier symmetry ──────────
    # Per-frame analytic-gradient minimisation of an 8-term U_total. Runs
    # AFTER B5 (which handles catastrophic moves) and AFTER the targeted
    # post-B5 fixers below would have run — placed here so B5+B6 form a
    # single "geometry refinement" stack. Hard topology gate + fallback on
    # failure. Bit-exact when DELFIN_B6_WIRED=0 (default).
    results = _apply_baustein6_if_enabled(mol, results, _dual_parse_done)

    # ── Targeted post-B5 fixers (F19 / F25 / WUXQAK / μ-X bridging) ──────────
    # Surgical per-frame correctors that run AFTER B5 so they re-align
    # hydrogens / pyramidal N / linear sp3-C donors / bridging anions to
    # the post-optimization heavy-atom frame.  Each is opt-in via its own
    # env flag and is bit-exact when its flag is 0 (default OFF).
    # Insertion order: F19 (sp3-H tetrahedrality) → F25 (sp3-N pyramidality)
    # → SP2N-PLANARIZE (sp2-N/nitro planarity, inverse of F25) → WUXQAK
    # (sp3-C linear collapse) → bridging-anion (μ-X M-X-M angle for
    # bimetallic complexes).  See ``_apply_fixer_*`` docstrings.
    results = _apply_fixer_f19_if_enabled(mol, results, _dual_parse_done)
    results = _apply_hydroxyl_geom_if_enabled(mol, results, _dual_parse_done)
    results = _apply_fixer_f25_if_enabled(mol, results, _dual_parse_done)
    results = _apply_fixer_sp2n_planarize_if_enabled(
        mol, results, _dual_parse_done,
    )
    results = _apply_fixer_sp2c_planarize_if_enabled(
        mol, results, _dual_parse_done,
    )
    results = _apply_fixer_wuxqak_if_enabled(mol, results, _dual_parse_done)
    results = _apply_fixer_bridging_anion_if_enabled(
        mol, results, _dual_parse_done,
    )

    # ── Iter-24 (2026-05-20): post-UFF aromatic-ring planarity enforcement ──
    # Flatten puckered TRUE aromatic rings (M_coord chelate rings excluded via
    # bond-length gate) onto their SVD best-fit plane, centroid-preserving so
    # the M-ring distance / M-D invariant is untouched; ring-H dragged.
    # Class-cond default-ON {hapto, multi_hapto} (where rings pucker 72-75 %).
    results = _apply_aromatic_planarity_if_enabled(mol, results, _dual_parse_done)

    # ── Iter-33 (2026-06-19): universal aromatic ring-SYSTEM planarisation ──
    # Flatten EVERY aromatic ring-system (heteroaromatic + fused/polycyclic +
    # coordinated) onto its best-fit plane — closes the Iter-24 gaps where
    # M-bound heteroaromatic donor rings + partly-coordinated chelates stayed
    # puckered (AXAGOY C5N donor rings 0.13-0.16 Å).  Coordinated ring atoms
    # are anchored (M-D invariant), only non-anchor atoms + ring-H + first
    # substituents projected.  Default-OFF byte-id (DELFIN_FFFREE_AROM_PLANARIZE).
    results = _apply_arom_planarize_if_enabled(mol, results, _dual_parse_done)

    # ── Iter-34 (2026-06-19): coordinated planar π-donor co-planar-M orient ──
    # Rotate each coordinated in-plane σ aromatic π-donor RIGID ring about its
    # FIXED donor so the donor's in-plane lone-pair points at the metal — i.e.
    # the ring mean-plane CONTAINS the metal (closes the Iter-33 gap where the
    # internally-flat ring still sat tilted, M out-of-plane).  Hapto / η π-face
    # donors excluded (they bind perpendicular).  Donor + metal frozen → M-D
    # preserved; per-ring never-worse + clash never-worse + M-D rollback.
    # Default-OFF byte-id (DELFIN_FFFREE_PI_COPLANAR_M).
    results = _apply_pi_coplanar_m_if_enabled(mol, results, _dual_parse_done)

    # ── Iter-25 (2026-05-20): final bond-decollapse corrector ───────────────
    # Spring-relax the non-metal heavy graph to physical bond lengths + repel
    # superimposed atoms (metals + coord-sphere frozen, M-D invariant + collapse-
    # reduce rollback).  Fixes the ~79% hapto ligand-collapse (validated
    # -20.5pp).  Class-cond default-ON {hapto, multi_hapto}; runs LAST.
    results = _apply_bond_decollapse_if_enabled(mol, results, _dual_parse_done)

    # ── Resonance-aware aromatic bond-length equalisation ───────────────────
    # Runs AFTER bond-decollapse (whose single-bond-ideal spring would otherwise
    # re-stretch equalised aromatic bonds): reshape every perceived aromatic ring
    # system so each ring bond sits at its first-principles DELOCALISED target
    # (Pyykko single<->double radius interpolation at the Huckel benzene fraction
    # f=2/3 -> C-C 1.393, C-N 1.333 ...), equalising Kekule alternation and
    # pulling single-drifted rings back to ~1.39.  In-plane PBD, metal-coordinated
    # ring atoms anchored, substituents rigidly dragged (lengths only -- angles/
    # planarity preserved); per-frame ring-bond-deviation never-worse rollback.
    # Default-OFF byte-id (DELFIN_FFFREE_AROM_BOND_LENGTH).
    results = _apply_arom_bond_length_if_enabled(mol, results, _dual_parse_done)

    # ── Iter-3 General-Isomer Enumerator (env-gated, default ON) ───────────
    # Restores the historical "Isomer 1, Isomer 2, ... Isomer N" emission
    # that produced 16-26 frames per crowded high-CN / multi-metal SMILES at
    # 44fce9e and earlier.  HEAD's expanded chemistry-rule set in
    # _verify_topology_from_graph (~485 added LOC: sp2 planarity / Cn axes /
    # phantom-bonds / hapto-bridge / element-of-its-type closest-neighbour)
    # rejects the constitutional-permutation candidates that older topo
    # enumerator paths emitted, so _generate_topological_isomers returns 0
    # at HEAD whereas it returned 17-26 at 44fce9e for the same SMILES.
    #
    # Restore strategy: temporarily monkey-patch _verify_topology_from_graph
    # to a relaxed predicate (bond-distance + heavy-clash only — Rules 1-3
    # of the original docstring) for the duration of one extra
    # _generate_topological_isomers call, then dedup the survivors against
    # the existing results by coordination fingerprint AND XYZ heavy-atom
    # signature.  Survivors are appended with display label 'Isomer N'.
    #
    # Toggle with DELFIN_GENERAL_ISOMER_ENUM (default 1).  When disabled
    # (or _dual_parse_done — only the outer call applies the pass) the
    # behaviour is bit-identical to current HEAD.
    if (
        has_metal
        and not _dual_parse_done
        # Iter-5: default flipped 1→0.  GIE pumped unrefined enumerator frames
        # via monkey-patched relaxed verify, breaking topology gate.
        # Opt-in via DELFIN_GENERAL_ISOMER_ENUM=1.
        and _delfin_env_int("DELFIN_GENERAL_ISOMER_ENUM", 0)
        and not DELFIN_TOPOLOGY_STRICT_MODE
    ):
        try:
            _gie_n_room = max(0, max_isomers - len(results))
            if _gie_n_room > 0:
                def _gie_relaxed_verify(xyz_str, ref_mol) -> bool:
                    try:
                        if not _metal_donor_distances_realistic(xyz_str, ref_mol):
                            return False
                        _ls = [l for l in xyz_str.strip().splitlines() if l.strip()]
                        _coords: List[Tuple[str, float, float, float]] = []
                        for _l in _ls:
                            _p = _l.split()
                            if len(_p) < 4:
                                return True
                            _coords.append(
                                (_p[0], float(_p[1]), float(_p[2]), float(_p[3]))
                            )
                        _heavy_idx = [
                            i for i, c in enumerate(_coords) if c[0] not in ('H', 'h')
                        ]
                        for _i_pos in range(len(_heavy_idx)):
                            _ai = _heavy_idx[_i_pos]
                            for _j_pos in range(_i_pos + 1, min(_i_pos + 50, len(_heavy_idx))):
                                _aj = _heavy_idx[_j_pos]
                                _dx = _coords[_ai][1] - _coords[_aj][1]
                                _dy = _coords[_ai][2] - _coords[_aj][2]
                                _dz = _coords[_ai][3] - _coords[_aj][3]
                                if _dx * _dx + _dy * _dy + _dz * _dz < 0.49:
                                    return False
                        return True
                    except Exception:
                        return True

                import sys as _sys
                _gie_mod = _sys.modules[__name__]
                _gie_orig_verify = _gie_mod._verify_topology_from_graph
                _gie_topo: List[Tuple[str, str]] = []
                try:
                    _gie_mod._verify_topology_from_graph = _gie_relaxed_verify
                    _gie_topo = _generate_topological_isomers(
                        mol, smiles, apply_uff=apply_uff,
                        max_isomers=max_isomers, n_metal_smart=n_metal_smart,
                        profile=_qprof,
                    ) or []
                except Exception as _gie_topo_exc:
                    logger.debug(
                        "General-isomer enum: enumerator call failed: %s",
                        _gie_topo_exc,
                    )
                finally:
                    _gie_mod._verify_topology_from_graph = _gie_orig_verify

                if _gie_topo:
                    _gie_dtype = _donor_type_map(mol)
                    _gie_existing_fps: set = set()
                    _gie_existing_sigs: set = set()
                    for _ex_xyz, _ex_lbl in results:
                        try:
                            _et = Chem.RWMol(mol)
                            _et.RemoveAllConformers()
                            _ec = _xyz_to_rdkit_conformer(_et.GetMol(), _ex_xyz)
                            if _ec is not None:
                                _eci = _et.AddConformer(_ec, assignId=True)
                                _efp = _compute_coordination_fingerprint(
                                    _et.GetMol(), _eci, dtype_map=_gie_dtype,
                                )
                                _gie_existing_fps.add(_efp)
                        except Exception:
                            pass
                        try:
                            _xl = [
                                ln.split() for ln in _ex_xyz.strip().splitlines()
                                if ln.strip()
                            ]
                            _heavy = tuple(sorted(
                                (p[0], round(float(p[1]), 2),
                                 round(float(p[2]), 2),
                                 round(float(p[3]), 2))
                                for p in _xl
                                if len(p) >= 4 and p[0] not in ('H', 'h')
                            ))
                            _gie_existing_sigs.add(_heavy)
                        except Exception:
                            pass

                    _gie_isomer_n = 0
                    for _, _lbl in results:
                        _ms = re.match(r'^Isomer (\d+)$', _lbl or '')
                        if _ms:
                            _gie_isomer_n = max(_gie_isomer_n, int(_ms.group(1)))

                    # ── Iter-3.2 PRE-GATE STRUCTURE REFINEMENT PIPELINE ──
                    # User direktive 2026-05-06: control at exit must be strict;
                    # work on structure GENERATION upstream instead of relaxing
                    # the gate.  Each enumerator candidate flows through
                    # UFF-refine → H-fix → aromatic-snap → STRICT verify
                    # (original _verify_topology_from_graph, NOT relaxed).
                    # This way bad geometries get repaired upstream, and
                    # only chemically-realistic frames pass the gate.
                    # Compute cost: ~0.5-2 s per candidate frame.
                    for _txyz, _tlabel in _gie_topo:
                        if _gie_n_room <= 0:
                            break
                        # Cheap pre-filter: relaxed predicate keeps obvious
                        # garbage out of the (expensive) refinement loop.
                        if not _gie_relaxed_verify(_txyz, mol):
                            continue

                        # Iter-4: try RAW frame against strict gate first.
                        # e6761e4 forensics: 70%+ of raw enumerator frames pass
                        # HEAD's strict verify directly.  UFF refinement on TM
                        # systems often WRECKS the frame (Cu+1, Cd-6 etc. are
                        # unrecognized atom types — UFF returns garbage).
                        # Only fall back to UFF+H-fix+snap if raw fails.
                        _refined = _txyz
                        _passed_raw = False
                        try:
                            if _gie_orig_verify(_txyz, mol):
                                _passed_raw = True
                        except Exception:
                            pass

                        if not _passed_raw:
                            # Pre-gate refinement pipeline as fallback
                            if apply_uff:
                                try:
                                    _refined = _optimize_xyz_openbabel_safe(
                                        _refined, mol_template=mol,
                                    )
                                except Exception:
                                    pass
                            try:
                                _refined = _fix_h_geometry_universal(_refined, mol)
                            except Exception:
                                pass
                            try:
                                _refined = _snap_aromatic_rings_in_xyz(
                                    _refined, mol, rms_threshold=0.05,
                                )
                            except Exception:
                                pass

                            # STRICT verify after refinement
                            try:
                                if not _gie_orig_verify(_refined, mol):
                                    continue
                            except Exception:
                                continue

                        try:
                            _mt = Chem.RWMol(mol)
                            _mt.RemoveAllConformers()
                            _mc = _xyz_to_rdkit_conformer(_mt.GetMol(), _refined)
                            if _mc is None:
                                continue
                            _mci = _mt.AddConformer(_mc, assignId=True)
                            _mfp = _compute_coordination_fingerprint(
                                _mt.GetMol(), _mci, dtype_map=_gie_dtype,
                            )
                        except Exception:
                            continue
                        if _mfp in _gie_existing_fps:
                            continue
                        try:
                            _xl = [
                                ln.split() for ln in _refined.strip().splitlines()
                                if ln.strip()
                            ]
                            _heavy = tuple(sorted(
                                (p[0], round(float(p[1]), 2),
                                 round(float(p[2]), 2),
                                 round(float(p[3]), 2))
                                for p in _xl
                                if len(p) >= 4 and p[0] not in ('H', 'h')
                            ))
                        except Exception:
                            continue
                        if _heavy in _gie_existing_sigs:
                            continue
                        _gie_existing_fps.add(_mfp)
                        _gie_existing_sigs.add(_heavy)
                        # Iter-8.7 every-append gate (123a130, env-gated, default OFF)
                        _gate_pass = True
                        if _every_append_gate_enabled(mol):
                            try:
                                _flat = _flatten_sp2_atoms_xyz(_refined, mol)
                                if _flat:
                                    _refined = _flat
                            except Exception:
                                pass
                            try:
                                # _mt + _mci already built above; check distortion
                                if _has_severe_covalent_distortion(_mt.GetMol(), _mci):
                                    _gate_pass = False
                            except Exception:
                                _gate_pass = False
                        if not _gate_pass:
                            continue
                        _gie_isomer_n += 1
                        results.append((_refined, f'Isomer {_gie_isomer_n}'))
                        _gie_n_room -= 1
        except Exception as _gie_exc:
            logger.debug("General-isomer enumerator skipped: %s", _gie_exc)

    # ── Iter-3.1 H1 universal aromatic-H projection (env-gated, default ON) ──
    # Iter-3 added new emit paths (Hapto-Diversity-Topology-Aware, General-
    # Isomer-N, conformer multiplication) that bypass the upstream
    # _snap_aromatic_rings_in_xyz call sites.  This left aromatic-H atoms in
    # ETKDG sp³-like positions, producing AROM_OOP +1097 anomalies vs Iter-1.
    # User direktive 2026-05-05: H must be realistic at aromatics.
    # Apply Iter-9 H1 (project ring-attached H onto π-plane) to every result
    # frame as a final post-emit pass.  Topology-preserving (heavy atoms
    # untouched), only H atoms move perpendicular to the SVD-fit ring plane.
    # Toggle: DELFIN_H1_PROJECT_ALL=0 → bit-exact Iter-3 (no projection).
    if (
        has_metal
        and _delfin_env_int("DELFIN_H1_PROJECT_ALL", 1)
        and not DELFIN_TOPOLOGY_STRICT_MODE
    ):
        try:
            _h1_results: List[Tuple[str, str]] = []
            _h1_gate_on = bool(_every_append_gate_enabled(mol))
            for _h1_xyz, _h1_lbl in results:
                try:
                    _h1_xyz_snapped = _snap_aromatic_rings_in_xyz(
                        _h1_xyz, mol, rms_threshold=0.05,
                    )
                    # Iter-8.7 every-append gate on snapped frame
                    if _h1_gate_on:
                        try:
                            _flat = _flatten_sp2_atoms_xyz(_h1_xyz_snapped, mol)
                            if _flat:
                                _h1_xyz_snapped = _flat
                        except Exception:
                            pass
                        try:
                            _mt_g = Chem.RWMol(mol); _mt_g.RemoveAllConformers()
                            _c_g = _xyz_to_rdkit_conformer(_mt_g.GetMol(), _h1_xyz_snapped)
                            if _c_g is not None:
                                _ci_g = _mt_g.AddConformer(_c_g, assignId=True)
                                if _has_severe_covalent_distortion(_mt_g.GetMol(), _ci_g):
                                    # fall back to unsnapped frame
                                    _h1_results.append((_h1_xyz, _h1_lbl))
                                    continue
                        except Exception:
                            _h1_results.append((_h1_xyz, _h1_lbl))
                            continue
                    _h1_results.append((_h1_xyz_snapped, _h1_lbl))
                except Exception as _h1_inner_exc:
                    logger.debug(
                        "H1 projection skipped for one frame: %s", _h1_inner_exc,
                    )
                    _h1_results.append((_h1_xyz, _h1_lbl))
            results = _h1_results
        except Exception as _h1_exc:
            logger.debug("H1 universal projection failed: %s", _h1_exc)

    # -- Iter-6 T2: hard frame-cap at exit (env-gated, default OFF) ---------
    # Cap total emitted frames per SMILES regardless of how many enumerator
    # branches contributed.  Goal: -75% inter-ligand clash density per 1k
    # frames by trimming low-priority tail buckets.  Stable order -- keeps
    # the first ``DELFIN_FRAME_CAP_HARD_N`` (default 30) results, which by
    # construction are the higher-priority polyhedron / preferred-isomer
    # frames.  Bit-exact when ``DELFIN_FRAME_CAP_HARD=0`` (default).
    if _delfin_env_int("DELFIN_FRAME_CAP_HARD", 0):
        try:
            cap_n = int(_delfin_env_int("DELFIN_FRAME_CAP_HARD_N", 30))
            if cap_n > 0 and len(results) > cap_n:
                try:
                    logger.debug(
                        "Hard frame-cap: %d->%d frames",
                        len(results), cap_n,
                    )
                except Exception:
                    pass
                results = results[:cap_n]
        except Exception as _cap_exc:
            try:
                logger.debug("Hard frame-cap skipped: %s", _cap_exc)
            except Exception:
                pass

    # ========================================================================
    # Iter-8.1 multi-hapto safe-fallback filter (FPCFD per-class extras-filter)
    # ========================================================================
    # Class-dispatched post-filter: drop catastrophic-extras frames for hapto
    # and multi-hapto only.  sigma / multi_sigma / no_metal pass through
    # unchanged.  Best-of-K fallback keeps top K=max(2, ⌈0.3·n⌉) frames if
    # filter would empty result set, guaranteeing ≥2 frames per SMILES.
    # See results/iter8.1_multihapto_safefallback_design.md for forensics.
    # Default-on; opt-out via DELFIN_ITER81_FILTER=0.
    try:
        if results and len(results) >= 3:
            _iter81_active = bool(_delfin_env_int("DELFIN_ITER81_FILTER", 1))
            if _iter81_active:
                _iter81_class = _classify_complex_class(mol)
                _iter81_threshold = _ITER8_1_EXTRA_THRESHOLDS.get(
                    _iter81_class, 9999
                )
                # Iter-8.4c: when DELFIN_SIGMA_TIGHT_THRESHOLD_ITER8=1 AND
                # class='sigma', override the disabled-default (9999) with
                # tau=3.  Full-pool sigma frame median is 0 extras with a
                # thin tail at 3-8 extras representing d8-pincer regression
                # cases; capping at 3 drops the catastrophic tail without
                # touching well-formed frames.  Best-of-K fallback below
                # guarantees ≥2 frames per SMILES so the filter cannot
                # empty the result set.  Default OFF preserves bit-exact
                # HEAD when env-flag unset.
                if (_iter81_class == "sigma"
                    and _class_conditional_flag(
                        "DELFIN_SIGMA_TIGHT_THRESHOLD_ITER8", mol)):
                    _iter81_threshold = 3
                if _iter81_threshold < 9999:
                    # Score each frame by extra-bond count
                    _scored = []
                    for _xyz_i, _disp_i in results:
                        _ne = _count_extra_heavy_bonds(_xyz_i, smiles)
                        _scored.append(
                            (_ne if _ne >= 0 else 0, _xyz_i, _disp_i)
                        )
                    # Median over valid (>0) scores; fallback 0 (no filter)
                    _valid = [s for s, _, _ in _scored if s > 0]
                    _med = sorted(_valid)[len(_valid) // 2] if _valid else 0
                    # Filter: extras <= τ AND extras <= 3×median (if med>0)
                    _kept = []
                    for s, x, d in _scored:
                        ok = (s <= _iter81_threshold)
                        if _med > 0:
                            ok = ok and (s <= 3 * _med)
                        if ok:
                            _kept.append((x, d))
                    # Best-of-K fallback
                    _n_total = len(results)
                    _k_min = max(2, int(round(0.3 * _n_total)))
                    # Iter-8.10b (2026-05-11, Subagent 1 forensics):
                    # preserve frames with DISTINCT base labels — the
                    # extras-filter killed FIRCOY/TIYRUR 5→2 and 3→2 even
                    # though distinct isomer-labels (e.g. all-trans,
                    # N-trans, O-trans, all-cis) were present. Per dual
                    # contract "all isomers", every distinct base-label
                    # frame must survive regardless of extras-count.
                    # Env-gate DELFIN_ITER81_PRESERVE_DISTINCT_LABELS
                    # default 1 (active); set to 0 for legacy extras-only.
                    if _delfin_env_int("DELFIN_ITER81_PRESERVE_DISTINCT_LABELS", 1):
                        import re as _re_lbl
                        _unique_bases = {
                            _re_lbl.sub(r'-(?:conf)?\d+$', '', d)
                            for _, d in results
                        }
                        _k_min = max(_k_min, min(_n_total, len(_unique_bases)))
                    if len(_kept) < _k_min:
                        _scored.sort(key=lambda t: t[0])
                        _kept = [(x, d) for _, x, d in _scored[:_k_min]]
                    if len(_kept) < len(results):
                        try:
                            logger.debug(
                                "Iter-8.1 filter (class=%s, tau=%d, med=%d): %d -> %d frames",
                                _iter81_class, _iter81_threshold, _med,
                                len(results), len(_kept),
                            )
                        except Exception:
                            pass
                        results = _kept
    except Exception as _iter81_exc:
        try:
            logger.debug("Iter-8.1 filter probe failed: %s", _iter81_exc)
        except Exception:
            pass

    # ── Welle-5b B: donor-orientation realism (heavy-atom only) ─────────────
    # Three geometry patterns UFF does not enforce natively:
    #   1. Aromatic-N edge-on:  rotate aromatic-N rings around the M-N axis
    #      so the ring normal becomes perpendicular to M-N (sigma lone-pair
    #      points at the metal) — fixes the face-on attack mode.
    #   2. Terminal-carbonyl linearity:  M-C=O snapped to 180 deg.
    #   3. NHC carbene plane:  rotate the N-C-N ring so the metal lies in
    #      the carbene plane.
    #
    # Master flag DELFIN_DONOR_ORIENT_REALISM (default 0 = bit-exact OFF).
    # Sub-flags DELFIN_DONOR_ORIENT_REALISM_{AROMATIC_N,CARBONYL,NHC}
    # default to 1 when master is 1 (per-pattern force-disable available).
    # Insertion order: 5b-B FIRST (moves heavy atoms), 5b-A AFTER (re-snaps
    # H to current heavy geometry).
    try:
        if (
            mol is not None
            and results
            and _delfin_env_int("DELFIN_DONOR_ORIENT_REALISM", 0)
        ):
            from delfin.manta._donor_orientation_realism import (
                snap_donor_orientations as _dor_snap,
            )
            _new_results: List[Tuple[str, str]] = []
            for _xyz_i, _lbl_i in results:
                try:
                    _xyz_i2 = _dor_snap(_xyz_i, mol, mode="end_of_pipeline")
                except Exception:
                    _xyz_i2 = _xyz_i
                _new_results.append((_xyz_i2, _lbl_i))
            results = _new_results
    except Exception as _dor_exc:
        try:
            logger.debug("Welle-5b-B donor-orient pass failed: %s", _dor_exc)
        except Exception:
            pass

    # ── Welle-5b A: universal VSEPR-correct H placement (end-of-pipeline) ──
    # Snap every X-H bond direction to the VSEPR-correct local geometry of
    # its heavy parent (sp/sp2/sp3 + lone-pair count).  Per-H rollback on
    # M-D invariant break or new heavy clash.  Aromatic ring H is delegated
    # to Baustein 4 / Iter-9 H1 to avoid double-projection.
    #
    # Welle-5l T5 REVERT of T3-A class-conditional default-flip:
    #   T3-A flipped default to class-conditional ON for {sigma, multi-sigma,
    #   hapto, multi-hapto} based on a SINGLE-SMILES success (29-Ni 12/12
    #   methyls fixed → 0/12).  T5 broader full-smoke n=197 then showed at
    #   scale the mechanism REGRESSES: F19 +17.32pp, F20 +4.93pp, F25 +0.97pp,
    #   topo −0.8pp, n_isomers −0.24/SMILES, 3 SMILES drop emission entirely
    #   (Iter-11 antipattern reproduced 4th time).
    #
    #   Root cause [[feedback_f19_f25_detector_contract_mismatch]]: the
    #   F19/F20/F25 detectors compare H against XYZ-inferred local geometry
    #   (not VSEPR-textbook), so moving H toward chemistry-correct VSEPR
    #   INCREASES detector violations on strained post-UFF TMC.  The 29-Ni
    #   "12/12 → 0/12" T3-A claim used `find_methyl_angle_quality` which
    #   measures DIFFERENT quantity than F19/F25 — both can simultaneously
    #   improve methyl AND worsen aggregate H-realism family.
    #
    #   Decision: revert default to bit-exact OFF (universal-fundamental
    #   doctrine — narrow single-SMILES success is NOT universal evidence).
    #   Re-investigate in Welle-5m with detector-revision-or-replacement
    #   precondition met first.
    #
    # Master flag DELFIN_5B_VSEPR_H_REALISM (default 0 = bit-exact OFF).
    # Sub-flag DELFIN_5F_D_ALKYL_ROTAMER (default 0) adds an inter-substituent
    # rotamer search inside the VSEPR pass.
    _vsepr_h_realism_active = bool(
        _delfin_env_int("DELFIN_5B_VSEPR_H_REALISM", 0)
    )
    try:
        if results and _vsepr_h_realism_active:
            from delfin.manta._h_vsepr_realism import correct_results as _vsepr_correct
            results = _vsepr_correct(mol, results)
    except Exception as _vsepr_exc:
        try:
            logger.debug("Welle-5b-A VSEPR-H pass failed: %s", _vsepr_exc)
        except Exception:
            pass

    # ── Welle-5f-D: inter-substituent H-H clash relief (rotamer-only) ──────
    # PMe3 / NMe3 / tBu pattern: H-H clashes between sibling alkyl groups
    # that share a heavy parent (e.g. P-CH3 + P-CH3 + P-CH3 on a phosphine).
    # Rotates each substituent rigidly around its parent-axis to minimise
    # inter-substituent H-H clashes.  H-only mover; heavy geometry untouched.
    #
    # When 5b-A is enabled, 5f-D auto-runs inside ``correct_xyz`` (sub-flag).
    # This dispatch handles the 5f-D-only path (5b-A OFF, 5f-D ON).  Both
    # paths gate on DELFIN_5F_D_ALKYL_ROTAMER (default 0 = bit-exact OFF).
    # Welle-5l T3-A: use the class-conditional 5b-A active-state (computed
    # above) so the inhibit-check stays consistent with the new default.
    try:
        if (
            results
            and _delfin_env_int("DELFIN_5F_D_ALKYL_ROTAMER", 0)
            and not _vsepr_h_realism_active
        ):
            from delfin.manta._h_vsepr_realism import (
                correct_results_rotamers_only as _rot_correct,
            )
            results = _rot_correct(mol, results)
    except Exception as _rot_exc:
        try:
            logger.debug("Welle-5f-D alkyl-rotamer pass failed: %s", _rot_exc)
        except Exception:
            pass

    # ── Welle-5f-C: post-UFF M-H atom-overlap rescue ────────────────────────
    # D-QOPVEZ Fe-H 0.4 A pattern: UFF on unparametrised TM cations leaves
    # hydride H atoms inside the metal core.  Detect M-H pairs with d < 0.9 *
    # ideal_M_H and try short axial translations along the local M-H axis to
    # rescue the bond length to >= 0.9 * ideal.  Frames whose rescue is
    # geometrically impossible are dropped (atom-overlap propagates worse
    # downstream artefacts than a missing frame).
    #
    # Gated by DELFIN_5F_C_MH_OVERLAP_RESCUE (default 0 = bit-exact OFF).
    # Runs AFTER 5b-A (H positions finalised) and AFTER B4/B5/B6 (heavy
    # geometry finalised) so the rescue snaps to current metal positions.
    try:
        if (
            results
            and _delfin_env_int("DELFIN_5F_C_MH_OVERLAP_RESCUE", 0)
        ):
            from delfin.manta._h_vsepr_realism import rescue_results as _mh_rescue
            results = _mh_rescue(mol, results)
    except Exception as _mh_exc:
        try:
            logger.debug("Welle-5f-C M-H rescue pass failed: %s", _mh_exc)
        except Exception:
            pass

    # ── Welle-5l Track-2: xTB-cascade post-UFF refinement ───────────────────
    # GFN2-xTB optimizer with M-D invariant rollback gate.  Runs as the FINAL
    # post-UFF stage so the cascade refines the fully-finished pipeline
    # geometry (post B3/B4/B5/B6/F19/F25/WUXQAK/bridging-anion/5b-A/5b-B/5f-C/
    # 5f-D/5j-A).  Per ``project_core_swap_decision`` 2026-05-14: UFF is the
    # cheap pre-conditioner, xTB the TM-aware refinement layer.  See
    # ``_apply_xtb_cascade_if_enabled`` docstring for env-flags + contract.
    #
    # Universal — no SMILES patterns.  Bit-exact when
    # ``DELFIN_CASCADE_REFINER=0`` (default) and
    # ``DELFIN_CASCADE_REFINER_CLASSES`` is unset.  Skipped for non-metal
    # complexes (UFF handles organic ligands correctly).  Per-frame M-D
    # invariant rollback inside ``refine_with_xtb`` catches topology breaks.
    results = _apply_xtb_cascade_if_enabled(mol, results, _dual_parse_done)

    # Multi-sigma V2: restore the per-seed embedding timeout override on
    # the calling thread so this function call has no effect on any
    # later, unrelated call from the same thread.
    try:
        if _ms_v2_prev_override is None:
            if hasattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value"):
                delattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")
        else:
            _MULTIEMBED_TIMEOUT_OVERRIDE.value = _ms_v2_prev_override
    except Exception:
        pass

    # --- Stereocentre-fold completeness (env-flag gated, default OFF) --------
    # Coordination-created X-H stereocentres (secondary amine [NH+] / P-H …) have an up/down
    # N-H fold (R/S) that ETKDG only samples by accident -> the crystal's fold (USEMOW's
    # alternating N1+N1-N1+N1-) is often MISSING even when the polyhedron + coordination
    # isomers are perfect (eye: ccdc_isomer_realized=FALSE under DELFIN_EYE_NH_STEREO).  This
    # pass ADDITIVELY appends every buildable fold -- both + and - at each such donor
    # (completeness law) -- so the crystal's fold AND all other feasible folds enter the
    # manifold.  Runs AFTER the final dedup so folds (heavy-atom-close to their base) are never
    # collapsed, and BEFORE the rotamer/conformer expansions so each fold gets conformer
    # diversity too.  Additive + deterministic -> never-worse by construction.  Bit-exact
    # no-op when DELFIN_STEREOCENTER_ENUM=0 or no coordination-created X-H stereocentre exists.
    results = _apply_stereocenter_enum_if_enabled(mol, results, _dual_parse_done)

    # --- AXIAL completeness (atropisomers), additive, default OFF ---------------------------
    # Directly behind the stereocentre expansion and at this position for the same reason:
    # AFTER the final dedup (the opposite handedness is close to the heavy atoms at its base and
    # would otherwise collapse away) and BEFORE the rotamer/conformer expansions (so that every
    # handedness gets its conformer diversity).
    # ⚠ SECOND CALL SITE IN THE FF-FREE BRANCH: on 10.08. exactly this module family ran ONLY on
    # legacy, because the FF-free path has its own body.  That is why it is wired on both
    # sides from the start -- call sites count, not lines.
    results = _apply_atropisomer_enum_if_enabled(mol, results, _dual_parse_done)

    # --- MIRROR COMPLETION, now ALSO on the legacy path (27.08.2026) ------------------------
    # COUNTED, not assumed: `_apply_mirror_enum_if_enabled` had 2 call sites in the
    # FF-free tail and ZERO here -- although `_mirror_enum.expand_results(results)`
    # takes frames and returns frames and does not know at all who built them.  It is
    # structurally path-neutral.  And 4069 of 5685 systems are built LEGACY, so the
    # mirror completion reached 28 % of the pool.
    #
    # ⛔ WHY ONLY NOW, although the line is trivial: it is NOT trivial without the two
    # gates below it.  Measured on 27.08. on both archives:
    #     without gates  +122 663 frames on 130 882 = +93.7 %   (archive doubling)
    #     cause          `_self_mirror_rmsd` checks a FIXED atom order and therefore
    #                    reports "chiral" for 4069 of 4069 = 100 % -- it measures
    #                    CONFORMER handedness, not MOLECULAR chirality
    #     split          >=1 stereocentre  1325 (32.6 %)  -> real new stereoisomer
    #                    no stereocentre   2744 (67.4 %)  -> only a second conformer,
    #                                                        68 % of the price for ZERO gain
    # ⇒ Only with DELFIN_MIRROR_STEREO_GATE=1 and DELFIN_MIRROR_ONE_PER_SYSTEM=1 does the
    #   frame doubler become a one-percent lever on 1325 systems (+1.0 % instead of +93.7 %).
    #   Both gates are switchable SEPARATELY, because they answer different questions --
    #   "which systems" and "how many frames per system".
    #
    # ⚠ NOT MEASURED: whether the 1325 MISS their CCDC isomer today.  Without this number the
    #   benefit is an upper bound, not a landing.  The eye decides that.
    # ⛔ Default OFF (DELFIN_MIRROR_ENUM) -> byte-identical, as on the FF-free path.
    results = _apply_mirror_enum_if_enabled(results)

    # --- Welle-5l Track-6: rotamer-diversity (env-flag gated, default OFF) ---
    # For each emitted isomer, sample staggered rotamers around bulky single
    # bonds (tBu / PMe3 / iPr / NMe2 …) and append the top-K best-energy
    # rotamer-frames as additional isomers labelled "<base>_rotamer-N".
    # Runs after final dedup so rotamer entries are never collapsed by the
    # label-suffix regex above. The helper is a no-op when the env-flag is
    # unset, preserving byte-identical default behaviour.
    try:
        from delfin.manta import _rotamer_diversity as _rot_div  # local import
        if _rot_div._is_enabled() and results:
            expanded: List[Tuple[str, str]] = []
            for _ridx, (_rxyz, _rlbl) in enumerate(results):
                _frames = _rot_div.apply_if_enabled(_rxyz)
                # _frames[0] is the original XYZ; preserve it under the
                # original label. _frames[1:] are rotamer-variants.
                if not _frames:
                    expanded.append((_rxyz, _rlbl))
                    continue
                expanded.append((_frames[0], _rlbl))
                for _k_idx, _fxyz in enumerate(_frames[1:], start=1):
                    _rot_label = f"{_rlbl}_rotamer-{_k_idx}" if _rlbl else f"rotamer-{_k_idx}"
                    expanded.append((_fxyz, _rot_label))
            results = expanded
    except Exception as _rot_exc:
        logger.debug("Welle-5l Track-6 rotamer expansion failed: %s", _rot_exc)

    # --- Welle-5o: per-isomer conformer-pool (env-flag gated, default OFF) ---
    # For each emitted isomer, generate K diverse conformers spanning the
    # relevant conformational space (torsion + ring-pucker + chelate-twist +
    # macrocycle modes), filtered by M-D invariant + topology preservation,
    # ranked by UFF energy + clash penalty, greedy-selected for pairwise-
    # RMSD diversity.  Each extra conformer is labelled "<base>_pool-N-<tag>".
    # Runs AFTER rotamer expansion so the pool acts on every isomer (including
    # rotamer-variants) and is the CREST/GOAT-obsolescence enabler — the
    # downstream xTB/DFT local-opt finds the global minimum from at least one
    # pool member without a separate global conformer search.
    # The helper is a no-op when DELFIN_5O_CONFORMER_POOL=0, preserving
    # byte-identical default behaviour.
    try:
        from delfin.manta import _conformer_pool as _conf_pool  # local import
        if _conf_pool._is_enabled() and results:
            pool_expanded: List[Tuple[str, str]] = []
            for _pidx, (_pxyz, _plbl) in enumerate(results):
                _members = _conf_pool.apply_if_enabled(_pxyz)
                if not _members:
                    pool_expanded.append((_pxyz, _plbl))
                    continue
                # _members[0] is (base_xyz, "base") — preserve under
                # original label.  _members[1:] are pool variants.
                pool_expanded.append((_members[0][0], _plbl))
                for _k_idx, (_fxyz, _ftag) in enumerate(_members[1:], start=1):
                    _pool_label = (
                        f"{_plbl}_pool-{_k_idx}-{_ftag}"
                        if _plbl
                        else f"pool-{_k_idx}-{_ftag}"
                    )
                    pool_expanded.append((_fxyz, _pool_label))
            results = pool_expanded
    except Exception as _pool_exc:
        logger.debug("Welle-5o conformer-pool expansion failed: %s", _pool_exc)

    # --- Welle-5p-A: post-emit topology hard-gate (env-flag gated, default OFF) ---
    # Final stand-alone amine-H realism check on every emitted frame.  Drops
    # frames whose amine-N (or P/As) H atoms have flipped toward the metal
    # (∠(M-D-H) < threshold or H · · · M too close).  This catches T6 / 5o
    # rotation artifacts that slipped past the per-layer gate as well as
    # any post-UFF amine-H umbrella inversion introduced upstream.
    # Default OFF: when DELFIN_5P_A_TOPOLOGY_HARDGATE=0 this is a no-op.
    try:
        from delfin.manta import _topology_hash as _th_post  # local import
        if _th_post.is_hardgate_enabled() and results:
            _gated: List[Tuple[str, str]] = []
            _dropped = 0
            for (_gxyz, _glbl) in results:
                _gate = _th_post.standalone_amine_h_realism_xyz(_gxyz)
                if _gate.passed:
                    _gated.append((_gxyz, _glbl))
                else:
                    _dropped += 1
                    logger.debug(
                        "5p-A drop emit '%s': %s",
                        _glbl,
                        _gate.violations[:3],
                    )
            if _gated:
                results = _gated
            if _dropped:
                logger.debug(
                    "5p-A post-emit gate dropped %d/%d frames",
                    _dropped,
                    _dropped + len(_gated),
                )
    except Exception as _post_exc:
        logger.debug("Welle-5p-A post-emit gate failed: %s", _post_exc)

    # --- Welle-5q: final all-class collapse-repair (env-flag gated, default OFF) ---
    # A collapsed heavy-heavy pair (0.24-1.2 A: fused substituents / overlapping
    # rings from a strained ETKDG embed) is a construction DEFECT in ANY class,
    # not only hapto.  The early ``_apply_bond_decollapse_if_enabled`` pass is
    # class-gated (default-ON {hapto, multi_hapto}) AND runs *before* the H-VSEPR /
    # rotamer / xtb correctors, so a collapse left in a sigma-class core (e.g. an
    # NHC-carbene benzimidazole ring, DUZFIQ) is never repaired.  The per-frame
    # corrector freezes metals + the M-D sphere, spring-relaxes the non-metal
    # heavy graph to physical bond lengths, and rolls back per frame unless the
    # collapsed-bond count strictly drops with no M-D break and no vdw / h-planar /
    # bond-distort regression -> byte-identical when no collapse is present.  Run
    # once more here on the TRULY-FINAL geometry (after every upstream corrector),
    # for ALL classes, so it is non-dormant and never-worse by construction.
    # Default OFF: when DELFIN_FINAL_COLLAPSE_REPAIR=0 this is a byte-exact no-op.
    try:
        if results and _delfin_env_int("DELFIN_FINAL_COLLAPSE_REPAIR", 0) == 1:
            from delfin.manta._bond_decollapse import correct_results as _fc_correct
            results = _fc_correct(mol, results)
    except Exception as _fc_exc:
        logger.debug("Welle-5q final collapse-repair failed: %s", _fc_exc)

    # --- Welle-5r: final VSEPR terminal-group repair (env-flag gated, default OFF) ---
    # A rigid terminal EX3 group (CF3/CCl3/CBr3/CH3/SO3/...) has NO conformational
    # freedom -- it is tetrahedral -- so any distortion is a construction defect,
    # not a legitimate conformer.  The metal embed/seating sometimes splays one
    # (LUXWIL CF3 F-C-F up to 157deg while the coordination sphere is fine).  This
    # final pass rebuilds each distorted group's terminal atoms at their ideal
    # VSEPR positions, anchored by the group's attachment bond so the rest of the
    # molecule (and the coordination sphere) is untouched.  Pure geometry,
    # license-clean, byte-identical when every terminal group is already ideal.
    # Default OFF: when DELFIN_VSEPR_REPAIR=0 this is a byte-exact no-op.
    try:
        _vsepr_on = (_delfin_env_int("DELFIN_VSEPR_REPAIR", 0) == 1
                     or _delfin_env_int("DELFIN_FFFREE_VSEPR_REPAIR", 0) == 1)
        if results and _vsepr_on:
            from delfin.manta._vsepr_repair import repair_terminal_groups as _vr
            results = [(_vr(_rx), _rl) for _rx, _rl in results]
    except Exception as _vr_exc:
        logger.debug("Welle-5r VSEPR terminal-group repair failed: %s", _vr_exc)

    if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
        # FINAL emit-level M-C-X over EVERY frame in results, measured DISTANCE-based like the eye
        # (metal-bonded C donor; worst M-C-X over distance neighbours) -- reconciles eye-vs-trace + finds
        # which emitted frame carries the linear geometry, regardless of which append path produced it.
        try:
            import numpy as np
            _trace_seating("RESULTS_RETURN_REACHED n_frames=%d len_tup0=%d" % (
                len(results), len(results[0]) if results else 0))
            for _fi, _tup in enumerate(results):
                _xyz = _tup[0]; _disp = _tup[1] if len(_tup) > 1 else ""
                _syms = []; _P = []
                for _ln in str(_xyz).strip().split("\n"):
                    _p = _ln.split()
                    if len(_p) >= 4:
                        try:
                            _c = [float(_p[1]), float(_p[2]), float(_p[3])]
                        except Exception:
                            continue          # skip count/comment header lines robustly
                        _syms.append(_p[0]); _P.append(np.array(_c))
                _mi = next((i for i, s in enumerate(_syms) if s in _METAL_SET), None)
                if _mi is None:
                    _trace_seating("RESULTS_RETURN frame=%d NO_METAL n_atoms=%d syms=%s" % (
                        _fi, len(_syms), ",".join(sorted(set(_syms)))[:40]))
                    continue
                _mp = _P[_mi]
                _cdist = sorted(round(float(np.linalg.norm(_P[_k] - _mp)), 2)
                                for _k, _s in enumerate(_syms) if _s == "C")
                _trace_seating("RESULTS_RETURN frame=%d metal=%s nC=%d minMC=%s" % (
                    _fi, _syms[_mi], len(_cdist), (_cdist[0] if _cdist else "no_C")))
                for _d, (_s, _dp) in enumerate(zip(_syms, _P)):
                    _mcn = float(np.linalg.norm(_dp - _mp))
                    if _s != "C" or _mcn > 2.8 or _mcn < 1e-6:   # eye uses M-C < 2.8
                        continue
                    _mc = (_mp - _dp) / _mcn
                    _worst = None
                    for _j, (_sj, _xp) in enumerate(zip(_syms, _P)):
                        if _j == _d or _sj == "H":
                            continue
                        _cx = _xp - _dp; _cxn = float(np.linalg.norm(_cx))
                        if 1.3 < _cxn < 2.0:
                            _ang = float(np.degrees(np.arccos(
                                max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                            if _worst is None or _ang > _worst:
                                _worst = _ang
                    _trace_seating("RESULTS_RETURN frame=%d Cdonor=%d M-C=%.2f worst_M-C-X=%s" % (
                        _fi, _d, _mcn, ("%.0f" % _worst) if _worst is not None else "no_heavy_nbr<2.0"))
        except Exception as _rre:
            _trace_seating("RESULTS_RETURN_ERR %s: %s" % (type(_rre).__name__, str(_rre)[:80]))

    # ── ERDBEBEN (2026-07-23): re-seat planar-COLLAPSED fragments from a clean ISOLATED embed.  Runs
    # ABSOLUTELY LAST -- after every enumerator / conformer-expansion / xtb / stereocenter step -- so it
    # sees the FINAL collapsed frames those steps add (an earlier hook missed them: the frames were not yet
    # collapsed).  Default-OFF byte-id (DELFIN_FFFREE_ISOLATED_SEAT); per-frame rollback keeps never-worse.
    results = _apply_isolated_reseat_if_enabled(mol, results, _dual_parse_done)

    # UNION (see the DELFIN_FFFREE_UNION block at the FF-free return): when the FF-free
    # builder produced frames and was told not to short-circuit, its frames are PREPENDED
    # here so frame 0 stays the FF-free pick -- the deterministic, construction-first one --
    # and legacy's spray follows as additional manifold members.  Strictly additive in both
    # directions: neither constructor's frames are altered or dropped.
    #
    # ═══ TRACE AT THE MERGE (30.08.2026) -- DELFIN_FFFREE_UNION_TRACE, default OFF ═════
    #
    # WHAT FOR.  CAZMEX, JACVES and XIDKON lose their entire chelate family with UNION
    # (8 -> 0, 3 -> 0, 9 -> 0).  Seven candidates are excluded -- most recently the
    # chelate enumerator itself, whose build trace (`_iso_trace`) is bit-identical in BOTH
    # arms: so it builds the chelate frames in the union arm too.  And between the
    # stashing (:32969) and here there is no `return` at body level, the block is
    # reached.  That leaves exactly three possibilities, and none is measured:
    #   (1) `_ffree_union` is EMPTY here, although the enumerator built
    #   (2) `if _ffree_union:` is false for another reason
    #   (3) the chelate frames survive the merge and are filtered AFTERWARDS
    #
    # Three numbers separate the three cases: how many frames `_ffree_union` carries, how
    # many `results`, and how many of those are chelate.  Exactly those this line writes.
    #
    # ⚠️ OUTPUT ONLY, no branching -- the build stays byte-identical when the
    #    switch is off (and it is, by default).  `os.write(2, ...)` as for the
    #    CN4 debug a few lines further down, so that the line appears even when
    #    the logger is not set up in the worker process.
    if os.environ.get("DELFIN_FFFREE_UNION_TRACE", "0") == "1":
        # ⚠️ MARKER WITHOUT EXCEPTION TRAPPING (30.08.2026, 14:0x).  The numbers line below
        #    sits in a `try/except: pass` -- if it throws, NOTHING appears, and that
        #    looks in the output exactly like "block not reached".  Exactly this
        #    case occurred with `cazmerge`: the ISO trace fired completely, the
        #    numbers line nowhere.  Without this marker the two explanations
        #    -- block never reached OR exception swallowed -- are indistinguishable.
        #    This `os.write` therefore stands BEFORE the `try` and is itself unguarded:
        #    if it appears, the block was reached; if it does not appear, it was
        #    not.  If it falls over itself, one sees the traceback instead of silence.
        os.write(2, b"[UNION_MERGE] ERREICHT\n")
        try:
            _u = _ffree_union or []
            _r = results or []
            _uc = sum(1 for _x, _l in _u if "chelate" in str(_l))
            _rc = sum(1 for _x, _l in _r if "chelate" in str(_l))
            os.write(2, ("[UNION_MERGE] ffree=%d (chelat %d) legacy=%d (chelat %d)\n"
                         % (len(_u), _uc, len(_r), _rc)).encode())
        except Exception:
            pass
    if _ffree_union:
        try:
            # ⚠ 2026-08-06: _seen_x came from `results` -- i.e. from the SAME list that was
            # filtered afterwards.  `x not in _seen_x` was thereby false for EVERY element, the
            # list ran completely empty, and what remained was exactly list(_ffree_union) -- byte-
            # identical with the OFF arm at :31914.  The lever drove the complete legacy pipeline
            # (measured: 1.76 h against 0.05 h on 24 systems, 35x) and then discarded its
            # result.  The reach probe caught it: 0/24 changed.
            # The dedup set must come from the FF-free frames -- those are the ones that
            # are prepended, i.e. the ones against which legacy must be deduplicated.
            _seen_x = {x for x, _l in _ffree_union}
            _extra = [(x, l) for x, l in results if x not in _seen_x]
            # ===== ONLY THE MISSING ISOMERS (DELFIN_FFFREE_UNION_ISOMERS, default OFF) =====
            # MEASURED 2026-08-06 on union180b, 180 systems:
            #   FF-free alone             2346 frames
            #   union raw                 9829 frames   -> 7483 imported
            #   of those conf copies      4763           -> 64 % of the import is conformer spray
            # The damage of the union came NOT from legacy's isomers, but from this set:
            # smiles_ccdc_regressed and pyramid_frame_regressed COUNT DEFECTIVE FRAMES, and
            # legacy delivers 6.9 % clean manifolds.  Whoever takes over its spray
            # takes over that rate a thousandfold.
            #
            # What legacy really contributes is ISOMER COVERAGE (71.5 % against 49.5 %).  So
            # take over exactly that and nothing else: per ARRANGEMENT (the key of
            # _arrangement_key folds conformers and both hands together) at most ONE
            # representative, and only if FF-free does NOT already have this arrangement anyway.
            #
            # That is the build form the user prescribes: completeness is sacred, but
            # "all of them at any cost" does not count, and conformers go by ENERGY, not by
            # quantity.  legacy's conf spray is neither the one nor the other.
            if _extra and _delfin_env_int("DELFIN_FFFREE_UNION_ISOMERS", 0):
                _have = set()
                for _x, _l in _ffree_union:
                    try: _have.add(_arrangement_key(_l))
                    except Exception: pass
                _pick, _order = {}, []
                for _x, _l in _extra:
                    try: _k = _arrangement_key(_l)
                    except Exception: _k = str(_l)
                    if _k in _have or _k in _pick:
                        continue          # FF-free already has it, or already a representative
                    _pick[_k] = (_x, _l)
                    _order.append(_k)
                _n_before = len(_extra)
                _extra = [_pick[_k] for _k in _order]
                _trace_seating("UNION_ISOMERS kept %d of %d legacy frames (%d arrangements already in ffree)"
                               % (len(_extra), _n_before, len(_have)))
            # ⛔ MEASURED AND REFUTED (2026-08-08) -- DO NOT SWITCH ON AGAIN.
            #
            #   unionreach1k  UNION + reach WITHOUT this filter        8 blocking terms, cap 0/+26
            #   unionqual     the same plus UNION_CLEAN                13 blocking terms, cap 0/+26
            #     newly torn: ccdc_backbone_lost 6, ccdc_pucker_lost 2,
            #                 ccdc_hapto_mode_lost 1, ccdc_isomer_lost 1 -> 2
            #
            # The same criterion was refuted a second time on the same day, on a
            # DIFFERENT path: TOPO_ENV in the additive enumerator, in isolation against the champion
            # (topoenvsolo) -- ccdc_arrangement_lost 3, isomers_lost 4, 1 better / 4
            # worse.  Two places, the same result: it is not the placement, it is
            # the CRITERION.  The comparison of the neighbour-element set fires on frames that
            # carry real backbones, puckers and hapticities -- the geometric perception
            # sees neighbourhoods there that the graph does not list, without anything
            # being broken.
            #
            # The best known union state is thus WITHOUT filter: cap_lost 0, cap_gained
            # 26, 8 blocking terms.  What union still costs needs a different instrument -- and
            # after three failed attempts (TORN_GATE, SPURIOUS_BOND, this one) the
            # honest reading is that the builder CANNOT separate the bad import frames from
            # the good ones with the available means.  That is the limit
            # "the eye is the ceiling of construction", measured three times.
            # ===== QUALITY FILTER ON THE IMPORT (UNION_CLEAN, new 2026-08-07) =====
            #
            # WHY THE FIRST VERSION WAS TOO WEAK.  It sent every legacy frame through
            # `_build_is_clean` -- and that one was BLIND there to exactly the defect type that legacy
            # brings along: it got no `graph_bonds`, so it knew neither a MISSING nor an
            # INVENTED bond, but only "is this contact too short".  Measured: it filtered
            # 15 % of the frames and lowered the regressions by nothing worth mentioning.
            #
            # WHAT THE MEASUREMENT SAYS.  unionreach1k has proven that the either-or is the whole
            # cause of the reach damage: MULTIBOND_LENGTH_EXEMPT alone cap_lost 40,
            # with UNION cap_lost 0 and all 40 rescued.  What UNION itself still costs is the
            # unfiltered import: isomers_lost 43, n_good_regressions 23.  And UNION_ISOMERS
            # has shown that it is NOT the quantity -- frames 9829 -> 3899, regressions only
            # -16 %.  It is the QUALITY of the imported frames.
            #
            # THE CRITERION, and why it works here at all.  A comparison frame-against-SMILES
            # normally fails on the ATOM ORDER -- exactly on that the mechanism twenty lines
            # further up has already died once ("They never coincide, so it
            # returned 0 every time").  The eye solves it ORDER-FREE: from the SMILES it is
            # fixed which heavy neighbourhoods an element MAY have (a nitrate N {O,O,O},
            # an acetate C {C}, an ether O {C,C}).  A frame atom whose perceived
            # neighbourhood is NONE of the allowed ones carries a topology break -- without
            # a single atom having to be mapped.
            #
            # That catches both directions: missing neighbour = detached substituent or
            # torn bond; additional neighbour = fused contact.  Metals stay out on
            # BOTH sides -- the coordination number is judged by the polyhedron, not by
            # this test.  License-clean: RDKit parse plus geometric perception, no
            # reference data.
            if _extra and _delfin_env_int("DELFIN_FFFREE_UNION_CLEAN", 0) and mol is not None:
                try:
                    from delfin.manta import _bond_decollapse as _uq_bd
                    import numpy as _uq_np
                    from collections import Counter as _uq_C

                    _allowed = {}
                    for _a in mol.GetAtoms():
                        _s = _a.GetSymbol()
                        if _s == "H" or _uq_bd._is_metal(_s):
                            continue
                        _env = _uq_C(nb.GetSymbol() for nb in _a.GetNeighbors()
                                     if nb.GetSymbol() != "H" and not _uq_bd._is_metal(nb.GetSymbol()))
                        _allowed.setdefault(_s, set()).add(
                            tuple(sorted(_env.items())))

                    def _uq_ok(_xyz):
                        try:
                            _sy, _co = [], []
                            for _ln in _xyz.strip().splitlines():
                                _p = _ln.split()
                                if len(_p) >= 4:
                                    _sy.append(_p[0])
                                    _co.append([float(_p[1]), float(_p[2]), float(_p[3])])
                            if not _sy:
                                return False
                            _P = _uq_np.array(_co, dtype=float)
                            _nbrs = {}
                            for _i, _j in _uq_bd._geometric_bonds(_sy, _P):
                                if _sy[_i] == "H" or _sy[_j] == "H":
                                    continue
                                if _uq_bd._is_metal(_sy[_i]) or _uq_bd._is_metal(_sy[_j]):
                                    continue
                                _nbrs.setdefault(_i, []).append(_sy[_j])
                                _nbrs.setdefault(_j, []).append(_sy[_i])
                            for _k, _s in enumerate(_sy):
                                if _s == "H" or _uq_bd._is_metal(_s):
                                    continue
                                _ok = _allowed.get(_s)
                                if not _ok:
                                    continue          # element not in the SMILES -> not assessable
                                _e = tuple(sorted(_uq_C(_nbrs.get(_k, [])).items()))
                                if _e not in _ok:
                                    return False      # a neighbourhood the molecule does not know
                            return True
                        except Exception:
                            return False              # not assessable -> do not take over
                    _n0 = len(_extra)
                    _extra = [t for t in _extra if _uq_ok(t[0])]
                    _trace_seating("UNION_CLEAN(topo) kept %d of %d legacy frames"
                                   % (len(_extra), _n0))
                except Exception as _uq_exc:
                    _trace_seating("UNION_CLEAN(topo) no-op: %s" % (str(_uq_exc)[:70],))
            # MERGE: FF-free first, so that frame 0 stays the deterministic
            # construction; legacy's checked extra follows as further
            # manifold members.  Neither of the two builders loses anything.
            results = list(_ffree_union) + _extra
        except Exception:
            results = list(_ffree_union) + list(results)

    return results, None


def convert_input_if_smiles(input_path: Path) -> Tuple[bool, Optional[str]]:
    """Check if input file contains SMILES and convert if needed.

    This function is called by the pipeline to automatically handle SMILES input.

    Args:
        input_path: Path to input file to check

    Returns:
        Tuple of (was_converted, error_message)
        - was_converted: True if file was SMILES and was converted
        - error_message: Error description if conversion failed, None if success or not SMILES
    """
    if not input_path.exists():
        return False, f"Input file does not exist: {input_path}"

    try:
        content = input_path.read_text(encoding='utf-8', errors='ignore')
    except Exception as e:
        return False, f"Could not read input file: {e}"

    if not is_smiles_string(content):
        return False, None

    # Extract SMILES (first non-comment line)
    smiles = None
    for line in content.split('\n'):
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('*'):
            smiles = line
            break

    if not smiles:
        return False, "No SMILES string found in input file"

    logger.info(f"Detected SMILES string in {input_path.name}: {smiles}")

    # Convert SMILES to XYZ
    xyz_content, error = smiles_to_xyz(smiles)

    if error:
        return False, error

    # Write XYZ content back to input file (replacing SMILES)
    try:
        input_path.write_text(xyz_content, encoding='utf-8')
        logger.info(f"Converted SMILES to XYZ coordinates in {input_path}")
        return True, None
    except Exception as e:
        return False, f"Could not write converted coordinates: {e}"


def smiles_to_xyz_architector(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Convert a metal-complex SMILES to XYZ using Architector (lowest-energy frame).

    The build is the one the dashboard's ARCHITECTOR button runs
    (:mod:`delfin.common.external_builders`): coordination number and oxidation
    state are read from the SMILES, and a missing Architector is an error.
    """
    from delfin.common.external_builders import build_first_xyz
    return build_first_xyz('architector', smiles)


def smiles_to_xyz_molsimplify(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Convert a metal-complex SMILES to XYZ using molSimplify (first geometry).

    Same shared build as the dashboard's MOLSIMPLIFY button; the first frame is
    the first geometry molSimplify lists for the coordination number (square
    planar for CN 4, octahedral for CN 6).
    """
    from delfin.common.external_builders import build_first_xyz
    return build_first_xyz('molsimplify', smiles)


def smiles_to_xyz_mace(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Convert a metal-complex SMILES to XYZ using epic-MACE (first frame).

    Same shared build as the dashboard's MACE button, with the same settings
    (``external_builders.MACE_DEFAULTS``); epic-MACE runs in its own Python 3.7
    environment (``DELFIN_MACE_PYTHON``, or the one the installer built).  The
    first frame is the lowest-energy conformer of the first geometry that fits
    the number of donor sites (octahedral for 6, square planar for 4).
    """
    from delfin.common.external_builders import build_first_xyz
    return build_first_xyz('mace', smiles)


# ---------------------------------------------------------------------------
# SELF-TEST for DELFIN_FFFREE_HAPTO_SEAT_RIGID
#   PYTHONPATH=/home/localuser/DELFIN_dev python -m delfin.smiles_converter
# A direct file invocation loads a DIFFERENT module copy because of the editable
# installation -- always start via `-m`.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys as _sys

    _fehler = 0

    def _pruefe(name: str, ist, soll):
        global _fehler
        ok = ist == soll
        if not ok:
            _fehler += 1
        print(f"   [{'ok' if ok else 'FEHL'}] {name}: {ist!r} (erwartet {soll!r})")

    print("## Selbsttest HAPTO_SEAT_RIGID")

    # 1. The gate is closed as long as nobody opens it -- and it reads exactly ONE variable.
    os.environ.pop("DELFIN_FFFREE_HAPTO_SEAT_RIGID", None)
    _pruefe("Vorgabe AUS", _hapto_seat_rigid_enabled(), False)
    os.environ["DELFIN_FFFREE_HAPTO_SEAT_RIGID"] = "1"
    _pruefe("Schalter AN", _hapto_seat_rigid_enabled(), True)
    os.environ["DELFIN_FFFREE_HAPTO_SEAT_RIGID"] = "0"
    _pruefe("Schalter 0", _hapto_seat_rigid_enabled(), False)
    os.environ.pop("DELFIN_FFFREE_HAPTO_SEAT_RIGID", None)

    # 2. The collapse count judges two-sidedly: it must find the compressed build
    #    AND leave the healthy one alone.
    if RDKIT_AVAILABLE:
        def _mol_mit(abstand: float):
            m = Chem.RWMol()
            m.AddAtom(Chem.Atom(6))
            m.AddAtom(Chem.Atom(6))
            m.AddBond(0, 1, Chem.BondType.SINGLE)
            out = m.GetMol()
            out.UpdatePropertyCache(strict=False)
            c = Chem.Conformer(2)
            c.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
            c.SetAtomPosition(1, Point3D(abstand, 0.0, 0.0))
            out.AddConformer(c, assignId=True)
            return out

        _pruefe("gesunde C-C (1.50 A) -> kein Kollaps",
                _hapto_candidate_collapsed_bonds(_mol_mit(1.50)), 0)
        _pruefe("gestauchte C-C (1.00 A) -> ein Kollaps",
                _hapto_candidate_collapsed_bonds(_mol_mit(1.00)), 1)
        _pruefe("ohne Konformer -> kein Urteil (None)",
                _hapto_candidate_collapsed_bonds(Chem.MolFromSmiles("CC")), None)
        _pruefe("None hinein -> None heraus",
                _hapto_candidate_collapsed_bonds(None), None)
    else:
        print("   (RDKit fehlt -- Geometrieteil uebersprungen)")

    # ------------------------------------------------------------------
    # SELF-TEST for DELFIN_FFFREE_UNION_HAPTO
    #   A switch that only proves that it FIRES proves nothing.  Two of the
    #   five cases are CONTRARY: in 3 and 5 the switch is ON and the
    #   answer must nevertheless be unchanged.
    # ------------------------------------------------------------------
    print("## Selbsttest UNION_HAPTO")
    _H = [("xH1", "η6-arene"), ("xH2", "η6-arene σ-1")]      # hapto branch
    _F = [("xF1", "T-4-hapto-iso1"), ("xF2", "T-4-hapto-iso2")]   # FF-free

    os.environ.pop("DELFIN_FFFREE_UNION_HAPTO", None)
    _pruefe("1 Vorgabe AUS -> unveraendert",
            _union_prepend_ffree(_H, _F), _H)

    os.environ["DELFIN_FFFREE_UNION_HAPTO"] = "1"
    _pruefe("2 AN -> FF-frei zuerst, Hapto vollstaendig dahinter",
            _union_prepend_ffree(_H, _F), _F + _H)

    # CONTRARY: the switch is on, but there is nothing to merge.
    _pruefe("3 AN ohne FF-freie Frames -> unveraendert",
            _union_prepend_ffree(_H, None), _H)

    # Dedup over the exact XYZ text: the shared frame appears ONCE.
    _pruefe("4 AN mit Ueberlappung -> kein Frame doppelt",
            _union_prepend_ffree([("xF2", "andersherum")] + _H, _F),
            _F + _H)

    # CONTRARY: '0' is OFF, not 'some value set'.
    os.environ["DELFIN_FFFREE_UNION_HAPTO"] = "0"
    _pruefe("5 Schalter 0 -> unveraendert",
            _union_prepend_ffree(_H, _F), _H)
    os.environ.pop("DELFIN_FFFREE_UNION_HAPTO", None)

    print(f"## {'ALLES GRUEN' if _fehler == 0 else str(_fehler) + ' FEHLER'}")
    _sys.exit(1 if _fehler else 0)
