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

from delfin.manta.chelate_templates import (
    _ITER84_SIGMA_CAPS_OVERRIDE,
    _SIGMA_CHELATE_CAPS_123A,
    _best_chelate_conformer_coords,
    _build_topology_template_mol,
    _build_topology_xyz_from_scratch,
    _chelate_conformer_candidates,
    _embed_fragment_procrustes,
)  # noqa: F401
# The module itself, for the one attribute _smiles_to_xyz_isomers_impl writes
# (_ITER84_SIGMA_CAPS_OVERRIDE); see the comment at that write.
from delfin.manta import chelate_templates as _chelate_templates

from delfin.manta.topo_isomers import (
    _build_topology_xyz,
    _build_topology_xyz_from_template,
    _emit_all_trans_by_type_arrangements,
    _emit_chelate_pucker_variants,
    _emit_nonmetal_ring_pucker_variants,
    _enumerate_hapto_sigma_isomers,
    _generate_topological_isomers,
    _verify_topology_from_graph,
)  # noqa: F401

from delfin.manta.binding_modes import (
    _find_linkage_alternatives,
    _generate_alternative_binding_modes,
    _generate_linkage_isomers,
    _rewire_linkage,
)  # noqa: F401

from delfin.manta.result_filters import (
    _apply_pi_coplanar_final,
    _apply_pi_inplane_final,
    _arrangement_key,
    _carbonyl_fix_filter,
    _clean_gate_filter,
    _conf_complete_filter,
    _coord_chirality_sign,
    _coord_integrity_filter,
    _coord_sphere_chirality_rmsd,
    _coord_sphere_donors,
    _enantiomer_mirror_filter,
    _filter_nonfinite_isomers,
    _gfnff_ensemble_rank_filter,
    _hapto_declash_filter,
    _heavy_dist_fp,
    _mirror_symmetrize_xyz,
    _mirror_xyz_coords,
    _permute_dedup_filter,
    _pointgroup_symmetrize_xyz,
    _proper_kabsch_rmsd,
    _rank_emitted_isomers,
    _sigma_declash_filter,
    _topology_gate_filter,
)  # noqa: F401
if os.environ.get("DELFIN_FFFREE_GEOM_IDEALS_REAL", "0") == "1":
    _GEOM_IDEAL_ANGLES.update(_GEOM_IDEAL_ANGLES_REAL)


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
    # MANTA split (2026-10): the override is a module global of
    # delfin.manta.chelate_templates, where _chelate_conformer_candidates reads
    # it.  The two writes in this function set that module's attribute, which
    # is exactly what the former ``global`` statement did while both lived here.
    _chelate_templates._ITER84_SIGMA_CAPS_OVERRIDE = None

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
                _chelate_templates._ITER84_SIGMA_CAPS_OVERRIDE = _SIGMA_CHELATE_CAPS_123A
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
                # MANTA split (2026-10): _verify_topology_from_graph and every
                # function that looks it up at call time (_generate_topological_isomers,
                # _build_topology_xyz_from_template, the pucker / all-trans emitters,
                # _enumerate_hapto_sigma_isomers) live in delfin.manta.topo_isomers,
                # so the temporary patch goes onto that module, as it went onto this
                # one while they were all here.
                _gie_mod = _sys.modules["delfin.manta.topo_isomers"]
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
