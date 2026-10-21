"""Plan for delfin/smiles_converter.py (line ranges of the monolith at 96c2a60a).

Usage: python tools/split_plan_smiles_converter.py <depgraph.json> <out_dir>

Writes <out_dir>/names/<module>.txt (one leading name per line) and
<out_dir>/plan.json, the step list split_extract.py consumes.
"""
import json
import sys
from pathlib import Path

S = Path(sys.argv[2])  # output directory
S.mkdir(parents=True, exist_ok=True)
g = json.load(open(sys.argv[1]))
by_line = []
for t in g["tops"]:
    if t["names"]:
        by_line.append((t["l0"], t["names"][0]))

# (module, doc, [(lo, hi), ...], extra_names, exclude_names)
GROUPS = [
    ("ml_tables", "Element sets, covalent radii, metal-ligand / metal-metal / metal-centroid bond-length tables, preferred coordination polyhedra per metal and the polyhedron catalogue of the MANTA constructor.",
     [(25, 1673)], [],
     ["_apply_h_placement_if_enabled", "_apply_mirror_enum_if_enabled", "_apply_me_bond_snap_if_enabled",
      "_HAPTO_QUICK_PREVIEW_CACHE", "_HybridHaptoFragment", "_HybridHaptoDecomposition", "_PrimaryOrganometalModule",
      "_all_polyhedra_codes", "_preferred_cn4_for"]),
    ("hapto_detect", "SMILES parsing, metal and hapto-group detection, complex-class classification and the hapto approximation of the MANTA constructor.",
     [(6233, 6578)], ["mol_from_smiles_rdkit"], []),
    ("converter_flags", "Switches, tunables, timeouts, seed schedule and quality profiles of the MANTA constructor (every module-level env read stays at import time, every helper reads at call time).",
     [(1675, 2007), (3349, 4145)], ["_all_polyhedra_codes", "_preferred_cn4_for", "_trace_seating"], []),
    ("metal_smiles", "Metal SMILES normalisation (dative bonds, donor hydrogens, organometallic carbons) of the MANTA constructor.",
     [(5754, 6232)], [], []),
    ("conformer_io", "Open Babel conformer generation, XYZ to RDKit conformer mapping and RDKit molecule to XYZ serialisation of the MANTA constructor.",
     [(5222, 5753)], ["_fix_zero_coord_hydrogens", "_mol_to_xyz", "_mol_to_xyz_conformer"], []),
    ("embed_timeout", "Timeout-guarded RDKit embedding (single and multi-conformer) of the MANTA constructor.",
     [(4147, 4304)], [], []),
    ("mol_prep", "Molecule preparation for embedding (cached), metal-donor distance rescale and robust multi-conformer embedding of the MANTA constructor.",
     [(16690, 17225)], [], ["mol_from_smiles_rdkit"]),
    ("topology_checks", "Graph-signature topology checks, metal connectivity verification and the geometry predicates of the MANTA constructor's gates.",
     [(17226, 17454), (17471, 18219), (18947, 19735)], ["_count_xyz_clashes"], ["_xyz_passes_final_geometry_checks"]),
    ("isomer_labels", "Donor typing, viable and bridging donors, ligand fragments, coordination fingerprints, canonical polyhedron forms and isomer labels of the MANTA constructor.",
     [(21573, 22646)], ["_is_viable_donor", "_find_bridging_donors", "_ligand_fragments"], []),
    ("hapto_scaffold", "Sphere-based constructive hapto scaffold: ansa bridges, shared eta rings, centroid bias, ring planarity, hapto geometry correction and BFS propagation in the MANTA constructor.",
     [(13189, 16689)], [], []),
    ("hybrid_fragments", "Hapto fragment extraction, embedding and rigid alignment for hybrid (hapto plus sigma) complexes in the MANTA constructor.",
     [(30, 65), (6579, 7447)], [], ["_align_hybrid_fragment_to_targets", "_append_hapto_preview_xyz", "_store_hapto_preview_candidate"]),
    ("secondary_metal_modules", "Secondary-metal coordination modules of multi-metal hapto complexes: geometry fits, pose optimisation and module assembly in the MANTA constructor.",
     [(7449, 11145)], ["_align_hybrid_fragment_to_targets", "_hybrid_bond_target_length"], []),
    ("hybrid_assembly", "Assembly and refinement of hybrid hapto complexes, donor pi-coplanarity, sequential multi-metal hapto building, hapto previews and final clash resolution in the MANTA constructor.",
     [(11146, 13188)], ["_append_hapto_preview_xyz", "_store_hapto_preview_candidate"], ["_hybrid_bond_target_length"]),
    ("geometry_quality", "Geometry quality scores, the final geometry checks, metal topology enforcement and isomer upper-bound estimates of the MANTA constructor.",
     [(19736, 20432), (20908, 21572)], ["_xyz_passes_final_geometry_checks"], []),
    ("coordination_enumerator", "Orbit enumeration of topological coordination isomers (Burnside / Polya) in the MANTA constructor.",
     [(22647, 23186)], [], []),
    ("ligand_placement", "Multi-metal scaffold, aromatic ring snapping and scaling, bridging donors, lone-pair tilt and ligand alignment of the MANTA constructor.",
     [(23735, 25431)], [], []),
    ("pre_uff_snap", "Pre-UFF metal-donor snap, topology gate switches, metalloid clamps and d8 square-planar flattening of the MANTA constructor.",
     [(25432, 25887)], [], []),
    ("uff_constraints", "Pyykko radii, donor sigma geometry and the UFF coordination constraints built from a template or from an XYZ in the MANTA constructor.",
     [(37876, 39209)], [], []),
    ("openbabel_optimize", "Geometric inter-ligand clash relief, template bond orders and the Open Babel UFF optimisation of the MANTA constructor.",
     [(39210, 40211)], [], []),
    ("hapto_candidates", "Hapto candidate topology check, quality scores, RDKit UFF refinement and best-candidate selection of the MANTA constructor.",
     [(20433, 20907)], [], []),
    ("stage_hooks", "The additive post-construction stages of the MANTA constructor, each behind its own DELFIN_FFFREE_* switch and delegating to a delfin.manta module.",
     [(853, 980), (2008, 3348)], [], []),
    ("embed_strategies", "Embedding strategies of the MANTA constructor: retries, the no-valence-check path, hydrogen geometry repair, clash-aware hydrogen placement, the manual metal embed and the unsanitised fallback.",
     [(4305, 5221)], ["_smiles_to_xyz_unsanitized_fallback"], ["_count_xyz_clashes"]),
    ("single_structure", "Single-structure SMILES to XYZ conversion of the MANTA constructor: quick hapto previews, the sanitised path and the legacy path.",
     [(36447, 37875)], ["_HAPTO_QUICK_PREVIEW_CACHE", "_SMILES_TOKEN_RE", "_SMILES_ATOM_TOKEN_RE", "_pin_metal_neighbour_hydrogens"],
     ["_smiles_to_xyz_unsanitized_fallback", "_fix_zero_coord_hydrogens", "_mol_to_xyz", "_mol_to_xyz_conformer"]),
    ("conformer_pools", "Template conformer ranking, organic conformer pools, TFD dedup, ring puckers and d8 SP-4 variants of the MANTA constructor.",
     [(26457, 26628), (29427, 29774), (29976, 30211)], [], []),
    ("chelate_templates", "Chelate conformer candidates (with the sigma cap override), Procrustes fragment embedding, the from-scratch topology builder and the topology template molecule of the MANTA constructor.",
     [(17455, 17470), (23187, 23734), (25888, 26177)], ["_build_topology_template_mol"], []),
    ("topo_isomers", "Topological isomer construction: the template builder, the graph topology verifier, the general topological isomer generator, hapto-sigma isomers, pucker variants and all-trans arrangements of the MANTA constructor.",
     [(26190, 26456), (26629, 28562)],
     ["_verify_topology_from_graph", "_enumerate_hapto_sigma_isomers", "_emit_chelate_pucker_variants",
      "_emit_nonmetal_ring_pucker_variants", "_emit_all_trans_by_type_arrangements"], ["_build_topology_template_mol"]),
    ("binding_modes", "Linkage isomers and alternative binding-site exploration of the MANTA constructor.",
     [(28563, 29137)], [], ["_is_viable_donor", "_find_bridging_donors", "_ligand_fragments"]),
    ("result_filters", "The frame filters applied to the finished manifold: non-finite, ranking, declash, carbonyl, topology gate, clean gate, pi planes, coordination integrity, conformer completion, GFN-FF rank, permutation dedup and mirror enumeration.",
     [(30967, 32360)], [], []),
]

plan = []
for mod, doc, ranges, extra, exclude in GROUPS:
    names = []
    for lo, hi in ranges:
        names += [n for (l, n) in by_line if lo <= l <= hi]
    names += extra
    names = [n for n in names if n not in exclude]
    # drop module header names (logger, imports) that never move
    names = [n for n in names if n not in ("logger",)]
    seen = set()
    names = [n for n in names if not (n in seen or seen.add(n))]
    (S / "names" ).mkdir(exist_ok=True)
    (S / "names" / f"{mod}.txt").write_text("\n".join(names) + "\n")
    plan.append({"module": f"delfin/manta/{mod}.py", "doc": doc,
                 "names_file": str(S / "names" / f"{mod}.txt")})
json.dump(plan, open(S / "plan.json", "w"), indent=1)
print("groups:", len(plan))
