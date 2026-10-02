"""Topological isomer construction: the template builder, the graph topology verifier, the general topological isomer generator, hapto-sigma isomers, pucker variants and all-trans arrangements of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import concurrent.futures
import math
import os
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

from delfin.common.logging import get_logger
from delfin.manta.chelate_templates import (
    _best_chelate_conformer_coords,
    _build_topology_xyz_from_scratch,
    _embed_fragment_procrustes,
)
from delfin.manta.conformer_io import (
    _mol_to_xyz_conformer,
    _xyz_to_rdkit_conformer,
)
from delfin.manta.conformer_pools import (
    _cp_basin_6ring,
    _frag_xyz_collapsed,
    _rank_template_conformers,
    _ring_canonical_snap_z,
)
from delfin.manta.converter_flags import (
    DELFIN_CHELATE_RANK_VARIANTS,
    DELFIN_MAX_PROCESS_WORKERS,
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
    DELFIN_TOPO_TEMPLATE_TOP_K,
    _delfin_env_float,
    _delfin_env_int,
    _deterministic_mode,
    _every_append_gate_enabled,
    _preferred_cn4_for,
    _resolve_quality_profile,
    _trace_seating,
)
from delfin.manta.coordination_enumerator import (
    _enumerate_topological_isomers,
)
from delfin.manta.geometry_quality import (
    _donor_h_points_at_metal,
    _ideal_polyhedron_angle_dev_per_metal,
)
from delfin.manta.hapto_detect import (
    _classify_complex_class,
    _find_hapto_groups,
)
from delfin.manta.isomer_labels import (
    _TOPO_GEOMETRY_VECTORS,
    _TOPO_TRANS_POSITIONS,
    _chelate_backbone_max_reach,
    _chelate_pairs,
    _classify_isomer_label,
    _compute_coordination_fingerprint,
    _donor_type_map,
    _extract_helicity_suffix,
    _find_bridging_donors,
    _label_from_canonical_form,
)
from delfin.manta.ligand_placement import (
    _align_and_orient_ligands,
    _build_multimetal_scaffold,
    _snap_aromatic_rings_in_xyz,
    _snap_bridging_donors_to_compromise,
)
from delfin.manta.ml_tables import (
    Chem,
    OPENBABEL_AVAILABLE,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_METAL_BOND_LENGTHS,
    _METAL_SET,
    _PREFERRED_CN5_GEOMETRY,
    _PREFERRED_CN6_GEOMETRY,
    _get_ml_bond_length,
)
from delfin.manta.mol_prep import (
    _prepare_mol_for_embedding,
    _rescale_metal_donor_distances,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel,
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.pre_uff_snap import (
    _D8_SQ_ISO_METALS,
    _cn5_enum_complete_enabled,
)
from delfin.manta.topology_checks import (
    _flatten_sp2_atoms_xyz,
    _has_atom_clash,
    _has_pi_ring_nonplanarity,
    _has_severe_covalent_distortion,
    _has_unphysical_metal_nonbonded_contact,
    _has_unphysical_oco_geometry,
    _metal_donor_distances_realistic,
)
from delfin.manta.uff_constraints import (
    _build_coordination_constraints_from_xyz,
)

logger = get_logger("delfin.smiles_converter")


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
