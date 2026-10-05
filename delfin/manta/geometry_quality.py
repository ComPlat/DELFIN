"""Geometry quality scores, the final geometry checks, metal topology enforcement and isomer upper-bound estimates of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    DELFIN_IDEAL_POLYHEDRON_MAX_DEV,
    DELFIN_SYMMETRY_WEIGHT,
    _class_conditional_flag,
    _delfin_env_int,
)
from delfin.manta.hapto_detect import (
    _find_hapto_groups,
)
from delfin.manta.isomer_labels import (
    _chelate_pairs,
    _donor_type_map,
)
from delfin.manta.ml_tables import (
    Chem,
    Point3D,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_METAL_BOND_LENGTHS,
    _METAL_SET,
)
from delfin.manta.topology_checks import (
    _has_pi_ring_nonplanarity,
    _has_severe_covalent_distortion,
    _has_unphysical_metal_nonbonded_contact,
    _has_unphysical_oco_geometry,
)

logger = get_logger("delfin.smiles_converter")


def _xyz_passes_final_geometry_checks(
    xyz_delfin: str,
    mol_template,
    skip_angle_check: bool = False,
) -> bool:
    """Final geometry sanity check for accepted XYZ outputs.

    Intended for relaxed-fallback candidates: keep only structures that still
    pass hard geometric plausibility checks when mapped back to the template
    molecular graph.

    When ``skip_angle_check`` is True, only covalent bond distortion is
    checked (not L-M-L angles). This is correct for topology-built
    structures whose donor positions are set by ideal-polyhedron vectors
    — their angles may not satisfy the sampling-oriented thresholds in
    ``_has_bad_geometry``.
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        # Ensure mol has explicit H so atom count matches XYZ (which
        # includes H atoms from the topology builder).
        mol_tmp = Chem.RWMol(mol_template)
        mol_tmp.RemoveAllConformers()
        conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_delfin)
        if conf is None:
            # XYZ has H but mol doesn't → add explicit H and retry.
            try:
                mol_h = Chem.AddHs(mol_template)
                mol_h = Chem.RWMol(mol_h)
                mol_h.RemoveAllConformers()
                conf = _xyz_to_rdkit_conformer(mol_h.GetMol(), xyz_delfin)
                if conf is None:
                    return True  # Can't validate → permissive
                cid = mol_h.AddConformer(conf, assignId=True)
                if _has_severe_covalent_distortion(mol_h.GetMol(), cid):
                    return False
                if _has_bad_geometry(mol_h.GetMol(), cid):
                    return False
                return True
            except Exception:
                return True  # Can't validate → permissive
        cid = mol_tmp.AddConformer(conf, assignId=True)
        if _has_severe_covalent_distortion(mol_tmp.GetMol(), cid):
            return False
        if not skip_angle_check and _has_bad_geometry(mol_tmp.GetMol(), cid):
            return False
        return True
    except Exception:
        return False


def _ml_distance_range(metal_symbol: str, donor_symbol: str) -> Tuple[float, float]:
    """Return (min_dist, max_dist) in Å for a metal-donor pair.

    Uses covalent radii sum with scaling factors:
    - min_dist = 0.7 × (r_M + r_D)
    - max_dist = 1.8 × (r_M + r_D)
    Falls back to (1.4, 3.5) if radii are unknown.
    """
    r_m = _COVALENT_RADII.get(metal_symbol)
    r_d = _COVALENT_RADII.get(donor_symbol)
    if r_m is not None and r_d is not None:
        r_sum = r_m + r_d
        return (0.7 * r_sum, 1.8 * r_sum)
    return (1.4, 3.5)


def _cn4_geometry_penalties(
    metal_pos,
    coord_positions: List[object],
    angles: List[float],
) -> Tuple[float, float]:
    """Return ``(tetra_pen, square_pen)`` for a 4-coordinate center."""
    tetra_pen = sum(abs(a - 109.5) for a in angles)
    square_targets = [90.0, 90.0, 90.0, 90.0, 180.0, 180.0]
    square_angle_pen = sum(
        abs(a - t) for a, t in zip(sorted(angles), square_targets)
    )

    # Square-planar candidates should also be planar. Add a donor-plane
    # penalty to distinguish flattened tetrahedral-like solutions.
    square_planarity_pen = 0.0
    try:
        import numpy as np

        pts = np.array(
            [[p.x, p.y, p.z] for p in coord_positions],
            dtype=float,
        )
        if pts.shape == (4, 3):
            centroid = pts.mean(axis=0)
            q = pts - centroid
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            normal = vh[-1]
            n_norm = float(np.linalg.norm(normal))
            if n_norm > 1e-12:
                normal = normal / n_norm
                donor_dev = np.abs(q @ normal)
                rms_donor = float(np.sqrt(np.mean(donor_dev * donor_dev)))
                mp = np.array([metal_pos.x, metal_pos.y, metal_pos.z], dtype=float)
                metal_dev = abs(float(np.dot(mp - centroid, normal)))
                square_planarity_pen = 35.0 * rms_donor + 45.0 * metal_dev
    except Exception:
        pass

    return tetra_pen, (square_angle_pen + square_planarity_pen)


def _preferred_cn4_geometry_score(
    metal_symbol: str,
    donor_symbols: List[str],
    metal_pos,
    coord_positions: List[object],
    angles: List[float],
) -> float:
    """Return a metal-aware geometry score for 4-coordinate centers."""
    tetra_pen, square_pen = _cn4_geometry_penalties(
        metal_pos,
        coord_positions,
        angles,
    )

    if metal_symbol in {'Pt', 'Pd'}:
        return min(square_pen, tetra_pen + 40.0)
    if metal_symbol == 'Au':
        return min(square_pen, tetra_pen + 28.0)
    if metal_symbol == 'Ni':
        n_n = donor_symbols.count('N')
        n_o = donor_symbols.count('O')
        n_p = donor_symbols.count('P')
        n_c = donor_symbols.count('C')
        # d8 Ni(II) with cyclometalated / N-rich donor sets strongly tends
        # toward square-planar arrangements.
        if n_c >= 1 or (n_n + n_o) >= 3 or n_p >= 2:
            return min(square_pen, tetra_pen + 22.0)
        return min(square_pen, tetra_pen + 10.0)
    if metal_symbol in {'Zn', 'Cd', 'Hg'}:
        return min(tetra_pen, square_pen + 22.0)
    if metal_symbol in {'Cu', 'Ag'}:
        return min(tetra_pen, square_pen + 10.0)

    return min(tetra_pen, square_pen)


# Ideal L-M-L angle target set per polyhedron NAME (deg).  Mirrors the per-CN best-of table inside
# _ideal_polyhedron_angle_dev_per_metal, but keyed by the geometry the frame was BUILT for (cf[0] / gn)
# so a frame can be scored vs its INTENDED polyhedron instead of the most-forgiving one.  Root motive
# (2026-07-22, AXOKED): the best-of measure lets TPR's 140-deg target absorb a distorted OCTAHEDRON
# (a rigid-scaffold Cl-Cl-ax that is ~49 deg off OH scores <30 because every angle finds SOME nearby
# ideal across OH+TPR).  Scoring a recovered arrangement vs the polyhedron it was ENUMERATED as (OH for
# an OH recovery) correctly flags it, while a genuine trig-prism (enumerated + built AS TPR) is scored
# vs TPR and kept -- so real prisms (AFAVOV/ADAVUZ) are untouched.
_GEOM_IDEAL_ANGLES: Dict[str, List[float]] = {
    'LIN': [180.0], 'TP': [120.0], 'TS': [90.0, 180.0],
    'TH': [109.47], 'SQ': [90.0, 180.0], 'SS': [90.0, 120.0, 180.0],
    'TBP': [90.0, 120.0, 180.0], 'SP': [90.0, 180.0],
    'OH': [90.0, 180.0], 'TPR': [76.0, 82.0, 140.0],
    'PBP': [72.0, 90.0, 144.0, 180.0],
    'SAP': [52.4, 73.1, 118.5, 143.1, 180.0],
    'DD': [62.2, 73.7, 117.4, 143.6, 180.0],
}


# ===== THE TARGET ANGLES AGAINST REALITY, NOT AGAINST THE BUILDER (2026-08-14) =====
# PRINCIPLE (user, 14.08.): "more realistic geometries must ALWAYS win."  Exactly that
# was violated, and in BOTH directions -- recomputed with harness/polyhedra_audit.py:
#
#   SP (square pyramid).  Here stood [90, 180].  The builder set 82.0/88.9/163.9, reality
#   is 100-105 apical-basal (VO(acac)2, [CuCl5]3-, [Ni(CN)5]3-).  With [90, 180]
#   the CORRECTED builder gets 25.0 degrees maximum deviation and the WRONG one only 16.1 --
#   the eye would have punished the right geometry and rewarded the wrong one.  With the
#   values below it turns around: 8.9 -> 0.3.  An eye that carries the old build error as
#   the ideal makes every repair unlandable.
#
#   SAP (square antiprism).  Here stand 52.4 and 180.0 degrees.  An antiprism has NO
#   antipodal pair -- 180 degrees cannot occur in it, nor can 52.4.  Since measurement is
#   against the NEAREST target value, surplus entries make the eye
#   LENIENT: a collapsed structure with a 180-degree pair would get deviation 0.
#   Computed for the equilateral antiprism: 74.9 (x16), 118.5 (x4), 141.6 (x8).
#
# ⚠ NOT touched, because not yet decided:
#   TPR  -- the prism angle depends on the height-to-width ratio, there is no "ideal"
#           prism.  Builder 70.0/90.4/131.6 against eye 76/82/140 (§2.8b).  ⚠ TPR6 is the
#           ONLY landed champion part -- nothing is changed here without a decision.
#   DD   -- also lists 180.0; whether a D2d dodecahedron has one is NOT recomputed.
#           Unchecked stays unchanged.
_GEOM_IDEAL_ANGLES_REAL: Dict[str, List[float]] = {
    'SP':  [87.0, 102.5, 155.0],      # C4v square pyramid, crystal values
    'SAP': [74.9, 118.5, 141.6],      # equilateral antiprism, computed
}


def _ideal_polyhedron_angle_dev_per_metal(
    mol, conf_id: int, only_geom: Optional[str] = None,
) -> Dict[int, float]:
    """Return {metal_idx: max-angle-deviation-deg} vs nearest ideal polyhedron.

    ``only_geom`` (default None -> unchanged best-of behaviour): when a polyhedron NAME is given and it
    is in ``_GEOM_IDEAL_ANGLES``, the metal's angles are scored vs THAT polyhedron's ideal targets only
    (no best-of).  Used by the FEAS_PREFERRED realism floor to score a recovered arrangement vs the
    polyhedron it was built for -- so a distorted octahedron cannot hide behind the trig-prism targets.

    Pure helper — does NOT reject anything, just measures how far each
    metal's L-M-L angles are from the closest textbook polyhedron for
    its coordination number.  Ideal angles by CN:

    - 2  → LIN (180)
    - 3  → TP / TS (120 / 90·2+180)
    - 4  → Td / SP  (109.47 / 90·4+180·2)
    - 5  → TBP / SPY (90·6+120·3+180 / 90·8+180·2)
    - 6  → Oh / TPR  (90·12+180·3 / 76·6+82·3+140·6)
    - 7  → PBP       (72·5+90·10+144·5+180)
    - 8  → SAP / DD  (52.4·4+73.1·4+118.5·4+143.1·4+180·2 / etc.)

    For each L-M-L pair we compute angle, measure distance to the nearest
    ideal target for the CN, and track the maximum for that metal.  The
    returned per-metal number in degrees is compared against
    ``DELFIN_IDEAL_POLYHEDRON_MAX_DEV``.  Used by ``_has_bad_geometry``
    (gate) and available for ranking / reporting.
    """
    out: Dict[int, float] = {}
    if not RDKIT_AVAILABLE:
        return out
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return out
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        metal_pos = conf.GetAtomPosition(metal_idx)
        neighbors = list(atom.GetNeighbors())
        cn = len(neighbors)
        if cn < 2:
            continue
        coord_positions = [conf.GetAtomPosition(nb.GetIdx()) for nb in neighbors]
        angles: List[float] = []
        for i in range(cn):
            for j in range(i + 1, cn):
                pa, pb = coord_positions[i], coord_positions[j]
                v1 = (pa.x - metal_pos.x, pa.y - metal_pos.y, pa.z - metal_pos.z)
                v2 = (pb.x - metal_pos.x, pb.y - metal_pos.y, pb.z - metal_pos.z)
                mag1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
                mag2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    continue
                dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angles.append(math.degrees(math.acos(cos_a)))
        if not angles:
            continue
        # INTENDED-GEOM scoring (only_geom): score vs the polyhedron the frame was BUILT for, no best-of.
        if only_geom and only_geom in _GEOM_IDEAL_ANGLES:
            ideals_list = [_GEOM_IDEAL_ANGLES[only_geom]]
        # Per-CN ideal angle targets, picked best-of over competing
        # polyhedra at the same CN.
        elif cn == 2:
            ideals_list = [[180.0]]
        elif cn == 3:
            ideals_list = [[120.0], [90.0, 180.0]]
        elif cn == 4:
            ideals_list = [[109.47], [90.0, 180.0]]
        elif cn == 5:
            ideals_list = [[90.0, 120.0, 180.0], [90.0, 180.0]]
        elif cn == 6:
            ideals_list = [[90.0, 180.0], [76.0, 82.0, 140.0]]
        elif cn == 7:
            ideals_list = [[72.0, 90.0, 144.0, 180.0]]
        elif cn == 8:
            ideals_list = [
                [52.4, 73.1, 118.5, 143.1, 180.0],
                [62.2, 73.7, 117.4, 143.6, 180.0],
            ]
        else:
            ideals_list = [[90.0, 180.0]]
        best_max_dev = float("inf")
        for ideals in ideals_list:
            max_dev = max(min(abs(a - t) for t in ideals) for a in angles)
            if max_dev < best_max_dev:
                best_max_dev = max_dev
        out[metal_idx] = best_max_dev
    return out


def _donor_h_points_at_metal(mol, conf_id: int, max_meh_deg: float = 80.0) -> bool:
    """FF-FREE geometric realism check -- MANTA-NATIVE (does NOT import WEDDELL: the construction filters
    with its OWN logic; the eye stays an independent judge).  Return True if ANY coordinating donor (N/O/S/
    Se/Te bearing an H) has an H pointing TOWARD the metal -- i.e. the smallest M-donor-H angle < max_meh_deg.

    A realistic coordinating N-H / O-H coordinates through its LONE PAIR (toward the metal), so its H points
    AWAY (M-D-H ~= 100-115 deg tetrahedral/pyramidal).  An H tilted toward the metal (small M-D-H) is
    physically unrealistic (the H crowds the metal, the lone pair points away).  Pure GEOMETRY (one angle) --
    no force field, no energy, DELFIN's own H is placed by construction so its position is meaningful."""
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            mp = conf.GetAtomPosition(atom.GetIdx())
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() not in ("N", "O", "S", "Se", "Te"):
                    continue
                dp = conf.GetAtomPosition(nb.GetIdx())
                for hb in nb.GetNeighbors():
                    if hb.GetSymbol() != "H":
                        continue
                    hp = conf.GetAtomPosition(hb.GetIdx())
                    v1 = (mp.x - dp.x, mp.y - dp.y, mp.z - dp.z)
                    v2 = (hp.x - dp.x, hp.y - dp.y, hp.z - dp.z)
                    n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
                    n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
                    if n1 < 1e-6 or n2 < 1e-6:
                        continue
                    cos_a = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)))
                    if math.degrees(math.acos(cos_a)) < max_meh_deg:
                        return True   # an H points toward the metal -> unrealistic donor orientation
        return False
    except Exception:
        return False


def _has_bad_geometry(mol, conf_id: int) -> bool:
    """Return True if the conformer has unrealistic metal-ligand geometry.

    Checks:
    1. Metal-ligand bond lengths must be within metal-specific distance range
       (derived from covalent radii, fallback to 1.4-3.5 Å)
    2. L-M-L angles: chelate bite angles (donors in the same chelate ring)
       may be as small as 40°; all other L-M-L pairs must be >=50°.
       This correctly allows 5- and 6-membered chelate ring bite angles
       (~55-70°) while still rejecting collapsed non-chelate geometries.
    3. Optional (env-gated): max single-angle deviation from the nearest
       ideal polyhedron must be ≤ ``DELFIN_IDEAL_POLYHEDRON_MAX_DEV``
       degrees.  Disabled when the env-var is 0.0 (default).
    """
    conf = mol.GetConformer(conf_id)
    if _has_unphysical_metal_nonbonded_contact(mol, conf_id):
        return True
    if _has_unphysical_oco_geometry(mol, conf_id):
        return True
    if _has_pi_ring_nonplanarity(mol, conf_id):
        return True
    if DELFIN_IDEAL_POLYHEDRON_MAX_DEV > 0.0:
        devs = _ideal_polyhedron_angle_dev_per_metal(mol, conf_id)
        if any(dev > DELFIN_IDEAL_POLYHEDRON_MAX_DEV for dev in devs.values()):
            return True
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        metal_pos = conf.GetAtomPosition(metal_idx)
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            continue

        metal_sym = atom.GetSymbol()
        nbr_indices = [nb.GetIdx() for nb in neighbors]
        coord_positions = []
        for nbr in neighbors:
            nbr_pos = conf.GetAtomPosition(nbr.GetIdx())
            dx = nbr_pos.x - metal_pos.x
            dy = nbr_pos.y - metal_pos.y
            dz = nbr_pos.z - metal_pos.z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            min_d, max_d = _ml_distance_range(metal_sym, nbr.GetSymbol())
            if dist < min_d or dist > max_d:
                return True
            coord_positions.append(nbr_pos)

        # Pre-compute chelate pairs (donors connected through non-metal path)
        chelate = _chelate_pairs(mol, metal_idx, nbr_indices)
        chelate_set = {frozenset(p) for p in chelate}
        n = len(coord_positions)
        n_pairs = n * (n - 1) // 2

        all_angles = []
        for i in range(n):
            for j in range(i + 1, n):
                pa, pb = coord_positions[i], coord_positions[j]
                v1 = (pa.x - metal_pos.x, pa.y - metal_pos.y, pa.z - metal_pos.z)
                v2 = (pb.x - metal_pos.x, pb.y - metal_pos.y, pb.z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    return True
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angle = math.degrees(math.acos(cos_a))
                all_angles.append((angle, frozenset([nbr_indices[i], nbr_indices[j]])))
                is_chelate_pair = frozenset([nbr_indices[i], nbr_indices[j]]) in chelate_set
                min_angle = 40.0 if is_chelate_pair else 50.0
                if angle < min_angle:
                    return True

        # NOTE: The old macrocyclic tetrahedral-rejection check (CN=4, all
        # pairs chelate, max angle < 115°) was removed.  It incorrectly
        # rejected valid tetrahedral structures (ideal 109.5°).  With the
        # BFS path-length cutoff in _chelate_pairs (max_path=4), macrocyclic
        # opposite-donor pairs are no longer marked as chelate, so the
        # "all pairs chelate" condition rarely triggers anyway.  The
        # _geometry_quality_score function already handles tetrahedral vs
        # square-planar ranking correctly.
    return False


def _geometry_quality_score(mol, conf_id: int) -> float:
    """Score how regular the metal coordination geometry is (lower = better).

    Measures:
    - Spread of M-L bond lengths (std-dev, ideally 0 for identical ligands)
    - Angular regularity:
      - 4-coordinate centers: best of tetrahedral (109.5 deg) OR
        square-planar (90/180 deg) targets
      - Other coordinations: deviation from nearest ideal (90/180 deg)
    """
    conf = mol.GetConformer(conf_id)
    total_penalty = 0.0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_pos = conf.GetAtomPosition(atom.GetIdx())
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            continue

        # Bond length uniformity
        dists = []
        coord_positions = []
        donor_symbols: List[str] = []
        for nbr in neighbors:
            nbr_pos = conf.GetAtomPosition(nbr.GetIdx())
            dx = nbr_pos.x - metal_pos.x
            dy = nbr_pos.y - metal_pos.y
            dz = nbr_pos.z - metal_pos.z
            dists.append(math.sqrt(dx*dx + dy*dy + dz*dz))
            coord_positions.append(nbr_pos)
            donor_symbols.append(nbr.GetSymbol())

        if dists:
            # Heterolept-aware bond-length penalty: group distances by
            # donor element (C, N, O, P, S, Cl, ...) and only penalise
            # dispersion WITHIN each group.  On homoleptic complexes
            # every donor shares one symbol so the sum collapses to
            # the original global std_d * 10 (unchanged behaviour).
            # On heteroleptic (e.g. Mn(CO)3(CO2R)(dppe) with C, O, P)
            # the correct ideals are 1.84 / 2.1 / 2.35 A — penalising
            # the cross-group spread would bias the collapse toward
            # structures that accidentally compressed Mn-P onto
            # Mn-C distance, which is chemically wrong.
            from collections import defaultdict as _dd
            _by_sym: Dict[str, List[float]] = _dd(list)
            for _d, _s in zip(dists, donor_symbols):
                _by_sym[_s].append(_d)
            for _grp in _by_sym.values():
                if len(_grp) < 2:
                    continue
                _mg = sum(_grp) / len(_grp)
                _sg = math.sqrt(sum((_x - _mg) ** 2 for _x in _grp) / len(_grp))
                total_penalty += _sg * 10  # weight bond-length spread

        # Angle regularity. For 4-coordinate centers, allow both tetrahedral
        # and square-planar patterns and keep the better one.
        angles = []
        for i in range(len(coord_positions)):
            for j in range(i + 1, len(coord_positions)):
                pa, pb = coord_positions[i], coord_positions[j]
                v1 = (pa.x - metal_pos.x, pa.y - metal_pos.y, pa.z - metal_pos.z)
                v2 = (pb.x - metal_pos.x, pb.y - metal_pos.y, pb.z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    continue
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angles.append(math.degrees(math.acos(cos_a)))

        if not angles:
            continue

        # Polyhedron angle penalties weighted x2 globally so ideal-
        # polyhedron adherence (Oh / TBP / Td / SAP / etc.) has
        # enough magnitude to dominate within-bucket ranking.
        # Multiplied by DELFIN_SYMMETRY_WEIGHT (default 1.0) so callers
        # can amplify the polyhedron-fidelity pressure without touching
        # code — raise to 2.0+ when quality trumps diversity.
        _POLY_W = 2.0 * DELFIN_SYMMETRY_WEIGHT
        if len(neighbors) == 2 and len(angles) == 1:
            # CN=2 (linear): ideal angle = 180°
            total_penalty += _POLY_W * abs(angles[0] - 180.0)
        elif len(neighbors) == 3 and len(angles) == 3:
            # CN=3: best of trigonal-planar (all 120°) vs T-shaped (90°,90°,180°)
            tp_pen = sum(abs(a - 120.0) for a in angles)
            ts_targets = sorted([90.0, 90.0, 180.0])
            ts_pen = sum(
                abs(a - t) for a, t in zip(sorted(angles), ts_targets)
            )
            total_penalty += _POLY_W * min(tp_pen, ts_pen)
        elif len(neighbors) == 4 and len(angles) == 6:
            total_penalty += _POLY_W * _preferred_cn4_geometry_score(
                atom.GetSymbol(),
                donor_symbols,
                metal_pos,
                coord_positions,
                angles,
            )
        elif len(neighbors) == 5:
            # CN=5: score against best of TBP or SP ideal angles.
            tbp_targets = sorted([90, 90, 90, 90, 90, 90, 120, 120, 120, 180])
            sp_targets = sorted([90, 90, 90, 90, 100, 100, 100, 100, 180, 180])
            sorted_a = sorted(angles)
            tbp_pen = sum(abs(a - t) for a, t in zip(sorted_a, tbp_targets))
            sp_pen = sum(abs(a - t) for a, t in zip(sorted_a, sp_targets))
            total_penalty += _POLY_W * min(tbp_pen, sp_pen)
        elif len(neighbors) == 6:
            # CN=6: score against best of Oh or TPR (trigonal prism).
            # Oh (Oh): 12x90 + 3x180 (15 pairs)
            # TPR (D3h): 6x76 (cap-basal) + 3x82 (basal-basal) + 6x140 (cap-cap)
            oh_targets = sorted([90] * 12 + [180] * 3)
            tpr_targets = sorted([76] * 6 + [82] * 3 + [140] * 6)
            sorted_a = sorted(angles)
            oh_pen = sum(abs(a - t) for a, t in zip(sorted_a, oh_targets))
            tpr_pen = sum(abs(a - t) for a, t in zip(sorted_a, tpr_targets))
            total_penalty += _POLY_W * min(oh_pen, tpr_pen)
        elif len(neighbors) == 7:
            # CN=7 (PBP, D5h): penalize each angle against nearest of 72°/90°/144°/180°
            ideal_7 = [72.0, 90.0, 144.0, 180.0]
            for a in angles:
                total_penalty += _POLY_W * min(abs(a - t) for t in ideal_7)
        elif len(neighbors) == 8:
            # CN=8: score against best of SAP (D4d) or DD (D2d).
            ideal_sap = [52.4, 73.1, 118.5, 143.1, 180.0]
            ideal_dd = [62.2, 73.7, 117.4, 143.6, 180.0]
            sap_pen = sum(min(abs(a - t) for t in ideal_sap) for a in angles)
            dd_pen = sum(min(abs(a - t) for t in ideal_dd) for a in angles)
            total_penalty += _POLY_W * min(sap_pen, dd_pen)
        else:
            # General: penalize distance from nearest ideal (90 or 180)
            for a in angles:
                dev_90 = abs(a - 90)
                dev_180 = abs(a - 180)
                total_penalty += min(dev_90, dev_180)

        # Symmetry bonus (additive penalty so smaller = better) —
        # rewards structures where donor-donor angles cluster tightly
        # around shared ideal values (90/109.5/120/180 deg).  For a
        # fully symmetric polyhedron every angle sits exactly on the
        # closest ideal, so per-bucket std-dev is near zero.  A
        # distorted structure with jagged angles has wide buckets
        # and a larger bonus value.  The contribution is intentionally
        # mild (weight 0.5) so it tips ties between otherwise-equal
        # candidates but does not override the primary polyhedron
        # penalty above.  Homoleptic symmetric complexes win over
        # equally-unfit distorted ones.
        if angles:
            _IDEAL_ANG = (72.0, 90.0, 109.5, 120.0, 144.0, 180.0)
            _buckets: Dict[float, List[float]] = {t: [] for t in _IDEAL_ANG}
            for a in angles:
                _nearest = min(_IDEAL_ANG, key=lambda t: abs(a - t))
                _buckets[_nearest].append(a)
            for _ideal, _vals in _buckets.items():
                if len(_vals) < 2:
                    continue
                _m = sum(_vals) / len(_vals)
                _sd = math.sqrt(sum((x - _m) ** 2 for x in _vals) / len(_vals))
                total_penalty += 5.0 * DELFIN_SYMMETRY_WEIGHT * _sd

        # Point-group / Cn-axis bonus: for each metal, test candidate
        # rotation axes (each principal axis of the donor point cloud
        # + the metal-centroid vector) at n = 6, 5, 4, 3, 2 and keep
        # the highest n where the donor positions are invariant under
        # rotation by 360/n within tolerance.  Higher n -> larger
        # bonus (negative penalty).  Captures Oh / Td / D3h / C4 / C3
        # symmetry of the local coordination sphere.  Weight 2.0 so
        # a detected C6 axis shaves ~12 points off the total penalty
        # (comparable to one 3-deg angle deviation improvement).
        try:
            if len(coord_positions) >= 2:
                import numpy as _np
                _pts = _np.array(
                    [(p.x - metal_pos.x, p.y - metal_pos.y, p.z - metal_pos.z)
                     for p in coord_positions],
                    dtype=float,
                )
                _max_order = 1
                _cands: List[_np.ndarray] = []
                _centroid = _pts.mean(axis=0)
                _cn = float(_np.linalg.norm(_centroid))
                if _cn > 1e-6:
                    _cands.append(_centroid / _cn)
                for _v in _pts:
                    _vn = float(_np.linalg.norm(_v))
                    if _vn > 1e-6:
                        _cands.append(_v / _vn)
                try:
                    _u, _s, _vh = _np.linalg.svd(_pts, full_matrices=False)
                    for _row in _vh:
                        _rn = float(_np.linalg.norm(_row))
                        if _rn > 1e-6:
                            _cands.append(_row / _rn)
                except Exception:
                    pass
                _tol_sq = 0.09
                for _axis in _cands:
                    for _n in (6, 5, 4, 3, 2):
                        _theta = 2.0 * math.pi / _n
                        _c = math.cos(_theta)
                        _sth = math.sin(_theta)
                        _K = _np.array([
                            [0.0, -_axis[2], _axis[1]],
                            [_axis[2], 0.0, -_axis[0]],
                            [-_axis[1], _axis[0], 0.0],
                        ])
                        _R = _np.eye(3) + _sth * _K + (1.0 - _c) * _K @ _K
                        _rot = _pts @ _R.T
                        _ok = True
                        for _rp in _rot:
                            _dists2 = ((_pts - _rp) ** 2).sum(axis=1)
                            if float(_dists2.min()) > _tol_sq:
                                _ok = False
                                break
                        if _ok and _n > _max_order:
                            _max_order = _n
                    if _max_order >= 6:
                        break
                if _max_order > 1:
                    total_penalty -= 5.0 * DELFIN_SYMMETRY_WEIGHT * (_max_order - 1)
        except Exception:
            pass

    # Inter-ligand soft close-contact penalty.  For every pair of
    # non-bonded non-metal heavy atoms, a graduated penalty kicks
    # in between r_cov_sum * 1.10 (Rule 10 hard threshold) and
    # r_cov_sum * 1.50 (vdW-like contact cutoff).  Rewards
    # structures where ligand fragments stay out of each other's
    # way even when technically allowed by Rule 10.  Catches cases
    # like NHC-methyl C 2.08 A from carbonyl-O — not a perceived
    # bond, but visually clashing.  Penalty per pair is scaled so
    # a single 2.0 A C-O contact adds ~2 pts, multiple clashes
    # compound quickly and move the structure down in ranking.
    try:
        if RDKIT_AVAILABLE:
            _CONTACT_SOFT_FRAC = 1.50
            _CONTACT_WEIGHT = 5.0
            _bonded_pairs: set = set()
            for _b in mol.GetBonds():
                _ii = _b.GetBeginAtom().GetIdx()
                _jj = _b.GetEndAtom().GetIdx()
                _bonded_pairs.add((min(_ii, _jj), max(_ii, _jj)))
            _heavy = [
                a for a in mol.GetAtoms()
                if a.GetAtomicNum() > 1 and a.GetSymbol() not in _METAL_SET
            ]
            _conf_local = mol.GetConformer(conf_id)
            for _i in range(len(_heavy)):
                _ai = _heavy[_i]
                _pi = _conf_local.GetAtomPosition(_ai.GetIdx())
                _ri = _COVALENT_RADII.get(_ai.GetSymbol())
                if _ri is None:
                    continue
                for _j in range(_i + 1, len(_heavy)):
                    _aj = _heavy[_j]
                    _pair = (
                        min(_ai.GetIdx(), _aj.GetIdx()),
                        max(_ai.GetIdx(), _aj.GetIdx()),
                    )
                    if _pair in _bonded_pairs:
                        continue
                    _rj = _COVALENT_RADII.get(_aj.GetSymbol())
                    if _rj is None:
                        continue
                    _pj = _conf_local.GetAtomPosition(_aj.GetIdx())
                    _dx = _pi.x - _pj.x
                    _dy = _pi.y - _pj.y
                    _dz = _pi.z - _pj.z
                    _d = math.sqrt(_dx * _dx + _dy * _dy + _dz * _dz)
                    _sum = _ri + _rj
                    _cut = _CONTACT_SOFT_FRAC * _sum
                    if _d < _cut:
                        total_penalty += _CONTACT_WEIGHT * (_cut - _d) / _sum
    except Exception:
        pass

    # Chelate-ring planarity bonus: aromatic rings that coordinate to a metal
    # should be flat.  Penalize RMSD of ring atoms from the best-fit plane.
    try:
        Chem.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
        metal_atom_indices = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET}
        if metal_atom_indices and ring_info and ring_info.NumRings() > 0:
            for ring in ring_info.AtomRings():
                # Check whether ≥2 ring atoms are direct neighbours of a metal
                ring_set = set(ring)
                n_coord_in_ring = sum(
                    1 for ridx in ring_set
                    if any(nbr.GetIdx() in metal_atom_indices
                           for nbr in mol.GetAtomWithIdx(ridx).GetNeighbors())
                )
                if n_coord_in_ring < 2:
                    continue  # Not a chelate ring

                # Collect ring-atom 3D positions
                positions = []
                for ridx in ring:
                    pos = conf.GetAtomPosition(ridx)
                    positions.append((pos.x, pos.y, pos.z))

                if len(positions) < 3:
                    continue

                # Fit a plane via centroid + SVD-like normal estimation using
                # cross products of consecutive edge vectors (rotation-invariant).
                cx = sum(p[0] for p in positions) / len(positions)
                cy = sum(p[1] for p in positions) / len(positions)
                cz = sum(p[2] for p in positions) / len(positions)
                vecs = [(p[0]-cx, p[1]-cy, p[2]-cz) for p in positions]

                # Accumulate a rough normal via cross products of consecutive vectors
                nx, ny, nz = 0.0, 0.0, 0.0
                n_v = len(vecs)
                for vi in range(n_v):
                    a_v = vecs[vi]
                    b_v = vecs[(vi + 1) % n_v]
                    nx += a_v[1]*b_v[2] - a_v[2]*b_v[1]
                    ny += a_v[2]*b_v[0] - a_v[0]*b_v[2]
                    nz += a_v[0]*b_v[1] - a_v[1]*b_v[0]
                n_mag = math.sqrt(nx*nx + ny*ny + nz*nz)
                if n_mag < 1e-8:
                    continue
                nx /= n_mag
                ny /= n_mag
                nz /= n_mag

                # RMSD of atoms from the plane
                deviations = [abs(v[0]*nx + v[1]*ny + v[2]*nz) for v in vecs]
                planarity_rmsd = math.sqrt(
                    sum(d*d for d in deviations) / len(deviations)
                )
                total_penalty += planarity_rmsd * 5.0
    except Exception:
        pass  # Ring-planarity scoring is optional; never block normal scoring

    return total_penalty


def _enforce_metal_topology(mol, conf_id: int = 0, min_nonbonded: float = 2.5):
    """Push non-bonded atoms away from metals to preserve SMILES topology.

    Operates as a final post-processing step: for every metal, any atom
    NOT bonded to it in the molecular graph but closer than *min_nonbonded*
    gets radially pushed outward.  Hapto ring atoms are moved as a rigid
    group to preserve ring geometry.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return
    try:
        import numpy as np
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, arr):
        conf.SetAtomPosition(i, Point3D(float(arr[0]), float(arr[1]),
                                         float(arr[2])))

    n = mol.GetNumAtoms()
    metal_indices = [i for i in range(n)
                     if mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET]
    if not metal_indices:
        return

    metal_bonded: Dict[int, set] = {}
    for mi in metal_indices:
        bonded = set()
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            bonded.add(nbr.GetIdx())
        metal_bonded[mi] = bonded

    hapto_groups = _find_hapto_groups(mol)
    hapto_of: Dict[int, int] = {}
    for gi, (_, grp) in enumerate(hapto_groups):
        for a in grp:
            hapto_of[a] = gi
    group_atoms: Dict[int, List[int]] = {}
    for gi, (_, grp) in enumerate(hapto_groups):
        group_atoms[gi] = list(grp)

    for _pass in range(10):
        any_moved = False
        for mi in metal_indices:
            mpos = _gp(mi)
            bonded = metal_bonded[mi]
            for ai in range(n):
                if ai == mi or ai in bonded:
                    continue
                if mol.GetAtomWithIdx(ai).GetSymbol() == 'H':
                    continue
                if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
                    continue
                apos = _gp(ai)
                d = float(np.linalg.norm(apos - mpos))
                if d >= min_nonbonded or d < 1e-8:
                    continue
                push_dir = (apos - mpos) / d
                push_dist = (min_nonbonded - d) * 1.1
                gi = hapto_of.get(ai)
                if gi is not None:
                    for ha in group_atoms[gi]:
                        _sp(ha, _gp(ha) + push_dir * push_dist)
                else:
                    _sp(ai, apos + push_dir * push_dist)
                any_moved = True
        if not any_moved:
            break

    # ---- P4 Layer 1: DELFIN_HAPTO_INTRALIG_FLOOR (BEGIN) -----------------
    # Universal non-degeneracy floor (default OFF → byte-identical).  The
    # multi-eta macrocycle placement bug assigns every linked eta-pair the
    # same [0,0,1] direction and the same combined_centroid, so >=2 eta
    # groups sharing one macrocyclic ring overlap EXACTLY (verified MEWCIA:
    # 4 pairs at d=0.000 A).  Catastrophic 0-A overlaps are invisible to the
    # metal-push loop above (it only repels atoms from METALS, not from each
    # other).  This pass separates ANY non-bonded heavy-atom pair closer
    # than ``min_intra`` (default 0.70 A) by pushing them apart symmetrically
    # — turning degenerate geometry into non-degenerate geometry for ALL
    # classes.  Deterministic: when d ~= 0 the separation axis is a stable
    # index-derived unit vector (no RNG, no clock).  Never emits non-finite
    # coordinates.
    if _delfin_env_int("DELFIN_HAPTO_INTRALIG_FLOOR", 0):
        min_intra = 0.70
        # Bonded-pair set (skip true chemical bonds — never separate them).
        bonded_pairs: set = set()
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            bonded_pairs.add((min(a1, a2), max(a1, a2)))
        heavy = [i for i in range(n)
                 if mol.GetAtomWithIdx(i).GetSymbol() != 'H'
                 and mol.GetAtomWithIdx(i).GetSymbol() not in _METAL_SET]

        def _det_axis(i, j):
            # Stable, deterministic unit vector from the atom-index pair
            # (used only when the two atoms are coincident).  No RNG / clock.
            h = (1.0 + (i * 131 + j * 17) % 97) / 98.0
            vec = np.array([
                np.sin(6.2831853 * h),
                np.cos(6.2831853 * h),
                np.sin(3.1415927 * h),
            ])
            nrm = float(np.linalg.norm(vec))
            return vec / nrm if nrm > 1e-9 else np.array([0.0, 0.0, 1.0])

        for _ipass in range(10):
            moved_any = False
            for _ix in range(len(heavy)):
                i = heavy[_ix]
                pi = _gp(i)
                for _jx in range(_ix + 1, len(heavy)):
                    j = heavy[_jx]
                    if (i, j) in bonded_pairs:
                        continue
                    pj = _gp(j)
                    delta = pj - pi
                    d = float(np.linalg.norm(delta))
                    if d >= min_intra:
                        continue
                    if d < 1e-6:
                        axis = _det_axis(i, j)
                    else:
                        axis = delta / d
                    push = (min_intra - d) / 2.0 + 1e-3
                    _sp(i, pi - axis * push)
                    _sp(j, pj + axis * push)
                    pi = _gp(i)
                    moved_any = True
            if not moved_any:
                break
    # ---- P4 Layer 1: DELFIN_HAPTO_INTRALIG_FLOOR (END) -------------------

    # ---- Wave-5 MULTIHAPTO_MM_ENFORCE (BEGIN) ----------------------------
    # The "push non-bonded apart" loop above protects topology by repelling
    # spurious near-contacts; it does NOT pull declared M-M sigma bonds back
    # to their ideal length when ETKDG/UFF stretched them.  Agent-3 Wave-5
    # forensics on the 22-SMILES multi-hapto class showed 17/22 SMILES carry
    # an Sn-Ir / Sn-Rh / Sn-Co main-group-TM sigma bond that ends up at
    # 3.4 A in the output -- above the topology-detector cutoff
    # (r_Sn + r_Ir + 0.45 = 3.25 A) so the M-M bond is dropped from the
    # parsed XYZ graph, causing topology_intact = False.
    #
    # MM_ENFORCE pulls each SMILES-declared (M, M) neighbour pair toward
    # the ideal distance from ``_METAL_METAL_BOND_LENGTHS`` (or covalent-
    # radii sum + 0.4 A fallback for unlisted hetero-pairs).  Symmetric
    # 50/50 move so neither metal drags its ligand frame disproportionally.
    # Universal: only fires on actual M-M edges (mono-hapto / sigma /
    # no_metal classes have none -- no-op).
    #
    # Phase-1.5 wire-in DISABLED (Welle-3 T7.2 revert, 2026-05-16): the
    # 50/50 metal shift translates each metal by ~half the M-M correction
    # but leaves the surrounding sigma-coord atoms (Cl/Me/Cp-ring) in
    # place, collapsing M-D distances to <1.0 A (verified D-COIRSN,
    # D-DIVPIJ, D-TENMIL: Sn-Cl 2.30->0.59 A).  On the 22-SMILES
    # multi_hapto set the wire-in produced 0/22 passes and -32% isomer
    # yield versus pre-wire-in.  Surgical fix: revert to opt-in only;
    # rigid-drag rewrite of ``_enforce_smiles_mm_distances`` is the
    # follow-up patch.
    #
    # Resolution precedence (post-revert):
    #   1. DELFIN_MULTIHAPTO_MM_ENFORCE_CLASSES set -> use that list.
    #   2. DELFIN_MULTIHAPTO_MM_ENFORCE=1 -> enabled on every class.
    #   3. env unset -> disabled everywhere (matches pre-Phase-1.5).
    if _class_conditional_flag(
        "DELFIN_MULTIHAPTO_MM_ENFORCE", mol, default=0,
    ):
        # Welle-5j Agent G (2026-05-17): ``DELFIN_5J_G_MM_RIGID_DRAG``
        # opt-in fixes the documented Sn-Cl drag bug (welle3 T7.2
        # forensics): when MM_ENFORCE shifts an Sn toward Ir to satisfy
        # the M-M target, the legacy metal-only shift leaves Sn's Cl
        # ligands at their original position, collapsing Sn-Cl to
        # 0.59 A.  ``rigid_drag=True`` translates each metal's
        # exclusive (non-bridging) sigma ligand frame with the metal
        # so Sn-Cl / Sn-Me bonds stay at their ETKDG / UFF length.
        # Default OFF (0) → bit-exact metal-only shift, matches the
        # legacy behaviour for any operator who has already set
        # ``DELFIN_MULTIHAPTO_MM_ENFORCE=1`` without the new flag.
        _mm_rigid = bool(_delfin_env_int("DELFIN_5J_G_MM_RIGID_DRAG", 0))
        try:
            _enforce_smiles_mm_distances(
                mol, conf, metal_indices, _gp, _sp, np,
                rigid_drag=_mm_rigid,
            )
        except Exception as _mm_exc:
            logger.debug("MULTIHAPTO_MM_ENFORCE skipped: %s", _mm_exc)


    # ---- Wave-5 MULTIHAPTO_MM_ENFORCE (END) ------------------------------


def _collect_metal_ligand_frame(mol, metal_idx: int,
                                metal_set_local) -> set:
    """Collect indices of atoms in a metal's exclusive sigma ligand frame.

    BFS outward from ``metal_idx`` through non-metal atoms.  An atom is
    in the frame iff its graph-shortest path back to ``metal_idx`` does
    NOT pass through any other metal atom (i.e. it "belongs" only to
    this metal).  Bridging atoms (shared between two metals) are
    excluded so the caller can choose to leave them in place or apply
    a 50/50 split independently.

    Stops at hapto-group boundary as a safety: the hapto-group rigid
    rotation is handled elsewhere (``_enforce_metal_topology`` rigid
    push uses ``group_atoms``); pulling hapto carbons here would
    deform the η-ring.  Hydrogens on frame atoms are included.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
    metal_idx : int
        Index of the metal whose ligand frame to collect.
    metal_set_local : set
        Snapshot of metal indices in the molecule (caller-provided so
        we don't recompute per metal).

    Returns
    -------
    set[int]
        Frame atom indices (excludes ``metal_idx`` itself, excludes
        any atom that has a different metal in its bonded neighbours).
    """
    frame: set = set()
    visited: set = set([metal_idx])
    # Direct sigma neighbours that are non-metal heavy atoms seed the BFS.
    seed_neighbours: list = []
    for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors():
        ni = nbr.GetIdx()
        if ni in metal_set_local:
            continue
        # Skip atoms that ALSO bond to a different metal (bridging /
        # shared ligand); we don't want to drag them with this metal.
        other_metal = False
        for nn in nbr.GetNeighbors():
            nni = nn.GetIdx()
            if nni != metal_idx and nni in metal_set_local:
                other_metal = True
                break
        if other_metal:
            continue
        seed_neighbours.append(ni)
        frame.add(ni)
        visited.add(ni)
    # BFS through bonded non-metal atoms, never crossing another metal
    # and never crossing back to ``metal_idx`` (loops are skipped by
    # ``visited`` membership).  Hydrogens come along for free at each
    # step (they have no further out-edges).
    queue = list(seed_neighbours)
    while queue:
        cur = queue.pop()
        for nn in mol.GetAtomWithIdx(cur).GetNeighbors():
            nni = nn.GetIdx()
            if nni in visited:
                continue
            if nni in metal_set_local:
                # Hit a different metal — stop branching here.
                continue
            # Reject atoms that bond to any OTHER metal (shared ligand
            # bridge); they stay in place / get handled by the other
            # metal's frame instead.
            other_metal = False
            for nnn in nn.GetNeighbors():
                nnnidx = nnn.GetIdx()
                if nnnidx == cur:
                    continue
                if nnnidx in metal_set_local and nnnidx != metal_idx:
                    other_metal = True
                    break
            if other_metal:
                continue
            visited.add(nni)
            frame.add(nni)
            queue.append(nni)
    return frame


def _enforce_smiles_mm_distances(
    mol,
    conf,
    metal_indices: List[int],
    _gp,
    _sp,
    np,
    tolerance: float = 0.20,
    rigid_drag: bool = False,
) -> int:
    """Pull SMILES-declared M-M sigma-bonded pairs toward ideal distance.

    For every graph-neighbour pair ``(m_i, m_j)`` where both atoms are in
    ``_METAL_SET``, look up the ideal bond length in
    ``_METAL_METAL_BOND_LENGTHS`` (keyed by ``frozenset({sym_i, sym_j})``);
    fall back to ``r_cov(sym_i) + r_cov(sym_j) + 0.4`` if the pair is not
    in the dict.  When the current Euclidean distance deviates from the
    ideal by more than ``tolerance`` (default 0.20 A), shift both metals
    symmetrically along the connecting axis so the new distance equals
    the ideal.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule whose conformer is being adjusted in place.
    conf : rdkit.Chem.Conformer
        Conformer whose atom positions are updated (mutated via ``_sp``).
    metal_indices : list[int]
        Indices of metal atoms in ``mol`` (matches ``_METAL_SET`` lookup).
    _gp, _sp : Callable[[int], np.ndarray] / Callable[[int, np.ndarray], None]
        Position get / set helpers closed over ``conf`` (provided by the
        caller, ``_enforce_metal_topology``).
    np : module
        ``numpy`` handle (passed in so this helper does not re-import).
    tolerance : float, optional
        Distance deviation (A) below which no adjustment is made.
    rigid_drag : bool, optional
        When True (Welle-5j Agent G ``DELFIN_5J_G_MM_RIGID_DRAG`` opt-in),
        each metal's exclusive sigma ligand frame translates with the
        metal so Sn-Cl / Sn-Me bonds do NOT break.  Per
        ``ITER-multihapto_wirein`` welle3 T7.2 the legacy metal-only
        shift collapses Sn-Cl from 2.30 A to 0.59 A on D-COIRSN /
        D-DIVPIJ / D-TENMIL.  Default False = legacy metal-only shift
        (bit-exact when caller does not pass ``rigid_drag=True``).

    Returns
    -------
    int
        Number of M-M pairs actually moved.
    """
    moved = 0
    seen: set = set()
    # Per-metal frame cache: only computed when rigid_drag is requested
    # AND we actually move that metal.  Empty dict acts as a sentinel.
    frame_cache: Dict[int, set] = {}
    metal_set_local: set = set(metal_indices) if rigid_drag else set()
    for mi in metal_indices:
        sym_i = mol.GetAtomWithIdx(mi).GetSymbol()
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nbr.GetIdx()
            if ni == mi:
                continue
            sym_j = nbr.GetSymbol()
            if sym_j not in _METAL_SET:
                continue
            pair = (min(mi, ni), max(mi, ni))
            if pair in seen:
                continue
            seen.add(pair)

            target = _METAL_METAL_BOND_LENGTHS.get(frozenset({sym_i, sym_j}))
            if target is None:
                r_i = _COVALENT_RADII.get(sym_i)
                r_j = _COVALENT_RADII.get(sym_j)
                if r_i is None or r_j is None:
                    continue
                target = r_i + r_j + 0.4

            pos_i = _gp(mi)
            pos_j = _gp(ni)
            vec = pos_j - pos_i
            cur = float(np.linalg.norm(vec))
            if cur < 1e-6:
                continue
            if abs(cur - target) <= tolerance:
                continue
            unit = vec / cur
            shift = (target - cur) * 0.5
            # Symmetric 50/50: m_i moves -shift along unit, m_j moves
            # +shift along unit, so the new separation equals target.
            delta_i = -unit * shift
            delta_j = unit * shift
            _sp(mi, pos_i + delta_i)
            _sp(ni, pos_j + delta_j)
            if rigid_drag:
                # Translate each metal's exclusive ligand frame by the
                # same delta so Sn-Cl / Sn-Me / Cp-Sn-Cl tetrahedra
                # follow rigidly.  Bridging atoms (those bonded to a
                # different metal) are deliberately excluded — they
                # would be double-counted between the two frames.
                if mi not in frame_cache:
                    frame_cache[mi] = _collect_metal_ligand_frame(
                        mol, mi, metal_set_local
                    )
                if ni not in frame_cache:
                    frame_cache[ni] = _collect_metal_ligand_frame(
                        mol, ni, metal_set_local
                    )
                for fi in frame_cache[mi]:
                    if fi == mi or fi == ni:
                        continue
                    _sp(fi, _gp(fi) + delta_i)
                for fj in frame_cache[ni]:
                    if fj == mi or fj == ni:
                        continue
                    _sp(fj, _gp(fj) + delta_j)
            moved += 1
    return moved


def _segment_distance_sq(p1, p2, p3, p4) -> float:
    """Squared minimum distance between 3D line segments P1-P2 and P3-P4.

    Each point is a tuple ``(x, y, z)``.  Uses the parametric approach:
    closest points on segment ``P1 + s*(P2-P1)`` and ``P3 + t*(P4-P3)``
    with ``s, t`` clamped to ``[0, 1]``.
    """
    d1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    d2 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])
    r = (p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2])

    a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]
    e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]
    f = d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2]

    EPS = 1e-12
    if a <= EPS and e <= EPS:
        return r[0] ** 2 + r[1] ** 2 + r[2] ** 2

    if a <= EPS:
        t = max(0.0, min(1.0, f / e))
        v = (r[0] - t * d2[0], r[1] - t * d2[1], r[2] - t * d2[2])
        return v[0] ** 2 + v[1] ** 2 + v[2] ** 2

    c = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2]
    if e <= EPS:
        s = max(0.0, min(1.0, -c / a))
        v = (r[0] + s * d1[0], r[1] + s * d1[1], r[2] + s * d1[2])
        return v[0] ** 2 + v[1] ** 2 + v[2] ** 2

    b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]
    denom = a * e - b * b

    if abs(denom) > EPS:
        s = max(0.0, min(1.0, (b * f - c * e) / denom))
    else:
        s = 0.0

    t = (b * s + f) / e
    if t < 0.0:
        t = 0.0
        s = max(0.0, min(1.0, -c / a))
    elif t > 1.0:
        t = 1.0
        s = max(0.0, min(1.0, (b - c) / a))

    v = (r[0] + s * d1[0] - t * d2[0],
         r[1] + s * d1[1] - t * d2[1],
         r[2] + s * d1[2] - t * d2[2])
    return v[0] ** 2 + v[1] ** 2 + v[2] ** 2


def _has_ligand_intertwining(mol, conf_id: int, threshold: float = 0.3) -> bool:
    """Return True if ligands are physically intertwined in the conformer.

    Decomposes the molecule into ligand fragments by removing metal centers,
    then checks if any heavy-atom bond from one fragment comes closer than
    *threshold* Å to any heavy-atom bond from another fragment (segment-to-
    segment distance).  This catches unphysical conformers where ligand arms
    pass through each other.
    """
    conf = mol.GetConformer(conf_id)
    metal_idxs = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET}
    if not metal_idxs:
        return False

    # Build adjacency graph for non-metal atoms (bonds to metals removed)
    non_metal = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in _METAL_SET}
    adj: Dict[int, set] = {i: set() for i in non_metal}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in metal_idxs or j in metal_idxs:
            continue
        adj[i].add(j)
        adj[j].add(i)

    # Connected components = ligand fragments
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

    if len(fragments) < 2:
        return False

    # Map atom index -> fragment index
    atom_frag: Dict[int, int] = {}
    for fi, frag in enumerate(fragments):
        for aidx in frag:
            atom_frag[aidx] = fi

    # Collect heavy-atom bonds per fragment (skip bonds involving H)
    frag_bonds: Dict[int, List[Tuple[int, int]]] = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in metal_idxs or j in metal_idxs:
            continue
        if mol.GetAtomWithIdx(i).GetAtomicNum() <= 1:
            continue
        if mol.GetAtomWithIdx(j).GetAtomicNum() <= 1:
            continue
        fi = atom_frag.get(i)
        if fi is None:
            continue
        frag_bonds.setdefault(fi, []).append((i, j))

    # Check inter-fragment bond-segment distances
    thresh_sq = threshold * threshold
    fkeys = sorted(frag_bonds)
    for a_idx in range(len(fkeys)):
        for b_idx in range(a_idx + 1, len(fkeys)):
            for i1, j1 in frag_bonds[fkeys[a_idx]]:
                p1 = conf.GetAtomPosition(i1)
                p2 = conf.GetAtomPosition(j1)
                pa, pb = (p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z)
                for i2, j2 in frag_bonds[fkeys[b_idx]]:
                    p3 = conf.GetAtomPosition(i2)
                    p4 = conf.GetAtomPosition(j2)
                    pc, pd = (p3.x, p3.y, p3.z), (p4.x, p4.y, p4.z)
                    if _segment_distance_sq(pa, pb, pc, pd) < thresh_sq:
                        return True
    return False


def _conformer_rmsd(mol, cid_a: int, cid_b: int) -> float:
    """Heavy-atom RMSD between two conformers after Kabsch alignment.

    Uses RDKit's ``AlignMol`` on a temporary copy so the original molecule
    coordinates are not modified.  Falls back to a sorted distance-matrix
    comparison (rotation-invariant) if alignment fails.
    """
    try:
        from rdkit.Chem import rdMolAlign
        heavy_map = [(i, i) for i, a in enumerate(mol.GetAtoms())
                     if a.GetAtomicNum() > 1]
        if not heavy_map:
            return float('inf')
        # AlignMol modifies probe coordinates → work on a copy
        probe = Chem.RWMol(mol)
        return rdMolAlign.AlignMol(probe, mol, prbCid=cid_a, refCid=cid_b,
                                   atomMap=heavy_map)
    except Exception:
        pass

    # Fallback: sorted distance-matrix comparison (rotation-invariant)
    try:
        conf_a = mol.GetConformer(cid_a)
        conf_b = mol.GetConformer(cid_b)
        heavy = [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1]
        if len(heavy) < 2:
            return float('inf')
        dists_a: List[float] = []
        dists_b: List[float] = []
        for i in range(len(heavy)):
            pa = conf_a.GetAtomPosition(heavy[i])
            pb = conf_b.GetAtomPosition(heavy[i])
            for j in range(i + 1, len(heavy)):
                qa = conf_a.GetAtomPosition(heavy[j])
                qb = conf_b.GetAtomPosition(heavy[j])
                dists_a.append(math.sqrt(
                    (pa.x - qa.x) ** 2 + (pa.y - qa.y) ** 2 + (pa.z - qa.z) ** 2))
                dists_b.append(math.sqrt(
                    (pb.x - qb.x) ** 2 + (pb.y - qb.y) ** 2 + (pb.z - qb.z) ** 2))
        dists_a.sort()
        dists_b.sort()
        return math.sqrt(
            sum((a - b) ** 2 for a, b in zip(dists_a, dists_b)) / len(dists_a))
    except Exception:
        return float('inf')


def _angle_class(pos_metal, pos_a, pos_b) -> str:
    """Classify the angle A-Metal-B as 'cis' (<135 deg) or 'trans'.

    Uses 135 deg as threshold (midpoint between ideal octahedral 90 and
    180 deg) for robust classification despite geometric noise from
    distance-geometry embedding.
    """
    v1 = (pos_a.x - pos_metal.x, pos_a.y - pos_metal.y, pos_a.z - pos_metal.z)
    v2 = (pos_b.x - pos_metal.x, pos_b.y - pos_metal.y, pos_b.z - pos_metal.z)
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    if mag1 < 1e-8 or mag2 < 1e-8:
        return 'cis'
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    angle_deg = math.degrees(math.acos(cos_angle))
    return 'cis' if angle_deg < 135 else 'trans'


def _estimate_isomer_upper_bound(
    mol, dtype_map: Optional[Dict[int, tuple]] = None
) -> Optional[int]:
    """Conservative upper bound on the number of distinct coordination isomers.

    Uses the standard Pólya counting simplification:  for a metal with
    coordination number CN and donor-type multiset :math:`k_1, k_2, \\dots`,
    the number of arrangements on a polyhedron of order :math:`|G|` is
    bounded above by

        CN! / (k_1! · k_2! · …) / max(1, CN // 2)

    The CN//2 divisor approximates the polyhedron rotation symmetry
    (Td: 12, D3h: 6, Oh: 24); this is intentionally generous so legitimate
    isomer counts pass.  For multi-metallic systems the per-metal bounds
    are multiplied and then clipped to an absolute ceiling.  When the
    bound cannot be computed (unusual CN, unclassifiable donors) returns
    ``None`` so callers skip enforcement.

    Intended use: log / soft-cap excess duplicates that the existing
    fingerprint dedup failed to collapse.  Not a hard reject — the
    symmetry-aware dedup and RMSD passes do the actual collapsing.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        if dtype_map is None:
            dtype_map = _donor_type_map(mol)
        import math as _m
        metals = [a for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        if not metals:
            return None
        per_metal: List[int] = []
        for m in metals:
            nbrs = [nb.GetIdx() for nb in m.GetNeighbors()]
            cn = len(nbrs)
            if cn < 2:
                per_metal.append(1)
                continue
            # Group donors by type
            type_counts: Dict[tuple, int] = {}
            for nb in nbrs:
                t = dtype_map.get(nb, (mol.GetAtomWithIdx(nb).GetSymbol(), frozenset()))
                type_counts[t] = type_counts.get(t, 0) + 1
            # CN! / product(k_i!)
            arrangements = _m.factorial(cn)
            for k in type_counts.values():
                arrangements //= _m.factorial(k)
            # Divide by approximate polyhedron order
            denom = max(1, cn // 2)
            bound = max(1, arrangements // denom)
            per_metal.append(bound)
        # Multi-metallic: multiply but cap
        product = 1
        for b in per_metal:
            product *= b
            if product > 500:
                product = 500
                break
        return int(product)
    except Exception:
        return None
