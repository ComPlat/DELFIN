"""Donor typing, viable and bridging donors, ligand fragments, coordination fingerprints, canonical polyhedron forms and isomer labels of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import re
from typing import Dict, FrozenSet, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_int,
)
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _METAL_SET,
)

logger = get_logger("delfin.smiles_converter")


def _donor_type_map(mol) -> Dict[int, tuple]:
    """Compute a donor type for each atom based on its chemical environment.

    Uses Morgan fingerprint bits at radius 2 to distinguish chemically
    inequivalent atoms of the same element (e.g. N atoms in a tridentate
    terpyridine vs. a bidentate phenylpyridine).

    The fingerprint is computed on a **charge-neutralised copy** of the
    template molecule so that SMILES variants which differ only in formal
    charge annotation (for example ``O#[C][Fe-2]([C]#O)(...)`` versus
    ``O#[C+][Fe-3]([C+]#O)(...)``) produce identical donor hashes.  The
    bonding pattern of every donor is preserved; only the annotation layer
    that encodes user bookkeeping (formal charge) is ignored.  Donor
    classes therefore track the chemical identity of the coordinating
    environment (CO vs carbene vs pyridine N) rather than the SMILES
    notation.

    Returns:
        Dict mapping atom index to ``(element_symbol, env_hash)`` where
        ``env_hash`` is a frozenset of Morgan bits unique to that
        chemical environment.  ``element_symbol`` is taken from the
        original (still charged) molecule, but the hash is computed on
        the neutralised copy.
    """
    result: Dict[int, tuple] = {}
    if not RDKIT_AVAILABLE:
        return result

    # Charge-neutralised working copy.  Atom ordering is preserved by
    # RDKit's RWMol copy constructor, so the atom index of every donor
    # in ``neutral_mol`` matches its index in the caller's ``mol``.
    try:
        neutral_mol = Chem.RWMol(mol)
        for _atom in neutral_mol.GetAtoms():
            try:
                _atom.SetFormalCharge(0)
            except Exception:
                pass
        try:
            neutral_mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    except Exception:
        neutral_mol = mol

    atom_bits: Dict[int, set] = {}
    try:
        Chem.FastFindRings(neutral_mol)
    except Exception:
        pass
    try:
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        gen.GetFingerprint(neutral_mol, additionalOutput=ao)
        for bit_id, entries in ao.GetBitInfoMap().items():
            for atom_idx, _radius in entries:
                if atom_idx not in atom_bits:
                    atom_bits[atom_idx] = set()
                atom_bits[atom_idx].add(bit_id)
    except Exception:
        # Fallback to legacy API
        try:
            from rdkit.Chem import rdMolDescriptors
            fp_info: Dict[int, list] = {}
            rdMolDescriptors.GetMorganFingerprint(neutral_mol, 2, bitInfo=fp_info)
            for bit_id, entries in fp_info.items():
                for atom_idx, _radius in entries:
                    if atom_idx not in atom_bits:
                        atom_bits[atom_idx] = set()
                    atom_bits[atom_idx].add(bit_id)
        except Exception:
            pass

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        env_hash = frozenset(atom_bits.get(idx, set()))
        result[idx] = (atom.GetSymbol(), env_hash)

    return result


def _compute_coordination_fingerprint(mol, conf_id: int, dtype_map: Optional[Dict[int, tuple]] = None) -> tuple:
    """Compute a hashable fingerprint describing the coordination geometry.

    Hybrid approach for robustness:
    1. **Trans pairs** (top floor(N/2) angles): which element pairs sit
       across the metal.  Robust against noise for mixed-element donors.
    2. **Same-element cis/trans pattern**: for each pair of identical-
       element donors, classify as cis (<135 deg) or trans.  This adds
       sensitivity for all-same-element complexes (e.g. Fe with 6 N)
       where trans-pair elements alone cannot distinguish isomers.
    3. **Detailed trans pairs**: uses Morgan-invariant-enriched donor types
       to distinguish same-element donors in different chemical environments
       (e.g. N in terpyridine vs. N in phenylpyridine of an MABC complex).

    Args:
        mol: RDKit molecule with conformers
        conf_id: Conformer ID to analyze
        dtype_map: Pre-computed donor type map from ``_donor_type_map(mol)``.
            If None, computed on the fly.
    """
    conf = mol.GetConformer(conf_id)
    fp_parts: List[tuple] = []

    # Use pre-computed donor types or compute on the fly
    if dtype_map is None:
        dtype_map = _donor_type_map(mol)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue

        metal_pos = conf.GetAtomPosition(atom.GetIdx())

        coord_atoms = sorted(
            ((nbr.GetSymbol(), nbr.GetIdx()) for nbr in atom.GetNeighbors()),
            key=lambda x: (x[0], x[1]),
        )
        n_coord = len(coord_atoms)
        if n_coord < 2:
            continue

        # Compute all pairwise angles
        # entries: (angle, symA, symB, idxA, idxB)
        angle_pairs: List[Tuple[float, str, str, int, int]] = []
        for i in range(n_coord):
            for j in range(i + 1, n_coord):
                sym_a, idx_a = coord_atoms[i]
                sym_b, idx_b = coord_atoms[j]
                v1 = (conf.GetAtomPosition(idx_a).x - metal_pos.x,
                      conf.GetAtomPosition(idx_a).y - metal_pos.y,
                      conf.GetAtomPosition(idx_a).z - metal_pos.z)
                v2 = (conf.GetAtomPosition(idx_b).x - metal_pos.x,
                      conf.GetAtomPosition(idx_b).y - metal_pos.y,
                      conf.GetAtomPosition(idx_b).z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    continue
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angle = math.degrees(math.acos(cos_a))
                angle_pairs.append((angle, sym_a, sym_b, idx_a, idx_b))

        # Part 1: trans pair elements via disjoint pairing that maximizes
        # total trans-angle sum (more robust than taking top N/2 raw angles).
        n_trans = n_coord // 2
        trans_pairs_raw: List[Tuple[str, str]] = []
        # Also track detailed trans pairs using donor types
        detailed_trans_raw: List[Tuple[tuple, tuple]] = []
        if n_coord % 2 == 0 and n_coord <= 8 and angle_pairs:
            idx_by_symbol = {idx: sym for sym, idx in coord_atoms}
            angle_map: Dict[Tuple[int, int], float] = {}
            for angle, _sa, _sb, ia, ib in angle_pairs:
                key = (ia, ib) if ia < ib else (ib, ia)
                angle_map[key] = angle

            donor_indices = tuple(sorted(idx_by_symbol.keys()))
            memo: Dict[Tuple[int, ...], Tuple[float, List[Tuple[int, int]]]] = {}

            def _best_pairing(rem: Tuple[int, ...]) -> Tuple[float, List[Tuple[int, int]]]:
                if len(rem) < 2:
                    return 0.0, []
                if rem in memo:
                    return memo[rem]

                first = rem[0]
                best_score = -1.0
                best_pairs: List[Tuple[int, int]] = []
                for k in range(1, len(rem)):
                    second = rem[k]
                    key = (first, second) if first < second else (second, first)
                    angle_val = angle_map.get(key, 0.0)
                    rest = rem[1:k] + rem[k + 1:]
                    sub_score, sub_pairs = _best_pairing(rest)
                    score = angle_val + sub_score
                    if score > best_score:
                        best_score = score
                        best_pairs = [(first, second)] + sub_pairs

                memo[rem] = (best_score, best_pairs)
                return memo[rem]

            _score, pair_idx = _best_pairing(donor_indices)
            for ia, ib in pair_idx[:n_trans]:
                trans_pairs_raw.append((idx_by_symbol[ia], idx_by_symbol[ib]))
                detailed_trans_raw.append((dtype_map.get(ia, (idx_by_symbol[ia],)),
                                           dtype_map.get(ib, (idx_by_symbol[ib],))))
        else:
            angle_pairs.sort(key=lambda x: -x[0])  # descending
            for _angle, sa, sb, ia, ib in angle_pairs[:n_trans]:
                trans_pairs_raw.append((sa, sb))
                detailed_trans_raw.append((dtype_map.get(ia, (sa,)),
                                           dtype_map.get(ib, (sb,))))

        trans_pairs = tuple(sorted(
            tuple(sorted((sa, sb))) for sa, sb in trans_pairs_raw
        ))
        # Normalize donor-type tuples for sortable comparison:
        # frozensets don't have total ordering, so convert to sorted tuples.
        def _norm_donor_type(dt: tuple) -> tuple:
            return tuple(
                tuple(sorted(x)) if isinstance(x, (frozenset, set)) else x
                for x in dt
            )
        detailed_trans = tuple(sorted(
            tuple(sorted((_norm_donor_type(da), _norm_donor_type(db))))
            for da, db in detailed_trans_raw
        ))

        # Part 2: pairwise cis/trans pattern for ALL same-element donor
        # pairs (not just when all donors share one element).  This adds
        # sensitivity for mixed complexes like MA2B2C2 where two different
        # "cis-A" arrangements produce identical trans-pair signatures but
        # differ in which same-element pairs are cis vs trans.
        same_elem_pattern: List[tuple] = []
        for angle, sa, sb, _ia, _ib in angle_pairs:
            if sa == sb:  # same element pair
                cls = 'cis' if angle < 135 else 'trans'
                same_elem_pattern.append((sa, cls))
        same_elem_pattern.sort()

        # Optional refined fingerprint (DELFIN_FP_REFINED=1, default 0):
        # additionally encode the multiset of CIS-pair Morgan-types so
        # arrangements that share trans-pair-types but differ in cis-
        # pair-types (e.g. DD A-trapezoid vs B-trapezoid type assignment,
        # or TBP axial-vs-equatorial donor placements) get distinct
        # fingerprints.  Defaults off because the existing dedup-pipeline
        # is iteration-order-sensitive and any fingerprint refinement
        # cascades into Cd-histidine-class regressions; safe to enable
        # opt-in for users who need finer stereochemistry distinction.
        if _delfin_env_int("DELFIN_FP_REFINED", 0):
            cis_detailed_raw: List[Tuple[tuple, tuple]] = []
            for angle, _sa, _sb, ia, ib in angle_pairs:
                if angle >= 135.0:
                    continue  # trans, already encoded above
                ta = _norm_donor_type(dtype_map.get(ia, (mol.GetAtomWithIdx(ia).GetSymbol(),)))
                tb = _norm_donor_type(dtype_map.get(ib, (mol.GetAtomWithIdx(ib).GetSymbol(),)))
                cis_detailed_raw.append(tuple(sorted((ta, tb))))
            cis_detailed = tuple(sorted(cis_detailed_raw))
            fp_parts.append((trans_pairs, tuple(same_elem_pattern), detailed_trans, cis_detailed))
        else:
            fp_parts.append((trans_pairs, tuple(same_elem_pattern), detailed_trans))

    return tuple(sorted(fp_parts))


def _classify_isomer_label(fingerprint: tuple, mol) -> str:
    """Translate a coordination fingerprint to a human-readable label.

    The fingerprint is ``((trans_pairs, same_elem_pattern, detailed_trans), ...)``
    per metal center.  Uses trans-pair analysis to classify all common
    octahedral and 4-coordinate patterns:

    6-coordinate (octahedral):
    - **MA6** / **MA5B**: single isomer → ``''``
    - **MA4B2**: cis / trans (minority pair B trans or not)
    - **MA3B3**: fac / mer
    - **MA2B2C2**: all-cis / all-trans / X-trans (X = element
      whose pair is trans while the other two are cis)
    - **MA3B2C**: fac/mer for A + cis/trans for B

    4-coordinate:
    - **MA2B2**: cis / trans
    """
    # Count coordinating atoms per element from the mol
    n_coord = 0
    elem_counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        for nbr in atom.GetNeighbors():
            sym = nbr.GetSymbol()
            elem_counts[sym] = elem_counts.get(sym, 0) + 1
            n_coord += 1

    count_signature = sorted(elem_counts.values())

    for metal_fp in fingerprint:
        trans_pairs, same_elem_pattern = metal_fp[0], metal_fp[1]

        # Count same-element trans pairs from trans_pairs
        same_trans: Dict[str, int] = {}
        for pair in trans_pairs:
            if pair[0] == pair[1]:
                same_trans[pair[0]] = same_trans.get(pair[0], 0) + 1

        # Supplement from same_elem_pattern ONLY when trans_pairs contain
        # exclusively same-element pairs (= all donors share one element).
        # In this case trans_pairs cannot distinguish isomers by element
        # and the cis/trans angle classification is the only differentiator.
        # For mixed-element complexes, trans_pairs already carry the needed
        # information; adding same_elem_pattern would double-count.
        has_hetero_trans = any(p[0] != p[1] for p in trans_pairs)
        if not has_hetero_trans:
            for sym, cls in same_elem_pattern:
                if cls == 'trans':
                    if sym not in same_trans:
                        same_trans[sym] = 0
                    same_trans[sym] += 1

        # --- 6-coordinate patterns ---
        if n_coord == 6:
            # MA6 or MA5B with homogeneous element set: usually a single isomer.
            # Exception: when all donors share one element symbol but have two
            # distinct Morgan-hash types (3+3), we can still classify fac/mer
            # using the detailed_trans information.  This handles tris-bidentate
            # complexes like Fe(citrate)3 where e.g. all donors are O but split
            # into alkoxo-O and carboxylato-O.
            if len(elem_counts) == 1 or count_signature == [1, 5]:
                if count_signature != [1, 5] and len(metal_fp) > 2:
                    detailed = metal_fp[2]  # tuple of (type_i, type_j) trans pairs
                    # Collect all donor types present
                    all_types = set()
                    for pair in detailed:
                        all_types.add(pair[0])
                        all_types.add(pair[1])
                    if len(all_types) == 2:
                        # Two distinct Morgan types, 3+3 split: fac or mer
                        n_same_type_trans = sum(1 for p in detailed if p[0] == p[1])
                        if n_same_type_trans == 0:
                            return 'fac'
                        return 'mer'
                return ''

            # MA4B2: cis/trans based on minority element (count==2).
            # When the majority 4-atom set splits into two Morgan-hash
            # subclasses (e.g. quinoline-N + imine-N for a hexadentate
            # bis-Schiff base), further distinguish the otherwise-generic
            # "cis" bucket by whether the majority subclasses stand trans
            # or cis to each other.  Without this, the X-twist / mer-mer
            # crossed topology collapses onto the same "cis" label as
            # truly all-cis arrangements and gets fingerprint-dedupped
            # before the caller ever sees it.
            if count_signature == [2, 4]:
                minority = [s for s, c in elem_counts.items() if c == 2][0]
                majority = [s for s, c in elem_counts.items() if c == 4][0]
                if same_trans.get(minority, 0) >= 1:
                    return 'trans'
                if len(metal_fp) > 2:
                    detailed = metal_fp[2]
                    same_sub = 0
                    cross_sub = 0
                    for pair in detailed:
                        if (
                            isinstance(pair[0], tuple) and isinstance(pair[1], tuple)
                            and pair[0][0] == majority and pair[1][0] == majority
                        ):
                            if pair[0][1] == pair[1][1]:
                                same_sub += 1
                            else:
                                cross_sub += 1
                    if cross_sub >= 1 and same_sub == 0:
                        return f'cis-cross-{majority}'
                    if same_sub >= 1 and cross_sub == 0:
                        return f'cis-same-{majority}'
                    if same_sub >= 1 and cross_sub >= 1:
                        return f'cis-mixed-{majority}'
                return 'cis'

            # MA3B3: fac/mer
            if count_signature == [3, 3]:
                for sym, count in elem_counts.items():
                    if count == 3:
                        if same_trans.get(sym, 0) == 0:
                            return 'fac'
                        return 'mer'

            # MA2B2C2: all-cis / all-trans / X-trans
            if count_signature == [2, 2, 2]:
                elems_with_2 = [s for s, c in elem_counts.items() if c == 2]
                trans_elems = [s for s in elems_with_2 if same_trans.get(s, 0) >= 1]
                n_trans = len(trans_elems)
                if n_trans == 0:
                    return 'all-cis'
                if n_trans == 3:
                    return 'all-trans'
                if n_trans == 1:
                    return f'{trans_elems[0]}-trans'
                # 2 trans pairs: label by the one that is cis
                cis_elems = [s for s in elems_with_2 if s not in trans_elems]
                if len(cis_elems) == 1:
                    return f'{cis_elems[0]}-cis'
                return ''

            # MA3B2C: fac/mer for A (count==3) + cis/trans for B (count==2)
            if count_signature == [1, 2, 3]:
                a_sym = [s for s, c in elem_counts.items() if c == 3][0]
                b_sym = [s for s, c in elem_counts.items() if c == 2][0]
                a_label = 'mer' if same_trans.get(a_sym, 0) >= 1 else 'fac'
                b_label = 'trans' if same_trans.get(b_sym, 0) >= 1 else 'cis'
                return f'{a_label}-{b_label}'

            # MA4BC: the majority element has 4, the other two have 1 each.
            # Classify by which pairs are trans.
            if count_signature == [1, 1, 4]:
                majority = [s for s, c in elem_counts.items() if c == 4][0]
                minorities = sorted(s for s, c in elem_counts.items() if c == 1)
                # Check if the two minorities are trans to each other
                minority_trans = any(
                    (sorted(p) == sorted(minorities)) for p in trans_pairs
                )
                n_maj_trans = same_trans.get(majority, 0)
                if minority_trans:
                    return f'{"/".join(minorities)}-trans'
                if n_maj_trans >= 2:
                    return f'{majority}-all-trans'
                return f'{"/".join(minorities)}-cis'

        # --- 2-coordinate patterns ---
        if n_coord == 2:
            return 'linear'

        # --- 3-coordinate patterns ---
        if n_coord == 3:
            if len(elem_counts) == 1:
                return ''  # MA3: single isomer (trigonal-planar)
            # MA2B: cis/trans only for T-shaped geometry (where trans pair exists)
            if count_signature == [1, 2] and trans_pairs:
                minority = [s for s, c in elem_counts.items() if c == 1][0]
                is_trans = any(minority in p for p in trans_pairs)
                if is_trans:
                    return 'trans'
                return 'cis'
            return ''

        # --- 4-coordinate patterns ---
        if n_coord == 4:
            if count_signature == [2, 2]:
                minority = sorted(elem_counts.keys())[0]
                if same_trans.get(minority, 0) >= 1:
                    return 'trans'
                return 'cis'

        # --- 5-coordinate patterns ---
        if n_coord == 5:
            if len(elem_counts) == 1:
                return ''  # MA5: single isomer
            if count_signature == [1, 4]:
                minority = [s for s, c in elem_counts.items() if c == 1][0]
                # Axial donors appear in a trans-pair; equatorial don't
                is_axial = any(minority in p for p in trans_pairs)
                return 'axial' if is_axial else 'equatorial'
            if count_signature == [2, 3]:
                minority = [s for s, c in elem_counts.items() if c == 2][0]
                n_trans = same_trans.get(minority, 0)
                return 'diaxial' if n_trans >= 1 else 'ax-eq'

        # --- 7-coordinate patterns ---
        if n_coord == 7:
            if len(elem_counts) == 1:
                return ''
            # For PBP: axial pair = first trans pair (only 1 real 180° pair)
            axial_pair = tuple(sorted(trans_pairs[0])) if trans_pairs else ()
            if axial_pair:
                if axial_pair[0] == axial_pair[1]:
                    return f'{axial_pair[0]}-{axial_pair[0]}-ax'
                return f'{axial_pair[0]}-{axial_pair[1]}-ax'
            return 'all-eq'

        # --- 8-coordinate patterns ---
        if n_coord == 8:
            if len(elem_counts) == 1:
                return ''  # MA8: single isomer
            # For mixed-element CN=8: count trans-pair types for differentiation
            n_total_same_trans = sum(same_trans.values())
            if count_signature == [4, 4]:
                # MA4B4: classify by number of same-element trans pairs
                if n_total_same_trans == 0:
                    return 'all-cis'
                if n_total_same_trans >= 4:
                    return 'all-trans'
                return f'{n_total_same_trans}-trans'
            return ''

        # --- 9-coordinate patterns ---
        if n_coord == 9:
            if len(elem_counts) == 1:
                return ''
            return ''

    return ''


# ---------------------------------------------------------------------------
# Topological isomer enumerator (Feature 1)
# ---------------------------------------------------------------------------

def _chelate_pairs(mol, metal_idx: int, donor_indices: List[int],
                    max_path: int = 5) -> List[FrozenSet]:
    """Find donor pairs connected through a short non-metal path (chelate constraints).

    BFS from each donor atom to every other donor, blocking the metal.
    A successful path with length <= *max_path* bonds means the two donors
    form a chelate ring small enough to enforce a *cis* constraint (up to a
    ``max_path + 1``-membered chelate ring counting the metal).

    Longer paths (e.g. opposite donors of a 14-membered macrocycle) are
    **not** marked as chelate because they can adopt trans arrangements.

    Args:
        max_path: Maximum number of bonds in the non-metal path between two
            donors to count as a chelate pair.  Default 5 corresponds to a
            7-membered chelate ring (donor–5 bridge atoms–donor + metal),
            capturing common diphosphine (dppp) and diamine chelates.

    Returns a list of frozensets {donor_i_idx, donor_j_idx}.
    """
    pairs: List[FrozenSet] = []
    n = len(donor_indices)
    for i in range(n):
        for j in range(i + 1, n):
            start = donor_indices[i]
            target = donor_indices[j]
            # BFS with distance tracking, blocking the metal atom
            visited = {metal_idx, start}
            queue: List[Tuple[int, int]] = [(start, 0)]  # (atom_idx, distance)
            found_dist = -1
            while queue and found_dist < 0:
                current, dist = queue.pop(0)
                if dist >= max_path:
                    # Cannot reach target within max_path from here
                    continue
                for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni == target:
                        found_dist = dist + 1
                        break
                    if ni not in visited:
                        visited.add(ni)
                        queue.append((ni, dist + 1))
            if 0 < found_dist <= max_path:
                pairs.append(frozenset([donor_indices[i], donor_indices[j]]))
    return pairs


def _chelate_backbone_max_reach(mol, a: int, b: int) -> Optional[float]:
    """First-principles UPPER BOUND on a chelate's donor-donor bite (Å): the sum of covalent bond
    lengths along the shortest LIGAND-backbone path a->b with ALL metal atoms removed.

    By the triangle inequality the straight-line donor-donor distance can NEVER exceed this contour
    length, for ANY conformer -- so this is a SAMPLING-INDEPENDENT, universal feasibility bound derived
    from FIRST PRINCIPLES (the molecular graph + element covalent radii), NOT from a template conformer
    that ETKDG happened to embed.  A short backbone physically cannot span a large bite (a 2-atom bridge
    cannot reach a trans-octahedral separation -> correctly infeasible); a long/flexible backbone can
    reach any bite up to its contour (Ir bis-tridentate CAN reach octahedral -> correctly feasible).

    Returns the contour length, or None if a and b are not connected once the metals are removed (not a
    genuine through-backbone chelate).  Never raises."""
    try:
        from collections import deque as _deque
        _metals = {at.GetIdx() for at in mol.GetAtoms() if at.GetSymbol() in _METAL_SET}
        if a in _metals or b in _metals:
            return None
        _prev = {a: -1}
        _dq = _deque([a])
        _found = False
        while _dq:
            _u = _dq.popleft()
            if _u == b:
                _found = True
                break
            for _nb in mol.GetAtomWithIdx(_u).GetNeighbors():
                _v = _nb.GetIdx()
                if _v in _metals or _v in _prev:
                    continue
                _prev[_v] = _u
                _dq.append(_v)
        if not _found:
            return None
        _pt = Chem.GetPeriodicTable()
        _total = 0.0
        _cur = b
        while _prev[_cur] != -1:
            _p = _prev[_cur]
            _total += (_pt.GetRcovalent(mol.GetAtomWithIdx(_cur).GetAtomicNum())
                       + _pt.GetRcovalent(mol.GetAtomWithIdx(_p).GetAtomicNum()))
            _cur = _p
        return _total
    except Exception:
        return None


def _canonical_oh(types: tuple) -> tuple:
    """Canonical form for octahedral (3 trans pairs: 0-1, 2-3, 4-5)."""
    pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
        tuple(sorted([types[4], types[5]])),
    ]))
    return ('OH', pairs)


def _canonical_sq(types: tuple) -> tuple:
    """Canonical form for square-planar (2 trans pairs: 0-1, 2-3)."""
    pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
    ]))
    return ('SQ', pairs)


def _canonical_tbp(types: tuple) -> tuple:
    """Canonical form for trigonal-bipyramidal (axial: 0,1; equatorial: 2,3,4)."""
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3], types[4]]))
    return ('TBP', axial, equatorial)


def _canonical_sp(types: tuple) -> tuple:
    """Canonical form for square-pyramidal (apical: 0; basal: 1,2,3,4; trans pairs 1-3, 2-4)."""
    basal_pairs = tuple(sorted([
        tuple(sorted([types[1], types[3]])),
        tuple(sorted([types[2], types[4]])),
    ]))
    return ('SP', types[0], basal_pairs)


def _canonical_pbp(types: tuple) -> tuple:
    """Canonical form for pentagonal-bipyramidal (axial: 0,1; equatorial: 2-6)."""
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3], types[4], types[5], types[6]]))
    return ('PBP', axial, equatorial)


def _canonical_th(types: tuple) -> tuple:
    """Canonical form for tetrahedral (all 4 sites equivalent, no trans pairs)."""
    return ('TH', tuple(sorted(types)))


def _canonical_lin(types: tuple) -> tuple:
    """Canonical form for linear (2 sites, 1 trans pair: 0-1)."""
    return ('LIN', tuple(sorted(types)))


def _canonical_tp(types: tuple) -> tuple:
    """Canonical form for trigonal-planar (3 equivalent sites, no trans pairs)."""
    return ('TP', tuple(sorted(types)))


def _canonical_ts(types: tuple) -> tuple:
    """Canonical form for T-shaped (trans pair: 0-1; unique: 2)."""
    trans = tuple(sorted([types[0], types[1]]))
    return ('TS', trans, types[2])


def _canonical_sap(types: tuple) -> tuple:
    """Canonical form for square-antiprismatic (D4d).

    Positions 0-3 = top square, 4-7 = bottom square (staggered by 45°).
    Trans pairs: 0-6, 1-7, 2-4, 3-5 (opposite vertices through center).
    """
    trans_pairs = tuple(sorted([
        tuple(sorted([types[0], types[6]])),
        tuple(sorted([types[1], types[7]])),
        tuple(sorted([types[2], types[4]])),
        tuple(sorted([types[3], types[5]])),
    ]))
    return ('SAP', trans_pairs)


def _canonical_dd(types: tuple) -> tuple:
    """Canonical form for dodecahedral (D2d, = two interpenetrating trapezoids).

    Positions 0-3 = A sites (larger trapezoid), 4-7 = B sites (smaller).
    Trans pairs: 0-2, 1-3 (A-site), 4-6, 5-7 (B-site).
    """
    a_pairs = tuple(sorted([
        tuple(sorted([types[0], types[2]])),
        tuple(sorted([types[1], types[3]])),
    ]))
    b_pairs = tuple(sorted([
        tuple(sorted([types[4], types[6]])),
        tuple(sorted([types[5], types[7]])),
    ]))
    return ('DD', a_pairs, b_pairs)


def _canonical_ss(types: tuple) -> tuple:
    """Canonical form for see-saw / C2v (CN=4).

    Positions 0,1 = axial (trans pair), 2,3 = equatorial (cis pair).
    """
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3]]))
    return ('SS', axial, equatorial)


def _canonical_tpr(types: tuple) -> tuple:
    """Canonical form for trigonal-prismatic (D3h, CN=6).

    Positions 0-2 = top triangle, 3-5 = bottom triangle (eclipsed).
    No 180° trans pairs exist in a trigonal prism.
    """
    top = tuple(sorted(types[:3]))
    bottom = tuple(sorted(types[3:6]))
    combined = tuple(sorted([top, bottom]))
    return ('TPR', combined)


def _canonical_coh(types: tuple) -> tuple:
    """Canonical form for capped octahedral (C3v, CN=7).

    Positions 0-5 = octahedral base (3 trans pairs: 0-1, 2-3, 4-5),
    position 6 = capping atom above one triangular face.
    """
    base_pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
        tuple(sorted([types[4], types[5]])),
    ]))
    return ('COH', base_pairs, types[6])


def _canonical_ttp(types: tuple) -> tuple:
    """Canonical form for tricapped trigonal-prismatic (D3h, CN=9).

    Positions 0-5 = prism vertices (top 0-2, bottom 3-5 eclipsed),
    positions 6-8 = equatorial caps.
    """
    prism = tuple(sorted(types[:6]))
    caps = tuple(sorted(types[6:9]))
    return ('TTP', prism, caps)


# ---------------------------------------------------------------------------
# Iter-2 Subagent D — CN10/11/12 polyhedra (lanthanide / actinide complexes).
# Env-gated: DELFIN_CN_HIGH_ENABLE=1 enables CN ≥ 10 enumeration.  Default
# OFF (=0) → bit-exact HEAD (existing dispatcher returns [] for cn ≥ 10).
#
# Polyhedra:
#   CN10  BCSAP (D4d, bicapped square-antiprism)
#         PAP   (D5d, pentagonal antiprism)
#   CN11  CPAP  (C5v, capped pentagonal antiprism)
#   CN12  ICOS  (Ih, icosahedron)
#         CUBO  (Oh, cuboctahedron)
#         HBP   (D6h, hexagonal prism — eclipsed)
#
# Pólya/Burnside counts validated for binary multisets [P^k X^(n-k)]
# against group-theoretic tables; see staged_patches/iter2D_cn10_11_12_design.md
# ---------------------------------------------------------------------------


def _canonical_bcsap(types: tuple) -> tuple:
    """Canonical form for bicapped square-antiprism (D4d, CN=10).

    Positions 0-3 = top square, 4-7 = bottom square (rotated 45°),
    8-9 = axial caps.  Trans pairs: (0,2),(1,3),(4,6),(5,7),(8,9).
    """
    sap_pairs = tuple(sorted([
        tuple(sorted([types[0], types[2]])),
        tuple(sorted([types[1], types[3]])),
        tuple(sorted([types[4], types[6]])),
        tuple(sorted([types[5], types[7]])),
    ]))
    caps = tuple(sorted([types[8], types[9]]))
    return ('BCSAP', sap_pairs, caps)


def _canonical_pap(types: tuple) -> tuple:
    """Canonical form for pentagonal antiprism (D5d, CN=10).

    Positions 0-4 = top pentagon, 5-9 = bottom pentagon (staggered 36°).
    Two rings interchangeable under σ_d.
    """
    top = tuple(sorted(types[0:5]))
    bot = tuple(sorted(types[5:10]))
    rings = tuple(sorted([top, bot]))
    return ('PAP', rings)


def _canonical_cpap(types: tuple) -> tuple:
    """Canonical form for capped pentagonal antiprism (C5v, CN=11).

    Positions 0-4 = top pentagon (cap-side), 5-9 = bottom pentagon,
    10 = axial cap.  Cap breaks σ_h → top and bottom NOT interchangeable.
    """
    top = tuple(sorted(types[0:5]))
    bot = tuple(sorted(types[5:10]))
    return ('CPAP', top, bot, types[10])


def _canonical_icos(types: tuple) -> tuple:
    """Canonical form for icosahedron (Ih, CN=12).

    All 12 vertices equivalent under Ih (vertex-transitive group).
    """
    return ('ICOS', tuple(sorted(types)))


def _canonical_cubo(types: tuple) -> tuple:
    """Canonical form for cuboctahedron (Oh, CN=12).

    Positions 0-3 = xy-plane square, 4-7 = xz-plane square,
    8-11 = yz-plane square.  Three squares interchangeable under Oh.
    """
    sq_xy = tuple(sorted(types[0:4]))
    sq_xz = tuple(sorted(types[4:8]))
    sq_yz = tuple(sorted(types[8:12]))
    sqs = tuple(sorted([sq_xy, sq_xz, sq_yz]))
    return ('CUBO', sqs)


def _canonical_hbp(types: tuple) -> tuple:
    """Canonical form for hexagonal prism (D6h, CN=12).

    Positions 0-5 = top hexagon, 6-11 = bottom hexagon (eclipsed).
    Two rings interchangeable under σ_h.
    """
    top = tuple(sorted(types[0:6]))
    bot = tuple(sorted(types[6:12]))
    rings = tuple(sorted([top, bot]))
    return ('HBP', rings)


def _label_from_canonical_form(cf: tuple) -> str:
    """Derive a display label directly from an enumerator canonical form.

    Re-classifying post-UFF geometry via :func:`_classify_isomer_label`
    loses information when UFF drifts an axial donor off the 180° line:
    a deliberate N-N-axial PBP arrangement can be mis-read as N-O-ax,
    collapsing three distinct constitutional isomers (N-N-ax, N-O-ax,
    O-O-ax) into two.  Labelling from the enumerator's canonical form
    preserves the intended topology regardless of UFF drift and lets
    the energy-based sort decide which survives.

    Returns an empty string for geometries where canonical-form labelling
    is not informative (e.g. homoleptic TH/TP), so the caller can fall
    back to classify-based labels.
    """
    if not cf:
        return ''
    geom = cf[0]

    # First, scan the canonical form to find elements that appear with more
    # than one Morgan-hash class (e.g. carbonyl-C and carbene-C are both
    # "C" but with class indices "0" and "1").  When an element has multiple
    # classes, the display label must keep the class digit; otherwise
    # distinct constitutional isomers collapse into one label.
    def _collect_labels(item):
        out: List[str] = []
        if isinstance(item, (tuple, list)):
            for x in item:
                out.extend(_collect_labels(x))
        elif isinstance(item, str) and re.match(r'[A-Za-z]', item or ''):
            out.append(item)
        return out

    _by_elem: Dict[str, set] = {}
    for _lbl in _collect_labels(cf[1:]):
        _m = re.match(r'([A-Za-z]+)(\d*)', _lbl)
        if _m:
            _by_elem.setdefault(_m.group(1), set()).add(_m.group(2) or '0')

    def _strip(symbol: str) -> str:
        # Donor labels may be Morgan-hash enriched (e.g. 'N0', 'O1').  Drop
        # the class digit only when its element has a single class in this
        # canonical form; otherwise keep it so that e.g. carbonyl-C ('C0')
        # and carbene-C ('C1') remain distinguishable in the label.
        if not symbol:
            return symbol
        m = re.match(r'([A-Za-z]+)(\d*)', str(symbol))
        if not m:
            return str(symbol)
        elem, digit = m.group(1), m.group(2) or '0'
        if elem in _by_elem and len(_by_elem[elem]) > 1:
            return f'{elem}{digit}'
        return elem

    def _pair(p: tuple) -> str:
        a, b = _strip(p[0]), _strip(p[1])
        return f'{a}-{b}' if a <= b else f'{b}-{a}'

    if geom in ('PBP', 'TBP', 'SS'):
        ax = cf[1]
        return f'{_pair(ax)}-ax'
    if geom == 'LIN':
        t = cf[1]
        a, b = _strip(t[0]), _strip(t[1])
        return '' if a == b else f'{a}-{b}'
    if geom == 'TS':
        ax, uniq = cf[1], cf[2]
        return f'{_pair(ax)}-ax/{_strip(uniq)}-eq'
    if geom == 'OH':
        pairs = cf[1]
        same = [p for p in pairs if p[0] == p[1]]
        if len(same) == 0:
            return 'all-cis'
        if len(same) == len(pairs):
            return 'all-trans'
        # Iter-8.6j (2026-05-11): join ALL same-pair elements (not just
        # first) so MA3BC2-type fingerprints don't collapse to the same
        # base label.  Previously OH-1 (Cl-Cl + N-N both trans) and OH-2
        # (Cl-Cl trans only) both returned "Cl-trans" → downstream label-
        # dedup at ~24492 merged them, dropping FIRCOY/DEQVIE/TIYRUR to 2
        # frames despite group-theory predicting 4-6.  Joining via '+'
        # makes labels distinct: "Cl-trans" vs "Cl+N-trans" → survive
        # base-collapse. Env-gate DELFIN_OH_LABEL_RICH default 1 (active).
        if _delfin_env_int('DELFIN_OH_LABEL_RICH', 1):
            same_elems = sorted({_strip(p[0]) for p in same})
            return f"{'+'.join(same_elems)}-trans"
        return f'{_strip(same[0][0])}-trans'
    if geom == 'SQ':
        pairs = cf[1]
        same = sum(1 for p in pairs if p[0] == p[1])
        return 'trans' if same >= 1 else 'cis'
    if geom == 'COH':
        base_pairs, cap = cf[1], cf[2]
        pair_sig = ','.join(_pair(p) for p in base_pairs)
        return f'cap-{_strip(cap)}/{pair_sig}'
    if geom == 'SP':
        # cf = ('SP', apical_type, basal_pair_tuple)
        apical, basal_pairs = cf[1], cf[2]
        pair_sig = ','.join(_pair(p) for p in basal_pairs)
        return f'ap-{_strip(apical)}/{pair_sig}'
    if geom == 'SAP':
        pairs = cf[1]
        same = sum(1 for p in pairs if p[0] == p[1])
        return f'{same}-trans' if same else 'all-cis'
    if geom == 'DD':
        a_pairs, b_pairs = cf[1], cf[2]
        same_a = sum(1 for p in a_pairs if p[0] == p[1])
        same_b = sum(1 for p in b_pairs if p[0] == p[1])
        return f'A{same_a}/B{same_b}-trans'
    if geom == 'TPR':
        combined = cf[1]
        top, bottom = combined[0], combined[1]
        return f'top-{"".join(_strip(s) for s in top)}/bot-{"".join(_strip(s) for s in bottom)}'
    if geom == 'TTP':
        prism, caps = cf[1], cf[2]
        cap_sig = ''.join(_strip(s) for s in caps)
        return f'caps-{cap_sig}'
    return ''


def _extract_helicity_suffix(cf: tuple) -> str:
    """Iter-2: extract helicity suffix from a canonical form, if present.

    Returns 'L', 'D', or '' (no helicity tag).  Looks for a ``('chir', X)``
    sub-tuple anywhere in cf — appended by ``helicity_aware_pairs`` when
    DELFIN_CHIRAL_ENUM=1 and ≥ 2 chelate pairs.
    """
    if not cf:
        return ''
    for item in cf:
        if (isinstance(item, tuple) and len(item) == 2
                and item[0] == 'chir' and item[1] in ('L', 'D')):
            return item[1]
    return ''


# Trans position pairs (indices into the geometry vector list) for each geometry
_TOPO_TRANS_POSITIONS: Dict[str, List[Tuple[int, int]]] = {
    'LIN': [(0, 1)],
    'TP':  [],
    'TS':  [(0, 1)],
    'OH':  [(0, 1), (2, 3), (4, 5)],
    'SQ':  [(0, 1), (2, 3)],
    'TH':  [],               # tetrahedral has no trans pairs
    'TBP': [(0, 1)],          # only axial-axial is strictly 180°
    'SP':  [(1, 3), (2, 4)],  # trans basal pairs
    'PBP': [(0, 1)],          # axial-axial
    'SAP': [(0, 6), (1, 7), (2, 4), (3, 5)],  # square-antiprismatic
    'DD':  [(0, 2), (1, 3), (4, 6), (5, 7)],  # dodecahedral
    'SS':  [(0, 1)],          # see-saw: axial pair only
    'TPR': [],                # trigonal-prismatic: no 180° pairs
    'COH': [(0, 1), (2, 3), (4, 5)],  # capped-oh: octahedral base trans pairs
    'TTP': [],                # tricapped trigonal prism: no 180° pairs
    # Iter-2D — CN10/11/12 polyhedra (env DELFIN_CN_HIGH_ENABLE=1)
    'BCSAP': [(0, 2), (1, 3), (4, 6), (5, 7), (8, 9)],            # CN10 D4d
    'PAP':   [(0, 7), (1, 8), (2, 9), (3, 5), (4, 6)],            # CN10 D5d
    'CPAP':  [],                                                  # CN11 C5v: cap unique
    'ICOS':  [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)],  # CN12 Ih (antipodal)
    'CUBO':  [(0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11)],  # CN12 Oh
    'HBP':   [(0, 3), (1, 4), (2, 5), (6, 9), (7, 10), (8, 11)],  # CN12 D6h
}


_TOPO_CANONICAL_FNS = {
    'LIN': _canonical_lin,
    'TP':  _canonical_tp,
    'TS':  _canonical_ts,
    'OH':  _canonical_oh,
    'SQ':  _canonical_sq,
    'TH':  _canonical_th,
    'TBP': _canonical_tbp,
    'SP':  _canonical_sp,
    'PBP': _canonical_pbp,
    'SAP': _canonical_sap,
    'DD':  _canonical_dd,
    'SS':  _canonical_ss,
    'TPR': _canonical_tpr,
    'COH': _canonical_coh,
    'TTP': _canonical_ttp,
    # Iter-2D — CN10/11/12 polyhedra
    'BCSAP': _canonical_bcsap,
    'PAP':   _canonical_pap,
    'CPAP':  _canonical_cpap,
    'ICOS':  _canonical_icos,
    'CUBO':  _canonical_cubo,
    'HBP':   _canonical_hbp,
}


# Idealized coordination vectors per geometry (bond length ~2 Å)
_TOPO_GEOMETRY_VECTORS: Dict[str, List[Tuple[float, float, float]]] = {
    'LIN': [(2, 0, 0), (-2, 0, 0)],
    'TP':  [(2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'TS':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0)],
    'OH':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)],
    'SQ':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0)],
    'TH':  [(1.155, 1.155, 1.155), (-1.155, -1.155, 1.155),
            (-1.155, 1.155, -1.155), (1.155, -1.155, -1.155)],
    'TBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'SP':  [(0, 0, 2.2), (2, 0, 0.4), (0, 2, 0.4), (-2, 0, 0.4), (0, -2, 0.4)],
    'PBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (0.618, 1.902, 0),
            (-1.618, 1.176, 0), (-1.618, -1.176, 0), (0.618, -1.902, 0)],
    # Square-antiprismatic (D4d): top square at z=+0.8, bottom at z=-0.8, rotated 45°
    'SAP': [(1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
            (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)],
    # Dodecahedral (D2d): two interpenetrating trapezoids
    'DD':  [(1.414, 0, 1), (0, 1.414, -1), (-1.414, 0, 1), (0, -1.414, -1),
            (0.9, 0.9, 0), (-0.9, 0.9, 0), (-0.9, -0.9, 0), (0.9, -0.9, 0)],
    # See-saw / C2v (CN=4): axial pair along z, equatorial pair in xz-plane
    'SS':  [(0, 0, 2), (0, 0, -2), (2, 0, 0.4), (-2, 0, 0.4)],
    # Trigonal-prismatic (D3h, CN=6): eclipsed top/bottom triangles
    'TPR': [(2, 0, 1), (-1, 1.732, 1), (-1, -1.732, 1),
            (2, 0, -1), (-1, 1.732, -1), (-1, -1.732, -1)],
    # Capped octahedral (C3v, CN=7): octahedral base + capping 7th
    'COH': [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2),
            (1.155, 1.155, 1.155)],
    # Tricapped trigonal prism (D3h, CN=9): 6 prism + 3 equatorial caps
    'TTP': [(1.633, 0, 1.155), (-0.816, 1.414, 1.155), (-0.816, -1.414, 1.155),
            (1.633, 0, -1.155), (-0.816, 1.414, -1.155), (-0.816, -1.414, -1.155),
            (1.0, 1.732, 0), (-2.0, 0, 0), (1.0, -1.732, 0)],
    # Iter-2D — CN10 BCSAP (D4d): SAP + 2 axial caps
    'BCSAP': [(1.414, 0, 0.6), (0, 1.414, 0.6), (-1.414, 0, 0.6), (0, -1.414, 0.6),
              (1, 1, -0.6), (-1, 1, -0.6), (-1, -1, -0.6), (1, -1, -0.6),
              (0, 0, 2.0), (0, 0, -2.0)],
    # Iter-2D — CN10 PAP (D5d): two staggered pentagons (36° offset)
    'PAP':   [(1.902, 0.000, 0.8), (0.588, 1.809, 0.8), (-1.539, 1.118, 0.8),
              (-1.539, -1.118, 0.8), (0.588, -1.809, 0.8),
              (1.539, 1.118, -0.8), (-0.588, 1.809, -0.8), (-1.902, 0.000, -0.8),
              (-0.588, -1.809, -0.8), (1.539, -1.118, -0.8)],
    # Iter-2D — CN11 CPAP (C5v): PAP + axial cap (top side)
    'CPAP':  [(1.902, 0.000, 0.6), (0.588, 1.809, 0.6), (-1.539, 1.118, 0.6),
              (-1.539, -1.118, 0.6), (0.588, -1.809, 0.6),
              (1.539, 1.118, -1.0), (-0.588, 1.809, -1.0), (-1.902, 0.000, -1.0),
              (-0.588, -1.809, -1.0), (1.539, -1.118, -1.0),
              (0, 0, 2.0)],
    # Iter-2D — CN12 ICOS (Ih): 12 golden-ratio vertices, scaled to ~2 Å
    # Vertices: (0,±1,±φ), (±1,±φ,0), (±φ,0,±1) with φ=(1+√5)/2 ≈ 1.618.
    # Numbering chosen so vertex i and i+6 are antipodal (i ↔ -i).
    'ICOS':  [(0, 1.051, 1.701), (1.051, 1.701, 0), (1.701, 0, 1.051),
              (0, -1.051, 1.701), (-1.051, 1.701, 0), (1.701, 0, -1.051),
              (0, -1.051, -1.701), (-1.051, -1.701, 0), (-1.701, 0, -1.051),
              (0, 1.051, -1.701), (1.051, -1.701, 0), (-1.701, 0, 1.051)],
    # Iter-2D — CN12 CUBO (Oh): 12 cuboctahedral vertices = (±1,±1,0)+perms × √2
    # 3 mutually-orthogonal squares (xy, xz, yz).  Numbering: i and i+2 (mod 4)
    # are antipodal within each square.
    'CUBO':  [(1.414, 1.414, 0), (-1.414, 1.414, 0),
              (-1.414, -1.414, 0), (1.414, -1.414, 0),    # xy square
              (1.414, 0, 1.414), (-1.414, 0, 1.414),
              (-1.414, 0, -1.414), (1.414, 0, -1.414),    # xz square
              (0, 1.414, 1.414), (0, -1.414, 1.414),
              (0, -1.414, -1.414), (0, 1.414, -1.414)],   # yz square
    # Iter-2D — CN12 HBP (D6h): hexagonal prism, eclipsed top/bottom hexagons.
    'HBP':   [(2.0, 0, 1.0), (1.0, 1.732, 1.0), (-1.0, 1.732, 1.0),
              (-2.0, 0, 1.0), (-1.0, -1.732, 1.0), (1.0, -1.732, 1.0),
              (2.0, 0, -1.0), (1.0, 1.732, -1.0), (-1.0, 1.732, -1.0),
              (-2.0, 0, -1.0), (-1.0, -1.732, -1.0), (1.0, -1.732, -1.0)],
}


# ---------------------------------------------------------------------------
# Alternative binding-site exploration (Feature 5)
# ---------------------------------------------------------------------------

def _is_viable_donor(mol, atom_idx: int, metal_bonded: set) -> bool:
    """Check if an atom is a viable donor (has lone pairs available for coordination).

    Rejects atoms that are:
    - Double-bonded to carbon without additional lone pairs (C=O carbonyl)
    - Already bonded to a metal
    - Terminal hydrogen
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    sym = atom.GetSymbol()

    # Already bonded to a metal → not an *alternative* donor
    if any(nbr.GetSymbol() in _METAL_SET for nbr in atom.GetNeighbors()):
        return False

    # Standard heteroatom donors: N, O, S, P
    if sym in ('N', 'O', 'S', 'P'):
        # Count double bonds to carbon — carbonyl O (C=O) has only one neighbor
        # and that bond is a double bond, so it has no remaining lone pair for sigma donation
        if sym == 'O':
            nbrs = list(atom.GetNeighbors())
            if len(nbrs) == 1 and nbrs[0].GetSymbol() == 'C':
                bond = mol.GetBondBetweenAtoms(atom_idx, nbrs[0].GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                    return False  # carbonyl oxygen — not a viable donor
        return True

    # Carbon donors: cyclometalated aromatic C (ppy-type) or NHC carbens
    if sym == 'C':
        # Aromatic C that could be cyclometalated
        if atom.GetIsAromatic():
            return True
        # NHC-type carbene: C bonded to ≥2 N atoms
        n_nitrogen_nbrs = sum(
            1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'N'
        )
        if n_nitrogen_nbrs >= 2:
            return True

    return False


def _find_bridging_donors(mol) -> List[Tuple[int, List[int]]]:
    """Detect bridging donor atoms bound to ≥2 metal centers.

    Returns list of ``(donor_idx, [metal_idx_1, metal_idx_2, …])``.
    """
    bridging: List[Tuple[int, List[int]]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METAL_SET:
            continue
        metal_nbrs = [
            nbr.GetIdx() for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() in _METAL_SET
        ]
        if len(metal_nbrs) >= 2:
            bridging.append((atom.GetIdx(), metal_nbrs))
    return bridging


def _ligand_fragments(mol, metal_idx: int) -> List[set]:
    """Decompose non-metal atoms into ligand fragments (connected components after removing metal).

    Bridging donors (atoms bound to ≥2 metals) are included in ALL
    fragments they connect to, so each metal center's fragment list is
    complete.
    """
    non_metal = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in _METAL_SET}
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

    # Bridging donors: ensure they appear in every fragment whose metal
    # they coordinate to.  This prevents incorrect fragment splitting
    # for μ-Cl, μ-OR, μ-oxo bridged complexes.
    bridging = _find_bridging_donors(mol)
    if bridging:
        for donor_idx, metal_list in bridging:
            for frag in fragments:
                # If any atom in frag is a neighbour of donor_idx (via adj),
                # the donor should be in this fragment.
                if donor_idx not in frag and frag & adj.get(donor_idx, set()):
                    frag.add(donor_idx)

    return fragments
