"""Pyykko radii, donor sigma geometry and the UFF coordination constraints built from a template or from an XYZ in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Any, Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_int,
)
from delfin.manta.hapto_detect import (
    _classify_complex_class,
    _find_hapto_groups,
)
from delfin.manta.isomer_labels import (
    _TOPO_GEOMETRY_VECTORS,
)
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_METAL_BOND_LENGTHS,
    _METAL_SET,
    _PREFERRED_CN6_GEOMETRY,
    _get_ml_bond_length,
)
from delfin.manta.pre_uff_snap import (
    _D8_SQ_ISO_METALS,
    _D8_SQ_METALS,
)

logger = get_logger("delfin.smiles_converter")


# --- donor geometry, calibrated on CCDC clean_v2 -----------------------------
# 307370 structures / 1433955 sigma-donor records (exactly one metal partner,
# non-hapto, H-complete), measured in full rather than sampled.  Two findings drive
# everything below.
#
# 1. ONE angle per element, and it follows the PERIODIC ROW, not the group and not
#    electronegativity: Si 104.5 / P 103.7 / S 104.8 span 1.1 deg across EN 1.90-2.58,
#    while P 103.7 -> As 102.1 -> Sb 99.8 tracks the row exactly.  The same number
#    falls out of two independent bins -- theta from 4-partner X-D-X and theta from
#    3-partner (angle sum / 3) agree to ~1 deg (S 104.8/103.8, As 102.1/103.1,
#    Sb 99.8/99.1) -- so it is one physical parameter, not two fitted ones.
#
# 2. The metal angle is NOT free.  Widening M-D-X closes X-D-X and vice versa:
#       phi(M-D-X) = 109.5 - 0.87 * (theta - 109.5)
#    reproduces the measured medians to <= 0.5 deg for C, N, Si, P, S, Ge, As, Sn, Sb.
#    One parameter per element; every other angle is geometrically implied.
_THETA_BY_PERIOD = {2: 109.5, 3: 104.0, 4: 101.5, 5: 100.5, 6: 99.0}


_GROUP_14 = frozenset((6, 14, 32, 50, 82))   # C  Si Ge Sn Pb


_GROUP_15 = frozenset((7, 15, 33, 51, 83))   # N  P  As Sb Bi


_GROUP_16 = frozenset((8, 16, 34, 52))       # O  S  Se Te


def _period_of(z: int) -> int:
    for _p, _hi in ((1, 2), (2, 10), (3, 18), (4, 36), (5, 54), (6, 86)):
        if z <= _hi:
            return _p
    return 7


# Pyykko & Atsumi 2009 published SINGLE, DOUBLE and TRIPLE covalent radii together; DELFIN already
# uses the single set (occupier.load_covalent_radii, source "pyykko2009").  Loading the other two
# turns the bond-order contraction from a global scale factor into an ATOM property, defined for
# every element instead of a hand-listed handful of pairs.
_PYYKKO_ORDER_ATTR = {2: "covalent_radius_pyykko_double", 3: "covalent_radius_pyykko_triple"}


_PYYKKO_ORDER_CACHE: Dict[Tuple[str, int], Optional[float]] = {}


def _pyykko_order_radius(sym: str, order: int) -> Optional[float]:
    """Pyykko 2009 double/triple covalent radius (A) for `sym`, or None if unavailable.

    Looked up per element and cached, so a build only ever queries the handful of elements it
    actually contains.  None makes every caller a no-op -- never-worse by construction if
    mendeleev is missing or has no value for that element.
    """
    key = (sym, order)
    if key in _PYYKKO_ORDER_CACHE:
        return _PYYKKO_ORDER_CACHE[key]
    val: Optional[float] = None
    attr = _PYYKKO_ORDER_ATTR.get(order)
    if attr is not None:
        try:
            from mendeleev import element  # type: ignore
            raw = getattr(element(sym), attr, None)
            if raw is not None:
                val = float(raw) / 100.0     # mendeleev stores pm
        except Exception:
            val = None
    _PYYKKO_ORDER_CACHE[key] = val
    return val


def _donor_sigma_geometry(atom):
    """Geometry target for a non-metal centre, counting the METAL as a sigma partner.

    2026-07-29.  Every hybridisation decision in the build path reads RDKit off a
    graph the metal was stripped from (`decompose.py` removes each M-D bond and
    re-sanitises the fragment), and `_build_uff_constraints_from_template` then drops
    bonded metals a second time when it collects `heavy_nbrs`.  A coordinated donor
    is therefore short one sigma partner, and `len(heavy_nbrs) < 3` skips it outright:

        R2N-M      R, R, M      -> 2 heavy non-metal -> no restraint at all, so the
                                   metal drifts out of the amido pi plane (15-23 deg
                                   measured on AFOMEO frame 4)
        pyridine-M C, C, M      -> 2                 -> metal leaves the ring plane
        M-NH2-R    H, H, R, M   -> 1                 -> no tetrahedral target either

    Counting the metal fixes both directions with one rule, atom- and bond-specific:

        >= 4 sigma partners            -> tetrahedral   (M-NR3, M-PR3, M-NH2-R)
        == 3 sigma partners, period 2  -> planar        (amido, aryl-C, carbene,
                                          pyridine-N: the metal lies IN the plane)
        == 3 sigma partners, period 3+ -> pyramidal     (phosphido, thioether,
                                          sulfoxide-S: inversion barrier wins)
        <= 2 sigma partners            -> no opinion    (alkoxide, thiolate, nitrile)

    Returns (mode, sigma_neighbour_indices, theta, phi); mode is
    "planar"    -> only the improper dihedral (the RING sets the individual angles,
                   see below), "pyramidal" -> theta among non-metal partners and phi
                   for every pair involving the metal, no improper, or
    None        -> no opinion, leave the caller's behaviour untouched.

    NOTE: a planar 3-partner donor must NOT be pushed to 120/120/120.  Only the SUM
    is invariant; the ring fixes the split.  Measured on N: flat 6-ring 118.1 internal
    / 120.8 to the metal, flat 5-ring 106.1 / 126.6, acyclic 118.3 / 120.6.  Forcing
    120 would bend every imidazole and pyrazole donor by ~14 deg.
    """
    sigma = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    if not any(n.GetSymbol() in _METAL_SET for n in sigma):
        return None, [], 0.0, 0.0
    z = atom.GetAtomicNum()
    theta = _THETA_BY_PERIOD.get(_period_of(z), 99.0)
    phi = 109.5 - 0.87 * (theta - 109.5)
    idx = sorted(n.GetIdx() for n in sigma)
    n_h = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
    n_h += atom.GetTotalNumHs()
    n_sigma = len(sigma) + n_h
    if n_sigma >= 4:
        # theta/phi collapse to 109.5/109.5 for period 2, i.e. plain tetrahedral.
        return "pyramidal", idx, theta, phi
    if n_sigma != 3:
        return None, idx, theta, phi
    if _period_of(z) == 2:
        # O is genuinely bimodal in the crystal -- 51 % planar, 26 % pyramidal, p10 of
        # the angle sum 332.9 deg.  Coordinated ether / alkoxide / aqua is flattened
        # but SOFT; forcing either target would be wrong, so decline to have an opinion.
        if z == 8:
            return None, idx, theta, phi
        return "planar", idx, theta, phi
    if z in _GROUP_14:
        # Three sigma partners on a heavy group-14 centre is a carbene analogue with an
        # M=E multiple bond, not a lone-pair donor: Si 96 %, Sn 79 % planar.
        return "planar", idx, theta, phi
    if z in _GROUP_15:
        # Genuinely bimodal (P: 50.4 % planar / 43.9 % pyramidal).  What separates them
        # is whether the donor sits in a coplanar ring: phosphinine 95.4 % planar,
        # a free phosphido pyramidal.
        return ("planar" if atom.GetIsAromatic() else "pyramidal"), idx, theta, phi
    # Group 16 stays pyramidal even inside an aromatic ring: S in a flat 5-ring is
    # 94.8 % PYRAMIDAL (sum 316.2 deg) and holds the metal 2.0 A off the ring plane.
    # "Conjugation flattens the donor" is measurably FALSE here (S 0.4 % vs 0.3 %).
    return "pyramidal", idx, theta, phi


def _build_uff_constraints_from_template(
    mol_template,
    xyz_delfin: Optional[str] = None,
) -> Optional[Dict]:
    """Build conservative OB-UFF constraints from the RDKit template graph.

    The goal is not to freeze the structure, but to keep fragile motifs
    chemically reasonable during UFF relaxation:
    - carboxyl-like O-C-O groups near trigonal-planar geometry
    - aromatic ring torsions near planarity
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return None

    constraints: Dict[str, list] = {
        "fix_atoms": [],
        "distances": [],
        "angles": [],
        "torsions": [],
    }
    seen_dist: set = set()
    seen_angle: set = set()
    seen_tors: set = set()

    # Optional source coordinates (same atom ordering as mol_template) to pick
    # the closest planar torsion target (0 or 180 deg).
    coords: Optional[List[Tuple[float, float, float]]] = None
    if xyz_delfin:
        try:
            lines = [l for l in xyz_delfin.splitlines() if l.strip()]
            if len(lines) == mol_template.GetNumAtoms():
                parsed: List[Tuple[float, float, float]] = []
                for line in lines:
                    parts = line.split()
                    if len(parts) < 4:
                        parsed = []
                        break
                    parsed.append((float(parts[1]), float(parts[2]), float(parts[3])))
                if len(parsed) == mol_template.GetNumAtoms():
                    coords = parsed
        except Exception:
            coords = None

    def _nearest_planar_target(a: int, b: int, c: int, d: int) -> float:
        if coords is None:
            return 0.0
        try:
            p1 = coords[a]
            p2 = coords[b]
            p3 = coords[c]
            p4 = coords[d]
            b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
            b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
            b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

            n1 = (
                b1[1] * b2[2] - b1[2] * b2[1],
                b1[2] * b2[0] - b1[0] * b2[2],
                b1[0] * b2[1] - b1[1] * b2[0],
            )
            n2 = (
                b2[1] * b3[2] - b2[2] * b3[1],
                b2[2] * b3[0] - b2[0] * b3[2],
                b2[0] * b3[1] - b2[1] * b3[0],
            )
            n1m = math.sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2])
            n2m = math.sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2])
            if n1m < 1e-10 or n2m < 1e-10:
                return 0.0
            n1u = (n1[0] / n1m, n1[1] / n1m, n1[2] / n1m)
            n2u = (n2[0] / n2m, n2[1] / n2m, n2[2] / n2m)
            dot = max(-1.0, min(1.0, n1u[0] * n2u[0] + n1u[1] * n2u[1] + n1u[2] * n2u[2]))
            dihedral = math.degrees(math.acos(dot))
            # Planar torsions are near 0 or 180.
            return 180.0 if abs(dihedral - 180.0) < abs(dihedral - 0.0) else 0.0
        except Exception:
            return 0.0

    try:
        # Carboxyl/carbonyl-like O-C-O units: keep broad C-O distances and ~120 deg angle.
        for atom in mol_template.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue
            o_neighbors = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() == 8 and n.GetSymbol() not in _METAL_SET
            ]
            if len(o_neighbors) != 2:
                continue

            c_idx = atom.GetIdx()
            o1_idx, o2_idx = sorted([o_neighbors[0].GetIdx(), o_neighbors[1].GetIdx()])

            angle_key = (o1_idx, c_idx, o2_idx)
            if angle_key not in seen_angle:
                seen_angle.add(angle_key)
                constraints["angles"].append((o1_idx, c_idx, o2_idx, 120.0))

            for o_idx in (o1_idx, o2_idx):
                dist_key = tuple(sorted((c_idx, o_idx)))
                if dist_key in seen_dist:
                    continue
                seen_dist.add(dist_key)
                bond = mol_template.GetBondBetweenAtoms(c_idx, o_idx)
                target = 1.27
                if bond is not None:
                    btype = bond.GetBondType()
                    if btype == Chem.BondType.DOUBLE:
                        target = 1.22
                    elif btype == Chem.BondType.SINGLE:
                        target = 1.31
                    elif btype == Chem.BondType.AROMATIC:
                        target = 1.27
                constraints["distances"].append((c_idx, o_idx, target))

            # Keep carboxyl group coplanar with its carbon backbone when possible.
            x_neighbors = sorted(
                n.GetIdx()
                for n in atom.GetNeighbors()
                if n.GetIdx() not in (o1_idx, o2_idx)
                and n.GetAtomicNum() > 1
                and n.GetSymbol() not in _METAL_SET
            )
            for x_idx in x_neighbors:
                x_atom = mol_template.GetAtomWithIdx(x_idx)
                anchor_candidates = sorted(
                    n.GetIdx()
                    for n in x_atom.GetNeighbors()
                    if n.GetIdx() != c_idx
                    and n.GetAtomicNum() > 1
                    and n.GetSymbol() not in _METAL_SET
                )
                if not anchor_candidates:
                    continue
                a_idx = anchor_candidates[0]
                for o_idx in (o1_idx, o2_idx):
                    key = (a_idx, x_idx, c_idx, o_idx)
                    rev = (o_idx, c_idx, x_idx, a_idx)
                    tors_key = key if key <= rev else rev
                    if tors_key in seen_tors:
                        continue
                    seen_tors.add(tors_key)
                    constraints["torsions"].append(
                        (a_idx, x_idx, c_idx, o_idx, _nearest_planar_target(a_idx, x_idx, c_idx, o_idx))
                    )

        # Aromatic ring planarity via torsion constraints on consecutive ring quartets.
        try:
            Chem.GetSymmSSSR(mol_template)
        except Exception:
            pass
        ring_info = mol_template.GetRingInfo()
        if ring_info is not None:
            for ring in ring_info.AtomRings():
                if len(ring) < 5 or len(ring) > 7:
                    continue
                if not all(mol_template.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    continue
                if any(mol_template.GetAtomWithIdx(i).GetSymbol() in _METAL_SET for i in ring):
                    continue

                n = len(ring)
                for i in range(n):
                    a = ring[(i - 1) % n]
                    b = ring[i]
                    c = ring[(i + 1) % n]
                    d = ring[(i + 2) % n]
                    if len({a, b, c, d}) < 4:
                        continue
                    key = (a, b, c, d)
                    rev = (d, c, b, a)
                    tors_key = key if key <= rev else rev
                    if tors_key in seen_tors:
                        continue
                    seen_tors.add(tors_key)
                    constraints["torsions"].append((a, b, c, d, _nearest_planar_target(a, b, c, d)))

        # MULTIBOND LENGTH (2026-07-29, env, default OFF): pin every heavy-heavy DOUBLE/TRIPLE bond
        # to the order-specific covalent-radii sum.  The carboxyl block above already does exactly
        # this -- for one atom type, with three hardcoded C-O numbers, and those numbers ARE the
        # Pyykko sums (C=O 0.67+0.57 = 1.24 vs its 1.22; C#O 0.60+0.53 = 1.13).  So the special case
        # is a hand-listed subset of a rule that holds for the whole periodic table.
        #
        # Measured gap it closes (CCDC clean_v2, 25.8 M bonds, signature-resolved bands): terminal
        # multiply-bonded heteroatoms on high-coordination centres are built at SINGLE-bond length --
        # ClO4- Cl-O +0.55 A (100 % outside the CSD band), N=N=N +0.295, S=O +0.20, coordinated C#N
        # +0.18 (472 systems), C#S +0.17 -- while P=O (+0.055) and N-O (-0.025) are already fine.
        # A wrong REFERENCE, not noise, so a reference is what it needs.
        #
        # Runs AFTER the carboxyl block on purpose: `seen_dist` makes that block win where it already
        # acts, so the flag only ever ADDS pins for bonds nothing constrained before.  Aromatic bonds
        # are left alone (measured offset +0.024 A -- not worth the blast radius), and M-L lengths
        # belong to the coordination path.
        if os.environ.get("DELFIN_FFFREE_MULTIBOND_LEN", "0") == "1":
            for _bond in mol_template.GetBonds():
                _bt = _bond.GetBondType()
                if _bt == Chem.BondType.TRIPLE:
                    _order = 3
                elif _bt == Chem.BondType.DOUBLE:
                    _order = 2
                else:
                    continue
                _a1, _a2 = _bond.GetBeginAtom(), _bond.GetEndAtom()
                if _a1.GetAtomicNum() <= 1 or _a2.GetAtomicNum() <= 1:
                    continue
                _s1, _s2 = _a1.GetSymbol(), _a2.GetSymbol()
                if _s1 in _METAL_SET or _s2 in _METAL_SET:
                    continue
                _r1 = _pyykko_order_radius(_s1, _order)
                _r2 = _pyykko_order_radius(_s2, _order)
                if _r1 is None or _r2 is None:
                    continue          # unknown element -> behave exactly as today
                _dk = tuple(sorted((_a1.GetIdx(), _a2.GetIdx())))
                if _dk in seen_dist:
                    continue
                seen_dist.add(_dk)
                constraints["distances"].append((_dk[0], _dk[1], round(_r1 + _r2, 3)))

        # Sp2 planarity at every 3-coordinate planar atom.  A ring-wise loop
        # alone cannot keep a fused-ring junction atom planar because the
        # two rings' torsion sets are independent and allow the shared atom
        # to fold between them.  Pinning one improper dihedral per sp2 atom
        # (A-X-B-C with X central) forces its three heavy neighbours to
        # stay coplanar with X for every ligand topology.
        #
        # Planarity alone however does not constrain which *angles* the
        # three neighbours make with X — a degenerate sp2 (e.g. nitrate
        # where two Os stack to close a 4-ring with the metal) can be
        # strictly coplanar yet have O-N-O angles of 0°/125°/126° instead
        # of 3 × 120°.  Adding pairwise angle constraints (120° for sp2,
        # 109.5° for sp3 X with ≥3 same-element neighbours) forces a
        # regular triangle/tetrahedron and generalises the carboxyl
        # O-C-O rule above to every rigid anion / cation (NO3-, SO4²-,
        # PO4³-, CO3²-, carbamate, urea, guanidinium, ...).
        _metal_sigma_count = (
            os.environ.get("DELFIN_FFFREE_METAL_SIGMA_COUNT", "0") == "1"
        )
        # MODE (2026-07-29): "full" pins a planar donor with an improper dihedral that includes the
        # metal; "angles" drops that dihedral and steers the metal with ANGLE targets instead.
        # full:1000 showed the improper fights the coordination polyhedron -- poly_cshm_vs_ccdc was
        # the worst-regressing axis, with coord_angle and graph_geom right behind, i.e. the damage
        # landed on angles AT THE METAL, not on the donor.  With a chelate, every donor demands the
        # metal in ITS plane while the polyhedron demands its own L-M-L angles; something gives, and
        # it was the polyhedron (ccdc_isomer_lost 2, topology_floor_ok false).
        _metal_sigma_mode = os.environ.get("DELFIN_FFFREE_METAL_SIGMA_MODE", "full")
        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() in _METAL_SET:
                continue
            if atom.GetAtomicNum() <= 1:
                continue
            is_sp2 = (
                atom.GetIsAromatic()
                or atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2
            )
            is_sp3_four_same = False
            if not is_sp2:
                hyb = atom.GetHybridization()
                if hyb == Chem.rdchem.HybridizationType.SP3:
                    _o_nbrs = [
                        n for n in atom.GetNeighbors()
                        if n.GetAtomicNum() == 8 and n.GetSymbol() not in _METAL_SET
                    ]
                    if len(_o_nbrs) >= 3:
                        is_sp3_four_same = True
            # METAL-SIGMA (2026-07-29, env, default OFF): count a bonded metal as a
            # sigma partner so a coordinated donor gets the geometry target its REAL
            # coordination number implies, instead of being dropped below for having
            # "too few" neighbours once the metal is subtracted.  See
            # _donor_sigma_geometry for the atom-specific rule.
            _geom = None
            _sigma_idx: List[int] = []
            _theta = _phi = 0.0
            if _metal_sigma_count:
                _geom, _sigma_idx, _theta, _phi = _donor_sigma_geometry(atom)
                if _geom is not None:
                    is_sp2 = _geom == "planar"
                    is_sp3_four_same = False
            if not (is_sp2 or is_sp3_four_same or _geom is not None):
                continue
            if _geom is not None:
                # sigma partners INCLUDING the metal -- the improper below then pins
                # the metal into the donor plane, which is the whole point.
                heavy_nbrs = _sigma_idx
            else:
                heavy_nbrs = sorted(
                    n.GetIdx()
                    for n in atom.GetNeighbors()
                    if n.GetAtomicNum() > 1 and n.GetSymbol() not in _METAL_SET
                )
            # Pairwise angles need two partners, the improper needs three.  The legacy
            # path keeps its 3-partner floor so it stays byte-identical.
            if len(heavy_nbrs) < (2 if _geom is not None else 3):
                continue
            x = atom.GetIdx()
            # In "angles" mode a planar donor is steered by angle targets alone -- no improper, so
            # nothing forces the metal into the donor plane against the polyhedron.
            # MODE "theta" (2026-07-29, third bisection step): constrain ONLY the donor's internal
            # X-D-X angle and emit nothing that involves the metal -- no improper, no phi.  Both
            # earlier modes damaged the polyhedron (poly_cshm_vs_ccdc was the worst-regressing axis
            # in BOTH, and "angles" was worse than "full": ccdc_isomer_lost 8 vs 2, isomers_lost 30
            # vs 9).  Dropping the improper therefore was not the answer; what the two modes still
            # SHARE is that they constrain angles AT the metal.  The coordination path already owns
            # M-D distances and D-M-D angles, so a second set of metal-involving targets competes
            # with it.  This mode shapes the LIGAND and leaves the metal to the path that owns it.
            _theta_only = (_geom is not None and _metal_sigma_mode == "theta")
            _planar_angles = (_geom == "planar" and _metal_sigma_mode == "angles")
            if is_sp2 and len(heavy_nbrs) >= 3 and not _planar_angles and not _theta_only:
                # Improper dihedral → planarity
                a, b, c = heavy_nbrs[:3]
                key = (a, x, b, c)
                rev = (c, b, x, a)
                tors_key = key if key <= rev else rev
                if tors_key not in seen_tors:
                    seen_tors.add(tors_key)
                    constraints["torsions"].append(
                        (a, x, b, c, _nearest_planar_target(a, x, b, c))
                    )
            _use_ring_angles = (_geom == "planar") and (_planar_angles or _theta_only)
            if _geom == "planar" and not _use_ring_angles:
                # "full" mode: deliberately NO angle constraints.  Only the angle SUM is invariant
                # for a planar donor; the ring sets the split, and pinning 120 would bend every
                # imidazole and pyrazole donor by ~14 deg.  The improper above delivers planarity.
                continue
            # Ring-aware targets for a planar donor -- pure polygon geometry, no fitted constant:
            # a planar n-ring has internal angle 180 - 360/n, and the two exocyclic angles split
            # what is left of 360.  6-ring -> 120 internal / 120 exocyclic, 5-ring -> 108 / 126.
            # Cross-check against the crystal: measured 118.0 / 120.75 and 106.0 / 126.5, so the
            # regular-polygon values land within ~2 deg without importing a single measured number.
            _ring_internal = _ring_exo = 0.0
            if _use_ring_angles:
                _rs = 0
                for _n in (3, 4, 5, 6, 7, 8):
                    if atom.IsInRingSize(_n):
                        _rs = _n
                        break
                if _rs >= 3:
                    _ring_internal = 180.0 - 360.0 / _rs
                    _ring_exo = (360.0 - _ring_internal) / 2.0
                else:
                    _ring_internal = _ring_exo = 120.0     # acyclic planar centre
            # Pairwise angle constraints on every heavy-heavy pair.
            target_angle = 120.0 if is_sp2 else 109.5
            for i in range(len(heavy_nbrs)):
                for j in range(i + 1, len(heavy_nbrs)):
                    ai, aj = heavy_nbrs[i], heavy_nbrs[j]
                    if _theta_only:
                        _has_m = (
                            mol_template.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET
                            or mol_template.GetAtomWithIdx(aj).GetSymbol() in _METAL_SET
                        )
                        if _has_m:
                            continue          # the coordination path owns everything at the metal
                        target_angle = _ring_internal if _use_ring_angles else _theta
                    elif _planar_angles:
                        # Internal iff BOTH partners share a ring bond with the donor.
                        _bi = mol_template.GetBondBetweenAtoms(x, ai)
                        _bj = mol_template.GetBondBetweenAtoms(x, aj)
                        _both_ring = bool(
                            _bi is not None and _bj is not None
                            and _bi.IsInRing() and _bj.IsInRing()
                        )
                        target_angle = _ring_internal if _both_ring else _ring_exo
                    elif _geom is not None:
                        # theta among non-metal partners, phi for any pair involving
                        # the metal -- one calibrated parameter per element, the rest
                        # geometrically implied (see _donor_sigma_geometry).
                        _has_m = (
                            mol_template.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET
                            or mol_template.GetAtomWithIdx(aj).GetSymbol() in _METAL_SET
                        )
                        target_angle = _phi if _has_m else _theta
                    angle_key = (ai, x, aj)
                    if angle_key in seen_angle:
                        continue
                    seen_angle.add(angle_key)
                    constraints["angles"].append((ai, x, aj, target_angle))

        # Metal-metal bond distance constraints.
        for bond in mol_template.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetSymbol() not in _METAL_SET or a2.GetSymbol() not in _METAL_SET:
                continue
            i1, i2 = a1.GetIdx(), a2.GetIdx()
            dist_key = tuple(sorted((i1, i2)))
            if dist_key in seen_dist:
                continue
            seen_dist.add(dist_key)
            mm_key = frozenset({a1.GetSymbol(), a2.GetSymbol()})
            target = _METAL_METAL_BOND_LENGTHS.get(mm_key)
            if target is None:
                r1 = _COVALENT_RADII.get(a1.GetSymbol())
                r2 = _COVALENT_RADII.get(a2.GetSymbol())
                if r1 is not None and r2 is not None:
                    target = r1 + r2 + 0.3
                else:
                    target = 2.5
            constraints["distances"].append((i1, i2, float(target)))
    except Exception:
        return None

    if not (
        constraints["fix_atoms"]
        or constraints["distances"]
        or constraints["angles"]
        or constraints["torsions"]
    ):
        return None
    return constraints


def _detect_chelate_donors(mol, metal_idx: int, donor_indices: List[int]) -> set:
    """Return the subset of ``donor_indices`` that share a ligand body
    with another donor of the same metal (i.e. chelate donors).

    A donor is chelate if, after removing the metal, another donor in
    ``donor_indices`` is still reachable through the molecular graph.
    Such donors are constrained by their ring closure and must NOT be
    frozen during UFF, otherwise UFF cannot relax ring internals and a
    ligand atom often gets dragged into the metal coordination sphere.
    """
    if not donor_indices:
        return set()
    donor_set = set(donor_indices)
    chelate: set = set()
    for start in donor_indices:
        if start in chelate:
            continue
        visited = {start, metal_idx}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in visited or ni == metal_idx:
                    continue
                visited.add(ni)
                stack.append(ni)
                if ni in donor_set and ni != start:
                    chelate.add(start)
                    chelate.add(ni)
    return chelate


def _build_coordination_constraints_from_xyz(
    mol_template,
    xyz_delfin: str,
    d8_trans=None,
    suppress_d8_sq: bool = False,
    force_d8_sq: bool = False,
    suppress_cn6_oh: bool = False,
    force_cn6_oh: bool = False,
) -> Optional[Dict]:
    """Auto-detect metal coordination from template graph and pin it during UFF.

    ``d8_trans`` (optional): the EXACT per-isomer trans donor-atom-index pairs
    from the topology enumerator (perm + _TOPO_TRANS_POSITIONS['SQ']).  When
    given and DELFIN_FFFREE_D8_SQ_ISO is on, the d8 square is imposed on THESE
    pairs (authoritative — no geometric guessing), which rescues strained
    macrocyclic frames the geometry alone cannot and cannot collapse isomers.
    Default None -> geometry fallback (byte-identical when the flag is off).

    ``suppress_d8_sq`` (default False): force the d8 SP-4 imposition OFF for this
    call even when DELFIN_FFFREE_D8_SQ_ISO is on.

    ``force_d8_sq`` (default False): impose the d8 SP-4 square for THIS call even when
    DELFIN_FFFREE_D8_SQ_ISO is off.  Used by the ADDITIVE d8 pass (DELFIN_FFFREE_D8_SQ_ADD)
    so the whole feature is ONE self-contained axis: the primary frame is the normal
    (tetrahedral) build and the SP-4 frame is added as a PURELY ADDITIVE sibling.  A bulky
    ligand that clashes under SP-4 keeps its valid tetrahedral primary (never a regression);
    a d8-no-valid system gains its valid SP-4 frame.  suppress overrides force.

    Unlike ``_build_coordination_uff_constraints`` (which needs explicit
    perm/geometry arguments from the topology enumerator), this function
    derives constraints directly from the template bond graph.  This
    makes it usable for ANY UFF call — sampling, linkage, alt-binding,
    hapto, etc.

    Only M-D DISTANCES are constrained (from the lookup table).  L-M-L
    angles are left free so UFF can optimize toward ideal polyhedral
    geometry.  This combines Principle 2 (distances from chemistry) with
    Principle 3 (UFF improves angles toward ideal).
    """
    base = _build_uff_constraints_from_template(mol_template, xyz_delfin=xyz_delfin)
    if base is None:
        base = {"fix_atoms": [], "distances": [], "angles": [], "torsions": []}

    # Baustein-5+6 Phase 3: collect metadata for opt-in soft-donor UFF mode.
    # Populated only when DELFIN_UFF_SOFT_DONORS=1 consumer reads this meta.
    # Format kept minimal so legacy callers that ignore the key are unaffected.
    _soft_meta: Dict[str, Any] = {
        "metal_indices": [],
        "donor_indices": [],     # monodentate (non-chelate) donors only
        "pairs": [],             # (m_idx, d_idx) tuples
        "class_label": "no_metal",
    }

    if not RDKIT_AVAILABLE or mol_template is None or not xyz_delfin:
        return base if any(base.values()) else None

    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return base if any(base.values()) else None
        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return base if any(base.values()) else None
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        seen_dist = {tuple(sorted((a, b))) for a, b, _t in base["distances"]}
        seen_angle = set()
        for a, b, c, _t in base["angles"]:
            seen_angle.add((a, b, c))
            seen_angle.add((c, b, a))

        try:
            _soft_meta["class_label"] = _classify_complex_class(mol_template)
        except Exception:
            _soft_meta["class_label"] = "no_metal"

        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = atom.GetIdx()
            m_sym = atom.GetSymbol()
            _soft_meta["metal_indices"].append(m_idx)

            donor_indices: List[int] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() <= 1:
                    continue
                donor_indices.append(nbr.GetIdx())

            # M-D distance constraints from lookup table.
            # Iter-8.6h (2026-05-11): env-flag DELFIN_UFF_NO_MD_CONSTRAINT
            # default 0.  When 1, M-D distance constraint is NOT added —
            # UFF optimizes M-D freely, allowing convergence to UFF
            # equilibrium (typically shorter than lookup-table values,
            # closer to 6efa34e/81f8a1f champion Fe-O 1.887 vs HEAD 1.999).
            # M-D invariant still protected by downstream
            # _verify_metal_connectivity (donor-drift catch).
            _skip_md_dist = bool(_delfin_env_int("DELFIN_UFF_NO_MD_CONSTRAINT", 0))
            if not _skip_md_dist:
                for d_idx in donor_indices:
                    key = tuple(sorted((m_idx, d_idx)))
                    if key in seen_dist:
                        continue
                    seen_dist.add(key)
                    d_sym = mol_template.GetAtomWithIdx(d_idx).GetSymbol()
                    if d_sym in _METAL_SET:
                        continue
                    bl = float(_get_ml_bond_length(m_sym, d_sym))
                    base["distances"].append((m_idx, d_idx, bl))

            # FREEZE metal + monodentate-donor atoms during UFF.
            # The topology builder placed donors at the ideal symmetric
            # _TOPO_GEOMETRY_VECTORS positions (Oh/Td/D4h/D3h/SAP/...).
            # Freezing M + donors keeps the high-symmetry polyhedron
            # exactly intact while still letting UFF relax ligand
            # internals (ring planarity, torsions, bond lengths beyond
            # the donor).
            #
            # CHELATE EXCEPTION: a donor that shares a ligand body with
            # another donor of the same metal must NOT be frozen.  The
            # topology builder's chelate-ring placement is geometric-
            # constraint-driven (ring closure), not pure ideal-polyhedron;
            # freezing leaves UFF unable to relax ring internals around
            # the fixed donor, and a ligand atom often gets dragged into
            # the metal coordination sphere, producing a spurious bond
            # and topology break.  Verified on D-IROWUK_4d_Rh_CN2
            # (tetrapyridyl): freeze causes max_dev 60°, no-freeze gives
            # max_dev 12.7° on Oh.
            chelate_donors = _detect_chelate_donors(
                mol_template, m_idx, donor_indices
            )
            if m_idx not in base["fix_atoms"]:
                base["fix_atoms"].append(m_idx)
            # Iter-8.6e (2026-05-11): env-flag controlled donor-relax.
            # When DELFIN_UFF_RELAX_DONORS=1, monodentate donors are NOT
            # frozen during UFF.  M-D distance constraints (added above)
            # keep bond lengths near the lookup-table values, but L-M-L
            # angles can relax based on donor-mixing / trans-influence
            # chemistry.  Validated: CESNIZ Ru CN6 4 -> 8 unique isomers,
            # CPOENR Re CN6 1 -> 4 unique isomers (single-SMILES per-frame
            # tests).
            #
            # Iter-8.6i DEFAULT-FLIP REVERTED 2026-05-11: smoke500 with
            # env=1 default showed sigma topo_pct -2.8pp regression
            # (84.9% -> 82.1%) with NO compensating hapto/multi-hapto
            # gain.  Default kept OFF (0); opt-in via env=1 for the
            # per-SMILES isomer recovery (CESNIZ/CPOENR pattern).
            _uff_relax_donors = bool(
                _delfin_env_int("DELFIN_UFF_RELAX_DONORS", 0)
            )
            # d8 SQUARE-PLANAR angle targets (env-gated DELFIN_FFFREE_D8_SQ_ANGLES): this function
            # leaves L-M-L angles FREE, and UFF's default for CN4 is tetrahedral -> a d8 CN4
            # (Pd/Pt/Ni/Au/Rh/Ir) built by ETKDG comes out tetrahedral/seesaw.  Give UFF explicit
            # 90/180 angle targets so it optimises toward SQUARE-planar, and do NOT freeze the donors
            # (below) so UFF can move them there.  This works WITH UFF -- it arranges bulky ligands
            # (no rigid-placement clash) and relaxes chelate backbones (handles the chelate majority
            # the geometric flatten could not).  Default OFF -> byte-identical.
            _d8_sq_angles = (os.environ.get("DELFIN_FFFREE_D8_SQ_ANGLES", "0") == "1"
                             and m_sym in _D8_SQ_METALS and len(donor_indices) == 4)
            if _d8_sq_angles:
                try:
                    import numpy as _np
                    _mp = _np.array(coords[m_idx], dtype=float)
                    _dv = {di: _np.array(coords[di], dtype=float) - _mp for di in donor_indices}
                    _dn = {di: (_dv[di] / (float(_np.linalg.norm(_dv[di])) + 1e-12)) for di in donor_indices}
                    _rem = list(donor_indices)
                    _trans_set = set()
                    while len(_rem) >= 2:                    # pair each donor with its most-opposite
                        a = _rem[0]
                        b = min(_rem[1:], key=lambda x: float(_np.dot(_dn[a], _dn[x])))
                        base["angles"].append((a, m_idx, b, 180.0))
                        _trans_set.add((a, b)); _trans_set.add((b, a))
                        _rem.remove(a); _rem.remove(b)
                    for _ii in range(len(donor_indices)):    # all remaining pairs are cis = 90
                        for _jj in range(_ii + 1, len(donor_indices)):
                            _a, _b = donor_indices[_ii], donor_indices[_jj]
                            if (_a, _b) not in _trans_set:
                                base["angles"].append((_a, m_idx, _b, 90.0))
                except Exception:
                    _d8_sq_angles = False
            # d8 SQUARE-PLANAR, PERM path (env-gated DELFIN_FFFREE_D8_SQ_ISO): when the topology
            # enumerator handed us the EXACT per-isomer trans pairs (perm + _TOPO_TRANS_POSITIONS['SQ']),
            # impose the square on THOSE pairs -- authoritative, no geometric guessing.  This rescues the
            # strained macrocyclic frames whose geometry alone is too ambiguous, and cannot collapse
            # isomers (each isomer's OWN trans is imposed).  Sets _d8_sq_angles so the geometry fallback
            # below is skipped (its `not _d8_sq_angles` guard).
            if (not _d8_sq_angles and d8_trans and not suppress_d8_sq
                    and (os.environ.get("DELFIN_FFFREE_D8_SQ_ISO", "0") == "1" or force_d8_sq)
                    and m_sym in _D8_SQ_ISO_METALS and len(donor_indices) == 4):
                try:
                    _dset = set(donor_indices)
                    _pt = [(int(_a), int(_b)) for (_a, _b) in d8_trans
                           if int(_a) in _dset and int(_b) in _dset]
                    _flat = [_x for _pr in _pt for _x in _pr]
                    if len(_pt) == 2 and len(set(_flat)) == 4:
                        _trans_set = set()
                        for (_a, _b) in _pt:
                            base["angles"].append((_a, m_idx, _b, 180.0))
                            _trans_set |= {(_a, _b), (_b, _a)}
                        for _ii in range(4):
                            for _jj in range(_ii + 1, 4):
                                _x, _y = donor_indices[_ii], donor_indices[_jj]
                                if (_x, _y) not in _trans_set:
                                    base["angles"].append((_x, m_idx, _y, 90.0))
                        _d8_sq_angles = True
                except Exception:
                    pass
            # d8 SQUARE-PLANAR, ISOMER-PRESERVING (env-gated DELFIN_FFFREE_D8_SQ_ISO, default-OFF).
            # ROOT fix for the D8_SQ_ANGLES isomer collapse: the plain pass pairs donors by geometric
            # most-opposite, which on an AMBIGUOUS tetrahedral ETKDG frame guesses an arbitrary trans and
            # forces cis and trans isomers into the SAME square (FONKOL/WIRFUA 2->1).  Instead: only impose
            # the square when the frame ALREADY shows a CLEAR trans pair (largest L-M-L angle > 135 deg).
            # That clear pair reflects the ISOMER's own arrangement (topology-builder SP-4 / seesaw carries
            # it), so reinforcing it -- plus the second trans pair by elimination -- preserves cis vs trans
            # and flattens seesaw -> square, WITHOUT guessing.  A purely tetrahedral frame (no angle > 135)
            # is left untouched (it is a conformer, not a distinct isomer) so nothing collapses.
            if (not _d8_sq_angles and not suppress_d8_sq
                    and (os.environ.get("DELFIN_FFFREE_D8_SQ_ISO", "0") == "1" or force_d8_sq)
                    and m_sym in _D8_SQ_ISO_METALS and len(donor_indices) == 4):
                try:
                    import numpy as _np
                    _mp = _np.array(coords[m_idx], dtype=float)
                    _dn = {}
                    for _di in donor_indices:
                        _v = _np.array(coords[_di], dtype=float) - _mp
                        _dn[_di] = _v / (float(_np.linalg.norm(_v)) + 1e-12)
                    def _ang(_a, _b):
                        _c = max(-1.0, min(1.0, float(_np.dot(_dn[_a], _dn[_b]))))
                        return math.degrees(math.acos(_c))
                    _d = list(donor_indices)
                    # Global 3-way matching: the ONLY 3 ways to split 4 donors into 2 trans
                    # pairs.  Pick the partition that best fits a square (2 trans ~180, 4 cis
                    # ~90).  This reads the ISOMER's arrangement from the WHOLE frame geometry
                    # (robust) instead of one greedy most-opposite pair (which collapsed cis/
                    # trans) or a single >135 deg pair (which missed strained macrocycles).
                    _parts = [((_d[0], _d[1]), (_d[2], _d[3])),
                              ((_d[0], _d[2]), (_d[1], _d[3])),
                              ((_d[0], _d[3]), (_d[1], _d[2]))]
                    _cand = []
                    for _p in _parts:
                        (_a, _b), (_c, _e) = _p
                        _t1, _t2 = _ang(_a, _b), _ang(_c, _e)
                        _cr = (_ang(_a, _c), _ang(_a, _e), _ang(_b, _c), _ang(_b, _e))
                        _dev = abs(180 - _t1) + abs(180 - _t2) + sum(abs(90 - _x) for _x in _cr)
                        _cand.append((_dev, min(_t1, _t2), _p))
                    _cand.sort(key=lambda t: t[0])
                    _best_dev, _min_trans, _bestp = _cand[0]
                    _runner_dev = _cand[1][0]
                    # Impose the square ONLY when this partition is CLEARLY a square: both its
                    # trans angles are open (>120 deg) AND it clearly beats the runner-up
                    # (margin >30 deg-sum).  A truly tetrahedral frame has no clear winner
                    # (small margin) -> left untouched (it is a conformer, not a distinct
                    # isomer), so nothing collapses.  A strained-but-square macrocycle frame
                    # DOES have a clear winner even at 120-135 deg -> rescued.
                    if _min_trans > 120.0 and (_runner_dev - _best_dev) > 30.0:
                        (_a, _b), (_c, _e) = _bestp
                        base["angles"].append((_a, m_idx, _b, 180.0))
                        base["angles"].append((_c, m_idx, _e, 180.0))
                        _trans_set = {(_a, _b), (_b, _a), (_c, _e), (_e, _c)}
                        for _ii in range(4):
                            for _jj in range(_ii + 1, 4):
                                _x, _y = donor_indices[_ii], donor_indices[_jj]
                                if (_x, _y) not in _trans_set:
                                    base["angles"].append((_x, m_idx, _y, 90.0))
                        _d8_sq_angles = True           # reuse the donor-free UFF path below (no freeze)
                except Exception:
                    pass
            # CN6 OCTAHEDRAL, PERM path (env-gated DELFIN_FFFREE_CN6_OH_ANGLES): the enumerator's EXACT
            # per-isomer trans pairs (perm + _TOPO_TRANS_POSITIONS['OH'], threaded in via d8_trans) impose
            # the octahedron on THIS isomer's own 3 trans axes -- no guessing -> preserves fac/mer/cis/trans
            # (fixes the isomer collapse the greedy fallback caused).  Sets _d8_sq_angles so the greedy
            # geometry fallback below is skipped.
            if (not _d8_sq_angles and d8_trans and not suppress_cn6_oh
                    and (os.environ.get("DELFIN_FFFREE_CN6_OH_ANGLES", "0") == "1" or force_cn6_oh)
                    and len(donor_indices) == 6
                    and _PREFERRED_CN6_GEOMETRY.get(m_sym, 'OH') == 'OH'):
                try:
                    _ds6 = set(donor_indices)
                    _cand6 = [(int(_a), int(_b)) for (_a, _b) in d8_trans
                              if int(_a) in _ds6 and int(_b) in _ds6]
                    _fl6 = [_x for _pr in _cand6 for _x in _pr]
                    if len(_cand6) == 3 and len(set(_fl6)) == 6:
                        _trans_set = set()
                        for (_a, _b) in _cand6:
                            base["angles"].append((_a, m_idx, _b, 180.0))
                            _trans_set |= {(_a, _b), (_b, _a)}
                        for _ii in range(6):
                            for _jj in range(_ii + 1, 6):
                                _x, _y = donor_indices[_ii], donor_indices[_jj]
                                if (_x, _y) not in _trans_set:
                                    base["angles"].append((_x, m_idx, _y, 90.0))
                        _d8_sq_angles = True
                except Exception:
                    pass
            # CN6 OCTAHEDRAL angle targets (env-gated DELFIN_FFFREE_CN6_OH_ANGLES, default-OFF).
            # Analogue of the d8 CN4 SP-4 pass for CN6: ETKDG/UFF can leave an OH-preferring CN6 metal
            # (Mo/W/Re/Ru/... -> _PREFERRED_CN6_GEOMETRY 'OH') as a TRIGONAL PRISM (whole-space measurement:
            # 30 systems TPR-6 built, OC-6 in the crystal -- the biggest single poly_match defect, and
            # UNLIKE d8 the frames are already VALID/topo-correct, just the wrong twist).  Give UFF octahedral
            # 90/180 targets so it relaxes the twist TPR -> OC.  Isomer-safe: the 3 trans pairs are taken
            # from THIS frame's own most-opposite donors (preserves fac/mer -- it's a twist correction, not
            # a donor-arrangement change), and only imposed when all 3 are clearly trans (>120 deg; a valid
            # TPR/OC frame sits at ~140-180, an ambiguous one does not -> skipped, nothing collapses).
            if (not _d8_sq_angles and not suppress_cn6_oh
                    and (os.environ.get("DELFIN_FFFREE_CN6_OH_ANGLES", "0") == "1" or force_cn6_oh)
                    and len(donor_indices) == 6
                    and _PREFERRED_CN6_GEOMETRY.get(m_sym, 'OH') == 'OH'):
                try:
                    import numpy as _np
                    _mp = _np.array(coords[m_idx], dtype=float)
                    _dn = {}
                    for _di in donor_indices:
                        _v = _np.array(coords[_di], dtype=float) - _mp
                        _dn[_di] = _v / (float(_np.linalg.norm(_v)) + 1e-12)
                    _rem = list(donor_indices)
                    _tp6 = []
                    _min_t6 = 999.0
                    while len(_rem) >= 2:               # greedy: pair each donor with its most-opposite
                        _a = _rem[0]
                        _b = min(_rem[1:], key=lambda x: float(_np.dot(_dn[_a], _dn[x])))
                        _cc = max(-1.0, min(1.0, float(_np.dot(_dn[_a], _dn[_b]))))
                        _min_t6 = min(_min_t6, math.degrees(math.acos(_cc)))
                        _tp6.append((_a, _b))
                        _rem.remove(_a); _rem.remove(_b)
                    if len(_tp6) == 3 and _min_t6 > 120.0:
                        _trans_set = set()
                        for (_a, _b) in _tp6:
                            base["angles"].append((_a, m_idx, _b, 180.0))
                            _trans_set |= {(_a, _b), (_b, _a)}
                        for _ii in range(6):
                            for _jj in range(_ii + 1, 6):
                                _x, _y = donor_indices[_ii], donor_indices[_jj]
                                if (_x, _y) not in _trans_set:
                                    base["angles"].append((_x, m_idx, _y, 90.0))
                        _d8_sq_angles = True           # reuse the donor-free UFF path below (no freeze)
                except Exception:
                    pass
            for d_idx in donor_indices:
                if d_idx in chelate_donors:
                    continue
                if _d8_sq_angles:
                    continue                                 # let UFF move donors to the square

                # Baustein-5+6 Phase 3: record monodentate M-D pair (used by
                # opt-in soft-donor UFF mode downstream).
                _soft_meta["donor_indices"].append(d_idx)
                _soft_meta["pairs"].append((m_idx, d_idx))
                if _uff_relax_donors:
                    # Skip freezing; rely on M-D distance constraint to
                    # preserve bond length while UFF relaxes the angle.
                    continue
                if d_idx not in base["fix_atoms"]:
                    base["fix_atoms"].append(d_idx)
                # sp3-C tetrahedral seating companion (DELFIN_FFFREE_SP3C_TET_SEAT): the seating now
                # places an sp3-C donor's substituents tetrahedrally (M-C-X ~109), but UFF with the donor
                # C frozen and its heavy substituent FREE swings that substituent back to linear (there is
                # no M-C-X angle term for the unparameterised metal).  Confirmed via DELFIN_TRACE_SEATING:
                # placement + orient give M-C-X 106, UFF then emits 178.  Freeze the sp3-C's heavy
                # first-shell substituent(s) too so the tetrahedral M-C-X survives UFF.  Byte-identical off.
                _dca = mol_template.GetAtomWithIdx(d_idx)
                if (os.environ.get("DELFIN_FFFREE_SP3C_TET_SEAT", "0") == "1"
                        and _dca.GetSymbol() == "C" and not _dca.GetIsAromatic()
                        and not _dca.IsInRing()):        # PENDANT sp3 alkyl donor only
                    for _nb in _dca.GetNeighbors():
                        if (_nb.GetAtomicNum() > 1 and _nb.GetIdx() != m_idx
                                and _nb.GetIdx() not in base["fix_atoms"]):
                            base["fix_atoms"].append(_nb.GetIdx())

        # Hapto protection: fix hapto metal atoms + all Cp ring atoms
        # during UFF. Without metal FF parameters, UFF treats the metal
        # as empty space, allowing ring atoms to drift out of plane.
        # Fixing these atoms preserves the hapto builder's geometry.
        try:
            hapto_groups = _find_hapto_groups(mol_template)
            fixed = set()
            for _hm, members in hapto_groups:
                fixed.add(_hm)  # fix the hapto metal
                for ci in members:
                    fixed.add(ci)  # fix all ring atoms
            for fi in sorted(fixed):
                if fi not in base["fix_atoms"]:
                    base["fix_atoms"].append(fi)
        except Exception:
            pass

        # Linear M-C≡O carbonyl constraint: M-C-O angle = 180°.
        # UFF without metal parameters bends CO ligands; this pins them.
        try:
            for atom in mol_template.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue
                # Is this a CO carbon? C bonded to metal AND to O via triple bond
                metal_nbr = None
                o_nbr = None
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() in _METAL_SET:
                        metal_nbr = nbr
                    elif nbr.GetAtomicNum() == 8:
                        bond = mol_template.GetBondBetweenAtoms(
                            atom.GetIdx(), nbr.GetIdx()
                        )
                        if bond is not None and bond.GetBondTypeAsDouble() >= 2.5:
                            # Triple bond to O
                            if len(list(nbr.GetNeighbors())) == 1:
                                o_nbr = nbr
                if metal_nbr is not None and o_nbr is not None:
                    m_idx = metal_nbr.GetIdx()
                    c_idx = atom.GetIdx()
                    o_idx = o_nbr.GetIdx()
                    akey = (m_idx, c_idx, o_idx)
                    if akey not in {(a, b, c) for a, b, c, _t in base["angles"]}:
                        base["angles"].append((m_idx, c_idx, o_idx, 180.0))
        except Exception:
            pass

        # Sp2 ring planarity constraints — extended beyond aromatic-flagged.
        # Pyridinium [N+] rings lose aromaticity flag but should stay planar.
        try:
            seen_tors = {
                tuple(sorted([a, b, c, d])): 1
                for a, b, c, d, _t in base["torsions"]
            }
            ring_info = mol_template.GetRingInfo()
            if ring_info is not None:
                for ring in ring_info.AtomRings():
                    if len(ring) < 5 or len(ring) > 7:
                        continue
                    n_sp2 = sum(
                        1 for ri in ring
                        if mol_template.GetAtomWithIdx(ri).GetIsAromatic()
                        or mol_template.GetAtomWithIdx(ri).GetHybridization()
                        == Chem.rdchem.HybridizationType.SP2
                    )
                    if n_sp2 < len(ring) * 0.6:
                        continue
                    if any(
                        mol_template.GetAtomWithIdx(ri).GetSymbol() in _METAL_SET
                        for ri in ring
                    ):
                        continue
                    n = len(ring)
                    for i in range(n):
                        a = ring[(i - 1) % n]
                        b = ring[i]
                        c = ring[(i + 1) % n]
                        d = ring[(i + 2) % n]
                        if len({a, b, c, d}) < 4:
                            continue
                        key = tuple(sorted([a, b, c, d]))
                        if key in seen_tors:
                            continue
                        seen_tors[key] = 1
                        base["torsions"].append((a, b, c, d, 0.0))
        except Exception:
            pass

        # Metallacycle planarity: RDKit's SSSR doesn't include metal-containing
        # rings (M-L bonds aren't typed as ring bonds). For conjugated chelates
        # like acac / salen / ppy, manually detect 5/6-membered metallacycles
        # where all non-metal atoms are sp2, and add planarity torsion
        # constraints that include the metal as the last ring vertex.  UFF
        # then pulls the metal back onto the ring plane — the "metal in
        # π-plane" rule that used to hold implicitly.
        #
        # sp² character is read from the bond graph (any bond of order
        # >= 1.5) rather than from RDKit's hybridisation / aromaticity
        # flags.  This is essential because ``_convert_metal_bonds_to_dative``
        # zeros the bond to the metal and sometimes drops aromatic flags
        # on the cyclometallated donor carbon, which used to silently
        # disable this torsion block for ppy-type ligands.
        def _is_sp2_graph_local(_atom):
            for _b in _atom.GetBonds():
                if (
                    _b.GetBondType() == Chem.BondType.AROMATIC
                    or _b.GetBondType() == Chem.BondType.DOUBLE
                    or _b.GetIsAromatic()
                    or _b.GetBondTypeAsDouble() >= 1.5
                ):
                    return True
            return False

        try:
            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    continue
                m_idx = atom.GetIdx()
                donors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                # For each pair of donors on this metal, find non-metal
                # path between them (chelate backbone).
                for i in range(len(donors)):
                    for j in range(i + 1, len(donors)):
                        d1, d2 = donors[i], donors[j]
                        # BFS from d1 to d2, blocking the metal
                        visited = {m_idx, d1}
                        prev: Dict[int, int] = {d1: -1}
                        queue = [d1]
                        found = False
                        while queue and not found:
                            cur = queue.pop(0)
                            if cur == d2:
                                found = True
                                break
                            for n in mol_template.GetAtomWithIdx(cur).GetNeighbors():
                                ni = n.GetIdx()
                                if ni in visited:
                                    continue
                                visited.add(ni)
                                prev[ni] = cur
                                queue.append(ni)
                                if ni == d2:
                                    found = True
                                    break
                        if not found or d2 not in prev:
                            continue
                        # Reconstruct path d1 → d2
                        path: List[int] = []
                        node = d2
                        while node != -1:
                            path.append(node)
                            node = prev.get(node, -1)
                        path.reverse()
                        # Only 3-5 atom paths (4-6 membered chelate including metal)
                        if len(path) < 3 or len(path) > 5:
                            continue
                        # Bond-graph sp² check — consistent with Rule 4 / 6
                        # in ``_verify_topology_from_graph`` and immune to
                        # sanitisation/dative-bond side effects.
                        all_sp2 = all(
                            _is_sp2_graph_local(mol_template.GetAtomWithIdx(pi))
                            for pi in path
                        )
                        if not all_sp2:
                            continue
                        # Add planarity torsion constraints including metal.
                        # Full cycle: metal - d1 - path... - d2 - metal
                        cycle = [m_idx] + path  # cycle is metal + path d1..d2
                        n = len(cycle)
                        for k in range(n):
                            a = cycle[(k - 1) % n]
                            b = cycle[k]
                            c = cycle[(k + 1) % n]
                            d = cycle[(k + 2) % n]
                            if len({a, b, c, d}) < 4:
                                continue
                            key = tuple(sorted([a, b, c, d]))
                            if key in seen_tors:
                                continue
                            seen_tors[key] = 1
                            base["torsions"].append((a, b, c, d, 0.0))

        except Exception:
            pass

        # Inter-fragment non-bonded repulsion.  Without this, the coordination
        # constraints can pin metal--donor distances while two ligand
        # fragments stay in bond-perception range of each other and UFF is
        # not free to push them apart.  For every pair of heavy atoms from
        # different non-metal fragments whose current distance is below
        # 1.25 * (r_cov_i + r_cov_j), add a distance target at
        # 1.40 * (r_cov_i + r_cov_j); OB's constraint system pushes them to
        # that target on the next minimisation step.  Pairs that are
        # already at or above the safe separation get no constraint, so
        # the rest of the geometry is not touched.
        try:
            non_metal_heavy = [
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() not in _METAL_SET and a.GetAtomicNum() > 1
            ]
            adj_nm: Dict[int, set] = {i: set() for i in non_metal_heavy}
            for bond in mol_template.GetBonds():
                bi = bond.GetBeginAtomIdx()
                bj = bond.GetEndAtomIdx()
                if bi in adj_nm and bj in adj_nm:
                    adj_nm[bi].add(bj)
                    adj_nm[bj].add(bi)
            visited: set = set()
            frag_id: Dict[int, int] = {}
            fid = 0
            for start in sorted(non_metal_heavy):
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
                for i_pos in range(len(non_metal_heavy)):
                    ai = non_metal_heavy[i_pos]
                    fi = frag_id[ai]
                    si = mol_template.GetAtomWithIdx(ai).GetSymbol()
                    ri = _COVALENT_RADII.get(si, 0.75)
                    xi, yi, zi = coords[ai]
                    for j_pos in range(i_pos + 1, len(non_metal_heavy)):
                        aj = non_metal_heavy[j_pos]
                        if frag_id[aj] == fi:
                            continue
                        sj = mol_template.GetAtomWithIdx(aj).GetSymbol()
                        rj = _COVALENT_RADII.get(sj, 0.75)
                        xj, yj, zj = coords[aj]
                        d = math.sqrt(
                            (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                        )
                        thresh = 1.25 * (ri + rj)
                        if d >= thresh:
                            continue
                        key = tuple(sorted((ai, aj)))
                        if key in seen_dist:
                            continue
                        seen_dist.add(key)
                        base["distances"].append(
                            (ai, aj, 1.40 * (ri + rj))
                        )
        except Exception:
            pass

    except Exception as exc:
        logger.debug("_build_coordination_constraints_from_xyz failed: %s", exc)

    # Attach soft-donor meta only when at least one metal is present.
    # Consumer (_optimize_xyz_openbabel) ignores this key unless
    # DELFIN_UFF_SOFT_DONORS=1.  Legacy callers that don't read it are
    # bit-exact unaffected.
    if _soft_meta["metal_indices"]:
        base["_soft_donor_meta"] = _soft_meta

    return base if any(base.values()) else None


def _build_coordination_uff_constraints(
    mol_template,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    xyz_delfin: Optional[str] = None,
) -> Optional[Dict]:
    """UFF constraints that preserve an idealized coordination polyhedron.

    Extends :func:`_build_uff_constraints_from_template` with hard
    metal-donor distance pins (from :func:`_get_ml_bond_length`) and
    donor-metal-donor angle targets derived from
    ``_TOPO_GEOMETRY_VECTORS[geometry]`` combined with ``perm``.

    These constraints keep octahedral/PBP/etc. geometry intact during
    UFF relaxation even when UFF lacks parameters for the metal.
    """
    base = _build_uff_constraints_from_template(
        mol_template, xyz_delfin=xyz_delfin
    )
    if base is None:
        base = {
            "fix_atoms": [],
            "distances": [],
            "angles": [],
            "torsions": [],
        }

    vectors = _TOPO_GEOMETRY_VECTORS.get(geometry)
    if not vectors or metal_idx is None:
        return base if any(base.values()) else None

    try:
        metal_sym = mol_template.GetAtomWithIdx(metal_idx).GetSymbol()
    except Exception:
        return base if any(base.values()) else None

    n_pos = len(vectors)
    if n_pos != len(perm) or n_pos != len(donor_atom_indices):
        return base if any(base.values()) else None

    # Unit vectors for each geometry position.
    unit_vecs: List[Tuple[float, float, float]] = []
    for vx, vy, vz in vectors:
        m = math.sqrt(vx * vx + vy * vy + vz * vz)
        if m < 1e-8:
            unit_vecs.append((0.0, 0.0, 0.0))
        else:
            unit_vecs.append((vx / m, vy / m, vz / m))

    # Track existing distance keys so we don't duplicate OCO pins.
    seen_dist = {tuple(sorted((a, b))) for a, b, _t in base["distances"]}

    # M-D distance pins.
    donor_pos_map: Dict[int, int] = {}
    for pos_idx, donor_list_idx in enumerate(perm):
        if donor_list_idx < 0 or donor_list_idx >= len(donor_atom_indices):
            continue
        donor_atom_idx = donor_atom_indices[donor_list_idx]
        donor_pos_map[donor_atom_idx] = pos_idx
        try:
            donor_sym = mol_template.GetAtomWithIdx(donor_atom_idx).GetSymbol()
        except Exception:
            continue
        bl = float(_get_ml_bond_length(metal_sym, donor_sym))
        key = tuple(sorted((metal_idx, donor_atom_idx)))
        if key in seen_dist:
            continue
        seen_dist.add(key)
        base["distances"].append((metal_idx, donor_atom_idx, bl))

    # L-M-L angle targets from ideal geometry vectors.
    seen_angle = {
        (a, b, c) for a, b, c, _t in base["angles"]
    }
    seen_angle.update(
        (c, b, a) for a, b, c, _t in base["angles"]
    )
    donors_with_pos = [
        d for d in donor_atom_indices if d in donor_pos_map
    ]
    for i in range(len(donors_with_pos)):
        for j in range(i + 1, len(donors_with_pos)):
            d_i = donors_with_pos[i]
            d_j = donors_with_pos[j]
            pi = donor_pos_map[d_i]
            pj = donor_pos_map[d_j]
            vi = unit_vecs[pi]
            vj = unit_vecs[pj]
            dot = vi[0] * vj[0] + vi[1] * vj[1] + vi[2] * vj[2]
            dot = max(-1.0, min(1.0, dot))
            ideal = math.degrees(math.acos(dot))
            key = (d_i, metal_idx, d_j)
            if key in seen_angle:
                continue
            seen_angle.add(key)
            seen_angle.add((d_j, metal_idx, d_i))
            base["angles"].append((d_i, metal_idx, d_j, ideal))

    if not (
        base["fix_atoms"]
        or base["distances"]
        or base["angles"]
        or base["torsions"]
    ):
        return None
    return base


# ---------------------------------------------------------------------------
# Bondi van-der-Waals radii (Å) — identical values to the inter-ligand-clash
# detector (find_inter_ligand_clash.py).  Used ONLY by the deterministic
# geometric clash-relief pass below so the floor it enforces matches the
# metric it is meant to improve.  Default 1.70 for anything missing.
# ---------------------------------------------------------------------------
_VDW_RADII_CLASH: Dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.11, "Ti": 1.95, "V": 1.92, "Cr": 1.89,
    "Mn": 1.97, "Fe": 1.94, "Co": 1.92, "Ni": 1.84, "Cu": 1.86, "Zn": 2.10,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.83, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18, "Mo": 2.17,
    "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43, "Ce": 2.42, "Pr": 2.40, "Nd": 2.39,
    "Pm": 2.38, "Sm": 2.36, "Eu": 2.35, "Gd": 2.34, "Tb": 2.33, "Dy": 2.31,
    "Ho": 2.30, "Er": 2.29, "Tm": 2.27, "Yb": 2.26, "Lu": 2.24, "Hf": 2.23,
    "Ta": 2.22, "W": 2.18, "Re": 2.16, "Os": 2.16, "Ir": 2.13, "Pt": 2.13,
    "Au": 2.14, "Hg": 2.23, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
}
