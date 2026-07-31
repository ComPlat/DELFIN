"""Analytic energy terms + gradients for Baustein 6 variational refinement.

Eight U-terms compose the total energy functional minimised by L-BFGS-B
in :mod:`delfin.manta._variational_refiner`:

* ``U_bond``      — bond-length penalty (covalent radii lookup)
* ``U_angle``     — bond-angle penalty per hybridization
* ``U_clash``     — one-sided vdW overlap (Bondi/Alvarez radii)
* ``U_topology``  — log-barrier on M-D bonds (HARD topology gate)
* ``U_A``         — Tier A coord-sphere pull to Hungarian-assigned slots
* ``U_B``         — Tier B Morgan equivalence (bond/angle means)
* ``U_C``         — Tier C per-fragment archetype point-group enforcement
* ``U_D``         — Tier D global molecular point-group enforcement

Each function returns ``(energy_value, gradient_array)`` where ``gradient_array``
has the same shape as ``coords`` (``[N, 3]``). Gradients are fully analytic
(no finite differencing) so L-BFGS-B can run with ``jac=True``.

All dependencies on :mod:`delfin.smiles_converter` (``_get_ml_bond_length``,
``_COVALENT_RADII``, ``_METAL_SET``) and :mod:`delfin.manta._vdw_radii` are lazily
imported inside helpers to keep the module light at import time and to avoid
the known circular-import risk with ``smiles_converter``.

Reference: ``iters/BAUSTEIN6_MASTERPLAN.md`` Section 3.3.
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Numerical constants
# ---------------------------------------------------------------------------

_EPS_COS = 1.0e-7      # clip for cos(theta) to avoid arccos NaN
_EPS_SIN = 1.0e-6      # treat angle as collinear / undefined below this
_EPS_DIST = 1.0e-8     # treat bond length below this as atom collapse
_EPS_BARRIER = 1.0e-4  # safety margin inside log-barrier domain

# Ideal angles (radians) per hybridization.
_THETA_SP = math.pi                # 180.0 deg
_THETA_SP2 = math.radians(120.0)   # 120.0 deg
_THETA_SP3 = math.radians(109.47)  # tetrahedral

# Covalent-bond-order scaling for ideal bond length d_ideal = scale * (r_i + r_j).
# Pyykkö-style; AROMATIC sits between SINGLE and DOUBLE.
#
# ⚠️ This is a GLOBAL factor applied to every element pair, and measurement says that is the defect.
# Pyykkö & Atsumi published single, DOUBLE and TRIPLE covalent radii per ELEMENT; collapsing all of
# them into one ratio per bond order is a simplification the data does not support.  Measured against
# CCDC clean_v2 (25.8 M bonds, signature-resolved bands): terminal multiply-bonded heteroatoms on
# high-coordination centres come out at SINGLE-bond length -- ClO4- Cl-O +0.55 A (100 % outside the
# band), N=N=N +0.295, S=O +0.20, coordinated C#N +0.18 across 472 systems, C#S +0.17 -- while P=O
# (+0.055) and N-O (-0.025) are fine.  A single ratio cannot be right for both groups at once.
# Kept as the FALLBACK for elements Pyykkö has no order-specific radius for; _pyykko_order_radius
# below is consulted first.
_BOND_ORDER_SCALE: Dict[float, float] = {
    1.0: 1.00,
    1.5: 0.93,
    2.0: 0.87,
    3.0: 0.78,
}

# Order-specific Pyykkö radii, looked up per element and cached, so a build only ever queries the
# handful of elements it contains.  DELFIN already uses the SINGLE set elsewhere
# (occupier.load_covalent_radii, source "pyykko2009"); the double/triple sets were published in the
# same work.  This turns the bond-order contraction from a global ratio into a property of the ATOM,
# defined across the periodic table instead of for a hand-listed set of pairs.
_PYYKKO_ORDER_ATTR: Dict[int, str] = {2: "covalent_radius_pyykko_double",
                                      3: "covalent_radius_pyykko_triple"}
_PYYKKO_ORDER_CACHE: Dict[Tuple[str, int], Optional[float]] = {}


def _pyykko_order_radius(sym: str, order: int) -> Optional[float]:
    """Pyykkö double/triple covalent radius (Å) for `sym`, or None when unavailable.

    None makes every caller fall back to the global-ratio path, so an element mendeleev has no
    order-specific value for behaves exactly as before -- never-worse by construction.
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
                val = float(raw) / 100.0          # mendeleev stores pm
        except Exception:
            val = None
    _PYYKKO_ORDER_CACHE[key] = val
    return val

_LAST_RESORT_BOND_LEN = 1.50   # Å — used if nothing else is known
_DEFAULT_VDW_FALLBACK = 1.80   # Å — also used for unknown elements


# ---------------------------------------------------------------------------
# Lazy-import helpers (avoid module-import cycles)
# ---------------------------------------------------------------------------

def _smiles_converter_metals():
    """Return ``_METAL_SET`` from :mod:`delfin.smiles_converter`."""
    try:
        from delfin.smiles_converter import _METAL_SET  # type: ignore
        return _METAL_SET
    except Exception:
        return frozenset()


def _smiles_converter_covalent_radii():
    """Return ``_COVALENT_RADII`` from :mod:`delfin.smiles_converter`."""
    try:
        from delfin.smiles_converter import _COVALENT_RADII  # type: ignore
        return _COVALENT_RADII
    except Exception:
        return {}


def _smiles_converter_ml_bondlen(metal_sym: str, donor_sym: str) -> float:
    """Return tabulated M-L bond length (Å), with safe fallback to 2.0 Å."""
    try:
        from delfin.smiles_converter import _get_ml_bond_length  # type: ignore
        return float(_get_ml_bond_length(metal_sym, donor_sym))
    except Exception:
        return 2.0


def _vdw_radius(symbol: str) -> float:
    """Return Bondi/Alvarez vdW radius for ``symbol``."""
    try:
        from delfin.manta._vdw_radii import get_vdw_radius  # type: ignore
        return float(get_vdw_radius(symbol))
    except Exception:
        return _DEFAULT_VDW_FALLBACK


# ---------------------------------------------------------------------------
# Internal mol introspection
# ---------------------------------------------------------------------------

def _symbols(mol) -> List[str]:
    """Element symbol per atom index."""
    return [a.GetSymbol() for a in mol.GetAtoms()]


def _is_metal_atom(mol, idx: int) -> bool:
    sym = mol.GetAtomWithIdx(idx).GetSymbol()
    return sym in _smiles_converter_metals()


def _enumerate_bonds(mol) -> List[Tuple[int, int, float]]:
    """Return list of ``(i, j, bond_order_double)`` for every bond in ``mol``."""
    out: List[Tuple[int, int, float]] = []
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        try:
            order = float(b.GetBondTypeAsDouble())
        except Exception:
            order = 1.0
        out.append((int(i), int(j), order))
    return out


def _ideal_bond_length(mol, i: int, j: int, order: float) -> float:
    """Return ideal bond length d_ideal for the bond ``(i, j)`` of order ``order``.

    Lookup order:
      1. M-L bond (one endpoint is a metal) → ``_get_ml_bond_length``.
      2. Otherwise covalent-radii sum × bond-order scale.
      3. Fallback ``_LAST_RESORT_BOND_LEN`` if radii unknown.
    """
    metals = _smiles_converter_metals()
    sym_i = mol.GetAtomWithIdx(i).GetSymbol()
    sym_j = mol.GetAtomWithIdx(j).GetSymbol()
    if sym_i in metals and sym_j not in metals:
        return _smiles_converter_ml_bondlen(sym_i, sym_j)
    if sym_j in metals and sym_i not in metals:
        return _smiles_converter_ml_bondlen(sym_j, sym_i)
    cov = _smiles_converter_covalent_radii()
    r_i = cov.get(sym_i)
    r_j = cov.get(sym_j)
    if r_i is None or r_j is None:
        return _LAST_RESORT_BOND_LEN
    # Order-specific radii first (a property of the ATOM), global ratio as fallback.
    _o = 3 if order >= 2.5 else (2 if order >= 1.75 else 0)
    if _o and os.environ.get("DELFIN_FFREE_PYYKKO_ORDER_RADII", "0") == "1":
        _ri = _pyykko_order_radius(sym_i, _o)
        _rj = _pyykko_order_radius(sym_j, _o)
        if _ri is not None and _rj is not None:
            return float(_ri + _rj)
    scale = _BOND_ORDER_SCALE.get(round(order * 2) / 2.0)  # snap to .5
    if scale is None:
        scale = 1.0
    return float(r_i + r_j) * scale


# --------------------------------------------------------------------------
# Signature-resolved ideal angle (measured, atom- and bond-specific)
# --------------------------------------------------------------------------
# WHY.  ``_hybridization_theta`` below reads RDKit's hybridization LABEL, which offers
# three numbers for the whole periodic table and is decided WITHOUT the metal.  Measured
# against the CSD that is wrong in two ways that matter:
#   * the X-D-X angle at a saturated donor follows the atom's PERIOD, not one 109.47 for
#     every element (2:109.5  3:104.0  4:101.5  5:100.5  6:99.0);
#   * a donor whose lone pair joins a neighbouring pi system is TRIGONAL PLANAR with the
#     metal IN that plane -- RDKit calls the same nitrogen sp3 and asks for 109.5.  That
#     second case is AFOMEO frame 4: R2N-M pulled out of the pi plane.
# Everything below keys off the ATOM and its BONDS (element, period, sigma-partner count
# WITH the metal counted, pi activity, ring membership and ring size).  No functional
# group is named anywhere -- "not every nitrogen is the same" is expressed as the bonding
# situation at that one nitrogen.
_THETA_BY_PERIOD: Dict[int, float] = {2: 109.5, 3: 104.0, 4: 101.5, 5: 100.5, 6: 99.0}
# Ring angles: crystal-measured where we have them, else the regular polygon
# (internal = 180 - 360/n, exocyclic = (360 - internal)/2).
_RING_INTERNAL_MEAS: Dict[int, float] = {5: 106.0, 6: 118.0}
_RING_EXOCYCLIC_MEAS: Dict[int, float] = {5: 126.5, 6: 120.75}
# M-D-X at a SATURATED donor: as the two substituents close down, the metal is pushed to a
# wider angle.  Measured slope, not fitted here.
_PHI_MD_SLOPE = 0.87
_THETA_PLANAR = 120.0
_THETA_LINEAR = 180.0


def _period_of(z: int) -> int:
    for p, hi in ((1, 2), (2, 10), (3, 18), (4, 36), (5, 54), (6, 86)):
        if z <= hi:
            return p
    return 7


def _ring_sizes_containing(mol, *idxs) -> List[int]:
    """Sizes of the rings that contain EVERY given atom index (smallest first)."""
    try:
        rings = mol.GetRingInfo().AtomRings()
    except Exception:
        return []
    out = [len(r) for r in rings if all(x in r for x in idxs)]
    return sorted(out)


def _is_pi_active(mol, idx: int) -> bool:
    """True if atom ``idx`` itself carries pi character (aromatic or a multiple bond)."""
    try:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic():
            return True
        for b in atom.GetBonds():
            if float(b.GetBondTypeAsDouble()) >= 1.5:
                return True
    except Exception:
        pass
    return False


def _has_conjugated_lone_pair(mol, idx: int) -> bool:
    """True if ``idx`` has a lone pair AND a neighbour with pi character.

    The lone pair then delocalises into that pi system and the centre flattens.  Decided
    from the element's group and the atom's own charge -- never from a group name.
    """
    try:
        atom = mol.GetAtomWithIdx(idx)
        z = atom.GetAtomicNum()
    except Exception:
        return False
    # Groups 15/16/17 across every period carry lone pairs; a negative formal charge adds
    # one to anything (a carbanion donor included), so no element is special-cased.
    has_lp = (z in (7, 8, 9, 15, 16, 17, 33, 34, 35, 51, 52, 53, 83, 84, 85)
              or atom.GetFormalCharge() < 0)
    if not has_lp:
        return False
    for nb in atom.GetNeighbors():
        if _is_pi_active(mol, int(nb.GetIdx())):
            return True
    return False


def _signature_theta(mol, i: int, j: int, k: int) -> float:
    """Ideal angle (radians) for the triple (i, j, k) from the bonding signature at ``j``.

    Order of the rules is chemistry, not convenience:
      1. ring-internal / ring-exocyclic -- a ring fixes the angle geometrically, and for a
         planar ring donor the metal simply IS the exocyclic substituent (pyridine N->M
         lands on ~121 deg, which the period rule would have got wrong);
      2. trigonal-planar centre -- 3 sigma partners (metal counted) plus pi character or a
         conjugated lone pair.  Applies to EVERY triple at ``j`` including M-D-X, which is
         precisely "metal in the pi plane";
      3. linear centre -- 2 sigma partners with a triple bond or two doubles;
      4. M-D-X at a saturated donor -- period rule through the measured slope;
      5. X-D-X at a saturated donor -- period rule.
    """
    metals = _smiles_converter_metals()
    try:
        a_j = mol.GetAtomWithIdx(j)
        sym_i = mol.GetAtomWithIdx(i).GetSymbol()
        sym_k = mol.GetAtomWithIdx(k).GetSymbol()
    except Exception:
        return _THETA_SP3
    i_is_m = sym_i in metals
    k_is_m = sym_k in metals

    # --- 1. rings ---------------------------------------------------------------
    # POLYGON ONLY UP TO SIX (DELFIN_FFREE_POLY6=1, default OFF -> byte-identical).
    #
    # 180 - 360/n is the internal angle of a REGULAR PLANAR polygon.  A seven- or
    # nine-membered ring is neither: it puckers, and its sides are not equal.  Measured
    # against the CSD angle table, the formula over-predicts by +16.3 deg at ring 7,
    # +11.5 at ring 8 and +10.0 above eight, while at ring 6 -- which carries 3.06 M of
    # the 4.67 M measured angles -- it is right to 5 deg.  A CH2 in a macrocycle is set by
    # its own hybridisation, not by the ring's perimeter.
    #
    # Judged alone against 1900 measured bins: 180 bins improve, 23 worsen; RMS 15.51 ->
    # 13.70, median |residual| 3.62 -> 2.88, bins inside their own sigma 39.6 -> 44.6 %.
    #
    # NOT changed, and measured as the wrong idea: letting a metal in the triple outrank
    # the ring rule.  That looked compelling -- M-D-X is the best branch in the law, 86.8 %
    # of its bins inside sigma, and the ring branch blocks it for every chelate donor --
    # but it scores 135 bins better against 188 worse.  A donor inside a chelate ring obeys
    # RING CLOSURE, not the period rule; the ring branch is right there and merely has the
    # wrong closure (a polygon instead of the actual geometry).  That case cannot be fixed
    # in a signature table at all -- it needs the 1,3 distances, i.e. the seating.
    _poly_max = 6 if os.environ.get("DELFIN_FFREE_POLY6", "0") == "1" else 10 ** 9
    common = _ring_sizes_containing(mol, i, j, k)
    if common and common[0] <= _poly_max:
        n = common[0]
        return math.radians(_RING_INTERNAL_MEAS.get(n, 180.0 - 360.0 / n))
    j_rings = _ring_sizes_containing(mol, j)
    if j_rings and j_rings[0] <= _poly_max:
        # Exocyclic only if exactly one of i/k shares a ring with j.
        in_ring = sum(1 for x in (i, k) if _ring_sizes_containing(mol, j, x))
        if in_ring == 1:
            n = j_rings[0]
            internal = _RING_INTERNAL_MEAS.get(n, 180.0 - 360.0 / n)
            return math.radians(
                _RING_EXOCYCLIC_MEAS.get(n, (360.0 - internal) / 2.0))

    # --- sigma partners, WITH the metal counted --------------------------------
    try:
        n_sigma = int(a_j.GetDegree())
        z_j = int(a_j.GetAtomicNum())
        max_order = max((float(b.GetBondTypeAsDouble()) for b in a_j.GetBonds()),
                        default=1.0)
        n_double = sum(1 for b in a_j.GetBonds()
                       if 1.75 <= float(b.GetBondTypeAsDouble()) < 2.5)
    except Exception:
        return _THETA_SP3

    # --- 2. trigonal planar (this is the AFOMEO case) --------------------------
    if n_sigma == 3 and (_is_pi_active(mol, j) or _has_conjugated_lone_pair(mol, j)):
        return math.radians(_THETA_PLANAR)

    # --- 3. linear -------------------------------------------------------------
    if n_sigma == 2 and (max_order >= 2.5 or n_double >= 2):
        return math.radians(_THETA_LINEAR)

    theta_base = _THETA_BY_PERIOD.get(_period_of(z_j), 109.5)

    # --- 4. M-D-X at a saturated donor ----------------------------------------
    if i_is_m != k_is_m:
        return math.radians(109.5 - _PHI_MD_SLOPE * (theta_base - 109.5))

    # --- 5. X-D-X at a saturated donor ----------------------------------------
    return math.radians(theta_base)


def _build_angle_targets(mol) -> List[Tuple[int, int, int, float]]:
    """``(i, j, k, theta_ideal)`` for every angle, resolved once per refinement.

    The caller memoises this in the per-refine ``sym_info`` dict, so the signature work
    is paid once instead of once per L-BFGS iteration (a 200-iteration run over a few
    thousand angles would otherwise dominate the runtime).
    """
    out = []
    for (i, j, k) in _enumerate_angles(mol):
        mb = _measured_angle_band(mol, i, j, k)
        if mb is not None:
            out.append((i, j, k, math.radians(mb[1])))     # measured p50 beats the rule
        else:
            out.append((i, j, k, _signature_theta(mol, i, j, k)))
    return out


def _hybridization_theta(mol, j: int) -> float:
    """Return ideal bond angle (radians) at atom ``j`` from RDKit hybridization."""
    try:
        from rdkit.Chem.rdchem import HybridizationType  # type: ignore
    except Exception:
        return _THETA_SP3  # safe default
    try:
        hyb = mol.GetAtomWithIdx(j).GetHybridization()
    except Exception:
        return _THETA_SP3
    if hyb == HybridizationType.SP:
        return _THETA_SP
    if hyb in (HybridizationType.SP2, HybridizationType.SP2D):  # rare alias
        return _THETA_SP2
    if hyb in (HybridizationType.SP3, HybridizationType.SP3D, HybridizationType.SP3D2):
        return _THETA_SP3
    # Unknown / unspecified → tetrahedral default.
    return _THETA_SP3


def _enumerate_angles(mol) -> List[Tuple[int, int, int]]:
    """Return all angle triples (i, j, k) with j the central atom.

    Iterates RDKit neighbour lists, skipping metal-centred angles (those are
    governed by ``U_A`` polyhedron pull rather than hybridization defaults).
    """
    metals = _smiles_converter_metals()
    triples: List[Tuple[int, int, int]] = []
    for atom in mol.GetAtoms():
        j = int(atom.GetIdx())
        if atom.GetSymbol() in metals:
            continue
        nbrs = [int(n.GetIdx()) for n in atom.GetNeighbors()]
        n = len(nbrs)
        if n < 2:
            continue
        for a in range(n):
            for b in range(a + 1, n):
                triples.append((nbrs[a], j, nbrs[b]))
    return triples


def _bonded_pairs_set(mol) -> set:
    """Return ``{(min(i,j), max(i,j))}`` for every covalent bond."""
    s = set()
    for b in mol.GetBonds():
        i = int(b.GetBeginAtomIdx())
        j = int(b.GetEndAtomIdx())
        if i > j:
            i, j = j, i
        s.add((i, j))
    return s


def _1_3_pairs_set(mol) -> set:
    """Return 1-3 pairs (atoms sharing a common neighbour). Excluded from clash."""
    s = set()
    for atom in mol.GetAtoms():
        nbrs = [int(n.GetIdx()) for n in atom.GetNeighbors()]
        m = len(nbrs)
        for a in range(m):
            for b in range(a + 1, m):
                i, j = nbrs[a], nbrs[b]
                if i > j:
                    i, j = j, i
                s.add((i, j))
    return s


def _enumerate_metal_donor_bonds(mol) -> List[Tuple[int, int]]:
    """Return list ``[(metal_idx, donor_idx)]`` for every M-L covalent bond."""
    metals = _smiles_converter_metals()
    pairs: List[Tuple[int, int]] = []
    for b in mol.GetBonds():
        i = int(b.GetBeginAtomIdx())
        j = int(b.GetEndAtomIdx())
        si = mol.GetAtomWithIdx(i).GetSymbol()
        sj = mol.GetAtomWithIdx(j).GetSymbol()
        if si in metals and sj not in metals:
            pairs.append((i, j))
        elif sj in metals and si not in metals:
            pairs.append((j, i))
    return pairs


# ---------------------------------------------------------------------------
# U_bond
# ---------------------------------------------------------------------------

_FLATBOTTOM_FRAC = 0.05   # half-width of the tolerated band, as a fraction of d_ideal

# Normal-distribution conversion: p90 - p10 = 2 * 1.2816 * sigma.  Kept next to the band
# readers so the gate in _variational_refiner and U_bond cannot drift apart on it.
_P10_P90_TO_SIGMA = 2.5631
_INVVAR_W_MAX = 25.0      # weight bound, i.e. sigma ratio 5 in either direction


# ---------------------------------------------------------------------------
# Bond signature — the SHARED contract between the functional and its calibration
# ---------------------------------------------------------------------------
# The measured band table is keyed by this function, and the functional looks its target
# up with this same function.  It therefore lives HERE, in the functional, and the
# offline measurement imports it -- so the key can never drift between the two sides.
# That is the whole reason not to port the eye's signature instead: two copies of a key
# function are two things that can disagree.
#
# Every field is a property of the ATOM and its BONDS -- element, degree with the metal
# counted, pi character, conjugated lone pair, smallest ring -- and of the bond itself.
# No functional group is named, and the same key carries N, S, Se, Sb, U without a branch.

def _atom_sig(mol, idx: int) -> str:
    a = mol.GetAtomWithIdx(idx)
    rings = _ring_sizes_containing(mol, idx)
    return "{},{},{},{},{}".format(
        a.GetSymbol(),
        int(a.GetDegree()),                       # metal counted: it is a neighbour
        1 if _is_pi_active(mol, idx) else 0,
        1 if _has_conjugated_lone_pair(mol, idx) else 0,
        rings[0] if rings else 0,
    )


def bond_signature_keys(mol, i: int, j: int) -> List[str]:
    """Keys for the bond (i, j), most specific first.

    The ladder is the fallback order for a sparse bin: drop the conjugated-lone-pair
    flag, then the ring size, then pi, then the degree.  A rare heavy-donor bond still
    lands on SOME measured band instead of falling back to a radii sum.
    """
    try:
        b = mol.GetBondBetweenAtoms(i, j)
        order = round(float(b.GetBondTypeAsDouble()) * 2) / 2.0
        bring = 1 if b.IsInRing() else 0
    except Exception:
        order, bring = 1.0, 0
    sa, sb = sorted((_atom_sig(mol, i), _atom_sig(mol, j)))
    fa, fb = sa.split(","), sb.split(",")
    el = sorted((fa[0], fb[0]))
    out = [
        f"{sa}|{sb}|{order}|{bring}",                                    # L0 full
        f"{','.join(fa[:4])}|{','.join(fb[:4])}|{order}|{bring}",         # L1 no ring size
        f"{','.join(fa[:3])}|{','.join(fb[:3])}|{order}|{bring}",         # L2 no lone pair
        f"{','.join(fa[:2])}|{','.join(fb[:2])}|{order}",                 # L3 element+degree
        f"{el[0]}|{el[1]}|{order}",                                       # L4 element pair
    ]
    return out


def angle_signature_keys(mol, i: int, j: int, k: int) -> List[str]:
    """Keys for the angle (i, j, k), j central, most specific first.

    Same contract as :func:`bond_signature_keys`: the offline measurement imports this,
    so key and lookup cannot drift.  The CENTRAL atom's signature leads -- it is what
    decides an angle -- and the two outer signatures are sorted so the key is symmetric.

    Measured need: rebuilding crystals with rule-based theta (period / pi / ring polygon)
    costs 0.20 A RMSD, the largest single contribution of any target class, while measured
    bond lengths cost 0.018 A.  The rule was a good first approximation; a measured p50
    per bin should beat it exactly as it did for lengths.
    """
    try:
        b_ij = mol.GetBondBetweenAtoms(i, j)
        b_jk = mol.GetBondBetweenAtoms(j, k)
        o_ij = round(float(b_ij.GetBondTypeAsDouble()) * 2) / 2.0
        o_jk = round(float(b_jk.GetBondTypeAsDouble()) * 2) / 2.0
    except Exception:
        o_ij = o_jk = 1.0
    ring = 1 if _ring_sizes_containing(mol, i, j, k) else 0
    sj = _atom_sig(mol, j)
    outer = sorted((f"{_atom_sig(mol, i)}~{o_ij}", f"{_atom_sig(mol, k)}~{o_jk}"))
    fj = sj.split(",")
    eo = sorted((mol.GetAtomWithIdx(i).GetSymbol(),
                 mol.GetAtomWithIdx(k).GetSymbol()))
    return [
        f"{sj}#{outer[0]}#{outer[1]}#{ring}",                             # L0 full
        f"{','.join(fj[:4])}#{outer[0]}#{outer[1]}#{ring}",                # L1 no ring size
        f"{sj}#{eo[0]}#{eo[1]}#{ring}",                                    # L2 outer = element
        f"{','.join(fj[:3])}#{eo[0]}#{eo[1]}#{ring}",                      # L3 + no lone pair
        f"{','.join(fj[:2])}#{eo[0]}#{eo[1]}",                             # L4 element+degree
    ]


def _evidence_keep(n_str: str, lo: float, mid: float, hi: float) -> bool:
    """May a bin measured on ``n`` structures with this width steer a build?

    ONE definition for every table, imported from polyhedra so the bond/angle readers here
    and the M-D / D-M-D reader there can never drift apart on what counts as evidence --
    that drift is exactly how the round-trip deactivation ended up sitting on one of three
    call sites.  Falls OPEN (keeps everything) if the helpers cannot be imported, so an
    import problem can never silently narrow the tables and be mistaken for a measurement.
    """
    try:
        from delfin.manta.polyhedra import (  # type: ignore
            _EVIDENCE_ON, _EVIDENCE_MIN_N, _EVIDENCE_MAX_REL,
        )
    except Exception:
        return True
    if not _EVIDENCE_ON():
        return True
    try:
        if int(float(n_str)) < _EVIDENCE_MIN_N():
            return False
    except ValueError:
        return True
    rel = (hi - lo) / abs(mid) if mid else 9.9
    return rel <= _EVIDENCE_MAX_REL()


_BOND_BAND_CACHE: Dict[str, Dict[str, Tuple[float, float, float]]] = {}
_ANGLE_BAND_CACHE: Dict[str, Dict[str, Tuple[float, float, float]]] = {}


def _angle_band_table() -> Dict[str, Tuple[float, float, float]]:
    """``key -> (p10, p50, p90)`` in DEGREES from ``DELFIN_FFREE_ANGLE_BANDS``.

    Same runtime-path split as the bond table: the measurement is CSD-derived and does
    not ship.  Unset -> empty -> every caller keeps the rule-based theta.
    """
    path = os.environ.get("DELFIN_FFREE_ANGLE_BANDS", "")
    if not path:
        return {}
    tbl = _ANGLE_BAND_CACHE.get(path)
    if tbl is not None:
        return tbl
    tbl = {}
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) >= 6:
                    try:
                        _lo, _md, _hi = float(p[3]), float(p[4]), float(p[5])
                        if not _evidence_keep(p[2], _lo, _md, _hi):
                            continue
                        tbl[p[1]] = (_lo, _md, _hi)
                    except ValueError:
                        continue
    except Exception:
        tbl = {}
    _ANGLE_BAND_CACHE[path] = tbl
    return tbl


def _measured_angle_band(mol, i: int, j: int, k: int):
    """``(p10, p50, p90)`` in degrees for this angle, walking the key ladder, or None."""
    tbl = _angle_band_table()
    if not tbl:
        return None
    for key in angle_signature_keys(mol, i, j, k):
        v = tbl.get(key)
        if v is not None:
            return v
    return None


def _bond_band_table() -> Dict[str, Tuple[float, float, float]]:
    """``key -> (p10, p50, p90)`` from the table ``DELFIN_FFREE_BOND_BANDS`` points at.

    The measurement is CSD-derived and therefore does NOT ship here; the path is supplied
    at runtime, the same split already used for the torsion L1 table.  Empty dict when
    unset, which makes every caller fall back to the radii-sum reference -- unchanged
    behaviour, never a crash.
    """
    path = os.environ.get("DELFIN_FFREE_BOND_BANDS", "")
    if not path:
        return {}
    tbl = _BOND_BAND_CACHE.get(path)
    if tbl is not None:
        return tbl
    tbl = {}
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                p = ln.rstrip("\n").split("\t")
                # level, key, n, p10, p50, p90
                if len(p) >= 6:
                    try:
                        _lo, _md, _hi = float(p[3]), float(p[4]), float(p[5])
                        if not _evidence_keep(p[2], _lo, _md, _hi):
                            continue
                        tbl[p[1]] = (_lo, _md, _hi)
                    except ValueError:
                        continue
    except Exception:
        tbl = {}
    _BOND_BAND_CACHE[path] = tbl
    return tbl


def _measured_bond_band(mol, i: int, j: int):
    """``(p10, p50, p90)`` for this bond from the measured table, walking the key ladder
    from most specific to least, or None when the table has nothing for it."""
    tbl = _bond_band_table()
    if not tbl:
        return None
    for k in bond_signature_keys(mol, i, j):
        v = tbl.get(k)
        if v is not None:
            return v
    return None


def _flat_bottom(d: float, target: float, half_width: float) -> Tuple[float, float]:
    """``(penalty, d(penalty)/dd)`` — zero inside ``target ± half_width``, quadratic
    outside.  C1 at the join.

    A harmonic term ``(d - d_ideal)^2`` has ONE zero-force point, so it is a restraint to
    a single length, not a statement about what is realistic.  Measured consequence: at
    107 real crystals U_bond pulls HARDER than on our own frames (ratio 0.79) -- our
    builder places bonds at exactly d_ideal, so we sit in the term's minimum while reality,
    which is a distribution, does not.  A flat bottom makes every length inside the
    tolerated band cost nothing, so a crystal is stationary by construction and only
    genuinely out-of-band geometry is pulled.
    """
    lo = target - half_width
    hi = target + half_width
    if d < lo:
        e = lo - d
        return e * e, -2.0 * e
    if d > hi:
        e = d - hi
        return e * e, 2.0 * e
    return 0.0, 0.0


def U_bond(coords: np.ndarray, mol, k_bond: float = 1000.0
           ) -> Tuple[float, np.ndarray]:
    """Bond-length penalty: Σ k_bond · (d - d_ideal)²  over covalent bonds.

    With ``DELFIN_FFREE_BOND_BAND=1`` the point restraint becomes a flat-bottomed band
    (see :func:`_flat_bottom`), which is what the crystal-stationarity measurement asks
    for.  Default OFF -> byte-identical.

    Skips M-L bonds (those are handled by ``U_topology`` and ``U_A``).
    """
    _band = os.environ.get("DELFIN_FFREE_BOND_BAND", "0") == "1"
    _invvar = os.environ.get("DELFIN_FFREE_INVVAR", "0") == "1"
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    metals = _smiles_converter_metals()
    if k_bond <= 0.0 or n == 0:
        return float(energy), grad

    # INVERSE-VARIANCE WEIGHTING (2026-07-30, default OFF, needs _band).
    #
    # k_bond is one global constant for every bond, which asserts that a C-C single bond
    # and a Ru-N bond are known equally well.  The measurement says otherwise: bin widths
    # differ by more than an order of magnitude.  Weighting by 1/sigma^2 states each
    # restraint in units of how well that bond is actually known.
    #
    # NORMALISED, not absolute.  sigma ~ 0.02 A would multiply every energy by ~2500 and
    # silently drown the other seven terms -- the weighting is meant to change the balance
    # WITHIN U_bond, not the balance BETWEEN terms.  Dividing by the molecule's own median
    # sigma keeps the total scale (and the meaning of k_bond) intact.
    #
    # And note what it couples to: the WIDTH OF THE BAND, a property of chemistry.  Never
    # the SIZE OF THE ERROR.  Weighting by how broken a site is would hand a torn ligand
    # enormous forces and tear it further -- the functional would be most destructive
    # exactly on the frames it exists to rescue.  Broken sites already pull harder without
    # that: the excursion enters quadratically, so 5 sigma is 100x the force of 0.5 sigma.
    _sig_ref = None
    if _band and _invvar:
        _sigs = []
        for (_i, _j, _o) in _enumerate_bonds(mol):
            if (mol.GetAtomWithIdx(_i).GetSymbol() in metals
                    or mol.GetAtomWithIdx(_j).GetSymbol() in metals):
                continue
            _b = _measured_bond_band(mol, _i, _j)
            if _b is None:
                continue
            _s = (_b[2] - _b[0]) / _P10_P90_TO_SIGMA
            if _s > 1.0e-9:
                _sigs.append(_s)
        if _sigs:
            _sig_ref = float(np.median(_sigs))

    for (i, j, order) in _enumerate_bonds(mol):
        sym_i = mol.GetAtomWithIdx(i).GetSymbol()
        sym_j = mol.GetAtomWithIdx(j).GetSymbol()
        if sym_i in metals or sym_j in metals:
            continue  # handled by U_topology / U_A
        diff = coords[i] - coords[j]
        d = float(np.linalg.norm(diff))
        if d < _EPS_DIST:
            continue
        d_id = _ideal_bond_length(mol, i, j, order)
        if _band:
            # MEASURED band first: centre p50, walls at p10 / p90 of the crystal
            # distribution for this bond's signature.  A flat bottom around the RADII-SUM
            # centre was measured to make discrimination WORSE (0.79 -> 0.59): the builder
            # places bonds at exactly that centre, so we sit at zero cost while reality,
            # being a distribution, does not.  The shape was never the problem; the centre
            # is.  Falls back to the radii sum when the table has no bin for this bond.
            _mb = _measured_bond_band(mol, i, j)
            _w = 1.0
            if _mb is not None:
                _p10, _p50, _p90 = _mb
                if d < _p10:
                    pen, dpen = (_p10 - d) ** 2, -2.0 * (_p10 - d)
                elif d > _p90:
                    pen, dpen = (d - _p90) ** 2, 2.0 * (d - _p90)
                else:
                    pen, dpen = 0.0, 0.0
                if _sig_ref is not None:
                    _s = (_p90 - _p10) / _P10_P90_TO_SIGMA
                    if _s > 1.0e-9:
                        # Bounded: one degenerate bin must not dominate the whole term.
                        _w = min(_INVVAR_W_MAX,
                                 max(1.0 / _INVVAR_W_MAX, (_sig_ref / _s) ** 2))
            else:
                pen, dpen = _flat_bottom(d, d_id, _FLATBOTTOM_FRAC * d_id)
            if pen == 0.0:
                continue
            energy += k_bond * pen * _w
            coef = k_bond * dpen * _w / d
        else:
            delta = d - d_id
            energy += k_bond * delta * delta
            coef = 2.0 * k_bond * delta / d
        grad[i] += coef * diff
        grad[j] -= coef * diff

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_angle
# ---------------------------------------------------------------------------

def U_angle(coords: np.ndarray, mol, k_angle: float = 100.0,
            targets: Optional[List[Tuple[int, int, int, float]]] = None
            ) -> Tuple[float, np.ndarray]:
    """Bond-angle penalty per hybridization: Σ k_angle · (θ - θ_ideal)².

    sp → 180°, sp2 → 120°, sp3 → 109.47°. Metal-centred triples skipped.

    ``targets`` optionally supplies ``(i, j, k, theta_ideal)`` from
    :func:`_build_angle_targets`, i.e. the measured signature angle instead of RDKit's
    three-value hybridization label.  Only the TARGET changes -- the analytic gradient
    below is a derivative of theta and is untouched by where theta_ideal came from.

    ``DELFIN_FFREE_ANGLE_BAND=1`` additionally replaces the harmonic well with a flat
    bottom between the measured p10 and p90 (default OFF -> byte-identical).

    WHY (found 2026-07-30 while building the inverse-variance weight).  The measured angle
    table was wired in for its CENTRE only: ``_signature_theta`` reads p50 and throws p10
    and p90 away, leaving a harmonic point restraint -- an infinitely narrow band.  That is
    exactly the shape the crystal-stationarity measurement rejected for bond lengths: a
    harmonic has ONE zero-force point, so a real crystal, being a distribution, is never
    stationary in it.  The argument transfers verbatim; it had simply never been applied
    here.  A measured centre with an invented width is still half an invention.
    """
    _aband = os.environ.get("DELFIN_FFREE_ANGLE_BAND", "0") == "1"
    _invvar = os.environ.get("DELFIN_FFREE_INVVAR", "0") == "1"
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_angle <= 0.0 or n == 0:
        return float(energy), grad

    _triples = (targets if targets is not None
                else [(i, j, k, None) for (i, j, k) in _enumerate_angles(mol)])

    # Same normalisation as U_bond: weights are RELATIVE within the term, so the balance
    # between the eight terms is untouched.  See the long note in U_bond.
    _sig_ref = None
    if _aband and _invvar:
        _sigs = []
        for _t in _triples:
            _b = _measured_angle_band(mol, _t[0], _t[1], _t[2])
            if _b is None:
                continue
            _s = (_b[2] - _b[0]) / _P10_P90_TO_SIGMA
            if _s > 1.0e-9:
                _sigs.append(_s)
        if _sigs:
            _sig_ref = float(np.median(_sigs))
    for (i, j, k, _theta_target) in _triples:
        x_ij = coords[i] - coords[j]
        x_kj = coords[k] - coords[j]
        norm_ij = float(np.linalg.norm(x_ij))
        norm_kj = float(np.linalg.norm(x_kj))
        if norm_ij < _EPS_DIST or norm_kj < _EPS_DIST:
            continue
        u = x_ij / norm_ij
        v = x_kj / norm_kj
        cos_t = float(np.dot(u, v))
        # Clip into safe arccos / sin window.
        cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        if sin_t < _EPS_SIN:
            continue
        theta = math.acos(cos_t)
        theta_id = (_theta_target if _theta_target is not None
                    else _hybridization_theta(mol, j))
        delta = theta - theta_id
        _w = 1.0
        if _aband:
            _ab = _measured_angle_band(mol, i, j, k)
            if _ab is not None:
                _lo = math.radians(_ab[0])
                _hi = math.radians(_ab[2])
                # Signed excursion past the nearer edge; exactly 0 inside the band, so a
                # crystal is stationary by construction.
                if theta < _lo:
                    delta = theta - _lo
                elif theta > _hi:
                    delta = theta - _hi
                else:
                    continue
                if _sig_ref is not None:
                    _s = (_ab[2] - _ab[0]) / _P10_P90_TO_SIGMA
                    if _s > 1.0e-9:
                        _w = min(_INVVAR_W_MAX,
                                 max(1.0 / _INVVAR_W_MAX, (_sig_ref / _s) ** 2))
        energy += k_angle * _w * delta * delta

        # ∂θ/∂x_i = -1/(sin θ · |x_ij|) · (v - cos θ · u)
        d_th_di = -(v - cos_t * u) / (norm_ij * sin_t)
        d_th_dk = -(u - cos_t * v) / (norm_kj * sin_t)
        d_th_dj = -(d_th_di + d_th_dk)
        coef = 2.0 * k_angle * _w * delta
        grad[i] += coef * d_th_di
        grad[j] += coef * d_th_dj
        grad[k] += coef * d_th_dk

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_torsion — planarity of conjugated bonds (the term the functional never had)
# ---------------------------------------------------------------------------
# U_total minimised EIGHT terms and not one of them was dihedral, so nothing in the
# functional held a conjugated pi system flat: bond, angle and clash are all satisfied by
# a twisted biaryl or a pyramidalised amide.  That is the "non-planar conjugated pi
# systems, especially with heteroatoms in the conjugation" defect, and it had no
# restoring force at all.
#
# The target is the eye's own measured law rather than a second, invented one:
#
#     thr(s) = 90 / (1 + exp((s - s50) / w))          [degrees]
#
# with ``s`` the bond's SHORTENING against that element pair's pure-single reference (pm).
# It is continuous in s and carries no group table -- amide and amidinate land on the same
# threshold at the same s, which is why the group table proved unnecessary.  Parameters are
# fits to the p99 of the pi-capable acyclic band (same provenance as _BOND_ORDER_SCALE's
# note above); only the fitted constants live here, never any structure data.
_TORSION_PI_LAW: Dict[str, Tuple[float, float]] = {
    "C-C": (9.00, 2.30), "C-N": (9.40, 2.00), "C-O": (6.05, 1.00),
    "C-S": (12.55, 2.85), "N-N": (5.85, 1.65),
}
_TORSION_PI_LAW_POOLED = (10.30, 1.95)
# An empty p orbital shortens a bond WITHOUT creating a rotation barrier, so boron is
# never judged by the law.
_TORSION_NO_LAW = frozenset(("B",))
_TORSION_SIN_MIN = 0.342     # sin 20 deg: a near-collinear X-A-B has no torsion plane
_TORSION_THR_MIN = 5.0       # clamp so 1/sin^2(thr) cannot blow up
_TORSION_L1_FALLBACK = 0.98  # pure-single stand-in = 0.98 * (r_cov_i + r_cov_j)
_TORSION_L1_CACHE: Dict[str, Dict[str, float]] = {}


def _torsion_pair_key(a: str, b: str) -> str:
    return "-".join(sorted((a, b)))


def _torsion_l1(sym_a: str, sym_b: str) -> float:
    """Pure-single reference length (Å) for the pair.

    The measured per-pair table lives in the private workspace and is NOT shipped here, so
    the default is the same covalent stand-in the calibration itself falls back to for an
    unseen pair.  ``DELFIN_FFREE_TORSION_L1`` may point at a two-column TSV
    (``pair<TAB>length``) to supply the measured values; without it the law still has the
    right SHAPE but a coarser s axis, which makes the restraint weaker, never wrong-signed.
    """
    # The measured bond-band table already CONTAINS this reference: the pure-single bin's
    # p50 for that element pair IS the pure-single length.  Preferring it closes the loop --
    # one measurement feeds both the length target and the torsion law -- and it fixes a
    # real miscalibration: with the covalent stand-in a biaryl's inter-ring bond looks
    # LONGER than its reference, so s comes out negative, thr -> 89 deg and the law reads
    # a conjugated bond as freely rotating.  Measured, that same bond is shortened.
    _bt = _bond_band_table()
    if _bt:
        _k = "|".join(sorted((sym_a, sym_b))) + "|1.0"
        _v = _bt.get(_k)
        if _v is not None:
            return float(_v[1])                       # p50 of the single-bond bin

    path = os.environ.get("DELFIN_FFREE_TORSION_L1", "")
    if path:
        tbl = _TORSION_L1_CACHE.get(path)
        if tbl is None:
            tbl = {}
            try:
                with open(path, "r") as fh:
                    for ln in fh:
                        parts = ln.split()
                        if len(parts) >= 2:
                            try:
                                tbl[parts[0]] = float(parts[1])
                            except ValueError:
                                continue
            except Exception:
                tbl = {}
            _TORSION_L1_CACHE[path] = tbl
        v = tbl.get(_torsion_pair_key(sym_a, sym_b))
        if v is not None:
            return float(v)
    cov = _smiles_converter_covalent_radii()
    r_a = cov.get(sym_a)
    r_b = cov.get(sym_b)
    if r_a is None or r_b is None:
        return _LAST_RESORT_BOND_LEN
    return _TORSION_L1_FALLBACK * float(r_a + r_b)


def _torsion_threshold_deg(sym_a: str, sym_b: str, d_actual: float) -> Optional[float]:
    """Allowed deviation from planarity (deg) for this bond, or None when the law is off."""
    if sym_a in _TORSION_NO_LAW or sym_b in _TORSION_NO_LAW:
        return None
    s_pm = (_torsion_l1(sym_a, sym_b) - float(d_actual)) * 100.0
    s50, w = _TORSION_PI_LAW.get(_torsion_pair_key(sym_a, sym_b),
                                 _TORSION_PI_LAW_POOLED)
    try:
        return 90.0 / (1.0 + math.exp((s_pm - s50) / w))
    except OverflowError:
        return 0.0 if s_pm > s50 else 90.0


def _build_torsion_targets(coords: np.ndarray, mol
                           ) -> List[Tuple[int, int, int, int, float]]:
    """``(i, j, k, l, weight)`` for every conjugation-restrained dihedral.

    Scope, decided per BOND from its own situation:
      * both ends must be pi-capable -- planar 3-coordinate, or at most 2 sigma partners.
        A tetrahedral end has no p orbital to conjugate through;
      * a ring bond is skipped UNLESS the ring is aromatic.  A saturated ring's torsions are
        already fixed by the ring itself, and restraining them is exactly how a puckered
        ring gets flattened -- a defect the eye has already caught us making;
      * metals are never an endpoint (the coordination sphere is U_A's job);
      * a near-collinear X-A-B substituent is dropped: it has no torsion plane.

    ``weight = 1/sin^2(thr)`` makes the energy read in units of the CRYSTAL-ALLOWED
    deviation -- at exactly the threshold a dihedral contributes 1.  The weight is measured,
    not tuned, and s is frozen from the input geometry so the target cannot drift with the
    coordinates and the gradient stays exact.
    """
    metals = _smiles_converter_metals()
    try:
        ring_info = mol.GetRingInfo()
    except Exception:
        return []
    out: List[Tuple[int, int, int, int, float]] = []

    def _pi_capable(idx: int) -> bool:
        a = mol.GetAtomWithIdx(idx)
        if a.GetDegree() <= 2:
            return True
        return a.GetDegree() == 3 and (_is_pi_active(mol, idx)
                                       or _has_conjugated_lone_pair(mol, idx))

    for bond in mol.GetBonds():
        j = int(bond.GetBeginAtomIdx())
        k = int(bond.GetEndAtomIdx())
        sym_j = mol.GetAtomWithIdx(j).GetSymbol()
        sym_k = mol.GetAtomWithIdx(k).GetSymbol()
        if sym_j in metals or sym_k in metals:
            continue
        if sym_j == "H" or sym_k == "H":
            continue
        try:
            if bond.IsInRing() and not bond.GetIsAromatic():
                continue
        except Exception:
            continue
        if not (_pi_capable(j) and _pi_capable(k)):
            continue
        d_jk = float(np.linalg.norm(coords[j] - coords[k]))
        thr = _torsion_threshold_deg(sym_j, sym_k, d_jk)
        if thr is None:
            continue
        thr = max(_TORSION_THR_MIN, min(90.0, thr))
        weight = 1.0 / (math.sin(math.radians(thr)) ** 2)
        # sin^2 of the tolerance: with DELFIN_FFREE_TORSION_TOL the term costs NOTHING
        # while the dihedral stays inside thr and only bites beyond it.  Measured need:
        # snapping conjugated torsions to exactly planar costs 0.20 A RMSD against real
        # crystals -- biaryls twist 20-40 deg, amides ~10 deg, and the eye's own law calls
        # thr(s) a TOLERANCE of 20-60 deg, not zero.  Forcing zero is over-planarising,
        # the same defect the eye already caught on ring puckers.
        s_thr = math.sin(math.radians(thr)) ** 2
        nbr_j = [int(n.GetIdx()) for n in mol.GetAtomWithIdx(j).GetNeighbors()
                 if int(n.GetIdx()) != k]
        nbr_k = [int(n.GetIdx()) for n in mol.GetAtomWithIdx(k).GetNeighbors()
                 if int(n.GetIdx()) != j]
        for i in nbr_j:
            for l in nbr_k:
                if i == l:
                    continue
                # Drop substituents with no torsion plane (near-collinear).
                b1 = coords[j] - coords[i]
                b2 = coords[k] - coords[j]
                b3 = coords[l] - coords[k]
                n1 = float(np.linalg.norm(np.cross(b1, b2)))
                n2 = float(np.linalg.norm(np.cross(b2, b3)))
                d1 = float(np.linalg.norm(b1)) * float(np.linalg.norm(b2))
                d2 = float(np.linalg.norm(b2)) * float(np.linalg.norm(b3))
                if d1 < _EPS_DIST or d2 < _EPS_DIST:
                    continue
                if n1 / d1 < _TORSION_SIN_MIN or n2 / d2 < _TORSION_SIN_MIN:
                    continue
                out.append((i, j, k, l, weight, s_thr))
    return out


def U_torsion(coords: np.ndarray, mol, k_torsion: float = 0.0,
              targets: Optional[List[Tuple[int, int, int, int, float]]] = None
              ) -> Tuple[float, np.ndarray]:
    """Planarity penalty on conjugated dihedrals: Σ k · w · sin²(φ).

    ``sin²φ`` has its minima at 0° and 180°, i.e. at BOTH planar arrangements, so the term
    flattens a twisted pi system without choosing cis or trans for it -- the isomer stays
    whatever the enumeration made it.  Gradient is analytic.
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0
    _tol_mode = os.environ.get("DELFIN_FFREE_TORSION_TOL", "0") == "1"
    if k_torsion <= 0.0 or coords.shape[0] == 0 or not targets:
        return float(energy), grad

    for _t in targets:
        i, j, k, l, weight = _t[:5]
        s_thr = _t[5] if len(_t) > 5 else 0.0
        b1 = coords[j] - coords[i]
        b2 = coords[k] - coords[j]
        b3 = coords[l] - coords[k]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        n1_sq = float(np.dot(n1, n1))
        n2_sq = float(np.dot(n2, n2))
        if n1_sq < _EPS_DIST or n2_sq < _EPS_DIST:
            continue
        # sin^2(phi) = 1 - cos^2(phi) with cos(phi) = (n1·n2)/(|n1||n2|).  Writing the
        # energy in n1/n2 instead of in phi removes the dihedral sign convention from the
        # problem entirely: the textbook dphi/dx formulas are stated for a b2 that points
        # the other way, and mixing the conventions still satisfies translation invariance
        # (the four terms sum to zero either way), so an inspection cannot catch it -- only
        # the FD check can, and it did, twice, at 2.1e-01 and 7.6e-02.  Differentiating the
        # cross products directly is mechanical and has no convention to get wrong.
        n1_len = math.sqrt(n1_sq)
        n2_len = math.sqrt(n2_sq)
        dot_n = float(np.dot(n1, n2))
        c = dot_n / (n1_len * n2_len)
        c = max(-1.0, min(1.0, c))
        s = 1.0 - c * c                       # = sin^2(phi), the planarity deviation
        if _tol_mode and s_thr > 0.0:
            # Zero cost inside the measured tolerance, quadratic in the EXCESS beyond it.
            # C1 at the join, and the derivative chains through s exactly as below.
            if s <= s_thr:
                continue
            ex = s - s_thr
            energy += k_torsion * weight * ex * ex
            dU_ds = 2.0 * k_torsion * weight * ex
        else:
            energy += k_torsion * weight * s
            dU_ds = k_torsion * weight
        # ds/dc = -2c, so dU/dc = dU_ds * (-2c)
        dU_dc = -2.0 * c * dU_ds
        inv = 1.0 / (n1_len * n2_len)
        # dc/dn1 and dc/dn2
        g = dU_dc * inv * (n2 - (dot_n / n1_sq) * n1)
        h = dU_dc * inv * (n1 - (dot_n / n2_sq) * n2)
        # Chain through n1 = b1 x b2 and n2 = b2 x b3 using a·(b x c) = b·(c x a).
        dU_db1 = np.cross(b2, g)
        dU_db2 = np.cross(g, b1) + np.cross(b3, h)
        dU_db3 = np.cross(h, b2)
        # b1 = r_j - r_i, b2 = r_k - r_j, b3 = r_l - r_k
        grad[i] -= dU_db1
        grad[j] += dU_db1 - dU_db2
        grad[k] += dU_db2 - dU_db3
        grad[l] += dU_db3

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_clash
# ---------------------------------------------------------------------------

def U_clash(coords: np.ndarray, mol, k_clash: float = 500.0,
            clash_factor: float = 0.85, h_h_factor: float = 0.75
            ) -> Tuple[float, np.ndarray]:
    """One-sided clash repulsion for non-bonded heavy pairs.

    Uses Bondi/Alvarez vdW sums from :mod:`delfin.manta._vdw_radii`. Metal/metal
    and 1-2, 1-3 connections are excluded. Penalty is C¹-continuous quadratic
    inside the overlap region; zero (with zero gradient) outside.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_clash <= 0.0 or n < 2:
        return float(energy), grad

    metals = _smiles_converter_metals()
    syms = _symbols(mol)
    bonded = _bonded_pairs_set(mol)
    angle13 = _1_3_pairs_set(mol)

    # Pre-tabulate per-element vdW once for speed.
    vdw_cache: Dict[str, float] = {}
    for s in set(syms):
        vdw_cache[s] = _vdw_radius(s)

    for i in range(n):
        si = syms[i]
        if si in metals:
            continue
        for j in range(i + 1, n):
            sj = syms[j]
            if sj in metals:
                continue
            if (i, j) in bonded or (i, j) in angle13:
                continue
            # Threshold scaling: H-H gets looser factor.
            if si == "H" and sj == "H":
                threshold = (vdw_cache[si] + vdw_cache[sj]) * h_h_factor
            else:
                threshold = (vdw_cache[si] + vdw_cache[sj]) * clash_factor
            diff = coords[j] - coords[i]
            d = float(np.linalg.norm(diff))
            if d < _EPS_DIST:
                # Atoms collapsed — push apart along arbitrary axis (x).
                overlap = threshold
                energy += k_clash * overlap * overlap
                push = np.array([2.0 * k_clash * overlap, 0.0, 0.0])
                grad[i] -= push
                grad[j] += push
                continue
            if d >= threshold:
                continue
            overlap = threshold - d
            energy += k_clash * overlap * overlap
            # ∂U/∂x_i = +2 k_c · overlap · (x_j - x_i) / d
            # (moving x_i toward x_j shortens d, U=k·(threshold-d)² grows
            #  → positive gradient in direction (x_j - x_i). L-BFGS-B
            #  moves -grad → pushes x_i AWAY from x_j.)
            push = (2.0 * k_clash * overlap / d) * diff  # = ∂U/∂x_i
            grad[i] += push
            grad[j] -= push

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_topology
# ---------------------------------------------------------------------------

# Width over which -log(t) is replaced by its 2nd-order Taylor continuation, as a
# fraction of the band width, so the relaxation scales with the band instead of being an
# absolute length that means something different for a 1.9 A and a 2.6 A bond.
_BARRIER_RELAX_FRAC = 0.05


def _relaxed_barrier_enabled() -> bool:
    return os.environ.get("DELFIN_FFREE_RELAXED_BARRIER", "0") == "1"


def _relaxed_log_barrier(t: float, t_rel: float) -> Tuple[float, float]:
    """``(-log(t), d/dt)`` for ``t >= t_rel``; below that, its 2nd-order Taylor at
    ``t_rel``.  C2-continuous at the join, convex, and finite for every t including
    t <= 0 -- so a frame that starts outside the band still sits on a real surface with
    a real gradient, and can be pulled back in instead of stranded on a plateau."""
    if t >= t_rel:
        return -math.log(t), -1.0 / t
    d = t - t_rel
    inv = 1.0 / t_rel
    return (-math.log(t_rel) - d * inv + 0.5 * d * d * inv * inv,
            -inv + d * inv * inv)


def U_topology(coords: np.ndarray, mol, k_topology: float = 10000.0,
               lo_frac: float = 0.85, hi_frac: float = 1.10
               ) -> Tuple[float, np.ndarray]:
    """Log-barrier on every M-D bond.

    For each M-D pair the bond length ``d`` must remain inside
    ``[lo_frac, hi_frac] · d_ideal``. Inside the open interval the barrier is

        barrier(d) = -log(d - lo) - log(hi - d)

    which diverges at the boundaries (preserving topology by construction).
    Outside the interval a large finite penalty with a repulsive linear
    gradient steers L-BFGS-B back into the safe domain (the algorithm needs
    finite values to make progress).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_topology <= 0.0 or n == 0:
        return float(energy), grad

    big_pen = 1.0e6  # placeholder energy when outside safe interval

    for (m_idx, d_idx) in _enumerate_metal_donor_bonds(mol):
        diff = coords[m_idx] - coords[d_idx]
        d = float(np.linalg.norm(diff))
        if d < _EPS_DIST:
            # Atoms collapsed — emit huge gradient to separate them.
            energy += k_topology * big_pen
            push = np.array([k_topology * big_pen, 0.0, 0.0])
            grad[m_idx] += push
            grad[d_idx] -= push
            continue
        sym_m = mol.GetAtomWithIdx(m_idx).GetSymbol()
        sym_d = mol.GetAtomWithIdx(d_idx).GetSymbol()
        d_id = _smiles_converter_ml_bondlen(sym_m, sym_d)
        lo = d_id * lo_frac
        hi = d_id * hi_frac

        # RELAXED BARRIER (DELFIN_FFREE_RELAXED_BARRIER, default OFF).
        # The branch below reports a CONSTANT energy outside the band together with a
        # huge gradient, i.e. g is not the derivative of f.  Measured in the self-test:
        # at 0.80x and at 1.25x ideal the energy is identical to the last decimal while
        # |g_an| = 9.0e9 against a true |g_fd| = 1.1e5 -- relative error 1.000.  Outside
        # its band this is not a potential energy surface at all, L-BFGS-B's Wolfe line
        # search cannot be satisfied, and the "minimisation" is not one.  It matters
        # because EVERY measured frame starts outside (topo_ok_input=False on all).
        # The standard remedy: continue -log(t) below a small width by its 2nd-order
        # Taylor expansion -- C2-continuous, convex, finite everywhere, and g == grad f.
        # The wall stops being infinite; topology is held by the (now relative) post-gate
        # instead, and the surface can pull an out-of-band frame smoothly back IN.
        # FLAT INTERIOR (DELFIN_FFREE_TOPO_BAND=1, default OFF).  The log barrier's
        # gradient vanishes ONLY at the band centre, i.e. at 0.975 x our reference M-D:
        # it is not a barrier, it is a stiff restraint to one length.  Measured at 107
        # real crystals it outweighs every other term ~1000x and scores OUR OWN FRAMES
        # FIVE TIMES BETTER THAN REALITY (RMS force 10786 vs 52833) -- because the builder
        # places M-D at exactly that reference, so we start in its minimum and crystals,
        # being a distribution, do not.  Packing forces explain excess force at a crystal;
        # they cannot explain a ratio below 1.  Flat bottom over the WHOLE band with
        # quadratic walls outside: zero force wherever the crystal actually lives, and the
        # hard topology guarantee sits in the post-gate, which is relative now anyway.
        if os.environ.get("DELFIN_FFREE_TOPO_BAND", "0") == "1":
            # MEASURED M-D band first, for exactly the reason it worked on U_bond: with the
            # radii/table reference the flat bottom is centred on OUR number, which the
            # builder already hits, so we sit at zero cost and reality does not.  Feeding
            # U_bond the measured centre took its crystal force from 195.85 to 16.03 and
            # its discrimination from 0.79 to 14.49; this is the same change on the term
            # that carries ~1000x the weight.
            _mb = _measured_bond_band(mol, m_idx, d_idx)
            if _mb is not None:
                _p10, _p50, _p90 = _mb
                if d < _p10:
                    pen, dpen = (_p10 - d) ** 2, -2.0 * (_p10 - d)
                elif d > _p90:
                    pen, dpen = (d - _p90) ** 2, 2.0 * (d - _p90)
                else:
                    pen, dpen = 0.0, 0.0
            else:
                _mid = 0.5 * (lo + hi)
                pen, dpen = _flat_bottom(d, _mid, 0.5 * (hi - lo))
            if pen != 0.0:
                energy += k_topology * pen
                coef = k_topology * dpen / d
                grad[m_idx] += coef * diff
                grad[d_idx] -= coef * diff
            continue

        if _relaxed_barrier_enabled():
            t_rel = _BARRIER_RELAX_FRAC * (hi - lo)
            b_lo, db_lo = _relaxed_log_barrier(d - lo, t_rel)
            b_hi, db_hi = _relaxed_log_barrier(hi - d, t_rel)
            energy += k_topology * (b_lo + b_hi)
            db_dd = db_lo - db_hi          # d(hi - d)/dd = -1
            coef = k_topology * db_dd / d
            grad[m_idx] += coef * diff
            grad[d_idx] -= coef * diff
            continue

        # Outside-domain handling: very large penalty + linear gradient
        # pointing back into the safe interval.
        if d <= lo + _EPS_BARRIER:
            energy += k_topology * big_pen
            # We need d to INCREASE → x_m should move AWAY from x_d.
            # ∂(d)/∂x_m = (x_m - x_d)/d, ∂(d)/∂x_d = -(...)
            # Use -∂U/∂d > 0 ⇒ force on m is along +(x_m - x_d)/d.
            coef = -k_topology * big_pen / max(d, _EPS_DIST)
            # coef·diff acts as +∂U/∂x_m; gradient = +coef·diff so the
            # minimiser moves m in direction +diff (away from donor).
            grad[m_idx] += coef * diff
            grad[d_idx] -= coef * diff
            continue
        if d >= hi - _EPS_BARRIER:
            energy += k_topology * big_pen
            # We need d to DECREASE → m moves TOWARD donor.
            coef = +k_topology * big_pen / max(d, _EPS_DIST)
            grad[m_idx] += coef * diff
            grad[d_idx] -= coef * diff
            continue

        # Smooth log barrier inside (lo, hi).
        b = -math.log(d - lo) - math.log(hi - d)
        energy += k_topology * b
        # ∂b/∂d = -1/(d - lo) + 1/(hi - d)
        db_dd = -1.0 / (d - lo) + 1.0 / (hi - d)
        coef = k_topology * db_dd / d
        grad[m_idx] += coef * diff
        grad[d_idx] -= coef * diff

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_A — Coord-sphere polyhedron pull (Tier A)
# ---------------------------------------------------------------------------

def U_A_coord_sphere(coords: np.ndarray, mol,
                     donor_targets: Dict[int, np.ndarray],
                     k_A: float = 100.0) -> Tuple[float, np.ndarray]:
    """Tier A: pull each donor toward its Hungarian-assigned slot.

    ``donor_targets[donor_idx]`` is a 3-vector giving the absolute target
    position (already in the working coordinate frame). The metal centres
    themselves are not pulled; their position is governed by the rest of
    the force field.
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_A <= 0.0 or not donor_targets:
        return float(energy), grad

    for donor_idx, target in donor_targets.items():
        if donor_idx < 0 or donor_idx >= coords.shape[0]:
            continue
        tgt = np.asarray(target, dtype=np.float64).reshape(3)
        delta = coords[donor_idx] - tgt
        energy += k_A * float(np.dot(delta, delta))
        grad[donor_idx] += 2.0 * k_A * delta

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_B — Local equivalence (Tier B)
# ---------------------------------------------------------------------------

def _angle_value(coords: np.ndarray, i: int, j: int, k: int) -> float:
    """Return the angle θ_ijk in radians, safely clipped."""
    x_ij = coords[i] - coords[j]
    x_kj = coords[k] - coords[j]
    nij = float(np.linalg.norm(x_ij))
    nkj = float(np.linalg.norm(x_kj))
    if nij < _EPS_DIST or nkj < _EPS_DIST:
        return _THETA_SP3
    cos_t = float(np.dot(x_ij, x_kj) / (nij * nkj))
    cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
    return math.acos(cos_t)


def U_B_equivalence(coords: np.ndarray, mol,
                    equiv_bond_pairs: List[List[Tuple[int, int]]],
                    equiv_angle_triples: List[List[Tuple[int, int, int]]],
                    k_B: float = 50.0) -> Tuple[float, np.ndarray]:
    """Tier B: chemically-equivalent bonds/angles pull toward their class mean.

    ``equiv_bond_pairs`` is a list of bond-classes, each a list of
    ``(i, j)`` index tuples.  ``equiv_angle_triples`` is the analogous list
    for angle triples ``(i, j, k)``.

    For tractable gradients the class mean is treated as a constant of the
    current iterate; L-BFGS-B will re-evaluate (mean,gradient) at every step
    so the mean updates implicitly between iterations.
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_B <= 0.0:
        return float(energy), grad

    # Bond equivalence classes.
    for cls in equiv_bond_pairs or []:
        if not cls or len(cls) < 2:
            continue
        ds: List[float] = []
        diffs: List[np.ndarray] = []
        ok: List[Tuple[int, int]] = []
        for (i, j) in cls:
            diff = coords[i] - coords[j]
            d = float(np.linalg.norm(diff))
            if d < _EPS_DIST:
                continue
            ds.append(d)
            diffs.append(diff)
            ok.append((i, j))
        if len(ds) < 2:
            continue
        mean = float(np.mean(ds))
        for (i, j), d, diff in zip(ok, ds, diffs):
            delta = d - mean
            energy += k_B * delta * delta
            coef = 2.0 * k_B * delta / d
            grad[i] += coef * diff
            grad[j] -= coef * diff

    # Angle equivalence classes.
    for cls in equiv_angle_triples or []:
        if not cls or len(cls) < 2:
            continue
        thetas: List[float] = []
        ok_triples: List[Tuple[int, int, int]] = []
        for (i, j, k) in cls:
            t = _angle_value(coords, i, j, k)
            thetas.append(t)
            ok_triples.append((i, j, k))
        if len(thetas) < 2:
            continue
        mean = float(np.mean(thetas))
        for (i, j, k), theta in zip(ok_triples, thetas):
            x_ij = coords[i] - coords[j]
            x_kj = coords[k] - coords[j]
            nij = float(np.linalg.norm(x_ij))
            nkj = float(np.linalg.norm(x_kj))
            if nij < _EPS_DIST or nkj < _EPS_DIST:
                continue
            u = x_ij / nij
            v = x_kj / nkj
            cos_t = float(np.dot(u, v))
            cos_t = max(-1.0 + _EPS_COS, min(1.0 - _EPS_COS, cos_t))
            sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
            if sin_t < _EPS_SIN:
                continue
            delta = theta - mean
            energy += k_B * delta * delta
            d_th_di = -(v - cos_t * u) / (nij * sin_t)
            d_th_dk = -(u - cos_t * v) / (nkj * sin_t)
            d_th_dj = -(d_th_di + d_th_dk)
            coef = 2.0 * k_B * delta
            grad[i] += coef * d_th_di
            grad[j] += coef * d_th_dj
            grad[k] += coef * d_th_dk

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_C — Per-fragment archetype point group (Tier C)
# ---------------------------------------------------------------------------

def U_C_fragment(coords: np.ndarray, mol, fragments: List[Dict],
                 k_C: float = 80.0) -> Tuple[float, np.ndarray]:
    """Tier C: per-fragment archetype point-group enforcement.

    Each ``fragments[k]`` is a dict carrying at least the keys::

        {
          "atoms": [i0, i1, ...],          # atom indices in the fragment
          "operations": [op_mat_1, ...],   # list of 3×3 numpy arrays (skip identity)
          "partners":  [{i: j, ...}, ...]  # parallel to operations: i → partner j under op
        }

    For every (op, partner-mapping) pair we minimise

        U += k_C · Σ_i ||x_{p(i)} - (centroid + op·(x_i - centroid))||²

    where ``centroid`` is the running fragment centroid. The chain-rule
    derivative w.r.t. the fragment centroid is distributed equally back to
    all fragment atoms (which keeps L-BFGS-B well-conditioned).
    """
    coords = np.asarray(coords, dtype=np.float64)
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_C <= 0.0 or not fragments:
        return float(energy), grad

    for frag in fragments:
        atoms = list(frag.get("atoms", []))
        ops = list(frag.get("operations", []))
        perms = list(frag.get("partners", []))
        if not atoms or not ops or not perms or len(ops) != len(perms):
            continue
        n_frag = len(atoms)
        if n_frag < 2:
            continue
        # Build sub-array of fragment coords once per fragment.
        frag_arr = coords[atoms]
        centroid = frag_arr.mean(axis=0)
        inv_n = 1.0 / float(n_frag)

        for op, perm in zip(ops, perms):
            if op is None or perm is None:
                continue
            op_mat = np.asarray(op, dtype=np.float64).reshape(3, 3)
            # Skip identity (it contributes zero penalty).
            if np.allclose(op_mat, np.eye(3), atol=1.0e-9):
                continue
            op_T = op_mat.T
            # Accumulator for the centroid back-distribution.
            centroid_force = np.zeros(3, dtype=np.float64)
            atom_set = set(atoms)
            for i, j in perm.items():
                if i not in atom_set or j not in atom_set:
                    continue
                xi_loc = coords[i] - centroid
                target = centroid + op_mat @ xi_loc
                delta = coords[j] - target
                if i == j:
                    # Atom is fixed by op (e.g. on a mirror plane / inversion centre).
                    # We still penalise residual displacement; ∂target/∂x_i = op,
                    # but x_j = x_i so the gradient combines.
                    energy += k_C * float(np.dot(delta, delta))
                    # ∂U/∂x_i = 2 k_C · (I - opᵀ) · delta   (j == i)
                    grad[i] += 2.0 * k_C * ((np.eye(3) - op_T) @ delta)
                    continue
                energy += k_C * float(np.dot(delta, delta))
                # ∂U/∂x_j = +2 k_C · delta
                grad[j] += 2.0 * k_C * delta
                # ∂U/∂x_i: target = centroid + op·(x_i - centroid)
                #          → ∂target/∂x_i = op
                # ∂U/∂x_i = -2 k_C · opᵀ · delta
                grad[i] -= 2.0 * k_C * (op_T @ delta)
                # Centroid contribution: ∂target/∂c = I - op
                # ⇒ centroid_force accumulates -2 k_C · (I - opᵀ) · delta
                centroid_force -= 2.0 * k_C * ((np.eye(3) - op_T) @ delta)
            # Distribute centroid gradient equally to all fragment atoms.
            if not np.allclose(centroid_force, 0.0):
                share = inv_n * centroid_force
                for a in atoms:
                    grad[a] += share

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_D — Global molecular point group (Tier D)
# ---------------------------------------------------------------------------

def U_D_global(coords: np.ndarray, mol, global_pg_ops: List[np.ndarray],
               atom_perms: Dict[int, Dict[int, int]],
               k_D: float = 30.0) -> Tuple[float, np.ndarray]:
    """Tier D: global molecular point-group enforcement.

    ``global_pg_ops`` is a list of 3×3 rotation/reflection matrices.
    ``atom_perms[op_index]`` is a dict ``{atom_idx: partner_atom_idx}``
    encoding how the operation permutes atoms (identity ops should be
    omitted or will be skipped automatically).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[0]
    grad = np.zeros_like(coords)
    energy = 0.0

    if k_D <= 0.0 or n == 0 or not global_pg_ops:
        return float(energy), grad

    centroid = coords.mean(axis=0)
    inv_n = 1.0 / float(n)

    for op_idx, op in enumerate(global_pg_ops):
        if op is None:
            continue
        op_mat = np.asarray(op, dtype=np.float64).reshape(3, 3)
        if np.allclose(op_mat, np.eye(3), atol=1.0e-9):
            continue
        perm = atom_perms.get(op_idx, {}) if atom_perms else {}
        if not perm:
            continue
        op_T = op_mat.T
        centroid_force = np.zeros(3, dtype=np.float64)
        for i, j in perm.items():
            if i < 0 or i >= n or j < 0 or j >= n:
                continue
            xi_loc = coords[i] - centroid
            target = centroid + op_mat @ xi_loc
            delta = coords[j] - target
            if i == j:
                # Atom fixed by op (lies on symmetry element).
                energy += k_D * float(np.dot(delta, delta))
                grad[i] += 2.0 * k_D * ((np.eye(3) - op_T) @ delta)
                continue
            energy += k_D * float(np.dot(delta, delta))
            grad[j] += 2.0 * k_D * delta
            grad[i] -= 2.0 * k_D * (op_T @ delta)
            centroid_force -= 2.0 * k_D * ((np.eye(3) - op_T) @ delta)
        if not np.allclose(centroid_force, 0.0):
            share = inv_n * centroid_force
            grad += share  # broadcast to every atom

    return float(energy), grad


# ---------------------------------------------------------------------------
# U_total
# ---------------------------------------------------------------------------

def U_total(coords: np.ndarray, mol, sym_info: Dict, params: Dict
            ) -> Tuple[float, np.ndarray]:
    """Aggregate all 8 terms with class-conditional k values.

    Parameters
    ----------
    coords : (N, 3) numpy array
        Current Cartesian coordinates (Å).
    mol : rdkit Mol
        Molecule with bonds/hybridization.
    sym_info : dict
        Pre-computed Tier A-D data. Recognised keys:
            ``donor_targets``        — dict[int, (3,)] for U_A
            ``equiv_bond_pairs``     — list[list[(i,j)]] for U_B
            ``equiv_angle_triples``  — list[list[(i,j,k)]] for U_B
            ``fragments``            — list[dict] for U_C
            ``global_ops``           — list[(3,3)] for U_D
            ``atom_perms``           — dict[int, dict[int,int]] for U_D
    params : dict
        Coefficients (default to spec values if missing):
            ``k_bond, k_angle, k_clash, k_topology, k_A, k_B, k_C, k_D``
        Plus optional ``clash_factor, h_h_factor, lo_frac, hi_frac``.

    Returns
    -------
    (energy_total, grad_total)
    """
    coords = np.asarray(coords, dtype=np.float64)
    sym_info = sym_info or {}
    params = params or {}

    k_bond = float(params.get("k_bond", 1000.0))
    k_angle = float(params.get("k_angle", 100.0))
    k_clash = float(params.get("k_clash", 500.0))
    k_topology = float(params.get("k_topology", 10000.0))
    k_A = float(params.get("k_A", 100.0))
    k_B = float(params.get("k_B", 50.0))
    k_C = float(params.get("k_C", 80.0))
    k_D = float(params.get("k_D", 30.0))
    clash_factor = float(params.get("clash_factor", 0.85))
    h_h_factor = float(params.get("h_h_factor", 0.75))
    lo_frac = float(params.get("lo_frac", 0.85))
    hi_frac = float(params.get("hi_frac", 1.10))

    e_total = 0.0
    g_total = np.zeros_like(coords)

    e, g = U_bond(coords, mol, k_bond=k_bond)
    e_total += e
    g_total += g

    # Signature-resolved angle targets (DELFIN_FFREE_THETA_SIGNATURE, default OFF ->
    # byte-identical).  Memoised in sym_info, which is created once per refinement, so the
    # signature work is paid once and not once per L-BFGS iteration; a plain dict also
    # means no global state and no id()-reuse hazard.
    _angle_targets = None
    if os.environ.get("DELFIN_FFREE_THETA_SIGNATURE", "0") == "1":
        _angle_targets = sym_info.get("_angle_targets")
        if _angle_targets is None:
            try:
                _angle_targets = _build_angle_targets(mol)
            except Exception:
                _angle_targets = []
            sym_info["_angle_targets"] = _angle_targets
        if not _angle_targets:
            _angle_targets = None
    e, g = U_angle(coords, mol, k_angle=k_angle, targets=_angle_targets)
    e_total += e
    g_total += g

    # Conjugation planarity (DELFIN_FFREE_TORSION, default OFF -> byte-identical).  The
    # scope and the per-bond weight are frozen from the INPUT geometry, memoised in the
    # per-refine sym_info, so the target cannot drift with the coordinates mid-minimisation
    # and the analytic gradient stays exact.
    k_torsion = float(params.get("k_torsion", 0.0))
    if k_torsion > 0.0 and os.environ.get("DELFIN_FFREE_TORSION", "0") == "1":
        _tors = sym_info.get("_torsion_targets")
        if _tors is None:
            try:
                _tors = _build_torsion_targets(coords, mol)
            except Exception:
                _tors = []
            sym_info["_torsion_targets"] = _tors
        if _tors:
            e, g = U_torsion(coords, mol, k_torsion=k_torsion, targets=_tors)
            e_total += e
            g_total += g

    e, g = U_clash(coords, mol, k_clash=k_clash,
                   clash_factor=clash_factor, h_h_factor=h_h_factor)
    e_total += e
    g_total += g

    e, g = U_topology(coords, mol, k_topology=k_topology,
                      lo_frac=lo_frac, hi_frac=hi_frac)
    e_total += e
    g_total += g

    e, g = U_A_coord_sphere(coords, mol,
                            sym_info.get("donor_targets", {}) or {},
                            k_A=k_A)
    e_total += e
    g_total += g

    e, g = U_B_equivalence(coords, mol,
                           sym_info.get("equiv_bond_pairs", []) or [],
                           sym_info.get("equiv_angle_triples", []) or [],
                           k_B=k_B)
    e_total += e
    g_total += g

    # U_C AND U_D REMOVED (2026-07-30), because they were MEASURED to do nothing.
    # Over 107 real crystals and our own frames the stationarity test reports:
    #     U_C   0.0000 / 0.0000
    #     U_D   0.0000 / 0.0000
    # Exactly zero force on BOTH sides.  A term whose gradient vanishes everywhere we
    # looked cannot change any outcome, so it is not a term, it is cost.
    #
    # And the reason is chemistry, not a fixable implementation detail: every complex
    # inspected came out as pg=C1.  Real coordination complexes have no global point group,
    # so U_D has no operation to apply; and a real ligand inside a complex is distorted
    # enough that the fragment-orbit guard (residual <= 0.35 A) rejects its automorphisms,
    # so U_C has no fragment to enforce.  That includes the automorphism rewrite made
    # earlier the same day -- the measurement says it changed nothing on real systems.
    #
    # Removing them also drops the Tier C/D half of _precompute_symmetry, which that
    # module's own docstring calls "one-shot, expensive".  So this is not merely dead
    # weight removed, it is time given back.
    #
    # Standing rule: corrections are scaffolding and are MEANT to fall; at equal outcome,
    # less code is strictly better.  k_C / k_D stay in the hyperparameter presets so an
    # ablation A/B can still be read against the old logs.
    #
    # DELFIN_FFREE_TIER_CD=1 puts them BACK.  A removal has to be confirmable, and it can
    # only be confirmed if both versions are runnable: the ablation arm sets this and, if
    # the byte partition comes out identical, the removal is proven to change nothing --
    # the strongest form of never-worse there is.  The flag exists to be deleted once that
    # measurement is in.
    if os.environ.get("DELFIN_FFREE_TIER_CD", "0") == "1":
        e, g = U_C_fragment(coords, mol,
                            sym_info.get("fragments", []) or [], k_C=k_C)
        e_total += e
        g_total += g
        e, g = U_D_global(coords, mol,
                          sym_info.get("global_ops", []) or [],
                          sym_info.get("atom_perms", {}) or {}, k_D=k_D)
        e_total += e
        g_total += g

    return float(e_total), g_total


# ---------------------------------------------------------------------------
# Self-test (smoke check + finite-difference gradient sanity)
# ---------------------------------------------------------------------------

def _fd_gradient(fn, coords: np.ndarray, eps: float = 1.0e-5) -> np.ndarray:
    """Central-difference numerical gradient of ``fn(coords) -> (energy, _)``."""
    g = np.zeros_like(coords)
    flat = coords.reshape(-1).copy()
    for k in range(flat.size):
        orig = flat[k]
        flat[k] = orig + eps
        e_p, _ = fn(flat.reshape(coords.shape))
        flat[k] = orig - eps
        e_m, _ = fn(flat.reshape(coords.shape))
        flat[k] = orig
        g.reshape(-1)[k] = (e_p - e_m) / (2.0 * eps)
    return g


if __name__ == "__main__":  # pragma: no cover
    # 4-atom synthetic Cu-N₂-O system (mock geometry; not chemistry-accurate).
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"RDKit required for self-test: {exc}")

    smi = "[Cu](N)(N)O"
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
                      dtype=np.float64)

    print(f"Self-test — {mol.GetNumAtoms()} atoms ({smi})\n")

    e, g = U_bond(coords, mol)
    print(f"U_bond      E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_angle(coords, mol)
    print(f"U_angle     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_clash(coords, mol)
    print(f"U_clash     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")
    e, g = U_topology(coords, mol)
    print(f"U_topology  E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Tier A: pull every non-metal donor (N, N, O) toward its current pos
    # plus a tiny offset so the term is non-trivial.
    metals = _smiles_converter_metals()
    donor_targets = {}
    for a in mol.GetAtoms():
        if a.GetSymbol() in metals:
            continue
        if not any(n.GetSymbol() in metals for n in a.GetNeighbors()):
            continue
        i = a.GetIdx()
        donor_targets[i] = coords[i] + np.array([0.05, 0.0, 0.0])
    e, g = U_A_coord_sphere(coords, mol, donor_targets)
    print(f"U_A         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}"
          f"  (n_donors={len(donor_targets)})")

    # Trivial Tier B: pair every M-D bond with itself (n=2) – use first two.
    md = _enumerate_metal_donor_bonds(mol)
    if len(md) >= 2:
        equiv_bonds = [[md[0], md[1]]]
    else:
        equiv_bonds = []
    e, g = U_B_equivalence(coords, mol, equiv_bonds, [])
    print(f"U_B         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}"
          f"  (n_pairs={len(equiv_bonds)})")

    # Trivial Tier C: single fragment with C2 rotation (identity-like partner).
    # We just exercise the code path; with partner=={i:i} U_C should give 0
    # (delta = 0 because target = x_i for identity rotation around centroid
    # only if op = I; choose op = -I to make it non-zero but symmetric).
    frag = {
        "atoms": list(range(mol.GetNumAtoms())),
        "operations": [-np.eye(3)],
        "partners": [{i: i for i in range(mol.GetNumAtoms())}],
    }
    e, g = U_C_fragment(coords, mol, [frag])
    print(f"U_C         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Trivial Tier D: identity-permutation under -I.
    e, g = U_D_global(coords, mol, [-np.eye(3)],
                      {0: {i: i for i in range(mol.GetNumAtoms())}})
    print(f"U_D         E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    sym_info = {
        "donor_targets": donor_targets,
        "equiv_bond_pairs": equiv_bonds,
        "equiv_angle_triples": [],
        "fragments": [frag],
        "global_ops": [],
        "atom_perms": {},
    }
    params = {}
    e, g = U_total(coords, mol, sym_info, params)
    print(f"U_total     E = {e:12.4f}   |g|_max = {np.max(np.abs(g)):.4e}")

    # Finite-difference cross-check (single coord, all terms aggregated).
    def total_fn(c):
        return U_total(c, mol, sym_info, params)

    g_an = g
    g_fd = _fd_gradient(total_fn, coords, eps=1.0e-5)
    diff = np.max(np.abs(g_an - g_fd))
    rel = diff / max(1.0e-9, float(np.max(np.abs(g_an))))
    print(f"\nFD gradient check on U_total:")
    print(f"  max|grad_an - grad_fd|  = {diff:.4e}")
    print(f"  relative                = {rel:.4e}")
    print(f"  PASS" if rel < 1.0e-3 else f"  FAIL (rel >= 1e-3)")

    # ------------------------------------------------------------------
    # IS THIS AN OPTIMISABLE SURFACE WHERE OUR FRAMES ACTUALLY LIVE?
    # The FD check above sits INSIDE U_topology's barrier domain.  Every measured frame
    # sits OUTSIDE it (topo_ok_input=False on all of them), and out there the barrier
    # reports a CONSTANT energy (k · big_pen) together with a huge gradient.  If f is flat
    # and g is not zero, g is not the derivative of f, L-BFGS-B's Wolfe line search cannot
    # be satisfied, and the "minimisation" is not one.  Test it where it matters.
    # ------------------------------------------------------------------
    print("\nU_topology OUTSIDE the barrier band (where our frames are):")
    _md = _enumerate_metal_donor_bonds(mol)
    if _md:
        _m_i, _d_i = _md[0]
        _sm = mol.GetAtomWithIdx(_m_i).GetSymbol()
        _sd = mol.GetAtomWithIdx(_d_i).GetSymbol()
        _ideal = _smiles_converter_ml_bondlen(_sm, _sd)
        for _frac, _where in ((0.99, "inside  (0.99 x ideal)"),
                              (0.80, "OUTSIDE (0.80 x ideal, below lo=0.85)"),
                              (1.25, "OUTSIDE (1.25 x ideal, above hi=1.10)")):
            _c = coords.copy()
            _v = _c[_d_i] - _c[_m_i]
            _n = float(np.linalg.norm(_v))
            if _n < 1e-9:
                continue
            _c[_d_i] = _c[_m_i] + (_v / _n) * (_ideal * _frac)

            def _topo_fn(cc):
                return U_topology(cc, mol)

            _e0, _g0 = U_topology(_c, mol)
            _gfd = _fd_gradient(_topo_fn, _c, eps=1.0e-6)
            _dd = float(np.max(np.abs(_g0 - _gfd)))
            _rr = _dd / max(1.0e-9, float(np.max(np.abs(_g0))))
            print(f"  {_where:38s} E = {_e0:14.4f}  "
                  f"|g_an| = {np.max(np.abs(_g0)):10.3e}  "
                  f"|g_fd| = {np.max(np.abs(_gfd)):10.3e}  "
                  f"rel = {_rr:9.3e}  {'PASS' if _rr < 1.0e-3 else 'FAIL'}")

    # ------------------------------------------------------------------
    # Tier C: the REAL call path.  detect_fragments() hands U_C_fragment namedtuples it
    # cannot read, which took the whole functional down to the fallback silently; the
    # orbit version must both survive the call and find symmetry without a pattern list.
    # ------------------------------------------------------------------
    print("\nTier C — fragment symmetry on the real call path:")
    from delfin.manta._fragment_archetypes import (detect_fragments,
                                                   detect_fragment_orbits)
    _fc = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1C(=O)[O-]"))
    AllChem.EmbedMolecule(_fc, randomSeed=0xC0DE)
    _fcc = np.array([list(_fc.GetConformer().GetAtomPosition(i))
                     for i in range(_fc.GetNumAtoms())], dtype=np.float64)
    try:
        U_C_fragment(_fcc, _fc, detect_fragments(_fc))
        print("  SMARTS path: returned (expected it to raise)")
    except Exception as _fx:
        print(f"  SMARTS path: {type(_fx).__name__} — "
              f"this is the silent demotion to _fallback_U_total")
    _orb = detect_fragment_orbits(_fc, _fcc)
    _oe, _og = U_C_fragment(_fcc, _fc, _orb)
    _nops = sum(len(f["operations"]) for f in _orb)
    _loose = sum(len(f["operations"])
                 for f in detect_fragment_orbits(_fc, _fcc, rms_tol=1.0e9))
    print(f"  orbit path : survives the call, E = {_oe:.4f}")
    print(f"  benzoate   : {_loose} graph automorphism(s), {_nops} accepted   "
          f"{'OK — the ring flip is NOT a geometric symmetry here (=O vs [O-] are '
             'distinct and off-axis), and the residual guard rejects it'
             if _nops == 0 else 'check: something asymmetric was accepted'}")
    # Positive control: a molecule that really carries the symmetry must keep it.
    for _smi, _what in (("c1ccccc1", "benzene"),
                        ("c1ccc(-c2ccccn2)nc1", "2,2'-bipyridine")):
        _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
        AllChem.EmbedMolecule(_m, randomSeed=0xC0DE)
        _mc = np.array([list(_m.GetConformer().GetAtomPosition(i))
                        for i in range(_m.GetNumAtoms())], dtype=np.float64)
        _o = detect_fragment_orbits(_m, _mc)
        _all = sum(len(f["operations"])
                   for f in detect_fragment_orbits(_m, _mc, rms_tol=1.0e9))
        _kept = sum(len(f["operations"]) for f in _o)
        print(f"  {_what:15s}: {_all:2d} automorphisms, {_kept:2d} geometric   "
              f"{'OK' if _kept > 0 else 'FAIL — real symmetry lost'}")

    # ------------------------------------------------------------------
    # U_torsion.  The dihedral gradient is the only new calculus in the functional, so
    # it gets its own finite-difference check on a real twisted conjugated system, plus
    # a negative control: a saturated ring must stay untouched.
    # ------------------------------------------------------------------
    print("\nU_torsion — conjugation planarity:")
    _bp = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1-c1ccccc1"))
    AllChem.EmbedMolecule(_bp, randomSeed=0xBEEF)
    _bc = np.array([list(_bp.GetConformer().GetAtomPosition(i))
                    for i in range(_bp.GetNumAtoms())], dtype=np.float64)
    _tt = _build_torsion_targets(_bc, _bp)
    _te, _tg = U_torsion(_bc, _bp, k_torsion=50.0, targets=_tt)
    print(f"  biphenyl:    {len(_tt):3d} restrained dihedrals   "
          f"E = {_te:10.4f}   |g|_max = {np.max(np.abs(_tg)):.4e}")

    def _tors_fn(c):
        return U_torsion(c, _bp, k_torsion=50.0, targets=_tt)

    _tg_fd = _fd_gradient(_tors_fn, _bc, eps=1.0e-6)
    _td = float(np.max(np.abs(_tg - _tg_fd)))
    _tr = _td / max(1.0e-9, float(np.max(np.abs(_tg))))
    print(f"  FD gradient: max|an-fd| = {_td:.4e}   relative = {_tr:.4e}   "
          f"{'PASS' if _tr < 1.0e-4 else 'FAIL'}")

    # Tolerance branch: a near-planar system must cost NOTHING, and a strongly twisted one
    # must still bite -- and its gradient must survive FD.  E = 0 alone proves nothing (it
    # is also what a dead code path returns), so twist the molecule and test there.
    _bc2 = _bc.copy()
    _ax = _bc2[:, 0] > float(np.median(_bc2[:, 0]))       # rotate one half about x
    _ang = math.radians(55.0)
    _ca, _sa = math.cos(_ang), math.sin(_ang)
    _bc2[_ax] = np.column_stack([
        _bc2[_ax][:, 0],
        _ca * _bc2[_ax][:, 1] - _sa * _bc2[_ax][:, 2],
        _sa * _bc2[_ax][:, 1] + _ca * _bc2[_ax][:, 2]])
    for _lbl, _tol in (("strict (planarity forced)", "0"), ("tolerance band", "1")):
        os.environ["DELFIN_FFREE_TORSION_TOL"] = _tol
        _e2, _g2 = U_torsion(_bc2, _bp, k_torsion=50.0, targets=_tt)

        def _tf(c, _t=_tt):
            return U_torsion(c, _bp, k_torsion=50.0, targets=_t)

        _fd2 = _fd_gradient(_tf, _bc2, eps=1.0e-6)
        _r2 = float(np.max(np.abs(_g2 - _fd2))) / max(1.0e-9,
                                                      float(np.max(np.abs(_g2))))
        print(f"  twisted 55deg, {_lbl:<26} E = {_e2:9.3f}   "
              f"FD rel = {_r2:.3e}   {'PASS' if _r2 < 1.0e-4 else 'FAIL'}")
    os.environ.pop("DELFIN_FFREE_TORSION_TOL", None)

    _cy = Chem.AddHs(Chem.MolFromSmiles("C1CCCCC1"))
    AllChem.EmbedMolecule(_cy, randomSeed=0xBEEF)
    _cc = np.array([list(_cy.GetConformer().GetAtomPosition(i))
                    for i in range(_cy.GetNumAtoms())], dtype=np.float64)
    _cn = len(_build_torsion_targets(_cc, _cy))
    print(f"  cyclohexane: {_cn:3d} restrained dihedrals   "
          f"{'OK — pucker untouched' if _cn == 0 else 'FAIL — would flatten the ring'}")

    # ------------------------------------------------------------------
    # Signature-angle cases.  These are the bonding situations the RDKit
    # hybridization label gets wrong; each names an ATOM and its BONDS, never a
    # functional group.  Run unconditionally so a regression shows up without a flag.
    # ------------------------------------------------------------------
    print("\nSignature angle (_signature_theta) — measured targets:")
    _cases = [
        # smiles, triple kind ("MDX" = metal-donor-X, "XDX" = X-donor-X), want deg, what
        ("[Cu][N+]1=CC=CC=C1", "MDX", 120.75,
         "pyridine N->M: planar 6-ring donor, metal IS the exocyclic substituent"),
        ("[Cu][N+]1=CC=CC=C1", "XDX", 118.0,
         "same N, ring-internal C-N-C: the 6-ring's measured internal angle"),
        ("C=C[N]([Cu])C=C", "MDX", 120.0,
         "amido R2N-M conjugated: trigonal planar, METAL IN THE PI PLANE (AFOMEO)"),
        ("C=C[N]([Cu])C=C", "XDX", 120.0,
         "same N, C-N-C: planar too -- one signature, every triple at that atom"),
        ("CN(C)[Cu]", "XDX", 109.5,
         "saturated amine C-N-C: period-2 rule, no pi neighbour"),
        ("CS(C)[Cu]", "XDX", 104.0,
         "period-3 donor, C-S-C closes to the period value (NOT 109.47)"),
        ("CS(C)[Cu]", "MDX", 114.28,
         "same S, M-S-C opens as the substituents close (measured slope 0.87)"),
        ("C[Se](C)[Cu]", "XDX", 101.5,
         "period-4 donor Se: the row keeps closing -- universal, not per-element"),
    ]
    for smi_c, _kind, expect, what in _cases:
        m = Chem.MolFromSmiles(smi_c)
        if m is None:
            print(f"  SKIP (unparsable): {smi_c}")
            continue
        m = Chem.AddHs(m)
        _mets = _smiles_converter_metals()
        # Pick the donor: the heavy atom bonded to the metal.
        _j = None
        for a in m.GetAtoms():
            if a.GetSymbol() in _mets:
                for nb in a.GetNeighbors():
                    if nb.GetSymbol() != "H":
                        _j = int(nb.GetIdx())
                        break
            if _j is not None:
                break
        if _j is None:
            print(f"  SKIP (no donor): {smi_c}")
            continue
        _nbrs = [int(n.GetIdx()) for n in m.GetAtomWithIdx(_j).GetNeighbors()
                 if m.GetAtomWithIdx(int(n.GetIdx())).GetSymbol() != "H"]
        _mi = [x for x in _nbrs if m.GetAtomWithIdx(x).GetSymbol() in _mets]
        _xi = [x for x in _nbrs if x not in _mi]
        if _kind == "MDX" and _mi and _xi:
            _i, _k = _mi[0], _xi[0]
        elif _kind == "XDX" and len(_xi) >= 2:
            _i, _k = _xi[0], _xi[1]
        else:
            print(f"  SKIP (no {_kind} triple): {smi_c}")
            continue
        got = math.degrees(_signature_theta(m, _i, _j, _k))
        ok = abs(got - expect) < 1.0
        print(f"  {'OK  ' if ok else 'FAIL'} {got:6.2f}deg (want {expect:6.2f})  {what}")

    # ----- inverse-variance weighting: does it move, and is its gradient right? -----
    #
    # Both halves are needed.  A lever that is byte-identical when off AND byte-identical
    # when on is dead code that looks alive; a lever that changes the energy but not the
    # gradient is worse -- it silently minimises a different function than it reports.
    # The first table for this test put every bond INSIDE its band, so the weight
    # multiplied zero and the term looked inert.  The bands below sit far outside, and
    # deliberately differ in WIDTH by a factor ~15 so the weights actually diverge.
    print("\nU_bond — inverse-variance weighting:")
    import tempfile as _tf
    with _tf.NamedTemporaryFile("w", suffix=".tsv", delete=False) as _fh:
        _tbl = _fh.name
        _fh.write("L4\tH|N|1.0\t9999\t0.50\t0.51\t0.52\n")   # sigma = 0.0078
        _fh.write("L4\tH|O|1.0\t9999\t0.40\t0.55\t0.70\n")   # sigma = 0.117
    _prev = {k: os.environ.get(k) for k in
             ("DELFIN_FFREE_BOND_BAND", "DELFIN_FFREE_BOND_BANDS", "DELFIN_FFREE_INVVAR")}
    try:
        os.environ["DELFIN_FFREE_BOND_BAND"] = "1"
        os.environ["DELFIN_FFREE_BOND_BANDS"] = _tbl
        os.environ.pop("DELFIN_FFREE_INVVAR", None)
        _e_off, _g_off = U_bond(coords, mol)
        os.environ["DELFIN_FFREE_INVVAR"] = "1"
        _e_on, _g_on = U_bond(coords, mol)
        _gfd = _fd_gradient(lambda c: U_bond(c, mol), coords, eps=1.0e-6)
        _d = float(np.max(np.abs(_g_on - _gfd)))
        _r = _d / max(1.0e-9, float(np.max(np.abs(_g_on))))
        _moved = abs(_e_on - _e_off) > 1.0e-6
        _gdiff = float(np.max(np.abs(_g_on - _g_off)))
        print(f"  E off = {_e_off:12.4f}   E on = {_e_on:12.4f}   "
              f"changed = {_moved}")
        print(f"  max|g_on - g_off| = {_gdiff:.4e}   (0 would mean the weight never "
              f"reached the gradient)")
        print(f"  FD on weighted U_bond: rel = {_r:.3e}   "
              f"{'PASS' if (_r < 1.0e-3 and _moved and _gdiff > 1.0e-9) else 'FAIL'}")
    finally:
        for _k, _v in _prev.items():
            if _v is None:
                os.environ.pop(_k, None)
            else:
                os.environ[_k] = _v

    # ----- U_angle flat bottom: the measured WIDTH, not just the measured centre -----
    # Key ladder level L4 is "{element},{degree}#{outer}#{outer}".  For [Cu](N)(N)O + Hs
    # the central atoms are N (degree 3) and O (degree 2); metal-centred triples are
    # skipped by the term itself.  The two bands differ in width by ~14x and both sit far
    # from the actual ~109 deg, so every angle is outside and the weights diverge.
    print("\nU_angle — measured band as a flat bottom:")
    with _tf.NamedTemporaryFile("w", suffix=".tsv", delete=False) as _fh:
        _atbl = _fh.name
        _fh.write("L4\tN,3#Cu#H\t9999\t10.0\t12.0\t14.0\n")     # sigma = 1.56 deg
        _fh.write("L4\tN,3#H#H\t9999\t10.0\t12.0\t14.0\n")
        _fh.write("L4\tO,2#Cu#H\t9999\t150.0\t160.0\t179.0\n")  # sigma = 11.3 deg
    _prev = {k: os.environ.get(k) for k in
             ("DELFIN_FFREE_ANGLE_BAND", "DELFIN_FFREE_ANGLE_BANDS", "DELFIN_FFREE_INVVAR")}
    try:
        os.environ.pop("DELFIN_FFREE_ANGLE_BAND", None)
        os.environ.pop("DELFIN_FFREE_INVVAR", None)
        os.environ["DELFIN_FFREE_ANGLE_BANDS"] = _atbl
        _ea_plain, _ = U_angle(coords, mol)
        os.environ["DELFIN_FFREE_ANGLE_BAND"] = "1"
        _ea_band, _ga_band = U_angle(coords, mol)
        os.environ["DELFIN_FFREE_INVVAR"] = "1"
        _ea_w, _ga_w = U_angle(coords, mol)
        _gfd = _fd_gradient(lambda c: U_angle(c, mol), coords, eps=1.0e-6)
        _r = (float(np.max(np.abs(_ga_w - _gfd)))
              / max(1.0e-9, float(np.max(np.abs(_ga_w)))))
        _ok = (abs(_ea_band - _ea_plain) > 1.0e-6      # band changed the shape
               and abs(_ea_w - _ea_band) > 1.0e-6      # weight changed it again
               and _r < 1.0e-3)                        # and the gradient followed
        print(f"  E harmonic = {_ea_plain:10.4f}   E band = {_ea_band:10.4f}   "
              f"E band+invvar = {_ea_w:10.4f}")
        print(f"  FD on weighted U_angle: rel = {_r:.3e}   "
              f"{'PASS' if _ok else 'FAIL'}")
    finally:
        for _k, _v in _prev.items():
            if _v is None:
                os.environ.pop(_k, None)
            else:
                os.environ[_k] = _v
