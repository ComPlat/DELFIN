"""Reference coordination polyhedra (unit vectors) + covalent radii for the
metal-FF-free builder (delfin.fffree).

Self-contained: ideal CN4/5/6 polyhedron vertex sets (tetrahedron, square planar,
trigonal bipyramid, square pyramid, octahedron, trigonal prism) and a covalent-radii
table.  Metal-donor distances use covalent-radii sums.
"""
from __future__ import annotations
import math
import os
import numpy as np


def _norm_rows(V: np.ndarray) -> np.ndarray:
    return V / np.linalg.norm(V, axis=1, keepdims=True)


def _ref_polyhedra():
    R = {}
    t = 1 / math.sqrt(3)
    # CN2 (iter-32f, DELFIN_FFFREE_CN_EXTEND): linear two-coordinate (180°), the
    # canonical d10 geometry for Cu(I)/Ag(I)/Au(I)/Hg(II) etc.  Two antipodal
    # vertices on the z-axis — the smallest non-trivial polyhedron.  Index order
    # 0=+z, 1=-z so the C2/inversion swap (0<->1) in _linear_group below is a real
    # geometric symmetry of this vertex set (same contract as the higher CNs).
    R[("CN2", "L-2 linear")] = np.array(
        [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]], float)
    # CN3 (iter-32c, User 2026-05-28 ADUMOD: Pd CN3 built as Td=109.5°/linear=180°
    # instead of correct SP-3 trigonal-planar 120° or d8 T-shape 90°/180°).
    R[("CN3", "SP-3 trigonal planar")] = np.array(
        [[1.0, 0.0, 0.0],
         [-0.5, math.sqrt(3) / 2, 0.0],
         [-0.5, -math.sqrt(3) / 2, 0.0]])
    R[("CN3", "T-3 T-shape")] = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    R[("CN4", "T-4 tetrahedron")] = np.array(
        [[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]])
    R[("CN4", "SP-4 square planar")] = np.array(
        [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
    R[("CN5", "TBP-5 trigonal bipyramid")] = np.array(
        [[0, 0, 1], [0, 0, -1], [1, 0, 0],
         [-0.5, math.sqrt(3) / 2, 0], [-0.5, -math.sqrt(3) / 2, 0]])
    # basal vertices in 90deg-cyclic order (45,135,225,315) so the proper-rotation
    # C4 in polya_isomer_count._spy_group (1->2->3->4) is a real geometric symmetry of
    # this vertex set -> chelate cis-edge enumeration & isomer dedup are consistent with
    # placement (the index space here IS the one assemble_from_config places into).
    R[("CN5", "SPY-5 square pyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [1, 1, 0.2], [-1, 1, 0.2], [-1, -1, 0.2], [1, -1, 0.2]], float))
    R[("CN6", "OC-6 octahedron")] = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    R[("CN6", "TPR-6 trigonal prism")] = _norm_rows(np.array(
        [[1, 0, 0.7], [-0.5, math.sqrt(3) / 2, 0.7], [-0.5, -math.sqrt(3) / 2, 0.7],
         [1, 0, -0.7], [-0.5, math.sqrt(3) / 2, -0.7], [-0.5, -math.sqrt(3) / 2, -0.7]]))
    # --- High-CN polyhedra (CN7-9). Vertex INDEX ORDER is chosen to match the
    # proper-rotation generators in polya_isomer_count (_pentagonal_bipyramid_group /
    # _square_antiprism_group / _tricapped_trigonal_prism_group) so isomer dedup is a
    # real geometric symmetry of this vertex set (same contract as SPY-5/TPR-6 above).
    # CN7 PB: idx 0,1 = axial (+z,-z); idx 2-6 = equatorial regular pentagon (0,72,...,288 deg).
    _pent = [[math.cos(2 * math.pi * k / 5), math.sin(2 * math.pi * k / 5), 0.0] for k in range(5)]
    R[("CN7", "PB-7 pentagonal bipyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [0, 0, -1]] + _pent, float))
    # CN8 square antiprism: idx 0-3 = top square (0,90,180,270 deg, +z); idx 4-7 = bottom
    # square (45,135,225,315 deg, -z) -- the 45deg stagger that defines the antiprism.
    _h8 = 0.62
    _top = [[math.cos(math.pi * k / 2), math.sin(math.pi * k / 2), _h8] for k in range(4)]
    _bot = [[math.cos(math.pi * (k + 0.5) / 2 + 0.0), math.sin(math.pi * (k + 0.5) / 2), -_h8] for k in range(4)]
    # bottom at 45,135,225,315: angle = 45 + 90k
    _bot = [[math.cos(math.radians(45 + 90 * k)), math.sin(math.radians(45 + 90 * k)), -_h8] for k in range(4)]
    R[("CN8", "SQAP-8 square antiprism")] = _norm_rows(np.array(_top + _bot, float))
    # CN9 tricapped trigonal prism: idx 0-2 = top triangle (0,120,240 deg, +z); idx 3-5 =
    # bottom triangle (eclipsed, -z); idx 6-8 = caps on the 3 rectangular faces (60,180,300 deg, z=0).
    _tri_t = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), 0.7] for k in range(3)]
    _tri_b = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), -0.7] for k in range(3)]
    _caps = [[1.3 * math.cos(math.radians(60 + 120 * k)), 1.3 * math.sin(math.radians(60 + 120 * k)), 0.0] for k in range(3)]
    R[("CN9", "TTP-9 tricapped trigonal prism")] = _norm_rows(np.array(_tri_t + _tri_b + _caps, float))
    return {k: _norm_rows(v) for k, v in R.items()}


REFS = _ref_polyhedra()

GEOM_BY_CN = {
    2: ["L-2 linear"],
    3: ["SP-3 trigonal planar", "T-3 T-shape"],
    4: ["T-4 tetrahedron", "SP-4 square planar"],
    5: ["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"],
    6: ["OC-6 octahedron", "TPR-6 trigonal prism"],
    7: ["PB-7 pentagonal bipyramid"],
    8: ["SQAP-8 square antiprism"],
    9: ["TTP-9 tricapped trigonal prism"],
}

# covalent radii (metal subset + donors); M-D = r(M) + r(D)
COV = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39, "Se": 1.20, "As": 1.19,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "La": 2.07,
}


def ref_vectors(geometry: str) -> np.ndarray:
    for (cn, shape), v in REFS.items():
        if shape == geometry:
            return v
    raise KeyError(geometry)


# --- Context-aware donor radii (env-gated; #305 root fix) -------------------
# Pyykkö & Atsumi single / double / triple-bond COVALENT radii (Å), OPEN data:
#   single  : Pyykkö & Atsumi, Chem. Eur. J. 2009, 15, 186  (table of r_cov^(1))
#   double  : Pyykkö & Atsumi, Chem. Eur. J. 2009, 15, 12770 (r_cov^(2))
#   triple  : Pyykkö, Riedel & Patzschke, Chem. Eur. J. 2005, 11, 3511 (r_cov^(3))
# These are PUBLISHED open reference values (no CCDC/CSD data).  The donor atom's
# effective radius is selected by its MAX heavy-neighbour bond order (single /
# aromatic≈double / double / triple).  The element-pair covalent sum (default,
# flag-OFF) is left BYTE-IDENTICAL.
PYYKKO_SINGLE = {
    "C": 0.75, "N": 0.71, "O": 0.63, "S": 1.03, "P": 1.11, "Se": 1.16,
    "F": 0.64, "Cl": 0.99, "Br": 1.14, "I": 1.33, "As": 1.21, "Te": 1.36,
}
PYYKKO_DOUBLE = {
    "C": 0.67, "N": 0.60, "O": 0.57, "S": 0.94, "P": 1.02, "Se": 1.07,
    "F": 0.59, "Cl": 0.95, "Br": 1.09, "I": 1.29, "As": 1.14, "Te": 1.30,
}
PYYKKO_TRIPLE = {
    "C": 0.60, "N": 0.54, "O": 0.53, "S": 0.95, "P": 0.94, "Se": 1.07,
    "As": 1.06, "Te": 1.21,
}

# Anionic / short-σ ligand-class shortening (Å), applied to the donor radius.
# A documented chemical rule (anionic σ-donors form notably shorter M–D bonds
# than neutral L donors), keyed ONLY on the local molecular graph — never on a
# specific SMILES.  Magnitudes are modest and physically reasonable; the
# classes (azide-N, cyanide-C/N, halide, hydroxide/alkoxide/oxo, amide/imide-N)
# are the canonical short anionic donors.  Open chemical knowledge, no CSD data.
_SIGMA_SHORTEN = {
    "azide": 0.12,      # terminal N of N–N=N (anionic, short Co–N≈1.94 vs py 2.16)
    "cyanide": 0.10,    # C/N of C≡N (strong-field σ-donor, short M–C/M–N)
    "halide": 0.08,     # F/Cl/Br/I anionic terminal halide
    "oxo_alkoxo": 0.12, # O with ≤1 heavy neighbour (oxo / hydroxide / alkoxide)
    "amide": 0.10,      # deprotonated (anionic) N donor
}


def _max_heavy_bond_order(atom) -> float:
    """Max bond order from ``atom`` to its HEAVY neighbours (H ignored).
    Aromatic bonds count as ~1.5.  Graph-only, deterministic."""
    from rdkit.Chem import BondType
    bo = 1.0
    for b in atom.GetBonds():
        other = b.GetOtherAtom(atom)
        if other.GetAtomicNum() == 1:
            continue
        bt = b.GetBondType()
        if b.GetIsAromatic() or bt == BondType.AROMATIC:
            v = 1.5
        elif bt == BondType.DOUBLE:
            v = 2.0
        elif bt == BondType.TRIPLE:
            v = 3.0
        else:
            v = 1.0
        if v > bo:
            bo = v
    return bo


def _pyykko_radius(sym: str, bond_order: float) -> float:
    """Bond-order-selected Pyykkö covalent radius (Å); fall back to COV[sym]."""
    if bond_order >= 2.5 and sym in PYYKKO_TRIPLE:
        return PYYKKO_TRIPLE[sym]
    if bond_order >= 1.5 and sym in PYYKKO_DOUBLE:     # aromatic≈1.5 -> double
        return PYYKKO_DOUBLE[sym]
    if sym in PYYKKO_SINGLE:
        return PYYKKO_SINGLE[sym]
    return COV.get(sym, 0.75)


def _sigma_shorten(atom, mol) -> float:
    """Detect a short anionic σ-donor class from the LOCAL graph around the donor
    ``atom`` and return the radius shortening (Å, >=0).  Graph-only, universal,
    deterministic.  Returns 0.0 for ordinary neutral L donors."""
    sym = atom.GetSymbol()
    heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
    # Halide donor: a lone / terminal F,Cl,Br,I bound to the metal.
    if sym in ("F", "Cl", "Br", "I"):
        return _SIGMA_SHORTEN["halide"]
    # Cyanide: donor C or N that is triple-bonded to a C/N partner (–C≡N / N≡C–).
    if sym in ("C", "N"):
        from rdkit.Chem import BondType
        for b in atom.GetBonds():
            o = b.GetOtherAtom(atom)
            if b.GetBondType() == BondType.TRIPLE and o.GetSymbol() in ("C", "N"):
                return _SIGMA_SHORTEN["cyanide"]
    # Azide terminal N: donor N bonded to an N that is itself bonded to a 3rd N
    # (the N–N=N chain), i.e. the coordinating end of an azide.
    if sym == "N":
        for nb in heavy:
            if nb.GetSymbol() == "N":
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() != atom.GetIdx() and nb2.GetSymbol() == "N":
                        return _SIGMA_SHORTEN["azide"]
    # Oxo / hydroxide / alkoxide: O donor with <=1 heavy neighbour.
    if sym == "O" and len(heavy) <= 1:
        return _SIGMA_SHORTEN["oxo_alkoxo"]
    # Amide / imide N: a deprotonated (formally anionic) N donor.
    if sym == "N" and atom.GetFormalCharge() < 0:
        return _SIGMA_SHORTEN["amide"]
    return 0.0


def md_distance(metal: str, donor: str, atom=None, mol=None) -> float:
    """Metal–donor placement distance (Å).

    DEFAULT (flag OFF or no donor context): the original element-pair covalent
    sum ``COV[metal] + COV[donor]`` — byte-identical to the historic behaviour.

    CONTEXT-AWARE (``DELFIN_FFFREE_MD_CONTEXT=1`` and ``atom`` given): NEUTRAL
    dative donors (pyridine / amine / imine / aromatic-N / ether / phosphine —
    the most common class) keep the HISTORIC single-bond covalent radius and so
    are byte-identical to the element-pair sum (a dative bond from a lone pair has
    single-bond length REGARDLESS of the donor's internal aromaticity/unsaturation
    — keying the M–D radius off the donor's own bond order erroneously shortened
    every aromatic/imine N donor).  Differentiation comes ONLY from shortening the
    genuinely-short anionic σ-donor classes detected from the local graph (azide /
    cyanide / halide / oxo-alkoxo / amide) — this is what distinguishes e.g.
    azide-N (short) from pyridine-N (unchanged) that the bare element sum cannot
    (#305 / GIXFIF).  Universal (graph-only, never SMILES-specific), deterministic."""
    if os.environ.get("DELFIN_FFFREE_MD_CONTEXT", "0") != "1" or atom is None:
        return COV.get(metal, 1.5) + COV.get(donor, 0.75)
    try:
        sym = atom.GetSymbol()
        # Base = historic single-bond covalent radius -> neutral dative donors
        # NEVER regress; only the short anionic σ-donor classes shorten.  (PYYKKO
        # bond-order radii + _pyykko_radius/_max_heavy_bond_order are retained above
        # as reference for a future explicit M=O / M≡N multiple-bond layer.)
        r_d = COV.get(sym, 0.75) - _sigma_shorten(atom, mol)
        md = COV.get(metal, 1.5) + r_d
        if not math.isfinite(md):
            md = COV.get(metal, 1.5) + COV.get(donor, 0.75)
    except Exception:
        md = COV.get(metal, 1.5) + COV.get(donor, 0.75)
    # Final guard: never non-finite or absurd.
    return float(min(4.0, max(0.8, md)))


def _kabsch_resid(P: np.ndarray, Q: np.ndarray) -> float:
    """Min mean-squared deviation aligning P onto Q by a proper rotation
    (Kabsch, determinant-corrected to forbid reflection)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = P @ R.T - Q
    return float((diff * diff).sum() / len(P))


def cshm(observed_vecs, geometry: str) -> float:
    """Continuous shape measure (0 = ideal, larger = worse) of observed donor
    unit-vectors against the ideal ``geometry`` polyhedron — scale-, rotation-
    and permutation-invariant (Kabsch over all vertex permutations, S =
    100·min_resid/obs_var).  Self-contained, deterministic; returns 0.0 for
    degenerate / size-mismatched input.  Used by the self-gate to reject
    catastrophic-coordination-shape outlier builds (#39)."""
    import itertools
    try:
        Q = ref_vectors(geometry)
    except KeyError:
        return 0.0
    P = np.asarray(observed_vecs, dtype=float)
    n = len(P)
    if n < 2 or n != len(Q):
        return 0.0

    def _unit_rms(a):
        rms = math.sqrt(float((a * a).sum()) / len(a))
        return a / rms if rms > 1e-9 else a

    P = _unit_rms(P)
    Q = _unit_rms(Q)
    obs_var = float((P * P).sum() / n)
    if obs_var < 1e-9:
        return 0.0
    best = float("inf")
    for perm in itertools.permutations(range(n)):
        r = _kabsch_resid(P, Q[list(perm)])
        if r < best:
            best = r
    return 100.0 * best / obs_var
