"""Reference coordination polyhedra (unit vectors) + covalent radii for the
metal-FF-free builder (delfin.manta).

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
    # CN3 (iter-32g, User 2026-06-19 ATENET: V/N3O complex missing the
    # trigonal-PYRAMIDAL CN3 isomer — 3 donors forming the base of a pyramid with
    # the metal at the apex ABOVE the donor plane, i.e. the "vacant tetrahedron" /
    # NH3 lone-pair geometry, distinct from the planar 120deg SP-3 and the 90/180deg
    # T-3).  C3v: C3 axis along +z, donors at azimuth 0/120/240, polar angle 68.16deg
    # from the -z axis -> all three M->donor vectors point into ONE hemisphere (z<0),
    # giving 107.0deg donor-M-donor angles (ammonia-like; between the 109.47deg ideal
    # vacant-tetrahedron and a flatter pyramid).  Index order 0,1,2 = the three base
    # vertices in 120deg-cyclic order so the proper-rotation C3 in
    # polya_isomer_count._trigonal_pyramidal_group (0->1->2->0) is a real geometric
    # symmetry of this vertex set (same contract as SP-3/T-3 above).  ADDITIVE only:
    # never replaces SP-3/T-3; reached solely via the additive enumerator under
    # DELFIN_FFFREE_CN3_PYRAMIDAL=1 (default OFF -> byte-identical).
    _tpy_theta = math.radians(68.158)
    _tpy_z = -math.cos(_tpy_theta)
    _tpy_s = math.sin(_tpy_theta)
    R[("CN3", "TPY-3 trigonal pyramidal")] = np.array(
        [[_tpy_s * math.cos(math.radians(120 * k)),
          _tpy_s * math.sin(math.radians(120 * k)), _tpy_z] for k in range(3)],
        float)
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


# Pyykkö & Atsumi (Chem. Eur. J. 2009, 15, 186) single-bond COVALENT radii (Å) for
# the heavy metalloid / post-transition DONOR elements that are ABSENT from ``COV``
# (so the bare ``COV.get(donor, 0.75)`` fallback gives an unphysically short M–D bond,
# e.g. Pd–Sb 2.14 instead of ~2.8 Å).  PUBLISHED open reference values (no CSD/CCDC
# data).  Used ONLY when DELFIN_FFFREE_METALLOID_DONOR=1 (default OFF -> byte-id, this
# table is never consulted).  Sb/As/Te/Se have published values; Bi/Ge/Sn/Pb included
# so every metalloid donor enabled by the resolver has a finite, realistic radius.
_METALLOID_COV = {
    "Sb": 1.40, "As": 1.21, "Bi": 1.51, "Te": 1.38, "Se": 1.16,
    "Ge": 1.20, "Sn": 1.39, "Pb": 1.46,
}


def _donor_cov(donor: str) -> float:
    """Donor covalent radius (Å) for the M–D distance.  Byte-identical default
    (``COV.get(donor, 0.75)``) unless DELFIN_FFFREE_METALLOID_DONOR=1 AND the donor is a
    heavy metalloid MISSING from COV — then use the published Pyykkö radius so the M–D
    bond has a realistic length instead of the 0.75 Å placeholder."""
    if donor in COV:
        return COV[donor]
    if os.environ.get("DELFIN_FFFREE_METALLOID_DONOR", "0") == "1" and donor in _METALLOID_COV:
        return _METALLOID_COV[donor]
    return COV.get(donor, 0.75)


_MD_BAND_CACHE: dict = {}


_TSV_CACHE: dict = {}


def _load_band_tsv(path: str) -> dict:
    """``key -> (p10, p50, p90)`` from any of the measured band TSVs, cached per path.

    One reader for every table (bonds, angles, metal-donor, D-M-D): they share the format
    ``level<TAB>key<TAB>n<TAB>p10<TAB>p50<TAB>p90`` because they come from one measurement
    pass.  Missing file or unreadable row -> empty, and every caller then keeps its
    historic behaviour rather than guessing.
    """
    if not path:
        return {}
    tbl = _TSV_CACHE.get(path)
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
                        tbl[p[1]] = (float(p[3]), float(p[4]), float(p[5]))
                    except ValueError:
                        continue
    except Exception:
        tbl = {}
    _TSV_CACHE[path] = tbl
    return tbl


def _md_band_table() -> dict:
    """``key -> (p10, p50, p90)`` from ``DELFIN_FFREE_MD_BANDS`` (CSD-derived, not shipped).

    Keys are what the offline measurement emitted: ``M,cn|donor_sig|chelate`` down to the
    bare ``M|donor``.  Unset -> empty -> every caller keeps the covalent-radii sum.
    """
    path = os.environ.get("DELFIN_FFREE_MD_BANDS", "")
    if not path:
        return {}
    tbl = _MD_BAND_CACHE.get(path)
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
                        tbl[p[1]] = (float(p[3]), float(p[4]), float(p[5]))
                    except ValueError:
                        continue
    except Exception:
        tbl = {}
    _MD_BAND_CACHE[path] = tbl
    return tbl


_CURRENT_CN = None


def set_current_cn(cn) -> None:
    """Record the coordination number of the complex currently being built.

    Threading ``cn`` through all twelve md_distance() call sites would touch half the
    build path for a value that is one number per complex and constant throughout it.
    One process builds one complex at a time (``_one.py`` pins DELFIN_MAX_PROCESS_WORKERS
    and the pool build runs single-threaded), so a module-level value is well defined here.
    It is only ever READ behind DELFIN_FFREE_MD_MEASURED, and if it is stale or unset the
    lookup simply falls back to the bare element pair -- never a wrong distance, only a
    less specific one.
    """
    global _CURRENT_CN
    try:
        _CURRENT_CN = int(cn) if cn else None
    except Exception:
        _CURRENT_CN = None


def _measured_md(metal: str, donor: str, cn=None, chelate=None):
    """Measured p50 for this metal-donor contact, most specific key first.

    The coordination number is what actually moves an M-D distance, and the measurement
    shows it cleanly: Cd-N runs 2.283 / 2.342 / 2.357 at CN 4 / 5 / 6, Cd-I 2.712 / 2.743 /
    2.837, monotonic for every donor type.  A CN-free p50 would flatten all of that into
    one number -- and worse, for Cu(II) it would place the Jahn-Teller axial bonds 0.47 A
    too short, because the element-pair band is 0.535 A wide and strongly asymmetric.

    Falls back down the ladder and finally to None, in which case the caller keeps the
    historic radii sum -- so a missing bin is never a wrong answer, only an unimproved one.
    """
    tbl = _md_band_table()
    if not tbl:
        return None
    keys = []
    if cn is not None:
        if chelate is not None:
            keys.append(f"{metal},{int(cn)}|{donor}|{int(bool(chelate))}")
        keys.append(f"{metal},{int(cn)}|{donor}")
    keys.append(f"{metal}|{donor}")
    for k in keys:
        row = tbl.get(k)
        if row is not None:
            return float(row[1])                                # p50
    return None


def md_distance(metal: str, donor: str, atom=None, mol=None, cn=None) -> float:
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
    # MEASURED M-D LENGTH (DELFIN_FFREE_MD_MEASURED=1, default OFF -> byte-identical).
    # A covalent-radii sum is not where crystals put a metal-donor bond, and this is a
    # SETTING, not an optimisation: the builder places the donor at this distance and
    # that is the end of it -- no weight, no barrier, no gate that can reject it.
    # Measured basis: over 107 crystals the same radii-sum reference made U_topology score
    # our own frames FIVE TIMES BETTER than reality (ratio 0.20); swapping it for the
    # measured band took the crystal force from 52833 to 67.39, a factor of 784.  Here the
    # same number is used one step earlier, where it costs nothing to be right.
    if os.environ.get("DELFIN_FFREE_MD_MEASURED", "0") == "1":
        _m = _measured_md(metal, donor, cn=(cn if cn else _CURRENT_CN))
        if _m is not None:
            return float(min(4.0, max(0.8, _m)))
    if os.environ.get("DELFIN_FFFREE_MD_CONTEXT", "0") != "1" or atom is None:
        return COV.get(metal, 1.5) + _donor_cov(donor)
    try:
        sym = atom.GetSymbol()
        # Base = historic single-bond covalent radius -> neutral dative donors
        # NEVER regress; only the short anionic σ-donor classes shorten.  (PYYKKO
        # bond-order radii + _pyykko_radius/_max_heavy_bond_order are retained above
        # as reference for a future explicit M=O / M≡N multiple-bond layer.)
        r_d = _donor_cov(sym) - _sigma_shorten(atom, mol)
        md = COV.get(metal, 1.5) + r_d
        if not math.isfinite(md):
            md = COV.get(metal, 1.5) + _donor_cov(donor)
    except Exception:
        md = COV.get(metal, 1.5) + _donor_cov(donor)
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
