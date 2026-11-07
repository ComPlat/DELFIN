"""Metal-donor target lengths from the bond-valence model (DELFIN_FFFREE_ML_TARGETS_BV, default OFF).

WHY.  Every metal-donor distance MANTA builds comes from one number per element pair: the
FF-free seating places a donor at the Cordero covalent-radius sum (polyhedra.md_distance), the
legacy path pins it to a hand-set table (ml_tables._get_ml_bond_length).  Neither knows the
oxidation state, the coordination number or what kind of donor it is, and those are what move
an M-D bond in a crystal: Co(III)-N (low spin) sits near 1.95 A, Co(II)-N (high spin) near
2.15 A, Cd(II)-N near 2.35 A.  Measured against the crystal structures, the best clean frame
of a typical system misses its crystal M-D lengths by 0.14-0.15 A (median of the maximum per
system), while two crystals of the same compound agree to 0.05 A in 88 % of the families.

THE MODEL.  Bond valence (I. D. Brown and D. Altermatt, Acta Cryst. B41 (1985) 244; N. E.
Brese and M. O'Keeffe, Acta Cryst. B47 (1991) 192):

    r_ij = R0_ij - b ln(s_ij),   b = 0.37 A,   sum_j s_ij = V (oxidation state of the metal)

The oxidation state comes from the ligand balance of the SMILES (_oxidation_state._ox_cut, the
same cut that the ionic O-H length and the spin rule use), and it is distributed over the CN
donors by the equal-valence rule, s_ij = V / CN.  R0 is either the literature value of the
element pair (M. O'Keeffe and N. E. Brese, J. Am. Chem. Soc. 113 (1991) 3226, data file
data/bv_okeeffe_brese_1991.tsv) or an effective R0 calibrated on crystal structures.

THE TABLE (one file, chosen by DELFIN_FFFREE_ML_TARGETS_BV_TABLE; unset -> the shipped table
data/ml_targets_cod.tsv, derived from the Crystallography Open Database, CC0; see
PROVENANCE.md).  Rows, most specific first:

    D  metal ox cn element class  n  median ...   direct median M-D at this oxidation state and CN
    R  metal ox *  element class  n  R0_eff ...   effective R0 = median(r + b ln(V/CN)), CN pooled
    R  metal ox *  element *      n  R0_eff ...   the same over all donor classes of the element
    R  metal *  *  element class  n  R0_eff ...   over all oxidation states (V carries them)
    R  metal *  *  element *      n  R0_eff ...

A donor takes the first level whose row exists with n >= DELFIN_FFFREE_ML_TARGETS_BV_MIN_N
(default 30); the literature R0 is the last level and is used only with
DELFIN_FFFREE_ML_TARGETS_BV_LITERATURE=1.  No level -> None, and the caller keeps its historic
value.  Oxidation state 0 has no bond-valence form (V = 0): only its direct rows apply.

DONOR CLASS = environment + charge type, read from the ligand as a closed-shell fragment
without the metal (the cut mol, or a decompose() ligand fragment):

    t0..t3   terminal donor (no heavy neighbour) with 0..3 H: halide/oxo/nitrido, hydroxo, aqua, ammine
    pi       donor in or on a pi system: aromatic (pyridine, imidazole, pyrazolate), imine,
             carbonyl, carboxylate, phenolate, nitrile, CO, singlet carbene
    sp3      saturated donor (amine, ether, alkoxide, thiolate, phosphane)
  + x / l    the ligand is anionic and its negative charge sits on the donor or on the donor's
             conjugated system (resonance-invariant: both acac O, both carboxylate O are x),
             otherwise l

Only mononuclear transition-metal complexes whose oxidation state the cut determines get a
target; hapto complexes, several metals, a nitrosyl or an unreadable valence -> None.

One complex is built at a time per process (as for polyhedra.set_current_cn), so the targets of
the complex being built are kept as module state: set_context(smiles) at the outermost build
call, clear_context() after it.  Default OFF -> set_context does nothing and every query
returns None -> byte-identical.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

from delfin.manta import _elements as _EL

FLAG = "DELFIN_FFFREE_ML_TARGETS_BV"
TABLE_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_TABLE"
MIN_N_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_MIN_N"
LITERATURE_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_LITERATURE"
# Where the targets act: "all" (default) = every caller; "seat" = only the FF-free seating of a
# donor (polyhedra.md_distance with the donor atom), so gate tolerances, legacy pins and every
# per-element question keep their historic value -- separates a better seat from a moved gate.
# "pins" = the seating plus the legacy UFF pins and pre-UFF snaps (ml_tables._pin_ml_bond_length); every
# other per-element question (topology checks, tolerance windows, clash limits) keeps its historic value.
SCOPE_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_SCOPE"
# "lengths" = every length a builder places, pins, snaps or scores against (all _get_ml_bond_length
# callers and the seating), while the topology checks and tolerance windows (ml_tables._check_ml_bond_length,
# polyhedra.md_distance without a donor atom) keep their historic value.  Measured 23.10.2026 on the 500
# bench: "pins" missed the legacy length setters outside uff_constraints/pre_uff_snap (soft-donor UFF,
# Open Babel constraints, post-optimiser, variational refiner, embed strategies) and lost 3.7 points K2.
_PIN_DEPTH = [0]           # > 0 while a pin or snap asks for its target length
_CHECK_DEPTH = [0]         # > 0 while a check asks for its comparison length
HISTORIC_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_HISTORIC"
LEVELS_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_LEVELS"
# The metallacycle of a chelate is embedded with its M-D radii at the targets (bounds pinned, donor-donor
# left to the ligand) instead of RDKit dative lengths; see embed_at_targets.
EMBED_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_EMBED"
ALL_LEVELS = ("N", "D", "R", "Re", "Ro", "Roe")
DEFAULT_LEVELS = ("N", "D", "R", "Re")
B_BV = 0.37
# A complex whose donors of one element get targets further apart than this is not given a
# per-element target: the callers that ask by element only (the legacy pins and the tolerance
# windows) then keep their historic value instead of one length for an oxo and a carboxylate O.
ELEMENT_SPREAD_MAX = 0.10

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
SHIPPED_TABLE = os.path.join(_DATA, "ml_targets_cod.tsv")
_LIT_FILE = os.path.join(_DATA, "bv_okeeffe_brese_1991.tsv")

_TABLE_CACHE: Dict[str, Dict[tuple, Tuple[int, float]]] = {}
_LIT: Optional[Dict[str, Tuple[float, float]]] = None
_CTX: Optional[dict] = None


def bv_targets_enabled() -> bool:
    return os.environ.get(FLAG, "0") == "1"


def embed_at_targets() -> bool:
    """True when a chelate of the complex being built embeds its metallacycle at the target radii.

    The seating (_orient_chelate_to_vertices) moves each donor radially onto its target.  The
    metallacycle is embedded with RDKit's own dative M-D length, close to the covalent sum the
    historic targets were; a bond-valence target 0.15-0.2 A away (Co(III)-O 2.04 -> 1.89, Zn(II)-N
    1.93 -> 2.14) then moves every donor against its own backbone, the donor-backbone bonds
    collapse, the self-gate rejects every isomer and the complex falls to legacy (CHELATE_EMPTY;
    AFEFIB, EGEPIR, FOJVIL, UDOMEG on dev300).  Embedding the ring at the very radii it is seated
    at leaves the rescale nothing to move: the bite comes out of the ring closure."""
    return (bv_targets_enabled() and _CTX is not None
            and os.environ.get(EMBED_ENV, "0") == "1")


def _min_n() -> int:
    try:
        return int(os.environ.get(MIN_N_ENV, "30"))
    except ValueError:
        return 30


# ---------------------------------------------------------------------------------------------
# literature R0
# ---------------------------------------------------------------------------------------------
def _literature_params() -> Dict[str, Tuple[float, float]]:
    global _LIT
    if _LIT is None:
        tbl = {}
        try:
            with open(_LIT_FILE) as fh:
                for ln in fh:
                    if ln.startswith("#") or ln.startswith("element"):
                        continue
                    p = ln.split()
                    if len(p) >= 3:
                        tbl[p[0]] = (float(p[1]), float(p[2]))
        except OSError:
            tbl = {}
        _LIT = tbl
    return _LIT


def r0_literature(metal: str, donor: str) -> Optional[float]:
    """O'Keeffe-Brese R0 of the pair, or None when an element is not tabulated."""
    p = _literature_params()
    if metal not in p or donor not in p:
        return None
    (ri, ci), (rj, cj) = p[metal], p[donor]
    return ri + rj - ri * rj * (math.sqrt(ci) - math.sqrt(cj)) ** 2 / (ci * ri + cj * rj)


# ---------------------------------------------------------------------------------------------
# donor classes
# ---------------------------------------------------------------------------------------------
def _ligand_fragment_atoms(atom) -> List:
    """Atoms of the connected fragment of ``atom`` (metals not crossed)."""
    seen = {atom.GetIdx()}
    stack = [atom]
    out = [atom]
    while stack:
        a = stack.pop()
        for n in a.GetNeighbors():
            if n.GetIdx() in seen or _EL.is_metal(n.GetSymbol()):
                continue
            seen.add(n.GetIdx())
            out.append(n)
            stack.append(n)
    return out


def _conjugated_system(atom) -> List:
    """The donor plus every atom reachable from it over conjugated bonds."""
    seen = {atom.GetIdx()}
    stack = [atom]
    out = [atom]
    while stack:
        a = stack.pop()
        for b in a.GetBonds():
            if not (b.GetIsConjugated() or b.GetIsAromatic()):
                continue
            n = b.GetOtherAtom(a)
            if n.GetIdx() in seen or _EL.is_metal(n.GetSymbol()):
                continue
            seen.add(n.GetIdx())
            out.append(n)
            stack.append(n)
    return out


def donor_class(atom) -> Optional[str]:
    """Class of a donor atom of a closed-shell, metal-free ligand fragment (see module doc).

    ``atom`` belongs to a mol whose aromaticity and conjugation flags are set.  Neighbours that
    are metals are ignored, so the function also reads a donor that still carries its metal bond.
    """
    try:
        if atom.GetAtomicNum() == 1:
            return None                                   # hydride: not a target of this model
        heavy = [n for n in atom.GetNeighbors()
                 if n.GetAtomicNum() > 1 and not _EL.is_metal(n.GetSymbol())]
        n_h = atom.GetTotalNumHs(includeNeighbors=True)
        bonds = [b for b in atom.GetBonds() if not _EL.is_metal(b.GetOtherAtom(atom).GetSymbol())]
        if not heavy:
            env = "t%d" % min(n_h, 3)
            conj = [atom]
        elif (atom.GetIsAromatic()
              or any(b.GetIsConjugated() or b.GetBondTypeAsDouble() > 1.0 for b in bonds)
              or any(n.GetIsAromatic() for n in heavy)
              or (atom.GetSymbol() == "C" and len(heavy) == 2 and n_h == 0
                  and len(bonds) == 2 and all(b.GetBondTypeAsDouble() == 1.0 for b in bonds))):
            # aromatic donor, donor in a conjugated / multiple bond, donor on an aromatic
            # ring, or a singlet carbene (two sigma bonds, no H): the SMILES cut and a ligand
            # perceived from coordinates differ in WHICH of these an atom is (aromatic NHC vs
            # carbene, pyridine vs Kekule imine), not in whether it is one -- one class
            env = "pi"
            conj = _conjugated_system(atom)
        else:
            env = "sp3"
            conj = [atom]
        q_frag = sum(a.GetFormalCharge() for a in _ligand_fragment_atoms(atom))
        neg = any(a.GetFormalCharge() < 0 for a in conj)
        return env + ("x" if (q_frag < 0 and neg) else "l")
    except Exception:
        return None


def _sanitized_copy(mol):
    from rdkit import Chem
    m = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(m)
        return m
    except Exception:
        pass
    m = Chem.Mol(mol)
    try:
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                         | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                         | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
                         catchErrors=True)
    except Exception:
        pass
    return m


def complex_donors(mol, metal_idx: int):
    """(metal, OS, [(donor idx, element, class)]) of a mononuclear complex, or None.

    Oxidation state and classes from the same cut of the SMILES as written (_ox_cut_full).
    None when the cut declines (doubt rule) or a donor cannot be classified.
    """
    from delfin.manta._oxidation_state import _ox_cut_full
    r = _ox_cut_full(mol, metal_idx)
    if r is None:
        return None
    os_, _d, cut, donors = r
    cm = _sanitized_copy(cut)
    out = []
    for i in donors:
        a = cm.GetAtomWithIdx(int(i))
        if a.GetAtomicNum() == 1:
            out.append((int(i), "H", None))
            continue
        c = donor_class(a)
        if c is None:
            return None
        out.append((int(i), a.GetSymbol(), c))
    return mol.GetAtomWithIdx(int(metal_idx)).GetSymbol(), int(os_), out


# ---------------------------------------------------------------------------------------------
# table
# ---------------------------------------------------------------------------------------------
def table_path() -> str:
    return os.environ.get(TABLE_ENV, "") or SHIPPED_TABLE


def load_table(path: str) -> Dict[tuple, Tuple[int, float]]:
    """``(kind, metal, ox, cn, element, class) -> (n, value)``; cn/class '*' where pooled.

    Kinds N (cn column 'cn/nx'), D, R, and P (bond-valence parameters: ('P', name) -> (0, value)).
    """
    tbl = _TABLE_CACHE.get(path)
    if tbl is not None:
        return tbl
    tbl = {}
    try:
        with open(path) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"):
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) < 8 or p[0] not in ("N", "D", "R", "P", "F", "J"):
                    continue
                try:
                    if p[0] == "P":
                        tbl[("P", p[1])] = (0, float(p[7]))
                    else:
                        tbl[(p[0], p[1], p[2], p[3], p[4], p[5])] = (int(p[6]), float(p[7]))
                except ValueError:
                    continue
    except OSError:
        tbl = {}
    _TABLE_CACHE[path] = tbl
    return tbl


def _levels():
    v = os.environ.get(LEVELS_ENV, "")
    return {x.strip() for x in v.split(",") if x.strip()} if v else set(DEFAULT_LEVELS)


def bv_params(table) -> Tuple[float, float]:
    """(b, rho) of the table's R rows; the literature (0.37, equal valence) when it names none."""
    b = table.get(("P", "b"), (0, B_BV))[1]
    rho = table.get(("P", "rho"), (0, 1.0))[1]
    return float(b), float(rho)


def n_x_donors(donors) -> int:
    """Number of X-type donors (class ...x, and hydrides) of a (idx, element, class) list."""
    return sum(1 for _i, el, c in donors if el == "H" or (c is not None and c.endswith("x")))


def donor_valence(ox: int, cn: int, nx: int, cls: Optional[str], rho: float) -> Optional[float]:
    """Bond valence of one donor: V w / sum w, w = rho for an X-type donor, 1 for an L-type one."""
    tot = rho * nx + (cn - nx)
    if ox <= 0 or tot <= 0:
        return None
    return ox * (rho if (cls is not None and cls.endswith("x")) else 1.0) / tot


def donor_target(table, metal: str, ox: int, cn: int, element: str, cls: Optional[str],
                 nx: Optional[int] = None, spin: Optional[str] = None):
    """(target A, level) of one donor, or (None, None).

    Levels, in order (DELFIN_FFFREE_ML_TARGETS_BV_LEVELS selects a subset; default N, D, R, Re):
      N    direct median at (ox, CN, nx, class); nx = number of X-type donors of the complex,
           so [ZnBr4]2- (nx 4) and ZnBr2L2 (nx 2) are different rows
      D    direct median at (ox, CN, class)
      R    effective R0 at (ox, class), r = R0 - b ln(s), s from the X/L valence split (table P rows)
      Re   effective R0 at (ox, element)
      Ro   effective R0 over all ox (class), Roe per element -- measured worse than the
           historic covalent sum (0.12-0.34 A median error), off unless named
      L    literature R0 (only with DELFIN_FFFREE_ML_TARGETS_BV_LITERATURE=1)
    """
    nmin = _min_n()
    use = _levels()
    if spin is not None and cls is not None:
        # the spin state's own rows first (phase 2, DELFIN_FFFREE_ML_TARGETS_BV_SPIN_BOTH)
        for kind, cnk, lev in (("N", f"{cn}/{nx}", "Ns"), ("D", str(cn), "Ds")):
            if kind == "N" and nx is None:
                continue
            row = table.get((kind, metal, str(ox), cnk, element, f"{cls}:{spin}"))
            if row is not None and row[0] >= nmin:
                return row[1], lev
    if "N" in use and nx is not None:
        row = table.get(("N", metal, str(ox), f"{cn}/{nx}", element, cls))
        if row is not None and row[0] >= nmin:
            return row[1], "N"
    if "D" in use:
        row = table.get(("D", metal, str(ox), str(cn), element, cls))
        if row is not None and row[0] >= nmin:
            return row[1], "D"
    if ox <= 0 or cn <= 0:
        return None, None
    b, rho = bv_params(table)
    s = donor_valence(ox, cn, nx if nx is not None else 0, cls, rho if nx is not None else 1.0)
    if s is None:
        return None, None
    lns = b * math.log(s)
    for o, c, lev in ((str(ox), cls, "R"), (str(ox), "*", "Re"), ("*", cls, "Ro"), ("*", "*", "Roe")):
        if lev not in use:
            continue
        row = table.get(("R", metal, o, "*", element, c))
        if row is not None and row[0] >= nmin:
            return row[1] - lns, lev
    if os.environ.get(LITERATURE_ENV, "0") == "1":
        r0 = r0_literature(metal, element)
        if r0 is not None:
            return r0 - B_BV * math.log(ox / float(cn)), "L"
    return None, None


def complex_targets(table, metal: str, ox: int, donors, spin: Optional[str] = None) -> Dict[int, Tuple[float, str]]:
    """``donor idx -> (target, level)`` for the donors that get one."""
    cn = len(donors)
    nx = n_x_donors(donors)
    out = {}
    for i, el, cls in donors:
        if el == "H" or cls is None:
            continue
        t, lev = donor_target(table, metal, ox, cn, el, cls, nx, spin)
        if t is not None:
            out[i] = (float(t), lev)
    return out


# ---------------------------------------------------------------------------------------------
# spin state (phase 2, DELFIN_FFFREE_ML_TARGETS_BV_SPIN_BOTH, default OFF)
# ---------------------------------------------------------------------------------------------
SPIN_BOTH_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_SPIN_BOTH"
# A spin state is "clear" when the crystals of this metal, oxidation state, CN and number of
# strong-field donors are at least this share in it; between the two shares both states are built.
SPIN_CLEAR_SHARE = 0.85
_3D_METALS = frozenset("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
_CPH = frozenset("C P As Sb Si Ge Sn B H".split())


def spin_both_enabled() -> bool:
    return os.environ.get(SPIN_BOTH_ENV, "0") == "1"


def spin_candidate(metal: str, ox: int, cn: int) -> bool:
    """3d metal with d4-d7 at CN >= 5: high and low spin both exist in octahedral fields."""
    from delfin.manta._oxidation_state import _OX_GROUP
    if metal not in _3D_METALS or metal not in _OX_GROUP:
        return False
    d = _OX_GROUP[metal] - int(ox)
    return 4 <= d <= 7 and int(cn) >= 5


def n_strong_field(pairs, n_hydride: int = 0) -> int:
    """Strong-field donors of (element, class) pairs: C/P/H-type and pi-accepting N (spectrochemical
    series, as in the spin rule of the xtb step); hydrides count as strong."""
    n = n_hydride
    for el, c in pairs:
        if el in _CPH or (el == "N" and c is not None and c.startswith("pi")):
            n += 1
    return n


def spin_by_rule(metal: str, ox: int, pairs, n_hydride: int = 0) -> str:
    """Ligand-field rule for an octahedral-type 3d d4-d7 centre (the spin rule of the xtb step,
    branch work/2026-09-29-spin-rule): Co(III) low spin; Mn(II) and Co(II) low spin only with a
    C/P/H-type majority (largest pairing energy / Delta_o ratio); otherwise low spin with a
    strong-field majority (spectrochemical series)."""
    n = len(pairs) + n_hydride
    if metal == "Co" and int(ox) == 3:
        return "ls"
    if metal in ("Mn", "Co") and int(ox) == 2:
        n_cph = n_hydride + sum(1 for el, _c in pairs if el in _CPH)
        return "ls" if 2 * n_cph > n else "hs"
    return "ls" if 2 * n_strong_field(pairs, n_hydride) > n else "hs"


def spins_for(table, metal: str, ox: int, donors) -> List[str]:
    """['ls'], ['hs'], ['ls', 'hs'] (both built) or [] (no spin choice / no measured share)."""
    cn = len(donors)
    if not spin_candidate(metal, ox, cn):
        return []
    ns = n_strong_field([(el, c) for _i, el, c in donors if el != "H"],
                        sum(1 for _i, el, _c in donors if el == "H"))
    row = table.get(("F", metal, str(ox), f"{cn}/{ns}", "*", "*"))
    if row is None or row[0] < _min_n():
        return []
    ls = row[1]
    if ls >= SPIN_CLEAR_SHARE:
        return ["ls"]
    if ls <= 1.0 - SPIN_CLEAR_SHARE:
        return ["hs"]
    return ["ls", "hs"]


# ---------------------------------------------------------------------------------------------
# context of the complex being built
# ---------------------------------------------------------------------------------------------
def _mol_from_smiles(smiles: str):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
    return m


def build_context(smiles: str, spin: Optional[str] = None) -> Optional[dict]:
    """Targets of one complex: ``{"metal", "ox", "cn", "by_class", "by_element", "spins", "spin"}``
    or None.  ``spin`` (ls/hs) selects that spin state's rows where the table has them."""
    from delfin.manta._oxidation_state import _OX_GROUP
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    mi = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP]
    if len(mi) != 1:
        return None
    cd = complex_donors(mol, mi[0])
    if cd is None:
        return None
    metal, ox, donors = cd
    table = load_table(table_path())
    spins = spins_for(table, metal, ox, donors) if spin_both_enabled() else []
    if spin is None and len(spins) == 1:
        spin = spins[0]                                  # a clear spin state: its own rows
    tg = complex_targets(table, metal, ox, donors, spin)
    # A donor the table has no row for keeps a historic length -- but the better of the two: the
    # legacy pair table (measured on held-out crystals 0.063 A median error) instead of the FF-free
    # covalent sum (0.126 A), so that switching a complex onto the FF-free path does not move its
    # unmeasured donors further from the crystal.  Not for heavy metalloid donors, whose legacy
    # values are the known overshoot METALLOID_MD_LEN replaces.  Off until measured in a build
    # (DELFIN_FFFREE_ML_TARGETS_BV_HISTORIC=1 turns it on;
    # default 0).
    if os.environ.get(HISTORIC_ENV, "0") == "1":
        from delfin.manta.ml_tables import _METAL_LIGAND_BOND_LENGTHS, _METALLOID_MD_DONORS
        for i, el, cls in donors:
            if i in tg or el == "H" or cls is None or el in _METALLOID_MD_DONORS:
                continue
            v = _METAL_LIGAND_BOND_LENGTHS.get((metal, el))
            if v is not None:
                tg[i] = (float(v), "H")
    if not tg:
        return None
    by_class: Dict[Tuple[str, str], List[float]] = {}
    for i, el, cls in donors:
        if i in tg:
            by_class.setdefault((el, cls), []).append(tg[i][0])
    by_class_t = {k: sum(v) / len(v) for k, v in by_class.items()}
    by_el: Dict[str, List[float]] = {}
    for (el, _c), v in by_class.items():
        by_el.setdefault(el, []).extend(v)
    n_el: Dict[str, int] = {}
    for _i, el, _cls in donors:
        n_el[el] = n_el.get(el, 0) + 1
    by_element = {}
    for el, v in by_el.items():
        # every donor of the element must have a target, and they must agree
        if len(v) == n_el.get(el, 0) and max(v) - min(v) <= ELEMENT_SPREAD_MAX:
            by_element[el] = sum(v) / len(v)
    return {"metal": metal, "ox": ox, "cn": len(donors), "by_class": by_class_t,
            "by_element": by_element, "levels": sorted({lv for _t, lv in tg.values()}),
            "spins": spins, "spin": spin}


def set_context(smiles, spin: Optional[str] = None) -> None:
    """Record the targets of the complex about to be built (no-op when the switch is off)."""
    global _CTX
    _CTX = None
    if not bv_targets_enabled() or not smiles:
        return
    try:
        _CTX = build_context(str(smiles), spin)
    except Exception:
        _CTX = None


def clear_context() -> None:
    global _CTX
    _CTX = None


def current_context() -> Optional[dict]:
    return _CTX


def _same_environment_target(ctx, donor: str, cls: Optional[str]) -> Optional[float]:
    """Target of the complex's donors of this element in the same environment (t0..t3, pi, sp3),
    whatever their x/l flag, when they agree within ELEMENT_SPREAD_MAX.

    The ligand fragments the builder works on are not the SMILES cut: decompose() hands a
    phenolate or quinolinate O over as a neutral, H-free atom, so the fragment reads 'pil' where
    the complex has 'pix'.  The class lookup then misses, and with a second O class in the
    complex (FOJVIL: methanol O 2.19, quinolinate O 2.06) the element lookup refuses as well --
    the donor kept its covalent-sum length 1.88 next to a quinolinate N at its table length 2.14.
    The environment is read from the bonds, which the fragment keeps; only the charge is lost."""
    if not cls:
        return None
    env = cls[:-1]
    vals = [t for (el, c), t in ctx["by_class"].items() if el == donor and c and c[:-1] == env]
    if not vals or max(vals) - min(vals) > ELEMENT_SPREAD_MAX:
        return None
    return float(sum(vals) / len(vals))


def target_for_donor(metal: str, donor: str, atom=None, source: str = "md") -> Optional[float]:
    """Target M-D length of the complex being built, or None (caller keeps its value).

    With ``atom`` (a donor of a ligand fragment), the target of its class; without, or when the
    class is not among the complex's donors, the per-element target if the complex has one.
    """
    ctx = _CTX
    if ctx is None or not bv_targets_enabled() or metal != ctx["metal"]:
        return None
    if atom is None:
        scope = os.environ.get(SCOPE_ENV, "all")
        if scope == "seat" or (scope == "pins" and _PIN_DEPTH[0] <= 0):
            return None
        if scope == "lengths" and _PIN_DEPTH[0] <= 0 and (source != "table" or _CHECK_DEPTH[0] > 0):
            return None
    if atom is not None:
        try:
            if atom.GetSymbol() == donor:
                cls = donor_class(atom)
                t = ctx["by_class"].get((donor, cls))
                if t is not None:
                    return float(t)
                t = _same_environment_target(ctx, donor, cls)
                if t is not None:
                    return t
        except Exception:
            pass
    t = ctx["by_element"].get(donor)
    return None if t is None else float(t)
