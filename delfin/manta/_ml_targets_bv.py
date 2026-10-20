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
SCOPE_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_SCOPE"
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
    """``(kind, metal, ox, cn, element, class) -> (n, value)``; cn/class '*' where pooled."""
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
                if len(p) < 8 or p[0] not in ("D", "R"):
                    continue
                try:
                    tbl[(p[0], p[1], p[2], p[3], p[4], p[5])] = (int(p[6]), float(p[7]))
                except ValueError:
                    continue
    except OSError:
        tbl = {}
    _TABLE_CACHE[path] = tbl
    return tbl


def donor_target(table, metal: str, ox: int, cn: int, element: str, cls: Optional[str]):
    """(target A, level) of one donor, or (None, None)."""
    nmin = _min_n()
    row = table.get(("D", metal, str(ox), str(cn), element, cls))
    if row is not None and row[0] >= nmin:
        return row[1], "D"
    if ox <= 0 or cn <= 0:
        return None, None
    lns = B_BV * math.log(ox / float(cn))
    # effective R0 at this oxidation state, then over all oxidation states of the metal (the
    # bond-valence model carries the oxidation-state dependence through V); per class, then
    # per element
    for o, c, lev in ((str(ox), cls, "R"), (str(ox), "*", "Re"), ("*", cls, "Ro"), ("*", "*", "Roe")):
        row = table.get(("R", metal, o, "*", element, c))
        if row is not None and row[0] >= nmin:
            return row[1] - lns, lev
    if os.environ.get(LITERATURE_ENV, "0") == "1":
        r0 = r0_literature(metal, element)
        if r0 is not None:
            return r0 - lns, "L"
    return None, None


def complex_targets(table, metal: str, ox: int, donors) -> Dict[int, Tuple[float, str]]:
    """``donor idx -> (target, level)`` for the donors that get one."""
    cn = len(donors)
    out = {}
    for i, el, cls in donors:
        if el == "H" or cls is None:
            continue
        t, lev = donor_target(table, metal, ox, cn, el, cls)
        if t is not None:
            out[i] = (float(t), lev)
    return out


# ---------------------------------------------------------------------------------------------
# context of the complex being built
# ---------------------------------------------------------------------------------------------
def _mol_from_smiles(smiles: str):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
    return m


def build_context(smiles: str) -> Optional[dict]:
    """Targets of one complex: ``{"metal", "ox", "cn", "by_class", "by_element"}`` or None."""
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
    tg = complex_targets(load_table(table_path()), metal, ox, donors)
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
            "by_element": by_element, "levels": sorted({lv for _t, lv in tg.values()})}


def set_context(smiles) -> None:
    """Record the targets of the complex about to be built (no-op when the switch is off)."""
    global _CTX
    _CTX = None
    if not bv_targets_enabled() or not smiles:
        return
    try:
        _CTX = build_context(str(smiles))
    except Exception:
        _CTX = None


def clear_context() -> None:
    global _CTX
    _CTX = None


def current_context() -> Optional[dict]:
    return _CTX


def target_for_donor(metal: str, donor: str, atom=None) -> Optional[float]:
    """Target M-D length of the complex being built, or None (caller keeps its value).

    With ``atom`` (a donor of a ligand fragment), the target of its class; without, or when the
    class is not among the complex's donors, the per-element target if the complex has one.
    """
    ctx = _CTX
    if ctx is None or not bv_targets_enabled() or metal != ctx["metal"]:
        return None
    if atom is None and os.environ.get(SCOPE_ENV, "all") == "seat":
        return None
    if atom is not None:
        try:
            if atom.GetSymbol() == donor:
                t = ctx["by_class"].get((donor, donor_class(atom)))
                if t is not None:
                    return float(t)
        except Exception:
            pass
    t = ctx["by_element"].get(donor)
    return None if t is None else float(t)
