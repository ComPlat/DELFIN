"""Derivation of the bond-valence metal-donor target tables read by ``_ml_targets_bv``.

Library for ``tools/derive_ml_targets.py``; nothing here runs during a build.  Two corpora, one
classification (``_ml_targets_bv.donor_class`` on a closed-shell, metal-free ligand mol):

* SMILES + crystal (CSD, bring your own licence; or any refcode;SMILES list with crystal xyz files):
  oxidation state and donor classes from the SMILES cut that MANTA uses at build time
  (``_ml_targets_bv.complex_donors``); the crystal is mapped onto the SMILES heavy-atom graph
  (VF2), and the donors of each Weisfeiler-Lehman class of the SMILES graph take the nearest
  crystal atoms of that class -- so a monodentate carboxylate gets its BOUND O, whichever of the
  two symmetry-equivalent O the graph map happened to pick.  A structure whose SMILES donor is not
  bonded in the crystal (> 1.25 x the covalent-radius sum) is refused, not averaged in.
* COD (CC0) coordinates only: donors within 1.25 x the covalent-radius sum of the metal; each
  ligand perceived from its coordinates (RDKit DetermineBondOrders; of the charges with a
  closed-shell Lewis structure the one whose bond orders fit the bond lengths); the complex
  charge from the CIF formula unit = k x complex + counter-ions/solvents (exact decomposition, one
  charge); oxidation state = complex charge - sum of ligand charges.  The table keeps a COD row
  only when an independent second reading agrees: the bond-valence sum with literature R0 within
  +-0.5 of that oxidation state (``consensus_ok``).

Table (format of ``_ml_targets_bv.load_table``):
  N  metal ox cn/nx element class   direct median at this CN and number of X-type donors nx
  D  metal ox cn    element class   direct median at this CN
  R  metal ox *     element class   effective R0 = median(d + b ln s), s from the valence split
  R  metal ox *     element *
  P  b / rho                        the bond-valence parameters the R rows were computed with
The valence split: s_j = V w_j / sum(w), w = rho for an X-type donor (anionic, hydride), 1 for an
L-type donor; rho = 1 is Brown's equal-valence rule.  ``fit_bv`` chooses (b, rho) on one half of
the structures and scores it on the other half, next to the literature (b = 0.37, rho = 1).
"""
from __future__ import annotations

import hashlib
import math
import random
import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from delfin.manta import _elements as EL
from delfin.manta import _ml_targets_bv as BV

DONOR_FAC = 1.25          # a donor is bonded below this multiple of the covalent-radius sum
BOND_FAC = 1.20           # non-metal heavy-atom bonds in a crystal
BOND_FAC_LOOSE = 1.30     # union graph of the monomorphism fallback
TRANS_DEG = 150.0
ROOT_RE = re.compile(r"^([A-Z]{6})")


# ------------------------------------------------------------------------------------------------
# identifiers and formulas
# ------------------------------------------------------------------------------------------------
def root_of(rid: str) -> str:
    """Refcode family (ABCDEF for ABCDEF, ABCDEF01, ...)."""
    m = ROOT_RE.match(rid)
    return m.group(1) if m else rid


def formula_key(counts) -> str:
    return " ".join(f"{e}{counts[e]}" for e in sorted(counts))


def smiles_formula(smi: str) -> Optional[str]:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None
    m.UpdatePropertyCache(strict=False)
    c = Counter()
    for a in m.GetAtoms():
        c[a.GetSymbol()] += 1
        c["H"] += a.GetTotalNumHs()
    if c["H"] == 0:
        del c["H"]
    return formula_key(c)


def load_exclusions(paths: Sequence[str]) -> Tuple[set, set]:
    """(refcode roots, formulas) of every system in the given id;SMILES / id|SMILES / id files."""
    roots, forms = set(), set()
    for p in paths:
        with open(p) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                sep = ";" if ";" in ln else ("|" if "|" in ln else None)
                rid, smi = (ln.split(sep, 1) if sep else (ln, None))
                roots.add(root_of(rid.strip()))
                if smi:
                    f = smiles_formula(smi.strip())
                    if f:
                        forms.add(f)
    return roots, forms


def read_first_xyz_frame(path) -> Optional[Tuple[List[str], np.ndarray]]:
    """First frame of an xyz file as (elements, coordinates)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        return None
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return None
    el, X = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        if len(p) < 4:
            return None
        el.append(EL.normalise(p[0]))
        X.append([float(p[1]), float(p[2]), float(p[3])])
    if len(el) != n:
        return None
    return el, np.array(X, float)


_RC = EL.covalent_radius


# ------------------------------------------------------------------------------------------------
# SMILES -> crystal map
# ------------------------------------------------------------------------------------------------
def _smiles_graph(mol):
    import networkx as nx
    G = nx.Graph()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1 and not EL.is_metal(a.GetSymbol()):
            G.add_node(a.GetIdx(), el=a.GetSymbol(), nh=a.GetTotalNumHs(includeNeighbors=True))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in G and j in G:
            G.add_edge(i, j)
    return G


def wl_classes(G) -> Dict[int, str]:
    import networkx as nx
    lab = {i: f"{d['el']}{d['nh']}" for i, d in G.nodes(data=True)}
    nx.set_node_attributes(G, lab, "lab")
    h = nx.weisfeiler_lehman_subgraph_hashes(G, node_attr="lab", iterations=4)
    return {i: (h[i][-1] if h[i] else lab[i]) for i in G}


def _crystal_heavy_graph(el, X, fac):
    import networkx as nx
    n = len(el)
    heavy = [i for i in range(n) if el[i] != "H" and not EL.is_metal(el[i])]
    r = np.array([_RC(e) for e in el])
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    nh = {i: 0 for i in heavy}
    metals = [i for i in range(n) if EL.is_metal(el[i])]
    for h in (i for i in range(n) if el[i] == "H"):
        if not heavy:
            break
        j = min(heavy, key=lambda k: D[h, k])
        if D[h, j] > 1.6 and metals and min(D[h, m] for m in metals) < D[h, j]:
            continue                                   # hydride
        nh[j] += 1
    G = nx.Graph()
    for i in heavy:
        G.add_node(i, el=el[i], nh=nh[i])
    for a in range(len(heavy)):
        for b in range(a + 1, len(heavy)):
            i, j = heavy[a], heavy[b]
            if D[i, j] < fac * (r[i] + r[j]):
                G.add_edge(i, j)
    return G


def map_smiles_to_crystal(mol, el, X):
    """({SMILES heavy non-metal idx: crystal idx}, how) or (None, reason)."""
    from networkx.algorithms import isomorphism as iso
    Gs = _smiles_graph(mol)
    Gx = _crystal_heavy_graph(el, X, BOND_FAC)
    if sorted(d["el"] for _i, d in Gs.nodes(data=True)) != sorted(d["el"] for _i, d in Gx.nodes(data=True)):
        return None, "composition"
    tries = [(Gx, lambda a, b: a["el"] == b["el"] and a["nh"] == b["nh"], False),
             (Gx, lambda a, b: a["el"] == b["el"], False),
             (_crystal_heavy_graph(el, X, BOND_FAC_LOOSE), lambda a, b: a["el"] == b["el"], True)]
    for Gc, nm, mono in tries:
        gm = iso.GraphMatcher(Gc, Gs, node_match=nm)
        ok = gm.subgraph_is_monomorphic() if mono else gm.is_isomorphic()
        if ok:
            return {s: x for x, s in gm.mapping.items()}, ("mono" if mono else "iso")
    return None, "graph"


def crystal_donor_distances(mol, metal_idx, donors, el, X):
    """((SMILES donor -> crystal idx, crystal metal idx, distances to the metal), how) or (None, why).

    Donors of one WL class take the nearest crystal atoms of that class (automorphism-proof)."""
    mp, how = map_smiles_to_crystal(mol, el, X)
    if mp is None:
        return None, how
    msym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    xm = [i for i, e in enumerate(el) if e == msym]
    if len(xm) != 1 or sum(1 for e in el if EL.is_metal(e)) != 1:
        return None, "metal"
    xm = xm[0]
    cls = wl_classes(_smiles_graph(mol))
    Dm = np.linalg.norm(X - X[xm], axis=1)
    by_cls = defaultdict(list)
    for i in donors:
        by_cls[cls[i]].append(i)
    out = {}
    for c, ds in by_cls.items():
        cand = sorted((mp[s] for s, cc in cls.items() if cc == c and s in mp), key=lambda x: Dm[x])
        if len(cand) < len(ds):
            return None, "class"
        for s, x in zip(sorted(ds), cand[:len(ds)]):
            out[s] = x
    for s, x in out.items():
        if Dm[x] > DONOR_FAC * (_RC(msym) + _RC(el[x])):
            return None, "donor_not_bonded"
    return (out, xm, Dm), how


def trans_partners(X, xm, donor_x):
    out = {}
    for a in donor_x:
        va = X[a] - X[xm]
        best, bang = None, TRANS_DEG
        for b in donor_x:
            if b == a:
                continue
            vb = X[b] - X[xm]
            c = float(va @ vb / (np.linalg.norm(va) * np.linalg.norm(vb)))
            ang = float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))
            if ang > bang:
                best, bang = b, ang
        out[a] = best
    return out


def record_from_smiles_crystal(rid: str, smi: str, el, X) -> dict:
    """One record (per-donor crystal M-D, keyed like the table) or {'st': reason}."""
    from delfin.manta._oxidation_state import _OX_GROUP
    mol = BV._mol_from_smiles(smi)
    if mol is None:
        return {"id": rid, "st": "smiles"}
    mis = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP]
    if len(mis) != 1 or any(EL.is_metal(a.GetSymbol()) and a.GetIdx() != mis[0] for a in mol.GetAtoms()):
        return {"id": rid, "st": "not_one_tm"}
    cd = BV.complex_donors(mol, mis[0])
    if cd is None:
        return {"id": rid, "st": "cut_none"}
    metal, ox, donors = cd
    heavy = [(i, e, c) for i, e, c in donors if e != "H"]
    got, how = crystal_donor_distances(mol, mis[0], [i for i, _e, _c in heavy], el, X)
    if got is None:
        return {"id": rid, "st": "map_" + how}
    s2x, xm, Dm = got
    tp = trans_partners(X, xm, list(s2x.values()))
    x2c = {s2x[i]: (e, c) for i, e, c in heavy}
    rec = []
    for i, e, c in heavy:
        x = s2x[i]
        t = tp.get(x)
        rec.append([e, c, round(float(Dm[x]), 4)] + (list(x2c[t]) if t is not None else [None, None]))
    pos = {s2x[i]: k for k, (i, _e, _c) in enumerate(heavy)}
    trans = [pos.get(tp.get(s2x[i]), -1) for i, _e, _c in heavy]
    return {"id": rid, "st": "ok", "src": "smiles", "root": root_of(rid), "metal": metal, "ox": ox,
            "cn": len(donors), "nx": BV.n_x_donors(donors),
            "formula": formula_key(Counter(el)), "donors": rec, "trans": trans}


# ------------------------------------------------------------------------------------------------
# COD: ligands, charges and oxidation state from coordinates
# ------------------------------------------------------------------------------------------------
_Q_ORDER = (0, -1, -2, -3, 1, -4, 2)
MISFIT_TIE = 0.005
MAX_ITER = 200000


def _bond_misfit(m, X, idxs) -> float:
    from delfin.manta import polyhedra as PLY
    s, n = 0.0, 0
    for b in m.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        e1, e2 = a1.GetSymbol(), a2.GetSymbol()
        d = float(np.linalg.norm(X[idxs[a1.GetIdx()]] - X[idxs[a2.GetIdx()]]))
        if b.GetIsAromatic():
            exp = 0.5 * (PLY._pyykko_radius(e1, 1) + PLY._pyykko_radius(e2, 1)
                         + PLY._pyykko_radius(e1, 2) + PLY._pyykko_radius(e2, 2))
        else:
            o = b.GetBondTypeAsDouble()
            exp = PLY._pyykko_radius(e1, o) + PLY._pyykko_radius(e2, o)
        s += (d - exp) ** 2
        n += 1
    return (s / n) ** 0.5 if n else 0.0


def perceive_fragment(el, X, idxs, bonds):
    """(charge, sanitized mol) of one ligand from its coordinates, or (None, None)."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    loc = {g: k for k, g in enumerate(idxs)}
    rw = Chem.RWMol()
    for g in idxs:
        a = Chem.Atom(el[g])
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    for i, j in bonds:
        if i in loc and j in loc:
            rw.AddBond(loc[i], loc[j], Chem.BondType.SINGLE)
    conf = Chem.Conformer(len(idxs))
    for g, k in loc.items():
        conf.SetAtomPosition(k, X[g].tolist())
    base = rw.GetMol()
    base.AddConformer(conf, assignId=True)
    cands = []
    for q in _Q_ORDER:
        m = Chem.Mol(base)
        try:
            rdDetermineBonds.DetermineBondOrders(m, charge=q, allowChargedFragments=True,
                                                 maxIterations=MAX_ITER, embedChiral=False)
            if sum(a.GetFormalCharge() for a in m.GetAtoms()) != q:
                continue
            if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
                continue
            Chem.SanitizeMol(m)
            cands.append((_bond_misfit(m, X, idxs), abs(q), q, m))
        except Exception:
            continue
    if not cands:
        return None, None
    best = min(c[0] for c in cands)
    near = sorted((c for c in cands if c[0] <= best + MISFIT_TIE), key=lambda c: (c[1], c[0]))
    return near[0][2], near[0][3]


def perceive_complex(el, X, mi):
    """(sum of ligand charges, [(xyz idx, element, class, distance)], n_hydride, status)."""
    n = len(el)
    rm = _RC(el[mi])
    D = np.linalg.norm(X - X[mi], axis=1)
    nonmetal = [i for i in range(n) if i != mi and not EL.is_metal(el[i])]
    if len(nonmetal) != n - 1:
        return None, None, 0, "multi_metal"
    r = np.array([_RC(e) for e in el])
    DD = np.linalg.norm(X[:, None] - X[None], axis=2)
    bonds = [(a, b) for ai, a in enumerate(nonmetal) for b in nonmetal[ai + 1:]
             if DD[a, b] < BOND_FAC * (r[a] + r[b])]
    nbr = {i: set() for i in nonmetal}
    for a, b in bonds:
        nbr[a].add(b)
        nbr[b].add(a)
    donors, hyd = [], []
    for i in nonmetal:
        if el[i] == "H":
            if not nbr[i] and D[i] < rm + 0.31 + 0.25:
                hyd.append(i)
            continue
        if D[i] < DONOR_FAC * (rm + r[i]):
            donors.append(i)
    if not donors:
        return None, None, 0, "no_donor"
    if any(el[i] in ("B", "Si") for i in donors):
        return None, None, 0, "agostic_b_si"
    dset = set(donors)
    if any(nbr[i] & dset for i in donors):
        return None, None, 0, "hapto"
    seen, frags = set(hyd), []
    for i in nonmetal:
        if i in seen:
            continue
        comp, st = [i], [i]
        seen.add(i)
        while st:
            a = st.pop()
            for b in nbr[a]:
                if b not in seen:
                    seen.add(b)
                    comp.append(b)
                    st.append(b)
        frags.append(sorted(comp))
    fid = {g: k for k, fr in enumerate(frags) for g in fr}
    hv = [i for i in nonmetal if el[i] != "H"]
    for ai, a in enumerate(hv):
        for b in hv[ai + 1:]:
            if fid.get(a) != fid.get(b) and DD[a, b] < 2.2:
                return None, None, 0, "fragment_split"
    qsum = -len(hyd)
    cls = {}
    for fr in frags:
        if not (set(fr) & dset):
            return None, None, 0, "spectator_fragment"
        q, m = perceive_fragment(el, X, fr, bonds)
        if q is None:
            return None, None, 0, "perception_failed"
        qsum += q
        loc = {g: k for k, g in enumerate(fr)}
        for g in fr:
            a = m.GetAtomWithIdx(loc[g])
            if (g not in dset and a.GetSymbol() == "C" and a.GetFormalCharge() < 0
                    and not a.GetIsAromatic() and not any(b.GetIsConjugated() for b in a.GetBonds())):
                return None, None, 0, "carbanion_nondonor"
        for g in fr:
            if g in dset:
                c = BV.donor_class(m.GetAtomWithIdx(loc[g]))
                if c is None:
                    return None, None, 0, "class_none"
                if c.startswith("t0") and el[g] in ("O", "N") and D[g] > 1.80:
                    return None, None, 0, "terminal_without_h"
                cls[g] = c
    return qsum, [(g, el[g], cls[g], float(D[g])) for g in donors], len(hyd), "ok"


# Counter-ions and solvents a CIF formula unit may carry next to the complex: (formula, charge).
_SPECIES = [
    ("Cl", -1), ("Br", -1), ("I", -1), ("F", -1), ("Cl O4", -1), ("B F4", -1), ("F6 P", -1),
    ("N O3", -1), ("F6 Sb", -1), ("As F6", -1), ("C F3 O3 S", -1), ("O4 S", -2), ("C24 H20 B", -1),
    ("C32 H12 B F24", -1), ("C N S", -1), ("N3", -1), ("C N", -1), ("Cl4 Zn", -2), ("Cl4 Co", -2),
    ("Na", 1), ("K", 1), ("Li", 1), ("Cs", 1), ("Rb", 1), ("H4 N", 1), ("C8 H20 N", 1),
    ("C16 H36 N", 1), ("C4 H12 N", 1), ("C24 H20 P", 1), ("C36 H30 N P2", 1), ("As C24 H20", 1),
    ("H2 O", 0), ("C H4 O", 0), ("C2 H6 O", 0), ("C2 H3 N", 0), ("C H2 Cl2", 0), ("C H Cl3", 0),
    ("C3 H6 O", 0), ("C3 H7 N O", 0), ("C2 H6 O S", 0), ("C4 H10 O", 0), ("C4 H8 O", 0),
    ("C7 H8", 0), ("C6 H6", 0), ("C6 H14", 0), ("C5 H12", 0), ("C4 H8 O2", 0), ("C Cl4", 0),
    ("C5 H5 N", 0), ("C3 H8 O", 0), ("C6 H12", 0), ("C2 H4 Cl2", 0), ("C6 H5 Cl", 0),
]


def _parse_species(s):
    c = Counter()
    for e, nn in re.findall(r"([A-Z][a-z]?)(\d*)", s.replace(" ", "")):
        c[e] += int(nn) if nn else 1
    return c


_SPEC = [(_parse_species(f), q) for f, q in _SPECIES]
_SPEC_ELEMENTS = set().union(*[set(sp) for sp, _q in _SPEC])


def _charges_of(res_t, start, memo, budget):
    key = (res_t, start)
    if key in memo:
        return memo[key]
    budget[0] -= 1
    if budget[0] < 0:
        raise RuntimeError("decomposition budget")
    if not res_t:
        return {0}
    res = dict(res_t)
    out = set()
    for k in range(start, len(_SPEC)):
        sp, q = _SPEC[k]
        if all(res.get(e, 0) >= n for e, n in sp.items()):
            nxt = dict(res)
            for e, n in sp.items():
                nxt[e] -= n
            sub = _charges_of(tuple(sorted((e, v) for e, v in nxt.items() if v)), k, memo, budget)
            out |= {q + s for s in sub}
            if len(out) > 1:
                break
    memo[key] = out
    return out


def _decompose_charges(res) -> set:
    if set(e for e, v in res.items() if v) - _SPEC_ELEMENTS:
        return set()
    try:
        return _charges_of(tuple(sorted((e, v) for e, v in res.items() if v)), 0, {}, [20000])
    except RuntimeError:
        return set()


def cif_formula_sum(path) -> Optional[Counter]:
    try:
        with open(path, errors="replace") as fh:
            txt = fh.read().splitlines()
    except OSError:
        return None
    for k, ln in enumerate(txt):
        m = re.match(r"^_chemical_formula_sum\s*(.*)$", ln.strip())
        if m:
            val = m.group(1).strip()
            j = k + 1
            while (not val or val == ";") and j < len(txt):
                val = txt[j].strip()
                j += 1
            val = val.strip("'\"; ")
            c = Counter()
            for e, nn in re.findall(r"([A-Z][a-z]?)\s*([0-9.]*)", val):
                try:
                    c[e] += float(nn) if nn else 1.0
                except ValueError:
                    return None
            return c
    return None


def complex_charge_from_formula(fs: Optional[Counter], counts: Counter) -> Optional[int]:
    """Charge of the complex: formula unit = k complexes + species of _SPECIES (one charge), or None."""
    if not fs:
        return None
    for k in (1, 2, 3, 4, 0.5):
        res2 = Counter()
        ok = True
        for e in set(fs) | set(counts):
            v = 2.0 * (fs.get(e, 0.0) / k - counts.get(e, 0))
            if v < -1e-3 or abs(v - round(v)) > 1e-3:
                ok = False
                break
            if round(v):
                res2[e] = int(round(v))
        if not ok:
            continue
        if not res2:
            return 0
        qs = _decompose_charges(res2)
        if len(qs) == 1:
            q2 = qs.pop()
            if q2 % 2 == 0:
                return -(q2 // 2)
        if len(qs) > 1:
            return None                        # two readings of the residual: refuse
    return None


def bvs_literature(metal: str, pairs) -> Optional[float]:
    """Bond-valence sum with the O'Keeffe-Brese R0, b = 0.37, over (element, distance) pairs."""
    s = 0.0
    for e, d in pairs:
        if e == "H":
            continue
        r0 = BV.r0_literature(metal, e)
        if r0 is None:
            return None
        s += math.exp((r0 - d) / 0.37)
    return round(s, 3)


def record_from_cod(cid: str, el, X, formula_sum: Optional[Counter], excl_forms=()) -> dict:
    from delfin.manta._oxidation_state import _OX_GROUP
    metals = [i for i, e in enumerate(el) if EL.is_metal(e)]
    if len(metals) != 1 or el[metals[0]] not in _OX_GROUP:
        return {"id": cid, "st": "not_one_tm"}
    counts = Counter(el)
    fk = formula_key(counts)
    if fk in excl_forms:
        return {"id": cid, "st": "excluded"}
    qc = complex_charge_from_formula(formula_sum, counts)
    if qc is None:
        return {"id": cid, "st": "charge_unproven"}
    mi = metals[0]
    qsum, donors, nhyd, st = perceive_complex(el, X, mi)
    if donors is None:
        return {"id": cid, "st": st}
    ox = qc - qsum
    msym = el[mi]
    if ox < 0 or ox > 8 or not (0 <= _OX_GROUP[msym] - ox <= 10):
        return {"id": cid, "st": "ox_unphysical"}
    tp = trans_partners(X, mi, [g for g, *_ in donors])
    cl = {g: (e, c) for g, e, c, _d in donors}
    rec = []
    for g, e, c, d in donors:
        t = tp.get(g)
        rec.append([e, c, round(d, 4)] + (list(cl[t]) if t is not None else [None, None]))
    nx = nhyd + sum(1 for _g, _e, c, _d in donors if c and c.endswith("x"))
    return {"id": cid, "st": "ok", "src": "cod", "root": cid, "metal": msym, "ox": ox, "q": qc,
            "cn": len(donors) + nhyd, "nx": nx, "formula": fk, "donors": rec,
            "bvs": bvs_literature(msym, [(e, d) for _g, e, _c, d in donors])}


def consensus_ok(r: dict, window: float) -> bool:
    """COD row kept only if the BVS (literature R0) agrees with the counted oxidation state."""
    return r.get("bvs") is not None and int(r["ox"]) > 0 and abs(r["bvs"] - int(r["ox"])) <= window


# ------------------------------------------------------------------------------------------------
# the table
# ------------------------------------------------------------------------------------------------
def _valences(ox: int, rec_donors, nx: int, cn: int, rho: float):
    """s_j per donor: V w_j / sum w, w = rho for X-type (class ...x; hydrides counted in nx)."""
    nl = cn - nx
    tot = rho * nx + nl
    if tot <= 0 or ox <= 0:
        return None
    return [ox * (rho if (c and c.endswith("x")) else 1.0) / tot for _e, c, *_rest in rec_donors]


def _stable_half(rid: str) -> int:
    return int(hashlib.sha1(rid.encode()).hexdigest()[:8], 16) & 1


def _r0_table(recs, b, rho, min_n):
    acc = defaultdict(list)
    for r in recs:
        ox, cn, nx = int(r["ox"]), int(r["cn"]), int(r.get("nx", 0))
        s = _valences(ox, r["donors"], nx, cn, rho)
        if s is None:
            continue
        for (e, c, d, *_t), sj in zip(r["donors"], s):
            if c is None:
                continue
            r0 = d + b * math.log(sj)
            acc[(r["metal"], ox, e, c)].append(r0)
            acc[(r["metal"], ox, e, "*")].append(r0)
    return {k: float(np.median(v)) for k, v in acc.items() if len(v) >= min_n}


def fit_bv(recs, b_grid=None, rho_grid=None, min_n=30):
    """Choose (b, rho) for the R levels: R0 on one half of the structures, |error| on the other."""
    b_grid = b_grid or [round(0.05 + 0.025 * i, 3) for i in range(17)]          # 0.05 .. 0.45
    rho_grid = rho_grid or [round(0.5 + 0.25 * i, 2) for i in range(15)]        # 0.5 .. 4.0
    A = [r for r in recs if _stable_half(r["id"]) == 0]
    B = [r for r in recs if _stable_half(r["id"]) == 1]

    def score(b, rho):
        t = _r0_table(A, b, rho, min_n)
        err = []
        for r in B:
            ox, cn, nx = int(r["ox"]), int(r["cn"]), int(r.get("nx", 0))
            s = _valences(ox, r["donors"], nx, cn, rho)
            if s is None:
                continue
            for (e, c, d, *_t), sj in zip(r["donors"], s):
                r0 = t.get((r["metal"], ox, e, c))
                if r0 is not None:
                    err.append(abs(r0 - b * math.log(sj) - d))
        return (float(np.median(err)) if err else 9.9), len(err)

    grid = {(b, rho): score(b, rho) for b in b_grid for rho in rho_grid}
    best = min(grid, key=lambda k: grid[k][0])
    lit = score(0.37, 1.0)
    return {"b": best[0], "rho": best[1], "score": grid[best][0], "n_scored": grid[best][1],
            "literature_score": lit[0], "literature_n": lit[1], "grid": grid}


def _boot_ci(groups, B=200, seed=20261020):
    if len(groups) < 10 or B <= 0:
        return float("nan"), float("nan")
    rnd = random.Random(seed)
    k = len(groups)
    ms = []
    for _ in range(B):
        v = []
        for _j in range(k):
            v.extend(groups[rnd.randrange(k)])
        ms.append(float(np.median(v)))
    ms.sort()
    return ms[int(0.025 * B)], ms[int(0.975 * B) - 1]


def build_table(recs, b: float, rho: float, min_n: int = 10, boot: int = 200):
    """Rows (kind, metal, ox, cn, element, class, n_struct, value, ci_lo, ci_hi, p25, p75, n_bonds)."""
    tbl = {k: defaultdict(lambda: defaultdict(list)) for k in ("N", "D", "R")}
    for r in recs:
        m, ox, cn, nx = r["metal"], int(r["ox"]), int(r["cn"]), int(r.get("nx", 0))
        s = _valences(ox, r["donors"], nx, cn, rho)
        for j, (e, c, d, *_t) in enumerate(r["donors"]):
            if c is None:
                continue
            tbl["N"][(m, str(ox), f"{cn}/{nx}", e, c)][r["id"]].append(d)
            tbl["D"][(m, str(ox), str(cn), e, c)][r["id"]].append(d)
            if s is not None:
                r0 = d + b * math.log(s[j])
                tbl["R"][(m, str(ox), "*", e, c)][r["id"]].append(r0)
                tbl["R"][(m, str(ox), "*", e, "*")][r["id"]].append(r0)
    rows = []
    for kind in ("N", "D", "R"):
        for key in sorted(tbl[kind]):
            g = list(tbl[kind][key].values())
            if len(g) < min_n:
                continue
            v = [x for gg in g for x in gg]
            lo, hi = _boot_ci(g, boot)
            rows.append((kind,) + key + (len(g), float(np.median(v)), lo, hi,
                                         float(np.percentile(v, 25)), float(np.percentile(v, 75)), len(v)))
    return rows


def write_target_table(path, rows, b, rho, header_lines: Sequence[str]):
    with open(path, "w") as fh:
        for h in header_lines:
            fh.write(f"# {h}\n")
        fh.write("# kind\tmetal\tox\tcn\telement\tclass\tn_struct\tvalue\tci_lo\tci_hi\tp25\tp75\tn_bonds\n")
        fh.write(f"P\tb\t-\t-\t-\t-\t0\t{b:.4f}\n")
        fh.write(f"P\trho\t-\t-\t-\t-\t0\t{rho:.4f}\n")
        for r in rows:
            fh.write("\t".join([r[0], r[1], r[2], r[3], r[4], r[5], str(r[6])]
                               + [f"{x:.4f}" for x in r[7:12]] + [str(r[12])]) + "\n")


# ------------------------------------------------------------------------------------------------
# spin labels (phase 2): two-component split of the per-structure mean M-D deviation
# ------------------------------------------------------------------------------------------------
SPIN_POSTERIOR = 0.90      # a structure is labelled only when its spin component is this certain
SPIN_MIN_SEPARATION = 0.08  # A between the two component means of a bin, else it is unimodal


def _two_gaussians(x, iters=200):
    """1-D two-component Gaussian mixture by EM: (means, sds, weights, posterior of component 1)."""
    x = np.asarray(x, float)
    mu = np.array([np.percentile(x, 15), np.percentile(x, 85)])
    sd = np.array([x.std() / 2 + 1e-3] * 2)
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        pdf = np.stack([w[k] / (sd[k] * math.sqrt(2 * math.pi))
                        * np.exp(-0.5 * ((x - mu[k]) / sd[k]) ** 2) for k in (0, 1)])
        post = pdf / np.maximum(pdf.sum(axis=0), 1e-300)
        nk = post.sum(axis=1)
        w = nk / nk.sum()
        mu = (post * x).sum(axis=1) / np.maximum(nk, 1e-12)
        sd = np.sqrt((post * (x - mu[:, None]) ** 2).sum(axis=1) / np.maximum(nk, 1e-12)) + 1e-3
    return mu, sd, w, post[1]


def label_spins(recs, min_n: int = 30):
    """Set rec['spin'] = 'ls' / 'hs' for spin-candidate structures (3d, d4-d7, CN >= 5).

    Per bin (metal, ox, CN, element, class) the per-structure mean distance is split into two
    Gaussians.  A bin is BIMODAL when the means lie >= SPIN_MIN_SEPARATION apart and the smaller
    component holds >= 10 %; there the lower component is low spin, the upper high spin, and a
    structure is labelled when its component is >= SPIN_POSTERIOR certain.  A structure whose
    donors sit only in unimodal bins takes the ligand-field rule (``BV.spin_by_rule``: Co(III)
    low spin; Mn(II)/Co(II) low spin only with a C/P/H majority; otherwise low spin with a
    strong-field majority).  Conflicting labels of one structure -> no label.
    Returns a per-(metal, ox) summary (aggregates)."""
    bins = defaultdict(lambda: defaultdict(list))
    cand = []
    for r in recs:
        if not BV.spin_candidate(r["metal"], int(r["ox"]), int(r["cn"])):
            continue
        cand.append(r)
        for e, c, d, *_t in r["donors"]:
            if c is not None:
                bins[(r["metal"], int(r["ox"]), int(r["cn"]), e, c)][r["id"]].append(d)
    votes = defaultdict(set)
    bimodal = set()
    for key, per in bins.items():
        if len(per) < min_n:
            continue
        ids = list(per)
        x = [float(np.mean(per[i])) for i in ids]
        mu, _sd, w, p_hi = _two_gaussians(x)
        if abs(mu[1] - mu[0]) < SPIN_MIN_SEPARATION or min(w) < 0.10:
            continue
        bimodal.add(key)
        hs = 1 if mu[1] > mu[0] else 0
        for i, ph in zip(ids, p_hi):
            p_hs = ph if hs == 1 else 1.0 - ph
            if p_hs >= SPIN_POSTERIOR:
                votes[i].add("hs")
            elif p_hs <= 1.0 - SPIN_POSTERIOR:
                votes[i].add("ls")
    summary = defaultdict(lambda: {"n": 0, "ls": 0, "hs": 0, "rule": 0, "bimodal_bins": 0})
    for key in bimodal:
        summary[(key[0], key[1])]["bimodal_bins"] += 1
    for r in cand:
        s = summary[(r["metal"], int(r["ox"]))]
        s["n"] += 1
        in_bimodal = any((r["metal"], int(r["ox"]), int(r["cn"]), e, c) in bimodal
                         for e, c, *_t in r["donors"])
        v = votes.get(r["id"], set())
        if in_bimodal:
            if len(v) == 1:
                r["spin"] = v.pop()
        else:
            nh = int(r["cn"]) - len(r["donors"])
            r["spin"] = BV.spin_by_rule(r["metal"], int(r["ox"]), [(e, c) for e, c, *_t in r["donors"]], nh)
            s["rule"] += 1
        if r.get("spin"):
            s[r["spin"]] += 1
    return dict(summary)


def spin_rows(recs, min_n: int = 10, boot: int = 200):
    """N/D rows per spin (class suffixed ':ls' / ':hs') and F rows (low-spin share per metal, ox,
    CN and number of strong-field donors) from spin-labelled records."""
    tbl = {k: defaultdict(lambda: defaultdict(list)) for k in ("N", "D")}
    frac = defaultdict(lambda: [0, 0])
    for r in recs:
        sp = r.get("spin")
        if sp is None:
            continue
        m, ox, cn, nx = r["metal"], int(r["ox"]), int(r["cn"]), int(r.get("nx", 0))
        ns = BV.n_strong_field([(e, c) for e, c, *_t in r["donors"]], cn - len(r["donors"]))
        frac[(m, str(ox), f"{cn}/{ns}")][0 if sp == "ls" else 1] += 1
        for e, c, d, *_t in r["donors"]:
            if c is None:
                continue
            tbl["N"][(m, str(ox), f"{cn}/{nx}", e, f"{c}:{sp}")][r["id"]].append(d)
            tbl["D"][(m, str(ox), str(cn), e, f"{c}:{sp}")][r["id"]].append(d)
    rows = []
    for kind in ("N", "D"):
        for key in sorted(tbl[kind]):
            g = list(tbl[kind][key].values())
            if len(g) < min_n:
                continue
            v = [x for gg in g for x in gg]
            lo, hi = _boot_ci(g, boot)
            rows.append((kind,) + key + (len(g), float(np.median(v)), lo, hi,
                                         float(np.percentile(v, 25)), float(np.percentile(v, 75)), len(v)))
    for key in sorted(frac):
        ls, hs = frac[key]
        if ls + hs >= min_n:
            rows.append(("F",) + key + ("*", "*", ls + hs, ls / (ls + hs), float("nan"), float("nan"),
                                         float("nan"), float("nan"), ls + hs))
    return rows


# ------------------------------------------------------------------------------------------------
# Jahn-Teller form of a crystal (counting only; the table rows per form come with the JT switch)
# ------------------------------------------------------------------------------------------------
JT_SPLIT = 0.15            # A: one trans pair this much longer/shorter than the other two -> distorted
JT_TIGHT = 0.08            # A: the other two pairs agree within this


def jt_form(r: dict, ref_median=None) -> Optional[str]:
    """'elongated', 'compressed', 'undistorted' or 'other' for a CN6 record with three trans pairs;
    None when the record has no three clean trans pairs.  Pair lengths are taken relative to the
    class median (ref_median[(metal, ox, cn, element, class)]) when given, so a pair of long donors
    (e.g. two Br among N) is not read as a distortion."""
    tr = r.get("trans")
    if int(r["cn"]) != 6 or not tr or len(r["donors"]) != 6:
        return None
    pairs, seen = [], set()
    for k, t in enumerate(tr):
        if t < 0 or k in seen or tr[t] != k:
            continue
        seen |= {k, t}
        dev = []
        for j in (k, t):
            e, c, d, *_x = r["donors"][j]
            ref = ref_median.get((r["metal"], int(r["ox"]), 6, e, c)) if ref_median else None
            dev.append(d - ref if ref is not None else d)
        pairs.append(sum(dev) / 2.0)
    if len(pairs) != 3:
        return None
    p = sorted(pairs)
    if p[2] - p[1] >= JT_SPLIT and p[1] - p[0] <= JT_TIGHT:
        return "elongated"
    if p[1] - p[0] >= JT_SPLIT and p[2] - p[1] <= JT_TIGHT:
        return "compressed"
    if p[2] - p[0] <= JT_TIGHT:
        return "undistorted"
    return "other"


def jt_rows(recs, metal: str = "Cu", ox: int = 2, min_n: int = 10, boot: int = 200):
    """Axial and equatorial M-D medians per donor element for the Jahn-Teller forms of one ion
    (CN6, three trans pairs): rows ('J', metal, ox, '6', element, 'elong_ax' | 'elong_eq' |
    'comp_ax' | 'comp_eq').  The axial pair of an elongated octahedron is its longest trans pair,
    of a compressed one its shortest (pair lengths relative to the class medians)."""
    sel = [r for r in recs if r["metal"] == metal and int(r["ox"]) == ox and int(r["cn"]) == 6]
    acc = defaultdict(list)
    for r in sel:
        for e, c, d, *_t in r["donors"]:
            acc[(e, c)].append(d)
    ref = {(metal, ox, 6, e, c): float(np.median(v)) for (e, c), v in acc.items() if len(v) >= 10}
    groups = defaultdict(lambda: defaultdict(list))
    for r in sel:
        form = jt_form(r, ref)
        if form not in ("elongated", "compressed"):
            continue
        tr = r["trans"]
        pairs = []
        for k, t in enumerate(tr):
            if t > k and tr[t] == k:
                dev = [r["donors"][j][2] - ref.get((metal, ox, 6, r["donors"][j][0], r["donors"][j][1]),
                                                     r["donors"][j][2]) for j in (k, t)]
                pairs.append((sum(dev) / 2.0, k, t))
        pairs.sort()
        ax = pairs[-1] if form == "elongated" else pairs[0]
        tag = "elong" if form == "elongated" else "comp"
        for j, (e, c, d, *_x) in enumerate(r["donors"]):
            pos = "ax" if j in (ax[1], ax[2]) else "eq"
            groups[(e, f"{tag}_{pos}")][r["id"]].append(d)
    rows = []
    for (e, lab) in sorted(groups):
        g = list(groups[(e, lab)].values())
        if len(g) < min_n:
            continue
        v = [x for gg in g for x in gg]
        lo, hi = _boot_ci(g, boot)
        rows.append(("J", metal, str(ox), "6", e, lab, len(g), float(np.median(v)), lo, hi,
                     float(np.percentile(v, 25)), float(np.percentile(v, 75)), len(v)))
    return rows


def jt_rows_sp5(recs, metal: str = "Cu", ox: int = 2, min_n: int = 10, boot: int = 200):
    """Apical and basal M-D medians of square-pyramidal CN5 crystals of one ion (the donor without a
    trans partner while the other four form two trans pairs): rows ('J', metal, ox, '5', element,
    'sp5_ap' | 'sp5_bas')."""
    groups = defaultdict(lambda: defaultdict(list))
    for r in recs:
        if r["metal"] != metal or int(r["ox"]) != ox or int(r["cn"]) != 5 or len(r["donors"]) != 5:
            continue
        tr = r.get("trans") or []
        lone = [k for k, t in enumerate(tr) if t < 0]
        paired = [k for k, t in enumerate(tr) if t >= 0 and tr[t] == k]
        if len(lone) != 1 or len(paired) != 4:
            continue
        for j, (e, c, d, *_x) in enumerate(r["donors"]):
            groups[(e, "sp5_ap" if j == lone[0] else "sp5_bas")][r["id"]].append(d)
    rows = []
    for (e, lab) in sorted(groups):
        g = list(groups[(e, lab)].values())
        if len(g) < min_n:
            continue
        v = [x for gg in g for x in gg]
        lo, hi = _boot_ci(g, boot)
        rows.append(("J", metal, str(ox), "5", e, lab, len(g), float(np.median(v)), lo, hi,
                     float(np.percentile(v, 25)), float(np.percentile(v, 75)), len(v)))
    return rows
