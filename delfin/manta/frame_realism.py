"""Crystal-free realism ranking of the emitted frames (DELFIN_FRAME_REALISM_RANK, default OFF).

THE RULE (user decision 29.10.2026, pre-registration addendum 1): the best, most realistic frames come
FIRST; nothing is dropped.  Until now frame 0 was whatever the construction produced first -- under
DELFIN_FFFREE_UNION the FF-free frames by construction ("so frame 0 stays the FF-free pick"), and the
legacy least-clash ranker reads only steric overlap, so 86 % of the frames tie and the enumeration order
decides.  Geometry audit (29.10.2026): frame 0 defective main 42.4 % / lc18 47.5 % while the best frame
improved -- the frame a user sees first got worse.

THE SCORE.  Every frame is mapped onto the SMILES graph (metal bonds cut; RDKit subgraph match on
(element, H count) labels), and every local unit gets a value and a class key from its chemistry, the
same quantities as the crystal-calibrated geometry audit:
  pi   ring_oop (conjugated ring out-of-plane), chel_twist (chelating biaryl torsion)
  mpl  donor_oop (M-D out of the sp2 donor plane), ring_mdist (metal from a donor ring's plane)
  dang lp_misdir (M-D against the ideal lone pair), sp_donor (linear donors), mdx_1sub (M-D=X)
  hyb  sp2_sum, sp3_dev, sp_lin
  con  m_contact (metal...non-donor), nb_clash (pairs >= 4 bonds apart)
  bnd  bond_len, ml_len
For a unit with class median m and p99 (p1 on the low side) q, the exceedance is e = (v - m) / (q - m)
(one- or two-sided as per quantity), and

    score = sum over units of max(0, e - 1)^2  +  0.001 * sum of e^2

-- zero while every unit sits inside the crystal p99, growing quadratically beyond it; the small second
term orders frames that are all inside.  Frames are sorted by (unmappable last, score, emitted index):
deterministic, a permutation (verified), labels travel with their frames.

THE TABLE holds only class statistics (n, median, p1, p99 per class, from crystal structures), never a
coordinate: DELFIN_FRAME_REALISM_TABLE=<json>, else the shipped COD table next to this module if present.
Without a table the pass is a no-op.  The CSD-derived table is built locally per licence and is never in
the repository; the shipped default must come from COD (CC0).

Off -> nothing here runs (byte-identical).
"""
from __future__ import annotations

import json
import math
import os
import threading
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

RANK_ENV = "DELFIN_FRAME_REALISM_RANK"
TABLE_ENV = "DELFIN_FRAME_REALISM_TABLE"
SHIPPED_TABLE = os.path.join(os.path.dirname(__file__), "data", "frame_realism_cod.json")

CAT = {"ring_oop": "pi", "chel_twist": "pi", "donor_oop": "mpl", "ring_mdist": "mpl", "lp_misdir": "dang",
       "sp_donor": "dang", "mdx_1sub": "dang", "sp2_sum": "hyb", "sp3_dev": "hyb", "sp_lin": "hyb",
       "m_contact": "con", "nb_clash": "con", "bond_len": "bnd", "ml_len": "bnd"}
LOW = {"m_contact", "nb_clash"}            # small = bad
TWO = {"mdx_1sub", "bond_len", "ml_len"}  # both sides

# Pyykko 2009 single-bond covalent radii and Bondi / Alvarez 2013 van der Waals radii (A), as in the
# eye's find_inter_ligand_clash tables (H vdW 1.10 as in the audit), so that the class values mean the
# same as in the table builder.
_RCOV = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Ne": 0.58, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15, "La": 2.07,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}
_VDW_TAB = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75,
    "Ar": 1.88, "K": 2.75, "Ca": 2.31, "Sc": 2.11, "Ti": 1.95, "V": 1.92, "Cr": 1.89, "Mn": 1.97,
    "Fe": 1.94, "Co": 1.92, "Ni": 1.84, "Cu": 1.86, "Zn": 2.10, "Ga": 1.87, "Ge": 2.11, "As": 1.85,
    "Se": 1.90, "Br": 1.83, "Kr": 2.02, "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18,
    "Mo": 2.17, "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18, "In": 1.93,
    "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16, "Cs": 3.43, "Ba": 2.68, "La": 2.43,
    "Ce": 2.42, "Pr": 2.40, "Nd": 2.39, "Pm": 2.38, "Sm": 2.36, "Eu": 2.35, "Gd": 2.34, "Tb": 2.33,
    "Dy": 2.31, "Ho": 2.30, "Er": 2.29, "Tm": 2.27, "Yb": 2.26, "Lu": 2.24, "Hf": 2.23, "Ta": 2.22,
    "W": 2.18, "Re": 2.16, "Os": 2.16, "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23, "Tl": 1.96,
    "Pb": 2.02, "Bi": 2.07,
}
_VDW_AUDIT = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "Si": 2.10, "Cl": 1.75,
              "Br": 1.85, "I": 1.98, "F": 1.47, "B": 1.92, "Se": 1.90, "As": 1.85, "Te": 2.06, "Ge": 2.11,
              "Sn": 2.17, "Sb": 2.06}
_METAL_Z = (set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81)) | set(range(89, 104))
            | {3, 11, 19, 37, 55, 4, 12, 20, 38, 56, 13, 31, 49, 81, 50, 82, 51, 83, 52, 84})


def rank_enabled() -> bool:
    return os.environ.get(RANK_ENV, "0") == "1"


def _fr_rcov(e: str) -> float:
    r = _RCOV.get(e)
    if r is None:
        try:
            from rdkit import Chem
            pt = Chem.GetPeriodicTable()
            r = pt.GetRcovalent(pt.GetAtomicNumber(e))
        except Exception:
            r = 1.5
    return float(r)


def _fr_vdw(e: str) -> float:
    return float(_VDW_AUDIT.get(e, _VDW_TAB.get(e, 2.0)))


def _fr_is_metal(e: str) -> bool:
    try:
        from rdkit import Chem
        return Chem.GetPeriodicTable().GetAtomicNumber(e) in _METAL_Z
    except Exception:
        return False


# ------------------------------------------------------------------------------------------- table
_TABLES: Dict[str, Optional[dict]] = {}


def table_path() -> str:
    p = os.environ.get(TABLE_ENV, "")
    return p if p else SHIPPED_TABLE


def load_table(path: Optional[str] = None) -> Optional[dict]:
    path = path or table_path()
    if path in _TABLES:
        return _TABLES[path]
    t = None
    try:
        with open(path) as fh:
            raw = json.load(fh)
        t = {k: {kk: (float(v["med"]), float(v["lo1"]), float(v["hi1"])) for kk, v in raw.get(k, {}).items()}
             for k in ("fine", "coarse", "all")}
    except Exception:
        t = None
    _TABLES[path] = t
    return t


def _fr_excess(T: dict, q: str, c: str, cc: str, v: float) -> float:
    t = T["fine"].get(f"{q}|{c}") or T["coarse"].get(f"{q}|{cc}") or T["all"].get(q)
    if t is None:
        return 0.0
    med, lo1, hi1 = t
    up = q not in LOW and (q not in TWO or v >= med)
    if up:
        e = (v - med) / max(hi1 - med, 1e-6)
    else:
        e = (med - v) / max(med - lo1, 1e-6)
    return max(e, 0.0)


# ---------------------------------------------------------------------------------- SMILES chemistry
def _fr_smiles_graph(smi: str) -> Optional[dict]:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None
    m.UpdatePropertyCache(strict=False)
    el = [a.GetSymbol() for a in m.GetAtoms()]
    nh = [a.GetTotalNumHs() for a in m.GetAtoms()]
    met = [_fr_is_metal(e) for e in el]
    bonds = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append((i, j, met[i] or met[j], b.GetBondTypeAsDouble()))
    return {"el": el, "nh": nh, "met": met, "bonds": bonds}


def _fr_chem(smi: str) -> Optional[dict]:
    """SMILES chemistry: hybridisation, aromaticity, rings, donors, chelate ring sizes, WL classes."""
    import networkx as nx
    from rdkit import Chem
    sg = _fr_smiles_graph(smi)
    if sg is None:
        return None
    m = Chem.MolFromSmiles(smi, sanitize=False)
    rw = Chem.RWMol(m)
    met = sg["met"]
    for b in list(m.GetBonds()):
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if met[i] != met[j]:
            d = j if met[i] else i
            rw.RemoveBond(i, j)
            a = rw.GetAtomWithIdx(d)
            a.SetFormalCharge(a.GetFormalCharge() - 1)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(sg["nh"][d])
    for a in rw.GetAtoms():
        if not met[a.GetIdx()]:
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(sg["nh"][a.GetIdx()])
    mol = rw.GetMol()
    ok = True
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        ok = False
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
                             catchErrors=True)
        except Exception:
            pass
    n = mol.GetNumAtoms()
    hyb, arom = [], []
    HY = Chem.rdchem.HybridizationType
    for a in mol.GetAtoms():
        i = a.GetIdx()
        arom.append(bool(a.GetIsAromatic()))
        if met[i] or sg["el"][i] == "H":
            hyb.append("x")
            continue
        orders = [b.GetBondTypeAsDouble() for b in a.GetBonds()]
        if ok:
            h = a.GetHybridization()
            h = "sp" if h == HY.SP else "sp2" if h == HY.SP2 else "sp3" if h == HY.SP3 else "x"
        else:
            if any(o == 3 for o in orders) or sum(o == 2 for o in orders) >= 2:
                h = "sp"
            elif any(o in (1.5, 2) for o in orders) or a.GetIsAromatic():
                h = "sp2"
            else:
                h = "sp3"
        hyb.append(h)
    multi = [any(b.GetBondTypeAsDouble() > 1.0 or b.GetIsAromatic() for b in a.GetBonds())
             for a in mol.GetAtoms()]
    conjx = [hyb[i] == "sp2" and not multi[i] and sg["el"][i] in ("N", "O", "S", "P") for i in range(n)]
    ri = mol.GetRingInfo()
    rings = [tuple(r) for r in ri.AtomRings()]
    minring = [0] * n
    for r in rings:
        for i in r:
            minring[i] = len(r) if minring[i] == 0 else min(minring[i], len(r))
    conj_rings = [r for r in rings if len(r) <= 8 and all(hyb[i] == "sp2" or arom[i] for i in r)]
    order = {}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        order[(min(i, j), max(i, j))] = 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()
    G = nx.Graph()
    for i, e in enumerate(sg["el"]):
        if e != "H":
            G.add_node(i)
    for i, j, mb, _ in sg["bonds"]:
        if i in G and j in G:
            G.add_edge(i, j)
    donors: Dict[int, set] = {}
    for i, j, mb, _ in sg["bonds"]:
        if mb and sg["el"][i] != "H" and sg["el"][j] != "H":
            for mm, d in ((i, j), (j, i)):
                if met[mm] and not met[d]:
                    donors.setdefault(mm, set()).add(d)
    Go = G.subgraph([i for i in G if not met[i]])
    chel = {}
    for mm, ds in donors.items():
        for d in ds:
            best = 0
            for d2 in ds:
                if d2 == d:
                    continue
                try:
                    L = nx.shortest_path_length(Go, d, d2)
                except Exception:
                    continue
                cand = L + 2
                best = cand if best == 0 else min(best, cand)
            chel[(mm, d)] = best
    # organic heavy graph: ligands, WL classes, (element, nH) labels
    Gs = nx.Graph()
    for i, e in enumerate(sg["el"]):
        if e != "H" and not met[i]:
            Gs.add_node(i, lab=f"{e}{sg['nh'][i]}")
    for i, j, mb, _ in sg["bonds"]:
        if not mb and i in Gs and j in Gs:
            Gs.add_edge(i, j)
    lig = {}
    for k, c in enumerate(sorted(nx.connected_components(Gs), key=min)):
        for i in c:
            lig[i] = k
    hsh = nx.weisfeiler_lehman_subgraph_hashes(Gs, node_attr="lab", iterations=4)
    cls = {i: hsh[i][-1] if hsh[i] else Gs.nodes[i]["lab"] for i in Gs}
    sp = dict(nx.all_pairs_shortest_path_length(G, cutoff=4))
    return {"sg": sg, "hyb": hyb, "arom": arom, "conjx": conjx, "rings": rings, "minring": minring,
            "conj_rings": conj_rings, "order": order, "G": G, "donors": donors, "chel": chel,
            "Gs_nodes": sorted(Gs.nodes), "Gs_edges": sorted(Gs.edges), "lig": lig, "cls": cls, "sp": sp}


# ------------------------------------------------------------------------------------------ mapping
def _fr_parse_xyz(xyz: str) -> Tuple[List[str], np.ndarray]:
    el, X = [], []
    for ln in str(xyz).splitlines():
        p = ln.split()
        if len(p) == 4:
            try:
                X.append([float(p[1]), float(p[2]), float(p[3])])
                el.append(p[0])
            except ValueError:
                continue
    return el, np.asarray(X, float).reshape(-1, 3)


def _fr_hparents(el, X) -> Dict[int, int]:
    n = len(el)
    nm = [i for i in range(n) if el[i] != "H" and not _fr_is_metal(el[i])]
    mets = [i for i in range(n) if _fr_is_metal(el[i])]
    hp: Dict[int, int] = {}
    if not nm:
        return hp
    Xn = X[nm]
    for h in range(n):
        if el[h] != "H":
            continue
        d = np.linalg.norm(Xn - X[h], axis=1)
        k = int(np.argmin(d))
        if d[k] > 1.6 and mets and min(np.linalg.norm(X[mm] - X[h]) for mm in mets) < d[k]:
            continue
        hp[h] = nm[k]
    return hp


def _fr_label_mol(labels, edges):
    """Plain mol of a labelled graph: element per atom, isotope = H count + 1 (or 0 = unlabelled), single bonds."""
    from rdkit import Chem
    t = Chem.RWMol()
    for e, h in labels:
        a = Chem.Atom(e)
        if h is not None:
            a.SetIsotope(int(h) + 1)
        a.SetNoImplicit(True)
        t.AddAtom(a)
    for a, b in edges:
        t.AddBond(int(a), int(b), Chem.BondType.SINGLE)
    m = t.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def _fr_canon(m):
    from rdkit import Chem
    return Chem.MolToSmiles(m), list(Chem.CanonicalRankAtoms(m, breakTies=True))


def _fr_map_frame(ch: dict, el, X, prev) -> Optional[Dict[int, int]]:
    """SMILES non-metal heavy index -> xyz index for one frame (prev reused if all its bonds hold).

    Isomorphism only, by canonical form: the frame's distance graph (1.20, then 1.25 and 1.30 x the covalent sum)
    must have the canonical SMILES of the labelled SMILES graph; equal canonical rank (ties broken) is then the
    same atom.  Polynomial, no backtracking search (a subgraph search on BEXVAC's menthyl/cyclohexyl phosphine
    ran for hours).  A frame whose bond graph is not the SMILES graph is unmappable and ranks last."""
    sg = ch["sg"]
    n = len(el)
    nm = [i for i in range(n) if el[i] != "H" and not _fr_is_metal(el[i])]
    if sorted(sg["el"][s] for s in ch["Gs_nodes"]) != sorted(el[i] for i in nm):
        return None
    r = np.array([_fr_rcov(e) for e in el])
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    lim = r[:, None] + r[None]
    if prev is not None:
        if all(D[prev[a], prev[b]] < 1.30 * lim[prev[a], prev[b]] for a, b in ch["Gs_edges"]):
            return prev
    hp = _fr_hparents(el, X)
    nhx = {i: 0 for i in nm}
    for h, p in hp.items():
        nhx[p] += 1
    lab_ok = (sorted((sg["el"][s], sg["nh"][s]) for s in ch["Gs_nodes"])
              == sorted((el[i], nhx[i]) for i in nm))
    key = "canon_h" if lab_ok else "canon_el"
    qorder = ch["Gs_nodes"]
    if key not in ch:
        pos = {s: k for k, s in enumerate(qorder)}
        qm = _fr_label_mol([(sg["el"][s], sg["nh"][s] if lab_ok else None) for s in qorder],
                           [(pos[a], pos[b]) for a, b in ch["Gs_edges"]])
        ch[key] = _fr_canon(qm)
    qsmi, qrank = ch[key]
    by_rank = {rk: qorder[k] for k, rk in enumerate(qrank)}
    sub_all = D[np.ix_(nm, nm)]
    lim_nm = lim[np.ix_(nm, nm)]
    labels = [(el[i], nhx[i] if lab_ok else None) for i in nm]
    for fac in (1.20, 1.25, 1.30):
        ia, ib = np.nonzero(np.triu(sub_all < fac * lim_nm, 1))
        if len(ia) != len(ch["Gs_edges"]):
            continue
        try:
            tsmi, trank = _fr_canon(_fr_label_mol(labels, zip(ia.tolist(), ib.tolist())))
        except Exception:
            continue
        if tsmi != qsmi:
            continue
        mp = {by_rank[trank[k]]: nm[k] for k in range(len(nm))}
        if all(D[mp[a], mp[b]] < fac * lim[mp[a], mp[b]] for a, b in ch["Gs_edges"]):
            return mp
    return None


def _fr_metal_map(sg, mp, el, X) -> Optional[Dict[int, int]]:
    smet = [i for i, m in enumerate(sg["met"]) if m]
    xmet = [i for i, e in enumerate(el) if _fr_is_metal(e)]
    nb = {m: [] for m in smet}
    for i, j, mb, _ in sg["bonds"]:
        if mb:
            for m, d in ((i, j), (j, i)):
                if m in nb and d in mp:
                    nb[m].append(mp[d])
    used, out = set(), {}
    for m in sorted(smet, key=lambda m: -len(nb[m])):
        cand = [x for x in xmet if el[x] == sg["el"][m] and x not in used]
        if not cand:
            return None
        best = (min(cand, key=lambda x: np.mean([np.linalg.norm(X[x] - X[d]) for d in nb[m]]))
                if nb[m] else cand[0])
        out[m] = best
        used.add(best)
    return out


# ------------------------------------------------------------------------------------------ geometry
def _fr_ang(u, v) -> Optional[float]:
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-6 or nv < 1e-6:
        return None
    return float(np.degrees(np.arccos(np.clip(u @ v / nu / nv, -1, 1))))


def _fr_plane(P):
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c)
    return c, vt[2]


def _fr_dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - (b0 @ b1n) * b1n
    w = b2 - (b2 @ b1n) * b1n
    return float(np.degrees(np.arctan2(np.cross(b1n, v) @ w, v @ w)))


def frame_units(ch: dict, el, X, mp, mm, hp) -> List[Tuple[str, str, str, float]]:
    """(quantity, fine class, coarse class, value) of every unit of one mapped frame (= the audit's units)."""
    sg = ch["sg"]
    cls = ch["cls"]
    hyb, arom, conjx = ch["hyb"], ch["arom"], ch["conjx"]
    out = []
    pos = {s: X[x] for s, x in mp.items()}
    for s, x in mm.items():
        pos[s] = X[x]
    Hof: Dict[int, list] = {}
    inv = {x: s for s, x in mp.items()}
    for h, p in hp.items():
        if p in inv:
            Hof.setdefault(inv[p], []).append(X[h])
    nbr = {s: [t for t in ch["G"].neighbors(s) if not sg["met"][t]] for s in mp}
    subs = {s: [pos[t] for t in nbr[s]] + Hof.get(s, []) for s in mp}
    nh_ok = {s: len(Hof.get(s, [])) == sg["nh"][s] for s in mp}
    fdon = {}
    for m, ds in ch["donors"].items():
        if m not in mm:
            continue
        want: Dict[str, list] = {}
        for d in sorted(ds):
            want.setdefault(cls.get(d), []).append(d)
        got = []
        for c, lst in want.items():
            pool = [s for s in sorted(mp) if cls.get(s) == c]
            pool.sort(key=lambda s: np.linalg.norm(pos[s] - pos[m]))
            got += pool[:len(lst)]
        fdon[m] = got
    chelsz = {}
    for m, got in fdon.items():
        for s in got:
            cands = [d for d in ch["donors"][m] if cls.get(d) == cls.get(s)]
            chelsz[(m, s)] = min(ch["chel"].get((m, d), 0) for d in cands) if cands else 0
    donor_rings = []
    for r in ch["conj_rings"]:
        P = np.array([pos[i] for i in r])
        c, nv = _fr_plane(P)
        dev = np.abs((P - c) @ nv)
        ar = all(arom[i] for i in r)
        out.append(("ring_oop", f"{len(r)}{'a' if ar else 'c'}", f"{len(r)}", float(dev.max())))
        for m, got in fdon.items():
            ds = [s for s in got if s in r]
            if ds:
                donor_rings.append((m, r))
                out.append(("ring_mdist", f"{len(r)}{'a' if ar else 'c'}-{sg['el'][ds[0]]}", f"{len(r)}",
                            float(abs((pos[m] - c) @ nv))))
    seen = set()
    for m, r1 in donor_rings:
        for m2, r2 in donor_rings:
            if m2 != m or r1 == r2:
                continue
            for a in r1:
                for b in nbr.get(a, []):
                    if b in r2 and a not in r2 and b not in r1:
                        key = (min(a, b), max(a, b))
                        if key in seen:
                            continue
                        seen.add(key)
                        pa = [t for t in nbr[a] if t in r1][0]
                        pb = [t for t in nbr[b] if t in r2][0]
                        t = _fr_dihedral(pos[pa], pos[a], pos[b], pos[pb])
                        out.append(("chel_twist", "biaryl", "biaryl", min(abs(t), 180 - abs(t))))
    for m, got in fdon.items():
        me = sg["el"][m]
        for d in got:
            de = sg["el"][d]
            v = pos[m] - pos[d]
            S = subs[d]
            csz = chelsz.get((m, d), 0)
            cz = "mono" if csz == 0 else str(min(csz, 8))
            h = hyb[d]
            if len(S) == 1 and len(nbr[d]) == 1 and hyb[nbr[d][0]] == "sp":
                h = "sp"
            ar = "a" if arom[d] else ""
            out.append(("ml_len", f"{me}-{de}{h}", f"{me}-{de}", float(np.linalg.norm(v))))
            if not nh_ok[d]:
                continue
            if h == "sp2" and len(S) == 2:
                c, nv = _fr_plane(np.array([pos[d]] + S))
                a0 = _fr_ang(v, nv)
                if a0 is not None:
                    out.append(("donor_oop", f"{de}{ar}2-{cz}", f"{de}2", abs(90 - a0)))
                u = [(p - pos[d]) / np.linalg.norm(p - pos[d]) for p in S]
                a = _fr_ang(v, -(u[0] + u[1]))
                if a is not None:
                    out.append(("lp_misdir", f"{de}{ar}sp2-{cz}", f"{de}sp2", a))
            elif h == "sp2" and len(S) == 1:
                x = nbr[d][0] if nbr[d] else None
                if x is not None:
                    a = _fr_ang(v, pos[x] - pos[d])
                    if a is not None:
                        out.append(("mdx_1sub", f"{de}={sg['el'][x]}-{cz}", f"{de}1", a))
                    Sx = subs[x]
                    if len(Sx) == 3:
                        c, nv = _fr_plane(np.array(Sx + [pos[x]]))
                        aa = _fr_ang(v, nv)
                        if aa is not None:
                            out.append(("donor_oop", f"{de}1-{cz}", f"{de}1", abs(90 - aa)))
            elif h == "sp3" and len(S) == 3:
                u = [(p - pos[d]) / np.linalg.norm(p - pos[d]) for p in S]
                a = _fr_ang(v, -(u[0] + u[1] + u[2]))
                if a is not None:
                    out.append(("lp_misdir", f"{de}sp3-{cz}", f"{de}sp3", a))
            elif h == "sp" and len(S) == 1:
                a = _fr_ang(v, S[0] - pos[d])
                if a is not None:
                    out.append(("sp_donor", f"{de}sp-{cz}", f"{de}sp", 180 - a))
    for s in sorted(mp):
        if not nh_ok[s]:
            continue
        S = subs[s]
        e, h = sg["el"][s], hyb[s]
        if h == "sp2" and len(S) == 3:
            u = [p - pos[s] for p in S]
            tot = sum(_fr_ang(u[a], u[b]) or 0 for a, b in ((0, 1), (0, 2), (1, 2)))
            kind = "a" if arom[s] else ("j" if conjx[s] else "d")
            out.append(("sp2_sum", f"{e}{kind}", e, 360 - tot))
        elif h == "sp3" and len(S) >= 3:
            u = [p - pos[s] for p in S]
            dv = max(abs((_fr_ang(u[a], u[b]) or 109.47) - 109.47)
                     for a in range(len(u)) for b in range(a + 1, len(u)))
            rr = ch["minring"][s]
            rk = "0" if rr == 0 else str(min(rr, 7))
            out.append(("sp3_dev", f"{e}r{rk}n{len(S)}", f"{e}r{rk}", dv))
        elif h == "sp" and len(S) == 2:
            a = _fr_ang(S[0] - pos[s], S[1] - pos[s])
            if a is not None:
                out.append(("sp_lin", e, "sp", 180 - a))
    for (i, j), o in sorted(ch["order"].items()):
        if (i in pos and j in pos and not sg["met"][i] and not sg["met"][j]
                and sg["el"][i] != "H" and sg["el"][j] != "H"):
            a, b = sorted((sg["el"][i], sg["el"][j]))
            d = np.linalg.norm(pos[i] - pos[j])
            out.append(("bond_len", f"{a}{b}{o:g}", f"{a}{b}", float(d / (_fr_rcov(a) + _fr_rcov(b)))))
    sp = ch["sp"]
    for m in sorted(mm):
        don = set(fdon.get(m, [])) | set(ch["donors"].get(m, set()))
        for s in sorted(mp):
            if s in don:
                continue
            t = sp.get(m, {}).get(s, 9)
            e = sg["el"][s]
            ek = e if e in ("C", "N", "O") else "X"
            d = np.linalg.norm(pos[s] - pos[m])
            out.append(("m_contact", f"t{min(t, 4)}{ek}", f"t{min(t, 4)}",
                        float(d / (_fr_rcov(sg["el"][m]) + _fr_rcov(e)))))
    hs = [s for s in sorted(pos) if sg["el"][s] != "H" and not sg["met"][s]]
    if len(hs) > 1:
        P = np.array([pos[s] for s in hs])
        D = np.linalg.norm(P[:, None] - P[None], axis=2)
        vd = np.array([_fr_vdw(sg["el"][s]) for s in hs])
        R = D / (vd[:, None] + vd[None])
        lig = ch["lig"]
        ia, ib = np.nonzero(np.triu(R <= 0.85, 1))
        for a, b in zip(ia.tolist(), ib.tolist()):
            sa, sb = hs[a], hs[b]
            t = sp.get(sa, {}).get(sb, 99)
            if t < 4:
                continue
            kind = "intra" if lig.get(sa) == lig.get(sb) else "inter"
            out.append(("nb_clash", f"{kind}{'4' if t == 4 else '5'}", kind, float(R[a, b])))
    return out


def frame_score(T: dict, units) -> float:
    s = soft = 0.0
    for q, c, cc, v in units:
        e = _fr_excess(T, q, c, cc, v)
        s += max(e - 1.0, 0.0) ** 2
        soft += e * e
    return s + 0.001 * soft


# ------------------------------------------------------------------------------------------ context
_CTX = threading.local()
_CHEM_CACHE: Dict[str, Optional[dict]] = {}


def set_context(smiles: Optional[str]) -> None:
    _CTX.smiles = smiles


def clear_context() -> None:
    _CTX.smiles = None


def current_smiles() -> Optional[str]:
    return getattr(_CTX, "smiles", None)


def _fr_chem_cached(smiles: str) -> Optional[dict]:
    if smiles not in _CHEM_CACHE:
        try:
            _CHEM_CACHE[smiles] = _fr_chem(smiles)
        except Exception:
            _CHEM_CACHE[smiles] = None
        if len(_CHEM_CACHE) > 64:
            _CHEM_CACHE.pop(next(iter(_CHEM_CACHE)))
    return _CHEM_CACHE[smiles]


def score_frames(xyzs: Sequence[str], smiles: Optional[str] = None) -> Optional[List[Optional[float]]]:
    """Realism score of every frame (None = not mappable onto the SMILES graph), or None when the pass
    cannot run (no table, no SMILES, unreadable SMILES)."""
    smiles = smiles or current_smiles()
    T = load_table()
    if not smiles or T is None:
        return None
    ch = _fr_chem_cached(smiles)
    if ch is None:
        return None
    out: List[Optional[float]] = []
    prev = None
    for xyz in xyzs:
        try:
            el, X = _fr_parse_xyz(xyz)
            mp = _fr_map_frame(ch, el, X, prev)
            if mp is None:
                out.append(None)
                continue
            prev = mp
            hp = _fr_hparents(el, X)
            best = None
            # THE MINIMUM OVER THE MOLECULE'S AUTOMORPHISMS: a symmetric ligand can be mapped onto the frame
            # in several ways (which of two equivalent N is bound, which C-N of an amidinate is the double
            # bond), and the bond-order, contact-distance and chelate-size classes depend on that choice.
            # The frame's realism is that of its best reading; every candidate is a graph isomorphism.
            for sig in _fr_automorphisms(ch):
                mps = {s: mp[sig[s]] for s in mp}
                mm = _fr_metal_map(ch["sg"], mps, el, X)
                if mm is None:
                    continue
                v = frame_score(T, frame_units(ch, el, X, mps, mm, hp))
                best = v if best is None or v < best else best
            out.append(round(best, 6) if best is not None else None)
        except Exception:
            out.append(None)
    return out


_AUTO_CAP = 16


_AUTO_STEPS = 200000


class _FrStepCap(Exception):
    pass


def _fr_automorphisms(ch: dict) -> List[Dict[int, int]]:
    """Up to _AUTO_CAP automorphisms of the organic SMILES graph ((element, nH) labels), identity first.

    VF2 (networkx) with a STEP budget (_AUTO_STEPS feasibility checks), not a time limit, so the result is the
    same on every machine; the search stops with what it has found when the budget is spent."""
    if "autos" in ch:
        return ch["autos"]
    import networkx as nx
    from networkx.algorithms import isomorphism as iso
    sg = ch["sg"]
    ident = {s: s for s in ch["Gs_nodes"]}
    autos = [ident]
    G = nx.Graph()
    for s in ch["Gs_nodes"]:
        G.add_node(s, lab=(sg["el"][s], sg["nh"][s]))
    G.add_edges_from(ch["Gs_edges"])
    steps = [0]

    def nm(a, b):
        steps[0] += 1
        if steps[0] > _AUTO_STEPS:
            raise _FrStepCap()
        return a["lab"] == b["lab"]
    try:
        for au in iso.GraphMatcher(G, G, node_match=nm).isomorphisms_iter():
            if au != ident:
                autos.append(dict(au))
            if len(autos) >= _AUTO_CAP:
                break
    except _FrStepCap:
        pass
    except Exception:
        autos = [ident]
    ch["autos"] = autos
    return autos


def order_by_realism(scores: Sequence[Optional[float]]) -> List[int]:
    """Indices best first: mappable before unmappable, then score, then the emitted index."""
    return sorted(range(len(scores)),
                  key=lambda i: (scores[i] is None, scores[i] if scores[i] is not None else 0.0, i))


def rank_isomers_realism(isomers, smiles: Optional[str] = None):
    """The emitted ensemble sorted best first (DELFIN_FRAME_REALISM_RANK).  A verified permutation: never
    drops, adds or alters a frame; any failure returns the input order."""
    if not isomers or not rank_enabled():
        return isomers
    try:
        xyzs = [it[0] if isinstance(it, (tuple, list)) else it for it in isomers]
        sc = score_frames(xyzs, smiles)
        if sc is None or len(sc) != len(isomers):
            return isomers
        order = order_by_realism(sc)
        if sorted(order) != list(range(len(isomers))):
            return isomers
        return [isomers[i] for i in order]
    except Exception:
        return isomers
