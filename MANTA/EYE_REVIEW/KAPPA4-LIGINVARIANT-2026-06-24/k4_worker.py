# single-structure worker: prints one JSON line with OFF+ON composite metrics
import sys; sys.path.insert(0, "/tmp/delfin_k4")
import os, json, numpy as np
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
for k, v in {"DELFIN_FFFREE_KAPPA4": "1", "DELFIN_FFFREE_HIGHCN": "1",
             "DELFIN_FFFREE_CHELATE_BACKBONE": "1", "DELFIN_FFFREE_PLANAR_MER": "1",
             "DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE": "1"}.items():
    os.environ[k] = v
from delfin.fffree import assemble_complex as AC
from delfin.fffree import converter_backend as CB
from delfin.fffree import polyhedra as PLY
from delfin.fffree import decompose as DEC
from delfin import _bond_decollapse as _bd
from rdkit import Chem
from rdkit.Chem import AllChem

CAP = []
_ofc = AC.assemble_from_config
def _wfc(metal, geom, config, ligands, **k):
    r = _ofc(metal, geom, config, ligands, **k)
    if r is not None and k.get("n_frames", 1) == 1:
        try:
            s, P, don = r
            CAP.append((geom, list(s), np.array(P, float), list(don)))
        except Exception:
            pass
    return r
AC.assemble_from_config = _wfc
CB.AC.assemble_from_config = _wfc

def cshm(geom, P, don):
    M = P[0]; v = [P[i] - M for i in don]
    try:
        return float(PLY.cshm(v, geom))
    except Exception:
        return float('nan')

def comp_donors(syms, P, donors):
    n = len(syms); bonds = _bd._geometric_bonds(syms, P); adj = {i: set() for i in range(n)}
    mi = next((i for i in range(n) if _bd._is_metal(syms[i])), 0)
    for i, j in bonds:
        if i == mi or j == mi:
            continue
        adj[i].add(j); adj[j].add(i)
    seen = set(); comps = []
    for s in range(n):
        if s == mi or syms[s] == "H" or s in seen:
            continue
        st = [s]; c = set()
        while st:
            x = st.pop()
            if x in c:
                continue
            c.add(x)
            for y in adj[x]:
                if syms[y] != "H" and y not in c:
                    st.append(y)
        seen |= c; comps.append(c)
    if not comps:
        return []
    big = max(comps, key=len)
    return [d for d in donors if d in big]

def sdd(P, idxs):
    return sorted(float(np.linalg.norm(P[a] - P[b])) for x, a in enumerate(idxs) for b in idxs[x + 1:])

def free_dd(smi):
    d = DEC.decompose(smi)
    if d is None:
        return None
    lg = max(d["ligands"], key=lambda l: l["denticity"])
    if lg["denticity"] < 3:
        return None
    mh = Chem.AddHs(lg["mol"])
    if AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE, useRandomCoords=False) != 0:
        if AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE, useRandomCoords=True) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mh)
    except Exception:
        pass
    Pf = np.array(mh.GetConformer().GetPositions(), float)
    return sdd(Pf, lg["donor_local_idxs"])

def run(smi, on):
    if on:
        os.environ["DELFIN_FFFREE_RIGID_LIGAND_SEAT"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_RIGID_LIGAND_SEAT", None)
    CAP.clear()
    try:
        r = CB._fffree_isomers(smi)
    except Exception:
        r = None
    o = {"n_iso": 0 if r is None else len(r)}
    if CAP:
        geom, s, P, don = CAP[0]
        o["geom"] = geom; o["cshm"] = cshm(geom, P, don)
        pd = comp_donors(s, P, don); o["built_dd"] = sdd(P, pd)
    return o

smi = sys.argv[1]
out = {"off": run(smi, False), "on": run(smi, True)}
fdd = free_dd(smi)
if fdd:
    out["free_dd"] = fdd
    for t in ("off", "on"):
        b = out[t].get("built_dd")
        if b and len(b) == len(fdd):
            out[t]["dd_rmsd"] = float(np.sqrt(np.mean((np.array(b) - np.array(fdd)) ** 2)))
print(json.dumps(out))
