# Curated adversarial eye package: OFF vs ON vs free-ligand donor cavity.
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
CAP = []
_ofc = AC.assemble_from_config
def _wfc(metal, geom, config, ligands, **k):
    r = _ofc(metal, geom, config, ligands, **k)
    if r is not None and k.get("n_frames", 1) == 1:
        try:
            s, P, don = r; CAP.append((geom, list(s), np.array(P, float), list(don)))
        except Exception: pass
    return r
AC.assemble_from_config = _wfc; CB.AC.assemble_from_config = _wfc
def cshm(geom, P, don):
    M = P[0]; v = [P[i] - M for i in don]
    try: return float(PLY.cshm(v, geom))
    except Exception: return float('nan')
def xyz(syms, P, comment):
    out = [str(len(syms)), comment]
    for s, p in zip(syms, P): out.append(f"{s:2} {p[0]:10.4f} {p[1]:10.4f} {p[2]:10.4f}")
    return "\n".join(out) + "\n"
def build(smi, on):
    if on: os.environ["DELFIN_FFFREE_RIGID_LIGAND_SEAT"] = "1"
    else: os.environ.pop("DELFIN_FFFREE_RIGID_LIGAND_SEAT", None)
    CAP.clear()
    try: r = CB._fffree_isomers(smi)
    except Exception: r = None
    if not CAP: return None
    geom, s, P, don = CAP[0]
    return geom, s, P, don, (0 if r is None else len(r))

sample = {d["ref"]: d for d in json.load(open("/tmp/rigid_sample.json"))}
data = {d["ref"]: d for d in json.load(open("/tmp/compare_result.json"))}
sample["ADELED"] = {"ref": "ADELED", "metal": "Pt", "cn": 4, "maxdent": 4,
    "smi": "CC1(C)C2=CNC3=[C]2[Pt-2]24[C]5=C(NC=C5C(C)(C)C5=[N+]2C(=C3C2=C(F)C(F)=C(F)C(F)=C2F)C=C5)C(C2=C(F)C(F)=C(F)C(F)=C2F)=C2C=CC1=[N+]24"}
CURATED = {
  "WIN_clear":  ["ADELED", "VIWKOC", "JEZNUB", "PIXGOS", "JEPDAJ", "BEDLIG"],
  "REGRESSION": ["HOKCUH", "SOTFIQ", "YOMLAN", "AJIDOM"],
  "LIGAND_FID": ["XELCUM", "JIRDOE", "GORJOO"],
}
OUT = "/tmp/delfin_k4/MANTA/EYE_REVIEW/KAPPA4-LIGINVARIANT-2026-06-24"
os.makedirs(OUT, exist_ok=True)
rm = ["# κ4 / rigid-polydentate LIGAND-INVARIANT cavity seating — adversarial eye package",
      "# Flag DELFIN_FFFREE_RIGID_LIGAND_SEAT (default OFF). 2026-06-24 iter36.",
      "#",
      "# OFF = current (per-donor rescale, polyhedron forced toward ideal).",
      "# ON  = cavity+rigid (ligand internal geometry preserved, polyhedron EMERGENT).",
      "# CShM = deviation from IDEAL polyhedron (LOWER=more ideal, but a rigid ligand's",
      "#        REAL crystal polyhedron is often distorted -> high CShM can be MORE physical).",
      "# dd   = built donor-donor cavity RMSD vs the FREE ligand (LOWER = ligand preserved).",
      "# QUESTION FOR THE EYE: which frame looks like the real coordination + intact ligand?",
      "#"]
for group, refs in CURATED.items():
    rm.append(f"\n## {group}")
    for ref in refs:
        if ref not in sample:
            rm.append(f"{ref}: (not in sample)"); continue
        e = sample[ref]; smi = e["smi"]
        for on, tag in ((False, "OFF"), (True, "ON")):
            b = build(smi, on)
            if b is None:
                rm.append(f"{ref}_{tag}: build None (legacy)"); continue
            geom, s, P, don, niso = b
            c = cshm(geom, P, don)
            open(f"{OUT}/{ref}_{tag}.xyz", "w").write(xyz(s, P, f"{ref} {tag} geom={geom} CShM={c:.2f} niso={niso}"))
        d = data.get(ref, {})
        def gv(t, key): return d.get(t, {}).get(key, "-")
        rm.append(f"{ref} {e['metal']} CN{e['cn']} κ{e['maxdent']}: "
                  f"CShM OFF={gv('off','cshm')} ON={gv('on','cshm')} | "
                  f"dd OFF={gv('off','dd_rmsd')} ON={gv('on','dd_rmsd')} | "
                  f"niso OFF={gv('off','n_iso')} ON={gv('on','n_iso')}")
open(f"{OUT}/README.md", "w").write("\n".join(str(x) for x in rm) + "\n")
print("eye package ->", OUT)
import glob
print("xyz files:", len(glob.glob(f"{OUT}/*.xyz")))
