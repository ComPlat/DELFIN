#!/usr/bin/env python
"""Tool-neutral spec (delfin.cluster_bench.prepare) -> Architector input / molSimplify inputs.

Runs in the Architector/molSimplify environment; no DELFIN import.  Importable (bench_worker.py uses it) and runnable as a CLI that
reports conversion statistics for a pool:

    bench_convert.py SPECS.jsonl LIST.txt  [--csv conversion_bench.csv]   (LIST: ID;SMILES or ID|SMILES)

Architector: {'core': {'metal', 'coreCN'}, 'ligands': [{'smiles', 'coordList'(0-based)}],
              'parameters': {'metal_ox', ...}}; metal_spin = Architector's own default rule
              (computed here, see bench_default_spin).
molSimplify: one input per molSimplify geometry of that CN (its coordinations.dict), ligands as
             SMILES strings + '-smicat' (1-based), '-ff uff -ffoption BA' (its standard
             recipe), '-spinmultiplicity' = Architector's default spin rule (same electronic
             assumption for both tools).  molSimplify's '-isomers' only works for ligands of its
             own dictionary (KeyError on SMILES), so it emits one structure per geometry.
molSimplify pitfall handled here: lig_load() treats ANY ligand string that is within
SequenceMatcher ratio > 0.6 of a ligands.dict name as that dictionary ligand ("typo" rescue),
e.g. methanol 'CO' -> carbonyl.  Colliding SMILES are rewritten to an equivalent spelling
(rooted / explicit-H) that does not collide; if none exists the system is a conversion failure.
"""
import argparse
import collections
import csv
import difflib
import json

# rdkit before openbabel/molSimplify: otherwise a foreign libstdc++ (LD_LIBRARY_PATH, e.g.
# /opt/orca/lib) gets loaded first and rdkit fails with CXXABI_1.3.15 not found.
from rdkit import Chem  # noqa: F401,E402

MS_GEOMS = {2: ["li"], 3: ["tpl"], 4: ["sqp", "thd"], 5: ["spy", "tbp"], 6: ["oct", "tpr"],
            7: ["pbp"], 8: ["sqap", "tdhd"]}

_LICORES = None


def bench_ms_licores():
    global _LICORES
    if _LICORES is None:
        from molSimplify.Scripts.io import getlicores
        _LICORES = getlicores()
    return _LICORES


def bench_ms_collides(s):
    keys = list(bench_ms_licores().keys())
    if s in keys:
        return True
    return max(difflib.SequenceMatcher(None, s, k).ratio() for k in keys) > 0.6


def bench_smiles_variants(smi, coord):
    """Equivalent spellings of a ligand SMILES with the donor positions re-mapped."""
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return
    for i, a in enumerate(m.GetAtoms()):
        a.SetIntProp("_o", i)
    Chem.Kekulize(m, clearAromaticFlags=True)
    # molSimplify drops the implicit H of very small SMILES ligands (water "O" -> bare O, "CO"
    # -> C,O): write them with explicit H first.
    for allH in ((True, False) if m.GetNumAtoms() <= 3 else (False, True)):
        for root in range(-1, m.GetNumAtoms()):
            s = Chem.MolToSmiles(m, kekuleSmiles=True, canonical=True, rootedAtAtom=root,
                                 allHsExplicit=allH)
            order = list(m.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
            pos = {old: i for i, old in enumerate(order)}
            yield s, sorted(pos[c] for c in coord)


def bench_default_spin(metal, ox):
    """Architector's own default: metal_charge_dict/metal_spin_dict, else mendeleev aufbau."""
    from architector import io_ptable
    import mendeleev
    if ox == io_ptable.metal_charge_dict.get(metal, 100):
        return int(io_ptable.metal_spin_dict[metal])
    # Architector 0.0.10 calls mendeleev.__dict__[metal] here, which mendeleev 1.3 no longer
    # provides (KeyError for Fe, Co, ... whenever ox != its reference ox).  Same rule, new API;
    # the result is passed to Architector explicitly as metal_spin.
    return int(mendeleev.element(metal).ec.ionize(ox).unpaired_electrons())


def bench_common_check(rec):
    if rec.get("status") != "ok":
        return "spec_not_ok"
    if any("radical_left" in l["flags"] for l in rec["ligands"]):
        return "radical_in_ligand"
    if not 0 <= rec["metal_ox"] <= 8:
        return "ox_out_of_range(%d)" % rec["metal_ox"]
    return None


def bench_ob_parses(smi):
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    m = ob.OBMol()
    return bool(conv.ReadString(m, smi)) and m.NumAtoms() > 0


def bench_to_architector(rec, params=None):
    """-> (input_dict, None) or (None, failure_category)."""
    bad = bench_common_check(rec)
    if bad:
        return None, bad
    from architector import io_ptable
    from architector.io_core import Geometries
    if rec["metal"] not in io_ptable.all_metals:
        return None, "arch_metal_unsupported"
    if rec["cn"] not in Geometries().cn_geo_dict:
        return None, "arch_cn_unsupported"
    ligs = []
    for l in rec["ligands"]:
        if not bench_ob_parses(l["smiles"]):
            return None, "ob_smiles_unparseable"
        ligs.append({"smiles": l["smiles"], "coordList": list(l["coordList"])})
    try:
        spin = bench_default_spin(rec["metal"], int(rec["metal_ox"]))
    except Exception:
        return None, "spin_default_failed"
    p = {"metal_ox": int(rec["metal_ox"]), "metal_spin": spin}
    p.update(params or {})
    return {"core": {"metal": rec["metal"], "coreCN": int(rec["cn"])},
            "ligands": ligs, "parameters": p}, None


def bench_to_molsimplify(rec):
    """-> (list of (geometry, input_dict without -rundir/-name), None) or (None, category)."""
    bad = bench_common_check(rec)
    if bad:
        return None, bad
    if rec["cn"] not in MS_GEOMS:
        return None, "ms_cn_unsupported"
    if max(l["denticity"] for l in rec["ligands"]) > 6:
        return None, "ms_denticity_gt6"
    smis, cats = [], []
    # molSimplify fills core sites in ligand order; a polydentate listed after monodentates runs
    # out of adjacent sites ("No more connecting points" -> KeyError in distgeom).  Its examples
    # list the highest denticity first -> stable sort by denticity, descending.
    for l in sorted(rec["ligands"], key=lambda l: -l["denticity"]):
        if not bench_ob_parses(l["smiles"]):
            return None, "ob_smiles_unparseable"
        chosen = None
        for s, c in bench_smiles_variants(l["smiles"], l["coordList"]):
            if not bench_ms_collides(s):
                chosen = (s, c)
                break
        if chosen is None:
            return None, "ms_ligand_name_collision"
        smis.append(chosen[0])
        cats.append(chosen[1])
    try:
        spin = bench_default_spin(rec["metal"], int(rec["metal_ox"]))
    except Exception:
        return None, "spin_default_failed"
    base = {"-core": rec["metal"],
            "-lig": ",".join(smis),
            "-ligocc": ",".join("1" for _ in smis),
            "-smicat": "[" + ",".join("[" + ",".join(str(c + 1) for c in cc) + "]"
                                      for cc in cats) + "]",
            "-coord": str(rec["cn"]),
            "-oxstate": str(int(rec["metal_ox"])),
            "-spinmultiplicity": str(spin + 1),
            "-ff": "uff", "-ffoption": "BA",
            "-keepHs": ",".join("yes" for _ in smis)}  # SMILES protonation is final; "auto" strips donor H (aqua -> hydroxo)
    out = []
    for g in MS_GEOMS[rec["cn"]]:
        d = dict(base)
        d["-geometry"] = g
        out.append((g, d))
    return out, None


def bench_load_specs(path):
    return {json.loads(l)["refcode"]: json.loads(l) for l in open(path)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("specs")
    ap.add_argument("pool")
    ap.add_argument("--csv", default="conversion_bench.csv")
    a = ap.parse_args()
    specs = bench_load_specs(a.specs)
    refs = [l.strip().split(";", 1)[0].split("|", 1)[0] for l in open(a.pool) if l.strip()]
    ca, cm = collections.Counter(), collections.Counter()
    rows = []
    for r in refs:
        rec = specs[r]
        ai, ae = bench_to_architector(rec)
        mi, me = bench_to_molsimplify(rec)
        ca[ae or "ok"] += 1
        cm[me or "ok"] += 1
        flags = sorted({f for l in rec.get("ligands", []) for f in l["flags"]})
        rows.append([r, rec.get("metal"), rec.get("cn"), rec.get("metal_ox"),
                     ae or "ok", me or "ok", len(mi) if mi else 0, "|".join(flags)])
    with open(a.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["refcode", "metal", "cn", "metal_ox", "architector", "molsimplify",
                    "ms_geometries", "ligand_flags"])
        w.writerows(rows)
    print("n =", len(refs))
    print("Architector:", dict(ca.most_common()))
    print("molSimplify:", dict(cm.most_common()))
    fl = collections.Counter(f for row in rows for f in row[7].split("|") if f)
    print("ligand flags (systems):", dict(fl))
