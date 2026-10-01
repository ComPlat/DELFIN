#!/usr/bin/env python
"""Tool-neutral spec (DELFIN split_complex_smiles, specs/*.jsonl) -> epic-MACE input.

Runs in the epic-MACE environment (python 3.7, rdkit 2020.09, epic-mace); no DELFIN import.  Importable (mace_worker.py
uses it) and runnable as the conversion census:

    mace_convert.py SPECS.jsonl LIST.txt [--geoms paper|extended] [--csv OUT.csv] [--check]

Input choice: the spec, i.e. the SAME metal / ligand SMILES / donor atoms (coordList) that
Architector and molSimplify get (bench_convert.py).  epic-MACE wants dative SMILES (donor -> metal,
donor atoms carrying atom-map numbers, one central atom); the spec's free-ligand SMILES are written
into exactly that form with the authors' own API (mace.ComplexFromLigands: ligand SMILES with
mapped donors + central-atom SMILES).  Dative SMILES written by other programs are NOT used:
for hapto systems they can differ in composition from the input SMILES (H on pi-bound C), and
they write a hapto ligand as n separate dative bonds, which MACE cannot read either (it needs a
centroid) -- see docs/CONSTRUCTION_BATCH.md.

Donor sites.  Donors of one ligand that are bonded to each other form one hapto group (eta-n,
n >= 2); every other donor is one site.  A hapto group becomes one MACE centroid dummy [*:1]
(epic-MACE >= "0.6.0" = GitHub master, docs/source/haptic_ligands.rst), encoded as the authors
document it:
  * 5- or 6-membered all-carbon ring ("anchor"): the centroid is bonded to ONE ring carbon; MACE
    expands it to the whole ring.  The anchor gets the free valence by dropping its formal charge
    or a radical, else one ring double bond at the anchor becomes single (the partner keeps a
    radical).  H counts are frozen, so the composition is exactly the spec's.
  * every other group ("star", eta2 alkene / eta3 allyl / eta4 diene / heteroatom rings ...): the
    centroid is bonded to every group atom and the bonds INSIDE the group are removed, as in the
    authors' eta2-ethylene "[*:4]([CH2])[CH2]" and eta3-allyl "[*:1]([CH2])([CH])[CH2]" examples.
The centroid is an X atom in MACE's xyz output; mace_worker.py drops it.

Geometry from the number of sites (not the spec's CN, which counts hapto atoms one by one):
  paper    (default; the geometries of the JCTC 2024 paper, v0.5.0):  6 -> OH,  4 -> SP
  extended (GitHub master adds TET/SPY/TBP/SAN):  6 -> OH,  5 -> SPY + TBP,  4 -> SP + TET,
           2 hapto centroids -> SAN.  Several geometries per system are all built (like
           molSimplify's sqp + thd).
Anything else is `not_expressible` (not a tool failure), with the reason.  The central atom is
written with the spec's oxidation state as formal charge when 0..8 (as in MACE's examples,
"[Ru+2]"); MACE uses it only for UFF typing.  Unlike Architector/molSimplify, MACE needs neither
the oxidation state nor a spin, so oxidation states outside 0..8 and radical ligands stay
expressible (recorded as flags).
"""
import argparse
import collections
import csv
import json
import sys

from rdkit import Chem

GEOMS = {"paper": {6: ["OH"], 4: ["SP"]},
         "extended": {6: ["OH"], 5: ["SPY", "TBP"], 4: ["SP", "TET"]}}

_DBLOCK_GROUP = {}
for _row in (["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"],
             ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"],
             ["Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]):
    for _i, _el in enumerate(_row):
        _DBLOCK_GROUP[_el] = _i + 3
_DBLOCK_GROUP["La"] = 3


def mace_d_count(metal, ox):
    g = _DBLOCK_GROUP.get(metal)
    if g is None or ox is None or not 0 <= ox <= g:
        return None
    return g - ox


def mace_donor_groups(mol, coord):
    """Connected components of the donor atoms in the ligand graph (sorted)."""
    coord = sorted(set(coord))
    cs = set(coord)
    seen, groups = set(), []
    for d in coord:
        if d in seen:
            continue
        comp, stack = [], [d]
        seen.add(d)
        while stack:
            a = stack.pop()
            comp.append(a)
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                j = n.GetIdx()
                if j in cs and j not in seen:
                    seen.add(j)
                    stack.append(j)
        groups.append(sorted(comp))
    return groups


def _freeze_h(mol):
    for a in mol.GetAtoms():
        a.SetNumExplicitHs(a.GetTotalNumHs())
        a.SetNoImplicit(True)


def _ring_of(mol, group):
    if len(group) not in (5, 6):
        return None
    gs = set(group)
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in group):
        return None
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == len(gs) and set(ring) == gs:
            return list(ring)
    return None


def _free_valence_for_anchor(rw, ring):
    """Pick the anchor of a carbon pi ring and give it one free valence -> anchor index."""
    rs = set(ring)
    # Adjacent ring carbons that both carry radicals (the spec writes some pi-bound arenes as
    # H-less [C] carbenes, which RDKit 2020.09 types SP3 -- MACE then does not see a pi ring)
    # are joined into ring double bonds: same atoms, same H, only the Lewis structure changes.
    n = len(ring)
    for k in range(n):
        i, j = ring[k], ring[(k + 1) % n]
        ai, aj = rw.GetAtomWithIdx(i), rw.GetAtomWithIdx(j)
        b = rw.GetBondBetweenAtoms(i, j)
        if (b is not None and b.GetBondType() == Chem.BondType.SINGLE
                and ai.GetNumRadicalElectrons() > 0 and aj.GetNumRadicalElectrons() > 0
                and not any(x.GetBondType() == Chem.BondType.DOUBLE
                            for x in list(ai.GetBonds()) + list(aj.GetBonds()))):
            b.SetBondType(Chem.BondType.DOUBLE)
            ai.SetNumRadicalElectrons(ai.GetNumRadicalElectrons() - 1)
            aj.SetNumRadicalElectrons(aj.GetNumRadicalElectrons() - 1)
    for i in ring:  # a radical ring carbon left without a double bond (SP3 for RDKit): anchor
        a = rw.GetAtomWithIdx(i)
        if a.GetNumRadicalElectrons() > 0 and not any(
                x.GetBondType() == Chem.BondType.DOUBLE for x in a.GetBonds()):
            a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() - 1)
            return i
    for i in ring:  # charged or radical ring carbon: neutralise / use the radical
        a = rw.GetAtomWithIdx(i)
        if a.GetFormalCharge() != 0 or a.GetNumRadicalElectrons() > 0:
            if a.GetFormalCharge() != 0:
                a.SetFormalCharge(0)
            else:
                a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() - 1)
            return i
    for i in ring:  # a ring double bond at the anchor becomes single, partner keeps a radical
        a = rw.GetAtomWithIdx(i)
        for b in a.GetBonds():
            j = b.GetOtherAtomIdx(i)
            if j in rs and b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetBondType(Chem.BondType.SINGLE)
                p = rw.GetAtomWithIdx(j)
                p.SetNumRadicalElectrons(p.GetNumRadicalElectrons() + 1)
                return i
    return ring[0]


def mace_ligand(lig):
    """Spec ligand -> (MACE ligand SMILES with mapped donors, [site info]) or (None, reason)."""
    m = Chem.MolFromSmiles(lig["smiles"])
    if m is None:
        return None, "rdkit2020_unparseable_ligand"
    try:
        Chem.Kekulize(m, clearAromaticFlags=True)
    except Exception:
        return None, "rdkit2020_kekulize_failed"
    _freeze_h(m)
    groups = mace_donor_groups(m, lig["coordList"])
    rw = Chem.RWMol(m)
    sites = []
    for g in groups:
        if len(g) == 1:
            rw.GetAtomWithIdx(g[0]).SetAtomMapNum(1)
            sites.append({"kind": "atom", "elem": m.GetAtomWithIdx(g[0]).GetSymbol()})
            continue
        ring = _ring_of(m, g)
        star = rw.AddAtom(Chem.Atom(0))
        rw.GetAtomWithIdx(star).SetAtomMapNum(1)
        if ring is not None:
            anc = _free_valence_for_anchor(rw, ring)
            rw.AddBond(star, anc, Chem.BondType.SINGLE)
            sites.append({"kind": "hapto", "eta": len(g), "enc": "anchor"})
        else:
            gs = set(g)
            for b in list(m.GetBonds()):
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                if i in gs and j in gs:
                    rw.RemoveBond(i, j)
            for i in g:
                rw.AddBond(star, i, Chem.BondType.SINGLE)
            sites.append({"kind": "hapto", "eta": len(g), "enc": "star"})
    mol = rw.GetMol()
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None, "sanitize_failed:%s" % str(e)[:60]
    return Chem.MolToSmiles(mol), sites


def mace_ca_smiles(metal, ox):
    if ox is None or not 0 <= ox <= 8 or ox == 0:
        return "[%s]" % metal
    return "[%s+%d]" % (metal, ox)


def mace_convert(rec, geoms="paper"):
    """-> (job dict, None) or (job or None, "not_expressible:<reason>").

    job = {"geoms": [...], "ligands": [...], "CA": "...", "info": {...}}"""
    if rec.get("status") != "ok":
        return None, "not_expressible:spec_%s" % (rec.get("bail") or rec.get("status"))
    ligs, sites = [], []
    for l in rec["ligands"]:
        s, info = mace_ligand(l)
        if s is None:
            return None, "not_expressible:%s" % info
        ligs.append(s)
        sites += info
    n_sites = len(sites)
    n_hapto = sum(1 for s in sites if s["kind"] == "hapto")
    ox = rec.get("metal_ox")
    d = mace_d_count(rec["metal"], ox)
    info = {"n_sites": n_sites, "n_hapto_sites": n_hapto,
            "hapto": [(s["eta"], s["enc"]) for s in sites if s["kind"] == "hapto"],
            "spec_cn": rec.get("cn"), "metal": rec["metal"], "metal_ox": ox, "d_count": d,
            "flags": sorted({f for l in rec["ligands"] for f in l["flags"]}),
            "ox_outside_0_8": ox is None or not 0 <= ox <= 8}
    if n_sites == 4:
        # MACE's paper geometries have no tetrahedron: CN4 goes to SP.  Record where that is
        # chemically unlikely (d10 / d0 centres with four donors are usually tetrahedral).
        info["cn4_tetrahedral_expected"] = d in (0, 10)
    g = list(GEOMS[geoms].get(n_sites, []))
    if geoms == "extended" and n_sites == 2 and n_hapto == 2:
        g = ["SAN"]
    if not g:  # the job (with geoms = []) is still returned for the census
        return ({"geoms": [], "ligands": ligs, "CA": mace_ca_smiles(rec["metal"], ox),
                 "info": info}, "not_expressible:sites_%d_no_%s_geometry" % (n_sites, geoms))
    return {"geoms": g, "ligands": ligs, "CA": mace_ca_smiles(rec["metal"], ox),
            "info": info}, None


def mace_check_job(job):
    """Build the Complex objects (no 3D) -> (None, info) or ("fail:<why>", info).

    Verifies that MACE reads every hapto centroid with the intended hapticity."""
    import mace
    want = sorted(e for e, _ in job["info"]["hapto"])
    out = {}
    for geom in job["geoms"]:
        try:
            X = mace.ComplexFromLigands(job["ligands"], job["CA"], geom)
        except Exception as e:
            return "fail:mace_init:%s:%s" % (type(e).__name__, str(e).splitlines()[0][:80]), out
        got = sorted(len(v) for v in getattr(X, "_haptic_DAs", {}).values())
        out[geom] = {"eta_detected": got}
        if got != want:
            out[geom]["eta_mismatch"] = True
    return None, out


def mace_load_specs(path, want=None):
    out = {}
    for l in open(path):
        r = json.loads(l)
        if want is None or r["refcode"] in want:
            out[r["refcode"]] = r
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("specs")
    ap.add_argument("pool", help="REFCODE;SMILES or REFCODE per line")
    ap.add_argument("--geoms", default="paper", choices=sorted(GEOMS))
    ap.add_argument("--csv", default=None)
    ap.add_argument("--check", action="store_true",
                    help="also build MACE's Complex objects (catches inputs MACE rejects)")
    a = ap.parse_args()
    refs = [l.strip().split(";", 1)[0].split("|", 1)[0] for l in open(a.pool) if l.strip()]
    specs = mace_load_specs(a.specs, set(refs))
    cnt, sub = collections.Counter(), collections.Counter()
    rows = []
    for i, r in enumerate(refs):
        rec = specs.get(r)
        if rec is None:
            cls, why, job = "not_expressible", "no_spec", None
        else:
            job, err = mace_convert(rec, a.geoms)
            if err:
                cls, why = "not_expressible", err.split(":", 1)[1]
            else:
                cls, why = "ok", ""
                if a.check:
                    bad, chk = mace_check_job(job)
                    if bad:
                        cls, why = "fail", bad[5:]
                    elif any(v.get("eta_mismatch") for v in chk.values()):
                        why = "eta_mismatch"
        cnt[cls] += 1
        key = why.split(":")[0] if why else ""
        if why.startswith("sites_"):
            key = why
        sub[(cls, key)] += 1
        inf = job["info"] if job else {}
        rows.append([r, rec.get("metal") if rec else "", rec.get("cn") if rec else "",
                     rec.get("metal_ox") if rec else "", bool(rec and rec.get("hapto")),
                     inf.get("n_sites", ""), inf.get("n_hapto_sites", ""),
                     "|".join(job["geoms"]) if job else "", cls, why,
                     inf.get("cn4_tetrahedral_expected", ""),
                     "|".join("eta%d-%s" % t for t in inf.get("hapto", []))])
        if (i + 1) % 2000 == 0:
            print("...", i + 1, dict(cnt), file=sys.stderr, flush=True)
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["refcode", "metal", "spec_cn", "metal_ox", "hapto", "n_sites",
                        "n_hapto_sites", "geoms", "class", "reason",
                        "cn4_tetrahedral_expected", "hapto_encoding"])
            w.writerows(rows)
    print("n =", len(refs), " geoms =", a.geoms, " check =", a.check)
    print("classes:", dict(cnt.most_common()))
    for (c, k), v in sorted(sub.items(), key=lambda t: (t[0][0], -t[1])):
        print("  %-16s %-40s %d" % (c, k or "-", v))
    hap = [r for r in rows if r[4]]
    print("hapto systems: %d, class ok: %d" % (len(hap), sum(1 for r in hap if r[8] == "ok")))
    cn4 = [r for r in rows if r[8] == "ok" and r[5] == 4]
    print("4-site ok: %d, of them tetrahedral expected (d0/d10): %d"
          % (len(cn4), sum(1 for r in cn4 if r[10] is True)))
