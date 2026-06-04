#!/usr/bin/env python3
"""Phase A+B analysis of b00f9a0 voll-pool isocoverage data."""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

JSONL = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL_isocoverage.jsonl")
OUT_DIR = Path("/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a6a0c001c3fc98a82/forensik")
OUT_DIR.mkdir(exist_ok=True)

# ---------- Phase A: load + per-SMILES coverage distribution ----------
records = []
with JSONL.open() as f:
    for line in f:
        line = line.strip()
        if line:
            records.append(json.loads(line))

print(f"Total isocov records: {len(records)}")

# Each record IS one SMILES (per-SMILES coverage)
print(f"Unique smiles_id values: {len(set(r.get('smiles_id', '') for r in records))}")

# Coverage histogram bins
bins = [(0, 0), (1, 25), (26, 50), (51, 75), (76, 99), (100, 100)]
hist = Counter()
for r in records:
    cov = r.get("coverage_pct", 0)
    if cov == 0:
        hist["0"] += 1
    elif cov == 100:
        hist["100"] += 1
    elif cov <= 25:
        hist["1-25"] += 1
    elif cov <= 50:
        hist["26-50"] += 1
    elif cov <= 75:
        hist["51-75"] += 1
    else:
        hist["76-99"] += 1

print("\n=== PHASE A: Per-SMILES Coverage Histogram ===")
for k in ["0", "1-25", "26-50", "51-75", "76-99", "100"]:
    n = hist.get(k, 0)
    pct = 100 * n / max(1, len(records))
    print(f"  {k:>6}%  : {n:>5}  ({pct:5.1f}%)")

# Coverage by class
print("\n=== Per-Class Coverage ===")
class_hist = defaultdict(lambda: Counter())
for r in records:
    cov = r.get("coverage_pct", 0)
    cls = r.get("class", "UNKNOWN")
    if cov == 0:
        class_hist[cls]["0"] += 1
    elif cov == 100:
        class_hist[cls]["100"] += 1
    else:
        class_hist[cls]["partial"] += 1
    class_hist[cls]["total"] += 1
    class_hist[cls]["sum_cov"] += cov

for cls, c in sorted(class_hist.items()):
    n = c["total"]
    mean = c["sum_cov"] / n
    print(f"  {cls:>15}: n={n:>5}  0%={c['0']:>4} ({100*c['0']/n:4.1f}%)  100%={c['100']:>4} ({100*c['100']/n:4.1f}%)  mean={mean:5.1f}%")

# ---------- Phase B: cluster 0% SMILES ----------
print("\n=== PHASE B: Clustering 0%-coverage SMILES ===")
zero_records = [r for r in records if r.get("coverage_pct", 0) == 0]
print(f"Total 0%-coverage records: {len(zero_records)}")

# Save for further analysis
with (OUT_DIR / "zero_coverage_records.jsonl").open("w") as f:
    for r in zero_records:
        f.write(json.dumps(r) + "\n")

# By metal
metal_counts = Counter(r.get("metal", "?") for r in zero_records)
print("\n  By metal:")
for m, n in metal_counts.most_common(20):
    print(f"    {m:>4}: {n}")

# Metal classification (3d/4d/5d/lanthanide/actinide)
PERIOD_3D = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
PERIOD_4D = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
PERIOD_5D = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}
LANTH = {"Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"}
ACT = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"}
MAIN = {"Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba", "Al", "Ga", "In", "Tl", "Si", "Ge", "Sn", "Pb", "Sb", "Bi"}

def metal_class(m):
    if m in PERIOD_3D: return "3d"
    if m in PERIOD_4D: return "4d"
    if m in PERIOD_5D: return "5d"
    if m in LANTH: return "4f-lanth"
    if m in ACT: return "5f-act"
    if m in MAIN: return "main-group"
    return "other"

metal_class_counts = Counter(metal_class(r.get("metal", "?")) for r in zero_records)
print("\n  By metal class:")
for mc, n in metal_class_counts.most_common():
    print(f"    {mc:>10}: {n}")

# By CN
cn_counts = Counter(r.get("cn", "?") for r in zero_records)
print("\n  By CN:")
for cn, n in sorted(cn_counts.items(), key=lambda x: (isinstance(x[0], str), x[0])):
    print(f"    CN{cn}: {n}")

# By polyhedron_class (parses denticity pattern)
print("\n  By polyhedron_class (top 25 patterns):")
poly_counts = Counter(r.get("polyhedron_class", "?") for r in zero_records)
for p, n in poly_counts.most_common(25):
    print(f"    {n:>4}  {p}")

# Extract dentate pattern (e.g. tridentate+bis-bidentate-NOO+NO+NO -> tridentate+bidentate*2)
def dentate_pattern(pc):
    if not pc: return "?"
    # split on +
    parts = pc.split("+")
    # extract leading word (or e.g. "tris-monodentate")
    out = []
    for p in parts:
        p = p.split("-")[0]
        # also handle bis/tris prefixes
        if p.startswith("bis"): p = "bis-bidentate" if "bidentate" in pc else p
        out.append(p)
    return "+".join(out[:3]) + ("+..." if len(parts) > 3 else "")

dpat = Counter()
for r in zero_records:
    pc = r.get("polyhedron_class", "")
    # simpler: count denticity words
    dents = re.findall(r"(monodentate|bidentate|tridentate|tetradentate|pentadentate|hexadentate)", pc)
    if not dents:
        dpat["NO_DENT"] += 1
    else:
        # canonical = sorted denticity tuple
        cnt = Counter(dents)
        key = tuple(sorted(cnt.items()))
        dpat[key] += 1

print("\n  By denticity composition (top 20):")
for k, n in dpat.most_common(20):
    print(f"    {n:>4}  {k}")

# Theoretical isomer counts in 0% records — was the polyhedron actually enumerated?
theory_cnts = [len(r.get("theoretical_isomers", [])) for r in zero_records]
print(f"\n  Theoretical_isomers count among 0%-coverage records:")
print(f"    n with theoretical=0: {sum(1 for t in theory_cnts if t == 0)}")
print(f"    n with theoretical=1: {sum(1 for t in theory_cnts if t == 1)}")
print(f"    n with theoretical>=2: {sum(1 for t in theory_cnts if t >= 2)}")
print(f"    n with theoretical>=10: {sum(1 for t in theory_cnts if t >= 10)}")
print(f"    n with theoretical>=50: {sum(1 for t in theory_cnts if t >= 50)}")
print(f"    max theoretical: {max(theory_cnts) if theory_cnts else 0}")

# CRITICAL: split by 0% records into:
#   (a) records that have theory >= 2 but observed 0 of them = TRUE coverage failure
#   (b) records that have theory == 0 = enumeration failure / polyhedron not in table
theory0 = [r for r in zero_records if len(r.get("theoretical_isomers", [])) == 0]
theory_ge2 = [r for r in zero_records if len(r.get("theoretical_isomers", [])) >= 2]
print(f"\n  Decomposition of 0%-coverage:")
print(f"    Theory=0 (enumeration broken, no isomers known): {len(theory0)} -- NOT real coverage issue")
print(f"    Theory>=2 + observed=0 (BUILDER missed all): {len(theory_ge2)} -- TRUE coverage failure")

# In theory0 group, why is theory=0? Likely polyhedron_class is not enumerable
print("\n  Top polyhedron_class for theory=0:")
poly0 = Counter(r.get("polyhedron_class", "?") for r in theory0)
for p, n in poly0.most_common(15):
    print(f"    {n:>4}  {p}")

# For theory >=2 + observed 0, what classes/CN/poly?
print("\n  Theory>=2 + observed=0 (TRUE coverage failures):")
print("    By class:")
for cls, n in Counter(r.get("class","?") for r in theory_ge2).most_common():
    print(f"      {cls:>12}: {n}")
print("    By CN:")
for cn, n in sorted(Counter(r.get("cn","?") for r in theory_ge2).items(), key=lambda x: (isinstance(x[0], str), x[0])):
    print(f"      CN{cn}: {n}")
print("    By metal class:")
for mc, n in Counter(metal_class(r.get("metal","?")) for r in theory_ge2).most_common():
    print(f"      {mc:>10}: {n}")
print("    Top polyhedron_class:")
for p, n in Counter(r.get("polyhedron_class","?") for r in theory_ge2).most_common(15):
    print(f"      {n:>4}  {p}")

# Save those true-failure records for deeper trace
with (OUT_DIR / "true_zero_failures.jsonl").open("w") as f:
    for r in theory_ge2:
        f.write(json.dumps(r) + "\n")

# Now also look at PARTIAL coverage records (1-99%) - what isomers are missing systematically?
partial = [r for r in records if 0 < r.get("coverage_pct", 0) < 100]
print(f"\n=== Partial-coverage records (1-99%): {len(partial)} ===")
print(f"  mean coverage: {sum(r.get('coverage_pct',0) for r in partial)/max(1,len(partial)):.1f}%")
# By class
for cls, n in Counter(r.get("class","?") for r in partial).most_common():
    cov = sum(r.get("coverage_pct",0) for r in partial if r.get("class","?")==cls) / max(1,n)
    print(f"    {cls:>12}: n={n}  mean={cov:.1f}%")
