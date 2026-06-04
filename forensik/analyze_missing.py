#!/usr/bin/env python3
"""Phase C/D: look at missing-isomer label distributions to find systematic gaps."""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

JSONL = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL_isocoverage.jsonl")
OUT_DIR = Path("/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a6a0c001c3fc98a82/forensik")

with JSONL.open() as f:
    records = [json.loads(line) for line in f if line.strip()]

# What kind of isomer labels appear most often as MISSING?
# Group by class and label-prefix to find systematic patterns.
print("=== Missing-isomer label analysis ===")

# Look at theoretical_geoms across all records
geom_dist_overall = Counter()
geom_dist_zero = Counter()
geom_dist_low = Counter()  # 0 < cov < 25
geom_dist_100 = Counter()
for r in records:
    geoms = tuple(sorted(r.get("theoretical_geoms", [])))
    cov = r.get("coverage_pct", 0)
    geom_dist_overall[geoms] += 1
    if cov == 0:
        geom_dist_zero[geoms] += 1
    elif cov <= 25:
        geom_dist_low[geoms] += 1
    elif cov == 100:
        geom_dist_100[geoms] += 1

print("\n  Theoretical geom-sets by coverage bucket:")
print(f"  {'geom-set':<30}{'total':>8}{'0%':>8}{'1-25%':>8}{'100%':>8}{'0%-rate':>10}")
for geoms, n in geom_dist_overall.most_common(25):
    n0 = geom_dist_zero.get(geoms, 0)
    nl = geom_dist_low.get(geoms, 0)
    n100 = geom_dist_100.get(geoms, 0)
    if n0 + nl < 5:
        continue
    rate = 100.0 * n0 / n
    print(f"  {str(geoms)[:30]:<30}{n:>8}{n0:>8}{nl:>8}{n100:>8}{rate:>9.1f}%")

# For 0% records, look at which isomer labels are most commonly missing
print("\n\n=== Most-missing isomer labels (in 0%-records) ===")
all_missing = Counter()
for r in records:
    if r.get("coverage_pct", 0) == 0 and len(r.get("theoretical_isomers", [])) >= 2:
        for lab in r.get("theoretical_isomers", []):
            all_missing[lab] += 1
print("  Top 30 missing labels:")
for lab, n in all_missing.most_common(30):
    print(f"    {n:>4}  {lab}")

# What's the structure of the missing labels?
# Pattern analysis: are these axial-cap (PBP), all-cis (TPR), etc?
prefix_dist = Counter()
for r in records:
    if r.get("coverage_pct", 0) == 0 and len(r.get("theoretical_isomers", [])) >= 2:
        for lab in r.get("theoretical_isomers", []):
            # extract prefix word
            m = re.match(r"([a-z]+|[A-Z][^/-]*)", lab)
            if m:
                prefix_dist[m.group(1)] += 1
print("\n  Top missing-label prefixes:")
for p, n in prefix_dist.most_common(20):
    print(f"    {n:>5}  {p!r}")

# Pattern-based: are these strings "cis/trans"-style 2-isomer simple?
two_iso_records = [r for r in records if len(r.get("theoretical_isomers", [])) == 2
                   and r.get("coverage_pct", 0) == 0]
print(f"\n\n=== 0%-coverage records with theory=2 (simple cis/trans only): {len(two_iso_records)} ===")
twopair = Counter()
for r in two_iso_records:
    ti = tuple(sorted(r.get("theoretical_isomers", [])))
    twopair[ti] += 1
print("  Top 15 pairs:")
for pair, n in twopair.most_common(15):
    print(f"    {n:>4}  {pair}")

# Detailed look at the bis-bidentate-NN+NN failure
# These should be M(en)2X2 type complexes with 3 expected isomers
print("\n\n=== bis-bidentate-NN+NN (Werner cis/trans) full detail ===")
recN = [r for r in records if r.get("polyhedron_class","") == "bis-bidentate-NN+NN" and r.get("coverage_pct",0)==0]
print(f"  Total 0%-coverage: {len(recN)}")
# What are their CNs?
print(f"  CNs: {Counter(r['cn'] for r in recN)}")
# Class?
print(f"  Class: {Counter(r['class'] for r in recN)}")
# Show 3 samples
for r in recN[:3]:
    print(f"\n    smiles_id: {r.get('smiles_id','?')}")
    print(f"    smiles: {r['smiles'][:130]}")
    print(f"    cn={r.get('cn','?')}  theory_geoms={r.get('theoretical_geoms',[])}")
    print(f"    theoretical_isomers: {r.get('theoretical_isomers',[])}")
    print(f"    observed: {[o.get('best_geom','?') for o in r.get('observed_isomers',[])]}")
    print(f"    obs labels: {[o.get('geometric_label','?') for o in r.get('observed_isomers',[])]}")
    print(f"    obs matched: {[o.get('matched_name','?') for o in r.get('observed_isomers',[])]}")

# Now look at partial-coverage records to see what is systematically missing
print("\n\n=== PARTIAL-coverage records: systematically-missing labels ===")
# Bin them by (class, polyhedron_class) and look at the top-missing labels
partial = [r for r in records if 0 < r.get("coverage_pct",0) < 100]
miss_per_class = defaultdict(Counter)
for r in partial:
    cls = r.get("class","?")
    for lab in r.get("named_missing", []):
        miss_per_class[cls][lab] += 1

for cls in ["sigma", "hapto", "multi-sigma", "multi-hapto"]:
    cc = miss_per_class.get(cls, Counter())
    if not cc:
        continue
    print(f"\n  {cls}: top 15 missing labels in partial records:")
    for lab, n in cc.most_common(15):
        print(f"    {n:>5}  {lab}")

# Look at "ax" vs "cap" splits for PBP/COH (CN7)
print("\n\n=== CN7 (PBP/COH) detail: ax vs cap label distribution ===")
cn7 = [r for r in records if r.get("cn",0) == 7]
print(f"  Total CN7 records: {len(cn7)}")
ax_obs = 0; cap_obs = 0; ax_miss = 0; cap_miss = 0
for r in cn7:
    for lab in r.get("named_observed", []):
        if "-ax" in lab: ax_obs += 1
        if lab.startswith("cap"): cap_obs += 1
    for lab in r.get("named_missing", []):
        if "-ax" in lab: ax_miss += 1
        if lab.startswith("cap"): cap_miss += 1
print(f"  ax labels observed: {ax_obs}  missing: {ax_miss}  -> miss-rate {100*ax_miss/max(1,ax_obs+ax_miss):.1f}%")
print(f"  cap labels observed: {cap_obs}  missing: {cap_miss}  -> miss-rate {100*cap_miss/max(1,cap_obs+cap_miss):.1f}%")
