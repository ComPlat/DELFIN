#!/usr/bin/env python3
"""Phase C: detailed examination of representative failure SMILES."""
import json
from collections import Counter, defaultdict
from pathlib import Path

OUT_DIR = Path("/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a6a0c001c3fc98a82/forensik")

with (OUT_DIR / "true_zero_failures.jsonl").open() as f:
    records = [json.loads(line) for line in f if line.strip()]

print(f"Records: {len(records)}")

# Group by CN x polyhedron_class
group = defaultdict(list)
for r in records:
    key = (r.get("cn","?"), r.get("polyhedron_class","?"))
    group[key].append(r)

# Show top groups + 1 representative SMILES each
print("\n=== Top failure groups (CN, polyhedron_class) with representative SMILES + theoretical_isomers count ===")
top = sorted(group.items(), key=lambda x: -len(x[1]))[:30]
for (cn, pc), recs in top:
    rep = recs[0]
    ti = rep.get("theoretical_isomers", [])
    obs = rep.get("named_observed", [])
    print(f"\n  [{len(recs):>3}] CN={cn}  poly={pc}")
    print(f"      SMILES: {rep['smiles'][:120]}")
    print(f"      smiles_id: {rep.get('smiles_id','?')}")
    print(f"      theoretical_geoms: {rep.get('theoretical_geoms',[])}")
    print(f"      theoretical_isomers ({len(ti)}): {ti[:6]}...")
    print(f"      observed: {len(obs)}; observed_named: {obs[:3]}")
    print(f"      named_missing_cn_mismatch: {len(rep.get('named_missing_cn_mismatch',[]))}")

# Also look at low-cn failures specifically
print("\n\n=== Detail: 'low-cn' failures (CN<4) ===")
lowcn = [r for r in records if r.get("polyhedron_class") == "low-cn"]
print(f"  Total: {len(lowcn)}")
cn_dist = Counter(r.get("cn","?") for r in lowcn)
print(f"  CN distribution: {dict(cn_dist)}")
for r in lowcn[:5]:
    ti = r.get("theoretical_isomers",[])
    print(f"    SMILES: {r['smiles'][:120]}")
    print(f"      cn={r.get('cn','?')}  class={r.get('class','?')}  metal={r.get('metal','?')}")
    print(f"      theoretical_isomers ({len(ti)}): {ti[:8]}")
    print(f"      observed: {len(r.get('named_observed',[]))} named")
    print()

# bis-monodentate-Cl+Cl etc
print("\n=== Detail: simple 2-monodentate-symmetric pairs (e.g. Cl+Cl) ===")
simple_pair = [r for r in records if "bis-monodentate" in r.get("polyhedron_class","") and "+" in r.get("polyhedron_class","")]
print(f"  Total bis-monodentate-X+Y: {len(simple_pair)}")
sp_pc = Counter(r.get("polyhedron_class","?") for r in simple_pair)
print(f"  Patterns: {dict(sp_pc.most_common(10))}")
# Pick one
for r in simple_pair[:5]:
    ti = r.get("theoretical_isomers",[])
    print(f"    SMILES: {r['smiles'][:120]}")
    print(f"      cn={r.get('cn','?')}  poly={r.get('polyhedron_class','?')}")
    print(f"      theoretical_isomers ({len(ti)}): {ti[:6]}")
    print(f"      observed: {len(r.get('named_observed',[]))} named")
    print()

# bis-bidentate-NN+NN
print("\n=== Detail: bis-bidentate-NN+NN failures ===")
bisN = [r for r in records if r.get("polyhedron_class","") == "bis-bidentate-NN+NN"]
print(f"  Total: {len(bisN)}")
for r in bisN[:5]:
    ti = r.get("theoretical_isomers",[])
    obs = r.get("named_observed",[])
    print(f"    SMILES: {r['smiles'][:120]}")
    print(f"      cn={r.get('cn','?')}  metal={r.get('metal','?')}")
    print(f"      theoretical_isomers ({len(ti)}): {ti}")
    print(f"      observed: {obs}")
    print()

# tris-bidentate-NN+NN+NN
print("\n=== Detail: tris-bidentate-NN+NN+NN failures ===")
trisN = [r for r in records if r.get("polyhedron_class","") == "tris-bidentate-NN+NN+NN"]
print(f"  Total: {len(trisN)}")
for r in trisN[:5]:
    ti = r.get("theoretical_isomers",[])
    obs = r.get("named_observed",[])
    print(f"    SMILES: {r['smiles'][:120]}")
    print(f"      cn={r.get('cn','?')}  metal={r.get('metal','?')}")
    print(f"      theoretical_isomers ({len(ti)}): {ti}")
    print(f"      observed: {obs}")
    print()

# tetradentate-NNNN
print("\n=== Detail: tetradentate-NNNN failures ===")
tt = [r for r in records if r.get("polyhedron_class","") == "tetradentate-NNNN"]
print(f"  Total: {len(tt)}")
for r in tt[:5]:
    ti = r.get("theoretical_isomers",[])
    obs = r.get("named_observed",[])
    print(f"    SMILES: {r['smiles'][:120]}")
    print(f"      cn={r.get('cn','?')}  metal={r.get('metal','?')}")
    print(f"      theoretical_isomers ({len(ti)}): {ti}")
    print(f"      observed: {obs}")
    print()
