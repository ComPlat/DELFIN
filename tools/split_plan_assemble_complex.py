"""Plan for delfin/manta/assemble_complex.py (line ranges of the file at 7a352904).

Usage: python tools/split_plan_assemble_complex.py <depgraph.json> <out_dir>

Writes <out_dir>/names/<module>.txt (one leading name per line) and
<out_dir>/plan.json, the step list split_extract.py consumes.
"""
import json
import sys
from pathlib import Path
S = Path(sys.argv[2])  # output directory
S.mkdir(parents=True, exist_ok=True)
g = json.load(open(sys.argv[1]))
by_line = [(t["l0"], t["names"][0]) for t in g["tops"] if t["names"] and t["kind"] not in ("Import", "ImportFrom")]
GROUPS = [
    ("assemble_donor_plane", "Donor-plane relaxation of the FF-free assembly: collapsed-bond checks, donor follow weights, arm order, riding hydrogens, beta scoring and the trilateration rescue switch.", [(647, 1177)], ["_bite_aware_targets", "_trilaterate_donor_targets"]),
    ("assemble_ligand_embed", "Ligand embedding for the FF-free assembly: rigid alignment, Kabsch, ring-bound tightening, metallacycle embedding, planar polydentate placement and rigid cavity conformers.", [(24, 646)]),
    ("assemble_orient", "Orientation of chelates onto polyhedron vertices, donor bend angles, sp-chain straightening, VSEPR reconstruction and diatomic donor orientation of the FF-free assembly.", [(1178, 1937)]),
    ("assemble_chelate", "Monodentate and chelate assembly of the FF-free constructor: clash relief, bite-aware and trilaterated donor targets, lone-pair orientation, bite closure and multichelate assembly.", [(1938, 2475)]),
    ("assemble_ligand_confs", "Ligand conformers of the FF-free constructor: cached MMFF-free relaxation, degenerate symmetrisation, clash count, torsion and joint declash frames, guarded refinement and sphere flex.", [(2476, 2926)]),
    ("assemble_fold_fp", "Fold fingerprints of the FF-free constructor: ring folds from blocks, the amplitude axis, metallacycle arms and complex RMSD for the dedup.", [(3117, 3646)]),
    ("assemble_hapto", "Hapto assembly of the FF-free constructor: eta ring placement, piano-stool leg tilt, ring spins, slip modes, puckers and the hapto ensemble.", [(4003, 4814)]),
    ("assemble_seat", "Seating passes of the FF-free constructor: global donor seat, OC-6 twist seat and ligand-DOF reseat.", [(4815, 5447)]),
    ("assemble_ensemble", "Heteroleptic assembly from mols, the heteroleptic ensemble, the constrained UFF relax and build_and_relax of the FF-free constructor.", [(2927, 3116), (3647, 4002)]),
]
plan = []
(S / "names").mkdir(exist_ok=True)
EXTRA_HOME = {"_bite_aware_targets": "assemble_donor_plane", "_trilaterate_donor_targets": "assemble_donor_plane"}
for mod, doc, ranges, *extra in GROUPS:
    names = [n for (l, n) in by_line if any(lo <= l <= hi for lo, hi in ranges) and EXTRA_HOME.get(n, mod) == mod]
    names += [n for n in (extra[0] if extra else []) if n not in names]
    (S / "names" / f"{mod}.txt").write_text("\n".join(names) + "\n")
    plan.append({"module": f"delfin/manta/{mod}.py", "doc": doc, "names_file": str(S / "names" / f"{mod}.txt")})
json.dump(plan, open(S / "plan.json", "w"), indent=1)
print("groups:", len(plan))
