"""Phase E: validate multi-metal assemble on 20 SMILES from master_v3_plus.

Run from worktree root:
    PYTHONHASHSEED=0 DELFIN_FFFREE_MULTI_METAL=1 python iters/MULTI_METAL_PHASE_B_validation.py

Outputs:
    iters/MULTI_METAL_PHASE_B_validation.json
"""
from __future__ import annotations

import json
import os
import re
import sys

from rdkit import Chem


METALS = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
          "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
          "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Sn"}


def parse_smiles(raw_line: str) -> str:
    """master_v3_plus.txt has lines like 'NAME|SMILES'."""
    if "|" in raw_line:
        return raw_line.split("|", 1)[1].strip()
    return raw_line.strip()


def count_metals(smi: str) -> int:
    bracket = re.findall(r"\[([A-Z][a-z]?)", smi)
    return sum(1 for a in bracket if a in METALS)


def main() -> int:
    master = "/home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt"
    if not os.path.exists(master):
        print(f"missing: {master}", file=sys.stderr)
        return 1
    multi: list[str] = []
    with open(master) as fp:
        for line in fp:
            smi = parse_smiles(line)
            if not smi:
                continue
            if count_metals(smi) >= 2:
                multi.append(smi)
            if len(multi) >= 60:
                break
    print(f"Found {len(multi)} multi-metal candidate SMILES; testing first 20")
    multi = multi[:20]

    os.environ["DELFIN_FFFREE_MULTI_METAL"] = "1"
    os.environ.setdefault("PYTHONHASHSEED", "0")

    from delfin.fffree import multi_metal_assemble as MMA

    rows: list[dict] = []
    for i, smi in enumerate(multi):
        rec: dict = {"i": i, "smiles": smi[:120]}
        try:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
        except Exception as exc:
            rec["error"] = f"parse:{exc}"
            rows.append(rec)
            continue
        if mol is None:
            rec["error"] = "rdkit_returned_none"
            rows.append(rec)
            continue
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                              | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                              | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception:
            pass
        try:
            g = MMA.build_metal_graph(mol)
        except Exception as exc:
            rec["error"] = f"graph:{exc}"
            rows.append(rec)
            continue
        rec["n_metals"] = len(g["metal_idxs"])
        rec["metal_syms"] = list(g["metal_syms"])
        rec["mm_declared"] = len(g["mm_edges_declared"])
        rec["mm_implied"] = len(g["mm_edges_implied"])
        rec["n_bridges"] = len(g["bridges"])
        try:
            sk = MMA.choose_skeleton(g)
        except Exception as exc:
            rec["error"] = f"skeleton:{exc}"
            rows.append(rec)
            continue
        rec["skeleton_kind"] = None if sk is None else sk["kind"]
        if sk is not None:
            rec["mm_distance_target"] = round(float(sk["mm_distance"]), 3)

        try:
            res = MMA.assemble_multi_metal(mol)
        except Exception as exc:
            rec["error"] = f"assemble:{exc}"
            rows.append(rec)
            continue
        if res is None:
            rec["built"] = False
            rec["status"] = "rollback_or_failed"
        else:
            syms, P, donors = res
            rec["built"] = True
            rec["n_atoms_out"] = int(P.shape[0])
            rec["n_donors"] = len(donors)
            rec["status"] = "ok"
            # M-M distance recovery: compute observed M-M for each metal pair
            mm_obs = []
            n_m = rec["n_metals"]
            for ii in range(n_m):
                for jj in range(ii + 1, n_m):
                    import numpy as _np
                    d = float(_np.linalg.norm(P[ii] - P[jj]))
                    mm_obs.append(round(d, 3))
            rec["mm_obs_all"] = mm_obs
        rows.append(rec)

    out_path = os.path.join(os.path.dirname(__file__),
                            "MULTI_METAL_PHASE_B_validation.json")
    with open(out_path, "w") as fp:
        json.dump({
            "n_total": len(rows),
            "n_built": sum(1 for r in rows if r.get("built")),
            "n_skeleton_chosen": sum(1 for r in rows if r.get("skeleton_kind")),
            "rows": rows,
        }, fp, indent=2, sort_keys=True)
    print(f"Wrote {out_path}")
    print(f"  total       = {len(rows)}")
    print(f"  skeleton OK = {sum(1 for r in rows if r.get('skeleton_kind'))}")
    print(f"  built OK    = {sum(1 for r in rows if r.get('built'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
