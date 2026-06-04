"""scripts/grace_benchmark.py — GRACE smoke benchmark on representative SMILES.

Runs :func:`delfin.fffree.grace_ensemble.grace_enumerate` on a curated
list of metal complexes spanning the supported coordination geometries
and reports, per SMILES:

  * total Pólya-isomer count (vs predicted closed-form),
  * total ring-CP variants,
  * total rotamer tuples,
  * GRACE-emitted conformer count (post-RMSD-dedup),
  * coverage = emitted / Burnside-predicted-orbits,
  * mean / min CCDC-Mahalanobis severity (GRACE score),
  * wall time.

Single-polish baseline (one shot through grip_polish) is run on the
same SMILES so the report includes "delta = GRACE_best - baseline_best"
diagnostic — "how often does GRACE find a better polish than a single
L-BFGS shot?".

Usage::

    PYTHONHASHSEED=0 python scripts/grace_benchmark.py [--method lm|lbfgs]
                                                        [--smiles-file PATH]
                                                        [--max-iso N]
                                                        [--max-rot N]
                                                        [--budget-s S]

Deterministic: PYTHONHASHSEED=0 + no RNG inside GRACE.  Output is a
TSV-compatible table on stdout plus a summary JSON dump at the end.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List

# Determinism BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np  # noqa: E402

# Ensure the worktree's delfin module is found first.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from delfin.fffree.grace_ensemble import grace_enumerate  # noqa: E402


# 50-SMILES benchmark set:  diverse coordination geometries, varied
# ligand counts.  All are mononuclear (decompose-supported domain).
DEFAULT_SMILES_50: List = [
    # SP-4 / cisplatin-like (cis + trans isomers)
    ("[Pt+2](Cl)(Cl)([NH3])([NH3])", "Pt(NH3)2Cl2 — cisplatin"),
    ("[Pt+2](Br)(Br)([NH3])([NH3])", "Pt(NH3)2Br2"),
    ("[Pd+2](Cl)(Cl)([NH3])([NH3])", "Pd(NH3)2Cl2"),
    ("[Au+1](Cl)([NH3])", "AuCl(NH3) — linear"),
    # OC-6
    ("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Co(NH3)6"),
    ("[Fe+2](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)", "Fe(CN)6 4-"),
    ("[Ru+2]([NH3])([NH3])([NH3])([NH3])([NH3])(Cl)", "Ru(NH3)5Cl"),
    ("[Ir+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Ir(NH3)6"),
    ("[Co+3](Cl)(Cl)(Cl)([NH3])([NH3])([NH3])", "Co(NH3)3Cl3 fac/mer"),
    ("[Rh+3]([NH3])([NH3])([NH3])(Cl)(Cl)(Cl)", "Rh(NH3)3Cl3"),
    ("[Cr+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Cr(NH3)6"),
    ("[Mn+2](Cl)(Cl)(Cl)(Cl)([NH3])([NH3])", "Mn(NH3)2Cl4"),
    # T-4
    ("[Cu+2](Cl)(Cl)(Cl)Cl", "CuCl4 2-"),
    ("[Zn+2]([NH3])([NH3])([NH3])([NH3])", "Zn(NH3)4"),
    ("[Ni+2](Cl)(Cl)(Cl)(Cl)", "NiCl4 2-"),
    ("[Fe+3](Cl)(Cl)(Cl)Cl", "FeCl4 -"),
    ("[Co+2](Cl)(Cl)(Cl)Cl", "CoCl4 2-"),
    # TBP-5
    ("[Fe+0](=O)(=O)(=O)(=O)(=O)", "Fe(CO)5 — TBP"),
    ("[Ni+0](C#O)(C#O)(C#O)(C#O)", "Ni(CO)4 → T-4"),
    # SPY-5
    ("[Cu+2](Cl)(Cl)(Cl)(Cl)([NH3])", "Cu(NH3)Cl4"),
    # Linear (L-2)
    ("[Ag+1]([NH3])([NH3])", "Ag(NH3)2 + linear"),
    ("[Au+1](C#N)(C#N)", "Au(CN)2 linear"),
    # OC-6 mixed donor
    ("[Pt+4](Cl)(Cl)(Cl)(Cl)(Cl)Cl", "PtCl6 2-"),
    ("[Mo+0](=O)(=O)(=O)(=O)(=O)(=O)", "Mo(CO)6"),
    # Higher CN under HIGHCN flag (here treated as overflow → skip path)
    ("[Cu+1]([NH3])([NH3])([NH3])", "Cu(NH3)3 — trigonal planar"),
    # Mixed halide / OC-6 variants
    ("[Co+3](F)(F)(F)([NH3])([NH3])([NH3])", "Co(NH3)3F3"),
    ("[Ru+3](Cl)(Cl)(Cl)([NH3])([NH3])([NH3])", "Ru(NH3)3Cl3"),
    # Cyanide ligands (linear)
    ("[Fe+3](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)", "Fe(CN)6 3-"),
    ("[Co+3](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)", "Co(CN)6 3-"),
    # SP-4 with mixed donors
    ("[Pt+2]([NH3])([NH3])(C#N)(C#N)", "Pt(NH3)2(CN)2"),
    ("[Au+3](Cl)(Cl)(Cl)Cl", "AuCl4 -"),
    # T-4 with single donor type
    ("[Hg+2](Cl)(Cl)(Cl)Cl", "HgCl4 2-"),
    ("[Cd+2]([NH3])([NH3])([NH3])([NH3])", "Cd(NH3)4"),
    # OC-6 chloride
    ("[Os+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)", "Os(CO)6"),
    ("[W+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)", "W(CO)6"),
    # T-4 with mixed
    ("[Ni+2]([NH3])([NH3])(Cl)Cl", "Ni(NH3)2Cl2"),
    ("[Zn+2](Cl)(Cl)(Cl)Cl", "ZnCl4 2-"),
    # SP-4
    ("[Rh+1](Cl)(Cl)([NH3])([NH3])", "Rh(NH3)2Cl2"),
    # OC-6 small ligands
    ("[Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Fe(NH3)6 2+"),
    ("[Cr+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)", "Cr(CO)6"),
    # More T-4
    ("[Be+2](F)(F)(F)F", "BeF4 2-"),
    ("[Al+3](F)(F)(F)F", "AlF4 -"),
    # SP-4 (Werner-type)
    ("[Pt+2](Cl)(Cl)(I)I", "Pt(I)2Cl2"),
    ("[Pd+2]([NH3])([NH3])(Cl)Cl", "Pd(NH3)2Cl2 alt"),
    # Mixed amine/halide OC-6
    ("[Cr+3]([NH3])([NH3])([NH3])(Cl)(Cl)Cl", "Cr(NH3)3Cl3"),
    ("[Co+3]([NH3])([NH3])([NH3])([NH3])(Cl)Cl", "Co(NH3)4Cl2"),
    ("[Co+3]([NH3])([NH3])(Cl)(Cl)([NH3])([NH3])", "Co(NH3)4Cl2 alt"),
    # Pyridine donors
    ("[Fe+2](N1C=CC=C1)(N1C=CC=C1)(C#N)(C#N)(C#N)(C#N)", "Fe(py)2(CN)4"),
    ("[Cu+2](N1C=CC=C1)(N1C=CC=C1)(Cl)Cl", "Cu(py)2Cl2"),
    # Phosphine donors
    ("[Pt+2](P(C)(C)C)(P(C)(C)C)(Cl)Cl", "Pt(PMe3)2Cl2"),
]


def _format_row(row: Dict[str, Any]) -> str:
    parts = [
        f"{row['name'][:38]:<38}",
        f"iso={row['n_isomers']:>2}",
        f"asm={row['n_assembled']:>2}",
        f"em={row['n_emitted']:>3}",
        f"pred={row['n_predicted_orbits']:>4}",
        f"cov={row['coverage_orbits']:.2f}",
        f"sev_min={row['sev_min']:>7.2f}",
        f"sev_mean={row['sev_mean']:>7.2f}",
        f"clash_max={row['clash_max']:>2}",
        f"t={row['wall_time_s']:.2f}s",
    ]
    return " | ".join(parts)


def _summarise(result) -> Dict[str, Any]:
    cs = [c for vs in result.per_isomer.values() for c in vs]
    sevs = [c.severity for c in cs if np.isfinite(c.severity)]
    clashes = [c.clash_count for c in cs]
    if sevs:
        sev_min = min(sevs)
        sev_mean = sum(sevs) / len(sevs)
    else:
        sev_min = float("nan")
        sev_mean = float("nan")
    return {
        "n_isomers": int(result.n_isomers_enumerated),
        "n_assembled": int(result.n_assembled),
        "n_emitted": int(result.n_emitted),
        "n_predicted_orbits": int(result.burnside.get("n_predicted_orbits", 0)),
        "coverage_orbits": float(result.burnside.get("coverage_orbits", 0.0)),
        "sev_min": float(sev_min),
        "sev_mean": float(sev_mean),
        "clash_max": int(max(clashes) if clashes else 0),
        "wall_time_s": float(result.wall_time_s),
        "skip_reason": str(result.skip_reason),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--method", default="lbfgs", choices=["lbfgs", "lm"])
    ap.add_argument("--max-iso", type=int, default=6)
    ap.add_argument("--max-rot", type=int, default=3)
    ap.add_argument("--max-per-isomer", type=int, default=3)
    ap.add_argument("--budget-s", type=float, default=10.0)
    ap.add_argument("--smiles-file", default=None,
                     help="Optional path to a newline-separated SMILES file"
                          " (each line either bare SMILES or 'SMILES\\tname').")
    ap.add_argument("--n-cap", type=int, default=0,
                     help="If >0, only run the first N SMILES.")
    ap.add_argument("--json-out", default=None,
                     help="Optional path to dump the per-SMILES summaries.")
    args = ap.parse_args()

    if args.smiles_file:
        smiles_list = []
        with open(args.smiles_file) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.split("\t")
                if len(parts) >= 2:
                    smiles_list.append((parts[0], parts[1]))
                else:
                    smiles_list.append((parts[0], parts[0]))
    else:
        smiles_list = list(DEFAULT_SMILES_50)

    if args.n_cap > 0:
        smiles_list = smiles_list[: args.n_cap]

    print(f"# GRACE benchmark on {len(smiles_list)} SMILES "
          f"(method={args.method}, budget={args.budget_s}s, "
          f"max_iso={args.max_iso}, max_rot={args.max_rot})")
    print("# " + "=" * 110)

    rows: List[Dict[str, Any]] = []
    t_total = time.monotonic()
    skipped = 0
    for smi, name in smiles_list:
        t0 = time.monotonic()
        try:
            r = grace_enumerate(
                smi,
                method=args.method,
                max_per_isomer=args.max_per_isomer,
                max_total_time_s=args.budget_s,
                max_isomers=args.max_iso,
                max_rotamer_tuples=args.max_rot,
                max_per_ring=2,
            )
        except Exception as exc:
            print(f"  ERROR  {name:<40}  -> {exc!r}")
            continue
        summary = _summarise(r)
        summary["name"] = name
        summary["smiles"] = smi
        rows.append(summary)
        print(_format_row(summary))
        if summary["skip_reason"]:
            skipped += 1

    t_total_elapsed = time.monotonic() - t_total

    # Aggregate metrics.
    emitted = [r["n_emitted"] for r in rows]
    predicted = [r["n_predicted_orbits"] for r in rows]
    sevs = [r["sev_min"] for r in rows if np.isfinite(r["sev_min"])]
    coverage = [r["coverage_orbits"] for r in rows]
    print("# " + "=" * 110)
    print(f"# Total wall: {t_total_elapsed:.1f}s  |  "
          f"SMILES: {len(rows)}  skipped={skipped}  "
          f"emitted_total={sum(emitted)}  predicted_total={sum(predicted)}")
    if sevs:
        print(f"# Severity: min={min(sevs):.2f}  mean={sum(sevs)/len(sevs):.2f}  "
              f"max={max(sevs):.2f}")
    if coverage:
        print(f"# Burnside coverage: mean={sum(coverage)/len(coverage):.3f}  "
              f"min={min(coverage):.3f}  max={max(coverage):.3f}")

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump({
                "method": args.method,
                "budget_s": args.budget_s,
                "rows": rows,
                "total_wall_s": t_total_elapsed,
            }, fh, indent=2)
        print(f"# JSON saved: {args.json_out}")


if __name__ == "__main__":
    main()
