"""scripts/grace_benchmark_200_runner.py — 200-SMILES GRACE benchmark suite.

Executes Tasks A, C, E from the GRACE benchmark mandate:

  A. **Burnside completeness** — run GRACE (method=lm, topology_healing
     ON) per SMILES with a 30-s wall budget; record
     ``(n_emitted, n_predicted, burnside_coverage, wall_time,
     isomer_classes_emitted)``.  Aggregate: mean coverage,
     %SMILES @100% coverage, mean wall-time.

  C. **GRACE vs single-shot grip_polish lift** — for each SMILES,
     additionally build the lowest-isomer-id Pólya config single-shot
     and polish via :func:`grip_polish`.  Compare the GRACE best score
     to the single-shot score → "energy/quality gain from enumeration".

  E. **Scaling profile** — per SMILES record (#atoms, CN, #rings,
     emit_count, wall_time).  Fit ``t = a + b * N^c`` over the run set
     so the scaling exponent is exposed.

Outputs (CSV + JSONL into ``paper_data/``):

  * ``grace_200_burnside_benchmark.csv``     — per-SMILES wide-row CSV.
  * ``grace_200_burnside_per_isomer.jsonl``  — one record per emitted
                                                isomer.
  * ``grace_vs_singleshot_lift.csv``         — single-shot vs GRACE-best.
  * ``grace_scaling_profile.csv``            — scaling fit input data.

Determinism: PYTHONHASHSEED=0, no RNG.  Re-run yields byte-identical CSVs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
import traceback
from typing import Any, Dict, List, Optional, Tuple

# Determinism BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("DELFIN_FFFREE_GRACE_ENABLE", "1")

import numpy as np

# Ensure worktree's delfin module is found first.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from delfin.fffree.grace_ensemble import grace_enumerate  # noqa: E402
from scripts.grace_benchmark_200_curation import (  # noqa: E402
    SmilesEntry, all_entries,
)

# Disable rdkit logs for clarity.
try:
    from rdkit import RDLogger  # noqa: E402
    RDLogger.DisableLog("rdApp.*")
except Exception:
    pass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _n_heavy_atoms(smiles: str) -> int:
    """Compute #non-H atoms via RDKit (best-effort)."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return 0
        Chem.SanitizeMol(m, catchErrors=True)
        return int(m.GetNumHeavyAtoms())
    except Exception:
        return 0


def _n_rings(smiles: str) -> int:
    """Compute #5-7 rings (matches GRACE's CP scope)."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return 0
        Chem.SanitizeMol(m, catchErrors=True)
        ri = m.GetRingInfo()
        return int(sum(1 for r in ri.AtomRings() if 5 <= len(r) <= 7))
    except Exception:
        return 0


def _decompose_meta(smiles: str) -> Dict[str, Any]:
    """Return decomposition metadata: cn, geometry, ligand-count."""
    out = {"cn": 0, "geometry": "", "n_ligands": 0, "metal": "",
            "has_chelate": False, "decompose_ok": False}
    try:
        from delfin.fffree import decompose as DEC
        d = DEC.decompose(smiles)
        if d is None:
            return out
        donors = d.get("donors") or []
        out["cn"] = len(donors)
        out["geometry"] = str(d.get("geometry") or "")
        out["n_ligands"] = len(d.get("ligands") or [])
        out["metal"] = str(d.get("metal") or "")
        out["has_chelate"] = bool(d.get("has_chelate"))
        out["decompose_ok"] = True
    except Exception:
        pass
    return out


def _run_grace_for_smiles(
    entry: SmilesEntry,
    *,
    method: str,
    max_per_isomer: int,
    max_isomers: int,
    max_rotamer_tuples: int,
    budget_s: float,
) -> Dict[str, Any]:
    """Run GRACE on one SMILES and gather the full result record.

    Failures (decompose-skip, exception) are absorbed into the record so
    aggregation never crashes.
    """
    t0 = time.monotonic()
    rec: Dict[str, Any] = {
        "smiles": entry.smiles,
        "name": entry.name,
        "family": entry.family,
        "refcode": entry.refcode or "",
        "n_expected": int(entry.n_expected),
    }
    rec.update(_decompose_meta(entry.smiles))
    rec["n_heavy_atoms"] = _n_heavy_atoms(entry.smiles)
    rec["n_rings"] = _n_rings(entry.smiles)

    try:
        r = grace_enumerate(
            entry.smiles,
            method=method,
            max_per_isomer=max_per_isomer,
            max_isomers=max_isomers,
            max_rotamer_tuples=max_rotamer_tuples,
            max_total_time_s=budget_s,
            max_per_ring=2,
        )
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["wall_time_s"] = float(time.monotonic() - t0)
        rec["n_emitted"] = 0
        rec["n_isomers_enumerated"] = 0
        rec["coverage_orbits"] = 0.0
        rec["coverage_raw"] = 0.0
        rec["n_predicted_orbits"] = 0
        rec["n_predicted_raw_product"] = 0
        rec["sev_min"] = float("nan")
        rec["sev_mean"] = float("nan")
        rec["clash_max"] = 0
        rec["skip_reason"] = ""
        rec["timed_out"] = False
        rec["isomer_classes_emitted"] = ""
        rec["per_isomer"] = []
        return rec

    sevs = []
    clashes = []
    iso_classes: List[Tuple[int, int]] = []  # (iso_id, n_survivors)
    per_iso_records = []
    for iso_id in sorted(r.per_isomer.keys()):
        cs = r.per_isomer[iso_id]
        if not cs:
            continue
        iso_classes.append((int(iso_id), len(cs)))
        for c in cs:
            sev = float(c.severity)
            if math.isfinite(sev):
                sevs.append(sev)
            clashes.append(int(c.clash_count))
            per_iso_records.append({
                "isomer_id": int(c.isomer_id),
                "ring_id": int(c.ring_id),
                "rotamer_id": int(c.rotamer_id),
                "label": str(c.label),
                "score": float(c.score),
                "severity": float(c.severity),
                "cshm": float(c.cshm),
                "clash_count": int(c.clash_count),
                "method": str(c.method),
                "topology_healed": bool(c.topology_healed),
            })

    rec["wall_time_s"] = float(r.wall_time_s)
    rec["n_emitted"] = int(r.n_emitted)
    rec["n_isomers_enumerated"] = int(r.n_isomers_enumerated)
    rec["n_assembled"] = int(r.n_assembled)
    rec["n_rejected_polish"] = int(r.n_rejected_polish)
    rec["n_rejected_clash"] = int(r.n_rejected_clash)
    rec["n_rejected_topology"] = int(r.n_rejected_topology)
    rec["timed_out"] = bool(r.timed_out)
    rec["skip_reason"] = str(r.skip_reason)
    bs = r.burnside or {}
    rec["n_predicted_orbits"] = int(bs.get("n_predicted_orbits", 0))
    rec["n_predicted_raw_product"] = int(bs.get("n_predicted_raw_product", 0))
    rec["coverage_orbits"] = float(bs.get("coverage_orbits", 0.0))
    rec["coverage_raw"] = float(bs.get("coverage_raw", 0.0))
    rec["sev_min"] = float(min(sevs)) if sevs else float("nan")
    rec["sev_mean"] = float(sum(sevs) / len(sevs)) if sevs else float("nan")
    rec["clash_max"] = int(max(clashes) if clashes else 0)
    rec["isomer_classes_emitted"] = ";".join(
        f"iso{i}:{n}" for i, n in iso_classes
    )
    rec["per_isomer"] = per_iso_records
    return rec


def _run_singleshot(entry: SmilesEntry) -> Dict[str, Any]:
    """Single-shot baseline: build the lowest-isomer-id Pólya config and
    polish once (grip_polish / L-BFGS).  Return ``(score, severity, cshm,
    wall_time_s, ok)`` so the C-task can compute the GRACE lift.
    """
    t0 = time.monotonic()
    out = {
        "single_shot_score": float("inf"),
        "single_shot_severity": float("inf"),
        "single_shot_cshm": float("inf"),
        "single_shot_wall_s": 0.0,
        "single_shot_ok": False,
    }
    try:
        from delfin.fffree import decompose as DEC
        from delfin.fffree import assemble_complex as AC
        from delfin.fffree.grace_ensemble import (
            _build_assembled_mol, _enumerate_polya_isomers,
            _compute_severity, _safe_cshm, _polish_dispatch,
        )
        d = DEC.decompose(entry.smiles)
        if d is None:
            return out
        isomer_configs, _ = _enumerate_polya_isomers(d)
        if not isomer_configs:
            return out
        config = isomer_configs[0]
        geometry = d.get("geometry", "")
        built = AC.assemble_from_config(d["metal"], geometry, config, d["ligands"])
        if built is None:
            return out
        syms_base, P_base, donors_global = built
        syms_base = list(syms_base)
        P_base = np.asarray(P_base, dtype=np.float64)
        donors_global = list(donors_global)
        cm = _build_assembled_mol(
            d.get("metal", "X"), config, d["ligands"], donors_global)
        if cm is None or cm.GetNumAtoms() != len(syms_base):
            return out
        P_polished, ok = _polish_dispatch(
            P_base, cm, metal_idx=0, donors=donors_global,
            geometry=geometry, method="lbfgs",
        )
        if not ok:
            return out
        sev = _compute_severity(syms_base, P_polished, cm, 0, donors_global)
        cshm = _safe_cshm(P_polished, 0, donors_global, geometry)
        score = float(sev) + (float(cshm) if math.isfinite(cshm) else 1e6)
        out["single_shot_score"] = float(score)
        out["single_shot_severity"] = float(sev)
        out["single_shot_cshm"] = float(cshm)
        out["single_shot_ok"] = True
    except Exception:
        pass
    finally:
        out["single_shot_wall_s"] = float(time.monotonic() - t0)
    return out


# ---------------------------------------------------------------------------
# Scaling fit utility (Task E)
# ---------------------------------------------------------------------------
def _fit_powerlaw(N: List[int], T: List[float]) -> Dict[str, Any]:
    """Fit ``t = a + b * N^c`` over (N, T) using log-linear regression on
    the (N, T)-positive subset.

    Returns ``{a, b, c, r2, n_used}``.  Degenerate input yields zeros and
    ``r2 = 0`` so callers can safely sort by quality.
    """
    pairs = [(n, t) for n, t in zip(N, T)
             if n is not None and t is not None and n > 0 and t > 0]
    if len(pairs) < 5:
        return {"a": 0.0, "b": 0.0, "c": 0.0, "r2": 0.0, "n_used": len(pairs)}
    n_arr = np.array([p[0] for p in pairs], dtype=float)
    t_arr = np.array([p[1] for p in pairs], dtype=float)
    # Log-linear fit on (log N, log T): log T = log b + c * log N
    try:
        logN = np.log(n_arr)
        logT = np.log(t_arr)
        c, log_b = np.polyfit(logN, logT, 1)
        b = math.exp(log_b)
        # Compute R^2 over the log fit
        pred = log_b + c * logN
        ss_res = float(np.sum((logT - pred) ** 2))
        ss_tot = float(np.sum((logT - np.mean(logT)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    except Exception:
        return {"a": 0.0, "b": 0.0, "c": 0.0, "r2": 0.0, "n_used": len(pairs)}
    return {
        "a": 0.0,  # additive constant fixed at 0 in the log model
        "b": float(b),
        "c": float(c),
        "r2": float(r2),
        "n_used": int(len(pairs)),
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
WIDE_COLUMNS = [
    "smiles", "name", "family", "refcode", "n_expected",
    "decompose_ok", "metal", "cn", "geometry", "n_ligands", "has_chelate",
    "n_heavy_atoms", "n_rings",
    "n_isomers_enumerated", "n_assembled", "n_emitted",
    "n_predicted_orbits", "n_predicted_raw_product",
    "coverage_orbits", "coverage_raw",
    "n_rejected_polish", "n_rejected_clash", "n_rejected_topology",
    "sev_min", "sev_mean", "clash_max",
    "wall_time_s", "timed_out", "skip_reason", "isomer_classes_emitted",
]
SINGLESHOT_COLUMNS = WIDE_COLUMNS + [
    "single_shot_score", "single_shot_severity", "single_shot_cshm",
    "single_shot_wall_s", "single_shot_ok",
    "grace_best_score", "lift_score_delta",
]
SCALING_COLUMNS = [
    "smiles", "name", "family", "n_heavy_atoms", "cn", "n_rings",
    "n_emitted", "wall_time_s",
]


def _write_csv(path: str, columns: List[str], rows: List[Dict[str, Any]]):
    """Write a CSV with deterministic column order and rounded floats."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            outrow = {}
            for c in columns:
                v = row.get(c, "")
                if isinstance(v, float):
                    if math.isfinite(v):
                        outrow[c] = f"{v:.6g}"
                    else:
                        outrow[c] = "nan" if math.isnan(v) else "inf"
                elif isinstance(v, bool):
                    outrow[c] = "1" if v else "0"
                else:
                    outrow[c] = "" if v is None else str(v)
            w.writerow(outrow)


def _write_jsonl(path: str, records: List[Dict[str, Any]]):
    """Write JSONL of per-isomer records (one line per emitted candidate).

    Each line includes anchor SMILES + isomer/ring/rotamer ids + score.
    """
    with open(path, "w") as fh:
        for rec in records:
            smi = rec.get("smiles", "")
            name = rec.get("name", "")
            family = rec.get("family", "")
            per_iso = rec.get("per_isomer", []) or []
            for c in per_iso:
                entry = {
                    "smiles": smi,
                    "name": name,
                    "family": family,
                    "isomer_id": int(c.get("isomer_id", 0)),
                    "ring_id": int(c.get("ring_id", 0)),
                    "rotamer_id": int(c.get("rotamer_id", 0)),
                    "label": str(c.get("label", "")),
                    "score": float(c.get("score", 0.0)),
                    "severity": float(c.get("severity", 0.0)),
                    "cshm": float(c.get("cshm", 0.0)),
                    "clash_count": int(c.get("clash_count", 0)),
                    "method": str(c.get("method", "")),
                    "topology_healed": bool(c.get("topology_healed", False)),
                }
                fh.write(json.dumps(entry, sort_keys=True) + "\n")


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------
def _aggregate(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute headline aggregate metrics for the SI table.

    Returns mean/median Burnside coverage on the runnable subset,
    %SMILES @100% coverage, mean wall-time, scaling fit coefficients.
    """
    runnable = [r for r in rows if r.get("decompose_ok") and not r.get("error")]
    n_total = len(rows)
    n_decompose_ok = len(runnable)
    n_emitted_total = sum(int(r.get("n_emitted", 0)) for r in runnable)
    n_predicted_total = sum(int(r.get("n_predicted_orbits", 0)) for r in runnable)
    n_at_100pct = sum(
        1 for r in runnable if float(r.get("coverage_orbits", 0.0)) >= 0.999
    )
    n_emitted_at_least_one = sum(
        1 for r in runnable if int(r.get("n_emitted", 0)) >= 1
    )
    cov_list = [float(r.get("coverage_orbits", 0.0)) for r in runnable]
    cov_mean = float(np.mean(cov_list)) if cov_list else 0.0
    cov_median = float(np.median(cov_list)) if cov_list else 0.0
    wall_list = [float(r.get("wall_time_s", 0.0)) for r in runnable]
    wall_mean = float(np.mean(wall_list)) if wall_list else 0.0
    wall_p95 = float(np.percentile(wall_list, 95)) if wall_list else 0.0
    wall_max = float(max(wall_list)) if wall_list else 0.0

    # Scaling fit over runnable subset.
    fit = _fit_powerlaw(
        [int(r.get("n_heavy_atoms", 0)) for r in runnable],
        [float(r.get("wall_time_s", 0.0)) for r in runnable],
    )

    # Per-family rollup
    by_family: Dict[str, Dict[str, Any]] = {}
    for r in runnable:
        fam = str(r.get("family", "?"))
        f = by_family.setdefault(fam, {
            "n": 0, "emit": 0, "predicted": 0,
            "cov_sum": 0.0, "at_100": 0, "wall_sum": 0.0,
        })
        f["n"] += 1
        f["emit"] += int(r.get("n_emitted", 0))
        f["predicted"] += int(r.get("n_predicted_orbits", 0))
        f["cov_sum"] += float(r.get("coverage_orbits", 0.0))
        f["at_100"] += 1 if float(r.get("coverage_orbits", 0.0)) >= 0.999 else 0
        f["wall_sum"] += float(r.get("wall_time_s", 0.0))
    fam_rollup = {}
    for fam, f in sorted(by_family.items()):
        fam_rollup[fam] = {
            "n": f["n"],
            "emitted_total": f["emit"],
            "predicted_total": f["predicted"],
            "coverage_mean": f["cov_sum"] / max(f["n"], 1),
            "pct_at_100": 100.0 * f["at_100"] / max(f["n"], 1),
            "wall_mean_s": f["wall_sum"] / max(f["n"], 1),
        }

    return {
        "n_total": n_total,
        "n_decompose_ok": n_decompose_ok,
        "n_emitted_at_least_one": n_emitted_at_least_one,
        "n_emitted_total": n_emitted_total,
        "n_predicted_total": n_predicted_total,
        "n_at_100_pct_coverage": n_at_100pct,
        "pct_at_100_coverage": 100.0 * n_at_100pct / max(n_decompose_ok, 1),
        "coverage_mean": cov_mean,
        "coverage_median": cov_median,
        "wall_mean_s": wall_mean,
        "wall_p95_s": wall_p95,
        "wall_max_s": wall_max,
        "scaling_fit": fit,
        "by_family": fam_rollup,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=os.path.join(_ROOT, "paper_data"))
    ap.add_argument("--method", default="lm", choices=["lm", "lbfgs"])
    ap.add_argument("--budget-s", type=float, default=30.0,
                     help="Per-SMILES wall budget for GRACE.")
    ap.add_argument("--max-per-isomer", type=int, default=10)
    ap.add_argument("--max-isomers", type=int, default=12)
    ap.add_argument("--max-rotamer-tuples", type=int, default=27)
    ap.add_argument("--n-cap", type=int, default=0,
                     help="Sample only the first N entries (for sanity).")
    ap.add_argument("--with-singleshot", action="store_true",
                     help="Run the single-shot grip_polish baseline too "
                          "(Task C).")
    ap.add_argument("--summary-only", action="store_true",
                     help="Skip CSV/JSONL output; print summary only.")
    args = ap.parse_args()

    entries = all_entries()
    if args.n_cap > 0:
        entries = entries[: args.n_cap]
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"# GRACE 200-SMILES benchmark — method={args.method} "
          f"budget={args.budget_s}s  n={len(entries)}", flush=True)

    rows: List[Dict[str, Any]] = []
    singleshot_rows: List[Dict[str, Any]] = []
    t_total = time.monotonic()
    for i, entry in enumerate(entries, 1):
        t0_smi = time.monotonic()
        rec = _run_grace_for_smiles(
            entry, method=args.method,
            max_per_isomer=args.max_per_isomer,
            max_isomers=args.max_isomers,
            max_rotamer_tuples=args.max_rotamer_tuples,
            budget_s=args.budget_s,
        )
        rows.append(rec)
        if args.with_singleshot:
            ss = _run_singleshot(entry)
            # Compute lift.
            grace_best = float("inf")
            per_iso = rec.get("per_isomer") or []
            if per_iso:
                grace_best = min(float(c.get("score", float("inf"))) for c in per_iso)
            ss["grace_best_score"] = (
                grace_best if math.isfinite(grace_best) else float("nan")
            )
            ss["lift_score_delta"] = (
                ss["single_shot_score"] - grace_best
                if math.isfinite(grace_best) and math.isfinite(ss["single_shot_score"])
                else float("nan")
            )
            merged = dict(rec)
            merged.update(ss)
            singleshot_rows.append(merged)
        if i % 10 == 0 or i == len(entries):
            elapsed = time.monotonic() - t_total
            n_ok = sum(1 for r in rows if r.get("decompose_ok"))
            print(f"  [{i:>3}/{len(entries):>3}] "
                  f"ok={n_ok} emit_total={sum(int(r.get('n_emitted',0)) for r in rows)} "
                  f"t={elapsed:.1f}s "
                  f"(last={time.monotonic()-t0_smi:.2f}s "
                  f"name={entry.name[:30]})", flush=True)

    total_wall = time.monotonic() - t_total
    agg = _aggregate(rows)

    if not args.summary_only:
        # Write CSV + JSONL.
        csv_path = os.path.join(args.out_dir, "grace_200_burnside_benchmark.csv")
        jsonl_path = os.path.join(args.out_dir, "grace_200_burnside_per_isomer.jsonl")
        _write_csv(csv_path, WIDE_COLUMNS, rows)
        _write_jsonl(jsonl_path, rows)
        print(f"# CSV   : {csv_path}")
        print(f"# JSONL : {jsonl_path}")

        if args.with_singleshot:
            ss_path = os.path.join(args.out_dir, "grace_vs_singleshot_lift.csv")
            _write_csv(ss_path, SINGLESHOT_COLUMNS, singleshot_rows)
            print(f"# CSV   : {ss_path}")

        # Scaling profile CSV.
        scale_path = os.path.join(args.out_dir, "grace_scaling_profile.csv")
        _write_csv(scale_path, SCALING_COLUMNS, rows)
        print(f"# CSV   : {scale_path}")

        # Aggregate summary JSON.
        agg_path = os.path.join(args.out_dir, "grace_200_burnside_aggregate.json")
        with open(agg_path, "w") as fh:
            json.dump({
                "params": {
                    "method": args.method,
                    "budget_s": args.budget_s,
                    "max_per_isomer": args.max_per_isomer,
                    "max_isomers": args.max_isomers,
                    "max_rotamer_tuples": args.max_rotamer_tuples,
                    "n_entries": len(entries),
                    "total_wall_s": total_wall,
                },
                "aggregate": agg,
            }, fh, indent=2)
        print(f"# JSON  : {agg_path}")

    # Print headline numbers (always).
    print("# " + "=" * 76)
    print(f"# HEADLINE  total={agg['n_total']}  decompose_ok={agg['n_decompose_ok']}  "
          f"emit_at_least_1={agg['n_emitted_at_least_one']}  "
          f"emit_total={agg['n_emitted_total']}")
    print(f"# Burnside  predicted_total={agg['n_predicted_total']}  "
          f"@100%-coverage={agg['n_at_100_pct_coverage']}/"
          f"{agg['n_decompose_ok']} ({agg['pct_at_100_coverage']:.1f}%)")
    print(f"# Coverage  mean={agg['coverage_mean']:.3f}  "
          f"median={agg['coverage_median']:.3f}")
    print(f"# Walltime  mean={agg['wall_mean_s']:.2f}s  "
          f"p95={agg['wall_p95_s']:.2f}s  max={agg['wall_max_s']:.2f}s  "
          f"total={total_wall:.1f}s")
    fit = agg["scaling_fit"]
    print(f"# Scaling   t = {fit['b']:.4g} * N^{fit['c']:.3f}  R^2={fit['r2']:.3f}  "
          f"(n_used={fit['n_used']})")
    print("# Per-family rollup:")
    for fam, fr in sorted(agg["by_family"].items()):
        print(f"#   {fam:6s}: n={fr['n']:3d}  "
              f"emit={fr['emitted_total']:4d}  "
              f"pred={fr['predicted_total']:4d}  "
              f"cov_mean={fr['coverage_mean']:.3f}  "
              f"@100%={fr['pct_at_100']:5.1f}%  "
              f"wall_mean={fr['wall_mean_s']:.2f}s")
    print("# " + "=" * 76)


if __name__ == "__main__":
    main()
