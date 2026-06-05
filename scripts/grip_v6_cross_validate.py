#!/usr/bin/env python3
"""Cross-validation: compare GRIP v5 (CCDC) and v6 (COD) fragment libraries.

For every fragment key present in BOTH libraries (pair_bond, triple_angle,
improper_pair, tm_categories, block-disagg), compute:

    diff_sigma = |mu_ccdc - mu_cod| / sqrt(0.5*(sigma_ccdc^2 + sigma_cod^2))

Classify each shared key:
    * both-agree     diff_sigma < 1
    * both-marginal  1 <= diff_sigma <= 2
    * both-disagree  diff_sigma > 2

Outputs:
    paper_data/grip_cod_vs_ccdc_agreement.json   summary statistics
    paper_data/grip_cod_vs_ccdc_agreement.csv    full per-key table

Usage::

    PYTHONHASHSEED=0 \\
    /home/qmchem_max/micromamba/envs/chemdarwin/bin/python \\
        scripts/grip_v6_cross_validate.py \\
        --ccdc reports/grip_lib_v5.npz \\
        --cod  reports/grip_lib_v6_cod.npz \\
        --out-dir paper_data
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

os.environ.setdefault("PYTHONHASHSEED", "0")

REPORTS = Path("/home/qmchem_max/agent_workspace/quality_framework/reports")
DEFAULT_CCDC = REPORTS / "grip_lib_v5.npz"
DEFAULT_COD = REPORTS / "grip_lib_v6_cod.npz"
DEFAULT_OUTDIR = Path("/home/qmchem_max/ComPlat/DELFIN/paper_data")

_AGREE_THRESHOLD = 1.0
_DISAGREE_THRESHOLD = 2.0


def _index_table(npz, keys_arr_name: str, mu_arr_name: str,
                 sigma_arr_name: str, n_arr_name: str
                 ) -> Dict[str, Tuple[float, float, int]]:
    """Index a (keys, mu, sigma, n) v5/v6 table into a dict."""
    if keys_arr_name not in npz.files:
        return {}
    keys = npz[keys_arr_name].tolist()
    mu = npz[mu_arr_name]
    sigma = npz[sigma_arr_name]
    n = npz[n_arr_name]
    out: Dict[str, Tuple[float, float, int]] = {}
    for i, k in enumerate(keys):
        m = float(mu[i])
        if not math.isfinite(m):
            continue
        s = float(sigma[i])
        nn = int(n[i])
        if nn < 5:
            continue
        out[str(k)] = (m, s, nn)
    return out


def _diff_sigma(mu1: float, sigma1: float, mu2: float, sigma2: float) -> float:
    pooled_var = 0.5 * (sigma1 * sigma1 + sigma2 * sigma2)
    pooled_sigma = max(math.sqrt(pooled_var), 1e-6)
    return abs(mu1 - mu2) / pooled_sigma


def _compare_table(
    ccdc_table: Dict[str, Tuple[float, float, int]],
    cod_table: Dict[str, Tuple[float, float, int]],
    category: str,
) -> Tuple[List[dict], dict]:
    """Compare two tables sharing the same key schema.

    Returns (per_key_rows, summary).
    """
    rows: List[dict] = []
    shared = sorted(set(ccdc_table.keys()) & set(cod_table.keys()))
    only_ccdc = sorted(set(ccdc_table.keys()) - set(cod_table.keys()))
    only_cod = sorted(set(cod_table.keys()) - set(ccdc_table.keys()))

    agree = marginal = disagree = 0
    diff_values: List[float] = []
    largest: List[Tuple[float, str, Tuple[float, float, int],
                        Tuple[float, float, int]]] = []
    for k in shared:
        mu_a, s_a, n_a = ccdc_table[k]
        mu_b, s_b, n_b = cod_table[k]
        d = _diff_sigma(mu_a, s_a, mu_b, s_b)
        if d < _AGREE_THRESHOLD:
            tag = "agree"
            agree += 1
        elif d > _DISAGREE_THRESHOLD:
            tag = "disagree"
            disagree += 1
        else:
            tag = "marginal"
            marginal += 1
        diff_values.append(d)
        largest.append((d, k, (mu_a, s_a, n_a), (mu_b, s_b, n_b)))
        rows.append({
            "category": category,
            "key": k,
            "mu_ccdc": round(mu_a, 6),
            "sigma_ccdc": round(s_a, 6),
            "n_ccdc": n_a,
            "mu_cod": round(mu_b, 6),
            "sigma_cod": round(s_b, 6),
            "n_cod": n_b,
            "diff_sigma": round(d, 4),
            "tag": tag,
        })

    largest.sort(key=lambda t: -t[0])
    top10 = [
        {
            "key": k,
            "diff_sigma": round(d, 3),
            "ccdc": {"mu": round(c[0], 3), "sigma": round(c[1], 3), "n": c[2]},
            "cod": {"mu": round(o[0], 3), "sigma": round(o[1], 3), "n": o[2]},
        }
        for (d, k, c, o) in largest[:10]
    ]
    summary = {
        "category": category,
        "n_shared": len(shared),
        "n_ccdc_only": len(only_ccdc),
        "n_cod_only": len(only_cod),
        "n_total_keys_seen": len(shared) + len(only_ccdc) + len(only_cod),
        "agree": agree,
        "marginal": marginal,
        "disagree": disagree,
        "pct_agree": round(100.0 * agree / max(len(shared), 1), 2),
        "pct_marginal": round(100.0 * marginal / max(len(shared), 1), 2),
        "pct_disagree": round(100.0 * disagree / max(len(shared), 1), 2),
        "median_diff_sigma": round(float(np.median(diff_values)), 4)
            if diff_values else None,
        "mean_diff_sigma": round(float(np.mean(diff_values)), 4)
            if diff_values else None,
        "p95_diff_sigma": round(float(np.percentile(diff_values, 95)), 4)
            if diff_values else None,
        "top10_disagreements": top10,
    }
    return rows, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ccdc", type=Path, default=DEFAULT_CCDC)
    parser.add_argument("--cod", type=Path, default=DEFAULT_COD)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    if not args.ccdc.exists():
        print(f"ERROR: CCDC library not found: {args.ccdc}", file=sys.stderr)
        return 1
    if not args.cod.exists():
        print(f"ERROR: COD library not found: {args.cod}", file=sys.stderr)
        return 1

    print(f"Loading CCDC v5 from {args.ccdc}")
    ccdc = np.load(args.ccdc, allow_pickle=True, mmap_mode="r")
    print(f"  version={int(ccdc['version'])}")
    print(f"Loading COD v6 from {args.cod}")
    cod = np.load(args.cod, allow_pickle=True, mmap_mode="r")
    print(f"  version={int(cod['version'])}")

    # Table-by-table cross-validation
    categories: List[Tuple[str, str, str, str, str]] = [
        # (label, keys_arr, mu_arr, sigma_arr, n_arr)
        ("pair_bond",        "pair_bond_keys",
         "pair_bond_mu", "pair_bond_sigma", "pair_bond_n"),
        ("triple_angle",     "triple_angle_keys",
         "triple_angle_mu", "triple_angle_sigma", "triple_angle_n"),
        ("improper_pair",    "improper_pair_keys",
         "improper_pair_mu", "improper_pair_sigma", "improper_pair_n"),
        ("tm_carbene",       "tm_carbene_keys",
         "tm_carbene_mu", "tm_carbene_sigma", "tm_carbene_n"),
        ("tm_hapto_eta2",    "tm_hapto_eta2_keys",
         "tm_hapto_eta2_mu", "tm_hapto_eta2_sigma", "tm_hapto_eta2_n"),
        ("tm_hapto_eta5",    "tm_hapto_eta5_keys",
         "tm_hapto_eta5_mu", "tm_hapto_eta5_sigma", "tm_hapto_eta5_n"),
        ("tm_hapto_eta6",    "tm_hapto_eta6_keys",
         "tm_hapto_eta6_mu", "tm_hapto_eta6_sigma", "tm_hapto_eta6_n"),
        ("tm_mu_bridge",     "tm_mu_bridge_keys",
         "tm_mu_bridge_mu", "tm_mu_bridge_sigma", "tm_mu_bridge_n"),
        ("tm_agostic",       "tm_agostic_keys",
         "tm_agostic_mu", "tm_agostic_sigma", "tm_agostic_n"),
        ("tm_ox_addition",   "tm_ox_addition_keys",
         "tm_ox_addition_mu", "tm_ox_addition_sigma", "tm_ox_addition_n"),
        ("tm_pair_bond_block", "tm_pair_bond_block_keys",
         "tm_pair_bond_block_mu", "tm_pair_bond_block_sigma",
         "tm_pair_bond_block_n"),
        ("tm_triple_angle_block", "tm_triple_angle_block_keys",
         "tm_triple_angle_block_mu", "tm_triple_angle_block_sigma",
         "tm_triple_angle_block_n"),
    ]

    all_rows: List[dict] = []
    per_category_summary: List[dict] = []
    totals = Counter()
    for (label, kname, muname, signame, nname) in categories:
        ccdc_table = _index_table(ccdc, kname, muname, signame, nname)
        cod_table = _index_table(cod, kname, muname, signame, nname)
        rows, summary = _compare_table(ccdc_table, cod_table, label)
        all_rows.extend(rows)
        per_category_summary.append(summary)
        totals["shared"] += summary["n_shared"]
        totals["ccdc_only"] += summary["n_ccdc_only"]
        totals["cod_only"] += summary["n_cod_only"]
        totals["agree"] += summary["agree"]
        totals["marginal"] += summary["marginal"]
        totals["disagree"] += summary["disagree"]
        print(
            f"  {label:24s} shared={summary['n_shared']:6d} "
            f"agree={summary['agree']:6d} marginal={summary['marginal']:6d} "
            f"disagree={summary['disagree']:6d} "
            f"({summary['pct_agree']:5.1f}% / {summary['pct_marginal']:5.1f}% / "
            f"{summary['pct_disagree']:5.1f}%)"
        )

    n_total_shared = totals["shared"]
    aggregate = {
        "schema_version": 1,
        "ccdc_path": str(args.ccdc),
        "cod_path": str(args.cod),
        "n_shared_total": totals["shared"],
        "n_ccdc_only_total": totals["ccdc_only"],
        "n_cod_only_total": totals["cod_only"],
        "agree_total": totals["agree"],
        "marginal_total": totals["marginal"],
        "disagree_total": totals["disagree"],
        "pct_agree_total": round(100.0 * totals["agree"] / max(n_total_shared, 1), 2),
        "pct_marginal_total": round(100.0 * totals["marginal"] / max(n_total_shared, 1), 2),
        "pct_disagree_total": round(100.0 * totals["disagree"] / max(n_total_shared, 1), 2),
        "agree_threshold_sigma": _AGREE_THRESHOLD,
        "disagree_threshold_sigma": _DISAGREE_THRESHOLD,
        "min_n_per_key": 5,
        "per_category": per_category_summary,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "grip_cod_vs_ccdc_agreement.json"
    csv_path = args.out_dir / "grip_cod_vs_ccdc_agreement.csv"
    with open(json_path, "w") as fh:
        json.dump(aggregate, fh, indent=2, sort_keys=True)
    print(f"\nWrote {json_path}")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "category", "key", "mu_ccdc", "sigma_ccdc", "n_ccdc",
            "mu_cod", "sigma_cod", "n_cod", "diff_sigma", "tag"])
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"Wrote {csv_path}  ({len(all_rows)} rows)")

    print("\n=== AGGREGATE ===")
    print(f"  shared keys     : {totals['shared']}")
    print(f"  ccdc-only keys  : {totals['ccdc_only']}")
    print(f"  cod-only keys   : {totals['cod_only']}")
    print(f"  agree           : {totals['agree']:6d}  "
          f"({aggregate['pct_agree_total']:.2f}%)")
    print(f"  marginal        : {totals['marginal']:6d}  "
          f"({aggregate['pct_marginal_total']:.2f}%)")
    print(f"  disagree        : {totals['disagree']:6d}  "
          f"({aggregate['pct_disagree_total']:.2f}%)")
    print(f"\nverdict: COD vs CCDC agreement = "
          f"{aggregate['pct_agree_total']:.1f}% "
          f"(target: >=90% would validate the COD lib for production merge)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
