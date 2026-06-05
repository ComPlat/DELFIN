#!/usr/bin/env python3
"""Detailed breakdown of v6 COD vs v5_full CCDC GRIP cross-validation.

Reads paper_data/grip_cod_vs_ccdc_agreement.csv and computes:
  * Per-fragment-class summary (bonds, angles, impropers, torsions, TM)
  * Per-metal-block summary (3d, 4d, 5d, f-block, main-group)
  * Top-20 worst disagreements + hypothesis for each
  * improper_pair forensics: how many disagreements driven by tiny ccdc-n vs cod-n?

Usage::
    PYTHONHASHSEED=0 \\
    /home/qmchem_max/micromamba/envs/delfin/bin/python \\
        scripts/grip_v6_cross_validate_breakdown.py \\
        --csv  paper_data/grip_cod_vs_ccdc_agreement.csv \\
        --json paper_data/grip_cod_vs_ccdc_agreement.json \\
        --out-md paper_data/grip_cross_validation_v5full_summary.md
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

os.environ.setdefault("PYTHONHASHSEED", "0")

# Periodic-table classification of metals.
# Source: standard chemistry — 3d = Sc..Zn, 4d = Y..Cd, 5d = Hf..Hg (+ La), f = Ce..Lu, Th..Lr
ELEMENTS_3D = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
ELEMENTS_4D = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
ELEMENTS_5D = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}
ELEMENTS_F = {
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr",
}
ELEMENTS_MG_METALS = {"Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
                       "Al", "Ga", "In", "Tl", "Sn", "Pb", "Sb", "Bi"}


def _classify_metal(elements: List[str]) -> str:
    """Classify a fragment by which metal block(s) it touches."""
    blocks = set()
    for e in elements:
        if e in ELEMENTS_3D:
            blocks.add("3d")
        elif e in ELEMENTS_4D:
            blocks.add("4d")
        elif e in ELEMENTS_5D:
            blocks.add("5d")
        elif e in ELEMENTS_F:
            blocks.add("f")
        elif e in ELEMENTS_MG_METALS:
            blocks.add("main_group_metal")
    if not blocks:
        return "non_metal"
    if len(blocks) == 1:
        return next(iter(blocks))
    return "mixed_metal"


def _extract_elements(key: str) -> List[str]:
    """Parse element symbols from a JSON-list-encoded fragment key.

    Keys look like ["Pd","sp2","Cl","sp3"] or ["C","sp2","C","sp2","C","sp2"].
    """
    # Match capitalized 1- or 2-char element symbols inside quotes.
    elements = re.findall(r'"([A-Z][a-z]?)"', key)
    return elements


def _hypothesize(row: dict, category: str) -> str:
    """Generate a plain-English hypothesis for why a specific row disagrees."""
    try:
        mu_a = float(row["mu_ccdc"])
        mu_b = float(row["mu_cod"])
        s_a = float(row["sigma_ccdc"])
        s_b = float(row["sigma_cod"])
        n_a = int(row["n_ccdc"])
        n_b = int(row["n_cod"])
    except (KeyError, ValueError):
        return "unparseable-row"

    notes = []
    # Small-n in either lib → noisy
    if n_a < 20 or n_b < 20:
        notes.append("low-n")
    # n imbalance > 10x
    if max(n_a, n_b) >= 10 * max(min(n_a, n_b), 1):
        notes.append("n-imbalance")
    # Very tight sigma in either → small denominator inflates diff_sigma
    if min(s_a, s_b) < 0.01:
        notes.append("tight-sigma")
    # Large absolute mu difference
    if category == "improper_pair":
        # Impropers are bond-pair angles around an atom; sensitive to bond detection
        notes.append("bond-detection-edge-case")
    if category == "pair_bond":
        delta = abs(mu_a - mu_b)
        if delta > 0.3:
            notes.append(f"|delta_mu|={delta:.2f}A")
    if category.startswith("tm_") or "_block" in category:
        if "n-imbalance" in notes or "low-n" in notes:
            notes.append("rare-metal-coord-class")
    if not notes:
        notes.append("genuine-database-disagreement")
    return ";".join(notes)


def load_rows(csv_path: Path) -> List[dict]:
    with open(csv_path, "r") as fh:
        r = csv.DictReader(fh)
        return list(r)


def per_class_summary(rows: List[dict]) -> List[dict]:
    by_class = defaultdict(lambda: Counter())
    for row in rows:
        cat = row["category"]
        tag = row["tag"]
        by_class[cat]["total"] += 1
        by_class[cat][tag] += 1
    out = []
    for cat in sorted(by_class.keys()):
        c = by_class[cat]
        total = c["total"]
        out.append({
            "category": cat,
            "shared": total,
            "agree": c["agree"],
            "marginal": c["marginal"],
            "disagree": c["disagree"],
            "pct_agree": round(100.0 * c["agree"] / max(total, 1), 2),
            "pct_marginal": round(100.0 * c["marginal"] / max(total, 1), 2),
            "pct_disagree": round(100.0 * c["disagree"] / max(total, 1), 2),
        })
    return out


def per_metal_block_summary(rows: List[dict]) -> List[dict]:
    by_block = defaultdict(lambda: Counter())
    for row in rows:
        elements = _extract_elements(row["key"])
        block = _classify_metal(elements)
        tag = row["tag"]
        by_block[block]["total"] += 1
        by_block[block][tag] += 1
    out = []
    for block in sorted(by_block.keys()):
        c = by_block[block]
        total = c["total"]
        out.append({
            "block": block,
            "shared": total,
            "agree": c["agree"],
            "marginal": c["marginal"],
            "disagree": c["disagree"],
            "pct_agree": round(100.0 * c["agree"] / max(total, 1), 2),
            "pct_marginal": round(100.0 * c["marginal"] / max(total, 1), 2),
            "pct_disagree": round(100.0 * c["disagree"] / max(total, 1), 2),
        })
    return out


def top_n_disagreements(rows: List[dict], n: int = 20) -> List[dict]:
    typed = []
    for row in rows:
        try:
            d = float(row["diff_sigma"])
        except (KeyError, ValueError):
            continue
        typed.append((d, row))
    typed.sort(key=lambda t: -t[0])
    return [
        {
            "rank": i + 1,
            "category": r["category"],
            "key": r["key"],
            "mu_ccdc": float(r["mu_ccdc"]),
            "sigma_ccdc": float(r["sigma_ccdc"]),
            "n_ccdc": int(r["n_ccdc"]),
            "mu_cod": float(r["mu_cod"]),
            "sigma_cod": float(r["sigma_cod"]),
            "n_cod": int(r["n_cod"]),
            "diff_sigma": float(r["diff_sigma"]),
            "hypothesis": _hypothesize(r, r["category"]),
        }
        for i, (_, r) in enumerate(typed[:n])
    ]


def improper_pair_forensics(rows: List[dict]) -> dict:
    """Validate the 'improper_pair 36% disagree' hypothesis from smoke."""
    imp = [r for r in rows if r["category"] == "improper_pair"]
    disagree = [r for r in imp if r["tag"] == "disagree"]
    if not disagree:
        return {"n_total": len(imp), "n_disagree": 0}
    low_n = [r for r in disagree if int(r["n_ccdc"]) < 20 or int(r["n_cod"]) < 20]
    tight_sigma = [r for r in disagree
                   if float(r["sigma_ccdc"]) < 0.01 or float(r["sigma_cod"]) < 0.01]
    n_imbal = [r for r in disagree
               if max(int(r["n_ccdc"]), int(r["n_cod"])) >=
                  10 * max(min(int(r["n_ccdc"]), int(r["n_cod"])), 1)]

    # Look at avg n on disagreements vs agreements
    agree = [r for r in imp if r["tag"] == "agree"]
    avg_n_ccdc_dis = sum(int(r["n_ccdc"]) for r in disagree) / max(len(disagree), 1)
    avg_n_cod_dis = sum(int(r["n_cod"]) for r in disagree) / max(len(disagree), 1)
    avg_n_ccdc_agr = sum(int(r["n_ccdc"]) for r in agree) / max(len(agree), 1)
    avg_n_cod_agr = sum(int(r["n_cod"]) for r in agree) / max(len(agree), 1)

    avg_sigma_ccdc_dis = sum(float(r["sigma_ccdc"]) for r in disagree) / max(len(disagree), 1)
    avg_sigma_cod_dis = sum(float(r["sigma_cod"]) for r in disagree) / max(len(disagree), 1)
    avg_sigma_ccdc_agr = sum(float(r["sigma_ccdc"]) for r in agree) / max(len(agree), 1)
    avg_sigma_cod_agr = sum(float(r["sigma_cod"]) for r in agree) / max(len(agree), 1)

    return {
        "n_total": len(imp),
        "n_disagree": len(disagree),
        "n_disagree_low_n": len(low_n),
        "n_disagree_tight_sigma": len(tight_sigma),
        "n_disagree_n_imbalance": len(n_imbal),
        "pct_disagree_low_n": round(100.0 * len(low_n) / max(len(disagree), 1), 2),
        "pct_disagree_tight_sigma": round(100.0 * len(tight_sigma) / max(len(disagree), 1), 2),
        "pct_disagree_n_imbalance": round(100.0 * len(n_imbal) / max(len(disagree), 1), 2),
        "avg_n_ccdc_agree": round(avg_n_ccdc_agr, 1),
        "avg_n_cod_agree": round(avg_n_cod_agr, 1),
        "avg_n_ccdc_disagree": round(avg_n_ccdc_dis, 1),
        "avg_n_cod_disagree": round(avg_n_cod_dis, 1),
        "avg_sigma_ccdc_agree": round(avg_sigma_ccdc_agr, 4),
        "avg_sigma_cod_agree": round(avg_sigma_cod_agr, 4),
        "avg_sigma_ccdc_disagree": round(avg_sigma_ccdc_dis, 4),
        "avg_sigma_cod_disagree": round(avg_sigma_cod_dis, 4),
    }


def smoke_vs_full_comparison(smoke_summary: dict, full_summary: dict) -> List[dict]:
    smoke_per_cat = {c["category"]: c for c in smoke_summary["per_category"]}
    full_per_cat = {c["category"]: c for c in full_summary["per_category"]}
    rows = []
    rows.append({
        "category": "TOTAL",
        "n_shared_smoke": smoke_summary["n_shared_total"],
        "n_shared_full": full_summary["n_shared_total"],
        "pct_agree_smoke": smoke_summary["pct_agree_total"],
        "pct_agree_full": full_summary["pct_agree_total"],
        "delta_agree_pp": round(full_summary["pct_agree_total"] -
                                smoke_summary["pct_agree_total"], 2),
        "pct_disagree_smoke": smoke_summary["pct_disagree_total"],
        "pct_disagree_full": full_summary["pct_disagree_total"],
        "delta_disagree_pp": round(full_summary["pct_disagree_total"] -
                                   smoke_summary["pct_disagree_total"], 2),
    })
    for cat in sorted(set(smoke_per_cat) | set(full_per_cat)):
        s = smoke_per_cat.get(cat, {})
        f = full_per_cat.get(cat, {})
        rows.append({
            "category": cat,
            "n_shared_smoke": s.get("n_shared", 0),
            "n_shared_full": f.get("n_shared", 0),
            "pct_agree_smoke": s.get("pct_agree", 0.0),
            "pct_agree_full": f.get("pct_agree", 0.0),
            "delta_agree_pp": round(f.get("pct_agree", 0.0) - s.get("pct_agree", 0.0), 2),
            "pct_disagree_smoke": s.get("pct_disagree", 0.0),
            "pct_disagree_full": f.get("pct_disagree", 0.0),
            "delta_disagree_pp": round(f.get("pct_disagree", 0.0)
                                       - s.get("pct_disagree", 0.0), 2),
        })
    return rows


def write_markdown(out_md: Path,
                   full_summary: dict,
                   smoke_summary: dict,
                   per_class: List[dict],
                   per_block: List[dict],
                   top20: List[dict],
                   improper_forensics: dict,
                   smoke_v_full: List[dict]) -> None:
    lines: List[str] = []
    lines.append("# GRIP v6 (COD) vs v5_full (CCDC) Cross-Validation Summary")
    lines.append("")
    lines.append(f"Date: 2026-06-05  ")
    lines.append(f"CCDC source: `{full_summary['ccdc_path']}`  ")
    lines.append(f"COD source: `{full_summary['cod_path']}`  ")
    lines.append(f"Agreement threshold: |mu_a - mu_b| / sqrt(0.5*(sigma_a^2 + sigma_b^2)) "
                 f"< {full_summary['agree_threshold_sigma']}  ")
    lines.append(f"Disagree threshold: > {full_summary['disagree_threshold_sigma']}  ")
    lines.append(f"Minimum n per key: {full_summary['min_n_per_key']}")
    lines.append("")

    # Headline numbers
    lines.append("## Headline")
    lines.append("")
    lines.append(f"- Total shared keys: **{full_summary['n_shared_total']:,}**")
    lines.append(f"- CCDC-only keys: {full_summary['n_ccdc_only_total']:,} "
                 "(features the COD lib does not see)")
    lines.append(f"- COD-only keys: {full_summary['n_cod_only_total']:,} "
                 "(features only COD encountered)")
    lines.append(f"- Both-agree: **{full_summary['agree_total']:,} "
                 f"({full_summary['pct_agree_total']:.2f}%)**")
    lines.append(f"- Both-marginal: {full_summary['marginal_total']:,} "
                 f"({full_summary['pct_marginal_total']:.2f}%)")
    lines.append(f"- Both-disagree: {full_summary['disagree_total']:,} "
                 f"({full_summary['pct_disagree_total']:.2f}%)")
    lines.append("")
    lines.append(f"Verdict: COD library reproduces **{full_summary['pct_agree_total']:.1f}%** "
                 "of CCDC-derived means within 1 pooled sigma.")
    lines.append("")

    # Smoke vs full
    lines.append("## Smoke (5k entries) vs Full (1.44M entries) v5 comparison")
    lines.append("")
    lines.append("| Category | n_shared smoke | n_shared full | agree% smoke | agree% full | "
                 "Delta agree pp | disagree% smoke | disagree% full | Delta disagree pp |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in smoke_v_full:
        lines.append(f"| {r['category']} | {r['n_shared_smoke']:,} | {r['n_shared_full']:,} | "
                     f"{r['pct_agree_smoke']:.2f} | {r['pct_agree_full']:.2f} | "
                     f"{r['delta_agree_pp']:+.2f} | "
                     f"{r['pct_disagree_smoke']:.2f} | {r['pct_disagree_full']:.2f} | "
                     f"{r['delta_disagree_pp']:+.2f} |")
    lines.append("")

    # Per-fragment-class
    lines.append("## Per-fragment-class breakdown (full v5)")
    lines.append("")
    lines.append("| Category | Shared | Agree | Marginal | Disagree | %Agree | %Marginal | %Disagree |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for c in per_class:
        lines.append(
            f"| {c['category']} | {c['shared']:,} | {c['agree']:,} | {c['marginal']:,} | "
            f"{c['disagree']:,} | {c['pct_agree']:.2f} | {c['pct_marginal']:.2f} | "
            f"{c['pct_disagree']:.2f} |"
        )
    lines.append("")

    # Per-metal-block
    lines.append("## Per-metal-block breakdown (full v5)")
    lines.append("")
    lines.append("Classification by all element symbols found in the fragment key. "
                 "`non_metal` = pure ligand fragment with no metal; "
                 "`mixed_metal` = key touches >=2 different blocks.")
    lines.append("")
    lines.append("| Block | Shared | Agree | Marginal | Disagree | %Agree | %Marginal | %Disagree |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for b in per_block:
        lines.append(
            f"| {b['block']} | {b['shared']:,} | {b['agree']:,} | {b['marginal']:,} | "
            f"{b['disagree']:,} | {b['pct_agree']:.2f} | {b['pct_marginal']:.2f} | "
            f"{b['pct_disagree']:.2f} |"
        )
    lines.append("")

    # improper_pair forensics
    lines.append("## improper_pair forensics (validating the smoke 36% disagree finding)")
    lines.append("")
    f = improper_forensics
    lines.append(f"- Smoke run reported **improper_pair 36% disagree** (177/491). "
                 f"Full v5 result: **{f['n_disagree']}/{f['n_total']} disagree "
                 f"= {100.0 * f['n_disagree'] / max(f['n_total'], 1):.2f}%**")
    lines.append(f"- Of those disagreeing keys: low-n (<20 in either lib) = "
                 f"{f['n_disagree_low_n']} ({f['pct_disagree_low_n']:.1f}%)")
    lines.append(f"- Of those disagreeing keys: tight sigma (<0.01 in either lib) = "
                 f"{f['n_disagree_tight_sigma']} ({f['pct_disagree_tight_sigma']:.1f}%)")
    lines.append(f"- Of those disagreeing keys: n imbalance (>=10x ratio) = "
                 f"{f['n_disagree_n_imbalance']} ({f['pct_disagree_n_imbalance']:.1f}%)")
    lines.append("")
    lines.append("Average sample size and sigma, agree vs disagree:")
    lines.append("")
    lines.append("| | n_ccdc | n_cod | sigma_ccdc | sigma_cod |")
    lines.append("|---|---:|---:|---:|---:|")
    lines.append(f"| agree | {f['avg_n_ccdc_agree']} | {f['avg_n_cod_agree']} | "
                 f"{f['avg_sigma_ccdc_agree']} | {f['avg_sigma_cod_agree']} |")
    lines.append(f"| disagree | {f['avg_n_ccdc_disagree']} | {f['avg_n_cod_disagree']} | "
                 f"{f['avg_sigma_ccdc_disagree']} | {f['avg_sigma_cod_disagree']} |")
    lines.append("")

    # Top 20 disagreements
    lines.append("## Top-20 worst disagreements (highest diff_sigma)")
    lines.append("")
    lines.append("| Rank | Category | Key | mu_ccdc | n_ccdc | mu_cod | n_cod | diff_sigma | Hypothesis |")
    lines.append("|---:|---|---|---:|---:|---:|---:|---:|---|")
    for t in top20:
        key_short = t["key"] if len(t["key"]) < 60 else t["key"][:57] + "..."
        lines.append(
            f"| {t['rank']} | {t['category']} | `{key_short}` | "
            f"{t['mu_ccdc']:.3f} | {t['n_ccdc']} | {t['mu_cod']:.3f} | {t['n_cod']} | "
            f"{t['diff_sigma']:.2f} | {t['hypothesis']} |"
        )
    lines.append("")

    out_md.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--smoke-json", type=Path, required=True,
                        help="Old smoke JSON for comparison")
    parser.add_argument("--out-md", type=Path, required=True)
    parser.add_argument("--out-breakdown-json", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.csv)
    full_summary = json.loads(args.json.read_text())
    smoke_summary = json.loads(args.smoke_json.read_text())

    per_class = per_class_summary(rows)
    per_block = per_metal_block_summary(rows)
    top20 = top_n_disagreements(rows, 20)
    improper_forensics = improper_pair_forensics(rows)
    smoke_v_full = smoke_vs_full_comparison(smoke_summary, full_summary)

    write_markdown(args.out_md, full_summary, smoke_summary,
                   per_class, per_block, top20, improper_forensics, smoke_v_full)
    print(f"Wrote {args.out_md}")

    breakdown = {
        "schema_version": 1,
        "date": "2026-06-05",
        "ccdc_path": full_summary["ccdc_path"],
        "cod_path": full_summary["cod_path"],
        "headline": {
            "n_shared_total": full_summary["n_shared_total"],
            "n_ccdc_only_total": full_summary["n_ccdc_only_total"],
            "n_cod_only_total": full_summary["n_cod_only_total"],
            "pct_agree_total": full_summary["pct_agree_total"],
            "pct_marginal_total": full_summary["pct_marginal_total"],
            "pct_disagree_total": full_summary["pct_disagree_total"],
        },
        "smoke_vs_full": smoke_v_full,
        "per_fragment_class": per_class,
        "per_metal_block": per_block,
        "improper_pair_forensics": improper_forensics,
        "top20_disagreements": top20,
    }
    with open(args.out_breakdown_json, "w") as fh:
        json.dump(breakdown, fh, indent=2, sort_keys=True)
    print(f"Wrote {args.out_breakdown_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
