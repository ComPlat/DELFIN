#!/usr/bin/env python3
"""Mission G1' Phase 3 — per-class aggregator over xrd_rmsd_*.jsonl files.

Combines per-builder JSONL outputs from xrd_rmsd_comparator.py into:
  - xrd_match_table.csv         (paper_data) — one row per (refcode, builder)
  - xrd_rmsd_per_class.{csv,md} — per-class mean/median heavy/bond/angle RMSD
  - xrd_rmsd_summary.json       — overall + per-class + headline claim

Classification axes (from CCDC index):
  - block  (d/f/p/s)
  - cn     (coordination number)
  - hapto  (organometallic flag from CCDC index; we coarsely tag η-class as
            organometallic_with_C_donor when CCDC labels organometallic)
  - bucket (composite: CN<n>-<block>; e.g. CN6-d, CN3-d)

Headline claim is computed as percentage improvement of best-DELFIN-builder
over UFF on the SHARED refcode subset (paired comparison, honest).

USAGE
    python scripts/xrd_rmsd_aggregate.py \
        --in paper_data/xrd_rmsd_uff.jsonl \
        --in paper_data/xrd_rmsd_delfin_aromatic.jsonl \
        --in paper_data/xrd_rmsd_delfin_fffree.jsonl \
        --out-dir paper_data
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def load_jsonl(path: Path) -> List[dict]:
    out: List[dict] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except Exception:
                continue
    return out


def class_bucket(rec: dict) -> str:
    cn = rec.get("cn")
    block = rec.get("block")
    if cn is None or block is None:
        return "unknown"
    return f"CN{cn}-{block}"


def hapto_class(rec: dict) -> str:
    """Coarse organometallic vs Werner classification using CCDC index."""
    is_om = rec.get("is_organometallic")
    if is_om is True:
        return "organometallic"
    if is_om is False:
        return "werner"
    return "unknown"


def coarse_cn(cn) -> str:
    if cn is None:
        return "CN?"
    try:
        n = int(cn)
    except Exception:
        return "CN?"
    if n <= 2:
        return "CN1-2"
    if n <= 4:
        return "CN3-4"
    if n <= 6:
        return "CN5-6"
    if n <= 8:
        return "CN7-8"
    return "CN9+"


def metric_summary(vals: List[float]) -> Dict[str, float]:
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    if not vals:
        return {"n": 0, "mean": float("nan"), "median": float("nan"),
                "min": float("nan"), "max": float("nan")}
    return {
        "n": len(vals),
        "mean": float(statistics.fmean(vals)),
        "median": float(statistics.median(vals)),
        "min": float(min(vals)),
        "max": float(max(vals)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inputs", action="append", required=True,
                    help="Per-builder JSONL (use multiple --in)")
    ap.add_argument("--out-dir", default="paper_data")
    ap.add_argument("--date-tag", default="2026_06_06")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load all
    by_label: Dict[str, Dict[str, dict]] = defaultdict(dict)
    all_recs: List[dict] = []
    for fp_str in args.inputs:
        fp = Path(fp_str)
        recs = load_jsonl(fp)
        for r in recs:
            label = r.get("label", "unknown")
            rc = r.get("refcode")
            if not rc:
                continue
            by_label[label][rc] = r
            all_recs.append(r)
        print(f"[aggregate] loaded {len(recs)} from {fp.name} (label={recs[0].get('label') if recs else '?'})")

    labels = sorted(by_label.keys())
    print(f"[aggregate] labels: {labels}")

    # ---- xrd_match_table.csv ----
    all_refcodes = sorted({r["refcode"] for r in all_recs})
    match_csv = out_dir / "xrd_match_table.csv"
    fields = ["refcode", "metal", "block", "cn", "geom", "is_organometallic",
              "n_atoms_ccdc", "n_atoms_raw"]
    # Per builder columns
    builder_fields: List[str] = []
    for lab in labels:
        for k in ("heavy_rmsd", "bond_rmsd", "angle_rmsd",
                  "n_match", "n_emitted", "candidate_file"):
            builder_fields.append(f"{lab}__{k}")
    with match_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields + builder_fields)
        for rc in all_refcodes:
            # Use first builder that has it for classification
            base = None
            for lab in labels:
                if rc in by_label[lab]:
                    base = by_label[lab][rc]
                    break
            if base is None:
                continue
            row = [rc, base.get("metal"), base.get("block"), base.get("cn"),
                   base.get("geom"), base.get("is_organometallic"),
                   base.get("n_atoms_ccdc"), base.get("n_atoms_raw")]
            for lab in labels:
                r = by_label[lab].get(rc, {})
                for k in ("heavy_rmsd", "bond_rmsd", "angle_rmsd",
                          "n_match", "n_emitted", "candidate_file"):
                    row.append(r.get(k, ""))
            w.writerow(row)
    print(f"[aggregate] wrote {match_csv}")

    # ---- per-class table ----
    # Group by cn-bucket + block, by builder
    class_keys = sorted({
        coarse_cn(r.get("cn"))
        for r in all_recs
    })
    block_keys = sorted({
        r.get("block") or "?"
        for r in all_recs
    })
    composite_buckets = sorted({
        f"{coarse_cn(r.get('cn'))} / {r.get('block') or '?'}-block"
        for r in all_recs
    })

    # For each composite bucket and label, gather heavy/bond/angle vals
    bucket_data: Dict[str, Dict[str, List[Tuple[float, float, float]]]] = (
        defaultdict(lambda: defaultdict(list))
    )
    for r in all_recs:
        bucket = f"{coarse_cn(r.get('cn'))} / {r.get('block') or '?'}-block"
        lab = r.get("label")
        h = r.get("heavy_rmsd")
        b = r.get("bond_rmsd")
        a = r.get("angle_rmsd")
        if lab is None:
            continue
        bucket_data[bucket][lab].append((h, b, a))

    # Also add organometallic/Werner buckets
    om_bucket_data: Dict[str, Dict[str, List[Tuple[float, float, float]]]] = (
        defaultdict(lambda: defaultdict(list))
    )
    for r in all_recs:
        bucket = hapto_class(r)
        lab = r.get("label")
        if lab is None:
            continue
        om_bucket_data[bucket][lab].append(
            (r.get("heavy_rmsd"), r.get("bond_rmsd"), r.get("angle_rmsd"))
        )

    # Write CSV
    perclass_csv = out_dir / f"xrd_rmsd_per_class_{args.date_tag}.csv"
    with perclass_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        header = ["bucket", "label", "n",
                  "heavy_rmsd_mean", "heavy_rmsd_median",
                  "bond_rmsd_mean", "bond_rmsd_median",
                  "angle_rmsd_mean", "angle_rmsd_median"]
        w.writerow(header)
        for bucket in sorted(bucket_data.keys()):
            for lab in labels:
                vals = bucket_data[bucket].get(lab, [])
                if not vals:
                    continue
                hs = [v[0] for v in vals]
                bs = [v[1] for v in vals]
                ang = [v[2] for v in vals]
                hsum = metric_summary(hs)
                bsum = metric_summary(bs)
                asum = metric_summary(ang)
                w.writerow([
                    bucket, lab, hsum["n"],
                    round(hsum["mean"], 4), round(hsum["median"], 4),
                    round(bsum["mean"], 4), round(bsum["median"], 4),
                    round(asum["mean"], 4), round(asum["median"], 4),
                ])
        # organometallic / werner
        for bucket in sorted(om_bucket_data.keys()):
            for lab in labels:
                vals = om_bucket_data[bucket].get(lab, [])
                if not vals:
                    continue
                hs = [v[0] for v in vals]
                bs = [v[1] for v in vals]
                ang = [v[2] for v in vals]
                hsum = metric_summary(hs)
                bsum = metric_summary(bs)
                asum = metric_summary(ang)
                w.writerow([
                    f"_chem_{bucket}", lab, hsum["n"],
                    round(hsum["mean"], 4), round(hsum["median"], 4),
                    round(bsum["mean"], 4), round(bsum["median"], 4),
                    round(asum["mean"], 4), round(asum["median"], 4),
                ])
    print(f"[aggregate] wrote {perclass_csv}")

    # ---- per-class MD table ----
    perclass_md = out_dir / f"xrd_rmsd_per_class_{args.date_tag}.md"
    with perclass_md.open("w") as fh:
        fh.write("# XRD-RMSD Per-Class Comparison (Mission G1')\n\n")
        fh.write("Comparison of DELFIN vs UFF outputs against CCDC XRD-resolved\n"
                 "structures (ground truth). Each cell reports mean (n) heavy-atom\n"
                 "Kabsch RMSD over best-emission-per-refcode within each pool.\n\n")
        fh.write("**CCDC ground-truth source**: "
                 "`agent_workspace/quality_framework/reference/ccdc_cleaned_xyz/` "
                 "(997 cleaned refcodes, counter-ions stripped, single-component, "
                 "disorder resolved).\n\n")
        fh.write("**Pool sources**:\n")
        for lab in labels:
            n_rc = len(by_label[lab])
            fh.write(f"- `{lab}` — {n_rc} refcodes scored\n")
        fh.write("\n")

        # Overall table
        fh.write("## Overall (all matched refcodes per builder)\n\n")
        fh.write("| builder | n | heavy-RMSD (Å) mean | heavy-RMSD median | "
                 "bond-RMSD (Å) mean | angle-RMSD (deg) mean |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for lab in labels:
            recs = list(by_label[lab].values())
            hs = [r["heavy_rmsd"] for r in recs]
            bs = [r["bond_rmsd"] for r in recs]
            ang = [r["angle_rmsd"] for r in recs]
            hs_s = metric_summary(hs)
            bs_s = metric_summary(bs)
            ang_s = metric_summary(ang)
            fh.write(f"| {lab} | {hs_s['n']} | "
                     f"{hs_s['mean']:.3f} | {hs_s['median']:.3f} | "
                     f"{bs_s['mean']:.3f} | {ang_s['mean']:.2f} |\n")

        # CN x block table
        fh.write("\n## Per CN-coarse × block bucket\n\n")
        fh.write("| bucket | builder | n | heavy-RMSD (Å) | "
                 "bond-RMSD (Å) | angle-RMSD (°) |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for bucket in sorted(bucket_data.keys()):
            for lab in labels:
                vals = bucket_data[bucket].get(lab, [])
                if not vals:
                    continue
                hs = [v[0] for v in vals]
                bs = [v[1] for v in vals]
                ang = [v[2] for v in vals]
                hs_s = metric_summary(hs)
                bs_s = metric_summary(bs)
                ang_s = metric_summary(ang)
                fh.write(f"| {bucket} | {lab} | {hs_s['n']} | "
                         f"{hs_s['mean']:.3f} | {bs_s['mean']:.3f} | "
                         f"{ang_s['mean']:.2f} |\n")

        # Organometallic vs Werner
        fh.write("\n## Organometallic vs Werner (CCDC `is_organometallic` flag)\n\n")
        fh.write("| class | builder | n | heavy-RMSD (Å) | "
                 "bond-RMSD (Å) | angle-RMSD (°) |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for bucket in sorted(om_bucket_data.keys()):
            for lab in labels:
                vals = om_bucket_data[bucket].get(lab, [])
                if not vals:
                    continue
                hs = [v[0] for v in vals]
                bs = [v[1] for v in vals]
                ang = [v[2] for v in vals]
                hs_s = metric_summary(hs)
                bs_s = metric_summary(bs)
                ang_s = metric_summary(ang)
                fh.write(f"| {bucket} | {lab} | {hs_s['n']} | "
                         f"{hs_s['mean']:.3f} | {bs_s['mean']:.3f} | "
                         f"{ang_s['mean']:.2f} |\n")

        # Paired-comparison headline (only refcodes covered by BOTH UFF and
        # any DELFIN pool — restricted subset for fairness)
        fh.write("\n## Paired comparison (shared refcodes — fair head-to-head)\n\n")
        uff_label = next((l for l in labels if "uff" in l.lower()), None)
        delfin_labels = [l for l in labels if l != uff_label]
        if uff_label and delfin_labels:
            uff_set = set(by_label[uff_label].keys())
            for dlab in delfin_labels:
                dset = set(by_label[dlab].keys())
                shared = sorted(uff_set & dset)
                if not shared:
                    continue
                u_h = [by_label[uff_label][rc]["heavy_rmsd"] for rc in shared]
                d_h = [by_label[dlab][rc]["heavy_rmsd"] for rc in shared]
                u_b = [by_label[uff_label][rc]["bond_rmsd"] for rc in shared]
                d_b = [by_label[dlab][rc]["bond_rmsd"] for rc in shared]
                wins = sum(1 for i in range(len(shared))
                           if d_h[i] < u_h[i])
                fh.write(f"### {dlab} vs {uff_label} (n_shared = {len(shared)})\n\n")
                fh.write(f"- {dlab} heavy-RMSD mean: "
                         f"{statistics.fmean(d_h):.3f} Å\n")
                fh.write(f"- {uff_label} heavy-RMSD mean: "
                         f"{statistics.fmean(u_h):.3f} Å\n")
                pct = (1.0 - statistics.fmean(d_h)
                       / statistics.fmean(u_h)) * 100
                fh.write(f"- **DELFIN improvement: {pct:+.1f}%** "
                         f"(per-refcode wins: {wins}/{len(shared)})\n")
                fh.write(f"- shared refcodes: "
                         f"{', '.join(shared[:20])}"
                         f"{'...' if len(shared) > 20 else ''}\n\n")
        else:
            fh.write("(no UFF or no DELFIN labels — paired comparison "
                     "not applicable)\n")

        # Paired comparison BETWEEN DELFIN variants
        if len(delfin_labels) >= 2:
            fh.write("\n## Paired comparison between DELFIN variants\n\n")
            for i, lab_a in enumerate(delfin_labels):
                for lab_b in delfin_labels[i + 1:]:
                    set_a = set(by_label[lab_a].keys())
                    set_b = set(by_label[lab_b].keys())
                    shared = sorted(set_a & set_b)
                    if not shared:
                        continue
                    h_a = [by_label[lab_a][rc]["heavy_rmsd"] for rc in shared]
                    h_b = [by_label[lab_b][rc]["heavy_rmsd"] for rc in shared]
                    wins_a = sum(1 for k in range(len(shared))
                                  if h_a[k] < h_b[k])
                    wins_b = sum(1 for k in range(len(shared))
                                  if h_b[k] < h_a[k])
                    fh.write(f"### {lab_a} vs {lab_b} "
                             f"(n_shared = {len(shared)})\n\n")
                    fh.write(f"- {lab_a} mean: "
                             f"{statistics.fmean(h_a):.3f} Å "
                             f"(wins: {wins_a})\n")
                    fh.write(f"- {lab_b} mean: "
                             f"{statistics.fmean(h_b):.3f} Å "
                             f"(wins: {wins_b})\n")
                    diff = statistics.fmean(h_a) - statistics.fmean(h_b)
                    fh.write(f"- delta: {diff:+.3f} Å "
                             f"({'A better' if diff < 0 else 'B better'})\n\n")

    print(f"[aggregate] wrote {perclass_md}")

    # ---- summary.json ----
    summary: Dict = {
        "date_tag": args.date_tag,
        "builders": {},
        "per_bucket": {},
        "per_chem_class": {},
        "paired_comparisons": {},
        "ccdc_ground_truth": (
            "agent_workspace/quality_framework/reference/ccdc_cleaned_xyz "
            "(997 cleaned refcodes)"
        ),
    }
    for lab in labels:
        recs = list(by_label[lab].values())
        hs = [r["heavy_rmsd"] for r in recs]
        bs = [r["bond_rmsd"] for r in recs]
        ang = [r["angle_rmsd"] for r in recs]
        summary["builders"][lab] = {
            "n_refcodes": len(recs),
            "heavy_rmsd": metric_summary(hs),
            "bond_rmsd": metric_summary(bs),
            "angle_rmsd": metric_summary(ang),
        }
    for bucket in sorted(bucket_data.keys()):
        summary["per_bucket"][bucket] = {}
        for lab in labels:
            vals = bucket_data[bucket].get(lab, [])
            if not vals:
                continue
            summary["per_bucket"][bucket][lab] = {
                "heavy_rmsd": metric_summary([v[0] for v in vals]),
                "bond_rmsd": metric_summary([v[1] for v in vals]),
                "angle_rmsd": metric_summary([v[2] for v in vals]),
            }
    for bucket in sorted(om_bucket_data.keys()):
        summary["per_chem_class"][bucket] = {}
        for lab in labels:
            vals = om_bucket_data[bucket].get(lab, [])
            if not vals:
                continue
            summary["per_chem_class"][bucket][lab] = {
                "heavy_rmsd": metric_summary([v[0] for v in vals]),
                "bond_rmsd": metric_summary([v[1] for v in vals]),
                "angle_rmsd": metric_summary([v[2] for v in vals]),
            }

    # Paired
    uff_label = next((l for l in labels if "uff" in l.lower()), None)
    delfin_labels = [l for l in labels if l != uff_label]
    if uff_label and delfin_labels:
        for dlab in delfin_labels:
            uff_set = set(by_label[uff_label].keys())
            dset = set(by_label[dlab].keys())
            shared = sorted(uff_set & dset)
            if not shared:
                continue
            u_h = [by_label[uff_label][rc]["heavy_rmsd"] for rc in shared]
            d_h = [by_label[dlab][rc]["heavy_rmsd"] for rc in shared]
            u_b = [by_label[uff_label][rc]["bond_rmsd"] for rc in shared]
            d_b = [by_label[dlab][rc]["bond_rmsd"] for rc in shared]
            wins = sum(1 for i in range(len(shared))
                       if d_h[i] < u_h[i])
            summary["paired_comparisons"][f"{dlab}_vs_{uff_label}"] = {
                "n_shared": len(shared),
                "shared_refcodes": shared,
                "uff_heavy_mean": float(statistics.fmean(u_h)),
                "delfin_heavy_mean": float(statistics.fmean(d_h)),
                "uff_bond_mean": float(statistics.fmean(u_b)),
                "delfin_bond_mean": float(statistics.fmean(d_b)),
                "delfin_wins_heavy_rmsd": wins,
                "delfin_pct_improvement_heavy": (
                    (1.0 - statistics.fmean(d_h)
                     / statistics.fmean(u_h)) * 100
                ),
            }

    summary_path = out_dir / "xrd_rmsd_summary.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[aggregate] wrote {summary_path}")


if __name__ == "__main__":
    main()
