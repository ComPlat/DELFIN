#!/usr/bin/env python3
"""G1 — paper-grade markdown report from JSON summary + pool-wide fallback."""
import argparse
import json
from pathlib import Path


def fmt(x, p=4):
    if x is None or x != x:
        return "n/a"
    return f"{x:.{p}f}"


def winner(uff, fff, heal):
    cands = []
    if uff == uff: cands.append(("UFF", uff))
    if fff == fff: cands.append(("DELFIN-fffree", fff))
    if heal == heal: cands.append(("DELFIN-heal", heal))
    if not cands:
        return "n/a"
    return min(cands, key=lambda kv: kv[1])[0]


def pct(uff, alt):
    if uff != uff or alt != alt or uff == 0:
        return "n/a"
    d = (uff - alt) / uff * 100.0
    sign = "+" if d > 0 else ""
    return f"{sign}{d:.1f}%"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--md", required=True)
    ap.add_argument("--fallback-fff", default=None)
    ap.add_argument("--fallback-heal", default=None)
    args = ap.parse_args()

    with open(args.json) as f:
        s = json.load(f)
    fall_fff = json.load(open(args.fallback_fff)) if args.fallback_fff else None
    fall_heal = json.load(open(args.fallback_heal)) if args.fallback_heal else None

    cfg = s["config"]
    o = s["overall"]
    h = s["honesty_signals"]
    lines = []
    L = lines.append

    L("# G1 — xtb / DFT Validation: DELFIN vs UFF on stratified sample")
    L("")
    L(f"- xtb method: `{cfg['xtb_method']}`")
    L(f"- Files requested: {cfg['total_files_requested']}  |  processed: {cfg['total_files_processed']}  |  frames: {cfg['total_frames']}  |  frames OK: {cfg['total_frames_ok']}")
    L(f"- Frames per file (max): {cfg['frames_per_file']}  |  prefer_differing: {cfg.get('prefer_differing', True)}")
    L(f"- Size filter: {cfg.get('min_atoms', 6)} <= atoms <= {cfg['max_atoms']}")
    L(f"- UFF pool:    `{cfg['uff_dir']}`")
    L(f"- fffree pool: `{cfg['fff_dir']}`")
    L(f"- heal pool:   `{cfg['heal_dir']}`")
    L("")
    L("Reference geometry is the xtb-GFN2 `--opt loose` relaxation starting from")
    L("the UFF baseline frame. Lower RMSD = closer to the DFT-level minimum.")
    L("")

    if fall_fff or fall_heal:
        L("## Whole-pool UFF-fallback fraction")
        L("")
        L("Across the entire common subset of UFF and DELFIN pools:")
        L("")
        L("| Pool | files | frames | identical to UFF | differing | mismatched |")
        L("|---|---|---|---|---|---|")
        if fall_fff:
            ov = fall_fff["overall"]
            L(f"| DELFIN-fffree | {ov['files']} | {ov['total']} | "
              f"{ov['frac_identical']*100:.1f}% | "
              f"{(1 - ov['frac_identical'] - ov['frac_mismatched'])*100:.1f}% | "
              f"{ov['frac_mismatched']*100:.1f}% |")
        if fall_heal:
            ov = fall_heal["overall"]
            L(f"| DELFIN-heal | {ov['files']} | {ov['total']} | "
              f"{ov['frac_identical']*100:.1f}% | "
              f"{(1 - ov['frac_identical'] - ov['frac_mismatched'])*100:.1f}% | "
              f"{ov['frac_mismatched']*100:.1f}% |")
        L("")
        if fall_fff:
            L("Pool-wide fffree fallback per stratum:")
            L("")
            L("| Class | files | frames | %identical to UFF | %differing |")
            L("|---|---|---|---|---|")
            for s_lbl, d in sorted(fall_fff["per_stratum"].items()):
                L(f"| {s_lbl} | {d['files']} | {d['total']} | "
                  f"{d['frac_identical']*100:.1f}% | "
                  f"{(1 - d['frac_identical'] - d['frac_mismatched'])*100:.1f}% |")
            L("")
        L("")

    L("## Honesty signals on the xtb sample")
    L("")
    L(f"- DELFIN-fffree frames present: **{o['n_fff_present']}** / {o['n_frames']}")
    L(f"  - identical to UFF baseline (<1e-3 Å): **{h['frac_fff_identical_to_uff']*100:.1f}%**")
    L(f"  - differing (constructive output): **{o['n_fff_differing']}** frames")
    L(f"- DELFIN-heal frames present: **{o['n_heal_present']}** / {o['n_frames']}")
    L(f"  - identical to UFF baseline: **{h['frac_heal_identical_to_uff']*100:.1f}%**")
    L(f"  - differing: **{o['n_heal_differing']}** frames")
    L("")
    L("Note: the xtb sample is biased toward differing frames (prefer_differing=1),")
    L("so the identical fraction here is lower than the pool-wide one above. The")
    L("'differing frames only' tables below are the correct comparator for")
    L("constructive logic.")
    L("")

    L("## Bond-length RMSD vs xtb (Å, mean) — per stratum, differing frames only")
    L("")
    L("| Class | n_fff_diff | UFF | DELFIN-fffree | n_heal_diff | DELFIN-heal | Best | UFF→fffree | UFF→heal |")
    L("|---|---|---|---|---|---|---|---|---|")
    for s_lbl, d in s["per_stratum"].items():
        w = winner(d["uff_bond_mean"], d["fff_bond_mean_diff"], d["heal_bond_mean_diff"])
        L(f"| {s_lbl} | {d['n_fff_differing']} | {fmt(d['uff_bond_mean'])} | "
          f"{fmt(d['fff_bond_mean_diff'])} | {d['n_heal_differing']} | "
          f"{fmt(d['heal_bond_mean_diff'])} | {w} | "
          f"{pct(d['uff_bond_mean'], d['fff_bond_mean_diff'])} | "
          f"{pct(d['uff_bond_mean'], d['heal_bond_mean_diff'])} |")
    w = winner(o["uff_bond_mean"], o["fff_bond_mean_diff"], o["heal_bond_mean_diff"])
    L(f"| **Total** | **{o['n_fff_differing']}** | **{fmt(o['uff_bond_mean'])}** | "
      f"**{fmt(o['fff_bond_mean_diff'])}** | **{o['n_heal_differing']}** | "
      f"**{fmt(o['heal_bond_mean_diff'])}** | **{w}** | "
      f"**{pct(o['uff_bond_mean'], o['fff_bond_mean_diff'])}** | "
      f"**{pct(o['uff_bond_mean'], o['heal_bond_mean_diff'])}** |")
    L("")

    L("## Bond-length RMSD vs xtb (Å, mean) — per stratum, all frames")
    L("")
    L("| Class | n | UFF | DELFIN-fffree | DELFIN-heal | Best | UFF→fffree | UFF→heal |")
    L("|---|---|---|---|---|---|---|---|")
    for s_lbl, d in s["per_stratum"].items():
        w = winner(d["uff_bond_mean"], d["fff_bond_mean_all"], d["heal_bond_mean_all"])
        L(f"| {s_lbl} | {d['n_frames']} | {fmt(d['uff_bond_mean'])} | "
          f"{fmt(d['fff_bond_mean_all'])} | {fmt(d['heal_bond_mean_all'])} | {w} | "
          f"{pct(d['uff_bond_mean'], d['fff_bond_mean_all'])} | "
          f"{pct(d['uff_bond_mean'], d['heal_bond_mean_all'])} |")
    w = winner(o["uff_bond_mean"], o["fff_bond_mean_all"], o["heal_bond_mean_all"])
    L(f"| **Total** | **{o['n_frames']}** | **{fmt(o['uff_bond_mean'])}** | "
      f"**{fmt(o['fff_bond_mean_all'])}** | **{fmt(o['heal_bond_mean_all'])}** | "
      f"**{w}** | **{pct(o['uff_bond_mean'], o['fff_bond_mean_all'])}** | "
      f"**{pct(o['uff_bond_mean'], o['heal_bond_mean_all'])}** |")
    L("")

    L("## Bond-angle RMSD vs xtb (deg, mean) — per stratum, differing frames only")
    L("")
    L("| Class | n_fff_diff | UFF | DELFIN-fffree | DELFIN-heal | Best |")
    L("|---|---|---|---|---|---|")
    for s_lbl, d in s["per_stratum"].items():
        w = winner(d["uff_angle_mean"], d["fff_angle_mean_diff"], d["heal_angle_mean_diff"])
        L(f"| {s_lbl} | {d['n_fff_differing']} | {fmt(d['uff_angle_mean'])} | "
          f"{fmt(d['fff_angle_mean_diff'])} | {fmt(d['heal_angle_mean_diff'])} | {w} |")
    w = winner(o["uff_angle_mean"], o["fff_angle_mean_diff"], o["heal_angle_mean_diff"])
    L(f"| **Total** | **{o['n_fff_differing']}** | **{fmt(o['uff_angle_mean'])}** | "
      f"**{fmt(o['fff_angle_mean_diff'])}** | **{fmt(o['heal_angle_mean_diff'])}** | **{w}** |")
    L("")

    L("## Heavy-atom Kabsch RMSD vs xtb (Å, mean) — per stratum, differing frames only")
    L("")
    L("| Class | n_fff_diff | UFF | DELFIN-fffree | DELFIN-heal | Best |")
    L("|---|---|---|---|---|---|")
    for s_lbl, d in s["per_stratum"].items():
        w = winner(d["uff_kabsch_mean"], d["fff_kabsch_mean_diff"], d["heal_kabsch_mean_diff"])
        L(f"| {s_lbl} | {d['n_fff_differing']} | {fmt(d['uff_kabsch_mean'])} | "
          f"{fmt(d['fff_kabsch_mean_diff'])} | {fmt(d['heal_kabsch_mean_diff'])} | {w} |")
    w = winner(o["uff_kabsch_mean"], o["fff_kabsch_mean_diff"], o["heal_kabsch_mean_diff"])
    L(f"| **Total** | **{o['n_fff_differing']}** | **{fmt(o['uff_kabsch_mean'])}** | "
      f"**{fmt(o['fff_kabsch_mean_diff'])}** | **{fmt(o['heal_kabsch_mean_diff'])}** | **{w}** |")
    L("")

    L("## Honest disclosure — classes where UFF beats DELFIN (differing frames)")
    L("")
    losers = []
    for s_lbl, d in s["per_stratum"].items():
        u = d["uff_bond_mean"]; f = d["fff_bond_mean_diff"]; h_ = d["heal_bond_mean_diff"]
        if u != u:
            continue
        if d["n_fff_differing"] == 0 and d["n_heal_differing"] == 0:
            continue
        beats_fff = (f != f) or (u < f)
        beats_heal = (h_ != h_) or (u < h_)
        if beats_fff and beats_heal:
            losers.append((s_lbl, u, f, h_, d["n_fff_differing"], d["n_heal_differing"]))
    if losers:
        for s_lbl, u, f_, h_, nf, nh in losers:
            L(f"- **{s_lbl}** UFF={fmt(u)}  beats  fffree={fmt(f_)} (n={nf}), heal={fmt(h_)} (n={nh})")
    else:
        L("- (none — DELFIN ties or wins every class where it produced differing output)")
    L("")

    L("## Paper claim (suggested)")
    L("")

    def paper_claim(method_label, key_diff, n_diff_key, kabsch_diff_key):
        u_b = o["uff_bond_mean"]
        f_b = o[key_diff]
        u_a = o["uff_angle_mean"]
        f_a = o[key_diff.replace("bond", "angle")]
        u_k = o["uff_kabsch_mean"]
        f_k = o[kabsch_diff_key]
        n = o[n_diff_key]
        L(f"### {method_label} (differing frames, n={n})")
        if f_b == f_b and u_b > 0:
            dpct = (u_b - f_b) / u_b * 100
            adj = "improvement toward" if dpct > 0 else "regression from"
            L(f"- bond RMSD: **{fmt(f_b)}** vs UFF **{fmt(u_b)}**  →  {dpct:+.1f}% {adj} DFT")
        if f_a == f_a and u_a > 0:
            dpct = (u_a - f_a) / u_a * 100
            adj = "improvement toward" if dpct > 0 else "regression from"
            L(f"- angle RMSD: **{fmt(f_a)}°** vs UFF **{fmt(u_a)}°**  →  {dpct:+.1f}% {adj} DFT")
        if f_k == f_k and u_k > 0:
            dpct = (u_k - f_k) / u_k * 100
            adj = "improvement toward" if dpct > 0 else "regression from"
            L(f"- heavy-atom Kabsch RMSD: **{fmt(f_k)}** vs UFF **{fmt(u_k)}**  →  {dpct:+.1f}% {adj} DFT")
        L("")

    paper_claim("DELFIN-fffree", "fff_bond_mean_diff",
                "n_fff_differing", "fff_kabsch_mean_diff")
    paper_claim("DELFIN-heal", "heal_bond_mean_diff",
                "n_heal_differing", "heal_kabsch_mean_diff")

    L("## Reproducibility")
    L("")
    L("```")
    L("bash scripts/xtb_validation_pipeline.sh")
    L("```")
    L("")
    L("Selection deterministic: sorted filenames, every-Nth file inside each")
    L("stratum (sqrt(class_count) weighting). Within each file, frame indices")
    L("where DELFIN differs from UFF are prioritized.")

    Path(args.md).write_text("\n".join(lines) + "\n")
    print(f"[G1] markdown report -> {args.md}")


if __name__ == "__main__":
    main()
