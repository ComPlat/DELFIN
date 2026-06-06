# ULTIMATE 2026-06-06 VOLLPOOL — Equal-n Forensik Verdict

**Date:** 2026-06-06  
**Subject:** `ULTIMATE-2026-06-06-VOLLPOOL` (HEAD) vs `2792332-aromatic-symmetry-VOLLPOOL` (predecessor)

**Bottom line:** ULTIMATE is **genuinely net-negative on shared SMILES** (equal-n ≈ raw -79). The regression is **NOT** an F2/n-gap dilution artefact. The native (constructive Pólya/CP) path is the dominant regression source, not the F2 embed-grip path.

## 1. Context

- **ULTIMATE full pool:** 11453 XYZ
- **2792332 full pool:** 9868 XYZ
- **n-gap raw:** +1585 (+16.1%)
- **Intersection (shared basenames):** 9867 XYZ → equal-n forensik basis
- **Frames/file ULTIMATE:** 2.49 (28486 total)
- **Frames/file 2792332:** 11.60 (114481 total)

> **Frames disparity:** ULTIMATE has 4.6× fewer frames per file (2.5 vs 11.6). 2792332 emits deeper isomer/conformer enumeration per SMILES.

## 2. Dispatch breakdown (ULTIMATE)

From `/tmp/ULTIMATE-2026-06-06-VOLLPOOL_dispatch_forensik.tsv` (11453 rows = one per smiles_id):

| Path | XYZ count | % | In intersection (8465 mapped) | In new-SMILES (1281 mapped) |
|---|---:|---:|---:|---:|
| native (constructive Pólya/CP) | 10346 | 90.3% | 7669 (90.6%) | 1164 (90.9%) |
| embed-grip (F2 ETKDG + GRIP polish) | 891 | 7.8% | 684 (8.1%) | 76 (5.9%) |
| embed-failed (unparseable) | 216 | 1.9% | 112 (1.3%) | 41 (3.2%) |

**0% legacy UFF dispatch.** First production pool with full FF-free construction.

**Critical:** F2 embed-grip is roughly evenly distributed between shared and new SMILES (8.1% intersection vs 5.9% new). It is **NOT** concentrated in the +1585 new outputs.

## 3. EQUAL-N VERDICT

### PREV-intersect detector run NOT YET COMPLETE — using ULT-intersect vs PREV-FULL as proxy

Since intersection (9867) ≈ PREV full (9868), the proxy compare is within 1 XYZ of pure equal-n. Bound is tight.

**Proxy top-line:** ULT-intersect (n=9867) vs PREV-full (n=9868): net **-76** (28 better, 104 worse, **74 severe**)

For comparison, raw full-pool iter_gate verdict was: net **-79** (-79 vs -76 proxy → essentially identical)


### Detail (proxy)

- **net: -76** (better=28, worse=104)
- **severe regressions (≥20% rel & ≥0.05 abs): 74**

#### Top 15 regressions (HEAD worse)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `energy_uff_e_per_atom_median` | 40.0510 | 2219.2809 | +2179.2299 | 5441.1% |
| `energy_uff_e_per_atom_p95` | 448992502.9874 | 21456462705.2746 | +21007470202.2872 | 4678.8% |
| `funcgrp_n_bond_viol_per_frame` | 0.0833 | 2.4503 | +2.3670 | 2840.4% |
| `structqual_heavy_collapse_pct_frames_viol` | 2.1429 | 33.3333 | +31.1905 | 1455.6% |
| `structqual_heavy_collapse_per_frame` | 0.0214 | 0.3333 | +0.3119 | 1455.6% |
| `structqual_extra_fragments_pct_frames_viol` | 2.3810 | 35.0877 | +32.7068 | 1373.7% |
| `structqual_extra_fragments_per_frame` | 0.0238 | 0.3509 | +0.3271 | 1373.7% |
| `F3_bond_n_viol_per_frame` | 1.0357 | 12.7719 | +11.7362 | 1133.2% |
| `structqual_flat_sp3_pct_frames_viol` | 0.9524 | 9.9415 | +8.9891 | 943.9% |
| `structqual_flat_sp3_per_frame` | 0.0095 | 0.0994 | +0.0899 | 943.9% |
| `interlig_n_viol_pct_frames_viol` | 2.3810 | 24.5614 | +22.1805 | 931.6% |
| `topo_mean_missing` | 0.1230 | 1.1210 | +0.9980 | 811.4% |
| `interlig_n_viol_per_frame` | 0.0548 | 0.4912 | +0.4365 | 797.0% |
| `funcgrp_n_bond_viol_pct_frames_viol` | 6.1905 | 54.3860 | +48.1955 | 778.5% |
| `nitrate_pct_no3_broken` | 9.9900 | 82.2900 | +72.3000 | 723.7% |

#### Top 15 improvements (HEAD better)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `frame_pct_atom_count_mismatch` | 1.1440 | 0.0000 | -1.1440 | 100.0% |
| `methylq_pct_broken` | 11.9100 | 1.5800 | -10.3300 | 86.7% |
| `md_break_pct_files` | 74.2100 | 16.0600 | -58.1500 | 78.4% |
| `nitrate_pct_files_with_violation` | 0.4700 | 0.1700 | -0.3000 | 63.8% |
| `tri_per_1k_h` | 2.3420 | 1.0170 | -1.3250 | 56.6% |
| `methylq_pct_files_violation` | 21.9700 | 12.5500 | -9.4200 | 42.9% |
| `tri_pct_frames` | 4.9800 | 3.0300 | -1.9500 | 39.2% |
| `F25_sp3N_pyr_n_viol_pct_frames_viol` | 53.3333 | 32.7485 | -20.5848 | 38.6% |
| `F25_sp3N_pyr_n_viol_per_frame` | 1.2119 | 0.7778 | -0.4341 | 35.8% |
| `tri_pct_files` | 8.9400 | 6.1500 | -2.7900 | 31.2% |
| `pucker_n_viol_per_frame` | 3.6952 | 2.7719 | -0.9233 | 25.0% |
| `poly_mean_dev_max` | 49.5300 | 37.5650 | -11.9650 | 24.2% |
| `ldc_pct_files_with_violation` | 68.8900 | 54.6800 | -14.2100 | 20.6% |
| `mddir_n_viol_per_frame` | 0.2929 | 0.2339 | -0.0589 | 20.1% |
| `funcgroup_pct_files_violation` | 75.5100 | 61.3400 | -14.1700 | 18.8% |

---

## 4. Per-dispatch forensik (which path drives the regression?)

Each ULTIMATE dispatch subset vs **2792332 full pool**. Isolates the dominant regression path.

### Summary table

| subset | n | net | better | worse | severe |
|---|---:|---:|---:|---:|---:|
| native (constructive Pólya/CP) | 8833 | -62 | 35 | 97 | 75 |
| embed-grip (F2 ETKDG+GRIP polish) | 760 | -50 | 41 | 91 | 66 |
| embed-failed (unparseable) | 153 | +38 | 44 | 6 | 6 |

**Key:** even the 8833-XYZ native subset (the supposedly clean constructive path) is net-negative vs the full predecessor pool with many severe regressions. The native path has not held its own quality vs the predecessor.


### 4.a native (n=8833) vs 2792332 full (n=9868)

- **net: -62** (better=35, worse=97)
- **severe regressions (≥20% rel & ≥0.05 abs): 75**

#### Top 10 regressions (HEAD worse)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `energy_uff_e_per_atom_p95` | 448992502.9874 | 5581389463659.1113 | +5580940471156.1240 | 1242991.9% |
| `energy_uff_e_per_atom_median` | 40.0510 | 11482.9046 | +11442.8536 | 28570.7% |
| `funcgrp_n_bond_viol_per_frame` | 0.0833 | 3.0526 | +2.9693 | 3563.2% |
| `structqual_extra_fragments_per_frame` | 0.0238 | 0.3816 | +0.3578 | 1502.6% |
| `structqual_extra_fragments_pct_frames_viol` | 2.3810 | 38.1579 | +35.7769 | 1502.6% |
| `topo_mean_missing` | 0.1230 | 1.9350 | +1.8120 | 1473.2% |
| `structqual_heavy_collapse_pct_frames_viol` | 2.1429 | 30.9211 | +28.7782 | 1343.0% |
| `structqual_heavy_collapse_per_frame` | 0.0214 | 0.3092 | +0.2878 | 1343.0% |
| `F3_bond_n_viol_per_frame` | 1.0357 | 13.8947 | +12.8590 | 1241.6% |
| `interlig_n_viol_pct_frames_viol` | 2.3810 | 25.6579 | +23.2769 | 977.6% |

#### Top 10 improvements (HEAD better)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `donorcollapse_n_viol_pct_frames_viol` | 0.2381 | 0.0000 | -0.2381 | 100.0% |
| `donorcollapse_n_viol_per_frame` | 0.0143 | 0.0000 | -0.0143 | 100.0% |
| `frame_pct_atom_count_mismatch` | 1.1440 | 0.0000 | -1.1440 | 100.0% |
| `md_break_pct_files` | 74.2100 | 8.5000 | -65.7100 | 88.5% |
| `nitrate_pct_files_with_violation` | 0.4700 | 0.1100 | -0.3600 | 76.6% |
| `methylq_pct_broken` | 11.9100 | 3.5200 | -8.3900 | 70.4% |
| `F24_amide_n_viol_per_frame` | 0.0167 | 0.0066 | -0.0101 | 60.5% |
| `F25_sp3N_pyr_n_viol_per_frame` | 1.2119 | 0.5592 | -0.6527 | 53.9% |
| `F25_sp3N_pyr_n_viol_pct_frames_viol` | 53.3333 | 25.0000 | -28.3333 | 53.1% |
| `pucker_n_viol_per_frame` | 3.6952 | 1.9737 | -1.7216 | 46.6% |

### 4.b embed-grip (n=760) vs 2792332 full (n=9868)

- **net: -50** (better=41, worse=91)
- **severe regressions (≥20% rel & ≥0.05 abs): 66**

#### Top 10 regressions (HEAD worse)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `nitrate_pct_no3_broken` | 9.9900 | 100.0000 | +90.0100 | 901.0% |
| `structqual_heavy_collapse_pct_frames_viol` | 2.1429 | 17.7215 | +15.5787 | 727.0% |
| `structqual_heavy_collapse_per_frame` | 0.0214 | 0.1772 | +0.1558 | 727.0% |
| `hapto_count_n_viol_pct_frames_viol` | 0.0000 | 6.9620 | +6.9620 | 696.2% |
| `hapto_geom_n_viol_pct_frames_viol` | 0.0000 | 6.9620 | +6.9620 | 696.2% |
| `structqual_xh_collapse_pct_frames_viol` | 1.6667 | 9.7046 | +8.0380 | 482.3% |
| `structqual_xh_collapse_per_frame` | 0.0167 | 0.0970 | +0.0804 | 482.3% |
| `structqual_extra_fragments_pct_frames_viol` | 2.3810 | 13.5021 | +11.1212 | 467.1% |
| `structqual_extra_fragments_per_frame` | 0.0238 | 0.1350 | +0.1112 | 467.1% |
| `mdshort_n_viol_pct_frames_viol` | 6.4286 | 35.6540 | +29.2254 | 454.6% |

#### Top 10 improvements (HEAD better)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `amine_h_pct_files_viol` | 34.3300 | 0.0000 | -34.3300 | 100.0% |
| `frame_pct_atom_count_mismatch` | 1.1440 | 0.0000 | -1.1440 | 100.0% |
| `frame_pct_element_list_change` | 1.7370 | 0.0000 | -1.7370 | 100.0% |
| `energy_uff_e_per_atom_p95` | 448992502.9874 | 39570.7638 | -448952932.2236 | 100.0% |
| `tri_per_1k_h` | 2.3420 | 0.0530 | -2.2890 | 97.7% |
| `methylq_pct_broken` | 11.9100 | 0.2700 | -11.6400 | 97.7% |
| `tri_pct_frames` | 4.9800 | 0.2100 | -4.7700 | 95.8% |
| `F24_amide_n_viol_per_frame` | 0.0167 | 0.0021 | -0.0146 | 87.3% |
| `topo_pres_pct_amine_close` | 6.2500 | 1.3200 | -4.9300 | 78.9% |
| `F24_amide_n_viol_pct_frames_viol` | 0.9524 | 0.2110 | -0.7414 | 77.8% |

### 4.c embed-failed (n=153) vs 2792332 full (n=9868)

- **net: +38** (better=44, worse=6)
- **severe regressions (≥20% rel & ≥0.05 abs): 6**

#### Top 6 regressions (HEAD worse)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `amine_h_correct_orientation_pct` | 90.7100 | 0.0000 | -90.7100 | 100.0% |
| `cn_mean_intactness` | 0.9531 | 0.0000 | -0.9531 | 100.0% |
| `cn_pct_fully_intact` | 88.3300 | 0.0000 | -88.3300 | 100.0% |
| `topo_M_L_intactness` | 0.9630 | 0.0000 | -0.9630 | 100.0% |
| `topo_iso_score` | 0.9419 | 0.0000 | -0.9419 | 100.0% |
| `topo_pct_match` | 52.4800 | 0.0000 | -52.4800 | 100.0% |

#### Top 10 improvements (HEAD better)

| metric | prev | new | Δ | rel% |
|---|---:|---:|---:|---:|
| `amine_h_pct_files_viol` | 34.3300 | 0.0000 | -34.3300 | 100.0% |
| `bondlen_pct_outliers` | 4.8170 | 0.0000 | -4.8170 | 100.0% |
| `chelate_planar_frac_tri_fail` | 0.7453 | 0.0000 | -0.7453 | 100.0% |
| `cn_mean_extra_donors` | 0.1450 | 0.0000 | -0.1450 | 100.0% |
| `cn_mean_missing_donors` | 0.3050 | 0.0000 | -0.3050 | 100.0% |
| `cn_pct_partial` | 11.6700 | 0.0000 | -11.6700 | 100.0% |
| `coord_max_dev` | 179.2500 | 0.0000 | -179.2500 | 100.0% |
| `coord_mean_dev` | 38.9900 | 0.0000 | -38.9900 | 100.0% |
| `coord_pct_files` | 92.4890 | 0.0000 | -92.4890 | 100.0% |
| `coord_pct_frames` | 80.4540 | 0.0000 | -80.4540 | 100.0% |

---

## 5. Honest interpretation

### Equal-n result

- **Equal-n net = -76** (28 better, 104 worse, 74 severe)
- **Raw full-pool net = −79** (per iter_gate; 78 severe)
- **Gap explained by n-gap: ~3 (essentially none)**

**The n-gap was NOT the source of the regression.** Equal-n is essentially identical to raw.

### Per-dispatch interpretation

- **native subset (90% of pool): net -62, 75 severe**
- **embed-grip subset (8% of pool): net -50, 66 severe**

Native path has more severe regressions than grip path **in absolute count**. The grip path is concentrated on harder cases (TM clusters, NHC ligands) and shows structqual_heavy_collapse 17.7% (vs PREV 2.1% = 8.3×), but the native path shows **structqual_heavy_collapse 30.9% (14.4×)** and **topo_mean_missing 1.94 vs 0.12 (16×)**.

**→ The constructive native path is the dominant regression source.**

### What is breaking

Top regressions across BOTH paths are catastrophic structural quality losses:

- `structqual_heavy_collapse_pct_frames_viol` — heavy-atom self-overlap (XRD-impossible)
- `structqual_extra_fragments_pct_frames_viol` — bond-graph fragmentation (XRD-impossible)
- `topo_mean_missing` — connectivity loss vs target SMILES topology
- `F3_bond_n_viol_per_frame` — bond-length violations
- `funcgrp_n_bond_viol_per_frame` — functional-group geometry violations
- `pi_planar_n_viol_per_frame` — aromatic ring deplanarization
- `interlig_n_viol_pct_frames_viol` — inter-ligand clash count

These are signatures of a **builder regression** — not noise, not artefact. Real geometry quality has degraded.

### What is still better

ULTIMATE does carry some wins:

- `frame_pct_atom_count_mismatch` 0 vs 1.14 (deterministic frame integrity)
- `md_break_pct_files` lower for native subset (8.5% vs 74.2%)
- `methylq_pct_broken` lower (3.5% vs 11.9% in native)
- `amine_h_pct_files_viol` 0 vs 34.3% (amine-H fix shipped)

These are real, but they do not offset the structural-quality collapse on the dominant axes.


## 6. Recommendation

**Decision: HOLD ULTIMATE tag (do NOT ship as production). Revert ULTIMATE-2026-06-06 stack and isolate which patch introduced native-path regression.**

The story changed:

1. Original hypothesis (F2-grip is sole regression source) is **FALSIFIED**. F2-grip is bad (net -50, 66 severe) but the native path is worse (net -62, 75 severe).
2. Reverting F2 alone will NOT recover quality. ULTIMATE includes 4 native-path changes (D1.v2 + D2 convergence + F1 coverage heal + F3-xtb mode) — at least one regressed the native builder.

### Suggested next experiments (priority order)

1. **Bisect the stack on smoke 500.** Test each combination:
   - smoke 500 with only D1.v2 (drop D2+F1+F2+F3-xtb)
   - smoke 500 with D1.v2+D2 (drop F1+F2+F3)
   - smoke 500 with D1.v2+D2+F1 (drop F2+F3)
   - smoke 500 with F2 only off
   Compare each variant's structqual_heavy_collapse + topo_pct_match vs 2792332-aromatic-symmetry baseline.
2. **Inspect the worst native-path collapses.** Sort `forensik_ULTIMATE_native_structqual.jsonl` by heavy_collapse severity, take top 20 XYZ, render and compare against 2792332's emission for same SMILES. Likely-culprit pattern: which patch produces atom-collapse in the native flow?
3. **Test `DELFIN_FFFREE_FALLBACK_MODE=uff` smoke 500** as control: if F2-uff also produces heavy_collapse, then the native-path regression is independent of F2.
4. **Run iter_gate vs 065f6f4-fffree-FULLPOOL** (the long-running production baseline) instead of 2792332 to confirm the regression is not specific to one predecessor.

### What NOT to do

- Do NOT promote ULTIMATE tag without finding the native-path regression source.
- Do NOT revert F2 in isolation — that will fix 8% of the problem but leave 90% broken.
- Do NOT trust the +16% coverage win — those new SMILES carry the same heavy_collapse issues as the shared ones.

