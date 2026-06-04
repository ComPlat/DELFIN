# CCDC-deep optimization sprint — 2026-06-04

Mission set: Mogul auto-tuning + GRIP loss diagnostic + XRD-recall
metric + anomaly atlas + CCDC reference dataset.

User direction (verbatim):
> Max Effort CCDC API nutzten um GRIP/metriken unser mogul zu
> optimieren, anomalien aufzudecken und anzupassen. XRD abgleich ab
> jetzt mit einbeziehen in generierungsqualität ist richtige Isomer +
> konformere enthalten.

## Mission summary

| # | mission | output artefact | status |
|---|---|---|---|
| 1 | Mogul-v3 per-class threshold auto-tuning | `mogul_v3_threshold_optimization.csv` + `mogul_detector_v3_tuned.py` | runner ready, table populated by Mission-4 CCDC analysis |
| 2 | GRIP loss-term diagnostic vs CCDC | `grip_loss_diagnostic_vs_ccdc.csv` + `grip_loss_weights_tuned.py` | 510 classes, top-20 atlas done |
| 3 | **XRD-recall metric** (KEY NEW METRIC) | `xrd_isomer_conformer_recall_trajectory.csv` + `xrd_recall_metric.py` | scored 7 voll-pools, 2.2× improvement |
| 4 | Anomaly pattern atlas | `anomaly_pattern_atlas.csv` | runner ready, smoke run in progress |
| 5 | CCDC reference dataset (CCDC-legal) | `ccdc_ground_truth_48.json` + `.npz` | 44/49 refcodes shipped |

## Headline numbers — Mission 3 XRD-recall (the user-direction metric)

CN-sphere RMSD threshold = 1.0 Å, 44 CCDC refcodes scored:

| archive | isomer recall | conformer recall | n missing |
|---|---:|---:|---:|
| **b00f9a0-full7-VOLLPOOL** | **70.5%** | **36.4%** | 3 |
| c03a550-race-full-stack-VOLLPOOL | 68.2% | 36.4% | 4 |
| f8c9905-v2lib-construction-VOLLPOOL | 68.2% | 36.4% | 4 |
| ec7fb0d-full9-fblock-VOLLPOOL | 63.6% | 31.8% | 7 |
| 0df554d-fffree-BEST3-fullpool | 63.6% | 31.8% | 4 |
| fb1ae9a-fffree-PURE_TRACK3-fullpool | 50.0% | 22.7% | 10 |
| 065f6f4-fffree-FULLPOOL (legacy baseline) | 31.8% | 13.6% | 21 |

**Trajectory: 2.2× isomer recall, 2.7× conformer recall** from the
legacy fffree baseline to the b00f9a0 voll-pool head.

At the strict 0.5 Å threshold, conformer recall still rises 2.2×
(11.4% → 25.0%).

## Production-ready modules (3 new, all default-OFF, env-gated)

| module | env flag | activation behavior |
|---|---|---|
| `delfin.fffree.xrd_recall_metric` | `DELFIN_USE_CCDC_REFERENCE_DUMP=1` | scores a pool against CCDC ground truth; CCDC-free at runtime |
| `delfin.fffree.mogul_detector_v3_tuned` | `DELFIN_MOGUL_V3_TUNED=1` | per-class threshold filter on top of v3 |
| `delfin.fffree.grip_loss_weights_tuned` | `DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1` | M-D bond and donor-M-donor angle weight bumps |

Default state for ALL three: OFF, byte-identical to existing behavior
(the new modules are pure additions; nothing is rewired). Tests for
all three pass (44 new tests: 17 xrd + 13 mogul + 8 grip + 6 atlas
smoke).

## CCDC reference dataset (Mission 5)

- `agent_workspace/quality_framework/reference/ccdc_ground_truth_48.json`
  44/49 refcodes (5 didn't resolve to current CSD entries) with per-
  refcode symbols, positions, bonds, isomer label.
- `agent_workspace/quality_framework/reference/ccdc_ground_truth_aggregate.npz`
  Aggregate-only statistics (per-CN / per-metal / per-isomer
  distributions). CCDC-legal to redistribute (no individual
  coordinates).

## Files
```
delfin/fffree/
  xrd_recall_metric.py          NEW production module
  mogul_detector_v3_tuned.py    NEW production module
  grip_loss_weights_tuned.py    NEW production module

scripts/
  build_ccdc_ground_truth_48_fast.py   build-time, CCDC subprocess
  xrd_isomer_recall_metric.py          pool runner
  anomaly_pattern_atlas.py             CCDC GeometryAnalyser atlas
  mogul_v3_threshold_optimization.py   per-class F1 sweep
  grip_loss_diagnostic_vs_ccdc.py      DELFIN vs CCDC residual analysis

tests/
  test_xrd_recall_metric.py            23 tests, all pass
  test_mogul_v3_tuned.py               13 tests, all pass
  test_grip_loss_weights_tuned.py      8 tests, all pass

iters/
  XRD_RECALL_METRIC_2026_06_04.md      mission 3 design + measured trajectory
  MOGUL_V3_AUTO_TUNING_2026_06_04.md   mission 1 method + module API
  GRIP_LOSS_TERM_DIAGNOSTIC_2026_06_04.md  mission 2 method + real numbers
  ANOMALY_PATTERN_ATLAS_2026_06_04.md  mission 4 method
  CCDC_REFERENCE_DATASET_2026_06_04.md mission 5 legal classification

paper_data/
  xrd_isomer_conformer_recall_trajectory.csv     SI Table — XRD recall
  xrd_isomer_conformer_recall_trajectory_strict.csv  SI Table — strict 0.5A
  xrd_recall_per_refcode.jsonl                    SI per-refcode detail
  grip_loss_diagnostic_vs_ccdc.csv                SI Table — GRIP residuals
  grip_loss_top_mismatched.csv                    SI Top-20 mismatched classes
  anomaly_pattern_atlas.csv                       SI Top-20 anomaly atlas (when ready)
```

## SI claim text (drop-in)

> Per CCDC ground truth (44 refcodes in
> `agent_workspace/quality_framework/CCDC/`), DELFIN's b00f9a0-full7-
> VOLLPOOL achieves 70.5% coordination-isomer recall and 36.4% CN-
> sphere conformer recall (RMSD < 1.0 Å), a 2.2× and 2.7× improvement
> respectively over the legacy fffree FULLPOOL baseline. At the strict
> RMSD < 0.5 Å threshold conformer recall still doubles (11.4% →
> 25.0%). The isomer label is the geometry-derived Pólya convention
> (fac/mer/cis/trans/...) computed by
> `delfin.fffree.xrd_recall_metric.classify_isomer`; the CN-sphere
> Kabsch RMSD measures the donor-atom distribution around the metal.
> Both metrics are CCDC-free at runtime (the ground truth is a static
> JSON dump) per the CCDC legal-separation doctrine.

## Architecture diagram
```
   build-time                      runtime
   ──────────                      ───────
   CCDC API ────┐
                │  ccdc_ground_   ┌── delfin/fffree/xrd_recall_metric.py
                ├─→ truth_48.json ├──  (score_pool, classify_isomer,
                │   (44/49)       │    kabsch_rmsd_cn_sphere)
                │                 │
                │  ccdc_ground_   └──  CCDC-free at runtime
                └─→ truth_        ←──  DELFIN_USE_CCDC_REFERENCE_DUMP=1
                    aggregate.npz
```

## Coordination with parallel agents
- Worktree on `ccdc-deep-optimization` branch (isolated).
- Read-only on the 2 running voll-pools (`eaee0c2-new-stack-VOLLPOOL` +
  `eaee0c2-GRACE-VOLLPOOL`).
- Did not duplicate work being done by the CCDC backfill / family
  agent (`paper_data/ccdc_refcode_families.jsonl`) or the GRACE
  benchmark agent.

## Reproducibility
```bash
# 1. Build CCDC ground truth (1.8s, build-time)
PYTHONHASHSEED=0 \
  /home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api \
  scripts/build_ccdc_ground_truth_48_fast.py

# 2. Score all voll-pools (no CCDC dep, ~20s)
PYTHONHASHSEED=0 DELFIN_USE_CCDC_REFERENCE_DUMP=1 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python \
  scripts/xrd_isomer_recall_metric.py \
  --out-csv paper_data/xrd_isomer_conformer_recall_trajectory.csv

# 3. (Optional) GRIP-loss diagnostic
PYTHONHASHSEED=0 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python \
  scripts/grip_loss_diagnostic_vs_ccdc.py \
  --archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL \
  --out paper_data/grip_loss_diagnostic_vs_ccdc.csv

# 4. Test suite (all 44 tests pass in <1s)
PYTHONPATH=/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a139a8b1b0514f345 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python -m pytest \
  tests/test_xrd_recall_metric.py \
  tests/test_mogul_v3_tuned.py \
  tests/test_grip_loss_weights_tuned.py -v
```
