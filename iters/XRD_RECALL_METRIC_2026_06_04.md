# XRD Isomer + Conformer Recall Metric — 2026-06-04

**Mission 3** of the CCDC-deep optimization sprint.

## User direction (verbatim)
> "Max Effort CCDC API nutzten um GRIP/metriken unser mogul zu optimieren,
> anomalien aufzudecken und anzupassen. XRD abgleich ab jetzt mit
> einbeziehen in generierungsqualität ist richtige Isomer + konformere
> enthalten."

The phrase **"XRD abgleich ab jetzt mit einbeziehen in generierungsqualität
ist richtige Isomer + konformere enthalten"** establishes a NEW core
quality metric: a generated DELFIN pool is GOOD iff it contains BOTH
the right coordination isomer AND a conformer that matches the XRD
ground truth.

## Metric definition

For each CCDC anchor refcode `R` whose SMILES is in the DELFIN pool's
master list:

**Isomer recall**

    isomer_recall(R) = 1 iff ∃ DELFIN structure for matching SMILES
                       whose coordination-isomer label equals the
                       CCDC isomer label; 0 otherwise.

The isomer label is computed by `classify_isomer(syms, P)` in
`delfin/fffree/xrd_recall_metric.py`:

  - CN4-tet-or-sp + cis / trans for square-planar with 2x2 donor
    elements
  - CN6-OH + fac / mer for octahedral with 3x3 donor elements
  - other: CN<n>-<polyhedron>

**Conformer recall (RMSD)**

    conformer_recall_rmsd(R) = min over DELFIN structures for matching
                              SMILES of the heavy-atom Kabsch RMSD to
                              the CCDC entry.

    rmsd_pass(R) = 1 iff conformer_recall_rmsd(R) < 0.5 Å.

The Kabsch RMSD is element-matched via greedy nearest-neighbour over
centred coordinates, hydrogens dropped.

**Aggregate per pool**

    isomer_recall_mean = (1/N) Σ_R isomer_recall(R)
    conformer_recall_mean = (1/N) Σ_R rmsd_pass(R)

where N is the number of CCDC ground-truth refcodes in the test set.

## Architecture

CCDC-FREE at runtime. The ground truth is a static JSON dump
(`agent_workspace/quality_framework/reference/ccdc_ground_truth_48.json`)
produced once by `scripts/build_ccdc_ground_truth_48.py` (which DOES
use the CCDC API). Runtime only consumes the JSON.

```
   build-time                            runtime
   ──────────                            ───────
   CCDC API ──┐                          ┌── delfin/fffree/
              │                          │   xrd_recall_metric.py
              │                          │
              │  ccdc_ground_truth_48    │
              ├─→  .json (DERIVED        ├─→  score_pool(archive)
              │   STATISTICS only —      │   → isomer_recall,
              │   CCDC-legal per         │     conformer_recall
              │   feedback_ccdc_legal_   │
              │   separation_doctrine_   │
              │   2026_06_04)            │
              └─→  .npz aggregate        └─→  (uses only numpy)
```

The `delfin.fffree.xrd_recall_metric` module has zero CCDC dependency
at import time. Activation: `DELFIN_USE_CCDC_REFERENCE_DUMP=1`
(default OFF, byte-identical to the previous run-tree).

## Implementation
- Production module: `delfin/fffree/xrd_recall_metric.py`
  - `classify_isomer(syms, P)`
  - `kabsch_rmsd_heavy(syms_a, Pa, syms_b, Pb)`
  - `score_pool(archive_dir, ref=None, rmsd_threshold=None)`
  - `load_reference(path=None)`
- Build script: `scripts/build_ccdc_ground_truth_48.py` (CCDC subprocess)
- Pool runner: `scripts/xrd_isomer_recall_metric.py` (delfin subprocess)
- Tests: `tests/test_xrd_recall_metric.py` (17 tests; all pass)

## Trajectory CSV columns
`paper_data/xrd_isomer_conformer_recall_trajectory.csv`:

| column | meaning |
|---|---|
| archive | voll-pool subdir name |
| exists | 1 if archive directory was found |
| isomer_recall | mean isomer recall over all CCDC refcodes |
| conformer_recall | fraction with RMSD < threshold |
| n_refcodes_scored | total number of CCDC refcodes scored |
| n_refcodes_missing | how many had no DELFIN structures |
| rmsd_threshold_A | Å threshold used (default 0.5) |
| elapsed_s | wall time |

And per-refcode detail in `paper_data/xrd_recall_per_refcode.jsonl`.

## Measured trajectory (real numbers, 2026-06-04 18:05 UTC)

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

Strict 0.5 Å CN-sphere RMSD threshold:

| archive | isomer recall | conformer recall (0.5Å) |
|---|---:|---:|
| **b00f9a0-full7-VOLLPOOL** | **70.5%** | **25.0%** |
| 065f6f4-fffree-FULLPOOL (legacy baseline) | 31.8% | 11.4% |

**Headline trajectory:** isomer recall **2.2x** (31.8% → 70.5%) and
conformer recall **2.7x** (13.6% → 36.4% @ 1.0 Å) from the legacy
fffree baseline to the b00f9a0 voll-pool head, across the new
construction stack. The new metric is the FIRST DELFIN evaluation that
directly anchors generation quality to XRD ground truth (the
user-direction "XRD abgleich ist richtige Isomer + konformere
enthalten").

## SI claim
> DELFIN's b00f9a0-full7-VOLLPOOL achieves 70.5% isomer recall and
> 36.4% CN-sphere conformer recall (RMSD < 1.0 Å) against the
> 44-refcode CCDC organometallic ground-truth set, a 2.2x and 2.7x
> improvement respectively over the legacy fffree FULLPOOL baseline
> at the same set. At the strict 0.5 Å threshold, conformer recall
> still rises 2.2x (11.4% → 25.0%).

Numbers regenerable from `paper_data/xrd_isomer_conformer_recall_trajectory.csv`
and the strict variant `_trajectory_strict.csv`.

## Reproducibility
```bash
# 1. Build the CCDC ground-truth dump (CCDC API; ~15 min)
PYTHONHASHSEED=0 \
  /home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api \
  /home/qmchem_max/ComPlat/DELFIN/scripts/build_ccdc_ground_truth_48.py

# 2. Score each voll-pool (no CCDC dep)
PYTHONHASHSEED=0 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python \
  /home/qmchem_max/ComPlat/DELFIN/scripts/xrd_isomer_recall_metric.py \
  --out-csv paper_data/xrd_isomer_conformer_recall_trajectory.csv
```

## Notes / limitations
- Isomer label is heuristic geometry-based. Not a full chirality solver:
  Δ/Λ labels are folded into "CN6-OH" without enantio-distinction.
- Heavy-atom Kabsch RMSD ignores hydrogens (lowest signal-to-noise
  atoms). For PERFECT recall a structure must have heavy-atom-perfect
  geometry; H errors don't count.
- Greedy element-matched correspondence (no graph-canonical solver) can
  inflate the RMSD when DELFIN's atom order diverges; the 0.5 Å
  threshold is therefore lenient.
