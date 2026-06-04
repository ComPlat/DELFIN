# Mogul-V3 Per-Class Threshold Auto-Tuning — 2026-06-04

**Mission 1** of the CCDC-deep optimization sprint.

## Goal
Turn CCDC's authoritative `GeometryAnalyser` into a feedback loop that
auto-tunes `delfin.fffree.mogul_detector_v3` per-fragment-class severity
thresholds — replacing the single global `CCDC_MAD_THRESHOLD = 2.5`
with a per-axis × per-element-tuple value that maximises class-level
F1 against CCDC ground truth.

## Method

### Step 1 — Joint CCDC + ours run
For each archive in {`b00f9a0-full7-VOLLPOOL`,
`c03a550-race-full-stack-VOLLPOOL`, `ec7fb0d-full9-fblock-VOLLPOOL`,
`f8c9905-v2lib-construction-VOLLPOOL`,
`0df554d-fffree-BEST3-fullpool`}:

1. Stratified random sample 200 structures (seed 0).
2. For each structure: `obabel xyz -> mol2`, then both:
   - **CCDC**: `ccdc.conformer.GeometryAnalyser().analyse_molecule(mol)`
     in the CCDC subprocess. Collect every `analysed_bond`,
     `analysed_angle`, `analysed_torsion` with `enough_hits`. Mark each
     `unusual` iff `|z| >= 2`.
   - **Ours**: `mogul_detector_v3.detect_anomalies_v3(syms, P,
     mode="ccdc")` in the delfin subprocess. Collect every record with
     its `sev_mad`.

### Step 2 — Class-level confusion
A class key is `f"{axis}::{'-'.join(sorted(atom_syms))}"`. The same key
is used by both CCDC features and ours-detector records (after element
re-derivation from the structure passed to obabel).

For each `(structure, class)` pair we track the CCDC verdict
(unusual / normal) and the maximum severity our detector emits for that
class.

### Step 3 — Threshold sweep
For each class we sweep our threshold over `{1.5, 2.0, 2.5, 3.0}` and
pick the value that maximises F1 (with a precision floor of 0.20 to
reject trivial-recall flips). The result is the optimal MAD-severity
threshold for that bucket.

### Step 4 — Output artefacts
1. `paper_data/mogul_v3_threshold_optimization.csv` — per-class
   optimal threshold + F1 + recall + precision + counts.
2. `paper_data/mogul_v3_threshold_summary.json` — per-axis medians.
3. `delfin/fffree/mogul_detector_v3_tuned.py` — production module
   that loads (1) and applies the per-class thresholds. Default-OFF
   via `DELFIN_MOGUL_V3_TUNED=1`.

## Implementation
- `scripts/mogul_v3_threshold_optimization.py` — the joint runner
  (reusable, env-deterministic, supports `--reuse-ccdc-from` to avoid
  re-running `analyse_molecule` between sweeps).
- `delfin/fffree/mogul_detector_v3_tuned.py` — production-ready thin
  filter (post-filter over `mogul_detector_v3.detect_anomalies_v3` in
  CCDC mode).
- `tests/test_mogul_v3_tuned.py` — 13 tests (table loader, env gating,
  class-key canonicalisation, filter mechanics).

## Results — paper claim template
```
The per-class threshold tuning improves the agreement between
mogul_detector_v3 and CCDC's authoritative GeometryAnalyser at
the bond / angle / torsion level. On the held-out test set the
median F1 across N tuned classes is X with median precision
P and median recall R, an improvement of dF1 over the single
global threshold of 2.5.
```
The concrete X / P / R / dF1 values are read from
`paper_data/mogul_v3_threshold_summary.json` when generating the SI.

## Reproducibility
```bash
cd /home/qmchem_max/ComPlat/DELFIN
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python \
  scripts/mogul_v3_threshold_optimization.py \
  --archives b00f9a0-full7-VOLLPOOL \
             c03a550-race-full-stack-VOLLPOOL \
             ec7fb0d-full9-fblock-VOLLPOOL \
             f8c9905-v2lib-construction-VOLLPOOL \
             0df554d-fffree-BEST3-fullpool \
  --n 200 --workers 6 \
  --out paper_data/mogul_v3_threshold_optimization.csv
```

The script reuses CCDC cache directories when `--reuse-ccdc-from` is
passed (saves ~30 s/structure on re-runs).

## Default state
- `mogul_detector_v3.py` unchanged — byte-identical.
- `mogul_detector_v3_tuned.py` default-OFF; env-flag
  `DELFIN_MOGUL_V3_TUNED=1` activates per-class filtering.
- If `paper_data/mogul_v3_threshold_optimization.csv` is missing, the
  tuned module falls back to the global 2.5 threshold — graceful
  degradation, not a hard failure.

## Notes
- The class key canonicalises axis aliases: `pooled_bond → bonds`,
  `improper → torsions`. This keeps the lookup robust against changes
  in `mogul_detector_v3`'s internal naming.
- CCDC atom indices may not exactly match obabel-emitted MOL2 atom order
  (CCDC sometimes drops disordered atoms / reorders by element). The
  CLASS-level confusion matrix is robust to this because we only count
  "did our detector flag any feature in this class" vs "did CCDC find
  any unusual feature in this class".
