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

## Results — smoke run (15 structures, b00f9a0-full7-VOLLPOOL)

2 classes tuned (small-N caveat — production needs n>=200 × 5+ archives):

| axis  | class           | optimal threshold     | recall | precision | F1   |
|---|---|---:|---:|---:|---:|
| bonds | bonds::C-C    | 1.5 (vs default 2.5) | 1.00 | 0.25 | 0.40 |
| angles| angles::C-C-C | 3.0 (vs default 2.5) | 1.00 | 0.33 | 0.50 |

**Reading:**
- **C-C bonds: CCDC is STRICTER than our v3** — default 2.5 misses
  unusual C-C bonds CCDC flags; lowering to 1.5 hits 100% recall.
- **C-C-C angles: CCDC is MORE LENIENT** — default 2.5 flags C-C-C
  angles CCDC considers normal; raising to 3.0 cuts false positives.

Consistent with chemical intuition: C-C bond lengths are tight (ring
+ chain distributions), C-C-C angles span sp/sp2/sp3 broadly.

## Production claim
> On the b00f9a0-full7-VOLLPOOL smoke set (15 structures), per-class
> threshold tuning identifies two classes where mogul_detector_v3's
> global threshold (2.5) is mis-calibrated: bonds::C-C should use 1.5
> (recall 1.0, precision 0.25), and angles::C-C-C should use 3.0
> (recall 1.0, precision 0.33). The production module
> `delfin.fffree.mogul_detector_v3_tuned` consumes
> `paper_data/mogul_v3_threshold_optimization.csv` and applies the
> per-class thresholds when `DELFIN_MOGUL_V3_TUNED=1`. A larger run
> (n>=200 × 5 archives) is needed for the full per-class table.

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
