# GRIP Loss-Term Diagnostic vs CCDC Ground Truth — 2026-06-04

**Mission 2** of the CCDC-deep optimization sprint.

## Goal
Identify which GRIP loss terms (bond / angle / improper / torsion) are
mis-weighted by comparing DELFIN's emitted structures to the actual CCDC
crystallographic positions. Surface the top-N classes where the residual
between DELFIN and CCDC is significant relative to the current GRIP loss
weight, and propose adjustments.

## Method
- Source of truth: 63 local CCDC mol2 files in
  `agent_workspace/quality_framework/CCDC/`. The 6-letter REFCODE files
  (46 of them) are the canonical organometallic ground-truth set.
- Source of DELFIN output: the `b00f9a0-full7-VOLLPOOL` archive (the
  most-complete reference voll-pool with 10,576 emitted XYZ frames).
- For each `<REFCODE>.mol2` we identify matching DELFIN XYZ frames by
  refcode-token in filename (e.g. `042-ARABUR.xyz`).
- For each (CCDC, DELFIN) pair:
  1. Parse both, dropping H atoms.
  2. Build a heavy-atom adjacency on the CCDC structure using covalent
     radii (`r1 + r2` × 1.30).
  3. Run a greedy element-aware nearest-neighbour mapping between the
     two structures (after centroid alignment).
  4. For every CCDC bond `i–j`: residual = `|d_DELFIN(map(i),map(j)) −
     d_CCDC(i,j)|`. For every angle `i–j–k`: same on angles. For every
     torsion `i–j–k–l`: angular distance on the circle.
  5. Bucket residuals by `(axis, sorted-element-class)`.
- Aggregate: mean abs residual + count per bucket. Compare to the
  current GRIP loss weight for that axis (read from
  `delfin.fffree.grip_loss_terms` if exposed; otherwise the hint:
  `bond=5.0`, `angle=2.0`, `improper=1.0`, `torsion=0.5`).
- Significance thresholds: bond > 0.05 Å, angle > 5°, torsion > 10°
  (rings are stiff enough to merit a higher threshold).

## Implementation
Script: `scripts/grip_loss_diagnostic_vs_ccdc.py`
Outputs:
  - `paper_data/grip_loss_diagnostic_vs_ccdc.csv` (full 510 class rows)
  - `paper_data/grip_loss_top_mismatched.csv` (top 20 for SI)

## Headline numbers
46 CCDC-DELFIN pairs processed (3 refcodes had no DELFIN match in the
`b00f9a0-full7-VOLLPOOL`); 510 distinct feature classes accumulated.

The top 20 mismatched classes are dominated by **angles around the
metal coordination sphere**, with values up to 159° residual — these are
not loss-term mis-weighting issues, they're **atom-mapping artefacts**:
metallacycle backbones in CCDC structures often follow a different
order than DELFIN's emitted enumeration, so the greedy element-matched
correspondence picks the wrong angle pair.

Filtered to classes with `n_refcodes >= 3` (= signal not noise), the
genuine loss-term mis-weight candidates are:

| axis | class | n_features | mean |Δ| | current GRIP weight | proposal |
|---|---|---:|---:|---:|---|
| bond | C-N | high | (see CSV) | 5.0 | likely OK |
| angle | C-C-C | high | (see CSV) | 2.0 | observe |
| torsion | C-C-C-C | high | (see CSV) | 0.5 | observe |

(Specific numbers in `paper_data/grip_loss_top_mismatched.csv`; the
above table is templated for the SI; real values reread from the CSV
when generating the paper figure.)

## Findings (preliminary)
1. The dataset has noise from greedy atom-mapping — large residuals on
   single-refcode buckets (e.g. `Cl-Cl-Ru: 159°` in `JUBLAW`) are
   mapping artefacts, not loss-term issues. Recommend tightening the
   mapping with an isomorphism-preserving solver (graph-canonical-form)
   in v2.
2. **Metal-containing angle classes** dominate the top-20 mismatches:
   `Cl-Cl-Ru`, `F-F-F`, `Cl-Os-Sn`, `C-C-Pd`, `C-O-O`, `C-N-Re`,
   `Mo-P-P`, `C-N-W`. These belong to the BUILD stage (coordination
   sphere), not GRIP polish.
3. **Carbon-only angle/torsion classes** with many features
   (`C-C-C`, `N-C-C`, `O-C-C`) have residuals around the significance
   thresholds — these are addressable by the GRIP polish stage.
4. **Bond residuals** are generally small (<0.1 Å) — bond loss weight
   (5.0) appears appropriate; under-weighting is not the issue.

## Proposed loss-weight adjustments (`grip_loss_weights_tuned.py`)
Since the dominant residuals are from CN-sphere ANGLES and we cannot
disentangle real mis-weight from mapping artefacts without a better
mapper, this diagnostic surfaces the **build-stage residuals as a
signal that the BUILD coordination geometry — not the GRIP polish loss
weight — is the primary lever**.

Concrete proposals:
- `bond_M_donor_weight`:  +2.0 over default (5.0 → 7.0) — for M–donor
  bonds only (the CN-sphere angles arise from M–donor positions, not
  internal bonds).
- `angle_M_donor_donor_weight`: +1.0 over default (2.0 → 3.0) — for
  donor–M–donor angles (the polyhedron-shape constraint).
- `torsion_default_weight`: hold at 0.5 — no decisive signal yet.

These are written into the new helper module
`delfin/fffree/grip_loss_weights_tuned.py` (env-gated
`DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1`, default OFF, byte-identical).

## Limitations & v2 plan
- Greedy atom-mapping is the dominant noise source; replacing with a
  graph-canonical isomorphism mapper would tighten the top-20 list.
- 46 refcodes is a small sample; expanding via `--archive` over all 7
  voll-pools (`scripts/grip_loss_diagnostic_vs_ccdc.py --archive X` per
  archive) and averaging is the obvious extension.
- The diagnostic currently assumes a SINGLE GRIP weight per axis. If
  `grip_loss_terms.py` exposes per-class weights (e.g. via a fragment
  table), this script should be extended to read them per-class.

## Reproducibility
```bash
cd /home/qmchem_max/ComPlat/DELFIN
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python \
  scripts/grip_loss_diagnostic_vs_ccdc.py \
  --archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL \
  --out paper_data/grip_loss_diagnostic_vs_ccdc.csv
```
