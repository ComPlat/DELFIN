# CCDC Reference Dataset — 2026-06-04

**Mission 5** of the CCDC-deep optimization sprint.

## Goal
Build a permanent CCDC-derived reference dataset that the production
runtime can use WITHOUT a CCDC license — the legal separation per
`[[feedback_ccdc_legal_separation_doctrine_2026_06_04]]`:

> CCDC is used by US AT BUILD TIME. Adopters / production runtime
> consume only the DERIVED ARTEFACTS (statistics, aggregates).

## Method
Script: `scripts/build_ccdc_ground_truth_48.py`

For each REFCODE in
`agent_workspace/quality_framework/CCDC/_VALID_CALIBRATION_SET.txt`
and `_HAPTO_TEST_SET.txt` (and any 6-letter `*.mol2` filename in
`agent_workspace/quality_framework/CCDC/`):

1. Load the entry via `ccdc.io.EntryReader('CSD').entry(refcode)`.
2. Extract heavy atoms + 3D positions + bond list + bond orders.
3. Identify the primary metal centre + coordination sphere via
   covalent-radii cutoffs.
4. Classify the polyhedron (CN + family) + heuristic isomer label
   (fac/mer/cis/trans/...) — same convention as DELFIN's Pólya output.
5. Run `ccdc.conformer.GeometryAnalyser().analyse_molecule(mol)` and
   record per-axis baselines (mean |z|, max |z|, p95 |z|).
6. Emit a per-refcode JSON record into the bundle.

## Two output artefacts (legal classification)

### A) Full per-refcode dump
`agent_workspace/quality_framework/reference/ccdc_ground_truth_48.json`

Contains per-refcode atomic positions, bonds, isomer labels, and Mogul
baseline statistics.

**Legal status:** internal-use only (build-time + research). NOT
shipped to adopters because individual atomic coordinates derived from
CSD entries fall under the CCDC redistribution policy.

This is the artefact consumed by `delfin.fffree.xrd_recall_metric`
during development / paper-figure generation.

### B) Aggregate-only NPZ
`agent_workspace/quality_framework/reference/ccdc_ground_truth_aggregate.npz`

Contains ONLY statistical aggregates:
- per-CN distribution
- per-metal-element distribution
- per-isomer-label distribution
- per-axis mean-|z| arrays (one float per refcode, anonymised)

**Legal status:** CCDC-derived but only **derived statistics** — no
individual coordinates leak. Shippable per CCDC's terms, which permit
aggregate / statistical re-distribution.

This is the artefact production runtime would import behind
`DELFIN_USE_CCDC_REFERENCE_DUMP=1`.

## Implementation notes
- The script lives in `scripts/`, NOT in `delfin/`, so the production
  package never imports CCDC.
- Output paths are explicit; the script does not modify `delfin/` and
  does not write to any existing artefact directory.
- Determinism: refcode list is sorted, output is dict-ordered, JSON
  emitted with `indent=2`.

## Reproducibility
```bash
PYTHONHASHSEED=0 \
  /home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api \
  /home/qmchem_max/ComPlat/DELFIN/scripts/build_ccdc_ground_truth_48.py
```

Expected runtime: ~10-20 min (CCDC `analyse_molecule` is the bottleneck:
~5-15 s per refcode × 66 refcodes).

## Aggregate NPZ schema
```
n_refcodes          (int32 scalar)
cn_keys             (int32 array, sorted unique CN)
cn_counts           (int32 array, per-CN count)
metal_keys          (object array, sorted unique metal symbols)
metal_counts        (int32 array)
iso_keys            (object array, sorted unique isomer labels)
iso_counts          (int32 array)
poly_keys           (object array, sorted unique polyhedron classes)
poly_counts         (int32 array)
bond_meanz_per_ref  (float32 array, one mean |z| per refcode with ≥5 features)
angle_meanz_per_ref (float32 array)
torsion_meanz_per_ref (float32 array)
ring_meanz_per_ref  (float32 array)
```

## Default state
- Build-time tool only. Runtime delfin/ stays CCDC-free.
- `DELFIN_USE_CCDC_REFERENCE_DUMP=1` activates consumption.
- Without the env flag, `delfin.fffree.xrd_recall_metric.load_reference`
  returns `None` (graceful degradation; the metric simply doesn't
  compute).
