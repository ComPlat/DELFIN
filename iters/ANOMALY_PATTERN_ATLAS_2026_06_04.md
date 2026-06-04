# Anomaly Pattern Atlas — 2026-06-04

**Mission 4** of the CCDC-deep optimization sprint.

## Goal
Catalogue the top-N most-anomalous fragment patterns in a target voll-
pool using CCDC's `GeometryAnalyser`, classify each into a fix-direction
(`build` / `polish` / `library`), and surface example refcodes so the
next iteration knows where to focus.

## Method
1. Stratified random sample of 500 XYZ structures from the target
   archive (eaee0c2-new-stack-VOLLPOOL or ec7fb0d-full9-fblock-VOLLPOOL
   when the new-stack pool is incomplete).
2. For each structure: convert to MOL2 (obabel), run
   `GeometryAnalyser.analyse_molecule(mol)`, collect every analysed
   feature (bond / angle / torsion / ring) with `enough_hits` and a
   defined `z_score`.
3. Bucket features by `(axis, fragment_label)` (the canonical CCDC
   fragment string).
4. For each bucket compute:
   - count of observations
   - mean abs(z)
   - max abs(z)
   - fraction unusual (|z| >= 2)
   - 5 example refcodes (extracted from the source XYZ filename)
5. Rank by mean abs(z), take top-20.
6. For each top-20 pattern, classify the suggested **fix direction**
   via regex over the fragment label:
   - `build`: M-D bond / eta-hapto pattern → build-time fix
   - `polish`: ring / aromatic / C-C / C-X / X-H → GRIP polish
   - `library`: explicit `L0` or `generalised` → expand COD fragment table

## Implementation
- `scripts/anomaly_pattern_atlas.py` — the runner. Reuses the CCDC
  subprocess pattern from `ccdc_mogul_validation.py`.
- Output CSV: `paper_data/anomaly_pattern_atlas.csv` (top-N)
- Optional JSONL: full ranked list (`--jsonl`)

## Top-20 atlas (template; populated when CCDC run completes)
| # | axis | fragment | mean|z| | n_obs | example refcodes | fix |
|---|---|---|---:|---:|---|---|
| 1 | bond | (TBD) | (TBD) | (TBD) | (TBD) | build |
| 2 | ... | ... | ... | ... | ... | ... |

(Concrete numbers in `paper_data/anomaly_pattern_atlas.csv` when the
job completes.)

## Reproducibility
```bash
cd /home/qmchem_max/ComPlat/DELFIN
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python \
  scripts/anomaly_pattern_atlas.py \
  --archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive/ec7fb0d-full9-fblock-VOLLPOOL \
  --n 500 --workers 6 \
  --out paper_data/anomaly_pattern_atlas.csv \
  --jsonl paper_data/anomaly_pattern_atlas_full.jsonl
```

## How to use the atlas in the next iteration
1. Filter rows where `fix_direction == build`. These are the structures
   the build stage is getting wrong. Match the example refcodes to the
   isomer / hapto / polyhedron audit and add specific test cases.
2. Filter rows where `fix_direction == polish`. These are addressable
   by tuning GRIP-loss weights (per the Mission-2 helper
   `grip_loss_weights_tuned.py`).
3. Filter rows where `fix_direction == library`. These indicate
   under-represented chemistry in the COD fragment library — schedule
   `metric_cod_fragment_lookup.py` re-build with expanded sources.

## Default state
- Read-only forensic tool.
- No production code changes (the atlas is consumed by humans /
  agents, not by the runtime).
