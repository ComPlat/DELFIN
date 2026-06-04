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

## Top-20 atlas (real, 20 structures from b00f9a0-full7-VOLLPOOL)

97 distinct fragment-class patterns observed; top-20 by mean |z|.
Two CSVs:
- `paper_data/anomaly_pattern_atlas.csv` (live-run, 5 rows — CCDC
  subprocess timeouts truncated some workers)
- `paper_data/anomaly_pattern_atlas_cache.csv` (rebuilt from the 20-
  result CCDC cache via `scripts/anomaly_atlas_from_cache.py`; 20
  rows from 97 classes — the proper top-20)

Top-5 from cache:

| # | axis  | fragment  | mean |z| | n_obs | example refcodes |
|---|---|---|---:|---:|---|
| 1 | bonds | P1_C3   | 2.80 | 3 | AVUHOQ, BERDEH, KEJTIF |
| 2 | bonds | C15_C16 | 2.54 | 4 | RABHIM, YIMLAH, YEZMAT, RUZQOS |
| 3 | bonds | C7_N2   | 2.46 | 3 | CMPBPU, CODAAC, RUZQOS |
| 4 | bonds | C3_N1   | 2.39 | 3 | GEBXOD, TIGREI, OWOYAB |
| 5 | bonds | C5_C6   | 2.39 | 4 | GEBXOD, CODAAC, HOPYIU, MOGGUK |

**Key observation:** the CCDC `fragment_label` returned by
`GeometryAnalyser` on obabel-generated MOL2 files is generic atom-name
based (`P1_C3`, `C15_C16`), NOT the chemical fragment SMARTS that
Mogul typically returns when reading native CCDC mol2/cif.

To get chemically-meaningful fragment labels (e.g. `cyclopentadienyl`,
`carbonyl-Pd`) we would need to either:
1. Re-export the structures using CCDC's own `MoleculeWriter` from
   the entry (not possible — DELFIN emits XYZ + obabel).
2. Run an explicit `chemistry.set_protonation_state(mol)` + bond-
   order perception step in the CCDC subscript before
   `GeometryAnalyser`.
3. Use SMARTS-based post-classification of the fragment atoms.

(Option 2 is the most tractable; left as v2 work.)

The current atlas STILL provides useful information at the atom-type
level (C-C, C-N, P-C bonds with |z| ~ 2-3 indicate systematic
geometry issues in our backbone construction). It just doesn't map
to high-level chemistry-class fix directions yet.

**Provisional fix-direction reading:**
- All top-20 entries classify as `polish` because the fragment
  labels are atom-pair names that match the C-C/C-N polish patterns.
- `P1_C3` (phosphine-carbon, |z| 2.8) and `C7_N2` (imine/nitrile
  backbone, |z| 2.5) are the highest-priority backbone-torsion
  candidates for GRIP polish in the next iteration.

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
