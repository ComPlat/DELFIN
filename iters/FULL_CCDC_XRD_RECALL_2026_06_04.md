# Full-CCDC XRD-Recall Metric — 2026-06-04

## User-Direktive

> "XRD abgleich ab jetzt mit einbeziehen in generierungsqualität ist
> richtige Isomer + konformere enthalten."  
> "48 CCDC ground truth ist unnötig wir haben jetzt alle ccdc stukturen"

## Branch

`full-ccdc-xrd-recall-metric` (worktree-isolated; no GUPPY/main merge).

## Scope vs Previous Work

The earlier `delfin/fffree/xrd_recall_metric.py` (main branch, 488 LOC)
operated against a hardcoded 48-refcode JSON dump
(`ccdc_ground_truth_48.json`). The user replaced that scope: we now use
the **full CSD 2026.1** (1,436,119 entries) via the Python API at build
time, with a ~500k-entry TMC subset shipped as static JSONL+NPZ at
runtime (CCDC-frei).

## Architecture

```
                                  build-time (CCDC API)
   CSD 2026.1                                 |
   1,436,119 entries  -->  build_ccdc_tmc_index.py
                                              |
                                              v
                            ccdc_tmc_index.jsonl  (~150-300 MB, ~500k records)
                            ccdc_tmc_index.npz    (compact lookup, ~10 MB)
                                              |
   smiles_master_v3_plus.txt    -->  build_per_smiles_ccdc_families.py
   (11458 SMILES)                             |
                                              v
                            per_smiles_ccdc_families.jsonl  (one line per SMILES)
                                              |
   ---------------------------------------------------------------- (legal cut)
                                              |
                                  runtime (CCDC-free)
                                              |
                            delfin.fffree.xrd_recall_metric_full
                                              |
                          compute_isomer_recall(emitted, table)
                          compute_conformer_recall(emitted, table, index, rmsd<0.5)
                          per_smiles_recall_report(...)
                          score_archive(archive_dir)
                                              |
                          _xrd_recall_report.jsonl + per-SMILES JSONL
                          (auto-emitted from pool_evaluator.py via env hook)
```

## Definition

For each SMILES with ≥1 CCDC family member:

```
isomer_recall(smiles) =
    |DELFIN-emitted-isomers ∩ CCDC-family-isomers (mod CN-skeleton)|
    / |CCDC-family-isomers|

conformer_recall(smiles, rmsd_threshold=0.5) =
    fraction of CCDC family refcodes for which SOME DELFIN emission has
    heavy-atom Kabsch RMSD < threshold (Å)
```

Global metric = mean over SMILES with ≥1 CCDC match (option:
weighted by family size).

CCDC-side isomer-class for a refcode is its coordination-shape label
(`octahedral_or_TP`, `tetrahedral_or_squareplanar`, `TBP_or_SPY`, …).
DELFIN-side isomer-class comes from `classify_isomer(syms, P)` which
returns `CN<n>-<polyhedron>` (e.g. `CN6-OH-fac`, `CN4-SP-cis`). The
matching uses CN equivalence + polyhedron-keyword substring matching
(`_isomer_match()`).

## Per-SMILES Match Lookup (Mission 2)

Two paths to a CCDC family per SMILES:

1. **Refcode-prefix**: master pool labels like `042-ARABUR|...` directly
   point at refcode `ARABUR`.
2. **Graph-signature**: `(metal, CN, donor_element_set, n_heavy_bucket)`
   computed regex-only from the SMILES, looked up in a pre-indexed
   `by_sig` map of the TMC index.

Method 1 is exact but covers only ~98 / 11458 SMILES. Method 2 catches
the rest at graph-isomorphism-by-signature level — coarser than full
substructure match, but orders-of-magnitude faster than
`SimilaritySearch` (which ~9 s / SMILES → 28 h on the full pool would be
prohibitive). The signature can be refined in v2 with InChIKey of the
metal-shell-substructure once we accept the build cost.

## Env-Flags (production runtime)

| Flag | Default | Effect |
| --- | --- | --- |
| `DELFIN_FFFREE_XRD_RECALL_TABLE_PATH` | (unset) | Path to per_smiles_ccdc_families.jsonl |
| `DELFIN_FFFREE_XRD_RECALL_INDEX_PATH` | (unset) | Path to ccdc_tmc_index.jsonl |
| `DELFIN_FFFREE_XRD_RECALL_RMSD_A` | 0.5 | Heavy-atom RMSD pass threshold (Å) |
| `DELFIN_FFFREE_XRD_RECALL_WEIGHTED` | "0" | "1" = weight global recall by family-size |
| `DELFIN_FFFREE_XRD_RECALL_MASTER_POOL` | smiles_master_v3_plus.txt | Master pool for label→smiles map |
| `DELFIN_QF_AUTO_XRD_RECALL` | (unset) | "1" = auto-run after pool emission |

All flags default to OFF byte-identical: importing
`delfin.fffree.xrd_recall_metric_full` is a no-op pure-python+numpy
operation.

## Build-Time Pipeline

```sh
# Mission 1 (60-100 min, single-threaded EntryReader iteration)
PYTHONHASHSEED=0 run_csd_python_api \
  agent_workspace/quality_framework/scripts/build_ccdc_tmc_index.py

# Mission 1b: pack compact NPZ
python agent_workspace/quality_framework/scripts/pack_ccdc_tmc_npz.py

# Mission 2 (~30s, in-memory join)
PYTHONHASHSEED=0 python \
  agent_workspace/quality_framework/scripts/build_per_smiles_ccdc_families.py
```

## Runtime Usage

```python
from delfin.fffree import xrd_recall_metric_full as XRDM

# Score a pool archive
res = XRDM.score_archive(
    "/path/to/archive",
    table_path="/path/to/per_smiles_ccdc_families.jsonl",
    index_path="/path/to/ccdc_tmc_index.jsonl",
)
# res = {'isomer_recall': 0.42, 'conformer_recall': 0.18,
#        'n_smiles_with_ccdc_family': 7341, ...}
```

Via env-flag from pool_evaluator:

```sh
DELFIN_QF_AUTO_XRD_RECALL=1 \
DELFIN_FFFREE_XRD_RECALL_TABLE_PATH=/path/to/families.jsonl \
DELFIN_FFFREE_XRD_RECALL_INDEX_PATH=/path/to/index.jsonl \
python pool_evaluator.py ...
# emits _xrd_recall_report.jsonl + _xrd_recall_per_smiles.jsonl in archive dir
```

## Tests

`tests/test_xrd_recall_metric_full.py` — 27 tests passing,
`PYTHONHASHSEED=0`, no CCDC dependency. Coverage:

- Loader smoke (missing path, env-var, cache hit, JSONL key resolution)
- Classifier fallback (octahedral, no-metal)
- Kabsch RMSD (identity, rotation-invariance)
- Recall computations (perfect, distorted, no-family, no-emission,
  unweighted vs weighted)
- Per-SMILES report structure
- Master-label-map sanitization + archive resolution
- CN-aware isomer matching across CCDC↔DELFIN label conventions

## Paper SI (draft text)

> **DELFIN-vs-CSD XRD-recall.** For each generated pool we measure
> two metrics against the full CSD 2026.1 (1,436,119 entries). A TMC
> subset (~500k d/f-block monomeric refcodes) is built once at build
> time via the CCDC Python API and shipped as static JSONL. For each
> input SMILES we precompute its CSD family — refcodes that share the
> `(metal, CN, donor-element-set, heavy-atom-bucket)` graph signature.
> **Isomer-recall** is the fraction of CCDC-family isomer-classes that
> appear in the DELFIN emission; **conformer-recall** at 0.5 Å is the
> fraction of CCDC-family refcodes for which some DELFIN emission has
> heavy-atom Kabsch-RMSD < 0.5 Å. DELFIN achieves **N%** isomer-recall
> and **M%** conformer-recall, averaged over **X** SMILES with ≥1 CCDC
> match. The metric is implementable as a runtime CCDC-free callable
> (`delfin.fffree.xrd_recall_metric_full`) and integrates directly into
> the pool evaluator as a post-pool hook.

(N, M, X to be filled by Mission 4 measurement.)

## Legal Separation

Doctrine: CCDC code/data ONLY at our build time. Runtime ships static
JSONL + NPZ; downstream adopters receive these files without any CCDC
import. The runtime module's `try: from … xrd_recall_metric import …`
falls back to a pure-numpy implementation when the sibling is absent.

## Status (2026-06-04 18:00)

- Mission 1 in progress (~100k / 1.44M, ~95 min ETA, 230 entries/s)
  → ~500k TMC refcodes expected on completion
- Mission 2 script committed, ready to run on completion of M1
- Mission 3 module + 27 pytest committed (8ffcfc7)
- Mission 4 ready to launch once M2 complete
- Mission 5 pool_evaluator hook committed in agent_workspace
- Mission 6 = this document

## Next Steps (post-completion)

1. Run Mission 1 to completion (~90 min) → `ccdc_tmc_index.jsonl`
2. Pack NPZ via `pack_ccdc_tmc_npz.py`
3. Run Mission 2 → `per_smiles_ccdc_families.jsonl`
4. Run Mission 4 measurement on b00f9a0, c03a550, ec7fb0d (and others)
5. Update Paper SI with N/M/X numbers
6. Open PR from `full-ccdc-xrd-recall-metric` for review (no auto-merge to GUPPY/main)
