# Experiment Manifest — 2026-06-06

All experiments deterministic (`PYTHONHASHSEED=0`, ETKDG seed=42).
Reproducible from launch scripts in this directory.

## Voll-pool runs

### ULTIMATE-2026-06-06-VOLLPOOL — ❌ FAIL net -79 vs 2792332

- **Launch script**: `voll-pool-ULTIMATE-pipeline.sh` (and `scripts/run_vollpool_ULTIMATE.sh` in repo)
- **GUPPY HEAD**: `82baa4c`
- **SMILES pool**: `pools/smiles_master_v3_plus.txt` (11458 SMILES)
- **XYZ archive**: `agent_workspace/quality_framework/xyz_archive/ULTIMATE-2026-06-06-VOLLPOOL/` (11453 XYZ)
- **Forensik dispatch TSV**: `ULTIMATE-2026-06-06-VOLLPOOL_dispatch_forensik.tsv` (11453 entries)
  - Dispatch breakdown: native=10346 (90.3%), embed-grip=891 (7.8%), embed-failed=216 (1.9%)
- **All-metrics**: `reports/all_metrics_ULTIMATE-2026-06-06-VOLLPOOL.json`
- **Pipeline log**: `voll-pool-ULTIMATE-pipeline.log`
- **iter_gate vs 2792332**: net **-79**, 78 severe, HEAL-FIRST verdict
- **Equal-n forensik subagent verdict**: `ULTIMATE_equal_n_verdict.md` (14KB) — proxy equal-n net -76, native subset -62
- **Per-dispatch quality split** (subagent forensik):
  - native subset (8833 XYZ): net -62, 75 severe
  - embed-grip subset (760 XYZ): net -50, 66 severe
  - embed-failed subset (153 XYZ): net +38 (trivial, all-zeros)
- **SMILES→dispatch mapping**: `xyz_to_dispatch.tsv` (9746 mapped of 11453 = 85%)

## Smokes — bisect series

### SMOKE-F1-OFF-2026-06-06 — F1 OFF control bisect

- **Launch script**: `smoke_F1_OFF.sh`
- **GUPPY HEAD**: `82baa4c` (but with `DELFIN_FFFREE_F1_COVERAGE` UNSET)
- **SMILES pool**: `/tmp/shared_smiles.txt` (550 SMILES)
- **XYZ archive**: `agent_workspace/quality_framework/xyz_archive/SMOKE-F1-OFF-2026-06-06/` (549 XYZ)
- **Forensik dispatch**: `SMOKE-F1-OFF-2026-06-06_dispatch_forensik.tsv`
  - Dispatch: native=150, embed-grip=330, embed-failed=70
- **All-metrics**: `reports/all_metrics_SMOKE-F1-OFF-2026-06-06.json`
- **Log**: `smoke_F1_OFF.log`
- **iter_gate vs ULTIMATE-intersect**: net **+46**, 19 severe regressions, HEAL-FIRST verdict
- **Verdict**: F1 IS primary regression source. F1-OFF beats ULTIMATE on internals.

### SMOKE-F1-HEAL-2026-06-06 — F1 + structqual gate v1

- **Launch script**: `smoke_F1_HEAL.sh`
- **GUPPY HEAD**: `82baa4c` + local F1-HEAL v1 patch (NOT committed, since reverted)
- **Patches** (now reverted):
  - `delfin/fffree/converter_backend.py`: `_f1_structqual_pass()` heavy_collapse/extra_fragments/topo_missing filter + `_F1_NO_LAST_RESORT` disable
- **Env**: `DELFIN_FFFREE_F1_STRUCTQUAL_GATE=1`, `DELFIN_FFFREE_F1_NO_LAST_RESORT=1`, `..._HEAVY_COLLAPSE_FACTOR=0.7`
- **SMILES pool**: `/tmp/shared_smiles.txt` (550)
- **XYZ archive**: 550 XYZ
- **Forensik dispatch**: `SMOKE-F1-HEAL-2026-06-06_dispatch_forensik.tsv`
  - Dispatch: native=149, embed-grip=329, embed-failed=72
- **Reason log**: `SMOKE-F1-HEAL-2026-06-06_reason.tsv` — 144 `f1_structqual_rejected` entries (ALL extra_fragments, ZERO heavy_collapse, ZERO topo_missing)
- **iter_gate vs ULTIMATE-intersect**: net **+36** (worse than F1-OFF's +46)
- **iter_gate vs F1-OFF**: net **-57** with 19 severe regressions
- **Verdict**: ❌ F1-HEAL v1 is NET-NEGATIVE vs F1-OFF. structqual reject removed valid hapto-π isomers (Cp/arene appear as separate component in RDKit bond graph). Patches reverted.

### SMOKE-F1F2-HEAL-2026-06-06 — F1-HEAL v2 + F2-mddir snap (interim)

- **Launch script**: `smoke_F1F2_HEAL.sh`
- **GUPPY HEAD**: `82baa4c` + local F1-HEAL v2 + F2-mddir patches (NOT committed, reverted before this smoke completed)
- **Code state during smoke**: workers cached pre-revert code (multiprocessing fork) — result validity needs careful interpretation
- **Patches** (now reverted):
  - F1-HEAL v2: hapto-π proximity whitelist (components within 2.5Å of metal not counted as "extra")
  - F2-mddir: M-D snap pre-GRIP polish in `embed_fallback.py` `_maybe_md_snap()` (chelate-aware subtree translation)
- **Env additions**: `DELFIN_FFFREE_F2_MD_SNAP=1`, `_LO=0.80`, `_HI=1.20`
- **XYZ archive**: 549-550 XYZ (in progress at archive time)
- **Forensik**: `SMOKE-F1F2-HEAL-2026-06-06_dispatch_forensik.tsv`, `SMOKE-F1F2-HEAL-2026-06-06_reason.tsv`
- **Status**: smoke completion + iter_gate pending at archive time

## Reverted patches — NOT in code base

The F1-HEAL v1, v2, and F2-mddir patches were reverted externally (system-marked intentional) after smoke runs. They remain documented above + in launch scripts + in smoke logs.

The structqual-gate approach was net-negative vs F1-OFF → wrong heal direction.

## Active subagent branches (in progress)

- `feat-phase2-grip-md-unfreeze-2026-06-06` — GRIP M+D unfreeze + CCDC-empirical priors (modifies `grip_polish.py`, new `grip_md_terms.py`)
- `feat-binding-mode-enum-2026-06-06` — coordination-mode enumeration wiring (modifies `coordination_mode_enum.py`, `converter_backend.py`)

## Reproducibility checklist (for any future run)

1. Set `PYTHONHASHSEED=0`
2. Use exact GUPPY commit hash from manifest entry
3. Use exact env-flag list from launch script
4. Use exact SMILES pool file
5. Use `--parallel`, `--timeout`, `--no-retry-on-zero` as in launch script
6. Compare via `iter_gate.py <label> --prev <prev-label>` for net-difference
7. Detector battery: `run_all_detectors.py <archive> --all --skip-xtb --cli-sample 50000 --parallel 16`

## Cross-references

- Memory: `feedback_ULTIMATE_voll_pool_FAIL_2026_06_06.md`
- Memory: `feedback_HONEST_xrd_rmsd_findings_2026_06_06.md`
- Memory: `project_F2_merged_with_polish_budget_2026_06_05.md`
- Memory: `project_night_autonomous_2026_06_05_to_06_06.md`
