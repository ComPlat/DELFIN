# Voll-Pool Code-State Map

**Purpose:** Map every voll-pool archive → exact code state. Enables retroactive
analysis "what code was in HEAD when this voll-pool ran?"

**Policy (active 2026-05-13):** commit + push BEFORE every voll-pool. Each archive
label maps 1:1 to a commit hash. See `feedback_commit_before_pool.md`.

Pre-policy iters (≤ 2026-05-13 morning) may have multiple voll-pools on the
same commit hash with different working-tree changes. For those, see the
patch column.

## Recovery commands

For commits in main branch history (column `Recovery` shows commit hash):
```bash
git checkout <hash> -- delfin/smiles_converter.py
# inspect, then restore current:
git checkout HEAD -- delfin/smiles_converter.py
```

For working-tree-only iters (column shows patch file):
```bash
git checkout 9c40b70 -- delfin/smiles_converter.py
git apply staged_patches/voll_pool_iters/<patch>.patch
# inspect, then restore current:
git checkout HEAD -- delfin/smiles_converter.py
```

## Full mapping

| Voll-pool archive | Date | Code commit | Patch file | Verdict |
|---|---|---|---|---|
| `HEAD-iter1-cand` | ~Apr | n/a (lost iter1) | — | early baseline |
| `HEAD-iter3-cand` / `iter3.1` / `iter3.2` | ~Apr | n/a (lost iter3.x) | — | early dev |
| `HEAD-iter4-cand` | ~Apr | n/a | — | early dev |
| `HEAD-iter5-cand` | ~Apr | n/a | — | early dev |
| `HEAD-iter6-cand` | ~Apr | n/a (env-gated patches OFF, baseline = 1e7eefe) | — | iter6 patches all dead-code |
| `HEAD-iter7-port1` / `-batch7flags` / `-selective5` | ~Apr | n/a | — | iter7 forensik |
| `HEAD-iter8.4a-port1` | ~May | n/a | — | iter8.4a smoke |
| `HEAD-iter8.4abc-vollpool` | 04-05 May | iter8.4abc commit (now bisect-confirmed regression) | — | REJECTED, root-cause of sigma drop |
| `HEAD-iter8.5a-vollpool` | 09-10 May | iter8.5a commit (e6761e4 port, class-dispatcher infra) | — | mostly OK, multi-σ 50% |
| `HEAD-iter8.6a-vollpool` | 10 May | iter8.6a commit (HAPTO_123A_PORT default-flip) | — | TOPO +16.62pp but hapto-intact catastrophe |
| `HEAD-iter8.6j-vollpool` / `-smoke500-baseline` | 11 May | iter8.6j commit (CN6-hetero enum depth) | — | OK baseline |
| `HEAD-iter8.10-vollpool` | 11 May | iter8.10 + iter8.10b commits | — | 11-May voll-pool baseline (used as "prev HEAD") |
| `HEAD-faf5f1b-port1` / `-vollpool` | 08 May | `faf5f1b` direct | — | multi-σ 48.5% champion (Wave-7 R analyzed) |
| `5dd520a-phase3a-vollpool` | 12 May AM | `5dd520a` (Phase 3B) | — | partial (2718 files) |
| `709c1a8-phase4d-vollpool` | 12 May PM | `709c1a8` (Phase 4D default-flip) | — | 11363 voll-pool, sigma +28.6pp smoke, multi-σ regression in voll-pool |
| `6ecd469-wave4b1-vollpool` | 12 May PM | `9c40b70` (Bundle 1 — committed AFTER pool ran) | — | net WIN: sigma 75.5% NEW HIGH, fail 9.4% NEW LOW |
| `9c40b70-wave5b2-vollpool` | 12 May night | `9c40b70` + Bundle 2 working-tree patches | `wave5b2_bundle2_9c40b70.patch` | REJECTED net loss 6w/9l/2t, M-D rollback over-aggressive |
| `9c40b70-wave6-fixA-vollpool` | 13 May AM | `9c40b70` + Fix A working-tree patch | `wave6_fixA_9c40b70.patch` | net ZERO per-class, 300s budget insufficient |
| `9c40b70-wave7b3-vollpool` | 13 May AM | `6c368de` (Wave-7 Bundle 3 — committed BEFORE pool finished) | — | RUNNING |

## Session 2026-05-13/14 — env-flag-configured pools

These pools share commits but differ by **env-flags** (the label encodes the
config). Recovery = `git checkout <commit>` + set the env-flags listed.

| Voll-pool archive | Commit | Env-flag config | JSONL | Verdict |
|---|---|---|---|---|
| `3943c2b-b5v2-vollpool` | `3943c2b` | `DELFIN_BAUSTEIN5=1` | `20260513_122130` | A: topo 47.64% (equal-n baseline) |
| `3943c2b-uffsoft-vollpool` | `3943c2b` | `DELFIN_UFF_SOFT_DONORS=1` | `20260513_122134` | B: topo 48.30% (+0.66pp, within ±1pp noise); bondlen real win |
| `fb4802f-b5v2-rigidH-vollpool` | `fb4802f` | `DELFIN_BAUSTEIN5=1 DELFIN_B5_RIGID_H=1` | `20260513_145341` | C: topo 46.55% (−1.09pp — rigid-H regresses) |
| `fb4802f-combined-vollpool` | `fb4802f` | `DELFIN_BAUSTEIN5=1 DELFIN_B5_RIGID_H=1 DELFIN_UFF_SOFT_DONORS=1` | `20260513_190750` | D: topo 46.21% (−1.44pp — "+6pp" headline was 300-file noise) |
| `fb4802f-b5v2regress-vollpool` | `fb4802f` | `DELFIN_BAUSTEIN5=1` (rigid-H default OFF) | `20260513_193052` | E: topo 45.97%; P-H-TRACK verified bit-exact default-OFF |
| `19821fb-champion-vollpool` | `19821fb` | `DELFIN_BAUSTEIN5=1 DELFIN_B5_RIGID_H=1 DELFIN_UFF_SOFT_DONORS=1 DELFIN_B5_PHASE_A5_SMART=1` (hard-C donor gate automatic) | `20260514_063038` | F: RUNNING — mitigated combined config |
| `19821fb-maxall-vollpool` | `19821fb` | all 16 Wave-1-7 flags ON (B5/rigidH/UFF-soft/PhaseA5-smart/multihapto-P1+P2/bridging-anion/CN5/CN3-CN4/high-CN/linear-CN2/mixed-metal/adaptive-timeout/multi-sigma-v2/class-aware-seeds) | `20260514_063039` | G: RUNNING — max-config |

**Key session finding:** the original C/D/E batteries scored only 300-file CLI
subsamples (±6pp noise). Re-scored with `--cli-sample 50000 --skip-xtb` →
`reports/all_metrics_<label>_FULL.json`. Equal-n verdict (63955 common frames):
raw patches do NOT improve topology; rigid-H + combined regress it. Only
UFF-soft is marginal-positive (within noise). F+G (mitigated configs) pending.

## Commit-graph (latest)

```
e69cc17 Archive: working-tree patches for past voll-pool iters on 9c40b70  ← THIS file added here
6c368de Wave-7 Bundle 3: champion code archeology — 3 evidence-based ports  ← wave7b3 code
9c40b70 Wave-4 Bundle 1: class-conditional Phase 4D defaults                ← wave4b1 code (post-pool commit)
6ecd469 Phase 6B: hapto BFS strict-fragmentation guard (env-gated default OFF)
8a18eea Phase 6A: multi-hapto centroid re-correction in CORRECT function
aad1f8c revert Phase 5B: unconditional _correct_hapto_geometry regressed multi-sigma
ceda1ec Phase 5B: remove secondary_metals gate on _correct_hapto_geometry
9c76e3d Phase 5A: universal H-bond-length normalization (env-gated default OFF)
709c1a8 Phase 4D: default-flip 3 env-flags from smoke500-validated wins   ← Phase4D voll-pool code
3734b50 Phase 4B: multi-sigma alt-binding wall-clock budget
606e799 Phase 4C: hash-seeded micro-jitter
a04ff14 Phase 4A: multi-hapto donor-drag guard
5dd520a Phase 3B: _class_conditional_flag helper                            ← Phase3A voll-pool code
c3c47cc Phase 3A: revert Iter-8.4abc sigma default-flips
...
```

## Pool result JSONLs (per-SMILES, in agent_workspace/quality_framework/results/)

| Archive | JSONL filename (chronological) |
|---|---|
| HEAD-faf5f1b-port1-vollpool | `smiles_master_v3__DELFIN__20260508_201624.jsonl` |
| HEAD-iter8.5a-vollpool | `smiles_master_v3__DELFIN__20260510_083239.jsonl` |
| HEAD-iter8.6a-vollpool | `smiles_master_v3__DELFIN__20260510_204726.jsonl` |
| HEAD-iter8.10-vollpool | `smiles_master_v3__DELFIN__20260511_132206.jsonl` |
| 5dd520a-phase3a-vollpool | `smiles_master_v3__DELFIN__20260512_132033.jsonl` (partial) |
| 709c1a8-phase4d-vollpool | `smiles_master_v3__DELFIN__20260512_155821.jsonl` |
| 6ecd469-wave4b1-vollpool | `smiles_master_v3__DELFIN__20260512_200413.jsonl` |
| 9c40b70-wave5b2-vollpool | `smiles_master_v3__DELFIN__20260512_235213.jsonl` |
| 9c40b70-wave6-fixA-vollpool | `smiles_master_v3__DELFIN__20260513_034338.jsonl` |
| 9c40b70-wave7b3-vollpool | `smiles_master_v3__DELFIN__20260513_071532.jsonl` (running) |

## Detector battery reports (in reports/)

`all_metrics_<archive>.json` per archive. Examples:
- `all_metrics_HEAD-iter8.10-vollpool.json`
- `all_metrics_709c1a8-phase4d-vollpool.json`
- `all_metrics_6ecd469-wave4b1-vollpool.json`
- `all_metrics_9c40b70-wave5b2-vollpool.json`
- `all_metrics_9c40b70-wave6-fixA-vollpool.json`
- `all_metrics_9c40b70-wave7b3-vollpool.json` (after wave7b3 completes)

## Wave-7 forensik reports (in ComPlat/DELFIN/iters/ and quality_framework/iters/)

- `WAVE4_A_md_break_forensik.md` — M-D break +7.71pp
- `WAVE4_B_hclash_forensik.md` — H-clash +11.44pp
- `WAVE4_C_topology_gate_design.md` — topology filter pre-return
- `WAVE4_D_multisigma_regression_forensik.md` — multi-σ -12.1pp
- `WAVE4_E_multihapto_class_forensik.md` — multi-hapto 0%
- `WAVE4_F_hapto_class_forensik.md` — hapto 17.6% stagnation
- `WAVE4_G_architecture_audit.md` — Säule 1 + router dormant
- `WAVE4_H_vdw_inter_clash_forensik.md` — vdW + inter clash
- `WAVE4_I_polyhedron_cshm_forensik.md` — polyhedron CSHM
- `WAVE4_J_bondorder_forensik.md` — bondord 0% aromatic (detector bug)
- `WAVE4_K_funcgrp_geom_forensik.md` — F3/F19/F20/F25
- `WAVE5_L_multisigma_regression.md` — wave4b1 multi-σ -12.1pp diagnosis
- `WAVE5_M_ligrealistic_diagnosis.md` — lig_real -6pp sampling noise
- `WAVE5_N_md_rollback_implementation.md` — P-MD-1 spec
- `WAVE5_O_multihapto_implementation.md` — SIMPLE_PATH + MM_ENFORCE spec
- `WAVE7_P_5687b3d_hapto_archeology.md` — HAPTO_NO_OB_WHEN_SIGMA finding (THE hapto breakthrough)
- `WAVE7_Q_81f8a1f_multihapto_archeology.md` — D-TIWWEB dead-code-fallback finding
- `WAVE7_R_faf5f1b_multisigma_archeology.md` — multi-σ runtime not chemistry
- `WAVE7_S_architecture_wirein.md` — Säule 1 + router wire-in spec
- `WAVE7_T_1e7eefe_topology_archeology.md` — debunked false champion
- `WAVE7_U_123a130_archeology.md` — sigma chelate-cap port default-flip
- `WAVE7_V_e6761e4_geometry_archeology.md` — sp3-C linear, σ-aryl OOP
- `WAVE7_W_cf1d480_merfac_archeology.md` — debunked false champion

## Per-iter journal (chronological event log)

`agent_workspace/quality_framework/results/iter_journal.md`

---

## Per-Class topology pass — voll-pool 11363 progression

```
run               sigma   hapto  multi-σ  multi-hapto  %fail  iso_ok
faf5f1b (08-May)  65.3%   17.6%   48.5%      0.0%      18.4%   47299
iter8.5a (10-May) 65.7%   18.1%   50.0%      0.0%      17.0%   23782 (partial)
iter8.6a (10-May) 49.9%   14.9%   30.0%      0.0%      37.6%   20166 (partial)
iter8.10 (11-May) 56.4%   17.6%   30.3%      0.0%      29.4%   43116
Phase4D  (12-May) 73.5%   17.4%   33.3%      0.0%      11.4%   56509
wave4b1  (12-May) 75.5%   17.6%   21.2%      0.0%       9.4%   58732  ← committed 9c40b70
wave5b2  (12-May) 73.5%   17.1%   27.3%      0.0%       9.4%   56105  ← rejected (Bundle 2)
wave6    (13-May) 75.4%   17.6%   21.2%      0.0%       9.5%   58652  ← reverted (Fix A)
wave7b3  (13-May) RUNNING (commit 6c368de) — expected hapto +N pp, sigma +N pp
```

## Rank-1 per metric (≥5k-file archives only, filtered)

Showing rank-1 winner among **voll-pool-sized archives** (excludes smoke/PROOF n=2 artifacts).

| Metric | Rank-1 winner | Value | HEAD-wave6 rank |
|---|---|---:|---:|
| per-frame topology % | `HEAD-iter8.6a-vollpool` | 67.19 | rank-16 |
| M-L intactness | `123a130` | 1.00 | rank-11 |
| M-D break files % | `HEAD-iter8.6a-vollpool` | 24.73 | rank-18 |
| h-clash files % | `81f8a1f` | 32.33 | rank-21 |
| h-clash per 1k H | `81f8a1f` | 9.79 | rank-22 |
| inter-clash files % | `e6761e4` | 57.00 | rank-11 |
| vdW clash % | `HEAD-iter8.10-vollpool` | 83.70 | rank-10 |
| CN intact % | `9c40b70-wave5b2-vollpool` | 99.08 | rank-5 |
| bondlen outlier % | `123a130` | 2.10 | rank-14 |
| F3 bond viol % | `123a130` | 28.79 | rank-16 |
| F19 H-Td viol % | `e6761e4` | 51.62 | rank-15 |
| **F20 H-planar viol %** | **`9c40b70-wave6-fixA-vollpool`** | **5.21** | **rank-1 ✓** |
| F25 sp3N pyr viol % | `6efa34e` | 33.15 | rank-20 |
| lig realistic % | `123a130` | 26.22 | rank-10 |
| polyhedron mean dev | `6efa34e` | 6.08 | rank-3 |
| CSHM max mean | `6efa34e` | 3.39 | rank-8 |
| extra bonds/frame | `123a130` | 0.06 | rank-19 |
| bondord correct % | `709c1a8-phase4d-vollpool` | 45.37 | rank-7 |
| stereo consistent % | `123a130` | 100.00 | tied |

**Caveats:**
- `HEAD-iter8.6a-vollpool` rank-1 on topo + M-D break is **detector-fooled** (per `feedback_8_6a_hapto_break.md`: hapto-intact actually 20.7%, false-pass via spurious chelate creation).
- All `topo_*` metrics are flagged **SUSPICIOUS / not CCDC-calibrated** per `feedback_metric_calibration_via_ccdc_zero.md`.
- `81f8a1f` rank-1 on h-clash is real per Wave-7 Q (champion-mechanism mineable).
- `123a130` rank-1 in 6 of 18 metrics (M-L intact, bondlen, F3, lig real, extra bonds, stereo) — Wave-7 U confirmed real champion, chelate-cap port code default-OFF in HEAD until Bundle 3.
- `e6761e4` rank-1 on inter-clash + F19 (Wave-7 V mixed champion).
- `6efa34e` rank-1 on F25 + polyhedron + CSHM.
- `HEAD-iter8.10-vollpool` rank-1 on vdW — prev HEAD baseline.

## Per-class topology pass rank-1 (real champion per class)

| Class | Rank-1 archive | Value | Source |
|---|---|---:|---|
| **sigma** | `9c40b70-wave6-fixA-vollpool` (≈ wave4b1 = HEAD committed) | **75.5%** | Bundle 1 win |
| **hapto** | `5687b3d` | ~32% (full archive) | Wave-7 P (HAPTO_NO_OB_WHEN_SIGMA=0 mechanism) |
| **multi-σ** | `faf5f1b` (08-May) | **48.5%** | Wave-7 R (runtime not chemistry — current 21.2%) |
| **multi-hapto** | `81f8a1f` (D-TIWWEB only) | 1/22 = 4.5% | Wave-7 Q (dead-code-fallback) |

**HEAD-wave6 per-class:** sigma rank-1 (75.5%), hapto rank-7-ish (17.6%), multi-σ rank-15-ish (21.2%), multi-hapto tied with most at 0%.

## Wave-7 forensik recommendations — current Bundle 3 covers

Bundle 3 (commit 6c368de, in voll-pool wave7b3 now) targets:
- **hapto rank-1** via Patch P 4A+4B (5687b3d mechanism port) — expected hapto 17.6% → 25-30%
- **sigma rank-1 maintain + extend** via Patch U 1 (123a130 chelate-cap default-ON) — expected -11.35pp gap close

Bundle 4 (deferred, evidence ready):
- **multi-hapto rank-1** via Patch Q (81f8a1f try-multiple-strategies)
- **multi-σ rank-1 attempt** via Patch R (300→540s budget)
- **e6761e4 sp3-C linear** via Patch V 5.3 (M-C-X 109.5° UFF constraint)
- **architecture activation** via Patch S (Säule 1 + ensemble_router wire-in)

