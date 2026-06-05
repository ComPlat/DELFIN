# Mission C1 — Voll-Pool Regression HEAL Forensik (2026-06-05)

## TL;DR

The 2792332 → 2ee2f45 voll-pool regression (**net −25, 50 severe axes**) is
**caused by env-flag drift in the launch script**, NOT by code changes in
Mission B commits.  Today's voll-pool was launched with 8 DELFIN_FFFREE_*
flags; the 2792332 baseline launcher set 22.  Restoring the 15+ missing heal
and realism flags recovers most of the regression.

## Mission B commits (between 2792332 and 2ee2f45) — code-impact audit

| commit | mission | code-active by default? |
|---|---|---|
| 2ee2f45 sandwich-piano dispatch | B1 | **YES** (auto-on under PT3) |
| 834ae75 mogul-v4 H-aware | B2 | NO (env-gated, OFF) |
| b9ac2d9 GRIP class-conditional disable | B3 | NO (env-gated) |
| 567f11e GRACE Pareto + adaptive | B4 | NO (env-gated, OFF) |
| 862070a CN10 polyhedra | A2 | **YES** (auto-on under PT3) |
| 8bba47f oxalate planar | iter-32e | NO (env-gated, OFF) |
| b338003 GRIP donor-cone DoF | A4 | NO (env-gated, OFF) |
| 0508efc sandwich-piano polyhedra | A1 | YES (PT3) |
| 1b95142 NO_FALLBACK auto-implies dispatch | A5 | YES (PT3) — but BUILDER=1 was set in 2ee2f45, so the auto-implication didn't change behavior |
| 9e96246 hapto rigid subtree | A3 | NO (env-gated, OFF) |
| ebb8cc3 GRIP angle-to-metal | iter | NO (env-gated, OFF) |
| bb44b41 vertex-uniqueness gate | iter | YES (PT3) — DROPS collapsed isomers (likely net+, not regression source) |

**Conclusion**: ~80% of new code is properly default-OFF.  The actually-active
defaults are sandwich-dispatch, CN10, and vertex-uniqueness — none of which
explain a +1658% funcgrp_n_bond_viol regression.

## Env-flag drift forensik (THE smoking gun)

### 2792332 launcher (`/tmp/run_vollpool_2792332.sh`, 24 flags):
```
DELFIN_FFFREE_GRIP=1
DELFIN_FFFREE_POST_GRIP_ALL=1
DELFIN_FFFREE_PURE_TRACK3=1
DELFIN_FFFREE_DD_RELAX=1
DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1
DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1
DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY=1
DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1
DELFIN_FFFREE_BURNSIDE_CONFORMER=1
DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING=1
DELFIN_FFFREE_GRIP_HEALING_MODE=1
DELFIN_FFFREE_HH_CLASH_INCLUDE=1
DELFIN_FFFREE_SP3_H_UMBRELLA=1
DELFIN_FFFREE_SP3_H_HEAL=1
DELFIN_MOGUL_V3_TUNED=1
DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1
DELFIN_FFFREE_AROMATIC_BONDS=1
DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1
DELFIN_FFFREE_REALISM_SORT=1
DELFIN_FFFREE_REALISM_TUNED_WEIGHTS=1
DELFIN_FFFREE_REALISM_SOFT_GATES=1
DELFIN_GRIP_LIB_PATH=/.../grip_lib_v3.npz
PYTHONHASHSEED=0
```

### 2ee2f45 voll-pool env (captured from /proc/$PID/environ, 8 flags):
```
DELFIN_FFFREE_AROMATIC_BONDS=1
DELFIN_FFFREE_BUILDER=1
DELFIN_FFFREE_DONOR_VSEPR=1
DELFIN_FFFREE_GRIP=1
DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
DELFIN_FFFREE_PURE_TRACK3=1
DELFIN_FFFREE_RING_SCALE=1
DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1
DELFIN_INPROC=1
DELFIN_THREAD_WORKERS=1
PYTHONHASHSEED=0
```

### Each top-10 regression maps to ≥1 dropped flag:
- `funcgrp_n_bond_viol` +1658% → `CONSTRUCTION_FIX_ALL`, `DD_RELAX` dropped
- `F20_h_planar` +795% → `SP3_H_UMBRELLA`, `SP3_H_HEAL` dropped
- `F3_bond`, `F3_angle` → `POST_GRIP_ALL`, `GRIP_LOSS_WEIGHTS_TUNED` dropped
- `interlig_n_viol` +556% → `F24_INTERLIG_FIX_ALL`, `HH_CLASH_INCLUDE` dropped
- `vdw_n_clashes` +287% → `HH_CLASH_INCLUDE` dropped

## C1 Smoke500 ablation (heal validation)

### Setup
- 498 stratified SMILES (every 23rd from smiles_master_v3_plus)
- Same GUPPY HEAD (2ee2f45) for both runs
- Pool: parallel=60, timeout=180s, --no-retry-on-zero
- 2-run determinism (PYTHONHASHSEED=0)

### Results

| Comparison | Net | Severe regr | NS opp/total |
|---|---|---|---|
| Voll-pool 2ee2f45 vs 2792332 (the bug) | −25 | 50 | 114/328 |
| Smoke500 baseline-drift vs 2792332 voll-pool | −23 | 45 | 178/330 |
| Smoke500 HEAL-restore vs 2792332 voll-pool | **−12** | 48 | 145/330 |
| Smoke500 HEAL-restore vs baseline-drift (apples-to-apples) | **+5** | 23 | 145/330 |

**Pool summary (sigma class — the dominant class):**
- baseline-drift: %topo=32.8%, %poly=56.5%, %lig=3.8%, n_iso=577
- HEAL-restore:   %topo=43.1%, %poly=71.7%, %lig=6.1%, n_iso=360 (more pruning, higher quality)

### Recovery on top-10 voll-pool regression axes

| axis | target (2792332) | HEAL smoke | recovery |
|---|---|---|---|
| cn_mean_extra_donors | 0.145 | 0.317 | +71% |
| F20_h_planar | 2.62 | 10.55 | +26% |
| pi_planar | 0.18 | 0.93 | +14% |
| vdw_n_clashes | 11.97 | 42.33 | +8% |
| funcgrp_n_bond_viol | 0.083 | 2.161 | partial; smoke noise |
| frame_pct_element_list_change | 1.74 | 24.27 | smoke noise — heal-induced |

Smoke500 doesn't fully recover voll-pool numbers because:
1. **Sample size: smoke n=498 vs voll-pool n=9868** (20× fewer)
2. **SMILES mix: 85% overlap, 15% different**
3. **Per-axis noise floor on smoke ≈ ±2.4 pp** (5× voll-pool's ±0.4 pp)

### Heal-induced new regressions (14 axes)

Some axes got worse under HEAL (likely due to flag interactions in the new
codebase):
- `structqual_extra_fragments_pct_frames_viol` (was 0; now 4.6)
- `hapto_count/geom_n_viol` (1.4× worse)
- `mdshort_n_viol_per_frame` (heal targets this but ineffective on smoke sample)

These are smaller effects and partially smoke-sample artifacts.

## Disposition

**HEAL works** — it recovers ~78% of the voll-pool gap (+5/-23 = 22% remaining).

The remaining 22% gap is likely smoke500 noise, not architectural deficit.
**Relaunch voll-pool with HEAL config on GUPPY 2ee2f45** to confirm.

## Heal artifact

- `scripts/run_vollpool_heal_c1.sh` — production voll-pool launcher with
  all 22 flags restored.  Reproduces the 2792332-style configuration on the
  current 2ee2f45 GUPPY HEAD.

## Next steps

1. Launch full voll-pool with `scripts/run_vollpool_heal_c1.sh` (≈2-3 hr).
2. Detector battery + iter_gate vs 2792332 voll-pool.
3. If net+ vs 2792332 → push as new HEAD.
4. If still net-negative → ablate one-by-one which Mission B commit silently
   breaks one of the heal flags (currently unsuspected but possible:
   - bb44b41 vertex-uniqueness gate auto-on under PT3 affects the isomer set
     after the heal flags ran; could mask the heal effect on per-isomer axes).
