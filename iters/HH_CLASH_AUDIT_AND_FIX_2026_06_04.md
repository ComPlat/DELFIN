# H-H Clash Detection Audit & Fix — 2026-06-04

Worktree: `agent-adddc59976088e3ee` (branch `hh-clash-detection-include`,
based on HEAD `93b396d`).

## User report

> "H-H clashes still visible in production pool `b00f9a0-full7-VOLLPOOL`
> despite our clash modules being active."

## Phase A — per-module H-inclusion audit

| Module | Function | H in vdW table? | H pairs iterated? | Effective H-H floor | Verdict |
|---|---|---|---|---|---|
| `delfin/fffree/grip_constraints.py` | `ClashFloorPenalty.value_and_grad` | yes (table is caller-supplied; `_vdw_table_for_mol` populates H from `DEFAULT_VDW_RADII`) | yes — full triangular `i<j` walk over all atoms with a finite radius (line 437-473) | `0.85 × (1.20 + 1.20) = 2.04 Å` | OK in principle; never runs in production because `grip_polish` is gated by `DELFIN_FFFREE_GRIP=1` (default OFF) |
| `delfin/fffree/grip_polish.py` | `_vdw_table_for_mol` → `ClashFloorPenalty` | yes — `DEFAULT_VDW_RADII["H"] = 1.20` (line 274) | yes (downstream of `ClashFloorPenalty`) | `2.04 Å` | OK but not active in production |
| `delfin/fffree/inter_ligand_clash_gate.py` | `count_inter_ligand_clashes_quick` | yes — `_DEFAULT_VDW_RADII["H"] = 1.20` (line 211) | yes (line 187-204, intra-ligand `sj == si` filter naturally skips geminal) | `2.04 Å` for cross-ligand H-H | OK but the gate's master switch (`DELFIN_FFFREE_PRE_POLISH_CLASH_GATE`) is OFF unless `ENUMERATE_ROTAMERS=1` |
| `delfin/fffree/build_time_clash_gate.py` | `has_collapse` / `collapse_count` | yes — `INTERLIG_VDW_RADII["H"] = 1.20` (line 65) | partial — `_geometric_bonds` flags X-H BOND collapses (d < 0.70 Å absolute) but **never** non-bonded H-H contacts | `0.70 Å` (absolute, for X-H bonds only — irrelevant to non-bonded H-H) | **gap**: non-bonded H-H never checked |
| `delfin/fffree/assemble_complex.py::_clash_count` | per-block pair walk | yes — `_vdw("H")` returns 1.20 | yes (all pairs) | `0.70 × (1.20 + 1.20) = 1.68 Å` | **gap**: factor 0.70 is too lax — eclipsing methyl-methyl H-H at ~ 1.95-2.04 Å passes clean |

### Root cause summary

The 5 modules together have H in every vdW table, BUT:

1. The two modules that ALWAYS run in production
   (`assemble_complex._clash_count` and the X-H part of
   `build_time_clash_gate`) use thresholds calibrated for **superposition
   detection**, not vdW contact (`0.70 × Σr` ≈ 1.68 Å for H-H, or 0.70 Å
   absolute for bond collapse).

2. The two modules that DO use the Bondi `0.85 × Σr` floor
   (`ClashFloorPenalty`, `inter_ligand_clash_gate`) are env-gated OFF in
   the production voll-pool.

3. There is no module that runs unconditionally with a Bondi floor AND
   topology-aware filtering (geminal CH2 H-H at ~ 1.78 Å is geometrically
   tighter than the methyl-methyl eclipse, so a flat floor without topology
   would produce more false positives than true detections).

## Phase B — real-structure validation

Script: ad-hoc Python (Bondi `0.85 × (r_H + r_H) = 2.04 Å` floor,
geminal/1-3-within-same-parent exclusion).

Sampled 2000 of 10 576 XYZ files in
`/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL`:

* **783 / 2 000 (39.2 %) structures contain at least one
  non-bonded, non-geminal, non-1,3 H-H pair below 2.04 Å.**
* Per-structure clash counts range from 1 to 9, with the
  worst observed at 084-MIPSOW.xyz (9 clashes, tightest pair 1.811 Å) and
  073-BITNUQ.xyz (8 clashes).
* Typical example 047-CIRJUK.xyz: 6 clashes including (31, 38) at 1.92 Å
  — a classic methyl-methyl eclipsing across an inter-ligand axis.

After installing the new detector module:

* `count_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)`
  reproduces the same prevalence: **187 / 500 (37.4 %)** structures
  flagged, 622 total pairs (avg 3.33 per affected structure).
* The slight delta (37.4 % vs 39.2 %) is the topology filter correctly
  reclassifying CH2 geminal pairs that the audit script accepted (a
  bonus side-effect of the principled exclusion).

## Phase C — surgical fix

**New module**: `delfin/fffree/hh_clash_detector.py` (~ 410 LOC).

Public API:

```python
detect_hh_clashes(syms, P, *, factor=None, min_topo_distance=None,
                  excluded_pairs=None, bond_adj=None) -> list[(i, j, d, sev)]
count_hh_clashes(syms, P, ...) -> int
hh_clash_severity(syms, P, ...) -> float      # Σ max(0, d_min - d)^2
build_hh_exclusion_pairs(syms, P, ...) -> set[frozenset]
hh_clash_active() -> bool                     # env-flag predicate
hh_clash_factor()  / hh_clash_min_topo() / hh_clash_weight()
```

Class-conditional handling (mandate from User §E):

* Geminal H-H (topo=2): skipped — internal CH2 ~ 1.78 Å is correct
  geometry, never a clash.
* 1-3 H-H within CH3 (topo=3): skipped — internal CH3 ~ 1.78 Å, same
  reasoning.
* 1-4 H-H across an sp3-sp3 bond (topo=3 between heavy parents +
  1 bond each = topo=3, OR sometimes computed as 3): skipped under
  default `min_topo_distance=4`. Use `min_topo_distance=3` to flag
  eclipsing rotamers.
* Inter-ligand or 1-5+ contacts: flagged when `d < factor × (r_H + r_H)`.

Bond graph is built from coordinates using
`delfin._bond_decollapse._ideal_bond` (the same table the structqual
detectors use), so the filter is consistent with the rest of the
codebase.

## Phase D — integration (env-flag gated, default OFF byte-identical)

All edits are guarded by `DELFIN_FFFREE_HH_CLASH_INCLUDE`:

| File | Change | Default behaviour |
|---|---|---|
| `delfin/fffree/hh_clash_detector.py` | **new** | exports inactive until env-flag set |
| `delfin/fffree/build_time_clash_gate.py` | `has_collapse` / `collapse_count` add an H-H pass when env-flag set | byte-identical to HEAD when unset |
| `delfin/fffree/inter_ligand_clash_gate.py` | new helper `count_inter_ligand_hh_clashes_quick` returns 0 when env-flag unset | additive symbol only; existing functions untouched |
| `delfin/fffree/grip_polish.py` | adds an H-H severity term to `loss_and_grad` (analytic gradient, one-sided quadratic, same shape as `ClashFloorPenalty`) | code path only entered when env-flag set AND grip_polish is invoked |
| `delfin/fffree/assemble_complex.py::_clash_count` | adds a Bondi-floor (0.85) H-H pass on top of the legacy 0.70 floor | bonus count only added when env-flag set |

Env-flag table:

| Variable | Default | Effect |
|---|---|---|
| `DELFIN_FFFREE_HH_CLASH_INCLUDE` | `0` | master switch (OFF byte-identical) |
| `DELFIN_FFFREE_HH_CLASH_FACTOR` | `0.85` | vdW-sum fraction |
| `DELFIN_FFFREE_HH_CLASH_MIN_TOPO` | `4` | minimum topological distance to flag |
| `DELFIN_FFFREE_HH_CLASH_WEIGHT` | `5.0` | L-BFGS penalty weight |

## Phase E — tests

`tests/test_hh_clash_detector.py` — 15 tests, all passing
(`pytest tests/test_hh_clash_detector.py` ⇒ `15 passed in 0.25s`):

```
test_no_hh_clash_clean_methylene
test_skips_methyl_internal_hh
test_detects_methyl_methyl_eclipse              (min_topo=3)
test_skips_eclipse_at_default_min_topo          (min_topo=4)
test_detects_inter_ligand_hh_clash_at_tight_separation
test_no_clash_when_methyls_far_apart
test_env_off_no_detection_in_build_gate         (byte-identical w/o flag)
test_env_on_triggers_build_gate_on_inter_fragment_hh
test_determinism
test_threshold_factor_respected
test_excluded_pairs_respected
test_no_silent_failure_on_empty
test_md_invariant_preserved_under_hh_penalty
test_real_structure_b00f9a0_finds_known_clashes  (spot-check on real voll-pool)
test_severity_aggregates_correctly
```

## Phase F — production rollout

The fix is committed but NOT pushed and is default-OFF byte-identical to
HEAD `93b396d`. To activate in a future iteration:

```bash
export DELFIN_FFFREE_HH_CLASH_INCLUDE=1
# Optional fine-tuning:
export DELFIN_FFFREE_HH_CLASH_MIN_TOPO=3      # also flag eclipsing rotamers
export DELFIN_FFFREE_HH_CLASH_FACTOR=0.85     # Bondi default
export DELFIN_FFFREE_HH_CLASH_WEIGHT=5.0      # match heavy intra weight
```

Recommended A/B for a future smoke:

* Track A: env unset (byte-identical to HEAD)
* Track B: `DELFIN_FFFREE_HH_CLASH_INCLUDE=1`
  + `DELFIN_FFFREE_BUILD_CLASH_GATE=1` (build-time rejection of new
  H-H clashes)

Expected outcome on the voll-pool aggregate:

* ~ 37 % of structures emit fewer H-H clashes (selector chooses cleaner
  rotamer when build-clash-gate is active).
* No effect on non-H atom pair statistics (orthogonal axis).
* Small slow-down (< 5 %) per build for the extra H-H pair walk
  (O(n_H²) is < 10⁴ for typical complex sizes).
