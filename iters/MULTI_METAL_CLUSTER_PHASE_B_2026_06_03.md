# Multi-Metal / Cluster construction path (Task #63 Phase B)

Date: 2026-06-03
Branch: `multi-metal-cluster-phase-b` (based on race-integration HEAD `3841b75`)
Env-flag: `DELFIN_FFFREE_MULTI_METAL=1` (also auto-enabled by `DELFIN_FFFREE_PURE_TRACK3=1`)
Default OFF: when flag unset the dispatch is bypassed -> single-metal path is byte-identical to `3841b75`.

## Phase A — design decision: build on prototype

After reading the existing assets:

| Existing asset | Lines | Status |
|---|---|---|
| `delfin/_cluster_scaffold.py` | 844 | PROTOTYPE T7.3 — strong skeleton enumeration (triangle/linear/tet/butterfly/oct/star), donor partition (terminal/edge/face), placement math, hapto skeleton. Not wired. |
| `delfin/fffree/multi_metal_polyhedra.py` | 224 | Clean vertex-math for M2/M3/M4-tet/M4-sq/M6-oct, M-M empirical bond table + covalent fallback. |
| `delfin/fffree/bridging_ligand.py` | 149 | μ²/μ³/μ⁴ detection on `mol` (≥2 metal bonds) + edge/face/apex placement math. |
| `delfin/fffree/decompose.py` | — | Already handles multi-metal by picking primary; doesn't construct the cluster. |
| `delfin/fffree/assemble_complex.py` | 1165 | Single-metal `assemble_from_config`; entry point. |
| `delfin/smiles_converter.py::_build_multimetal_scaffold` | — | Legacy dimer-only path inside the non-fffree converter. |

**Decision: BUILD ON.** Reasons:
* `multi_metal_polyhedra.py` + `bridging_ligand.py` + `_cluster_scaffold.py` together cover ~90% of the math; what is missing is the *orchestration*.
* Rewriting would re-invent the empirical M-M tables, the polyhedron formulas and the donor-partition logic without gain.
* Production-safety: leaving the prototype untouched protects the running race-pool.

What is added:
* `delfin/fffree/multi_metal_assemble.py` (new orchestrator, ~470 LOC)
* `delfin/fffree/cluster_scaffold_polyhedra.py` (new facade wrapping `_cluster_scaffold` PROTOTYPE for use from the fffree namespace)
* `delfin/fffree/assemble_complex.py` (env-gated dispatch only, ~22 LOC added at top of `assemble_from_config`)
* `tests/test_multi_metal_assemble.py` (~12 tests)

## Phase B — env-gated dispatcher in `assemble_complex.py`

`assemble_from_config(metal, geometry, config, ligands, refine)` checks:

1. `DELFIN_FFFREE_MULTI_METAL=1` (or `DELFIN_FFFREE_PURE_TRACK3=1`)
2. AND `len(ligands_containing_metal_atoms) >= 1` OR a `multi_metal_spec` is attached to `config`.

When both true and a multi-metal spec is detected, dispatch to `multi_metal_assemble.assemble_multi_metal(...)`; otherwise fall through to the current single-metal path **byte-identical**.

The dispatcher is read-only on input objects and has a safe try/except — any failure falls back to the legacy single-metal path; the caller observes no change.

## Phase C — `multi_metal_assemble.py` architecture

```
assemble_multi_metal(mol, ligands, metals, config) -> (syms, P, donors) | None
        │
        ├── 1. Build metal-connectivity graph from `mol` (M-M bonds + μ-bridged pairs)
        │
        ├── 2. Pick metal-skeleton shape:
        │       n=2 -> linear dimer       (multi_metal_polyhedra.m2_dimer)
        │       n=3 -> triangle / linear  (multi_metal_polyhedra.m3_triangle / _cluster_scaffold._skeleton_linear_m3)
        │       n=4 -> tetrahedron / square (multi_metal_polyhedra.m4_tetrahedron / m4_square_planar / _cluster_scaffold._skeleton_butterfly_m4)
        │       n=5 -> star_m5            (_cluster_scaffold._skeleton_star_m5)
        │       n=6 -> octahedron         (multi_metal_polyhedra.m6_octahedron)
        │       n>6 -> generic (centroid + lex-sorted)
        │       Choice is graph-driven: declared M-M-edges (full / partial)
        │
        ├── 3. Lex-sort metal indices for determinism; map metal_idx -> skeleton vertex
        │
        ├── 4. Place metals at scaled vertex coords (mm_distance from empirical table)
        │
        ├── 5. Place μ-bridging atoms (detect via bridging_ligand.detect_bridging_ligands):
        │      μ² -> edge midpoint  + perpendicular offset (M-X chemistry)
        │      μ³ -> face centroid  + face-normal offset
        │      μ⁴ -> centroid of 4 metals
        │      M-X distance = covalent-radius sum (clamped to [1.7, 2.6] Å)
        │
        ├── 6. Per-metal coordination sphere for non-bridging donors:
        │      Each metal still gets its own polyhedron orientation (radial from cluster centroid).
        │      Terminal-donor radial axis = (M_pos - cluster_centroid).
        │      Reuse `_cluster_scaffold.place_skeleton_with_ligands` for the terminal layer.
        │
        ├── 7. Hard-rollback gates:
        │      - per-bond M-D distance ±0.05 Å of target (single-metal contract)
        │      - per M-M ±0.10 Å of target (multi-metal contract — looser)
        │      - failure -> return None, caller falls back to single-metal
        │
        └── 8. Return (syms, P, donors) tuple identical in shape to `assemble_from_config`.
```

Determinism: PYTHONHASHSEED=0, lex-sort metal indices, lex-sort donors per metal, lex-sort skeleton edges.

## Phase D — tests `tests/test_multi_metal_assemble.py`

| # | Test | What it asserts |
|---|---|---|
| 1 | `test_env_off_byte_identical` | Without env var, dispatcher returns `None` (legacy fall-through). |
| 2 | `test_metal_connectivity_graph_dimer` | `[Cu](mu-Cl)Cl]2` -> 2 metals, 1 M-M edge implied via Cl bridge. |
| 3 | `test_skeleton_choice_dimer` | n=2 -> M2_dimer vertices, correct mm_distance. |
| 4 | `test_skeleton_choice_triangle` | n=3 + 3 M-M edges -> triangle. |
| 5 | `test_skeleton_choice_tetrahedron` | n=4 + ≥6 M-M edges -> tetrahedron. |
| 6 | `test_skeleton_choice_octahedron` | n=6 + 12 M-M edges -> octahedron. |
| 7 | `test_mu2_bridge_geometry` | μ²-bridge atom equidistant from both metals, ±0.05 Å of target. |
| 8 | `test_mu3_face_cap` | μ³-cap above face centroid, equidistant from 3 metals. |
| 9 | `test_mm_distance_lookup` | `m_m_distance('Mo','Mo')` -> 2.20 Å (paddlewheel). |
| 10 | `test_determinism_same_seed_same_output` | Two runs of `assemble_multi_metal` produce identical coords. |
| 11 | `test_dispatch_falls_back_single_metal` | Single-metal SMILES never enters the multi-metal path. |
| 12 | `test_hard_rollback_on_violation` | Forced bad input -> returns None (caller safe). |

## Phase E — validation on test SMILES

20 multi-metal SMILES from `smiles_master_v3_plus.txt` (out of 134 found total).
Run with `PYTHONHASHSEED=0 DELFIN_FFFREE_MULTI_METAL=1`.

Results (see `iters/MULTI_METAL_PHASE_B_validation.json`):

|                              | count | note |
|---|---:|---|
| Total SMILES tested          | 20    | |
| Skeleton chosen              | 20    | every multi-metal mol matched a polyhedron |
| Built (passed hard-rollback) | 18    | 90% pass rate |
| Hard-rollback rejected       |  2    | safe fall-through to single-metal path |

Skeleton distribution:
* M2_dimer       : 17
* M3_linear      :  1
* M4_butterfly   :  1
* M6_octahedron  :  1

M-M distance recovery: every built mol has observed M-M = target within tolerance
(e.g. dimers: 2.51, 2.812, 2.764 A matching their empirical COD-p50 values exactly).

The two hard-rollback cases are heteronuclear pairs (Sc+Fe, Sn+Co) whose
empirical pair-distances do not match the regular-polyhedron geometry the
prototype uses for n>=3.  This is **correct** behaviour: the orchestrator
refuses to ship a mis-scaled M-M, the caller silently falls back to the
single-metal pipeline, no regression.

## Production safety

- Working tree limited to this worktree (`agent-afc9b0709e1b1be4c`).
- NO push, NO edits outside.
- All new code is behind `DELFIN_FFFREE_MULTI_METAL=1` — when unset, every code-path is identical to `3841b75`.
- Hard-rollback: any deviation from M-D ±0.05 Å or M-M ±0.10 Å -> caller falls back to single-metal.
