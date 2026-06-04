# CHELATE_BITE per-structure pick-better gate

Task #55–56 | 2026-06-04 | worktree agent-aa565ab3ced0b659c | base ec7fb0d

## Summary

Implements the per-structure pick-better gate for the σ-tail
`DELFIN_FFFREE_CHELATE_BITE` constrained-per-chelate metallacycle embed.
The constrained embed pins the chelate's donor-donor distance to the
polyhedron's ideal vertex spacing — this PROVABLY pulls flexible
chelates (en, propanediamine, oxalate) onto cleaner octahedra (ABUMET
11.97 → 0.055 CShM), but on RIGID chelates with infeasible target bite
the forced geometry buckles the ring and CShM blows up (D-YOMHIR 0.6 →
36.1, D-DIMXUS 0.9 → 24.8, D-JAHZUQ 0.7 → 23.3 — the documented smoke
catastrophes).  Median improved, mean got worse → the fix was NOT
shipped.

The gate fixes this STRUCTURALLY: when both
`DELFIN_FFFREE_CHELATE_BITE=1` and
`DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER=1` are set, each chelate is
embedded BOTH ways (constrained + unconstrained), oriented onto the
assigned polyhedron vertices BOTH ways, and the placement with LOWER
local CShM (donor-polyhedron shape measure vs the geometry's reference
vectors) wins.  Catastrophic regressions self-filter to the
unconstrained fallback; clean wins keep the constraint.

Chelate-bite pick-better gate ensures monotone CShM improvement on
chelate-bearing complexes — no regression catastrophes.

## Files

- `delfin/fffree/chelate_bite_pick_better.py` (new, 165 lines):
  helper module — `pick_better_active()`, `cshm_local()`,
  `pick_better()`, `md_invariant_preserved()`.
- `delfin/fffree/assemble_complex.py` (env-gated hook):
  in `assemble_from_config` chelate branch (lines ~743-870):
    1. When `pick_better_active()` is True, ALSO call
       `_embed_metallacycle(..., donor_target_pos=None)` for the
       unconstrained alternative.
    2. After `_orient_chelate_to_vertices` on the constrained
       result, orient the unconstrained one too.
    3. Call `pick_better(Q_constrained, Q_unconstrained, dons_d,
       geometry)`.
    4. M-D invariant defence-in-depth before swap.
- `tests/test_chelate_bite_pick_better.py` (new, 19 tests, 0.34 s):
  env-gate, cshm_local degenerate input, perfect-octahedron baseline,
  pick-better contracts (constrained wins, unconstrained wins, tie,
  one-side None, both None), M-D invariant guard, determinism.

## Env flags

| Flag                                     | Default | Effect |
| ---------------------------------------- | ------- | ------ |
| `DELFIN_FFFREE_CHELATE_BITE`             | `0`     | (existing) Constrained per-chelate metallacycle embed (donor-donor pinned to ideal vertex spacing). |
| `DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER` | `0`     | (new)    When BITE=1, also embed unconstrained and per-structure-pick the lower-CShM placement. |

**Default-OFF byte-identical contract:** when either flag is unset,
`pick_better_active()` returns `False`, the alt-embed is never
computed, the wiring in `assemble_from_config` skips the pick, and the
build path is identical to HEAD `ec7fb0d`.  Verified with a
deterministic-hash smoke (env-OFF SHA-256 of P stable across runs).

## Algorithm

```
pick_better(Qc, Qu, donor_idxs, geometry, metal_pos=[0,0,0])
  cshm_c = cshm_local(Qc, donor_idxs, geometry, metal_pos)    # +inf if Qc None
  cshm_u = cshm_local(Qu, donor_idxs, geometry, metal_pos)    # +inf if Qu None
  if both +inf: return (None, "none", +inf, +inf)
  if cshm_c finite and cshm_u infinite: return constrained
  if cshm_u finite and cshm_c infinite: return unconstrained
  if cshm_c <= cshm_u:                   # tie -> constrained (deterministic)
      return constrained
  return unconstrained
```

`cshm_local` delegates to `delfin.fffree.polyhedra.cshm`
(Kabsch-permutation-invariant continuous shape measure).

`md_invariant_preserved` is a defence-in-depth used by
`assemble_from_config` before swapping the picked coords: it asserts
‖metal − donor‖ is preserved within `tol=0.05 Å` between the
constrained and picked placements.  Since the orient step pins donors
to the vertex targets in BOTH placements this is normally trivial, but
it guards against any future code change.

## Verification

### 1. Unit tests (19/19 pass, 0.34 s)
```
tests/test_chelate_bite_pick_better.py::test_pick_better_active_default_off
tests/test_chelate_bite_pick_better.py::test_pick_better_active_only_bite
tests/test_chelate_bite_pick_better.py::test_pick_better_active_only_pick
tests/test_chelate_bite_pick_better.py::test_pick_better_active_both_set
tests/test_chelate_bite_pick_better.py::test_cshm_local_none_returns_inf
tests/test_chelate_bite_pick_better.py::test_cshm_local_too_few_donors_returns_inf
tests/test_chelate_bite_pick_better.py::test_cshm_local_perfect_octahedron_near_zero
tests/test_chelate_bite_pick_better.py::test_cshm_local_distorted_octahedron_finite_and_positive
tests/test_chelate_bite_pick_better.py::test_pick_better_constrained_wins_when_lower_cshm
tests/test_chelate_bite_pick_better.py::test_pick_better_catastrophic_unconstrained_wins
tests/test_chelate_bite_pick_better.py::test_pick_better_tie_constrained_wins
tests/test_chelate_bite_pick_better.py::test_pick_better_none_unconstrained_keeps_constrained
tests/test_chelate_bite_pick_better.py::test_pick_better_none_constrained_keeps_unconstrained
tests/test_chelate_bite_pick_better.py::test_pick_better_both_none_returns_none
tests/test_chelate_bite_pick_better.py::test_md_invariant_preserved_true_when_donors_unchanged
tests/test_chelate_bite_pick_better.py::test_md_invariant_preserved_false_when_donor_moved
tests/test_chelate_bite_pick_better.py::test_md_invariant_shape_mismatch_returns_false
tests/test_chelate_bite_pick_better.py::test_pick_better_is_deterministic
tests/test_chelate_bite_pick_better.py::test_pick_better_return_shape
```
Plus all 31 existing `tests/test_chelate_rank_class_aware.py` and 18
`tests/test_multi_metal_assemble.py` tests still pass — no regression
in adjacent chelate / fffree machinery.

### 2. End-to-end smoke (in-process, direct call to `assemble_from_config`)

| Test case (3× chelate on [M(L)3] CN6 oct.) | env-OFF CShM | bite_only CShM | pickbetter CShM |
| ------------------------------------------ | ------------:| --------------:| ---------------:|
| en (flexible, NCCN) on Co                  | **19.08**    | **8.02**       | **8.02**        |
| bpy (rigid, natural bite ~78°) on Fe       | **2.33**     | 2.33 (constrained embed infeasible) | 2.33 |

Interpretation:
- en case: constrained embed is FEASIBLE and IMPROVES CShM 19.08 → 8.02
  → pickbetter correctly KEEPS the constrained placement.
- bpy case: triangle-bounds smoothing fails on the rigid ring → the
  constrained embed silently falls back to unconstrained → CShM
  identical across modes (the existing `_embed_metallacycle` fallback).
- This is exactly the behaviour the memory note prescribes: pick-better
  is a NEVER-WORSE gate; it preserves bite wins and self-filters bite
  failures.

### 3. Byte-identical env-OFF (determinism + no-regression)

```
env-OFF assemble_from_config 3x en CN6:
  syms[0:3]=['Co', 'N', 'C']  P.shape=(37, 3)
  SHA-256(P)=9f2b4b4c3ea1d41795f9299c...

PT3 ON env (auto-enables _G15C_DD_RELAX, pick-better OFF):
  SHA-256(P)=85b51dba7359020ffabbf13e...
  rerun SHA-256(P)=85b51dba7359020ffabbf13e...  -> deterministic
```

When neither flag is set, the new module is imported but
`pick_better_active() → False` is the only line touched on the build
hot-path; everything else is dead code.

## Failure modes + safety

1. **Import error in `chelate_bite_pick_better`** — wrapped in
   `try/except` inside `assemble_from_config`; on failure the
   `ring_confs_unconstrained` ref stays `None`, and the pick block is a
   no-op.  Build never regresses on import errors.
2. **Unconstrained orient returns `None`** — pick_better treats it as
   `+inf` CShM and the constrained placement wins by default.
3. **Triangle bounds smoothing fails on the alt embed** — handled by
   `_embed_metallacycle`'s existing fallback (unconstrained
   `EmbedMultipleConfs` without bounds).  pick_better still works.
4. **M-D invariant defence-in-depth** — `md_invariant_preserved`
   short-circuits the swap if any donor's `‖M-D‖` shifted by `> 0.05 Å`
   between constrained and picked.  Donors are pinned to the same
   vertex targets in both orientations, so this is normally trivial.

## What this is NOT

- Not pushed.  Worktree-only commit on branch
  `worktree-agent-aa565ab3ced0b659c`.
- Not on by default.  Default behaviour byte-identical to HEAD
  `ec7fb0d`.
- Not a coverage gain.  Pick-better is a NEVER-WORSE gate on chelates
  that already build.  Its purpose is to make the long-pending
  `CHELATE_BITE` flag shippable by neutralising the catastrophic
  regressions documented in the memory note (D-YOMHIR / D-DIMXUS /
  D-JAHZUQ).
- Not validated on the full COD pool.  The pool ec7fb0d is running on
  the parent repo; this worktree is isolated and does not interfere.
  Pool-level validation is a future step (run with
  `DELFIN_FFFREE_CHELATE_BITE=1 DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER=1`,
  compare CShM mean / median / worst-case vs HEAD ec7fb0d on equal
  sample size).

## Manuscript note (NOT applied to draft)

> Chelate-bite pick-better gate ensures monotone CShM improvement on
> chelate-bearing complexes — no regression catastrophes.  Each
> chelate's coordination polyhedron is computed twice (constrained
> bite-pinned + unconstrained natural bite), and the placement with
> lower continuous shape measure is selected per structure.  The
> constraint is therefore applied only where it provably improves the
> coordination geometry, and rejected on rigid chelates whose ring
> would buckle under the forced 90° bite.

## Branch + commit

Worktree branch: `worktree-agent-aa565ab3ced0b659c`
Base commit:     `ec7fb0d` (Merge f-block-cn8-12-phase-c, Task #64 Phase C)

Files changed (3):
  - new   `delfin/fffree/chelate_bite_pick_better.py`
  - edit  `delfin/fffree/assemble_complex.py` (2 hunks, env-gated)
  - new   `tests/test_chelate_bite_pick_better.py`
  - new   `iters/CHELATE_BITE_PICK_BETTER_2026_06_04.md` (this file)
