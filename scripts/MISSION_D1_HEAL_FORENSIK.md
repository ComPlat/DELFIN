# Mission D1 — Voll-Pool Regression HEAL Forensik (2026-06-05)

## TL;DR

The 2792332 → 2ee2f45 voll-pool regression (net −26, 51 severe axes) was
**NOT** caused by env-flag drift (Mission C1's hypothesis was wrong).  It was
caused by **commit 1b95142 (Mission A5)** which made
`DELFIN_FFFREE_PURE_TRACK3=1` auto-imply the `NO_FALLBACK` semantic.  All
production voll-pools running PT3 with the 2792332 22-flag stack thereafter
lost 81% of XYZ emissions to silent zero-isomer returns: the fffree backend
fails on ~80% of metal-SMILES (per A5's own forensik), and PT3 alone now
blocks the legacy fallback that would have rescued them.

The C1 heal (which restored the 22-flag stack on 8654d8f) made no
difference because the regression is in the dispatch code itself, not in
the launcher.

The D1 heal decouples the dispatch gate (BUILDER) from the post-fail block
gate (NO_FALLBACK), restoring the 2792332-byte-equivalent behaviour:
`PT3=1 alone → fall back to legacy on fffree fail`.

## Forensik (bisect)

Smoke10 of `head -10 smiles_master_v3_plus.txt` with the EXACT 2792332
22-flag environment, varying only the GUPPY HEAD commit:

| commit | non-empty XYZ / 10 | verdict |
|---|---|---|
| 2792332 (baseline) | 6 | GOOD |
| bb44b41 (vertex-uniqueness) | 6 | GOOD |
| 9e96246 (hapto rigid subtree) | 6 | GOOD |
| **1b95142 (A5 NO_FALLBACK auto-imply)** | **1** | **BAD ← culprit** |
| 0508efc (sandwich-piano) | 1 | BAD (inherits A5) |
| 8bba47f (oxalate) | 1 | BAD |
| 8654d8f (current HEAD, pre-heal) | 1 | BAD |
| 8654d8f + D1 heal | 7 | GOOD (matches/exceeds 2792332) |

`9e96246` → `1b95142` is the failure transition.  Inspecting the diff of
`1b95142` shows the entire change is the dispatch logic in
`delfin/smiles_converter.py` line 26342:

```python
# A5 (BROKE production)
_no_fallback = (
    _delfin_env_int("DELFIN_FFFREE_PURE_TRACK3", 0)
    or _delfin_env_int("DELFIN_FFFREE_NO_FALLBACK", 0)
)
_fffree_dispatch = (
    _delfin_env_int("DELFIN_FFFREE_BUILDER", 0) or _no_fallback
)
```

Whenever `_no_fallback` is true and fffree fails, the function returns
`[], None` — no isomers, no error.  PT3=1 in production voll-pools
silently activates this.

## Voll-pool A1 forensik (per-SMILES failure)

```
results/smiles_master_v3_plus__DELFIN__20260605_183956.jsonl
                                       (heal-c1, 8654d8f, broken)
  total              11452 SMILES
  zero-isomer        9296  (81.2%)
  with error         0     <-- silent, no error captured
  retried            0
  recovered          0

results/smiles_master_v3_plus__DELFIN__20260604_220901.jsonl
                                       (2792332, baseline, healthy)
  total              11456 SMILES
  zero-isomer        1670  (14.6%)   <-- mostly timeouts
  with error         1670
    1587 timeout
      62 RuntimeError
       8 All strategies failed
       8 Failed to generate 3D coordinates for SMILES
       3 exit 127
       1 Error converting SMILES to XYZ
       1 exit -6
  retried            88
  recovered          20
```

XYZ archive size verification:

```
2792332-aromatic-symmetry-VOLLPOOL/   non-empty 9787 / 9868  (99.2 % usable)
2ee2f45-heal-c1-VOLLPOOL/             non-empty 2156 / 11456 (18.8 % usable)
```

The `n_files_total` in the C1 forensik report (11452) was misleading
because `benchmark_runner.py:693` **creates the XYZ file immediately on
entry** (`open("w")`), so a 0-isomer build still leaves a 0-byte file
behind.  `compare_pools.py` counted those.

## D1 heal

```python
# D1 (RESTORED production)
_no_fallback = _delfin_env_int("DELFIN_FFFREE_NO_FALLBACK", 0)
_fffree_dispatch = (
    _delfin_env_int("DELFIN_FFFREE_BUILDER", 0)
    or _delfin_env_int("DELFIN_FFFREE_PURE_TRACK3", 0)
    or _no_fallback
)
```

Behaviour matrix:

| flags set | dispatch fffree? | fall back to legacy on fail? |
|---|---|---|
| (none) | NO (skip block) | n/a |
| BUILDER=1 | YES | YES |
| PT3=1 | YES | YES (matches 2792332) |
| NO_FALLBACK=1 | YES | NO (strict) |
| PT3=1 + NO_FALLBACK=1 | YES | NO (strict, paper measurement mode) |

The contract A5's own commit message advocated for ("keep NO_FALLBACK
default OFF for production, set ON for honest fffree quality
measurements") is now actually implemented.

## Test status

`tests/test_no_fallback_pt3.py`: all 11 tests pass.  No assert change
needed; only the module docstring was updated.  The
`test_pure_track3_alone_implicitly_enables_dispatch` test continues to
pass because PT3=1 still triggers the dispatch — only the post-fail
return-`[]` is now gated on the explicit `NO_FALLBACK` flag.

## Smoke100 validation (8654d8f + D1 heal)

```
head -100 smiles_master_v3_plus.txt, 22-flag stack, parallel=32
non-empty XYZ:  82 / 100  (82.0 %)
zero-byte XYZ:   4
no archive:     14  (timeouts / outright failures)
```

Compare: at the same HEAD without the heal, smoke10 gave 1/10 emission.
The D1 heal brings emission back into the 2792332 regime.

## Voll-pool launch

`scripts/run_vollpool_D1.sh D1-heal-VOLLPOOL` — exact 2792332 22-flag
stack on GUPPY HEAD with D1 heal applied.  INPROC=1 fast-path enabled.

Partial-pool emission (interim measurement at n~450 SMILES):

```
files emitted     :  ~450
non-empty XYZ     :  ~441  (97.3 %)
zero-byte XYZ     :  ~9    (2.0 %)
silent zero rate  :  < 5 %  (vs 81 % pre-heal)
```

The healed pipeline produces real XYZ at the same rate as 2792332's
baseline.  Per-SMILES build cost is what slows the pool — the same
3-4h voll-pool wall-clock that 2792332 took.

Full pool completion + iter_gate vs `2792332-aromatic-symmetry-VOLLPOOL`
pending.  Detectors run via:

```
scripts/MISSION_D1_HEAL_FORENSIK.md  (this file)
scripts/run_vollpool_D1.sh           (launcher)
/tmp/D1_post_pool_eval.sh            (battery + iter_gate driver)
```

## Why this should beat 2792332 (not just match)

2792332's voll-pool was built almost entirely on the legacy/UFF
fallback path (per A5 forensik on the same flag stack: only ~20 % of
metal-SMILES use fffree-native, the rest fall through to legacy).
The D1 heal restores that same legacy fallback, so on those SMILES
the structures will be byte-equivalent.

But the ~20 % that DO build natively now get the post-2792332
improvements that are auto-active under PT3:

| Mission B/A commit | improvement | active under PT3? |
|---|---|---|
| bb44b41 vertex-uniqueness | drops collapsed isomers | YES |
| 0508efc SANDWICH-10 + PIANO-STOOL-8 | new polyhedra | YES |
| 862070a CN10 polyhedra (BICAP/CSAP/SAP) | CN10 coverage | YES |
| 2ee2f45 sandwich-piano dispatch + cone topo | improved hapto | YES |
| 8bba47f oxalate planar projector | improved internals | YES (autogated) |

The net should be net+ vs 2792332 because the legacy-path SMILES are
preserved (no regression) and the fffree-native SMILES are improved.

## Author / Provenance

Heal branch: `mission-D1-heal` (pushed to origin).
Heal commit: `24849e3` — hmaximilian <hmaximilian496@gmail.com>.
No Co-Authored-By trailers.
