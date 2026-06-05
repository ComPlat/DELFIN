# Comprehensive CCDC Extraction + GRIP Lookup Fix (2026-06-05)

## Mission

User directive 2026-06-05: "alles bitte fixen und aufnehmen sprich jegliche
fragmente aus ccdc extrahieren" — fix the GRIP lookup *and* extract every
fragment from CCDC.

## Root cause of the legacy GRIP miss rate

`grip_lib_v3.npz` has 1,142,066 master keys. ~16k of them mention Iridium,
~20k mention Platinum, etc. Yet a query for ``lookup_bond("Ir","sp3","C","sp2")``
returned ``None``.

After tracing the chain in ``delfin/fffree/grip_mogul_lookup.py``:

* Lookup builds the key
  ``["Ir","sp3", [["C",1,-1,"sp2"]], []]`` (orientation 1) and
  ``["C","sp2", [["Ir",1,-1,"sp3"]], []]`` (orientation 2) with single
  neighbour.
* The library *does not* store metal-centred fragments at all (the build
  script's ``_infer_hyb`` returns ``"*"`` for metals and centres with
  ``hyb="*"`` are explicitly skipped in v3 ``centre_keys`` emission).
* The library *does* store C-centred fragments with Ir among their
  neighbours — but in keys like
  ``["C","sp3", [["C",..],["C",..],["C",..],["Ir",1,-1,"*"]], []]`` with the
  C atom's *full* neighbour list.
* No fallback level in v3's `_fallback_levels` reduces the neighbour list to
  ``[["Ir",..]]`` alone. The two never meet.

### Why a naive "neighbour-superset" aggregation is *wrong*

Pooling these C-centred keys and reporting their `bond_mu` *does not* give
the Ir-C distance — the stored `bond_mu` is the median over *every*
bond distance at that centre (C-C + C-H + C-Ir = ~1.47 Å, dominated by
the more numerous C-C/C-H bonds). The Ir-C bond is buried.

### The actual fix: pair-resolved tables (v4)

The library has to record bond statistics *per ordered pair of
(Z, hyb)* directly. We add three new tables to the npz, populated per
actual bond/angle/improper rather than aggregated over each centre's
full neighbour list:

* ``pair_bond_keys[]`` — key ``[Z_lo, hyb_lo, Z_hi, hyb_hi]``, sorted.
* ``triple_angle_keys[]`` — key ``[Z_c, hyb_c, Z_l, Z_r]`` with Z_l<=Z_r.
* ``improper_pair_keys[]`` — key ``[Z_c, hyb_c, sorted(Z_a,Z_b,Z_c)]``.

Each table carries ``mu``, ``sigma`` (1.4826 * MAD), ``n``. The keys
support hybridisation wildcards ``"*"`` so metals (always ``"*"``) match
correctly.

## Deliverables

### 1. Fixed lookup `delfin/fffree/grip_mogul_lookup.py`

Adds:

* Env flag ``DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK`` (default OFF —
  byte-identical to legacy when unset).
* New v4-aware fields on ``GripLibrary``: ``has_pair_tables``,
  ``n_pair_bond``, ``n_triple_angle``, ``n_improper_pair``.
* New methods ``_v4_lookup_bond``, ``_v4_lookup_angle``,
  ``_v4_lookup_improper`` — these consult the pair tables only and apply
  a 4-step hyb-wildcard fallback chain on each side.
* ``_tm_lookup_bond/angle/improper`` are wired to the v4 methods and
  return ``None`` against v1/v2/v3 libraries (the env flag is harmless
  on those — no chemically wrong injection).
* The legacy ``lookup_bond/angle/improper`` chain runs *first*; the v4
  fallback runs *only* when the legacy chain misses and the env flag is
  set.

### 2. New extraction script `scripts/grip_build_mogul_lib_v4_ccdc.py`

Same CCDC-native iteration as v3 (forkserver pool, EntryReader per
worker, identifier-sorted indices for determinism). Additions:

* Emits all v3 fragment-centred keys exactly as before (backward
  compatible).
* Adds per-bond, per-angle, per-improper emission to the new pair
  tables. Metal centres are included for impropers (the v3 path skipped
  them).
* Determinism preserved: PYTHONHASHSEED=0, np/random seed 42,
  identifier lex-sort, no RNG in aggregation, single-threaded BLAS.

### 3. New library `grip_lib_v4.npz`

Smoke build (500 CSD entries) — already validated end-to-end:

```
Pair bond keys: 196 (n>=5: 121)
Triple angle:   454 (n>=5: 283)
Improper pair:  127 (n>=5: 63)
Lookups now return chemically correct values:
  Ir-C(sp2):  2.003 A (real ~2.0 A)   [was None in v3]
  Pd-Cl:      2.308 A (real 2.30 A)   [unchanged, terminal-Cl path]
  Cu-N(sp3):  1.979 A (real ~2.0 A)   [was None in v3]
  Cu-Cl:      2.256 A (real ~2.25 A)  [was None in v3]
```

Full build (1,436,119 entries, 48 workers) is running in the
background; expected pair_bond keys ~10000-50000, triple_angle keys
~30000-100000.

### 4. Tests `tests/test_grip_lib_v4_pair.py` (8 new)

* ``test_v4_pair_tables_loaded`` — has_pair_tables True on v4
* ``test_tm_fallback_default_off`` — byte-identical legacy when env unset
* ``test_tm_fallback_on_returns_pair_bond`` — Ir-C hits v4 pair table
* ``test_tm_fallback_hyb_wildcard`` — wildcard fallback through hyb chain
* ``test_tm_fallback_min_n_filters`` — min_n contract honoured
* ``test_tm_fallback_angle_centered_at_metal`` — C-Ir-C angle hits
* ``test_tm_fallback_improper_metal_centre`` — improper at Ir hits
* ``test_v3_library_returns_none_with_flag`` — critical: v3 library does
  NOT inject chemically wrong values when the env flag is on

All 8 pass.

### 5. Regression check on existing GRIP tests

```
266 passed, 1 failed (pre-existing: test_real_v2_library_smoke_if_present
checks lookup_bond("C","sp3","C","sp3") which exposes the same
single-neighbour-key architectural limitation -- not caused by this work)
```

## File changes (paths absolute)

* `delfin/fffree/grip_mogul_lookup.py` — TM-aware fallback (+~200 lines),
  v4 pair-table loaders, env-gated default-off integration.
* `tests/test_grip_lib_v4_pair.py` — new, 8 tests.
* `scripts/grip_build_mogul_lib_v4_ccdc.py` — new, full CCDC pair-table
  extraction (in agent_workspace/quality_framework, not in DELFIN).
* `reports/grip_lib_v4.npz` — new artefact (in agent_workspace, full
  build running in background).
* `iters/CCDC_FULL_EXTRACTION_2026_06_05.md` — this report.

## Activation in production

To opt-in to v4 priors in the running pipeline:

```bash
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v4.npz
export DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK=1
```

Without either env var, the system is byte-identical to today. Both
must be set for the v4 metal-organic bond priors to take effect.

## Constraints met

* Worktree isolation (`comprehensive-ccdc-extraction` branch in
  `.claude/worktrees/`).
* No push to GUPPY/main; new branch only.
* Default OFF byte-identical preserved (env flag).
* PYTHONHASHSEED=0 set at module/script import.
* CCDC API only at build-time; runtime uses static npz.
* Determinism preserved end-to-end.
