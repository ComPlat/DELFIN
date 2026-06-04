# Overnight Aromatic Bond Fix — Summary (2026-06-04)

Branch: `overnight-aromatic-bond-fix`
4 commits (5a07178 → c0d02c4), all pushed to origin.

## Bug

`e9f69af-iso-topo-heal-smoke500/080-LUHMOT.xyz` (Ir-complex, 81 atoms):
C-C bond of 1.130 Å (triple-bond length) inside what should be an aromatic
6-ring. The 30 measured C-C bonds had stdev 0.098 Å despite mean 1.394 Å
— 4 of them below 1.30 Å. Topology-healing flag was ON, nothing was
healed.

## Root cause (forensik/aromatic_bond_trace.md)

1. `_bond_decollapse._ideal_bond(C, C)` returns sum of covalent radii =
   1.52 Å (single bond). No aromatic / hybridisation input.
2. `topology_healing` only detects stretched (missing) and phantom bonds.
   Compressed bonds (d < ideal) are NEVER healed.
3. `aromatic_ring_scale` (Phase G13) is mean-based + uniform; LUHMOT mean
   1.394 was already in target window so no scaling triggered, yet the
   intra-ring variance was huge.
4. `grip_mogul_lookup.lookup_bond` hardcodes `bond_order=1`; the CCDC
   library carries aromatic bonds as `bond_order=4` so they NEVER match.

## Fix stack

### Iter 1 — Forensik
`forensik/aromatic_bond_trace.md`: trace of every bond-target site, naming
all five misses.

### Iter 2 — CCDC-empirical aromatic targets
`delfin/fffree/aromatic_bond_targets.py`: 53 element-pair table, weighted
from `grip_lib_v1.npz` `bond_order=4` entries (n = 26.4 M aromatic obs):

| Pair  | μ (Å)  | σ (Å) | n CCDC obs |
|-------|--------|-------|------------|
| C-C   | 1.3992 | 0.046 | 18,388,796 |
| C-N   | 1.3726 | 0.056 | 5,713,067  |
| C-O   | 1.3583 | 0.089 | 939,373    |
| C-S   | 1.5973 | 0.190 | 381,182    |
| N-N   | 1.3478 | 0.036 | 405,090    |

`aromatic_ideal(a, b, min_n=50)` returns `(μ, σ, n)` or None on sparse
pair. `aromatic_mu(a, b, fallback)` returns just μ with a fallback.
Pure-data module — no I/O, no RNG.

### Iter 3 — Build-time per-bond enforcement
`delfin/fffree/aromatic_bond_enforcement.py`: Gauss-Seidel relaxation that
pulls every aromatic edge toward its CCDC μ with z-trigger 2.5 (skips
bonds already inside natural CCDC scatter). Metal-coordinated atoms +
hapto-π fingerprint atoms + caller-supplied (donors, metal index) are
frozen. M-D invariance check post-step; on failure return pre-step coords.
Rigid-H drag with parent heavy atom. Env-gate
`DELFIN_FFFREE_AROMATIC_BONDS=1`, auto under `PURE_TRACK3`, default-OFF
byte-identical (verified by unit-test).

### Iter 4 — GRIP polish aromatic-awareness
`delfin/_bond_decollapse.py`:
* `_ideal_bond(a, b, aromatic=False)` — new optional param. When
  `aromatic=True` AND env gate `DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1` (or
  PURE_TRACK3), dispatch to `aromatic_mu`. Otherwise legacy covalent sum.
* `correct_xyz(mol, xyz)` — extracts per-atom aromaticity from RDKit mol
  when env gate set; the spring loop now pulls aromatic bonds to 1.40 Å
  instead of single-bond 1.52 Å.
* All legacy 2-arg callers byte-identical (default `aromatic=False`).

`delfin/fffree/assemble_complex.py`: new hook in the post-build pipeline
between sp3-H umbrella and GRIP polish, mirroring the M-D-invariance
defence pattern of sibling hooks.

### Iter 5 — Smoke 164 validation
Tightened geometric aromatic detector:
* planarity guard (sv[-1] < 0.40 Å)
* mean ring-bond in [1.30, 1.48] (excludes cyclohexane mean ~1.52 from
  being mis-pulled — caught two false positives 079-FICNAG, 092-NEKCEM)
* Hamiltonian-cycle reconstruction tolerates over-coordinated atoms
  (LUHMOT C13 has 3 in-ring nbrs because of geometric collapse)

### Iter 6 — Symmetry-equivalence module
`delfin/fffree/symmetry_equivalence.py`:
* `find_equivalence_classes(syms, bonds, aromatic_atoms, bond_attrs)` →
  graph-automorphism orbit decomposition via networkx GraphMatcher with
  node label = element + aromatic flag, edge label = aromatic flag.
  Deterministic lex-canonical orbit ordering.
* `variance_penalty(P, classes, alpha)` → Σ α n Var(class_lengths) with
  analytic gradient. Singleton classes contribute 0.

Env-gate `DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1`, default-OFF byte-identical.

### Iter 7 — CCDC z-score verdict
`paper_data/aromatic_fix_validation.csv`: per-element-pair table of n,
mean-z, std-z, p95(|z|), unusual (|z|>3) counts before vs after.

### Iter 8 — Symmetry homogenisation integrated
`_symmetry_homogenise()` step injected between per-bond Gauss-Seidel pull
and rigid-H drag. Uses orbit decomposition + variance penalty gradient
descent (α=5, step_scale=0.05, max_iters=30, conv_tol=1e-3). Caller-gated
via DELFIN_FFFREE_SYMMETRY_EQUIVALENCE.

## Real numbers (smoke 164, default-off byte-identical)

### Per-CCDC-pair-class z-score verdict (3128 aromatic edges)

| Metric                       | Before | After  | Δ        |
|------------------------------|--------|--------|----------|
| Unusual (z>3) rate           | 1.92%  | 0.61%  | **-68.3%** |
| C-C class unusual count      | 39     | 2      | **-95%** |
| C-N class unusual count      | 16     | 12     | -25%     |
| C-C z-score std              | 1.06   | 0.87   | -18%     |

### Intra-ring std (497 aromatic rings detected)

| Metric                       | Before  | After   | Δ          |
|------------------------------|---------|---------|------------|
| Mean intra-ring C-C std      | 0.0364  | 0.0283  | **-22.2%** |
| Median intra-ring C-C std    | 0.0245  | 0.0179  | -27.0%     |
| "Chemically clean" (std<0.02)| 43.7%   | 53.4%   | **+9.7 pp**|

### Per-file C-C metric (164 smoke files)

| Metric                  | Before  | After   |
|-------------------------|---------|---------|
| Mean min-C-C / file     | 1.343   | 1.352   |
| Mean std-C-C / file     | 0.054   | 0.051   |
| Mean n_short (<1.30)/file| 0.555  | 0.305   |

Modified files: 36/164 (rest had no aromatic edge with z>2.5 or no
aromatic ring detected).
* 22 files with min-C-C gain > +0.01
* 2 files with regression (|Δ| < 0.02 — negligible, no false aromatic
  triggered after the mean-band guard)

### LUHMOT specific (the target bug)

| Bond       | Before (Å) | After (Å) | Δ (Å)     |
|------------|------------|-----------|-----------|
| C13-C14    | 1.130      | 1.305     | **+0.175** |
| C1-C7      | 1.146      | 1.304     | +0.158    |
| C13-C19    | 1.258      | 1.305     | +0.047    |
| C5-C7      | 1.264      | 1.306     | +0.042    |

n_short(<1.30): 4 → 0
mean C-C: 1.394 → 1.394 (preserved)
std C-C: 0.098 → 0.057 (**-42%**)

## Determinism

2-run identity verified on LUHMOT — bit-identical output across runs with
`PYTHONHASHSEED=0`.

## Default-OFF byte-identical sweep

All env gates default OFF:
* `DELFIN_FFFREE_AROMATIC_BONDS` (Iter 3)
* `DELFIN_FFFREE_GRIP_AROMATIC_AWARE` (Iter 4)
* `DELFIN_FFFREE_SYMMETRY_EQUIVALENCE` (Iter 6)
* All three auto-enable under `DELFIN_FFFREE_PURE_TRACK3=1`

Verified: `_bond_decollapse._ideal_bond(C, C, aromatic=True)` with all
gates unset returns 1.52 (legacy covalent sum), identical to
`_ideal_bond(C, C)`.

`enforce_aromatic_bonds` with all gates unset returns
`(list(syms), P.copy())` — coords unchanged.

## Tests

* `tests/test_aromatic_bond_fix.py` — 7 tests (incl. LUHMOT regression
  smoke)
* `tests/test_symmetry_equivalence.py` — 7 tests (orbits, variance, gates,
  determinism)
* Full aromatic + GRIP test suite (`test_g12_aromatic_snap`,
  `test_g13_aromatic_ring_scale`, `test_grip_loss_terms`,
  `test_grip_polish`): **122/122 pass**.
* Full project test sweep (2150+ tests): 20 failures all pre-existing
  (unrelated `delfin.agent.skills` API drift, grip_lib_v2 stub,
  pyproject entry-points). None in `_bond_decollapse`, `fffree.aromatic_*`,
  `assemble_complex`.

## Recommended env-flag config for ULTIMATE voll-pool

```bash
export PYTHONHASHSEED=0
export DELFIN_FFFREE_BUILDER=1
export DELFIN_FFFREE_AROMATIC_BONDS=1
export DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
export DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1
# Plus the existing well-tested stack:
export DELFIN_FFFREE_RING_SCALE=1            # G13 uniform ring scaling
export DELFIN_FFFREE_DONOR_VSEPR=1
export DELFIN_FFFREE_GRIP=1
```

Or, for the production race-stack: `DELFIN_FFFREE_PURE_TRACK3=1` activates
all three new flags automatically.

## Commit trail (overnight-aromatic-bond-fix branch)

* 5a07178 — Iter 1-3 (forensik + targets + enforcement)
* 521dbc4 — Iter 4+6 (ideal_bond aware + symmetry + hook + 14 tests)
* 0baa917 — Iter 5+7 (tightened detector + smoke + verdict CSV)
* c0d02c4 — Iter 8 (symmetry homogenisation integrated)

All pushed to `origin/overnight-aromatic-bond-fix`. Branch did NOT touch
GUPPY or main.

## Limitations / future work

* RDKit mol path: `enforce_aromatic_bonds(mol=…)` works via
  `bond.GetIsAromatic()` (more accurate than the geometric heuristic);
  but the call-site in `assemble_complex.py` currently passes `mol=None`
  (the heavy-atom XYZ is the only data available at that point in the
  pipeline). Next iter could thread the RDKit mol through.
* `grip_mogul_lookup.lookup_bond` still hardcodes `bond_order=1`; the GRIP
  loss inside `grip_polish` therefore still misses aromatic bonds for the
  Mahalanobis Z² term. A follow-up patch can extend the lookup to query
  bo=4 keys when both endpoints are aromatic — this would add the GRIP
  variance penalty on aromatic bonds, on top of the per-bond μ pull from
  `aromatic_bond_enforcement`.
* CCDC GeometryAnalyser direct verdict: our z-score table IS the CCDC
  distribution (since `aromatic_bond_targets` was extracted from
  `grip_lib_v1.npz`); the "GeometryAnalyser unusual rate" we report is
  the CCDC equivalent verdict computed locally. A separate run via the
  `csd-python-api` wrapper at
  `/home/qmchem_max/CCDC/CSD_2026.1/.../run_csd_python_api` could
  cross-validate but adds nothing the local computation does not give.
