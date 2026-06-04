# Aromatic Bond Trace — Forensik (Iteration 1)

## Test case
`e9f69af-iso-topo-heal-smoke500/080-LUHMOT.xyz` (Ir-complex, 81 atoms, Isomer 5).

Measured C-C distribution (heavy graph, d<1.65 Å threshold):
- n = 30 bonds
- mean = 1.394 Å (looks aromatic-ish)
- **std = 0.098 Å** (catastrophic variance)
- min = 1.130 Å (C13-C14, triple-bond length, broken aromatic)
- 10 bonds below 1.36 Å

10 worst:
```
C13-C14 = 1.130   (in 6-ring with 16,17,18,19)
C 1-C 7 = 1.146
C13-C19 = 1.258
C 5-C 7 = 1.264
C 3-C 4 = 1.310
C33-C34 = 1.312
C34-C35 = 1.333
C 4-C 5 = 1.346
C28-C29 = 1.350
C40-C41 = 1.352
```

Neighbourhood: atom 13 (C) has phantom 1.13 Å contact to C14 plus ring-distance 1.95-2.20 Å contacts to C16/17/18 (atoms that should be ring-bonded at 1.40 Å). The aromatic 6-ring 13-14-16-17-18-19 has internally collapsed during build/embed.

## Where bond-length targets are SET

### (1) `_bond_decollapse._ideal_bond(a, b)` — naïve sum-of-covalent
File `delfin/_bond_decollapse.py:68-69`:
```python
def _ideal_bond(a, b):
    return _COV.get(a, 0.9) + _COV.get(b, 0.9)
```
`_COV['C']=0.76` → C-C = 1.52 Å (single-bond target). NEVER 1.40 (aromatic). NO bond-type input.

Consumers of `_ideal_bond`:
- `_count_collapsed` (factor 0.82)
- `_md_snapshot`, `_md_broken` (factor 1.30 = M-D detection)
- `_geometric_bonds` (factor 1.30)
- `_count_vdw_clashes`
- `_count_bond_distort` (band [-0.25, +0.085])
- `correct_xyz` (the spring relaxer pulls bonded pairs toward `_ideal_bond`, so aromatic C-C is pulled toward 1.52, OVER-stretching aromatic rings).

### (2) `grip_fragment_detect.py:367` — IS aromatic-aware
```python
aromatic = _atom_is_aromatic(atom_a) and _atom_is_aromatic(atom_b)
hit = lib.lookup_bond(z1, hyb1, z2, hyb2, rmin, aromatic, min_n=min_n)
```
This is correct. But:
- only runs when `DELFIN_FFFREE_GRIP=1`
- a bond is SKIPPED if its endpoints are in `frozen_or_donor` (metal/donor protection)
- a bond is SKIPPED if either endpoint is in `hapto_set` (η-protection)
- a bond is SKIPPED if the library has no hit for that (z1,hyb1,z2,hyb2,rmin,aromatic) tuple

So aromatic rings INSIDE the chelate / on donor-shell-1 paths get NO GRIP cure.

### (3) `aromatic_ring_scale.py` (Phase G13) — uniform RING scale
Only changes scale if `|mean_bond - 1.40| > 0.05` AND clamped 0.85..1.05 (max 5% expansion).
- LUHMOT 6-ring 13-14-16-17-18-19 mean = (1.130 + 1.434 + ring-edges-1.95-2.20) → way mixed. Probably outside the 5/6-ring "aromatic" detection because mean is in [1.30, 1.65]. Even if it triggered, uniform scaling does NOT reduce intra-ring variance.

### (4) `topology_healing.py` — NO compressed-bond detection
Only `detect_phantom_bonds` (spatially bonded but topologically NOT) and `detect_missing_bonds` (topologically bonded but stretched > MISSING_FACTOR × ideal, default 1.5). Compressed bonds (d < ideal) are NEVER healed.

## Root cause

The LUHMOT 1.130 Å C-C bond is an **embed/orient artifact**: during DG bounds embed of the chelate, ring carbons project onto a near-degenerate configuration, then chelate orient + Kabsch rotation + radial donor placement leaves the ring squashed. The current cure stack does not catch this:

1. `_bond_decollapse.correct_xyz` pulls C-C toward 1.52 (single-bond ideal) — would over-stretch aromatic rings AND its accept-if-better gate rejects partial fixes inside an aromatic ring (see `feedback_refine_aromatic_pairwise_plateau`).
2. `aromatic_ring_scale` is **mean-based** and **uniform** — does nothing for per-bond variance, and skips rings whose mean is in [1.30, 1.65].
3. `topology_healing` ignores compressed bonds.
4. `grip_polish` would help but is gated by env, and ring atoms inside chelate may be filtered out.

## What to fix (Iteration 2-6 plan)

### Iter 2: Empirical CCDC aromatic ideals → `aromatic_bond_targets.py`
A static element-pair table of aromatic ideal bond lengths (C-C 1.395, C-N 1.337, C-O 1.362, C-S 1.713, N-N 1.345, etc.), keyed on RDKit aromaticity, replacing naïve covalent-sum for any aromatic edge. Built from CCDC fragment statistics or canonical reference values (Allen 1987 / IUCr).

### Iter 3: Build-time aromatic ring REWRITE → `aromatic_bond_enforcement.py`
- Detect every aromatic ring (RDKit `GetIsAromatic` on **both** atoms of a heavy-heavy edge)
- For 5/6-rings: replace ring coordinates with an ideal regular polygon at the in-ring mean-COD ideal C-C 1.395 Å, oriented to best-fit the current ring centroid+normal (so M-D and ligand placement preserved)
- Substituent and H drag rigidly with parent ring atom
- Skip rings any atom of which is metal-coordinated (preserves hapto-π)
- env-gate `DELFIN_FFFREE_AROMATIC_BONDS=1`, default-OFF byte-identical

### Iter 4: GRIP loss aromatic-aware fallback in `_bond_decollapse._ideal_bond`
- Add an `aromatic=False` parameter; when True, dispatch to `aromatic_bond_targets`
- Update `_bond_decollapse.correct_xyz` to pass `aromatic=…` derived from a passed-in RDKit mol (when available, else legacy single-bond ideal). Env-gate `DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1`, default-OFF byte-identical.

### Iter 5: Smoke 200 validation per change
Measure intra-ring C-C std on aromatic rings (mean over all built complexes).

### Iter 6: Symmetry-equivalence GLOBAL extension
Graph-automorphism (networkx) → bond/angle equivalence classes → variance penalty (α·std(class)) in GRIP loss.

### Iter 7: CCDC GeometryAnalyser verdict on smoke output

### Iter 8+: iterate until aromatic-C-C std < 0.02 Å mean across smoke output.

## Severity of LUHMOT (paper SI snippet)

C13-C14 = 1.130 Å is a clear bug: shorter than a C≡C triple bond (1.20 Å), inside what should be an aromatic 6-ring. ALL build-time / refinement passes ran (topology_healing flag was ON), yet none caught it. This trace explains why.
