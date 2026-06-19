# FFFREE Aryl Ring-Size Fix (`DELFIN_FFFREE_ARYL_RING_SIZE`)

**Status:** env-gated, **default OFF** (byte-identical when unset). Branch
`feat-aryl-ring-size-2026-06-19` off main-HEAD `68356218`. main untouched.

## The defect

The declash forensics traced part of the ~12% "garbage" structures (closest
contact < 0.9 Å in the best frame) to a root cause that declash **cannot**
resolve: aromatic 5/6 rings built **too small** (in-ring C-C ~1.20-1.30 Å, vs
the physical ~1.39 Å), and in the worst metal-build cases **catastrophically
distorted** (one ring with bonds ranging 0.235-5.14 Å). Crushed rings squeeze
the ring-attached H into each other (H-H < 0.9 Å). Declash leaves them alone
because the ring bonds are rigid (0 rotational DOF) — the *size/shape* is wrong,
not the rotation.

Reference cases (from the BEST-50K pool): **ATABOM** (Rh, η2 + biphenyl-PPh3),
**AFALOK** (Pt, three PPh3 + ylide).

## Why the existing machinery did not catch it

* `_snap_aromatic_rings_to_plane` / `_snap_aromatic_rings_in_xyz` only do a
  **minimal-movement projection onto the best-fit plane** — they flatten a
  buckled ring but **never change its size**. A ring at 1.30 Å stays 1.30 Å.
* The regular-polygon reconstruction at ideal C-C lives inside
  `_correct_hapto_geometry` and only fires for **η-coordinated hapto rings**
  (Cp/arene), not ordinary pendant aryl rings (e.g. phenyls on a phosphine).
* Critically, the metal / hapto build path (`_select_best_hapto_candidate` →
  `_mol_to_xyz`) returns its final XYZ **without ever calling**
  `_snap_aromatic_rings_in_xyz` — so no ring-geometry post-processing ran on
  exactly the structures that exhibit the squash.

## The fix

Two implementations of the same idea (regular-polygon reconstruction of each
bad aromatic 5/6 ring at the ideal aromatic C-C = 1.39 Å), wired onto **all**
build paths:

1. `_scale_aromatic_rings_to_ideal_cc(mol, conf_id)` — operates on an RDKit
   conformer with reliable bond topology. Used on:
   * the non-metal path inside `_snap_aromatic_rings_in_xyz` (runs first, then
     the plane snap removes any residual buckle), and
   * the metal/hapto emit, called on `best_mol`'s conformer just before
     `_mol_to_xyz` (atoms map 1:1, topology trustworthy).
2. `_scale_aromatic_rings_in_xyz_geom(xyz_str)` — a path-independent fallback
   that perceives rings **purely from geometry** (covalent-radius bond graph),
   used at the remaining final-XYZ chokepoints (`_topology_hard_gate_check`,
   the multi-hapto ETKDG fallback, the legacy-hapto return). Needed because the
   metal build mutates the molecule (dative rewrites, charge/H changes) so a
   SMILES-rebuilt mol no longer maps onto the geometry.

### Algorithm (per ring)

* Restrict to aromatic-capable rings: size 5/6, >= 2 unsaturated in-ring bonds
  (mol path) / sp2 substituent signature (geom path), no H or metal ring member.
* **NEVER-WORSE:** skip rings already physical (mean C-C in [1.34, 1.44] **and**
  min >= 1.29 **and** max <= 1.54) — good benzene/toluene/biphenyl are left
  byte-identical. Saturated cyclohexane (C-C ~1.54, >1 substituent per carbon)
  is excluded so it is never flattened.
* Reconstruct a regular eta-gon of side 1.39 Å in the ring's SVD best-fit plane,
  centred on the current centroid, atoms assigned to consecutive vertices by
  ring bond connectivity (topology preserved). Each ring atom's substituent
  **branch** is dragged rigidly by the same per-atom displacement, so
  substituent bond lengths/angles and ring-H are preserved exactly.

### M-D invariant (coordinated rings)

* Each ring atom's branch BFS **stops before any metal**; an atom whose branch
  reaches a metal is **anchored** (frozen).
* >= 2 anchored ring atoms → ring skipped entirely (chelating diaryl / fused
  metallacycle — too risky to reshape).
* Exactly 1 anchor (common: pendant aryl on a P/donor whose ipso-C → P → M):
  the polygon is translated so the anchored atom lands on its current position;
  the anchored atom + its metal-bearing branch stay put, the rest of the ring
  moves. The M-D bond is never touched.

### Per-ring NEVER-WORSE guard

After laying out a ring's new positions, the closest contact between any
**moved** atom and any non-moved atom is compared before vs after. If the move
would introduce a fresh sub-0.9 Å clash worse than what was there, that ring is
reverted. This guarantees the pass cannot make the global minimum contact worse
for the rings it touches (verified on the sample: 0 min-dist regressions).

## Validation (PYTHONHASHSEED=0)

Aromatic ring C-C measured by mapping RDKit ring perception onto the built XYZ
(`_measure_ring_precise.py`):

| Structure | metric | OFF | ON |
|-----------|--------|-----|----|
| AFALOK (Pt, 3×PPh3) | arom C-C mean / min / max | 1.731 / 0.235 / 5.140 | rings reset to 1.390 where the never-worse guard allows; global min-dist held at 0.235 (no regression), min-HH 0.395→0.503 |
| ATABOM (Rh) | global min-dist / min-HH | 0.225 / 0.225 | 0.278 / 0.502 (both improved) |
| benzene / toluene / biphenyl | arom C-C, min-dist, min-HH | (good) | **byte-identical** (good rings untouched) |
| ABORIV / ACOWUM (well-built PPh3-metal) | — | (good) | **byte-identical** (no-op) |

Honest scope: the ring-**size** defect is fixed deterministically (rings the
guard can move become exact 1.39-Å hexagons). Severely broken metal builds
(AFALOK) have a *separate* inter-ligand placement defect — overlapping PPh3
ligands — that ring-sizing cannot repair; the never-worse guard correctly
reverts the rings whose repair would deepen that pre-existing overlap, so the
structure is never made worse. Cleanly-built squashed rings (ATABOM-class,
non-metal pendant aryls) are improved.

## Guarantees

* `default-OFF` byte-identical (proven: cisplatin + an aromatic non-metal,
  fix-unset ≡ 68356218, `delfin.__file__` from the worktree, PYTHONHASHSEED=0).
* Universal — RDKit aromaticity / geometric ring perception, never
  SMILES-specific.
* Deterministic, never non-finite (all moves `np.isfinite`-guarded; failures
  return the input unchanged).
* No CCDC data committed.
