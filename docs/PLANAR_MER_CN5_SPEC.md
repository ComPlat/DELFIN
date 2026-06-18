# Rigid Planar Tridentate — CN5 (TBP-5 / SPY-5) Meridional Geometry Fix

Flag: `DELFIN_FFFREE_PLANAR_MER_CN5` (default `"0"` -> byte-identical, branch never entered)

## Problem

The planar-tridentate -> meridional fix (`DELFIN_FFFREE_PLANAR_MER`, commit `e2c00c7c`)
solved the OC-6 class (ADERAE: fac -> mer, outer-outer N-Cd-N 141.6deg). For CN5
(TBP-5 / SPY-5, ANUCOE class) the rigid planar tridentate gets the correct
*meridional vertex ASSIGNMENT* (via `_meridional_triples` + `enumerate_chelate_configs`),
but the *placement geometry* stays FOLDED:

- `assemble_from_config` -> `_embed_metallacycle` embeds the chelating conformer
  WITHOUT pinning the donor-donor distances to the assigned vertices (the bite
  constraint is gated behind a SEPARATE flag `DELFIN_FFFREE_CHELATE_BITE`).
- `_orient_chelate_to_vertices` then rigid-Kabsch-rotates that conformer onto the
  target vertex directions and only rescales each donor's RADIUS to the M-D ideal —
  it KEEPS each donor's embedded DIRECTION, i.e. the embed's natural folded bite.
- A single rigid rotation cannot OPEN a ~110deg bite onto the meridional TBP axis.

ANUCOE baseline (V, CN5, TBP-5, O,N,O Schiff-base tridentate + 2 monodentate O):
- outer-outer O-V-O angle = **110.6deg** (should be meridional ~150-165deg)
- CShM vs TBP-5 = **12.696** (poor)
- M+3-donor plane RMS high; the two monodentate O crowd the folded tridentate
  (donor-donor angles 61.8deg / 65.7deg).

ANUCOE crystal truth (CCDC index, internal validation only):
- 5-coordinate, donors 4xO + 1xN, intermediate TBP/SPY.
- meridional outer-outer pair O-V-N = **150.2deg** (the tridentate spans ~150deg).

## Fix

When `DELFIN_FFFREE_PLANAR_MER_CN5=1` AND the ligand is tagged `rigid_planar`
(by the existing `_is_rigid_planar_tridentate`, only reachable under
`DELFIN_FFFREE_PLANAR_MER=1`) AND denticity 3 AND the geometry is CN5
(TBP-5 / SPY-5): the metallacycle embed is forced through the distance-geometry
BITE CONSTRAINT (the same machinery as `DELFIN_FFFREE_CHELATE_BITE`, but activated
independently for this scoped class). The DG bounds pin:
  - each M-donor distance to the ideal `polyhedra.md_distance(metal, donor)`
  - each donor-donor distance to the ideal meridional-vertex separation
    (`polyhedra.ref_vectors` for TBP-5 / SPY-5 scaled by `md_distance`)
so the embedded conformer ALREADY carries the meridional bite. The rigid Kabsch
orient then only chooses the orientation; the bite is correct by construction.
The rigid conjugated backbone is dragged along by the DG embed (it is one
quasi-rigid unit), and the M-D target distances are exactly the ideal so the
M-D invariant is preserved (each donor at `md_distance`).

SPY-5: `_meridional_triples` is extended so the SPY-5 trans-basal arc (apical
central + trans-basal outer pair, outer-outer ~163.9deg, coplanar through the
metal) is admitted (TBP-5 already returns axial-axial+equatorial triples at
180deg/90deg). Universal, derived geometrically from the ideal vertex set; no
SMILES / refcode keying.

## Never-worse self-gate (`_build_config_never_worse`)

A blanket forced bite is NOT uniformly better: it opens the bite toward the ideal
for most rigid-planar CN5 ligands, but for some it would DISTORT the coordination
shape WORSE (measured: QOFLOP CShM 0.48 -> 10.6, AGEZUK 2.94 -> 12.7 under a blanket
force).  So `_fffree_chelate_isomers` routes each rigid-planar-CN5 config through
`_build_config_never_worse`, which builds the frame BOTH ways — `planar_bite=True`
(meridional bite-constrained embed) and `planar_bite=False` (the historic folded
seating) — and KEEPS the one with the lower FULL-SHELL CShM (`polyhedra.cshm` over
the coordination donors vs the target geometry).  The meridional build wins where it
helps; the folded build survives where it would not (strict never-worse).

`assemble_from_config(planar_bite=...)`: a per-call override of the env gate (None =
use `DELFIN_FFFREE_PLANAR_MER_CN5`).  The never-worse caller passes True/False to get
both builds; the env gate alone is used elsewhere.  The M-D distances are pinned to
the ideal (`md_distance`) in the bite-constrained embed, preserving the M-D invariant
within the DG tolerance (+-0.05A bounds).

Validated never-worse (14 reference rigid-planar CN5 builders, CShM OFF -> ON):
ANUCOE 12.7->2.1, AXOKAZ 5.86->0.67, FIPJIX 16.1->3.8, BOJZIJ 13.3->9.1, ACITIP
15.6->0.5, ADEQUX 4.7->2.3 (improved); QOFLOP 0.48->0.48, AGEZUK 2.94->2.94, NUBDAR
6.53->6.53, CAZTEC 11.6->11.6, ABUBEI 2.01->2.01, AFUTEB 3.36->3.36 (folded kept,
unchanged); AMUKUP None->8.4, OJIGUK None->7.3 (newly build).  EVERY case <= OFF.

## Validation (ANUCOE, flags
`DELFIN_FFFREE_BUILDER=1 DELFIN_FFFREE_CN_EXTEND=1 DELFIN_FFFREE_PLANAR_MER=1
DELFIN_FFFREE_PLANAR_MER_CN5=1`, `PYTHONHASHSEED=0`):

- outer-outer O-V-O: 110.6deg -> **153.8deg** (meridional)
- CShM vs TBP-5: 12.696 -> **2.098**
- M+3-donor plane RMS -> 0.199A
- central N at ~90deg to both outer donors (correct meridional central seating)
- no donor-donor clashes (min 2.757A)

## Byte-identity

Default OFF -> `_planar_mer_cn5_enabled()` returns False -> the forced bite is
never applied and the SPY-5 trans-basal triples are never emitted (the
`_meridional_triples` extension is itself gated). cisplatin / Pt(en)Cl2 / ADERAE
(OC-6) builds are byte-identical with the flag unset vs the parent commit
`e2c00c7c`. The CN5 fix never touches the OC-6 path (geometry scoped to CN5).
