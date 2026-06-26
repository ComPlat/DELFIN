# SPEC — Universal Aromatic Ring-System Planarisation (Iter-33, 2026-06-19)

Branch: `feat-aromatic-planarity-2026-06-19` (off main-HEAD `68356218`)
Flag:   `DELFIN_FFFREE_AROM_PLANARIZE` (default-OFF, byte-identical when unset)
Module: `delfin/_arom_planarize.py`
Wired:  `delfin/smiles_converter.py` →
        `_apply_arom_planarize_if_enabled(mol, results, dual_parse_done)`,
        called in the post-emit chain right after the Iter-24
        `_apply_aromatic_planarity_if_enabled` and before bond-decollapse.

## The defect
(Hetero)aromatic + fused/polycyclic π-ring systems buckle out of plane.
Forensic anchor — CONSOL-archive `AXAGOY` (Os, 49 atoms, phenanthroline-type
N,N-chelate + 2 pendant phenyls), frame 0:

| ring | elems | coord | OOP (Å) |
|------|-------|-------|---------|
| C5N donor 1 | CCCCCN | yes | 0.1289 |
| C5N donor 2 | CCCCCN | yes | 0.1561 |
| fused central C6 | CCCCCC | yes | 0.1327 |
| pendant phenyl A | CCCCCC | no | 0.0628 |
| pendant phenyl B | CCCCCC | no | 0.0858 |

Flat-aromatic target < 0.05 Å.

## Why heteroaromatics / fused rings stayed bent (root cause)
The existing Iter-24 `_aromatic_ring_flattener`:
1. **Skips every coordinated ring** — `pendant_rings = [r for r in rings if
   not _is_coordinated(r)]`. The worst offenders are exactly the M-bound
   heteroaromatic donor rings → never touched.
2. **Fuses only pendant rings** — a polycyclic chelate that is partly
   coordinated is never flattened as one unit.
3. **Default-OFF + class-gated to {hapto, multi_hapto}** — a σ-class complex
   like AXAGOY never reaches it.

The detector itself is NOT element-biased: it accepts ring N/O/S (it found
both C5N rings). The miss was the coordinated/fused exclusion, not hetero
detection.

## The fix
- Geometric aromatic-ring perception (5/6-membered C/N/O/S ring, mean
  intra-ring heavy bond < 1.46 Å). Universal, no SMILES/refcode/element-list
  specialisation; robust to atom-order drift and to exotic metal valences
  (Os-4) that defeat RDKit sanitisation.
- Fuse ALL aromatic rings sharing an atom into one polycyclic π-system; flatten
  each as a single rigid unit onto its common best-fit plane.
- **M-D preservation**: any ring-system atom within M-D bonding distance
  (< 2.60 Å) of a metal is an *anchor* and is NEVER moved. The best-fit plane
  is constrained to pass through the anchor(s); only NON-anchor atoms are
  projected onto it. With anchors fixed and the metal untouched, every
  M-anchor distance is preserved exactly → M-D invariant holds by construction
  (plus a hard `_md_invariant_violated` rollback as belt-and-suspenders).
  Pure-pendant (no-anchor) systems use centroid-preserving projection
  (M-centroid invariant), as Iter-24 did.
- Ring-attached H + first substituent heavy atoms dragged rigidly onto the
  new plane (projected along the plane normal — keeps in-plane radial, removes
  only the out-of-plane wiggle). Metals are never dragged.

## Guards (NEVER-WORSE + determinism + finiteness)
- Per-system: keep only if that system's worst ring OOP strictly decreased.
- Frame-level: global worst aromatic-ring OOP must not rise.
- Hard M-D / M-H invariant rollback (tol 0.05 Å).
- `np.all(np.isfinite)` guard before emit.
- No RNG, no hash dependence → deterministic (PYTHONHASHSEED-independent).
- Bit-exact to input when flag unset / no qualifying ring / already flat.

## Validation (measured this session)
Run the fix directly on the puckered CONSOL-archive geometries (fix as a pure
post-pass on the emitted XYZ). AXAGOY frame 0 worst-ring OOP:
  before 0.156 Å → after < 0.05 Å (see report).
byte-id OFF proven: cisplatin + a non-aromatic complex, flag unset, output
identical to `68356218`.
