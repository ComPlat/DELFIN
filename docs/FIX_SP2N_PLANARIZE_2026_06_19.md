# SP2N-PLANARIZE — universal sp2-N / nitro planarisation fixer (2026-06-19)

Branch: `feat-nitro-sp2n-2026-06-19` (off `main` HEAD `68356218`).
Flag: `DELFIN_FFFREE_SP2N_PLANARIZE` (default `0` = OFF, byte-identical).

## Problem (observed)
The AVIDAM (CONSOL archive) nitro group was severely distorted in the built
complex:
- N-O = 1.43 / 1.12 Å (target ~1.22 each, symmetric)
- O-N-O = 107° (target ~125°)
- N **1.00 Å out of the C-O-O plane** (pyramidalised; nitro N is planar sp2)

This is the signature of an sp2 nitrogen treated as sp3 / pushed out of plane
during complex assembly.

## Part A — confirm test: is MMFF (conformer_enum) the cause?
Tested on the experimental branch `feat-conf-nommff-proof-2026-06-18`
(toggle `DELFIN_FFFREE_CONF_NOMMFF`) at the **free-ligand** level
(`conformer_enum.enumerate_conformers`), MMFF-ON vs MMFF-OFF:

| Ligand SMILES                            | MMFF ON N-O / O-N-O / oop | MMFF OFF N-O / O-N-O / oop |
|------------------------------------------|---------------------------|----------------------------|
| `c1ccccc1N(=O)=O` (nitrobenzene)         | 1.24/1.24, 124.5°, 0.0    | 1.24/1.24, 124.5°, 0.0     |
| `O=C(NN=Cc1ccccc1N(=O)=O)c1ccncc1`       | 1.237/1.237, 125.8°, 0.002| 1.237/1.238, 125.6°, 0.002 |
| `OC(=O)c1cccc(c1)N(=O)=O`                | 1.239/1.239, 124.6°, 0.0  | 1.239/1.239, 124.6°, 0.0   |

**Verdict: MMFF is NOT the cause.** The free-ligand nitro is perfect in BOTH
modes. The defect therefore arises in the **complex-assembly** stage
(metal-sphere construction + global heavy-atom relax), not in the
conformer-enumeration MMFF. The fix must run as a post-assembly geometry
repair, so it is wired as a per-frame post-B5 fixer (after F25).

## Why `_fix_sp3_n_pyramidality` (F25) does not catch it
F25 is the **inverse** fixer:
- it targets **sp3** N that is too **flat** (3-angle sum > 348°) and makes it
  **more pyramidal** (target 328°);
- it only displaces **H** substituents — a nitro N has no H, so it is skipped
  with `no_H_substituent`;
- it filters to RDKit `SP3` hybridisation; nitro N is `SP2`, so it is excluded
  up front;
- the nitro defect needs the opposite correction (flatten an over-pyramidal
  sp2 N), which F25 cannot produce.

A new, dedicated fixer was therefore required.

## Fix
`delfin/_fix_sp2n_planarize.py` :: `planarize_sp2_nitrogen(xyz, mol, ...)`.

Detection (universal, graph-driven, never SMILES-specific):
- **nitro**: N bonded to exactly 2 O + 1 C, total degree 3, not in a ring,
  not metal-bonded, and neither O bonded to a metal (coordinated nitrite is
  left to the coordination path).
- (optional, `DELFIN_FIX_SP2N_AMIDE_IMINE=1`) acyclic RDKit-sp2 amide/imine N.

Construction for each triggered nitro (oop > 0.20 Å, |ΔN-O| > 0.12 Å, or
|O-N-O − 125°| > 12°):
1. Plane through {C, O1, O2}; normal = (O1-C) × (O2-C).
2. Keep C and the C-N bond length; project the C→N direction into the plane;
   place N along it (N becomes coplanar).
3. Place O1, O2 symmetrically in that plane at `no_target` (1.22 Å) from N,
   separated by `ono_target` (125°), with the O1 side chosen to minimise
   motion.
Only {N, O1, O2} move; the carbon skeleton stays rigid.

Safety / doctrine:
- M-D safe: metal-bonded N and metal-bonded O are skipped entirely.
- Never-worse: per-group rollback if planarity does not improve or a NEW
  non-bonded heavy clash appears.
- Deterministic: closed-form, no RNG, no FF call.
- Never non-finite: guarded; rollback to input on any NaN/Inf.
- Bit-exact OFF: the dispatch helper returns the input untouched when the flag
  is 0.

Env-flags:
```
DELFIN_FFFREE_SP2N_PLANARIZE=0    # master, default OFF (byte-id)
DELFIN_FIX_SP2N_OOP_DEG_A=0.20    # out-of-plane trigger (Å)
DELFIN_FIX_SP2N_NO_A=1.22         # target N-O length (Å)
DELFIN_FIX_SP2N_ONO_DEG=125.0     # target O-N-O angle (°)
DELFIN_FIX_SP2N_AMIDE_IMINE=0     # extend to amide/imine N
DELFIN_FFFREE_SP2N_PLANARIZE_CLASSES=  # optional class allow-list
```

Wiring: `delfin/smiles_converter.py`
`_apply_fixer_sp2n_planarize_if_enabled(mol, results, dual_parse_done)`,
called after `_apply_fixer_f25_if_enabled` in the post-B5 fixer chain.

## Validation
- Synthetic broken nitromethane (N 0.81 Å oop, 1.43/1.12 Å, 82.9°) →
  after fix: 1.220/1.220 Å, 125.0°, oop 0.000, C-N length preserved.
- byte-id OFF: see commit message / report.
- pytest fffree: green.

## Status
default-OFF, byte-id when OFF. NOT flipped. NOT pushed. `main` untouched.
