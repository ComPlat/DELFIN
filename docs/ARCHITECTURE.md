# DELFIN Architecture — FF-free TMC Pipeline

## Two-Path Construction

DELFIN combines **chemistry-aware constructive build** for known classes with
**universal embed fallback** for edge cases. Both paths support a
**user-selectable polish stage** so adopters can choose the quality/speed
tradeoff that suits their use case.

User-vision (2026-06-06): the constructive path targets ≥ ETKDG/RDKit
coverage; both paths together = 100 % coverage with paper-grade quality.

```
SMILES
  │
  ├─→ Universal Classification (ALWAYS runs)
  │   - Topological isomer enumeration (Pólya-complete)
  │   - Linkage/spin/charge isomer enumeration
  │   - Labels: ["fac", "mer", "Δ", "Λ", ...]
  │
  ├─→ Per isomer-label: 3D-Build dispatch
  │
  │   Path 1 — Constructive (chemistry-aware)
  │   ┌─────────────────────────────────────────┐
  │   │ Polyhedron vertex assembly              │
  │   │   (OC-6, TBP-5, SAP-9, sandwich, etc.)  │
  │   │ Cremer-Pople ring puckering             │
  │   │ Aromatic-plane + symmetry homogenisation│
  │   │                                         │
  │   │ → GRIP-polish (FF-free)                 │
  │   │ → GRACE conformer ensemble              │
  │   └─────────────────────────────────────────┘
  │       Coverage target: ≥ ETKDG/RDKit (F1 expansion)
  │
  │   OR (when class unknown)
  │
  │   Path 2 — Universal Embed
  │   ┌─────────────────────────────────────────┐
  │   │ RDKit ETKDGv3 distance-geometry embed   │
  │   │   (deterministic seed=42, no random)    │
  │   │                                         │
  │   │ → Selectable post-process (F3):         │
  │   │   uff / grip / xtb / none / both / all  │
  │   └─────────────────────────────────────────┘
  │       Coverage: catches every edge case (F2/F3)
  │
  └─→ Output: Pólya-complete isomers, CP-complete conformers,
              FF-free + user-selectable polish, deterministic
```

## FF-free Polish — Mathematics

UFF target: `E_UFF = Σ k(r - r_uff_ideal)²` with parameters from 1992 gas-phase fits.

DELFIN GRIP target: `E_GRIP = Σ z²` where `z = (r - μ_CCDC) / σ_CCDC` from 1.7M crystal observations.

Same Hooke-Spring functional form, but:
- **GRIP-targets are crystal data**, not gas-phase parameters
- **Per-class disaggregation**: 3d/4d/5d/f-block separate
- **TM-categories**: Carbene/η²/η⁵/η⁶/μ-Bridge first-class
- **Cleaning applied**: counter-ions/solvent/disorder stripped pre-extraction
- **Multimodal handling**: Gaussian Mixture for bimodal torsions

## Dual-Source Library (v5 CCDC + v6 COD)

```
grip_lib_v5.npz (CCDC, 188 MB)
  - 1,436,119 entries processed (83.58% extraction success)
  - 2,449,647 master keys
  - 7 TM-categories: carbene, hapto-η², η⁵, η⁶, μ-bridge, agostic, ox-addition
  - Block-disagg: 3d/4d/5d/f separately aggregated

grip_lib_v6_cod.npz (COD, 18 MB)
  - 50,813 / 51,192 entries (99.26% extraction)
  - Adopter-friendly: no CCDC licence needed at runtime
  - 83.75% cross-validation agreement with v5

MergedLibrary
  - N-weighted closed-form mean+sigma merge
  - Provenance flags: both-agree / both-marginal / both-disagree / ccdc-only / cod-only
```

## Fallback Mode (F2 + F3) — 6 user-selectable polish stages

`DELFIN_FFFREE_FALLBACK_MODE` env-flag:

| Mode | Pipeline | Story |
|---|---|---|
| `uff` (DEFAULT, unset) | RDKit ETKDG + UFF minimise | backward compatibility |
| `grip` | RDKit ETKDG + FF-free GRIP polish | FF-free production paper claim |
| `xtb` (F3) | RDKit ETKDG + GFN2-xTB relaxation | semi-empirical reference |
| `none` | Return `[]` on fffree-native fail | paper-grade "no fallback" measurement |
| `both` | grip + uff (concatenated) | A/B comparison |
| `all` (F3) | grip + uff + xtb per conformer | full A/B/C comparison for the paper |

Behaviour notes:
- Mode-setting (any non-`uff` value) also implicitly enables fffree dispatch,
  so adopters only need to flip one flag.
- Default unset = `"uff"` is **byte-identical** to the pre-F2 legacy pipeline.
- `xtb` mode probes
  `/home/qmchem_max/micromamba/envs/delfin/bin/xtb`, then
  `/home/qmchem_max/micromamba/bin/xtb`, then `$PATH`.
  When no xtb binary is found, the polish is **silently skipped** and the
  output carries an `-raw` suffix so the caller can detect partial polish.
- `all` mode emits 3× outputs per conformer (one each grip/uff/xtb branch);
  the legacy UFF concatenation is **not** added on top of `all` so UFF is
  not duplicated.
- xtb runs are pinned single-threaded (`OMP_NUM_THREADS=1`) and
  `--norestart` so the byte-identity contract holds across parallel
  workers.
- Per-call xtb wall-clock cap: `DELFIN_FFFREE_XTB_TIMEOUT` seconds
  (default 120, clamped to `[5, 1800]`).

Paper claim (F3): *DELFIN supports user-selectable post-processing —
`uff` / `grip` / `xtb` / `none` — over the same construction path.
The polish stage is a single env-flag for end-users.*

## Determinism

- `PYTHONHASHSEED=0` mandatory
- ETKDG seed=42 fixed, `useRandomCoords=False`
- L-BFGS Mersenne-Twister seeded per restart
- Sorted iteration (no dict-order dependence)
- Float64 throughout
- 2-run byte-identical verified per merge

## Coverage Expansion Plan (F1)

Per failure-class, add: polyhedron + Pólya-group + dispatcher + tests.

| Class | Polyhedra needed | Status |
|---|---|---|
| CN3-9 standard | All implemented | ✓ |
| f-block CN8-12 | Phase C done | ✓ |
| Sandwich, piano-stool, half-sandwich | A7 done | ✓ |
| CN10 BICAP/CSAP/SAP non-f-block | A2 done | ✓ |
| Multi-metal μ-bridge | F1 Bucket 1 | in progress |
| Mixed hapto+σ | F1 Bucket 2 | in progress |
| Spiro/bridged chelates | F1 Bucket 5 | queued |
| Carbene Fischer/Schrock | F1 Bucket 6 | queued |
| Exotic donors Se/Te/As/Sb | F1 Bucket 7 | queued |
