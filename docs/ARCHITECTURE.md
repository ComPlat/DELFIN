# DELFIN Architecture — FF-free TMC Pipeline

## Two-Path Construction

DELFIN combines **chemistry-aware constructive build** for known classes with
**universal embed fallback** for edge cases. Both paths apply the same
**FF-free Mahalanobis polish** on CCDC/COD crystal statistics.

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
  │   Constructive Path (chemistry-aware)
  │   ┌─────────────────────────────────────────┐
  │   │ Polyhedron vertex assembly              │
  │   │   (OC-6, TBP-5, SAP-9, sandwich, etc.)  │
  │   │ Cremer-Pople ring puckering             │
  │   │ Aromatic-plane + symmetry homogenisation│
  │   │                                         │
  │   │ → GRIP-polish (FF-free)                 │
  │   │ → GRACE conformer ensemble              │
  │   └─────────────────────────────────────────┘
  │       Coverage: 20% → 70%+ (F1 expansion)
  │
  │   OR (when class unknown)
  │
  │   Universal Embed Path
  │   ┌─────────────────────────────────────────┐
  │   │ RDKit ETKDGv3 distance-geometry embed   │
  │   │   (deterministic seed=42, no random)    │
  │   │                                         │
  │   │ → GRIP-polish (FF-free, same as above)  │
  │   └─────────────────────────────────────────┘
  │       Coverage: catches edge cases (F2)
  │
  └─→ Output: Pólya-complete isomers, CP-complete conformers,
              all FF-free, all deterministic
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

## Fallback Mode (F2)

`DELFIN_FFFREE_FALLBACK_MODE` env-flag:
- `grip` (DEFAULT, FF-free): ETKDG embed + GRIP-polish
- `uff`: legacy UFF (backward compatible)
- `none`: return empty (NO_FALLBACK measurement)
- `both`: emit both grip-polished AND uff outputs (comparison)

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
