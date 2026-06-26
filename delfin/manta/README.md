<p align="center">
  <img src="../logo/MANTA_readme_demo.gif" alt="MANTA logo animation" width="820">
</p>

<p align="center">
  <a href="../../LICENSE"><img src="https://img.shields.io/badge/License-LGPL--3.0--or--later-blue.svg?style=for-the-badge" alt="License"></a>
  <img src="https://img.shields.io/badge/data-license--clean-success.svg?style=for-the-badge" alt="License-clean">
  <img src="https://img.shields.io/badge/determinism-byte--identical-blue.svg?style=for-the-badge" alt="Deterministic">
  <img src="https://img.shields.io/badge/part%20of-DELFIN-blueviolet.svg?style=for-the-badge" alt="Part of DELFIN">
</p>

# MANTA

**MANTA** is the coordination-structure engine inside [DELFIN](../../README.md). From a transition-metal **SMILES**, it **constructs** the provably-complete set of coordination isomers and expands each into conformers — then energy-ranks the ensemble. It is **deterministic**, uses **no force field at the metal**, and uses **no licensed crystallographic data at runtime**.

> **Construct, don't search.** CREST / GOAT explore the conformational space *stochastically* (metadynamics, genetic search). MANTA *enumerates* the coordination manifold from group theory and ideal geometry, then relaxes and ranks. The two are complementary: MANTA gives a deterministic, complete, license-clean starting manifold; a global optimizer can refine it.

## What it does

| Stage | How |
|-------|-----|
| **Coordination isomers** | Burnside–Pólya enumeration over the donor set on ideal polyhedra — provably complete (verified against 24/24 textbook cases: *fac/mer*, Λ/Δ, *cis/trans*, …) |
| **Geometry** | Donors placed on ideal polyhedron vertices (Pyykkö covalent radii) + VSEPR; metallacycle backbones embedded by distance geometry — **no force field at the metal** |
| **Conformers** | Hindered-rotation / ring-pucker enumeration per isomer, RMSD-deduplicated |
| **Ranking** | Optional GFN2-xTB (or GFN-FF) energy ranking of the emitted ensemble (via the `xtb` CLI) |

Every result is **deterministic** — same SMILES in, byte-identical manifold out, independent of process or hash seed.

## Usage

The user-facing entry point lives in the DELFIN facade (a dedicated `delfin-manta` CLI is on the roadmap):

```python
from delfin.smiles_converter import smiles_to_xyz_isomers

# returns the constructed, deduplicated coordination-isomer manifold as XYZ blocks
isomers = smiles_to_xyz_isomers("[Co](N)(N)(N)(N)(Cl)Cl")
```

Optional GFN2-xTB energy ranking of the ensemble (requires `xtb` on `PATH`):

```bash
DELFIN_FFFREE_GFNFF_RANK=1 DELFIN_CONF_RANK_METHOD=gfn2 python your_script.py
```

All new behaviour is **flag-gated and default-OFF** (byte-identical to the prior pipeline unless explicitly enabled).

## Package layout

`delfin/manta/` holds the engine (68 modules). Key entry points:

| Module | Role |
|--------|------|
| `polya_isomer_count.py` | Burnside–Pólya isomer enumeration (the completeness theorem) |
| `polyhedra.py` | ideal coordination polyhedra + vertex sets |
| `decompose.py` | SMILES → coordination decomposition (donors, denticity, hapticity) |
| `assemble_complex.py` | places donors on vertices, orients chelates, seats the metal |
| `conformer_enum.py` / `conformer_complete.py` | per-isomer conformer construction |
| `_gfnff_rank.py` | GFN2 / GFN-FF energy ranking via the `xtb` CLI |

## Scope & limitations (read this)

MANTA is honest about its boundaries:

- **Starting geometries, not production geometries.** Output is a *constructive* geometry (≈ UFF-quality), correct in topology/coordination but **not** xTB/DFT-accurate. Always relax with xTB/DFT before any production calculation.
- **Recall is gas-phase-bounded.** MANTA constructs the gas-phase-feasible manifold. Crystal-packing-distorted poses are, by construction, not reachable — this is a shared physical ceiling, not an enumeration gap.
- **Build-rate boundaries.** A minority of systems (some strained metallacycles whose distance-geometry embed collapses, certain multi-hapto / cage motifs) do not build FF-free and fall back or are skipped — documented, not hidden.
- **Completeness is over the modelled coordination model** (denticity, hapticity, ideal polyhedra). Exotic bonding outside that model is out of scope.

## License & data provenance

LGPL-3.0-or-later. The runtime uses **only open priors** (Pyykkö covalent radii, open-database-derived ideal distances) — **no CCDC/CSD/GRIP data ships or is read at runtime**. A mechanical license guard (`scripts/license_guard.py` + CI test) enforces this.
