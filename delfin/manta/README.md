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

**MANTA** is the deterministic constructor inside [DELFIN](../../README.md). From a transition-metal **SMILES** it **enumerates** the coordination isomers as symmetry orbits and **builds** each of them with the metal held on its coordination polyhedron, adding ring-pucker, pose and mirror siblings. The result is the **isomer–conformer ensemble** of the input: the finite set of frames emitted for one SMILES (by analogy with the conformer–rotamer ensemble of CREST); a single structure is a frame. On this path **no force field acts on a free metal centre**, and **no licensed crystallographic data is read at runtime**.

> **Construct, don't search.** CREST / GOAT explore the conformational space *stochastically* (metadynamics, genetic search). MANTA *constructs* the ensemble from group theory and ideal geometry. The two are complementary: MANTA gives a deterministic, license-clean set of starting frames, with the enumerated isomer count as an explicit denominator; a global optimizer can refine any of them.

## What it does

| Stage | How |
|-------|-----|
| **Coordination isomers** | The admissible placements of the donor arms on the vertices of an ideal polyhedron are reduced to orbits under the rotation group of the polyhedron and the interchange of identical ligands (Burnside–Pólya counting); mirror images are distinct orbits. The enumerator agrees with the orbit count and with brute force for monodentate ligands and for bidentates whose donors differ in element or are symmetry-equivalent (`tests/test_burnside_enumerator.py`, `tests/test_prescribed_isomer_enumerator.py`); for bidentates with two different donors of one element it can undercount, and for tri- and higher-dentate ligands it can under- or overcount |
| **Geometry** | Donors placed on the ideal vertices at covalent-radius metal–donor distances (Cordero; Pyykkö–Atsumi for heavy main-group donors), each ligand turned so that its lone pair points at the metal; a chelate is embedded as a metallacycle and rotated onto its vertices; the ligand periphery is then relieved with metal and donors frozen |
| **Siblings** | Ring puckers (Cremer–Pople targets, relaxed with UFF while metal and donors are fixed points), pose, lone-pair and re-embedding siblings and a mirror frame, each gated and de-duplicated. The ensemble is a union of these sets, not a product; only the isomer part has a closed-form count, the sibling sets are heuristic and capped, and no torsion grid is enumerated |
| **Self-gate** | Every frame must pass a check against bond collapse, atom overlap, shell intrusion and a continuous-shape limit before it is emitted |
| **Ranking** | Optional GFN2-xTB (or GFN-FF) energy ranking of the emitted ensemble (via the `xtb` CLI); off by default |

Inputs outside the constructor's scope (it covers one metal centre with coordination number 4, 5 or 6 and ligands below a size cap per donor arm), or whose isomers all fail the self-gate, are built entirely by a **legacy path** that minimises UFF over the whole complex with the metal as a typed, fixed atom. The two paths are either-or per system; their outputs are not combined.

Construction uses fixed seeds and total orders in enumeration, selection and de-duplication, and the champion configuration also seeds the random nudges Open Babel's force fields apply at near-zero contacts (previously drawn from a clock-seeded generator), so the same SMILES on the same commit and environment gives **byte-identical frames**. Across machines, byte identity also needs the same C math library: different glibc versions round a few transcendental functions differently in the last digit at rare arguments, which can change a de-duplication or ranking decision; the [reference environment](../../README.md#reference-environment-for-reproducible-manta-builds) names the one we build in.

## Usage

MANTA has its own command, `delfin-manta`, and is reachable from a DELFIN run
with `smiles_converter=MANTA` in CONTROL.txt. From Python:

```python
from delfin.smiles_converter import smiles_to_xyz_isomers

# returns the isomer–conformer ensemble (one XYZ block per frame)
isomers = smiles_to_xyz_isomers("[Co](N)(N)(N)(N)(Cl)Cl")
```

Optional GFN2-xTB energy ranking of the ensemble (requires `xtb` on `PATH`):

```bash
DELFIN_FFFREE_GFNFF_RANK=1 DELFIN_CONF_RANK_METHOD=gfn2 python your_script.py
```

All new behaviour is **flag-gated and default-OFF** (byte-identical to the prior pipeline unless explicitly enabled).

## Package layout

`delfin/manta/` holds the constructor (about 120 modules). Key entry points:

| Module | Role |
|--------|------|
| `polya_isomer_count.py` | Burnside–Pólya orbit count of the coordination isomers |
| `polyhedra.py` | ideal coordination polyhedra + vertex sets |
| `decompose.py` | SMILES → coordination decomposition (donors, denticity, hapticity) |
| `assemble_complex.py` | places donors on vertices, orients chelates, seats the metal |
| `_ring_pucker.py` | ring-pucker siblings (Cremer–Pople targets) |
| `conformer_enum.py` / `conformer_complete.py` | ligand-conformer enumeration and the conformer-completion pass; not used by the shipped configuration (the completion pass is behind a switch that is off) |
| `_gfnff_rank.py` | GFN2 / GFN-FF energy ranking via the `xtb` CLI |

## Scope & limitations (read this)

MANTA is honest about its boundaries:

- **Starting geometries, not production geometries.** Output is a *constructive* geometry (≈ UFF-quality), correct in topology/coordination but **not** xTB/DFT-accurate. Always relax with xTB/DFT before any production calculation.
- **Recall is gas-phase-bounded.** MANTA constructs a gas-phase ensemble. Crystal-packing-distorted poses are, by construction, not reachable — this is a shared physical ceiling, not an enumeration gap.
- **Build-rate boundaries.** A minority of systems (some strained metallacycles whose distance-geometry embed collapses, certain multi-hapto / cage motifs) do not build FF-free and fall back or are skipped — documented, not hidden.
- **Enumeration is over the modelled coordination model** (denticity, hapticity, ideal polyhedra). Only the enumerated orbits that build and pass the self-gate are realised, so the realised isomers are bounded by the enumeration, not guaranteed equal to it; how much of the enumerated set a run realises is measured, not asserted. Exotic bonding outside that model is out of scope.

## License & data provenance

LGPL-3.0-or-later.

MANTA itself ships **no CCDC/CSD structures, coordinates or fragment libraries and no CSD-derived tables** and never calls the CSD API. Covalent radii are Cordero et al. (2008), with the Pyykkö–Atsumi radii for multiple bonds and heavy main-group donors; the five delocalised aromatic bond lengths (`delfin/fffree/aromatic_bond_targets.py`) and the donor-plane tolerance `_BETA_BAND` (5.2°, `assemble_donor_plane.py`) are our own statistics over the Crystallography Open Database (COD, CC0); further COD bond-length medians serve an optional refinement path that is off by default; the metal–ligand distance tables of the legacy path are hand-set typical values, not a database statistic. Every value and its source: [`PROVENANCE.md`](PROVENANCE.md).

A few construction paths that are off in the shipped configuration were calibrated on structure statistics that are not distributed; they read their values from tables named by environment variables and have no value without them. Third-party libraries carry published structure-derived knowledge of their own: RDKit's ETKDG, for example, uses experimental torsion preferences derived from the CSD (Riniker & Landrum, 2015).

Where every non-trivial number comes from is listed in [PROVENANCE.md](PROVENANCE.md). A mechanical license guard (`scripts/license_guard.py` + CI test) blocks CCDC/CSD imports and data files in the tree.
