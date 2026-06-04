<<<<<<< HEAD
"""DELFIN fffree — FF-free TMC structure generation pipeline.

Universal, deterministic, mathematically complete 3D structure generation
for transition metal complexes from SMILES strings. No force field is used
in the construction pipeline.

Core architecture (post 2026-05-30/31 night-sprint, 15 new modules):

  INPUT (SMILES)
  │
  ├─ oxidation_state    : metal + oxidation state inference
  ├─ spin_states        : multiplicity enumeration per d-count
  ├─ linkage_isomers    : ambidentate κ-binding choices
  ├─ solvent_counterion : solvent/counterion strip
  │
  ├─ polya_isomer_count : Pólya-Burnside isomer enumeration
  ├─ decompose          : donor classification (σ/π/hapto/μ-bridging)
  │
  ├─ polyhedra          : CN3-9 coordination polyhedra
  ├─ multi_metal_polyhedra : M2/M3/M4/M6 clusters
  ├─ fblock_polyhedra   : CN7-12 for Ln/An
  ├─ metal_sphere_builder  : vertex placement
  │
  ├─ conformer_enum     : ETKDG + torsion grid
  ├─ ring_pucker        : Cremer-Pople universal N>=4
  ├─ ring_pucker_integration : RDKit Mol -> variants
  ├─ lema_saavedra      : 2-stage robust ring construction
  ├─ macrocycle         : porphyrin NSD + calixarene
  ├─ backbone_torsion   : coupled gauche-anti pairs
  ├─ bridging_ligand    : mu2/mu3/mu4 placement
  ├─ assemble_complex   : orient + place donor on vertex
  │
  ├─ refine             : FF-free defect-loss gradient descent
  ├─ donor_slide        : Thomson-sphere donor relaxation
  ├─ cod_ideals         : COD-empirical bond targets
  │
  ├─ symmetry           : point-group detection
  └─ burnside           : orbit-counting completeness proof

  OUTPUT: XYZ x (isomer x conformer x linkage x spin) x multiplicity

Mathematical foundations cited:
  - Cremer & Pople JACS 97, 1354 (1975)              -- ring puckering
  - Lema-Saavedra & Fernandez-Ramos JCP 164, 074109 (2026)  -- 2-stage CP
  - Burnside 1897 / Cauchy-Frobenius                  -- orbit counting
  - Jentzen, Shelnutt, Smith JACS 1997               -- NSD porphyrin basis
  - Gutsche 1989                                      -- calixarene conformations
  - Pophristic & Goodman Nature 411, 565 (2001)      -- pentane-effect
  - Bethe 1929 / Van Vleck                            -- crystal field theory
  - Riniker & Landrum JCIM 2015                       -- ETKDG distance geometry

Env-flag DELFIN_FFFREE_BUILDER=1 enables the FF-free path in smiles_converter.
Env-flag DELFIN_FFFREE_PURE_TRACK3=1 enables the FULL universal FF-free
architecture (auto-enables RING_PUCKER, SYMMETRY, BURNSIDE, MACROCYCLE,
MULTI_METAL, FBLOCK, SPIN_STATES, LINKAGE, SOLVENT, BRIDGING, BACKBONE,
OXSTATE, COD_BONDS, RIGID_H_DRAG). Default OFF, byte-identical when unset.
"""

__version__ = "0.3.0"  # Bumped for 2026-05-30/31 night-sprint architecture
