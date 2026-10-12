# Where MANTA's numbers come from

The numerical values the MANTA constructor uses fall into six classes. This file names the class and the source of every table and constant that sets a geometry; the many gate thresholds are summarised as one row.

| Class | Meaning |
|---|---|
| **Literature** | Published values, cited below |
| **COD** | Our own statistics on the Crystallography Open Database (CC0; Gražulis et al., *Nucleic Acids Res.* 40, D420, 2012) |
| **Theory** | Ideal geometry: VSEPR angles, regular polyhedra and polygons |
| **Hand-set** | Chosen by hand: round thresholds, typical distances from textbook ranges, some adjusted to individual published structures. Several thresholds were afterwards *checked* against crystal structures (for example, that a gate stays silent on clean structures); checking does not make them derived values |
| **CSD-derived, to be replaced** | Our own statistics on CSD crystal structures. Exactly two values of the shipped configuration are in this class; see below |
| **Not shipped** | Values used by paths that are off in the shipped configuration and that were calibrated on structure statistics we do not distribute. The code reads them from a file named by an environment variable and has no value without it |

## Values of the shipped configuration

| Value | Where | Class and source |
|---|---|---|
| Covalent radii (single bond) | `ml_tables._COVALENT_RADII`, `polyhedra.COV`, `_elements.COV_R` and their copies | **Literature**: Cordero et al., *Dalton Trans.* 2008, 2832. Low-spin values for Mn/Fe/Co in `ml_tables` and `_elements`; the mean of low and high spin for Mn/Fe/Co in `polyhedra.COV` |
| Covalent radii for double and triple bonds and for heavy main-group donors | `polyhedra.PYYKKO_SINGLE/DOUBLE/TRIPLE`, `polyhedra._METALLOID_COV`, `_energy_terms`/`uff_constraints._pyykko_order_radius` (via mendeleev) | **Literature**: Pyykkö & Atsumi, *Chem. Eur. J.* 15, 186 and 12770 (2009); Pyykkö, Riedel & Patzschke, *Chem. Eur. J.* 11, 3511 (2005). `_METALLOID_COV` mixes Pyykkö and Cordero values. Known deviation: the Te double-bond radius stands at 1.30 Å, the published value is 1.28 Å |
| Bond-order scale factors 1.00 / 0.93 / 0.87 / 0.78 | `_energy_terms._BOND_ORDER_SCALE` | **Hand-set** (fallback where no order-specific radius exists) |
| Ideal coordination polyhedra, VSEPR angles (109.47°, 120°, 180°), regular-polygon ring angles | `polyhedra`, `_energy_terms`, `uff_constraints` | **Theory**. The SPY-5 basal offset (0.2, or −0.3135 under its switch) is **hand-set** to reproduce the textbook angles of C4v square pyramids |
| Legacy-path metal–ligand distances, M–H, As/Sb/Te/B/Se donors, M–M and metal–centroid distances, row offsets | `ml_tables._METAL_LIGAND_BOND_LENGTHS`, `_METAL_HYDRIDE_BOND_LENGTHS`, `_METAL_LIGAND_BOND_LENGTHS_T33`, `_METAL_METAL_BOND_LENGTHS`, `_HAPTO_CENTROID_DISTANCES`, `_METAL_ROW_OFFSET` | **Hand-set**: typical distances from textbook ranges and Cordero radii sums, some adjusted to individual published structures (for example Ir(ppy)2(acac), [Rh(ppy)2Cl(MeCN)], ferrocene, uranocene). They were set before the CSD statistics used elsewhere in this project existed and are not database averages. The M–H and As/Sb/Te/B/Se tables sit behind switches that are off |
| van der Waals radii | `_vdw_radii.py` and copies | **Literature** for the common main-group elements: Bondi, *J. Phys. Chem.* 68, 441 (1964). The values for main-group elements Bondi lacks and all metal values are of unverified origin and count as **hand-set**; the metal values do not reproduce Alvarez (2013), which the module header cites |
| Gate and repair thresholds (collapse floors 0.82 / 0.78 / 0.70, the X–H floor 0.55, H···H floor 1.50 Å, clash factors 0.65 / 0.85, tetra-volume 0.40 Å³, bond floor 0.65, aromaticity cut 1.46 Å and similar) | `_bond_decollapse`, `converter_backend`, `converter_flags`, `_h_placement`, `hapto_declash`, `sphere_flex`, `_arom_criterion` and others | **Hand-set**, several checked against clean crystal structures for silence |
| Interpolation of aromatic radii (2/3 of the way from single to double) | `_arom_bond_length._AROM_FRACTION` | **Hand-set** on Literature radii |
| **Delocalised aromatic bond lengths** C–C, C–N, C–O, C–S, N–N | `delfin/fffree/aromatic_bond_targets.py` (`AROMATIC_BOND_TARGETS`), used by the aromatic-bond seat | **CSD-derived, to be replaced**: medians of our own statistics on CSD crystal structures (bond order "aromatic", ligand bonds in metal complexes). The same method applied to the COD reproduces them to within 0.012 Å, C–C identically; the COD values will replace them |
| **Donor-plane tolerance** 4.8° | `assemble_donor_plane._BETA_BAND`, used by the ring-pucker, lone-pair and β siblings | **CSD-derived, to be replaced**: the upper out-of-plane angle of the metal at conjugated donors in clean CSD crystal structures. It will be re-measured on the COD and replaced |

## Values used only by paths that are off in the shipped configuration

| Value | Where | Class |
|---|---|---|
| COD bond-length medians per element pair | `cod_ideals._COD_P50` (switch `DELFIN_FFFREE_COD_BONDS`) | **COD** |
| Per-metal M–C target lengths | `ml_tables._mc_len_crystal_median` (switch `DELFIN_FFFREE_MC_LEN_CRYSTAL`) | **Not shipped**: `DELFIN_FFFREE_MC_LEN_TABLE` |
| Per-period donor angles, measured 5- and 6-ring angles, M–D–X slope | `_energy_terms._period_theta`, `_ring_internal_deg`, `_ring_exocyclic_deg`, `_phi_md_slope`, used by `_energy_terms._signature_theta` (switch `DELFIN_FFREE_THETA_SIGNATURE`) and `uff_constraints._donor_sigma_geometry` (switch `DELFIN_FFFREE_METAL_SIGMA_COUNT`) | **Not shipped**: `DELFIN_FFREE_THETA_TABLE`. Without it the tetrahedral angle, the regular polygon and a zero slope apply |
| Element-normalised aromaticity threshold | `_arom_criterion._arom_ratio_max` (switch `DELFIN_FFFREE_AROM_CRITERION_RADII`) | **Not shipped**: `DELFIN_FFFREE_AROM_RATIO_TABLE`. Without it the absolute 1.46 Å gate decides |
| Terminal M=E multiple-bond lengths | `ml_tables._ml_me_band` | **Not shipped**: `DELFIN_ML_ME_BANDS` |
| Measured bond, angle, metal–donor and donor–metal–donor bands, inverse-variance weights | `polyhedra`, `_energy_terms`, `assemble_chelate` | **Not shipped**: `DELFIN_FFREE_BOND_BANDS`, `DELFIN_FFREE_ANGLE_BANDS`, `DELFIN_FFREE_MD_BANDS`, `DELFIN_FFREE_DMD_BANDS`, `DELFIN_FFREE_INVVAR` |
| Shape budgets for the ligand swing | `joint_declash` | **Not shipped**: `DELFIN_FFFREE_LIGAND_SWING_CSHM`, `DELFIN_FFFREE_LIGAND_SWING_MAX_START_CSHM` (default 0 = off) |

The table files share one simple format (`key<TAB>value` for the scalar tables, see `_unshipped_tables.py`; the band tables are described in `polyhedra._load_band_tsv`). The shipped configuration sets none of these variables.

## Third-party software

Libraries MANTA calls carry their own parameters, some of them derived from structure databases:

- **RDKit** ETKDG conformer embedding uses experimental torsion preferences derived from the CSD (Riniker & Landrum, *J. Chem. Inf. Model.* 55, 2562, 2015). UFF parameters: Rappé et al., *J. Am. Chem. Soc.* 114, 10024 (1992).
- **Open Babel** builds 3D coordinates from its own fragment and ring templates and supplies UFF for the legacy path.
- **xtb** (GFN2-xTB / GFN-FF), used only for the optional energy ranking.

## Identifiers in comments and tests

Comments and tests name individual structures by their CSD refcodes, as public identifiers of the examples a fix was found on. No measured value of an entry and no coordinates are stored with them.
