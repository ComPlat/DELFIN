# PROVENANCE — rc2_fe_co5

- System: Fe(CO)5, iron pentacarbonyl. Target state: singlet (S=0),
  trigonal-bipyramidal D3h. 18-electron Fe(0) d8 low-spin closed shell.
- Method: r2SCAN-3c (r2SCAN meta-GGA + D4 + gCP, def2-mTZVPP basis),
  TightOpt + NumFreq. ORCA 6.1.1. Input keywords: `! r2SCAN-3c TightOpt
  NumFreq` with `%output jsonpropfile True end`, `%pal nprocs 2 end`,
  `%maxcore 1000`; charge/multiplicity `* xyz 0 1`.
- Job id: 8172015 (submitted by the operator through DELFIN's own
  submit path, mode=orca, 2026-10-31; time limit 01:30:00).
- Wall time vs estimate: job wall-clock 00:15:44 vs my estimate 30 min
  (the operator raised the limit from my 30 min to 01:30:00 because
  NumFreq needs 66 displaced gradients — the estimate was too generous
  for this run; actual ORCA TOTAL RUN TIME 13 min 26 s, .out line 5352).
- Core-hours: CPU Utilized 00:40:38 (0.68 core-hours).
- CPU efficiency: 64.57% of 01:02:56 core-walltime (delfin_8172015.out
  JOB FEEDBACK line 85). Memory: 454.81 MB, 22.74% of 1.95 GB.
- Spin-state evidence (computed state and citation in the same
  sentence): the computed state is an RKS closed-shell singlet
  (Multiplicity 1 throughout the .out, no UHF/S2 block, SCF converged),
  which is the state expected for Fe(CO)5 because Fe(0) is d8 and five
  2-electron CO donors close the shell at the 18-electron count, the
  ground-state structure of gas-phase Fe(CO)5 being the D3h
  trigonal bipyramid as reviewed by Portius and Buhl,
  doi.org/10.1021/acs.organomet.9b00559 (Organometallics 2019); its
  triplet excited-state chemistry belongs to 266 nm photochemistry,
  doi.org/10.1021/jacs.4c07523 (JACS 2024), not to the ground state.
  A higher spin state is therefore not plausible for the ground state;
  both DOIs verified through DELFIN's web tools this session.
- Relativistic treatment: none (3d metal — DELFIN's `_should_use_rel`,
  delfin/utils.py:203-217, is False for Fe; r2SCAN-3c's def2-mTZVPP is
  all-electron for Fe; the .out contains no ECP block).
- Physics gates (from rc2_fe_co5.out): SCF converged, geometry
  converged, TERMINATED NORMALLY, ZERO imaginary modes (grep count 0 of
  "imaginary mode"; NumFreq ran to completion), Gibbs G = -1830.41077 Eh
  (final), ZPE = 0.04261 Eh, single point -1830.41644 Eh.
- Geometry: optimized Fe-C_ax = 1.8017 A, Fe-C_eq = 1.7976 A,
  Fe-C-O linear, equatorial C-Fe-C 120 deg — D3h symmetry kept to
  numerical precision (rc2_fe_co5.xyz, final ORCA coordinates, energy
  -1830.416435 Eh).
- Files committed as fixture: .inp, .out (trimmed only above 2 MB;
  geometry, SCF/convergence, spin, frequency and final-energy blocks
  kept), .gbw, .json/.property.json, PROVENANCE.md.
