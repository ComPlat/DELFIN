# PROVENANCE — rc2_ru_o4

- System: RuO4, ruthenium tetroxide. Target state: singlet (S=0),
  tetrahedral Td, Ru(VIII) d0.
- Method: r2SCAN-3c (RKS), TightOpt + NumFreq, ORCA 6.1.1.
  Input keywords: `! r2SCAN-3c TightOpt NumFreq` with
  `%output jsonpropfile True end`, `%pal nprocs 2 end`,
  `%maxcore 1000`; charge/multiplicity `* xyz 0 1`.
- Job id: 8172021 (submitted by the operator through DELFIN's own
  submit path, mode=orca, 2026-10-31; time limit 00:45:00).
- Wall time vs estimate: job wall-clock 00:04:32 vs my estimate 20 min
  (actual ORCA TOTAL RUN TIME 2 min 14 s, .out line 3551; the estimate
  was generous — NumFreq for 5 atoms is 12 displaced gradients).
- Core-hours: CPU Utilized 00:07:14 (0.12 core-hours).
- CPU efficiency: 39.89% of 00:18:08 core-walltime
  (delfin_8172021.out JOB FEEDBACK line 85).
- Spin-state evidence (computed state and citation in the same
  sentence): the computed state is an RKS closed-shell singlet
  (Multiplicity 1, .out line 3124), which is the only state possible
  for RuO4 because Ru is in the rare +8 oxidation state — a d0 metal —
  so the valence shell holds zero unpaired electrons and no alternative
  spin state exists for this electron count; the molecular-oxide
  character and +8 state of RuO4 are described in
  doi.org/10.1021/acs.jpcc.4c08816 (JPCC 2025), verified through
  DELFIN's web tools this session. A higher spin state is not merely
  implausible, it is impossible for d0.
- Relativistic treatment (what DELFIN chose and why):
  DELFIN's own rule `_should_use_rel` (delfin/utils.py:203-217)
  classifies Ru as a 4d element for which a relativistic treatment is
  suggested, and its set_main_basisset path (utils.py:254-260) would
  attach a ZORA-def2-SVP main basis. That rule was NOT applied here
  because this run went through ORCA's composite keyword r2SCAN-3c,
  which fixes its own basis and dispersion (api.py:2518-2536
  documents a composite as self-contained) and for Ru uses the
  def2-ECP/def2-mTZVPP pair: the .out confirms "Group 1, Type Ru ECP
  Def2-ECP (replacing 28 core electrons, lmax=3)" (line 381) and the
  basis def2-mTZVPP (line 209). Scalar relativity is thus carried by
  the energy-consistent def2 ECP itself (quasi-relativistic
  derivation, Andrae et al., Theor. Chim. Acta 77, 123-141 (1990),
  doi.org/10.1007/BF01114537, for the Y-Cd/Hf-Hg ecp-28/46 families
  Ru falls in; the citation is printed by ORCA itself in the .out at
  line 195). Stacking ZORA on top of an ECP basis would double-count
  relativistic effects, so no ZORA keyword appears in the .inp —
  verified: the .inp contains no %basis block and no ZORA.
- Physics gates (from rc2_ru_o4.out): SCF converged, geometry
  converged, TERMINATED NORMALLY, ZERO imaginary modes (grep count 0
  of "imaginary mode"), Gibbs G = -395.79978 Eh, ZPE = 0.01268 Eh,
  single point -395.78564 Eh.
- Geometry: perfect tetrahedron, all four Ru-O = 1.701 A
  (sqrt(3) x 0.9820 A per final .xyz; O-Ru-O = 109.47 deg within
  numerical precision; Td kept).
- Files committed as fixture: .inp, .out (trimmed only above 2 MB),
  .gbw, .json/.property.json, PROVENANCE.md.
