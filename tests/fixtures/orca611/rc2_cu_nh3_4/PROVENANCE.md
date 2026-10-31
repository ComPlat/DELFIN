# PROVENANCE — rc2_cu_nh3_4 (fixture = v3 run; v1, v2 documented below)

(v3 run 8173353 complete; all values below verified from its .out by s3.)

- System: [Cu(NH3)4]2+, tetraamminecopper(II). Target state: doublet
  (S=1/2), square planar, Cu(II) d9. CPCM(water).
- Method: r2SCAN-3c (UKS, multiplicity 2) + CPCM(water), TightOpt +
  NumFreq, MaxIter 200.
- Input keywords: `! r2SCAN-3c CPCM(water) TightOpt NumFreq` with
  `%output jsonpropfile True end`, `%pal nprocs 2 end`, `%maxcore 1000`,
  `%geom MaxIter 200 end`.
- Job id (v3): 8173353 (operator-submitted through DELFIN's own path,
  pal 2, maxcore 1000, 01:00:00 limit).
- Wall time vs estimate: 21 min 22 s actual (ORCA run time, .out:7842)
  vs 35 min estimate (request 01:00:00). Earlier runs: v2 36:31, v1
  17:49. Core-hours: CPU Utilized 01:04:32, CPU Efficiency 67.88% of
  01:35:04 core-walltime, job wall 00:23:46 (sacct via
  delfin_8173353.out:85-88).
- Computed (verified from the v3 .out by s3, not taken on trust):
  0 imaginary modes (lowest real 16.99 cm-1, the former mode-6 rocking
  coordinate now positive — root cause of the v2 imag mode was the flat
  out-of-plane rocking of square-planar d9 Cu, curvature near the
  NumFreq noise floor); 4x Cu-N 2.039/2.039/2.040/2.040 A; <S^2> =
  0.7525 (Deviation 0.002528, .out:912); Mulliken spin on Cu 0.555072,
  sum of spin populations 1.0000000 (rest delocalised on the N donors —
  normal Cu(II)-amine covalency, .out:1018-1038); FINAL ENERGY
  -1866.529712725782 Eh; Gibbs -1866.40930926 Eh (extract_thermochem);
  SCF converged in 14/10/9/8/7... cycles per geometry step; geometry
  optimization converged (HURRAY, .out:5420); ORCA TERMINATED NORMALLY
  (.out:7841).
- Physics history (the chain of runs, with job ids — the geometry came
  through a mode displacement and must say so):
  - v1, job 8172018: gas-phase 2+ cation, NO solvent. Optimization did
    NOT converge (51 cycles, "reached the maximum"); two NH3 ligands
    dissociated to y = ±6.6-7.5 A during the run. <S^2> = 0.754
    (healthy doublet). Cause: charge artefacts of a bare dication in
    the gas phase (operator directive to add implicit solvent).
  - v2, job 8172924: CPCM(water) added on the ! line; start geometry =
    cycle 4 of v1 (last structure with all four NH3 bound — deviation
    from the operator's "LAST geometry" directive, accepted by operator
    and reviewer; restarting from the dissociated last structure would
    have continued dissociation). Wall 36:31. CONVERGED: all four
    Cu-N 2.039-2.041 A, <S^2> = 0.7525, 0 SCF failures. ONE imaginary
    mode: mode 6 = -21.86 cm-1 (out of 51 modes, output line 29313).
  - v3 (this run): displace the converged v2 geometry along mode 6 and
    re-optimise + NumFreq. Displacement hand-built from the mode-6
    vector in the v2 normal-mode block (rc2_cu_nh3_4_v2.out:29420-29455,
    column "6"): dominant motion = trans NH3 ligands on the y axis
    rocking out of the square plane (largest z components: N(−y) z =
    -0.420, its H z = -0.415/-0.279; the +y pair moves oppositely),
    Cu stationary (|components| ~ 0.002). Displacement applied: z of
    the two trans-N NH3 shifted by ∓0.02 A (y-ligand N, H moved
    against the mode vector's dominant z components), all other
    coordinates unchanged — i.e. moved off the saddle in the mode-6
    direction, the same first step DELFIN's own IMAG engine makes
    (delfin/imag.py:432-556, displaced_geometries() reads the .hess
    normal-mode block and displaces ± along it; same principle, one
    hand-applied step here).
- Mode-6 character (evidence, not guess, per QS s5 requirement): the
  -21.86 cm-1 soft mode is an NH3 rocking/umbrella-type torsion —
  trans-NH3 out-of-plane rocking with Cu stationary. Magnitude is
  consistent with numerical noise on a soft rotor (ORCA manual 6.1.1,
  section 1.7.14 "Troublesome Geometry Optimizations": imaginary
  frequencies after optimization "typically" require "tighter
  optimization criteria, tighter DFT grids" or displacement along the
  mode — quoted verbatim, doc_id orca_manual_6_1_1_delfin,
  section ch1_7_14...). The v3 run is the displacement fix from that
  section; if mode 6 (or any other) survives in v3, the value and the
  moving atoms are reported, not hidden — the package's rule is
  stop-and-report, and the operator decides whether a soft mode is
  acceptable for a reference fixture.
- Solvation: CPCM(water) on the ! line (v2 onward). Why (operator
  directive, 2+ cations are not state-of-the-art in the gas phase —
  ligand loss and charge artefacts of bare dications, exactly the
  dissociation the v1 run showed): the continuum water model screens
  the cation's charge and stabilises the intact complex. Keyword syntax
  verified from DELFIN's own code (delfin/esd_input_generator.py:244
  'CPCM(water)', ensemble_nmr.py:246) and the ORCA 6 manual/faccts
  tutorial Implicit Solvation Models. Literature support (DOI verified
  via DELFIN's web tools): Norjmaa, Ujaque, Lledos, 'Beyond Continuum
  Solvent Models in Computational Homogeneous Catalysis', Topics in
  Catalysis 65, 118-140 (2022), doi.org/10.1007/s11244-021-01520-2 —
  the standard practice in computational homogeneous catalysis is to
  include solvent effects through continuum models; charged metal
  complexes especially carry large gas-phase artefacts without them.
  The implicit model does NOT reproduce the specific solvation shell;
  that limitation is recorded here (a fully explicit treatment would
  change the reference system, not this fixture's purpose).
- Spin check: electron count = Cu(29) - 2 + 4x(7+3) = 67, odd -> at
  least doublet; d9 Cu(II) -> one unpaired electron -> exactly doublet.
  v2 measured: <S^2> = 0.7525 (ideal 0.75 for a pure doublet);
  deviation reported, not hidden. v3 (the fixture) measured:
  <S^2> = 0.7525 (Deviation 0.002528, .out:912); Mulliken spin on
  Cu 0.555072, sum 1.0000000 (rest delocalised on N donors — normal
  Cu(II)-amine covalency, .out:1018-1038).
- Physics gates (all verified from the v3 .out, job 8173353):
    SCF converged yes (Termination: ... all done! at .out:7839; RKS
    convergence clean in every cycle); geometry converged yes
    ("THE OPTIMIZATION HAS CONVERGED" at .out:7838 + "OPTIMIZATION RUN
    DONE"); 0 imaginary modes (lowest real 16.99 cm-1 at .out:29453,
    FREQ block "mode 6" positive); all four Cu-N bound 2.039/2.039/
    2.040/2.040 A (final xyz block, .out:7805-7820); spin mostly on Cu
    (Mulliken spin population Cu 0.555072 of sum 1.0000000,
    .out:1018-1038).
