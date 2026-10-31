# PROVENANCE — butane relaxed dihedral scan, r2SCAN-3c (RC1 fixture)

## What this fixture is
A real ORCA 6.1.0 run (not constructed): relaxed surface scan of the central
C1-C2-C3-C4 dihedral of n-butane, 24 points at 15 deg steps covering the full
circle (Scan D 0 1 2 3 = 180, 525, 24; values mod 360 give the 0-360 profile),
started from the GOAT global minimum (anti). All 24 constrained optimizations
converged.

## Job and run facts
- Fixture run: job 8173642 (v5), ORCA "6.1.0 - RELEASE" (banner line in the .out),
  24/24 HURRAY, 0 "did not converge", ORCA TERMINATED NORMALLY. Wall time 23 min
  (relaxscanact grading: operator report; .out TOTAL RUN TIME line agrees).
- Resources: SLURM pal 2 cores, %maxcore 2000 MB, time limit 05:00:00, actual
  ~23 min -> core-hours 2 x 0.39 h = 0.77 core-h. CPU efficiency line: taken from
  the SLURM job accounting (delfin_8173642.out); see the job folder
  .gate/practice/rc1_butane_scan_r2scan3c_v5/ for the raw log.
- Estimate stated before the run: ~1.5 min/step -> 40-60 min; measured 23 min.
- Submitted through DELFIN's submit path by the operator
  (delfin.api.submit_calculation, mode=orca, folder
  .gate/practice/rc1_butane_scan_r2scan3c_v5).

## Physics (from rc1_butane_scan_r2scan3c_v5.relaxscanact.dat, cross-checked
## against the .out per-step FINAL SINGLE POINT ENERGY lines, max diff < 1e-6 Eh)
- anti (180 deg, step 1) = global minimum: -158.39247808 Eh.
- gauche grid points 300/420 deg: -158.39152313 / -158.39152809 Eh
  -> gauche-anti = 0.599 / 0.596 kcal/mol at the grid points. This sits at the
  lower edge of the expected 0.6-0.9 kcal/mol band and is CONSISTENT with the
  literature value ~0.86 kcal/mol (Wiberg & Murcko, JACS 1988, 110, 8029-8038,
  DOI 10.1021/ja00232a012, trans/gauche 0.86 kcal/mol, ZPE-corrected), because a
  15-deg grid point cannot sit exactly in the true gauche minimum between grid
  points and therefore slightly under-reports the gap. No band was tuned for it.
- syn (360 deg, step 13) = global maximum: -158.38354658 Eh -> syn-anti =
  5.605 kcal/mol.
- eclipsed (240 deg): -158.38714903 Eh -> 3.344 kcal/mol above anti.
- Profile symmetric about 360 deg to ~1e-5 Eh (steps k and 26-k agree).
- Ordering holds: syn > eclipsed > gauche > anti at every pair checked.

## Provenance chain (v1 -> v5; every failed run kept, never rerun into an old folder)
- v1 (never submitted): draft with an explicit { D 0 1 2 3 C } constraints block
  next to the Scan block - refused before submission (the Scan block itself
  constrains the coordinate; ORCA 6.1.1 manual section 4.2.1/4.2.6).
- v2, job 8171987: FAILED - ORCA 6.1.1 rejected "NoFreq" as a simple keyword
  ("UNRECOGNIZED OR DUPLICATED KEYWORD(S) ... NOFREQ", .out lines 182-186).
  Exit code 1, nothing salvageable.
- v3, job 8172117: FAILED - scan start 0 deg demanded a 180 deg jump from the
  anti start geometry; "GSTEP: could not impose initial constraints"
  (.out lines ~293-298). ORCA 6.1.1 manual section 4.2.2 Tip: the scan start
  value should not be far from the initial structure.
- v4, job 8172492: FAILED - 24/25 steps ran, but the step at 540 deg (= 180
  mod 360, the anti wrap point) did not converge; the dihedral target sits on
  the +/-180 deg wrap so the optimiser chases a discontinuity; +1.9 Eh energy
  collapse mid-step; SHARK/PROPERTIES abort downstream. Steps 1-24 data intact
  but the fixture rule requires EVERY point converged.
- v5, job 8173642: SUCCESS - scan 180, 525, 24 (24 distinct angles, no
  duplicate end point), %geom MaxIter 100 retained from the v4 diagnosis
  (harmless), everything else identical to v4. THIS run is the fixture.

## Files
- rc1_butane_scan_r2scan3c_v5.inp - ORCA input (verbatim from the run).
- rc1_butane_scan_r2scan3c_v5.out - TRIMMED (the full run .out is 6.0 MB, above
  the 2 MB fixture limit; the trim keeps head/banner/input echo, per step the
  last convergence block + FINAL SINGLE POINT ENERGY, and the trailing scan
  summary + normal-termination lines; every step keeps its convergence block;
  trim markers "[[trimmed: ...]]" mark the cuts; the full .out lives in the
  practice folder .gate/practice/rc1_butane_scan_r2scan3c_v5/, not committed).
- rc1_butane_scan_r2scan3c_v5.property.json - ORCA property output, unmodified.
- rc1_butane_scan_r2scan3c_v5.relaxscanact.dat - ORCA scan table (angle, energy
  in Eh per step), unmodified.
- No .gbw committed (restart files excluded by the package rules).
