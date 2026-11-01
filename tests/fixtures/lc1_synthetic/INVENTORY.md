# LC1 Inventory — licence-safe test data (nacht-s1, 2026-11-01)

Task: no ORCA-generated data in the public repository (EULA_ORCA_2025 §3 c–g:
DATA may be shared only for academic purposes; database sharing only
non-commercial with EULA notice; software using DATA must be under the EULA;
ML on DATA only academic). DELFIN is a public LGPL-3.0 repo — ORCA output
files must not be committed.

## Classification evidence

Markers searched over the whole tracked `tests/fixtures/` tree:
`ORCA TERMINATED NORMALLY`, `* O   R   C   A *` banner, `Program Version 6`,
`Host name:`, `CONSTRUCTED`, `SYNTHETIC`, plus `git log -1 -- <group>` per
file group (origin commit given per row below).

Classification used only file/git evidence, never a test's docstring
(reviewer nacht-s2 pre-registered this; confirmed: the
`test_a_method_reviewer_checks_the_method.py:27` docstring says "synthetic"
but the file content is what decides).

## (A) ORCA-generated — must leave the repo (after operator confirms step 2)

Tree `tests/fixtures/orca611/` — 46 tracked files. Real ORCA 6.1.x output
(trimmed, banner host replaced by `fixture-host` in commit d7f08a958) except
the marked subset below:

| fixture | files | evidence | git origin |
|---|---|---|---|
| h2o_optfreq | .out/.json/.property.json/.property.txt/.inp | banner :2, `Program Version 6.1.1` :54, `Host name:` :175; **.gbw 1,010 KB binary** (unparsable, cannot prove provenance — treated as A) | d7f08a958 (host anonymisation) |
| h2o_tddft | .out/.json/.property.json/.property.txt/.inp/.gbw | same banner/host lines | d7f08a958 |
| o2_triplet | .out/.json/.property.json/.property.txt/.inp/.gbw | same | d7f08a958 |
| scf_not_converged | all 6 incl. .gbw | same | d7f08a958 |
| rc1_butane_scan_r2scan3c | .out/.inp/.property.json/.relaxscanact.dat | real cluster run per PROVENANCE.md (job 8172015, 27k lines) | 824df5ed1 |
| rc2_cu_nh3_4, rc2_fe_co5, rc2_ru_o4 | .out/.inp/.property.json each | real cluster runs per their PROVENANCE.md | 6d6fe645b |
| sh1_cu_nh3_imagmode | .out/.inp/.hess | real wave-23 failure copy per PROVENANCE.md (banner + host line) | 063a6b558 |
| sh1_gstep_constraints, sh1_input_error_nofreq, sh1_scan_maxcycles_shark | .out each | banner + host line, wave-23 failure copies | 063a6b558 |

Tree `tests/fixtures/orca_experience/` — 8 files, run1–4 .out/.inp:
RECLASSIFIED A -> B after operator review. Evidence per file (all four
.out added NEW in origin commit 1c10a5ad7, `git diff 1c10a5ad7~1
1c10a5ad7` shows `new file mode` for all 8; no host/user/Working-dir
banner block, which every real ORCA 6.x output carries — compare
orca611 where Host name: is at :175, Working dir. at :177):
- run1.out (26 lines pre-marker): wall time "123.45s 6.0" — synthetic
  decimal, and the trailing "6.0" is not valid ORCA wall-time format.
- run2.out (23 lines): wall 456.78s synthetic; "SCF CONVERGED AFTER 5
  CYCLES" immediately followed by "SCF NOT CONVERGED AFTER 124 CYCLES"
  — physically impossible sequence in one real run.
- run3.out (18 lines): wall 900.00s exactly round; geometry cycles 1
  then 7 with nothing between.
- run4.out (13 lines): wall 50.00s round; no sections at all.
Conclusion: constructed, not cut from a real run. They now carry the
SYNTHETIC first-line marker (guard rule family 2 exception).
Git origin: 1c10a5ad7.

## (B) hand-written / synthetic — stays

- `orca611/butane_dlpno_ccsdt_pvtz` (.out/.json/.property.json/.inp):
  first lines say `Starting time: CONSTRUCTED fixture, no real run`,
  `ORCA terminated normally (CONSTRUCTED fixture -- replace by a cluster
  run)`. BUT it carries the full `* O R C A *` banner and
  `****ORCA TERMINATED NORMALLY****` (:43) — it is hand-written in ORCA
  format. Legal class B, licence-safe; under my new guard rule family it
  would be flagged by the ORCA-signature scan and needs the SYNTHETIC
  first-line marker added (planned, not a data problem; parser impact must
  be gated — see reviewer note, orca611 fixtures are read by
  test_opi_orca_dlpno via parse_orca_output).
  Git origin: 5b6cb3f81.
- `orca611/butane_goat_gfn2` (.out 3022 lines, .json, .property.json,
  .inp, 2 .xyz): mixed — header + most sections look real ORCA, but :2996
  `CONSTRUCTED SECTION (the real run was killed before it printed this)` and
  :3018 `CONSTRUCTED normally-`. .json says `Date: CONSTRUCTED fixture`.
  Classified A-for-removal (mixed provenance, indistinguishable real
  sections present; .finalensemble.xyz/.globalminimum.xyz written by
  ORCA/xTB tooling). Safe: move to private store. Git origin: 5b6cb3f81.
- `method_reviewer/water_b3lyp`, `water_pbe0` (.out only, no .inp): ORCA
  6.0.1-shaped, real wall time `0 days 0 hours 0 minutes 41 seconds 212
  msec` (:59) — real-run signature, NO synthetic marker in file. Docstring
  claims synthetic; evidence does not. Classified A: move to private store,
  replace with synthetic fixtures (parser needs are small — the method
  reviewer checks banner/method lines only). Git origin: 87869f912.
- `science_workspace/run_a..d.out`: xtb-shaped synthetic ("xtb output
  (fixture)" header, no xtb real output format markers). xtb is LGPL;
  synthetic anyway. Class B, stays. Git origin: see `git log -1 --
  tests/fixtures/science_workspace` (same workspace commit family as
  behavior_workspace).
- `behavior_workspace/xtb_sample.out`: xtb synthetic. Class B stays.
- `energy_profile/xtb_*.out.gz` (6 files): xtb outputs, gz. Need one check
  each for ORCA banner before final ruling; listed provisionally B pending
  `zcat`-style inspection via a gate probe.
- `orca_inputs/**` (.inp only, 26 files): ORCA *input* files are hand-written
  by DELFIN/users, not ORCA-generated DATA. Class B. EULA §3 speaks of
  "data generated with SOFTWARE" (§3 c, d) — an input file is written by
  the user into ORCA, not generated by it, so §3 does not restrict it;
  inputs stay.
- control_templates, bench_runs, choice_offer_corpus, benchmark_baseline,
  chemistry/user_project workspaces: not QM output. Class B.

## (C) other licensed programs

- xtb: LGPL (github.com/grimme-lab/xtb, LGPL-3.0). Its output may ship.
- CREST: LGPL-3.0. Not present as fixtures.
- DELFIN does not drive Gaussian/Turbomole/Q-Chem/VASP/Molpro as fixture
  sources (no such outputs found under tests/). Rule family added anyway
  (task step 4) so future use is blocked.

## Consumers (tests that read class-A fixtures)

tests/test_opi_equivalence_q04.py, q05, q06, test_opi_orca_dlpno.py,
test_opi_orca_goat.py, test_opi_orca_sp_dft.py, test_orca_out_index.py,
test_rc1_adversarial_fixture.py, test_rc1_butane_scan_fixture.py,
test_rc2_adversarial_fixtures.py, test_rc2_tm_fixtures.py,
test_sh1_recovery_classification.py, test_ms1_spectra_checks.py,
test_a_method_reviewer_checks_the_method.py, test_opi_orca_read.py.

Plan per phase: phase 2 red test for guard rule family (tests/test_zz_ first,
then tests/test_lc1_*), extend scripts/license_guard.py (note: it currently
scans only tracked .py, license_guard.py:165-167 — non-py scan is the core
gap, reviewer point 1 confirmed), LICENSE_GUARD.md policy, private-store
listing for operator (step 2), then synthetic replacement fixtures.
