# PROVENANCE — sh1 fixtures (package SH1, wave 24)

All four fixtures are trimmed copies of real ORCA 6.1.x outputs from the
wave-23 failures the operator copied into `.gate/evidence/` for this
package. Anonymisation applied to every file: the banner host name was
replaced with `fixture-host` (the placeholder `tests/
test_orca_fixtures_carry_no_host_identity.py` enforces). No home paths,
user names or other host identities in any file.

## sh1_input_error_nofreq — v2.out (complete copy, 189 lines)

- Run: relaxed-scan step input refused at parse time.
- What ORCA wrote: `INPUT ERROR` + `UNRECOGNIZED OR DUPLICATED
  KEYWORD(S) IN SIMPLE INPUT LINE` + the keyword `NOFREQ`, then the
  source trace `orca_main/main_input_keywordline.cpp, line 13872`.
- Why NOFREQ is not an ORCA keyword: the frequency control is `Freq`
  (with e.g. `NoFreqs` not existing either); the correct way to skip a
  frequency job is not to request one. Verified against the ORCA 6.1.1
  manual index (doc search, "simple keywords", Freq entry).
- Recovery class today: OrcaErrorType.INPUT_ERROR
  (delfin/orca_recovery.py:112-125, pattern "INPUT ERROR" +
  "UNRECOGNIZED OR DUPLICATED KEYWORD"). Correctly NOT healed: an
  input error goes to the user. This fixture is the CONTROL case —
  the classification must keep firing, no change may break it.

## sh1_gstep_constraints — v3.out (complete copy, 299 lines)

- Run: relaxed scan (dihedral) whose start geometry sits 180 degrees
  away from the requested starting constraint.
- What ORCA wrote (line 296-298): `GSTEP: could not impose initial
  constraints. Please check your input and try again!` from
  `orca_tools/Tool-GStep/qcgstep.cpp, line 521`.
- Recovery class today: UNKNOWN (no pattern contains GSTEP or
  "could not impose"). Gap F1: the run is retried generically instead
  of being told to the user with the real reason.

## sh1_scan_maxcycles_shark — v4.out trimmed (1050 lines)

- Run: relaxed scan along a butane dihedral with r2SCAN-3c TightOpt
  TightSCF (pal 2). The last scan step sits on the +-180 degree
  dihedral wrap; the constrained optimisation inside that step runs to
  its maximum cycles; the SCF energy collapses inside that step; the
  PROPERTIES stage (SHARK) then cannot read its input file and the
  whole run aborts.
- What ORCA wrote: line ~158241 of the original: `WARNING !!! The
  optimization did not converge but reached the maximum number of
  optimization cycles. Please check your results very carefully.`,
  then `Storing optimized geometry ... failed!` x3, then
  `ERROR (SHARK): Failed to read input file (...SHARKINP.tmp)` twice
  and `ORCA finished by error termination in PROPERTIES`.
- Trim: the first 700 lines (header, input echo, method banner, first
  scan steps) plus the last 350 lines (the failure block). The middle
  (26 intermediate scan steps' nearly identical iteration tables) is
  not in the committed file; nothing a parser reads was cut. The .inp
  and the scan data files are not committed (the failure is in the
  .out alone).
- Recovery class today: UNKNOWN. Gap F2: neither the out-of-cycles
  warning inside a scan step nor the SHARK/PROPERTIES abort is
  classified; the real cause (unconverged scan step) is not surfaced.

## sh1_cu_nh3_imagmode — cu_v2 (out 30185 lines, inp 522, hess 1363)

- System: [Cu(NH3)4]2+ tetraamminecopper(II), doublet, r2SCAN-3c +
  CPCM(water), TightOpt + NumFreq. Same run family as the committed
  `rc2_cu_nh3_4` fixture (whose PROVENANCE documents the whole v1/v2/v3
  story in detail).
- What ORCA wrote: a CONVERGED optimisation ("THE OPTIMIZATION HAS
  CONVERGED") whose frequencies show exactly ONE imaginary mode,
  `-21.86 cm**-1  ***imaginary mode***` at .out:29313, mode index 6
  (NH3 rocking character, documented in rc2_cu_nh3_4/PROVENANCE.md
  lines 61-73 with the ORCA manual section 1.7.14 reference).
- Why this fixture: IMAG's default noise floor is -50 cm-1
  (delfin/imag.py:85), a measured value (control_validator.py:1227-
  1234: 71 archived structures, 60 had modes in [-1,-50] cm-1 whose
  displaced single points sat within +-1 mEh of the saddle). -21.86
  cm-1 sits inside that noise band, so with the DEFAULT no key IMAG
  skips the mode by design. The red tests pin BOTH sides: with the
  default the mode is skipped and not warned about; with an explicit
  allow_imaginary_freq (e.g. 20) IMAG acts on the real frequency
  block. The committed .hess carries the real `$vibrational_frequencies`
  block and the real displacements, so the test exercises the true
  reader, not a typed dict.
- QS ruling 2026-10-31 (nacht-s5, recorded in .gate/HANDOFF.md): the
  noise floor stays; the fix is the route for an explicit value.
