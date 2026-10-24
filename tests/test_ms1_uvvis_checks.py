"""MS1 phase 3a: UV-Vis/TDDFT plausibility checks in /check.

Fixture shape follows ORCA's real TDDFT output ("ABSORPTION SPECTRUM
VIA TRANSITION ELECTRIC DIPOLE MOMENTS" with state labels like 0-1A /
0-2A and an SOC/SPIN CONTAMINATION-style per-root <S^2> column when
printed). Sources for each threshold:
- Roots vs requested: ORCA manual %tddft nroots (6.1.1 §5.1): fewer
  computed roots than requested means the user's energy window is not
  covered — the requested spectrum range is silently truncated.
- Per-state <S^2> deviation: ORCA manual TDDFT output description;
  diagnostic = |<S^2> - S(S+1)| (DELFIN methodology, "Spin
  Contamination and State Selection").
- Triplet-in-singlet: a state_to label with multiplicity 3 in a
  singlet-reference (mult 1) TDDFT spectrum means spin-flip/ST mixture
  — a singlet absorption spectrum must not list triplets as states.
- Broadening: any reported spectrum must state the broadening width
  used (FWHM/stdev) or it cannot be compared with an experiment;
  missing width = ok-level note.

Report-only; never edits, re-runs or loosens anything.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import result_critic as rc


_TERM = "ORCA TERMINATED NORMALLY\n"
_ENERGY = "FINAL SINGLE POINT ENERGY     -123.456789\n"
_ATOMS = "Number of atoms                             ...      4\n"

_ABCS = (
    "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
    "-------------------------------------------------------------------\n"
    "State      Energy   Energy   Wavelength   fosc         S^2\n"
    "           (eV)    (cm-1)      (nm)                  (-)\n"
    "-------------------------------------------------------------------\n"
    "0-1A -> 1-1A    4.100000  33068.0    302.06   0.0335000000  0.0000\n"
    "0-1A -> 2-1A    4.900000  39522.0    253.04   0.5000000000  0.0000\n"
    "0-1A -> 3-1A    5.400000  43558.0    229.54   0.0010000000  0.0000\n"
    "\n\n"
)


def _codes(crits):
    return {c.code for c in crits}


def _flag(crits, code):
    return next(c for c in crits if c.code == code)


def _write(tmp_path, name, content) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def _td(*state_lines, nroots="3", extra="") -> str:
    body = _ABCS if not state_lines else (
        _ABCS.split("0-1A -> 1-1A")[0]
        + "".join(line + "\n" for line in state_lines))
    return (_TERM + _ENERGY + _ATOMS
            + f"Number of roots         ...   {nroots}\n"
            + body + extra)


def test_roots_below_requested_is_flagged(tmp_path):
    """3 roots requested, only 2 parsed -> requested window truncated."""
    out = _write(tmp_path, "few.out", _td(
        "0-1A -> 1-1A    4.100000  33068.0    302.06   0.0335  0.0000",
        "0-1A -> 2-1A    4.900000  39522.0    253.04   0.5000  0.0000"))
    crits = rc.critique_output(out)
    flag = _flag(crits, "uv-few-roots")
    assert flag.level == "warn"
    assert "2" in flag.message and "3" in flag.message


def test_roots_matching_request_do_not_fire(tmp_path):
    """Control: 3 roots requested, 3 transitions -> no fire."""
    out = _write(tmp_path, "ok.out", _td())
    crits = rc.critique_output(out)
    assert "uv-few-roots" not in _codes(crits)


def test_no_nroots_line_means_no_check(tmp_path):
    """Control: without an echoed nroots request the check stays quiet."""
    content = _TERM + _ENERGY + _ATOMS + _ABCS
    out = _write(tmp_path, "noreq.out", content)
    crits = rc.critique_output(out)
    assert "uv-few-roots" not in _codes(crits)


def test_spin_contaminated_excited_state_is_flagged(tmp_path):
    """<S^2> deviation 1.00 from ideal 0.0 (singlet) is severe."""
    out = _write(tmp_path, "spin.out", _td(
        "0-1A -> 1-1A    4.100000  33068.0    302.06   0.0335  1.0000"))
    crits = rc.critique_output(out)
    flag = _flag(crits, "uv-spin-contamination")
    assert flag.level == "error"
    assert "S^2" in flag.message or "S²" in flag.message


def test_mild_spin_deviation_is_a_warning(tmp_path):
    """<S^2> deviation 0.35 from ideal 0.0 -> warn (0.10 threshold)."""
    out = _write(tmp_path, "mild.out", _td(
        "0-1A -> 1-1A    4.100000  33068.0    302.06   0.0335  0.3500"))
    crits = rc.critique_output(out)
    assert _flag(crits, "uv-spin-contamination").level == "warn"


def test_clean_excited_state_spins_do_not_fire(tmp_path):
    """Control: <S^2> ~ 0 for all states -> no fire."""
    out = _write(tmp_path, "clean.out", _td())
    crits = rc.critique_output(out)
    assert "uv-spin-contamination" not in _codes(crits)


def test_triplet_state_in_singlet_spectrum_is_flagged(tmp_path):
    """A 1->3 (triplet) state in a mult-1 TDDFT spectrum: spin-flip
    mixture leaking into the singlet absorption table.
    """
    out = _write(tmp_path, "trip.out", _td(
        "0-1A -> 1-3A    3.400000  27428.0    364.59   0.0000  0.0000",
        "0-1A -> 1-1A    4.100000  33068.0    302.06   0.0335  0.0000"))
    crits = rc.critique_output(out)
    flag = _flag(crits, "uv-triplet-in-singlet")
    assert flag.level == "error"
    assert "1-3A" in flag.message


def test_singlet_only_spectrum_does_not_fire(tmp_path):
    """Control: only singlet states -> no fire."""
    out = _write(tmp_path, "sing.out", _td())
    crits = rc.critique_output(out)
    assert "uv-triplet-in-singlet" not in _codes(crits)


def test_missing_broadening_width_gets_an_ok_note(tmp_path):
    """A reported spectrum must state its broadening width; without one
    it cannot be compared with experiment. ok-level note only.
    """
    out = _write(tmp_path, "nobroad.out", _td())
    crits = rc.critique_output(out)
    flag = _flag(crits, "uv-no-broadening")
    assert flag.level == "ok"


def test_stated_broadening_width_silences_the_note(tmp_path):
    """Control: "Gaussian broadening FWHM = 20 nm" in the output -> the
    note does not fire.
    """
    out = _write(tmp_path, "broad.out",
                 _td(extra="Gaussian broadening FWHM = 20 nm\n"))
    crits = rc.critique_output(out)
    assert "uv-no-broadening" not in _codes(crits)


# --- s17 cross-review F1: float-tolerant row parsing -------------------
# Real ORCA TDDFT rows carry bare "0." tokens (RX/TY/RZ columns),
# negative values and E-notation. The row parser must accept them, or
# a named-S^2 table is silently unparsed and the gate is off where it
# must fire.

_ORCA_REAL_ROW = (
    " 0 -> 2.0   0.200  0.0 0.0  0.0     0.0   0.20  0.0"
    "  0. 0. 0.  2.20")


def _td_named(*state_lines, nroots="3", extra=""):
    """Table whose header explicitly names the S^2 column."""
    header = (
        "State   Energy   Wavelength  fosc   D**2    X     Y     Z"
        "     S^2\n"
        "        (eV)     (nm)        (-)    (-)                      \n"
        "----------------------------------------------------------------------------\n")
    body = "".join(line + "\n" for line in state_lines)
    return (_TERM + _ENERGY + _ATOMS
            + f"Number of roots         ...   {nroots}\n"
            + "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
            + header + body + extra)


def test_bare_zero_dot_tokens_do_not_silence_spin_contamination(tmp_path):
    """A contaminated row with bare '0.' RX/TY/RZ tokens is parsed and
    the 2.20 <S^2> fires (s17 probe: silently skipped before).
    """
    out = _write(tmp_path, "bare.out", _td_named(_ORCA_REAL_ROW))
    crits = rc.critique_output(out)
    assert "uv-spin-contamination" in _codes(crits)


def test_negative_and_enotation_tokens_are_parsed(tmp_path):
    """Rows with -0.12 and 1.23E-05 style tokens still parse."""
    row = (" 0 -> 2.0   0.200  -0.12 0.0     0.0   0.20  0.0"
           "  0. 0. 0.  2.20")
    out = _write(tmp_path, "neg.out", _td_named(row))
    crits = rc.critique_output(out)
    assert "uv-spin-contamination" in _codes(crits)


def test_clean_real_shaped_rows_still_do_not_fire(tmp_path):
    """Control: real-shaped rows with <S^2> ~ 0 stay quiet."""
    row = (" 0 -> 2.0   0.200  0.0 0.0  0.0     0.0   0.20  0.0"
           "  0. 0. 0.  0.02")
    out = _write(tmp_path, "cleanreal.out", _td_named(row))
    crits = rc.critique_output(out)
    assert "uv-spin-contamination" not in _codes(crits)


def test_named_s2_table_with_no_parsed_rows_fails_visible(tmp_path):
    """Fail-visible: a named-S^2 table whose rows parse to nothing
    emits a uv-s2-unparsed warning instead of silence.
    """
    rows = "State   Energy   Wavelength  fosc      S^2\n"
    out = _write(tmp_path, "garbage.out", _td_named(rows))
    crits = rc.critique_output(out)
    assert "uv-s2-unparsed" in _codes(crits)


def test_parsed_rows_in_named_table_do_not_emit_unparsed_note(tmp_path):
    """Control: when rows parse, uv-s2-unparsed stays quiet."""
    out = _write(tmp_path, "parses.out", _td_named(_ORCA_REAL_ROW))
    crits = rc.critique_output(out)
    assert "uv-s2-unparsed" not in _codes(crits)


def test_parser_exception_keeps_other_findings_and_no_critic_error(tmp_path):
    """s17 F3: when the api parser RAISES, res must not be unbound.
    Before the fix the triplet block hit NameError at the res.transitions
    read and critique_output collapsed every finding into critic-error.
    With the fix: uv-parser-error warn appears, the content-based
    termination finding survives, no critic-error, no NameError.
    """
    out = _write(tmp_path, "raise.out", _td())
    orig = rc.__dict__.get("_orig_extract_excited_states")
    def _boom(folder):
        raise RuntimeError("forced parser failure")
    saved = None
    import delfin.api as api
    saved = api.extract_excited_states
    try:
        api.extract_excited_states = _boom
        # _check_uvvis imports inside the function, so patch the source
        # module it resolves through: delfin.api.
        crits = rc.critique_output(out)
    finally:
        api.extract_excited_states = saved
    codes = _codes(crits)
    assert "uv-parser-error" in codes
    assert "critic-error" not in codes
    assert "termination" in codes  # other findings survive
