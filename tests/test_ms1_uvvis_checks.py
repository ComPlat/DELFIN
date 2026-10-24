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
