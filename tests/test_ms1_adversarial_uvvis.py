"""MS1 adversarial phase-3a tests (reviewer nacht-s4).

F2 — the triplet-in-singlet check (result_critic.py:401-410) fires on
ANY triplet state label, without reading the reference multiplicity.
A triplet-reference calculation (Multiplicity ... 3 in the echo, ORCA
manual 6.1.1 %tddft: TDA/TDDFT on a triplet reference lists triplet
excitations legitimately) would be flagged "error" — a false positive
on physically correct output.

F3 — the per-state <S^2> line regex (result_critic.py:383-385) accepts
only the builder's invented 6-column shape. Real recorded ORCA TDDFT
tables (tests/test_the_spectrum_says_which_line_is_bright.py) have the
extra D**2 X Y Z columns after fosc; with an S^2 column printed after
them the regex matches nothing and spin contamination is silently
never flagged on a real-shape table.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import result_critic as rc


def _codes(crits):
    return {c.code for c in crits}


def _write(tmp_path, name, content) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


_TABLE_HDR = (
    "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
    "-------------------------------------------------------------------\n"
    "State      Energy   Energy   Wavelength   fosc         S^2\n"
    "           (eV)    (cm-1)      (nm)\n"
    "-------------------------------------------------------------------\n"
)


def _triplet_ref_out() -> str:
    """Triplet-reference TDDFT output: multiplicity 3 in the echo, and
    triplet excitation rows (0-3A -> 1-3A) — physical, not a mixture.
    """
    return (
        "ORCA TERMINATED NORMALLY\n"
        "FINAL SINGLE POINT ENERGY     -129.654321\n"
        "Multiplicity                     ...    3\n"
        "Number of atoms                             ...      4\n"
        + _TABLE_HDR
        + "0-3A  ->  1-3A    2.200000     17744.2     563.6     0.010000\n"
        + "0-3A  ->  2-3A    2.650000     21373.7     467.9     0.412300\n"
        + "\n\n"
    )


def test_f2_triplet_reference_spectrum_is_not_flagged(tmp_path):
    """A triplet-reference calculation legitimately lists triplet states.
    Must not be flagged as 'triplet in singlet'."""
    out = _write(tmp_path, "tripref.out", _triplet_ref_out())
    crits = rc.critique_output(out)
    assert "uv-triplet-in-singlet" not in _codes(crits)


def _real_shape_out() -> str:
    """Real recorded ORCA shape: fosc column followed by D**2 X Y Z, then
    an S^2 column (as ORCA prints with higher TDDFT print levels).
    """
    return (
        "ORCA TERMINATED NORMALLY\n"
        "FINAL SINGLE POINT ENERGY     -123.456789\n"
        "Number of atoms                             ...      4\n"
        "Number of roots         ...   3\n"
        + _TABLE_HDR
        + "0-1A  ->  1-1A    4.100000     33068.0     302.06"
        + "     0.0335    0.500  0.1  0.2  0.0     1.0000\n"
        + "0-1A  ->  2-1A    4.900000     39522.0     253.04"
        + "     0.5000    0.400  0.1  0.2  0.0     0.0000\n"
        + "0-1A  ->  3-1A    5.400000     43558.0     229.54"
        + "     0.0010    0.300  0.1  0.2  0.0     0.0000\n"
        + "\n\n"
    )


def test_f3_spin_contamination_detected_in_real_shape_table(tmp_path):
    """S^2 deviation 1.0 in a table whose rows carry the real D**2 X Y Z
    columns before S^2 — must still be flagged."""
    out = _write(tmp_path, "realshape.out", _real_shape_out())
    crits = rc.critique_output(out)
    assert "uv-spin-contamination" in _codes(crits)
