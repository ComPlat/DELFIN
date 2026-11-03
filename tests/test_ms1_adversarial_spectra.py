"""MS1 adversarial review tests (reviewer nacht-s4).

Each test attacks an implemented MS1 check from a shape the builder's
own tests do not cover. Threshold sources are cited where a test pins
a number.

F1 — the linear-molecule marker ("linear molecule" in the output text)
cannot match real ORCA output: the manual (6.1.1 section 4.7
"Thermochemistry") shows ORCA reports linearity via the thermochemistry
block ("Point Group: ..., Symmetry Number: ..." and "Rotational
constants in cm-1: ..."), never a literal "linear molecule" phrase. A
linear molecule (e.g. CO2, N2, HCN) with the CORRECT 3N-5 modes must
not be flagged as a mode-count error.
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


def _linear_shape(n_real_modes: int) -> str:
    """Real-shape linear-molecule freq output (3 atoms, e.g. HCN).

    Follows the manual 6.1.1 section 4.7 output example: linearity shows
    through a zero first rotational constant (a linear rotor has one
    zero constant) — NOT a "linear molecule" phrase.
    """
    term = "ORCA TERMINATED NORMALLY\n"
    energy = "FINAL SINGLE POINT ENERGY     -92.912345678\n"
    atoms = "Number of atoms                             ...      3\n"
    freqs = _FREQ_HDR = "VIBRATIONAL FREQUENCIES\n-----------------------\n"
    freqs += "   0:         0.00 cm**-1\n"
    freqs += "   1:         0.00 cm**-1\n"
    for i in range(n_real_modes):
        freqs += f"   {2 + i}:    {700.0 + 300.0 * i:9.2f} cm**-1\n"
    thermo = (
        "Point Group: C*v, Symmetry Number: 1\n"
        "Rotational constants in cm-1: 0.000000 1.512345 1.512345\n"
    )
    return term + energy + atoms + freqs + "\n" + thermo


def test_f1_linear_molecule_correct_3n_minus_5_is_not_flagged(tmp_path):
    """HCN-like, 3 atoms, linear: 3N-5 = 4 real modes — correct output.

    On the builder's implementation the marker never matches, the check
    expects 3N-6 = 3, sees 4 real modes and wrongly reports an error.
    """
    out = _write(tmp_path, "hcn.out", _linear_shape(4))
    crits = rc.critique_output(out)
    assert "freq-count-mismatch" not in _codes(crits)
