"""SC1 red tests: result_claims.check — every number in an answer must
have a source.

An agent's ANSWER is the product here, not the calculation: a number
stated with a unit must come from an output that is part of the
evidence, an energy must name method/basis/solvent, a number whose
source the result critic flags (imaginary mode, unconverged SCF, spin
contamination) must carry that warning, and only a clearly marked
estimate is exempt from provenance. Structure, not word lists.
"""

from __future__ import annotations

import textwrap

from delfin.agent import result_claims as rcl


def _calc_evidence(text: str, flags=None) -> dict:
    return {"kind": "calc", "ref": "calc/example",
            "text": text, "flags": flags or []}


GOOD_OUT = textwrap.dedent("""\
    ORCA TERMINATED NORMALLY
    Total Charge              :     0
    Electronic Energy        -1234.567890 Eh
    FINAL SINGLE POINT ENERGY -1234.567890
    HURRAY! THE OPTIMIZATION HAS CONVERGED.
    No imaginary modes found.
    ...
""")


def test_made_up_energy_is_flagged():
    ans = "The reaction barrier is -15.3 kcal/mol."
    report = rcl.check(ans, evidence=[])
    assert rcl.worst_level(report) == "error"
    codes = {f.code for f in report}
    assert "claim-no-source" in codes


def test_number_from_evidence_passes():
    ans = ("The electronic energy is -1234.567890 Eh "
           "(B3LYP/def2-SVP, solvent water).")
    ev = [_calc_evidence(
        "Electronic Energy        -1234.567890 Eh\n"
        "HURRAY! THE OPTIMIZATION HAS CONVERGED.\n"
        "ORCA TERMINATED NORMALLY\n")]
    report = rcl.check(ans, evidence=ev)
    assert rcl.worst_level(report) == "ok", report


def test_energy_without_method_basis_solvent_is_flagged():
    ans = "The electronic energy is -1234.567890 Eh."
    ev = [_calc_evidence("Electronic Energy -1234.567890 Eh\n"
                         "ORCA TERMINATED NORMALLY\n")]
    report = rcl.check(ans, evidence=ev)
    codes = {f.code for f in report}
    assert "claim-no-method" in codes


def test_estimate_is_exempt_from_provenance():
    ans = "A barrier of about 15 kcal/mol is expected (rough estimate)."
    report = rcl.check(ans, evidence=[])
    codes = {f.code for f in report}
    assert "claim-no-source" not in codes


def test_flagged_source_carries_the_warning():
    ans = "The electronic energy is -1234.567890 Eh."
    ev = [_calc_evidence(
        "Electronic Energy -1234.567890 Eh\n"
        "ORCA TERMINATED NORMALLY\n",
        flags=["min-has-imag"])]
    report = rcl.check(ans, evidence=ev)
    codes = {f.code for f in report}
    assert "claim-flagged-source" in codes


def test_flagged_source_passes_when_answer_warrants_it():
    ans = ("The electronic energy is -1234.567890 Eh — NOTE: the "
           "frequency check found an imaginary mode; treat this "
           "structure as a saddle point, not a minimum.")
    ev = [_calc_evidence(
        "Electronic Energy -1234.567890 Eh\n"
        "ORCA TERMINATED NORMALLY\n",
        flags=["min-has-imag"])]
    report = rcl.check(ans, evidence=ev)
    codes = {f.code for f in report}
    assert "claim-flagged-source" not in codes


def test_never_raises_and_empty_answer_is_ok():
    assert rcl.check("", evidence=[]) == []
    assert rcl.check("Plain prose with no numbers at all.",
                     evidence=[]) == []
