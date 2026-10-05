"""An equation the answer states is checked against its own arithmetic.

Welle 11, package B, phase 4 (narrowed with a reason -- see the docstring
of ``scan_for_wrong_arithmetic``). Wave 10 (s17): the interpretation said
"too positive" while the difference the cited sources give is negative.

WHY THE BROAD VERSION IS IMPOSSIBLE, stated so it does not get retried:
checking the sign of every derived difference against the observation
pool cannot work. The pool is an unordered bag of floats with no pair
identity, and for ANY two pool numbers a and b both ordered differences
(a-b and b-a) exist, so every sign of a matched magnitude is compatible
with SOME ordered difference. A sign check there is a no-op by
construction. The wave-10 defect is only detectable where the claim
NAMES ITS OPERANDS: "0.087844 - 0.031522 = -0.056322" states A, B and C,
and then the arithmetic -- sign included -- is checkable exactly.

So the check that exists: an equation "A - B = C" or "A + B = C" written
in the answer must be true. It fires only when the text itself asserts
an identity that is false, which is what the wave-10 answer did (it
wrote the operands and then a result whose sign contradicted them).
"""

from __future__ import annotations

from delfin.agent.verify_guard import scan_for_wrong_arithmetic


def test_a_sign_flipped_subtraction_is_flagged():
    """The wave-10 shape: operands stated, result sign impossible."""
    flags = scan_for_wrong_arithmetic(
        "ΔEST = 0.031522 - 0.087844 = 0.056322 Hartree.")
    assert len(flags) == 1, flags
    assert flags[0].claimed == "0.056322"


def test_a_correct_subtraction_passes():
    flags = scan_for_wrong_arithmetic(
        "ΔEST = 0.031522 - 0.087844 = -0.056322 Hartree.")
    assert flags == []


def test_a_correct_addition_passes():
    flags = scan_for_wrong_arithmetic("Summe: 1.5 + 2.25 = 3.75")
    assert flags == []


def test_a_wrong_addition_is_flagged():
    flags = scan_for_wrong_arithmetic("Summe: 1.5 + 2.25 = 4.75")
    assert len(flags) == 1, flags


def test_sign_flipped_with_unicode_minus():
    flags = scan_for_wrong_arithmetic(
        "ΔE = 0.087844 − 0.031522 = −0.056322 Hartree.")
    assert len(flags) == 1, flags


def test_inside_a_fenced_block_is_not_checked():
    """Code blocks are not claims -- the same rule every scanner here
    follows."""
    text = "```\n0.031522 - 0.087844 = 0.056322\n```"
    assert scan_for_wrong_arithmetic(text) == []


def test_prose_without_an_equation_is_not_checked():
    assert scan_for_wrong_arithmetic(
        "Die Differenz beträgt -2.31 eV.") == []


def test_rounded_results_within_printed_precision_pass():
    """The result is held to the precision it prints, like every
    quantity claim here: 0.0219 stands for 0.0218707..."""
    flags = scan_for_wrong_arithmetic(
        "Gap: -1705.197735 - (-1705.219605) = 0.021871 Hartree.")
    assert flags == [], flags


def test_small_differences_near_zero_err_toward_silence():
    """A disagreement inside the absolute tolerance (5e-3, the module's
    own floor) is not a claim this check makes: a guard that has to
    guess should guess toward silence."""
    flags = scan_for_wrong_arithmetic("x: 0.001 - 0.002 = -0.003")
    assert flags == [], flags


def test_energy_scale_equation_is_checked():
    """The real magnitudes: two Hartree energies, wrong-sign result. The
    equation as answers actually write it -- operand minus a parenthesised
    negative operand."""
    flags = scan_for_wrong_arithmetic(
        "E(S1) - E(T1) = -1705.197735 - (-1705.219605) = -0.021871")
    assert len(flags) == 1, flags
    # And the correctly-signed form passes:
    ok = scan_for_wrong_arithmetic(
        "E(S1) - E(T1) = -1705.197735 - (-1705.219605) = 0.021871")
    assert ok == [], ok


def test_the_message_names_the_operands_and_both_results():
    flags = scan_for_wrong_arithmetic(
        "ΔEST = 0.031522 - 0.087844 = 0.056322 Hartree.")
    msg = flags[0].message()
    assert "0.031522" in msg and "0.087844" in msg
    assert "-0.056322" in msg and "0.056322" in msg
