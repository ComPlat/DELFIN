"""Adversarial pins on scan_for_wrong_arithmetic (phase 4, 0ba71610).

The phase-4 check verifies a stated equation "A +- B = C" against its own
operands -- the only place the wave-10 sign flip (s17) is detectable.
These pins are the independent-reviewer contract: wrong sums and sign
errors are flagged, correct equations are not, and equations that merely
LIVE inside code blocks or tables (a scanner could naively trip on them)
are never flagged. Every case here is correct behaviour on the current
code -- they go RED if a future change over-fires (flagging code) or
under-fires (missing a wrong sign).

Reference: delfin/agent/verify_guard.py --
  ArithmeticFlag (:603), _ARITHMETIC_EQ_RE (:632),
  _read_decimal (:638) reads under both separator conventions,
  scan_for_wrong_arithmetic (:695, public entry).
"""

from __future__ import annotations

from delfin.agent.verify_guard import scan_for_wrong_arithmetic


def _eqs(text: str) -> list[str]:
    return [f.equation for f in scan_for_wrong_arithmetic(text)]


def test_a_wrong_sum_is_flagged():
    # 1.5 + 2.25 is 3.75, not 4.75 -- wrong magnitude.
    eqs = _eqs("the total is 1.5 + 2.25 = 4.75 kcal/mol")
    assert any("4.75" in e for e in eqs), eqs
    assert any("=" in e for e in eqs)


def test_a_wrong_difference_is_flagged():
    # 0.9 - 0.3 = 0.6, not 0.2. (Operands must be UNAMBIGUOUS: a
    # three-digit fractional tail like "0.062" is read both as a decimal
    # and as a grouped thousand, and the guard stays silent on it.)
    eqs = _eqs("0.9 - 0.3 = 0.2")
    assert eqs, eqs


def test_a_sign_error_is_flagged():
    # The operands give -0.056322, the claim says +0.056322: sign flipped.
    # This is the wave-10 s17 shape -- operands are IN the claim, so the
    # scanner can and must catch it.
    eqs = _eqs("the gap is 0.031522 - 0.087844 = 0.056322 eV")
    assert eqs, eqs


def test_a_correct_equation_is_not_flagged():
    eqs = _eqs("emission is 3.0 - 2.1 = 0.9 eV")
    assert eqs == [], eqs


def test_a_correct_sum_with_unicode_minus_is_not_flagged():
    # The unicode minus is normalised; a correct equation stays silent.
    eqs = _eqs("total = 2.5 + 1.25 = 3.75")
    assert eqs == [], eqs


def test_a_parenthesised_negative_operand_is_not_flagged_when_right():
    # E(S1) - (-1705.219605) as a correct statement: no flag.
    eqs = _eqs("E = 100.0 - (-1705.219605) = 1805.219605")
    assert eqs == [], eqs


def test_an_equation_inside_a_fenced_code_block_is_not_flagged():
    # A scanner that strips claim regions must not trip on code that
    # happens to contain "a - b = c".
    text = (
        "the function is:\n"
        "```\n"
        "def f():\n"
        "    return 0.062 - 0.031 = 0.041  # illustrative, not a claim\n"
        "```\n"
    )
    eqs = _eqs(text)
    assert eqs == [], eqs


def test_an_equation_inside_a_markdown_table_cell_is_not_flagged():
    # A table row "a | b | c" is not "a = b"; only real equations are.
    text = (
        "| metric | value |\n"
        "|--------|-------|\n"
        "| gap    | 3.90 - 2.31 = 1.59 |\n"
    )
    eqs = _eqs(text)
    assert eqs == [], eqs
