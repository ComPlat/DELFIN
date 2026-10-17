"""Adversarial tests for data_checks.unit_consistency (package D1, reviewer s31).

Every case here is an ALLOW-direction attack: a value the builder's own
tests did NOT feed, where a wrong answer looks plausible. All must pass
against 3e35a5d2; a red one is a finding.
"""
import math

import pytest

dchecks = pytest.importorskip("delfin.agent.data_checks")


def test_underscore_thousands_read_as_german_decimal_is_undecided_not_number():
    """1_234 must not silently become 1234 nor a unit string."""
    result = dchecks.unit_consistency(["1_234", "1.5 kcal/mol"])
    # It is not a readable number and not the same unit column: refuse.
    assert result["consistent"] is False


def test_whitespace_unit_still_detected():
    """'12.5  kcal/mol' (double space) must still carry the unit."""
    result = dchecks.unit_consistency(["12.5  kcal/mol", "13.1 kcal/mol"])
    assert result["consistent"] is True
    assert result["units"] == ["kcal/mol"]


def test_case_sensitive_units_are_different_readings():
    """'KJ/MOL' vs 'kJ/mol' collide as distinct units — do not fold case."""
    result = dchecks.unit_consistency(["12.5 kJ/mol", "13.1 KJ/MOL"])
    assert result["consistent"] is False
    assert set(result["units"]) == {"kJ/mol", "KJ/MOL"}


def test_bare_float_object_has_no_unit():
    """A real float carries no unit; a float-typed column stays consistent."""
    result = dchecks.unit_consistency([1.5, 2.5, 3.0])
    assert result["consistent"] is True
    assert result["units"] == []


def test_bool_is_not_silently_a_number():
    """bool must not pass the isinstance(number) branch as 1.0/0.0."""
    result = dchecks.unit_consistency([True, 1.5])
    assert result["consistent"] is False


def test_negative_thousands_is_undecided():
    """-1.234 has both readings (minus 1.234 / minus 1234) — undecided."""
    result = dchecks.unit_consistency(["-1.234", "2.5"])
    assert result["undecided"] is True
    assert result["consistent"] is False


def test_three_group_thousands_is_undecided():
    """1.234.567 (dot-thousands, two groups) is also a both-readings case."""
    result = dchecks.unit_consistency(["1.234.567", "2.5"])
    assert result["undecided"] is True or result["unreadable"] != []
    assert result["consistent"] is False


def test_percent_and_degree_units():
    """'%' and '°' are units; '45 %' and '45 %' must agree, not collide."""
    result = dchecks.unit_consistency(["45 %", "46 %", "47 %"])
    assert result["consistent"] is True
    assert result["units"] == ["%"]


def test_none_and_na_are_missing_not_unreadable():
    """None/''/'n/a' are missing values — never unreadable, never a unit."""
    result = dchecks.unit_consistency([None, "", "n/a", "1.5 kcal/mol"])
    assert result["unreadable"] == []
    assert result["consistent"] is True


def test_truly_unreadable_cell_is_named_by_index():
    """'1,5 mol L' (comma decimal, unknown unit) is unreadable at index 0."""
    result = dchecks.unit_consistency(["1,5 mol L", "1.5 kcal/mol"])
    assert result["unreadable"] == [0]
    assert result["consistent"] is False


def test_zero_and_exponent_notation():
    """1e-3 and 0 parse; scientific notation is not a unit string."""
    result = dchecks.unit_consistency(["1e-3 kcal/mol", "0 kcal/mol"])
    assert result["consistent"] is True


def test_unicode_minus_is_not_a_number():
    """U+2212 MINUS is not ASCII '-' — a hostile look-alike must not parse."""
    result = dchecks.unit_consistency(["−1.5 kcal/mol", "1.5 kcal/mol"])
    # Either unreadable (refusal) or consistent — but never silently −1.5==1.5.
    assert result["unreadable"] != [] or result["consistent"] is False
