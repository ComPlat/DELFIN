"""Red controls for reviewer s31's phase-2 findings F1, F3, F4, F5, F6.

Run against the unchanged phase-2 commit 3e35a5d2: all must be RED
there and green after the phase-2b fix.
"""
from __future__ import annotations

import pytest

import delfin.agent.data_checks as dchecks


def test_f1_underscore_thousands_is_undecided_not_a_number():
    """float() accepts '1_234' silently — must mark undecided instead."""
    result = dchecks.unit_consistency(["1_234", "1.5 kcal/mol"])
    assert result["consistent"] is False, result
    assert result["undecided"] is True, result


def test_f3_missing_key_column_is_a_refusal():
    """A typo'd key reports zero duplicates today; must refuse."""
    with pytest.raises(dchecks.DataCheckError):
        dchecks.find_duplicates(
            [{"a": 1, "b": 2}, {"a": 1, "b": 2}], key=["nonexistent"])


def test_f4_outlier_rejects_nan_and_inf():
    for bad in (float("nan"), float("inf")):
        with pytest.raises(dchecks.DataCheckError):
            dchecks.flag_outliers([1.0, 2.0, 3.0, bad])


def test_f5_outlier_rejects_non_positive_k():
    for bad in (-1.5, 0.0):
        with pytest.raises(dchecks.DataCheckError):
            dchecks.flag_outliers([1.0, 2.0, 3.0, 4.0], k=bad)


def test_f6_summarize_rejects_nan_and_bool():
    with pytest.raises(dchecks.DataCheckError):
        dchecks.summarize([1.0, float("nan")])
    with pytest.raises(dchecks.DataCheckError):
        dchecks.summarize([1.0, True])


def test_f2_bare_plus_unit_column_is_mixed_shape_not_consistent():
    """QS ruling: a bare number does not inherit the column unit."""
    result = dchecks.unit_consistency(["13.1 kcal/mol", "8.5"])
    assert result["consistent"] is False, result
    assert result["mixed_shape"] is True, result
    assert result["n_bare"] == 1, result
    assert result["undecided"] is False, result
    # a pure bare column is not mixed_shape
    pure = dchecks.unit_consistency(["1.5", "2.0"])
    assert pure["consistent"] is True and pure["mixed_shape"] is False
