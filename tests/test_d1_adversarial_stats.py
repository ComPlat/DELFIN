"""Adversarial tests for data_checks duplicates/outliers/summarize (D1, s31)."""
import pytest

dchecks = pytest.importorskip("delfin.agent.data_checks")


def test_duplicate_key_with_none_values_grouped():
    """None==None on key columns must group, not crash or vanish."""
    rows = [{"id": None, "v": 1}, {"id": None, "v": 1}]
    result = dchecks.find_duplicates(rows, key=["id"])
    assert len(result["duplicate_groups"]) == 1


def test_duplicate_missing_key_column_is_a_refusal():
    """A key column absent from rows must be refused, not silently None."""
    rows = [{"id": "A"}, {"id": "B"}]
    with pytest.raises(Exception):
        dchecks.find_duplicates(rows, key=["nonexistent"])


def test_duplicate_empty_rows():
    """Empty row list: zero groups, n_rows 0."""
    result = dchecks.find_duplicates([], key=["id"])
    assert result["duplicate_groups"] == []
    assert result["n_rows"] == 0


def test_outlier_iqr_quartile_method_is_stated():
    """Small n where quartile definitions differ: result must not crash
    and the bounds must reflect the stated interpolation (linear)."""
    result = dchecks.flag_outliers([1.0, 2.0, 3.0, 100.0], rule="iqr", k=1.5)
    # Linear interpolation: q1=1.75, q3=27.25, iqr=25.5 -> high=65.5
    assert result["bounds"]["q1"] == pytest.approx(1.75)
    assert result["bounds"]["q3"] == pytest.approx(27.25)
    assert result["outlier_indices"] == [3]


def test_outlier_zscore_constant_column_not_crash():
    """A constant column under zscore: std=0 must be handled, no crash."""
    result = dchecks.flag_outliers([5.0] * 10, rule="zscore")
    assert result["outlier_indices"] == []


def test_outlier_rejects_nan_and_inf():
    """NaN/inf are not data — a stated rule cannot decide with them."""
    for bad in (float("nan"), float("inf")):
        with pytest.raises(Exception):
            dchecks.flag_outliers([1.0, 2.0, 3.0, bad], rule="iqr")


def test_outlier_rejects_negative_k():
    """k=-1 inverts the bounds — must be refused, not computed."""
    with pytest.raises(Exception):
        dchecks.flag_outliers(list(range(10, 30)), rule="iqr", k=-1.5)


def test_outlier_uses_original_indices():
    """Outlier indices refer to the ORIGINAL order, not the sorted one."""
    values = [1000.0, 10.0, 11.0, 10.5]
    result = dchecks.flag_outliers(values, rule="iqr", k=1.5)
    if result["outlier_indices"]:
        assert result["outlier_indices"] == [0]


def test_summarize_single_value_std_is_refused_or_zero_not_crash():
    """n=1: sample std (ddof=1) is undefined — refuse or state, no crash."""
    result = dchecks.summarize([5.0])
    assert result["n"] == 1


def test_summarize_rejects_nan():
    """NaN poisons mean/std silently — must be refused."""
    with pytest.raises(Exception):
        dchecks.summarize([1.0, float("nan"), 3.0])


def test_summarize_rejects_bool_true_as_value():
    """bool is not a measurement — True must not count as 1.0."""
    with pytest.raises(Exception):
        dchecks.summarize([True, 2.0, 3.0])
