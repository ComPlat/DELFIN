"""Red controls for delfin/agent/data_checks.py (package D1).

Each test states the failure it guards. All cases must be red on the
previous commit (module absent) and green after phase 2.
"""
import math

import pytest

dchecks = pytest.importorskip("delfin.agent.data_checks")


def test_unit_consistency_flags_mixed_units():
    """A column mixing 'kcal/mol' and 'kJ/mol' must be flagged, not summed."""
    values = ["12.5 kcal/mol", "52.3 kJ/mol", "13.1 kcal/mol"]
    result = dchecks.unit_consistency(values)
    assert result["consistent"] is False
    assert set(result["units"]) == {"kcal/mol", "kJ/mol"}
    # A same-unit column passes.
    ok = dchecks.unit_consistency(["1.5 kcal/mol", "2.0 kcal/mol"])
    assert ok["consistent"] is True


def test_unit_consistency_undecided_decimal_is_refused_not_guessed():
    """'8.986' cannot be 8.986 or 8986 — the helper must refuse."""
    result = dchecks.unit_consistency(["8.986", "12.34", "9.876"])
    assert result["consistent"] is False
    assert result["undecided"] is True


def test_find_duplicates_reports_repeated_blocks():
    """A duplicated block of rows must be reported with row indices."""
    rows = [
        {"id": "A", "v": 1},
        {"id": "B", "v": 2},
        {"id": "A", "v": 1},  # duplicate
        {"id": "C", "v": 3},
    ]
    result = dchecks.find_duplicates(rows, key=["id", "v"])
    assert result["duplicate_groups"] == [
        {"key": {"id": "A", "v": 1}, "rows": [0, 2]}
    ]
    clean = dchecks.find_duplicates(rows[:2] + [{"id": "C", "v": 3}],
                                    key=["id", "v"])
    assert clean["duplicate_groups"] == []


def test_flag_outliers_by_stated_rule():
    """An outlier must be flagged by the rule given, never by judgement."""
    values = [10.0] * 20 + [1000.0]
    result = dchecks.flag_outliers(values, rule="iqr", k=1.5)
    assert result["outlier_indices"] == [20]
    # The rule is carried back so the report can state it.
    assert result["rule"] == "iqr"
    # Inliers survive a mild rule.
    ok = dchecks.flag_outliers(list(range(10, 30)), rule="iqr", k=1.5)
    assert ok["outlier_indices"] == []


def test_flag_outliers_unknown_rule_refused():
    """An unknown outlier rule is a refusal, not a silent default."""
    with pytest.raises(ValueError):
        dchecks.flag_outliers([1.0, 2.0], rule="vibes")


def test_summary_carries_n_and_spread():
    """A summary must carry n and a spread measure — not just a mean."""
    result = dchecks.summarize([1.0, 2.0, 3.0, 4.0])
    assert result["n"] == 4
    assert result["mean"] == pytest.approx(2.5)
    assert result["std"] == pytest.approx(math.sqrt(5.0 / 3.0))
    assert result["median"] == 2.5
    assert result["min"] == 1.0 and result["max"] == 4.0


def test_summary_on_empty_input_refused():
    """A summary of nothing is a refusal — n=0 must never be silently a 0."""
    with pytest.raises(ValueError):
        dchecks.summarize([])
