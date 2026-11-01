"""The ISC-vs-reference comparison is a product check, not a script.

The referee logic lives in delfin/common/esd_numerics.py next to the
reference value it compares against: isc_level_matches() gates the
comparison on the run's level (functional, dispersion, basis, solvent,
linewidth) and isc_reference_comparison() produces the verdict record
(PASS/FAIL/NOT GRADED with number, unit, source, reason) that any
consumer — a grader script, a report, a test — can rely on. Before the
logic landed there, only the uncommitted .gate grader knew these rules
(review finding F5; operator ruling 2026-10-31: product checks belong
in owned files).
"""

from __future__ import annotations

from delfin.common.esd_numerics import (
    ISC_REFERENCE,
    isc_level_matches,
    isc_reference_comparison,
)


class TestLevelGate:
    def test_matched_level_passes(self):
        control = {
            "functional": "PBE0",
            "disp_corr": "D4",
            "main_basisset": "def2-SVP",
            "solvent": "",
        }
        ok, reason = isc_level_matches(control, "Homogeneous linewidth is: 50 cm-1")
        assert ok is True
        assert "match" in reason.lower()

    def test_functional_mismatch_is_named(self):
        control = {
            "functional": "B3LYP",
            "disp_corr": "D4",
            "main_basisset": "def2-SVP",
            "solvent": "",
        }
        ok, reason = isc_level_matches(control, None)
        assert ok is False
        assert "functional=B3LYP" in reason

    def test_solvent_mismatch_is_named(self):
        control = {
            "functional": "PBE0",
            "disp_corr": "D4",
            "main_basisset": "def2-SVP",
            "solvent": "acetone",
        }
        ok, reason = isc_level_matches(control, None)
        assert ok is False
        assert "solvent" in reason

    def test_linewidth_mismatch_is_named(self):
        control = {
            "functional": "PBE0",
            "disp_corr": "D4",
            "main_basisset": "def2-SVP",
            "solvent": "",
        }
        ok, reason = isc_level_matches(control, "Homogeneous linewidth is: 120 cm-1")
        assert ok is False
        assert "LINEW=120" in reason

    def test_linewidth_absent_from_output_does_not_fail(self):
        control = {
            "functional": "PBE0",
            "disp_corr": "D4",
            "main_basisset": "def2-SVP",
            "solvent": "",
        }
        ok, reason = isc_level_matches(control, None)
        assert ok is True


class TestReferenceComparison:
    def test_matched_rate_in_band_passes(self):
        record = isc_reference_comparison(
            rate=2.5e7, problem=None, level_ok=True
        )
        assert record["verdict"] == "PASS"
        assert record["unit"] == "ratio"
        assert "2.46e+07" in record["value"] or "2.46" in record["value"]

    def test_rate_outside_band_fails(self):
        record = isc_reference_comparison(rate=5.0e8, problem=None, level_ok=True)
        assert record["verdict"] == "FAIL"

    def test_level_mismatch_is_not_graded(self):
        record = isc_reference_comparison(
            rate=2.5e7, problem=None, level_ok=False, level_reason="functional=B3LYP"
        )
        assert record["verdict"] == "NOT GRADED"
        assert "functional=B3LYP" in record["value"]

    def test_missing_rate_is_not_graded(self):
        record = isc_reference_comparison(rate=None, problem=None, level_ok=True)
        assert record["verdict"] == "NOT GRADED"

    def test_problem_rate_is_never_compared(self):
        record = isc_reference_comparison(
            rate=2.5e7,
            problem="WARNING: negative rates are unphysical",
            level_ok=True,
        )
        assert record["verdict"] == "NOT GRADED"
        assert "negative rates" in record["value"]

    def test_reference_is_the_in_tree_value(self):
        assert ISC_REFERENCE["rate"] == 2.46e7
        assert ISC_REFERENCE["source"] == (
            "delfin/common/esd_numerics.py docstring "
            "(ORCA 6.1.1 formaldehyde, TD-DFT, LINEW 50, MAXTIME auto)"
        )
