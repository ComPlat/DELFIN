"""Adversarial tests for OCCUPIER package OC1 (reviewer nacht-s4).

These test the REVIEWER's view of the OCCUPIER code paths: the Gibbs/FSPE
energy extractors and the sequence-reference parsers. They are written
adversarially: they encode what a correct extractor should do, not what the
code happens to do.
"""
import pytest

from delfin.occupier import (
    _parse_dependency_indices,
    _resolve_primary_source,
)
from delfin.occupier_flat_extraction import _parse_energy


class TestGibbsExtractionLastLine:
    """F1 (s3 finding, independently verified): _parse_energy scans FORWARD
    for the Gibbs line, so with several Gibbs blocks in one file it returns
    the FIRST, not the last. find_FSPE (occupier.py:986) scans bottom-up and
    returns the LAST occurrence. An ORCA output that carries more than one
    thermochemistry block (restart appended, or a completed FREQ run whose
    earlier geometry steps printed partial blocks) must yield the last one.
    """

    def test_gibbs_returns_last_occurrence(self, tmp_path):
        out = tmp_path / "output.out"
        out.write_text(
            "ORCA TERMINATED NORMALLY\n"
            "Some header\n"
            "Final Gibbs free energy         ...  -100.000000000 eh\n"
            "FINAL SINGLE POINT ENERGY      -99.000000000 eh\n"
            "Second thermochemistry block after a restart:\n"
            "Final Gibbs free energy         ...  -101.500000000 eh\n",
            encoding="utf-8",
        )
        got = _parse_energy(out, use_gibbs=True)
        # The last Gibbs line is -101.5; a forward scan wrongly returns -100.0.
        assert got == pytest.approx(-101.5), (
            f"_parse_energy returned {got}; the last Gibbs line is -101.5 — "
            "the extractor scans forward instead of bottom-up"
        )

    def test_gibbs_bad_token_does_not_abort(self, tmp_path):
        """One unparseable Gibbs token must not end the search: the last block
        still holds a valid number and a wrong 'None' silently drops the FoB
        from the winner comparison.
        """
        out = tmp_path / "output2.out"
        out.write_text(
            "ORCA TERMINATED NORMALLY\n"
            "Final Gibbs free energy         ...  not-a-number eh\n"
            "Final Gibbs free energy         ...   -50.250000000 eh\n",
            encoding="utf-8",
        )
        got = _parse_energy(out, use_gibbs=True)
        assert got == pytest.approx(-50.25), (
            f"_parse_energy returned {got}; a malformed token in an earlier "
            "block aborted the scan instead of skipping it"
        )
