"""Adversarial tests: OCCUPIER contamination trigger (module level)."""
import pytest

from delfin.occupier import _check_contamination_trigger


class TestContaminationTrigger:
    """_check_contamination_trigger decides whether an extra BS(m,1) run
    is added after the main sequence. Adversarial cases: disabled flag,
    BS winner, no results, threshold edge."""

    def test_disabled_returns_none(self):
        seq = [{"index": 1, "m": 3}]
        assert _check_contamination_trigger(
            seq, {1: -10.0},
            {"contamination_bs_trigger": "no"},
            freq_enabled=False, apm="2",
        ) is None

    def test_no_results_returns_none(self):
        seq = [{"index": 1, "m": 3}]
        assert _check_contamination_trigger(
            seq, {},
            {"contamination_bs_trigger": "yes"},
            freq_enabled=False, apm="2",
        ) is None

    def test_bs_winner_not_retriggered(self):
        seq = [{"index": 1, "m": 3}, {"index": 2, "m": 2, "BS": "2,1"}]
        assert _check_contamination_trigger(
            seq, {1: -10.0, 2: -11.0},
            {"contamination_bs_trigger": "yes"},
            freq_enabled=False, apm="2",
        ) is None

    def test_clean_winner_no_trigger(self, tmp_path, monkeypatch):
        # Winner contamination below threshold -> no extra job.
        seq = [{"index": 1, "m": 3}]
        monkeypatch.chdir(tmp_path)
        (tmp_path / "output.out").write_text(
            "ORCA TERMINATED NORMALLY\n", encoding="utf-8")
        entry = _check_contamination_trigger(
            seq, {1: -10.0},
            {"contamination_bs_trigger": "yes",
             "contamination_bs_threshold": "0.5"},
            freq_enabled=False, apm="2",
        )
        assert entry is None

    def test_threshold_accepts_fractional_string(self):
        # A fractional threshold string must parse; a garbage value falls
        # back to 0.5 without raising.
        seq = [{"index": 1, "m": 3}]
        # No output files exist -> trigger cannot read contamination and
        # must return None, not raise.
        for cfg in ({"contamination_bs_threshold": "0.75"},
                    {"contamination_bs_threshold": "garbage"}):
            assert _check_contamination_trigger(
                seq, {1: -10.0},
                {"contamination_bs_trigger": "yes", **cfg},
                freq_enabled=False, apm="2",
            ) is None
