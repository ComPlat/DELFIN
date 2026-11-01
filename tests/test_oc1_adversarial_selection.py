"""Adversarial tests: OCCUPIER sequence-reference parsers and selection."""
import pytest

from delfin.occupier import (
    _parse_dependency_indices,
    _resolve_primary_source,
)
from delfin.reporting.occupier_selection import OccupierSelector


class TestResolvePrimarySource:
    def test_zero_means_geometry_start(self):
        # 0 means "start from the input geometry" — must be preserved, not
        # silently turned into a fallback.
        assert _resolve_primary_source(0, fallback=5) == 0

    def test_list_picks_first_nonnegative(self):
        assert _resolve_primary_source([0, 3], fallback=1) == 0

    def test_all_negative_list_falls_back(self):
        assert _resolve_primary_source([-1, -2], fallback=4) == 4

    def test_none_falls_back(self):
        assert _resolve_primary_source(None, fallback=7) == 7

    def test_string_with_spaces(self):
        assert _resolve_primary_source(" 3 ", fallback=0) == 3

    def test_garbage_string_falls_back(self):
        assert _resolve_primary_source("Fe", fallback=2) == 2


class TestParseDependencyIndices:
    def test_none_empty_zero(self):
        assert _parse_dependency_indices(None) == set()
        assert _parse_dependency_indices("") == set()
        assert _parse_dependency_indices(0) == set()

    def test_separated_string(self):
        assert _parse_dependency_indices("2;3|4") == {2, 3, 4}

    def test_zero_dependency_ignored(self):
        # 0 is the geometry start, not a FoB index; it must not enter deps.
        assert _parse_dependency_indices("0,2") == {2}

    def test_negative_ignored(self):
        assert _parse_dependency_indices("-1,2") == {2}


class TestOccupierSelectorEdgeCases:
    """Selection must not crash or misrank on degenerate inputs."""

    def _sel(self, sequence, dev_cache=None, **cfg):
        return OccupierSelector(
            {"occupier_selection": "tolerance", **cfg},
            sequence,
            dev_cache or {},
        )

    def test_all_energies_none(self):
        seq = [{"index": 1, "m": 3}, {"index": 2, "m": 5}]
        assert self._sel(seq).select_candidates([None, None]) == (None, None)

    def test_empty_sequence(self):
        assert self._sel([]).select_candidates([]) == (None, None)

    def test_picks_minimum_energy(self):
        seq = [{"index": 1, "m": 3}, {"index": 2, "m": 5}]
        idx, e = self._sel(seq).select_candidates([-10.0, -10.5])
        assert idx == 2 and e == pytest.approx(-10.5)

    def test_dev_max_filters_out_contaminated_candidate(self):
        seq = [{"index": 1, "m": 3}, {"index": 2, "m": 5, "BS": "2,1"}]
        dev = {1: 0.01, 2: 2.0}
        # Without dev_max the contaminated lower energy wins; with dev_max=1.0
        # it is excluded and FoB 1 is selected.
        idx, _ = self._sel(seq, dev_cache=dev, dev_max=1.0).select_candidates(
            [-10.5, -10.0]
        )
        assert idx == 1
