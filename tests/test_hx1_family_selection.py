"""HX1 Phase 1 red test: weakest-family selection with a ceiling skip.

Contract under test: delfin/agent/harness_loop.pick_weakest_family
- picks the family with the LOWEST score among those with headroom;
- skips families at or above the ceiling (classic families scored
  91-95 in the 31.10. E2 baseline must not be picked);
- raises LookupError when NO family has headroom;
- raises ValueError on an empty score mapping.
"""
from __future__ import annotations

import pytest

from delfin.agent.harness_loop import pick_weakest_family


def test_picks_lowest_scoring_family_with_headroom() -> None:
    scores = {"generalist": 55.0, "workflow": 62.0, "esd": 93.0}
    assert pick_weakest_family(scores) == "generalist"


def test_family_at_ceiling_is_skipped() -> None:
    scores = {"esd": 93.0, "occupier": 91.0, "workflow": 70.0}
    assert pick_weakest_family(scores, ceiling=90.0) == "workflow"


def test_family_exactly_at_ceiling_is_also_skipped() -> None:
    scores = {"esd": 90.0, "workflow": 88.0}
    assert pick_weakest_family(scores, ceiling=90.0) == "workflow"


def test_no_headroom_raises_lookup_error() -> None:
    with pytest.raises(LookupError):
        pick_weakest_family({"esd": 95.0, "im": 92.0}, ceiling=90.0)


def test_empty_scores_raises_value_error() -> None:
    with pytest.raises(ValueError):
        pick_weakest_family({})
