"""HX1 Phase 2 red test: noise-aware contrast between probe and
reference runs.

Contract under test: delfin/agent/harness_loop.contrast_supported
- a probe vs reference gain is supported only when the mean quality
  gain exceeds the noise band from bench_watch (max(8, 1.5*spread) —
  the same threshold the nightly watcher uses for regressions);
- a gain INSIDE the noise band is recorded as no-contrast, never as
  a gain (the brief's second red test; bench_watch: identical runs
  flip 6/48 tasks, stdev ~20);
- the returned verdict names the band that was applied, so the report
  can cite it (QS requirement 1: "a contrast that names the stdev
  band from compare_run, not just a point estimate");
- paired samples are required: reference and probe must have the
  same task ids in the same order (matched conditions), otherwise
  ValueError — unmatched samples cannot give a contrast.
"""
from __future__ import annotations

import pytest

from delfin.agent.harness_loop import contrast_supported

# Per-task paired quality samples.  Reference noise stdev ~20 → band
# = max(8, 1.5*20) = 30.
_TASKS = ["t1", "t2", "t3", "t4"]
_REFERENCE = [40.0, 45.0, 35.0, 50.0]     # mean 42.5
_NOISY_PROBE = [45.0, 50.0, 40.0, 55.0]   # +5.0 mean gain: INSIDE band
_BIG_PROBE = [70.0, 75.0, 65.0, 80.0]     # +30.0 mean gain: AT the band


def test_gain_inside_noise_band_is_not_a_gain() -> None:
    verdict = contrast_supported(_TASKS, _REFERENCE, _NOISY_PROBE,
                                 spread=20.0)
    assert verdict["supported"] is False
    assert verdict["band"] == pytest.approx(30.0)
    assert verdict["mean_gain"] == pytest.approx(5.0)


def test_gain_at_the_band_is_not_supported() -> None:
    # Strictly greater required: a gain equal to the band is not
    # evidence (the watcher itself flags only above the threshold).
    verdict = contrast_supported(_TASKS, _REFERENCE, _BIG_PROBE,
                                 spread=20.0)
    assert verdict["supported"] is False


def test_gain_clearly_above_the_band_is_supported() -> None:
    probe = [75.0, 80.0, 70.0, 85.0]  # +35.0 > 30.0
    verdict = contrast_supported(_TASKS, _REFERENCE, probe, spread=20.0)
    assert verdict["supported"] is True
    assert verdict["band"] == pytest.approx(30.0)
    assert verdict["mean_gain"] == pytest.approx(35.0)


def test_band_uses_bench_watch_floor_for_quiet_tasks() -> None:
    # spread=0 → band = QUALITY_DROP_FLOOR (8), not 0.
    verdict = contrast_supported(["t1"], [50.0], [60.0], spread=0.0)
    assert verdict["band"] == pytest.approx(8.0)
    assert verdict["supported"] is True


def test_unmatched_tasks_are_refused() -> None:
    with pytest.raises(ValueError):
        contrast_supported(["t1", "t2"], [50.0, 60.0], [55.0, 65.0, 70.0],
                           spread=10.0)


def test_no_reference_samples_are_refused() -> None:
    with pytest.raises(ValueError):
        contrast_supported([], [], [], spread=10.0)
