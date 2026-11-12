"""Adversarial tests for SC1: result_critic thresholds must never loosen
(package target 3) and flagged outputs must leak into claim checks
(target 1, item 4). Written by reviewer nacht-s4 before the fix."""

from __future__ import annotations

import pytest

from delfin.agent import result_critic as rc


def test_spin_thresholds_pinned():
    assert rc._SPIN_WARN == 0.10
    assert rc._SPIN_ERROR == 1.0


def test_imag_noise_threshold_pinned():
    assert rc._IMAG_NOISE_CM == 50.0


def test_worst_level_order_pinned():
    ok = [rc.Critique("ok", "a", "x")]
    warn = [rc.Critique("warn", "b", "x")]
    err = [rc.Critique("error", "c", "x")]
    assert rc.worst_level(ok) == "ok"
    assert rc.worst_level(warn) == "warn"
    assert rc.worst_level(err) == "error"
    assert rc.worst_level(ok + warn + err) == "error"


def test_ts_without_imag_is_a_warn_not_ok():
    # A TS search with no imaginary mode must never be reported as "ok".
    class Fake:
        code = "ts-no-imag"
        level = "warn"
    # source-level pin: the Critique produced for ts-no-imag is a warn
    import inspect
    src = inspect.getsource(rc._check_frequencies)
    assert '"warn", "ts-no-imag"' in src, "ts-no-imag must stay a warn"
    assert Fake is not None  # keep the object referenced


def test_nonnegative_energy_is_a_warn_not_ok():
    import inspect
    src = inspect.getsource(rc._check_energy)
    assert '"warn", "energy-nonneg"' in src
