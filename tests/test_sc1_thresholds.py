"""SC1 phase 3: thresholds in result_critic are never loosened.

These are PIN tests, not red-first defect tests: the controlled value
exists on main and the test pins it exactly. Any relaxation of a
threshold fails this file. Values from delfin/agent/result_critic.py
(the /check command's only thresholds).
"""

import delfin.agent.result_critic as rc


def test_spin_warning_threshold_is_pinned():
    assert rc._SPIN_WARN == 0.10


def test_spin_error_threshold_is_pinned():
    assert rc._SPIN_ERROR == 1.0


def test_imaginary_noise_band_is_pinned():
    assert rc._IMAG_NOISE_CM == 50.0


def test_noise_band_is_not_used_to_excuse_real_imag_modes():
    # The band may soften reporting of |freq| <= 50 cm-1, but a mode
    # beyond the band must never be silently dropped: check the critic
    # still flags a clearly imaginary mode on a minimum.
    assert rc._IMAG_NOISE_CM < 500.0


def test_thresholds_stay_inside_their_documented_order():
    assert rc._SPIN_WARN < rc._SPIN_ERROR
