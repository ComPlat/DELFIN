"""Red control for phase 1: SLURM time-limit parsing must not invent values.

A time string with more than three colon-separated fields is not a SLURM
time and must be refused, not silently read as minutes. Pinned by the
reviewer's adversarial test first (nacht-s14, F1), green only after the
fix.
"""
from __future__ import annotations

import pytest

from delfin.slurm_submit import time_limit_seconds


def test_four_field_time_is_refused_not_read_as_minutes():
    with pytest.raises(ValueError):
        time_limit_seconds("1:2:3:4")


def test_known_forms_still_parse():
    assert time_limit_seconds("48:00:00") == 48 * 3600
    assert time_limit_seconds("2-00:00:00") == 2 * 86400
    assert time_limit_seconds("90") == 90 * 60
    # SLURM's 2-field form is MM:SS (cluster_utils._parse_walltime: "without
    # a day prefix it is MM:SS"); D-HH:MM only after a day count.
    assert time_limit_seconds("01:30") == 90
    assert time_limit_seconds("1-02") == 26 * 3600
