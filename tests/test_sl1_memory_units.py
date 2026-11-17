"""Red control, SL1 phase 2: memory units parsed without a 512x overcount.

Reviewer F6 (nacht-s14, batch 2): _parse_memory_string puts 'B' in the
megabyte set, so '512B' reads as 512 MB. Bytes are their own unit.
"""
from __future__ import annotations

import pytest

from delfin.cluster_utils import _parse_memory_string


@pytest.mark.parametrize("text,expected_mb", [
    ("8000", 8000),        # bare number: MB by convention
    ("8000MB", 8000),
    ("8000M", 8000),
    ("8G", 8 * 1024),
    ("8GB", 8 * 1024),
    ("512KB", 0),          # < 1 MB truncates
    ("2T", 2 * 1024 * 1024),
])
def test_unit_table(text, expected_mb):
    assert _parse_memory_string(text) == expected_mb


def test_bytes_are_not_megabytes():
    # '512B' is half a kilobyte, never a 512 MB allocation.
    assert _parse_memory_string("512B") == 0
