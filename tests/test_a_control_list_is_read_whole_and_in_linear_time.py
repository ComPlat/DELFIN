"""A multi-line list in a CONTROL text is read in linear time.

MEASURED 2026-11-23 under Python 3.11, the interpreter DELFIN and CI run
on. ``config._collect_literal_list`` built its buffer with ``buffer +=
line``; there every append may copy the text so far, and an open list of
a megabyte took 0.19 s at 256 KB, 0.74 s at 512 KB and 3.05 s at 1 MB --
four times per doubling. Python 3.13 happened to append in place, so the
same code measured linear there, and the bounded-time test of the outside-
input parsers failed only on CI, under coverage (6.3 s against 5 s).

Whether an append copies depends on whether the allocator can grow the
block where it lies, so the quadratic case does not show on every run.
Growth is therefore judged as a RATIO over a wide span -- eight times the
lines, best of three each: linear gives about x8, quadratic about x64 --
not as a wall-clock bound, which a slow runner breaks alike for both.
"""
from __future__ import annotations

import time

from delfin import config


def _best_of_three(n_lines: int) -> float:
    text = "orca_settings=[\n" + "  1,\n" * n_lines + "]\n"
    best = float("inf")
    for _ in range(3):
        started = time.perf_counter()
        config.parse_control_text(text)
        best = min(best, time.perf_counter() - started)
    return best


def test_a_long_list_is_read_in_time_linear_in_its_length():
    small, large = _best_of_three(50_000), _best_of_three(400_000)
    ratio = large / max(small, 1e-6)
    assert ratio < 24, (f"x{ratio:.1f} for eight times the lines "
                        f"({small:.3f} s -> {large:.3f} s): quadratic")


def test_the_list_still_reads_as_before():
    text = "orca_settings=[\n  'a',  # the first\n  'b',\n]\n"
    assert config.parse_control_text(text).get("orca_settings") == ["a", "b"]
