"""One long line made the grounding scanner take minutes.

``_LOC_PATH`` ends in a greedy component that then has to find the
extension dot. On a long run of word characters containing no dot, the
engine gives back one character at a time from every starting position --
quadratic. Six patterns use it, and the scanner runs on every final
answer and every delegate report.

Measured before the bound, on one unbroken line:

    20 kB    6.6 s
    50 kB   46.5 s
   100 kB   minutes

Ordinary multi-line text of the same size cost 16 ms, which is why this
never showed up in normal use. It needed a pasted blob, a minified dump,
a base64 payload -- input nobody anticipates, where the only symptom is a
turn that stops responding and no error anywhere.

Bounding each path component to 200 characters makes the work per
position constant. That is far past any real path, so the set of strings
matched is unchanged; the tests below check the detection still holds as
well as the time.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import verify_guard as vg


_BUDGET_S = 2.0


def _scan(text: str):
    return vg.scan_for_ungrounded_location_claims(
        text, observed_files=frozenset(), max_flags=6)


@pytest.mark.parametrize("size", [20_000, 100_000])
def test_one_long_line_is_scanned_quickly(size):
    start = time.time()
    _scan("x" * size)
    assert time.time() - start < _BUDGET_S


def test_a_long_line_of_path_like_text_is_also_quick():
    """Word characters, slashes and dots, but never a valid extension --
    the shape that maximises the backtracking."""
    start = time.time()
    _scan("delfin/agent/" + ("abc.def/" * 2000))
    assert time.time() - start < _BUDGET_S


def test_ordinary_text_is_unaffected():
    start = time.time()
    _scan("a finding in delfin/agent/api_client.py:1234\n" * 500)
    assert time.time() - start < _BUDGET_S


# ---------------------------------------------------------------------------
# ...and it still finds what it is for
# ---------------------------------------------------------------------------

def test_a_nested_path_with_a_line_is_still_flagged():
    flags = _scan("The limit is set in delfin/agent/api_client.py:1234.")
    assert [(f.path, f.line) for f in flags] == [
        ("delfin/agent/api_client.py", 1234)]


def test_a_bare_filename_with_a_line_is_still_flagged():
    flags = _scan("See engine.py:99 for the loop.")
    assert [(f.path, f.line) for f in flags] == [("engine.py", 99)]


def test_a_dotted_filename_still_matches():
    flags = _scan("Config lives in delfin/agent/settings.local.json:12.")
    assert [f.path for f in flags] == ["delfin/agent/settings.local.json"]


def test_a_hyphenated_directory_still_matches():
    flags = _scan("Written in my-pkg/sub-dir/mod_a.py:7.")
    assert [f.path for f in flags] == ["my-pkg/sub-dir/mod_a.py"]


def test_a_dotted_module_name_is_still_not_a_path():
    """The extension set exists so `numpy.linalg` never reads as a file."""
    assert _scan("numpy.linalg.norm is used at 12.") == []


def test_an_unknown_extension_is_still_not_a_path():
    assert _scan("See archive.tar:12 for the bundle.") == []
