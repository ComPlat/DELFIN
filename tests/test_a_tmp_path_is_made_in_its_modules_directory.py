"""``tmp_path`` is created in one numbered directory per module.

The layout contract of the override in ``tests/conftest.py``:
``<basetemp>/t<k>/<pytest's name><n>``. The base directory gains one
entry per module instead of one per test, which is what keeps the
creation of a ``tmp_path`` from scanning every earlier test's directory
(O(N^2) over a 29 000-test run).
"""

from __future__ import annotations

import re

import pytest

_SEEN: dict = {}


def test_the_directory_sits_one_level_below_the_base(tmp_path,
                                                     tmp_path_factory):
    base = tmp_path_factory.getbasetemp().resolve()
    here = tmp_path.resolve()
    assert here.parent.parent == base
    assert re.fullmatch(r"t\d+", here.parent.name)
    assert re.fullmatch(
        r"test_the_directory_sits_one_le\d+", here.name), here.name
    assert here.is_dir() and not any(here.iterdir())


@pytest.mark.parametrize("which", ["first", "second"])
def test_tests_of_one_module_share_its_directory(tmp_path, which):
    _SEEN[which] = tmp_path
    if which == "second":
        assert _SEEN["first"].parent == tmp_path.parent
        assert _SEEN["first"] != tmp_path


def test_the_factory_still_numbers_in_the_base(tmp_path_factory):
    made = tmp_path_factory.mktemp("factory_made")
    assert made.parent == tmp_path_factory.getbasetemp()
    assert re.fullmatch(r"factory_made\d+", made.name)
