"""Property-based tests keep their cache out of the checkout.

Input: the test session's environment and hypothesis's storage location,
and the checkout-leak guard in tests/conftest.py.

Output:
- HYPOTHESIS_STORAGE_DIRECTORY is set for the session (and for every
  process a test starts) to a directory outside the checkout;
- hypothesis itself, where installed, resolves its storage there, not to
  ``<cwd>/.hypothesis``;
- ``.hypothesis`` is no longer build noise to the checkout guard, so a
  cache written into the checkout fails the run.

Measured before the change: one run of
tests/test_outside_input_parsers_hold_on_any_text.py from the checkout
root left ``.hypothesis/`` (unicode_data/, constants/, a .gitignore)
beside pyproject.toml, which the guard skipped by name.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _outside_checkout(path: Path) -> bool:
    try:
        path.resolve().relative_to(REPO.resolve())
    except ValueError:
        return True
    return False


def test_the_session_points_hypothesis_outside_the_checkout():
    where = os.environ.get("HYPOTHESIS_STORAGE_DIRECTORY", "")
    assert where, "HYPOTHESIS_STORAGE_DIRECTORY is not set for the test session"
    assert _outside_checkout(Path(where)), where


def test_hypothesis_resolves_its_storage_there():
    pytest.importorskip("hypothesis")
    from hypothesis.configuration import storage_directory

    path = storage_directory("unicode_data", intent_to_write=False).path
    assert _outside_checkout(path), path


def test_a_cache_in_the_checkout_is_not_build_noise_to_the_guard():
    text = (REPO / "tests" / "conftest.py").read_text(encoding="utf-8")
    noise = re.search(r"_CHECKOUT_NOISE = \((.*?)\)\n", text, re.S).group(1)
    assert '".hypothesis"' not in noise
