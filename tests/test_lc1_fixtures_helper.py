"""LC1 helpers: access to the private fixture store for real ORCA output.

Imported by tests/conftest.py (re-exported there) and directly by the
consumer test files, so `from tests.conftest import ...` (a second,
plain-module import of conftest) is not needed. See LICENSE_GUARD.md,
rule family 2: real licensed-program output lives outside the repo,
addressed through the DELFIN_PRIVATE_FIXTURES environment variable;
tests skip when it is absent and never fail in CI.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

PRIVATE_FIXTURES_ENV = "DELFIN_PRIVATE_FIXTURES"


def private_fixtures_path(rel: str = "") -> Path:
    """Path inside the private fixture store, or raise LookupError."""
    root = os.environ.get(PRIVATE_FIXTURES_ENV, "").strip()
    if not root:
        raise LookupError(
            f"{PRIVATE_FIXTURES_ENV} not set: real ORCA fixtures live in the "
            "private store outside the repository (LICENSE_GUARD.md, rule "
            "family 2); set the variable to the store root to run this test."
        )
    p = Path(root) / rel
    if not p.exists():
        raise LookupError(f"{PRIVATE_FIXTURES_ENV} set but {p} is missing")
    return p


def private_fixtures_or_skip(rel: str) -> Path:
    """Path inside the private store, or SKIP (never fail in CI).

    Called at import time by several consumer modules, where the skip
    must cover the whole module; pytest raises a collection ERROR
    instead of a skip without allow_module_level=True. Inside a test the
    flag is harmless (the module is already collected).
    """
    try:
        return private_fixtures_path(rel)
    except LookupError as exc:
        pytest.skip(str(exc), allow_module_level=True)


def lc1_fixture_dir(rel: str) -> Path:
    """Directory of a real-ORCA fixture: the committed legacy path while it
    still exists, else the private store (SKIP without the env var)."""
    legacy = Path(rel)
    if legacy.is_dir() and any(legacy.iterdir()):
        return legacy
    return private_fixtures_or_skip(rel)
