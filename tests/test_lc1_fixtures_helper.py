"""LC1 helpers: access to the private fixture store for real ORCA output.

Imported by tests/conftest.py (re-exported there) and directly by the
consumer test files, so `from tests.conftest import ...` (a second,
plain-module import of conftest) is not needed. See LICENSE_GUARD.md,
rule family 2: real licensed-program output lives outside the repo,
addressed through the DELFIN_PRIVATE_FIXTURES environment variable (or
the paths.fixtures_dir settings entry); tests skip when no store is
configured and never fail in CI. Files are hash-checked against the
store manifest (delfin/fixtures_store.py); a modified file is refused,
never silently used.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.fixtures_store import private_fixtures_path as _store_path
from delfin.fixtures_store import store_root as _store_root

PRIVATE_FIXTURES_ENV = "DELFIN_PRIVATE_FIXTURES"


def private_fixtures_path(rel: str = "") -> Path:
    """Hash-checked path inside the private fixture store, or raise LookupError."""
    return _store_path(rel)


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


def private_fixtures_store_root() -> Path:
    """The store root, or raise LookupError (env var, then settings)."""
    return _store_root()
