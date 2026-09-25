"""The record of a kept session names the way back in, token included.

It is therefore never on disk in a form anyone else on the machine
could read -- not the directory, not the file, and not the moment
between writing and setting the mode.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.dashboard import session as S

KID = "aaaa1111-0000-4000-8000-000000000001"


@pytest.fixture(autouse=True)
def clean(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    S._reset_for_tests()
    yield
    S._reset_for_tests()


def _mode(p) -> int:
    return stat.S_IMODE(os.stat(p).st_mode)


@pytest.mark.parametrize("umask", [0o022, 0o002, 0o000])
def test_directory_and_record_are_owner_only_whatever_the_umask(tmp_path, umask):
    old = os.umask(umask)
    try:
        path = S.write_record("uc3n990-ab12", kid=KID, root=str(tmp_path / "kept"))
    finally:
        os.umask(old)
    assert path, S._last_write_error
    assert _mode(Path(path).parent) == 0o700
    assert _mode(path) == 0o600


def test_an_older_wider_directory_is_narrowed(tmp_path):
    kept = tmp_path / "kept"
    kept.mkdir(mode=0o755)
    os.chmod(kept, 0o755)
    assert _mode(kept) == 0o755
    S.write_record("uc3n990-ab12", kid=KID, root=str(kept))
    assert _mode(kept) == 0o700


def test_nothing_readable_is_left_beside_the_record(tmp_path, monkeypatch):
    kept = tmp_path / "kept"
    S.write_record("uc3n990-ab12", kid=KID, root=str(kept))
    leftovers = [p for p in kept.iterdir() if p.suffix == ".tmp" or p.name.startswith(".")]
    assert not leftovers


def test_a_record_that_cannot_be_written_leaves_no_temp_file(tmp_path, monkeypatch):
    kept = tmp_path / "kept"

    def boom(*a, **k):
        raise OSError("disk says no")

    monkeypatch.setattr(os, "replace", boom)
    path = S.write_record("uc3n990-ab12", kid=KID, root=str(kept))
    assert path == ""
    assert "disk says no" in S._last_write_error
    assert not any(kept.iterdir())
