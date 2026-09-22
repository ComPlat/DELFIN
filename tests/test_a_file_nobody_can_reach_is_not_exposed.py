"""A 0644 file inside a 0700 directory is not readable by anybody else.

The doctor reported, on a developer machine:

    OPENAI_API_KEY sits in ~/.bashrc:135, which group, others can read
    (mode 0644)

True about the file, and wrong about the world. The home directory above
it is ``drwx------``: on Unix a file is reached only by walking the
directories above it, and each of those needs its execute bit for the
asker. Nobody but the owner may enter, so nobody but the owner can open
the file, whatever its own mode says. The nightly backup preserves the
same 0700, so the copy is as closed as the original.

A 0700 home is the ordinary configuration, which makes this a warning
that fires on almost every host while naming a hole that is shut. That is
worse than saying nothing: a security check people learn to scroll past
has stopped being a security check.

So the question is reachability, not the last mode bits — and it is asked
of the whole path.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import credentials as C


LINE = 'export OPENAI_API_KEY="sk-not-a-real-key-0123456789abcdef"\n'


def _home_with_key(tmp_path, *, home_mode: int, file_mode: int):
    home = tmp_path / "home"
    home.mkdir()
    rc = home / ".bashrc"
    rc.write_text(LINE, encoding="utf-8")
    os.chmod(rc, file_mode)
    os.chmod(home, home_mode)
    return home, rc


def _found(home):
    return C.keys_in_readable_shell_files(("OPENAI_API_KEY",), home=home)


# -- the report -------------------------------------------------------------

def test_a_readable_file_in_a_private_directory_is_not_reported(tmp_path):
    """The case from the machine: 0644 file, 0700 home."""
    home, _rc = _home_with_key(tmp_path, home_mode=0o700, file_mode=0o644)
    try:
        assert _found(home) == [], (
            "nobody can enter the directory, so nobody can open the file")
    finally:
        os.chmod(home, 0o700)


def test_a_group_readable_file_in_a_private_directory_is_not_reported(tmp_path):
    home, _rc = _home_with_key(tmp_path, home_mode=0o700, file_mode=0o640)
    assert _found(home) == []


# -- what must still be caught ---------------------------------------------

def test_a_readable_file_in_a_traversable_directory_is_reported(a_path_others_can_walk):
    """The real hole, which this must not stop finding."""
    tmp_path = a_path_others_can_walk
    home, rc = _home_with_key(tmp_path, home_mode=0o755, file_mode=0o644)
    rows = _found(home)
    assert len(rows) == 1, rows
    row = rows[0]
    assert row["name"] == "OPENAI_API_KEY"
    assert row["file"] == str(rc)
    assert "others" in row["readable_by"]


def test_the_value_is_never_carried(a_path_others_can_walk):
    tmp_path = a_path_others_can_walk
    home, _rc = _home_with_key(tmp_path, home_mode=0o755, file_mode=0o644)
    blob = repr(_found(home))
    assert "sk-not-a-real-key" not in blob, blob


def test_a_private_file_in_an_open_directory_is_not_reported(tmp_path):
    """Already true before, and it stays true."""
    home, _rc = _home_with_key(tmp_path, home_mode=0o755, file_mode=0o600)
    assert _found(home) == []


def test_a_group_readable_file_under_a_group_traversable_directory_is_reported(
        a_path_others_can_walk):
    tmp_path = a_path_others_can_walk
    home, _rc = _home_with_key(tmp_path, home_mode=0o750, file_mode=0o640)
    rows = _found(home)
    assert len(rows) == 1, rows
    assert "group" in rows[0]["readable_by"]


# -- and it never falls over ------------------------------------------------

def test_a_path_that_cannot_be_stat_ed_is_quiet(tmp_path):
    assert C.keys_in_readable_shell_files(
        ("OPENAI_API_KEY",), home=tmp_path / "nothing-here") == []


def test_it_asks_about_the_whole_path_not_only_the_parent(tmp_path):
    """One closed directory anywhere above is enough to shut the door."""
    outer = tmp_path / "outer"
    outer.mkdir()
    home = outer / "inner"
    home.mkdir()
    rc = home / ".bashrc"
    rc.write_text(LINE, encoding="utf-8")
    os.chmod(rc, 0o644)
    os.chmod(home, 0o755)
    os.chmod(outer, 0o700)          # shut two levels up
    try:
        assert _found(home) == []
    finally:
        os.chmod(outer, 0o755)
