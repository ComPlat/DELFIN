"""The chain index is a state file like the others.

Input: chain_index.record() with and without DELFIN_STATE set.

Before: the index was always written to ~/.delfin, so DELFIN_STATE (the
override every other state writer honours) did not move it, and its
directory was created under the process umask.
Now it lives under state_paths.state_root(), and the directory and file
are owner-only.
"""
from __future__ import annotations

import os
import stat

from delfin.tools import chain_index


def _mode(path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def test_delfin_state_moves_the_index(tmp_path, monkeypatch):
    state = tmp_path / "state"
    monkeypatch.setenv("DELFIN_STATE", str(state))
    assert chain_index.record("c1", tmp_path / "c1")
    assert chain_index.index_path() == state / chain_index.INDEX_FILENAME
    assert [e["name"] for e in chain_index.load()] == ["c1"]


def test_without_override_it_stays_in_the_home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    monkeypatch.delenv("DELFIN_STATE", raising=False)
    monkeypatch.setattr("pathlib.Path.home", lambda: home)
    assert chain_index.index_path() == home / ".delfin" / "chain_index.jsonl"


def test_index_and_its_directory_are_owner_only(tmp_path, monkeypatch):
    state = tmp_path / "state"
    monkeypatch.setenv("DELFIN_STATE", str(state))
    old = os.umask(0o002)
    try:
        assert chain_index.record("c1", tmp_path / "c1")
    finally:
        os.umask(old)
    assert _mode(state) == 0o700
    assert _mode(chain_index.index_path()) == 0o600
