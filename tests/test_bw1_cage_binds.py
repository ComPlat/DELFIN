"""BW1 unit tests: cage_binds.ensure_argv_bind_targets.

Runs everywhere (no bwrap): the module walks a finished cage argv and
creates every missing bind target with the owner's mode before bwrap runs.
Pins, per QS ruling: mask targets are created (0700 under ~/.delfin), never
degraded to a -try form; file masks get a file destination; the command
tail is untouched; a non-createable target is left alone (bwrap's own
refusal, diagnosable) and the call never raises."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.agent import cage_binds as CB


@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    return h


def test_missing_approval_codes_dir_is_created_0700(home):
    argv = ["bwrap", "--ro-bind", "/", "/", "--tmpfs",
            str(home / ".delfin" / "approval_codes"), "true"]
    out = CB.ensure_argv_bind_targets(argv)
    d = home / ".delfin" / "approval_codes"
    assert d.is_dir()
    assert stat.S_IMODE(d.stat().st_mode) == 0o700
    assert out is argv  # the argv is unchanged, only the filesystem acted


def test_missing_parent_chain_is_created(home):
    argv = ["bwrap", "--tmpfs", str(home / ".delfin" / "trust_store"), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert (home / ".delfin" / "trust_store").is_dir()


def test_existing_target_is_not_touched(home):
    d = home / ".aws"
    d.mkdir(mode=0o750)
    argv = ["bwrap", "--tmpfs", str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert stat.S_IMODE(d.stat().st_mode) == 0o750


def test_target_outside_delfin_gets_0755(tmp_path):
    d = tmp_path / "somewhere" / "deep"
    argv = ["bwrap", "--bind", str(d), str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert d.is_dir()
    assert stat.S_IMODE(d.stat().st_mode) == 0o755


def test_try_bind_targets_are_created_too(home):
    # -try forms are part of BIND_FLAGS: creating the target never weakens
    # them (a created target means the bind happens instead of being
    # skipped -- strictly tighter).
    d = home / ".delfin" / "messenger_pairing"
    argv = ["bwrap", "--ro-bind-try", str(d), str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert d.is_dir()


def test_file_mask_destination_is_created_as_a_file(home):
    f = home / ".delfin" / "credentials.json"
    argv = ["bwrap", "--ro-bind", "/dev/null", str(f), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert f.is_file()
    assert stat.S_IMODE(f.stat().st_mode) == 0o700


def test_command_tail_is_untouched():
    argv = ["bwrap", "--ro-bind", "/", "/", "--", "sh", "-c", "true",
            "--tmpfs", "not-a-flag"]
    before = list(argv)
    CB.ensure_argv_bind_targets(argv)
    assert argv == before


def test_non_creatable_target_does_not_raise(tmp_path, monkeypatch):
    locked = tmp_path / "locked"
    locked.mkdir()
    os.chmod(locked, 0o555)
    argv = ["bwrap", "--tmpfs", str(locked / "sub" / "deeper"), "x"]
    out = CB.ensure_argv_bind_targets(argv)
    assert out == argv
    os.chmod(locked, 0o755)  # so tmp cleanup works


def test_relative_operands_are_ignored():
    argv = ["bwrap", "--ro-bind", "rel", "paths", "true"]
    before = list(argv)
    CB.ensure_argv_bind_targets(argv)
    assert argv == before


def test_dev_bind_target_created(home):
    d = home / "run" / "user"
    argv = ["bwrap", "--dev-bind", "/dev", str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert d.is_dir()
    assert stat.S_IMODE(d.stat().st_mode) == 0o755


def test_home_door_target_gets_0700(home):
    # s4's live repro: ~/.ssh is a door target of _process_cage_options
    # (plain --tmpfs, build-time isdir). A created door is owner-only.
    argv = ["bwrap", "--tmpfs", str(home / ".ssh"), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert (home / ".ssh").is_dir()
    assert stat.S_IMODE((home / ".ssh").stat().st_mode) == 0o700


def test_masks_have_no_try_form():
    # The degradation QS rejected: there is no --tmpfs-try in the walked
    # vocabulary, and DELFIN's argv must never contain one -- a skipped
    # try-mask re-exposes the path it was to hide.
    for flag in CB.BIND_FLAGS:
        assert not flag.endswith("-try") or flag in (
            "--ro-bind-try", "--bind-try")
    assert "--tmpfs-try" not in CB.BIND_FLAGS
    assert CB.MASK_FLAGS == ("--tmpfs",)


def test_try_bind_missing_source_is_not_created(home):
    # s4 F1: bwrap's contract for a -try bind is SKIP when the source is
    # missing. Creating that source would turn the skip into a bind of a
    # fresh empty dir; the destination is still ensured.
    src = home / "state" / "missing_src"
    dst = home / "state" / "dst"
    argv = ["bwrap", "--ro-bind-try", str(src), str(dst), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert not src.exists()      # the skip semantics for src are kept
    assert dst.is_dir()          # the target side is still made to exist


def test_plain_bind_missing_source_is_created(home):
    # Contrast (the DELFIN form, src==dst): a plain bind with a missing
    # source gets it created -- bwrap would refuse, not skip.
    src = home / "state" / "missing_plain"
    argv = ["bwrap", "--ro-bind", str(src), str(src), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert src.is_dir()
