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


def test_plain_bind_source_is_never_created(home):
    # Operator ruling, round 2: the ensure creates ONLY mask targets.
    # A missing plain-bind source is bwrap's own diagnosable refusal; an
    # empty dir made in its place would silently run the command against
    # an empty workspace/undo store.
    src = home / "state" / "missing_workspace"
    argv = ["bwrap", "--ro-bind", str(src), str(src), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert not src.exists()


def test_try_bind_source_is_never_created(home):
    # Same rule for -try forms: bwrap skips a missing source by contract,
    # and creating it would turn the skip into a bind of a fresh empty
    # dir. Masks have no -try form, so nothing -try is created at all.
    d = home / ".delfin" / "messenger_pairing"
    argv = ["bwrap", "--ro-bind-try", str(d), str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert not d.exists()


def test_file_mask_destination_is_created_as_a_file(home):
    f = home / ".delfin" / "credentials.json"
    argv = ["bwrap", "--ro-bind", "/dev/null", str(f), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert f.is_file()
    assert stat.S_IMODE(f.stat().st_mode) == 0o600


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


def test_dev_bind_operand_is_never_created(home):
    # --dev-bind is a plain bind (mask rule does not apply): a missing
    # operand stays missing -- bwrap's own refusal is the correct outcome.
    d = home / "run" / "user"
    argv = ["bwrap", "--dev-bind", "/dev", str(d), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert not d.exists()


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


def test_file_mask_destination_parents_are_created(home):
    # The file-mask rule creates the DESTINATION (a file) -- its parent
    # chain comes with it; the /dev/null source itself is of course not
    # touched.
    f = home / "private" / "deep" / "mask"
    argv = ["bwrap", "--ro-bind", "/dev/null", str(f), "x"]
    CB.ensure_argv_bind_targets(argv)
    assert f.is_file()


def test_a_mask_target_outside_home_is_owner_only_too(tmp_path, monkeypatch):
    """No created mask target is group- or world-readable (CodeQL
    py/overly-permissive-file): outside HOME as well as inside it."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    target = tmp_path / "elsewhere" / "masked"
    CB.ensure_argv_bind_targets(["bwrap", "--tmpfs", str(target), "--", "true"])
    assert stat.S_IMODE(target.stat().st_mode) == 0o700
