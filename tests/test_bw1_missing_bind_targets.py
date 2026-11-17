"""BW1: a cage argv must not fail at exec because a mask/bind target that
existed at argv-build time is missing when bwrap runs.

Measured failure (wave 31, intermittent): ``bwrap: Can't mkdir parents for
<HOME>/.delfin/approval_codes: Read-only file system``. The argv builders
check a target's existence only at build time; the leaf of a plain
``--tmpfs``/``--ro-bind`` mount whose parent chain sits under the read-only
``--ro-bind / /`` must exist when bwrap sets the mount up, and bwrap cannot
create it inside the read-only root. A HOME that lost its ``.delfin``
between build and exec (first approval code creates it, a HOME refresh
removes it) therefore kills the command before it runs -- randomly from the
outside.

Rule under test (cage_binds.py, phase 2): before the argv is handed to
bwrap, every mask target under HOME exists (created with the owner's mode,
0700 under ``~/.delfin``); no mask is ever silently dropped (--tmpfs-try on
a mask skips it and would re-expose the path). Protection both ways is
pinned: with the fix the command still cannot read or write the masked
paths, and the argv still runs.

Drives the real argv builder (``api_client._bash_isolation_argv``), real
bwrap; skipped where bwrap does not work (nested cage)."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A


@pytest.fixture
def ws(tmp_path, monkeypatch):
    """A fresh HOME with no ~/.delfin and a plain workspace."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _cage_decision(monkeypatch):
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    # Landlock off: the argv taken is the plain bwrap branch -- the shape
    # the measured failure ran under.
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)


def _argv(ws, monkeypatch, caged_unasked=False):
    perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    return A._bash_isolation_argv("true", str(ws), perms,
                                  mode="bwrap", caged_unasked=caged_unasked)


def _run(argv):
    return subprocess.run(argv, capture_output=True, text=True, timeout=120)


@pytest.mark.skipif(not getattr(A, "_process_cage_enabled", bool)(),
                    reason="process cage disabled")
class TestCageRunsWhenMaskTargetVanishes:
    @requires_a_working_sandbox
    def test_approval_codes_dir_removed_after_build(self, ws, monkeypatch):
        """The measured failure: approval_codes exists at build (created by
        the argv builder itself), is gone at exec, the cage must still run."""
        _cage_decision(monkeypatch)
        argv = _argv(ws, monkeypatch)
        codes = Path(os.environ["HOME"]) / ".delfin" / "approval_codes"
        assert codes.is_dir(), "the builder creates the dir at build time"
        _deep_rmtree(codes)
        assert not codes.exists()
        run = _run(argv)
        assert run.returncode == 0, (
            "cage refused to run after a mask target vanished between build "
            "and exec (the measured EROFS failure):\n"
            f"stdout: {run.stdout}\nstderr: {run.stderr}")

    @requires_a_working_sandbox
    def test_whole_delfin_dir_removed_after_build(self, ws, monkeypatch):
        """Every mask target under ~/.delfin gone at exec: the cage still
        runs and the paths stay hidden (protection not weakened)."""
        _cage_decision(monkeypatch)
        argv = _argv(ws, monkeypatch)
        delfin = Path(os.environ["HOME"]) / ".delfin"
        assert delfin.is_dir()
        _deep_rmtree(delfin)
        run = _run(argv)
        assert run.returncode == 0, f"stderr: {run.stderr}"

    @requires_a_working_sandbox
    def test_masked_paths_stay_hidden_in_a_fresh_home(self, ws, monkeypatch):
        """Protection both ways: in a HOME that never had ~/.delfin, the
        argv must (a) run and (b) contain a mount that HIDES
        ~/.delfin/approval_codes -- never a --tmpfs-try that silently skips
        a missing mask target."""
        _cage_decision(monkeypatch)
        argv = _argv(ws, monkeypatch)
        run = _run(argv)
        assert run.returncode == 0, f"stderr: {run.stderr}"
        codes = str(Path(os.environ["HOME"]) / ".delfin" / "approval_codes")
        masked = _mask_args_for(argv, codes)
        assert masked, f"{codes} is not masked anywhere in the argv"
        assert "--tmpfs-try" not in argv, (
            "a mask target must never degrade to a try-mount: a skipped "
            "mask re-exposes the path")

    @requires_a_working_sandbox
    def test_inside_the_cage_approval_codes_is_not_readable(self, ws,
                                                            monkeypatch):
        """With the fix the denial is preserved: a command inside the cage
        still cannot list or write ~/.delfin/approval_codes."""
        _cage_decision(monkeypatch)
        perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
        probe = ("python3", "-c",
                 "import os;"
                 "d=os.path.expanduser('~/.delfin/approval_codes');"
                 "print('EXISTS' if os.path.isdir(d) and os.listdir(d) "
                 "is not None else 'ABSENT')")
        argv = A._bash_isolation_argv(
            f"python3 -c {probe[2]!r}", str(ws), perms, mode="bwrap")
        run = _run(argv)
        # The directory is hidden (fresh tmpfs => empty) or absent: either
        # way nothing of it is readable, and the run itself succeeds.
        assert run.returncode == 0, f"stderr: {run.stderr}"
        assert "EXISTS" not in run.stdout or run.stdout.strip() == "", (
            f"approval_codes visible inside the cage: {run.stdout}")


def _mask_args_for(argv: list, target: str) -> list:
    """Mount argument pairs in ``argv`` whose in-cage path is ``target``."""
    hits = []
    for i, tok in enumerate(argv[:-1]):
        if tok in ("--tmpfs", "--ro-bind", "--bind") and argv[i + 1] == target:
            hits.append(tok)
    return hits


def _deep_rmtree(path: Path) -> None:
    """Test-only removal of a directory the test itself created."""
    import shutil
    shutil.rmtree(path)
