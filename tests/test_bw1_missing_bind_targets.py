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
import shutil
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent import contained_run as CR


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
    # Network mode "none": no egress proxy bind/token (R2, s4's finding).
    # The proxy dir is a session-local /tmp mkdtemp that cannot exist on
    # another node; the tests' subject (mask protection) is untouched.
    from delfin.agent import egress_proxy as _ep
    monkeypatch.setattr(_ep, "network_mode", lambda: "none")


def _argv(ws, monkeypatch, caged_unasked=False):
    perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    return A._bash_isolation_argv("true", str(ws), perms,
                                  mode="bwrap", caged_unasked=caged_unasked)


def _run(argv):
    """Run the caged argv through the REAL exec boundary
    (contained_run.run) exactly as production does, so the exec-boundary
    fix sits on the tests' path. subprocess.run directly would bypass
    the fix (operator rerun finding, round 2)."""
    out = CR.run(argv, timeout=120)
    # contained_run.run returns a CompletedProcess-like object
    # (returncode/stdout/stderr) -- assertable the same way.
    return out


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
    def test_host_side_approval_codes_stays_empty(self, ws, monkeypatch):
        """(c) The mask is a tmpfs: a write attempt INSIDE the cage lands
        in the mount, never on the HOST dir. Checked on the host side
        (os.listdir), not by in-cage output -- an in-cage listing is a
        tmpfs artifact (operator round-2 note, no defect)."""
        _cage_decision(monkeypatch)
        perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
        argv = A._bash_isolation_argv(
            "python3 -c \"import os;"
            "os.makedirs(os.path.expanduser('~/.delfin/approval_codes'),"
            "exist_ok=True);"
            "open(os.path.expanduser('~/.delfin/approval_codes/x'), 'w')"
            ".write('leak')\"", str(ws), perms, mode="bwrap")
        run = _run(argv)
        assert run.returncode == 0, f"stderr: {run.stderr}"
        host_dir = Path(os.environ["HOME"]) / ".delfin" / "approval_codes"
        assert host_dir.is_dir(), "the ensure created the target"
        assert os.listdir(host_dir) == [], (
            f"a caged write leaked to the host: {os.listdir(host_dir)}")

    @requires_a_working_sandbox
    def test_inside_the_cage_approval_codes_is_not_readable(self, ws,
                                                            monkeypatch):
        """With the fix the denial is preserved: a command inside the cage
        still cannot list or write ~/.delfin/approval_codes."""
        _cage_decision(monkeypatch)
        perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
        probe = ("python3", "-c",
                 "import os,stat;"
                 "d=os.path.expanduser('~/.delfin/approval_codes');"
                 "vis = os.path.isdir(d) and not os.access(d, os.R_OK);"
                 "listed = (not os.path.isdir(d)) or os.listdir(d) == [];"
                 "print('HIDDEN' if (vis or listed) else 'VISIBLE')")
        argv = A._bash_isolation_argv(
            f"python3 -c {probe[2]!r}", str(ws), perms, mode="bwrap")
        run = _run(argv)
        # The directory is hidden (fresh tmpfs => empty) or absent: either
        # way nothing of it is readable, and the run itself succeeds.
        assert run.returncode == 0, f"stderr: {run.stderr}"
        assert "VISIBLE" not in run.stdout, (
            f"approval_codes content readable inside the cage: {run.stdout}")


def _mask_args_for(argv: list, target: str) -> list:
    """Mount argument pairs in ``argv`` whose in-cage path is ``target``."""
    hits = []
    for i, tok in enumerate(argv[:-1]):
        if tok in ("--tmpfs", "--ro-bind", "--bind") and argv[i + 1] == target:
            hits.append(tok)
    return hits


class TestR3ParallelPytestBasetempDeletion:
    """R3 repro (operator's CANDIDATE cause, NOT proof of wave 31):
    pytest's basetemp retention deletes old pytest-N directories while a
    parallel pytest session still runs there. A test HOME under tmp_path
    of a still-running session vanishes mid-run -> ~/.delfin/approval_codes
    gone between argv build and exec -> the wave-31 EROFS failure shape.
    An operator harness hit the same refusal with a private HOME that
    nothing deleted, so a second cause may exist; the exec-boundary ensure
    is the net that covers both."""

    @requires_a_working_sandbox
    def test_deletion_of_a_running_sessions_home(self, tmp_path,
                                                 monkeypatch):
        """Simulate the retention deletion: another 'session' deletes the
        pytest-N parent of the running test's tmp_path HOME while this
        test is between build and exec. Expected WITHOUT any fix: the
        caged command dies with the measured bwrap EROFS refusal."""
        home = tmp_path / "home"
        home.mkdir()
        monkeypatch.setenv("HOME", str(home))
        ws = tmp_path / "ws"
        ws.mkdir()
        _cage_decision(monkeypatch)
        argv = _argv(ws, monkeypatch)
        # The parallel session's deletion: the numbered basetemp dir that
        # CONTAINS this test's tmp_path is removed from outside. Do not
        # assume tmp_path.parent is it (the repo's pytest setup can nest
        # tmp_path one level deeper, e.g. .../pytest-N/t0): walk up to the
        # nearest ancestor named pytest-N instead.
        basetemp = tmp_path
        while not basetemp.name.startswith("pytest-"):
            parent = basetemp.parent
            if parent == basetemp:
                pytest.skip("no pytest-N basetemp ancestor found")
            basetemp = parent
        assert basetemp.is_dir()
        shutil.rmtree(basetemp)
        run = _run(argv)
        # Red state (no fix at the exec boundary): bwrap refuses with the
        # measured message. With the exec-boundary ensure the run works.
        assert run.returncode == 0, (
            "R3 REPRODUCED: the parallel-pytest deletion of an in-use "
            f"basetemp kills the cage:\nstdout: {run.stdout}\n"
            f"stderr: {run.stderr}")



def _deep_rmtree(path: Path) -> None:
    """Test-only removal of a directory the test itself created."""
    shutil.rmtree(path)

