"""BW1 red test for the EXEC-BOUNDARY fix (px pair contained_run.py): the
one place that starts a caged bwrap argv must call
cage_binds.ensure_argv_bind_targets right before Popen.

Per the operator's round-2 verdict: the build-time api_client wiring is
DROPPED -- only the exec boundary can observe the state bwrap actually
sees. Without the contained_run change, a mask target that vanishes
between argv build and exec reaches bwrap and bwrap dies with
``Can't create file ...: Read-only file system`` -- the measured defect
(and the R3 basetemp scenario). With the change, the target exists again
at exec.

Real call site (``contained_run.run``), real bwrap; skips where bwrap
cannot work (nested cage -- the operator runs these). Red on the branch
without the applied pair; green once the operator builds it."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import contained_run as CR


@requires_a_working_sandbox
def test_vanished_approval_codes_does_not_kill_the_cage(tmp_path,
                                                        monkeypatch):
    """The measured defect, through the REAL exec site: remove
    ~/.delfin/approval_codes between build and exec; the cage must still
    run the command (the exec-boundary ensure re-creates the target)."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = tmp_path / "ws"
    ws.mkdir()
    argv, _ = _argv(ws, monkeypatch)
    codes = home / ".delfin" / "approval_codes"
    if codes.is_dir():
        shutil.rmtree(codes)
    run = _run_caged(argv, monkeypatch, ws)
    assert run.returncode == 0, (
        "the cage died on a mask target that vanished between build and "
        f"exec (the measured EROFS defect):\nstderr: {run.stderr}")


def _argv(ws, monkeypatch):
    """Build the bwrap argv the way the real builder does (the same
    decision points _cage_decision pins in the sibling test file)."""
    from delfin.agent import api_client as A
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    from delfin.agent import egress_proxy as _ep
    monkeypatch.setattr(_ep, "network_mode", lambda: "none")
    perms = A.KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    return A._bash_isolation_argv("true", str(ws), perms, mode="bwrap"), None


def _run_caged(argv, monkeypatch, ws):
    """Run the argv through the REAL exec boundary (contained_run.run),
    not a bare Popen -- that is the call site the px patch wires.
    ``argv[0]`` is "bwrap", which the patch's guard matches."""
    rc_out = CR.run(argv, cwd=str(ws), timeout=120)
    return rc_out
