"""_run_login_shell must not be able to recurse.

A login shell's module initialization can itself spawn `bash -lc`
(~75 nested `bash -lc module -t avail` on the login node, 2026-09-26,
after a gate run of the analysis-probe tests). The depth guard passes a
marker environment to the probe shell; the marker is inherited by every
child, so a nested probe (this test simulates the recursive side by
calling _run_login_shell with the marker already set) refuses to start
instead of adding another level.
"""

from __future__ import annotations

import os
import subprocess
from unittest import mock

from delfin.system_tools import _run_login_shell


def test_a_nested_login_shell_probe_refuses_to_start(monkeypatch):
    """With the marker in the environment (what every child of a probe
    shell sees), no new bash is started at all."""
    monkeypatch.setenv("_DELFIN_LOGIN_SHELL_PROBE", "1")
    calls = []
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: calls.append(a) or mock.DEFAULT, raising=True)
    result = _run_login_shell("module -t avail")
    assert result is None, "the nested probe returned a result anyway"
    assert not calls, "the nested probe still started a subprocess"


def test_the_first_probe_marks_the_child_environment(monkeypatch):
    """The first probe must pass a copy of the environment carrying the
    marker, so that anything the login shell starts sees it."""
    seen_env = {}

    def fake_run(cmd, **kwargs):
        seen_env.update(kwargs.get("env") or {})
        return subprocess.CompletedProcess(cmd, 0, "stub", "")

    monkeypatch.delenv("_DELFIN_LOGIN_SHELL_PROBE", raising=False)
    monkeypatch.setattr(subprocess, "run", fake_run)
    _run_login_shell("module -t avail")
    assert seen_env.get("_DELFIN_LOGIN_SHELL_PROBE") == "1", (
        "the probe child does not carry the recursion marker"
    )
    # And the parent process env stays untouched (no leak upward).
    assert "_DELFIN_LOGIN_SHELL_PROBE" not in os.environ
