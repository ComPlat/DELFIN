"""The login-shell probe always finds the system tools.

Input: ``system_tools._run_login_shell`` started from a process whose
PATH lacks the system directories (a test gate's reduced tool directory).

Before: the probe shell got that PATH unchanged. Its profile called a
tool it could not find, bash ran ``command_not_found_handle``, and the
handler prints through ``$(gettext ...)`` -- gettext was not found
either, so the handler called itself in a new subshell, again and again:
815 nested ``bash -lc module -t avail`` from one probe (2026-11-14).
Now /usr/bin and /bin are appended when missing; the caller's own
directories keep their order and come first.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from delfin import system_tools
from delfin.system_tools import _run_login_shell


class _FakePopen:
    seen: dict = {}

    def __init__(self, cmd, **kwargs):
        _FakePopen.seen = dict(kwargs.get("env") or {})
        self.args, self.returncode, self.pid = cmd, 0, -1

    def communicate(self, timeout=None):
        return "", ""


@pytest.fixture()
def fake_popen(monkeypatch):
    monkeypatch.delenv("_DELFIN_LOGIN_SHELL_PROBE", raising=False)
    monkeypatch.setattr(subprocess, "Popen", _FakePopen)
    _FakePopen.seen = {}
    return _FakePopen


def test_a_reduced_path_gets_the_system_dirs(monkeypatch, fake_popen):
    monkeypatch.setenv("PATH", "/gate/cipath")
    _run_login_shell("module -t avail")
    parts = fake_popen.seen["PATH"].split(os.pathsep)
    assert parts[0] == "/gate/cipath"
    assert "/usr/bin" in parts and "/bin" in parts


def test_a_full_path_stays_as_it_was(monkeypatch, fake_popen):
    full = os.pathsep.join(["/opt/x/bin", "/usr/bin", "/bin"])
    monkeypatch.setenv("PATH", full)
    _run_login_shell("module -t avail")
    assert fake_popen.seen["PATH"] == full


def test_an_empty_path_is_not_left_empty(monkeypatch, fake_popen):
    monkeypatch.setenv("PATH", "")
    _run_login_shell("module -t avail")
    assert fake_popen.seen["PATH"] == os.pathsep.join(system_tools._SYSTEM_DIRS)


def test_the_callers_path_is_not_changed(monkeypatch, fake_popen):
    monkeypatch.setenv("PATH", "/gate/cipath")
    _run_login_shell("module -t avail")
    assert os.environ["PATH"] == "/gate/cipath"
