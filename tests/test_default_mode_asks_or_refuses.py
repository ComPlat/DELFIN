"""What ``default`` does at the gate, attended and with nobody to ask.

``default`` is now the start posture when neither ``--permission-mode``
nor a settings ``default_mode`` names one, for a terminal session as for
``chat -p`` and ``run``. These tests pin what that posture lets through,
through the real executor and permission object:

* attended (a confirm callback bound -- the terminal broker, or the
  approvals broker of ``--supervised``): a write and a shell command that
  is not on the auto-allow list are both put to the callback first, and a
  "no" leaves the disk untouched;
* nobody to ask (no callback -- ``-p`` / ``run`` without
  ``--supervised``): a shell command off the auto-allow list is refused
  and not run. A write inside the workspace -- a file tool, or a shell
  write such as ``touch`` -- is ALLOWED, by the design documented at
  ``_gate_write_path`` (writing has no allow list, so refusing it would
  leave every unattended run unable to produce anything), and pinned here
  so that a change to it is a decision rather than an accident;
* ``plan`` still refuses both.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from delfin.agent import api_client as _A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture(autouse=True)
def _uncaged(monkeypatch):
    """These tests pin default mode's question for an off-list command.
    Where bubblewrap and Landlock hold it, such a command runs unasked
    (_caged_command_may_run_unasked, pinned in its own test file); without
    this the verdict would depend on the host."""
    monkeypatch.setattr(_A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)


class _Asker:
    """Bound-method callback, as the brokers hand over; always says no."""

    last_timed_out = False

    def __init__(self):
        self.asked: list[str] = []

    def callback(self, tool_name, args, preview):
        self.asked.append(tool_name)
        return False


def _perms(ws: Path, mode: str, asker: "_Asker | None") -> KitToolPermissions:
    return KitToolPermissions(
        workspace=ws, mode=mode,
        confirm_callback=asker.callback if asker else None)


def _ws(tmp_path) -> Path:
    ws = tmp_path / "project"
    ws.mkdir()
    return ws


_SHELL_WRITE = "touch made-by-shell.txt"


def test_attended_default_asks_before_a_write(tmp_path):
    ws = _ws(tmp_path)
    asker = _Asker()
    out = _doc_executor.execute(
        "write_file", {"path": str(ws / "a.txt"), "content": "x"},
        _perms(ws, "default", asker))
    assert asker.asked == ["write_file"], out
    assert not (ws / "a.txt").exists()


def test_attended_default_asks_before_a_command(tmp_path):
    ws = _ws(tmp_path)
    asker = _Asker()
    out = _doc_executor.execute(
        "bash", {"command": _SHELL_WRITE}, _perms(ws, "default", asker))
    assert asker.asked, out
    assert not (ws / "made-by-shell.txt").exists()


@pytest.mark.parametrize("command", [
    "rm seed.txt", "chmod +x seed.txt", "curl -s http://127.0.0.1:1/",
    "git push"])
def test_attended_default_asks_before_a_command_off_the_allow_list(
        tmp_path, command):
    ws = _ws(tmp_path)
    (ws / "seed.txt").write_text("s")
    asker = _Asker()
    out = str(_doc_executor.execute(
        "bash", {"command": command}, _perms(ws, "default", asker)))
    assert asker.asked, out
    assert "exit_code" not in out, "a refused command must not have run"
    assert (ws / "seed.txt").exists()


@pytest.mark.parametrize("command", [
    "rm seed.txt", "chmod +x seed.txt", "curl -s http://127.0.0.1:1/",
    "git push"])
def test_unattended_default_refuses_a_command_off_the_allow_list(
        tmp_path, command):
    ws = _ws(tmp_path)
    (ws / "seed.txt").write_text("s")
    out = str(_doc_executor.execute(
        "bash", {"command": command}, _perms(ws, "default", None)))
    assert "exit_code" not in out, out
    assert "auto-allow" in out, out
    assert (ws / "seed.txt").exists()
    assert os.access(ws / "seed.txt", os.X_OK) is False


def test_unattended_default_writes_inside_the_workspace(tmp_path):
    """The documented exception: a headless default run can produce files
    inside the workspace roots, by a file tool or a shell write (measured:
    touch, mv, a redirect and mkdir run; rm, chmod, curl and git push are
    refused)."""
    ws = _ws(tmp_path)
    _doc_executor.execute(
        "write_file", {"path": str(ws / "a.txt"), "content": "x"},
        _perms(ws, "default", None))
    _doc_executor.execute(
        "bash", {"command": _SHELL_WRITE}, _perms(ws, "default", None))
    assert (ws / "a.txt").read_text() == "x"
    assert (ws / "made-by-shell.txt").exists()


def test_unattended_default_refuses_a_write_outside_the_workspace(tmp_path):
    ws = _ws(tmp_path)
    outside = tmp_path / "elsewhere.txt"
    _doc_executor.execute(
        "write_file", {"path": str(outside), "content": "x"},
        _perms(ws, "default", None))
    assert not outside.exists()


def test_plan_refuses_a_write_and_a_command_even_attended(tmp_path):
    ws = _ws(tmp_path)
    asker = _Asker()
    perms = _perms(ws, "plan", asker)
    _doc_executor.execute(
        "write_file", {"path": str(ws / "a.txt"), "content": "x"}, perms)
    _doc_executor.execute("bash", {"command": _SHELL_WRITE}, perms)
    assert not (ws / "a.txt").exists()
    assert not (ws / "made-by-shell.txt").exists()
