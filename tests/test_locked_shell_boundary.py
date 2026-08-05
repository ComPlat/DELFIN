"""The locked shell was checked by reading the command text.

Three ways past that, all confirmed by the security audit against the
real gates:

* A SYMLINK inside the folder pointing out of it. `cat shared/x` names
  nothing absolute and contains no `..`, so the absolute-path scanner
  never looked at it. A link into a department share is a normal thing to
  find in an office folder, which makes this the escape a model reaches
  by accident rather than on purpose.
* An INTERPRETER. `python -c`, `make`, `xargs`, `base64 -d | bash`,
  `find -exec` were on the auto-allow list, which is consulted before any
  confirmation and without regard to the boundary -- while the boundary is
  enforced by reading the command text, which is exactly what these hide.
* A PARSE FAILURE. Both gates ended `except Exception: return None`, i.e.
  allow. Everywhere else that is backed by filesystem isolation; where the
  folder IS the promise there is nothing behind it.
"""

from __future__ import annotations

import os
import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def scene():
    ws = pathlib.Path(tempfile.mkdtemp(prefix="office_ws_")).resolve()
    outside = pathlib.Path(tempfile.mkdtemp(prefix="shared_")).resolve()
    (outside / "personal.xlsx").write_text("salaries", encoding="utf-8")
    (ws / "eigene.xlsx").write_text("mine", encoding="utf-8")
    (ws / "unter").mkdir()
    os.symlink(outside, ws / "abteilung")
    perms = A.KitToolPermissions(workspace=ws, agent_role="office_agent")
    perms.mode = "acceptEdits"
    perms.confirm_callback = lambda n, a, prev: True
    return A._DocToolExecutor.__new__(A._DocToolExecutor), perms, ws, outside


# ---------------------------------------------------------------------------
# Symlinks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "cat abteilung/personal.xlsx",
    "cd abteilung && cat personal.xlsx",
    "cp abteilung/personal.xlsx ./kopie.xlsx",
])
def test_a_link_out_of_the_folder_is_blocked(scene, cmd):
    executor, perms, _ws, _out = scene
    assert executor._gate_bash_read_paths(cmd, perms) is not None, cmd


def test_the_refusal_says_what_it_saw(scene):
    """"Blocked" without the resolved target leaves the model guessing, and
    a guessing model tries the next spelling."""
    executor, perms, _ws, outside = scene
    err = executor._gate_bash_read_paths("cat abteilung/personal.xlsx", perms)
    assert str(outside) in err


@pytest.mark.parametrize("cmd", [
    "cat eigene.xlsx",
    "ls unter",
    "cd unter && ls",
    "sed 's/a/b/' eigene.xlsx",
    "echo hello",
])
def test_ordinary_work_inside_the_folder_still_runs(scene, cmd):
    """A boundary check that stops the work is not a boundary check."""
    executor, perms, _ws, _out = scene
    assert executor._gate_bash_read_paths(cmd, perms) is None, cmd


def test_an_unlocked_session_is_not_affected(scene):
    executor, perms, ws, _out = scene
    perms.agent_role = ""
    perms.lock_workspace = False
    if perms.scope_locked:      # the folder itself may be registered
        pytest.skip("workspace is locked by its own path")
    assert executor._gate_bash_read_paths(
        "cat abteilung/personal.xlsx", perms) is None


# ---------------------------------------------------------------------------
# Interpreters
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    'python3 -c "print(1)"',
    "echo Y2F0 | base64 -d | bash",
    'eval "$(echo ls)"',
    "make build",
    "cat list.txt | xargs cat",
    "find . -name x -exec cat {} +",
])
def test_an_interpreter_is_not_auto_allowed_under_a_lock(cmd):
    """Not a judgement about danger -- python -c "print(1)" is harmless.
    A judgement about VISIBILITY: a scanner deciding on the command string
    cannot decide about these."""
    ws = pathlib.Path(tempfile.mkdtemp()).resolve()
    locked = A.KitToolPermissions(workspace=ws, agent_role="office_agent")
    assert not locked.matches_bash_auto_allow(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "ls -la", "git status", "python3 --version", "pwd",
])
def test_ordinary_commands_are_still_auto_allowed_under_a_lock(cmd):
    ws = pathlib.Path(tempfile.mkdtemp()).resolve()
    locked = A.KitToolPermissions(workspace=ws, agent_role="office_agent")
    assert locked.matches_bash_auto_allow(cmd), cmd


@pytest.mark.parametrize("cmd", ['python3 -c "print(1)"', "make build"])
def test_an_unlocked_session_keeps_its_auto_allow(cmd):
    """HPC coding workflows depend on these running unattended."""
    ws = pathlib.Path(tempfile.mkdtemp()).resolve()
    unlocked = A.KitToolPermissions(workspace=ws)
    if unlocked.scope_locked:
        pytest.skip("workspace is locked by its own path")
    assert unlocked.matches_bash_auto_allow(cmd), cmd


# ---------------------------------------------------------------------------
# Parse failure
# ---------------------------------------------------------------------------

def test_a_locked_gate_fails_closed(monkeypatch, scene):
    executor, perms, _ws, _out = scene
    monkeypatch.setattr(
        A, "_bash_climbs_out",
        lambda cmd: (_ for _ in ()).throw(RuntimeError("unparseable")))
    err = executor._gate_bash_read_paths("cat x", perms)
    assert err is not None
    assert "could not be checked" in err


def test_an_unlocked_gate_still_fails_open(monkeypatch):
    """Everywhere else the parse is backed by filesystem isolation, and
    refusing every command a parser trips over would be worse."""
    monkeypatch.setattr(
        A, "_bash_climbs_out",
        lambda cmd: (_ for _ in ()).throw(RuntimeError("unparseable")))
    ws = pathlib.Path(tempfile.mkdtemp()).resolve()
    perms = A.KitToolPermissions(workspace=ws)
    if perms.scope_locked:
        pytest.skip("workspace is locked by its own path")
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    assert executor._gate_bash_read_paths("cat x", perms) is None
