"""Office mode works in one folder and nowhere else.

The ordinary sandbox is a boundary with a door: a read outside it is
offered to the user for confirmation, and a directory can be granted
while the session runs. That is right when the agent works on a code
base whose dependencies live elsewhere. It is wrong for a folder of
administrative documents, where the answer to "may I read the rest of
the machine" is decided in advance.

A locked scope has no door. Outside is refused, not asked about; the
workspace cannot be widened; and the refusal holds in every permission
mode, including bypassPermissions, and against a confirm callback that
answers yes.

Filesystem isolation (bwrap) is the real containment and is forced on
for these sessions where it works. The path rules tested here are what
stands in for it where it does not.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from delfin.agent.api_client import (
    _SCOPE_LOCKED_ROLES,
    _DocToolExecutor,
    KitToolPermissions,
    _bash_climbs_out,
    _bash_paths_outside,
)


@pytest.fixture
def scene():
    """An OFFICE folder with a sibling the agent must never reach."""
    root = Path(tempfile.mkdtemp())
    office = root / "OFFICE"
    office.mkdir()
    (office / "liste.csv").write_text("A,B\n1,2\n", encoding="utf-8")
    outside = root / "privat"
    outside.mkdir()
    (outside / "geheim.txt").write_text("SECRET\n", encoding="utf-8")
    return office, outside


def _perms(office: Path, mode: str = "bypassPermissions") -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(office))
    perms.agent_role = "office_agent"
    perms.mode = mode
    # A callback that says yes to everything: the lock must not depend on
    # the user noticing a dialog.
    perms.confirm_callback = lambda *a: True
    return perms


# ---------------------------------------------------------------------------
# The lock follows the role
# ---------------------------------------------------------------------------

def test_the_office_role_is_locked_without_anyone_setting_a_flag(scene):
    office, _ = scene
    assert "office_agent" in _SCOPE_LOCKED_ROLES
    assert _perms(office).scope_locked is True


def test_other_roles_are_not_locked(scene):
    office, _ = scene
    perms = KitToolPermissions(workspace=str(office))
    perms.agent_role = "solo_agent"
    assert perms.scope_locked is False
    assert KitToolPermissions(workspace=str(office)).scope_locked is False


def test_the_lock_can_also_be_asked_for_directly(scene):
    office, _ = scene
    perms = KitToolPermissions(workspace=str(office))
    perms.lock_workspace = True
    assert perms.scope_locked is True


def test_a_locked_session_reaches_nothing_but_its_workspace(scene):
    office, _ = scene
    perms = _perms(office)
    perms.read_only_workspace_dirs = (Path("/tmp"),)
    perms.confirm_write_dirs = (Path("/var"),)
    perms.extra_workspace_dirs = (Path("/usr"),)
    assert perms.all_readable_roots() == (office,)
    assert perms.all_workspace_roots() == (office,)


def test_the_workspace_cannot_be_widened_at_runtime(scene):
    office, outside = scene
    with pytest.raises(ValueError) as exc:
        _perms(office).add_extra_dir(outside)
    assert "locked" in str(exc.value)


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------

def test_reading_inside_works(scene):
    office, _ = scene
    out = _DocToolExecutor().execute(
        "read_file", {"path": "liste.csv"}, _perms(office))
    assert "A,B" in out


def test_reading_outside_is_refused_not_offered(scene):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute(
        "read_file", {"path": str(outside / "geheim.txt")}, _perms(office)))
    assert "outside" in out["error"]
    assert "no confirmation can grant it" in out["error"]
    assert "SECRET" not in out["error"]


def test_a_relative_escape_is_refused(scene):
    office, _ = scene
    out = json.loads(_DocToolExecutor().execute(
        "read_file", {"path": "../privat/geheim.txt"}, _perms(office)))
    assert "error" in out


def test_the_document_reader_is_bound_by_the_same_rule(scene):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute(
        "read_document", {"path": str(outside / "geheim.txt")},
        _perms(office)))
    assert "outside" in out["error"]


def test_the_table_comparison_is_bound_by_the_same_rule(scene):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute("compare_tables", {
        "left": str(outside / "a.csv"), "right": "liste.csv", "key": "A",
    }, _perms(office)))
    assert "outside" in out["error"]


# ---------------------------------------------------------------------------
# Writes
# ---------------------------------------------------------------------------

def test_writing_outside_is_refused(scene):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute(
        "write_file", {"path": str(outside / "neu.txt"), "content": "x"},
        _perms(office)))
    assert "error" in out
    assert not (outside / "neu.txt").exists()


def test_writing_inside_works(scene):
    office, _ = scene
    _DocToolExecutor().execute(
        "write_file", {"path": "bericht.txt", "content": "ok"},
        _perms(office))
    assert (office / "bericht.txt").read_text() == "ok"


# ---------------------------------------------------------------------------
# The shell
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("command", [
    "cat {outside}/geheim.txt",
    "ls -la {outside}",
    "cp {outside}/geheim.txt .",
    "python3 -c \"print(open('{outside}/geheim.txt').read())\"",
])
def test_shell_commands_naming_an_outside_path_are_blocked(scene, command):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute(
        "bash", {"command": command.format(outside=outside)}, _perms(office)))
    assert "outside" in out["error"]
    assert "SECRET" not in json.dumps(out)


@pytest.mark.parametrize("command", [
    "cd .. && cat privat/geheim.txt",
    "cat ../privat/geheim.txt",
    "cp ../privat/geheim.txt .",
])
def test_a_relative_climb_out_is_blocked(scene, command):
    """No absolute path appears in these, so the path scanner never sees
    them — this is what the '..' rule is for."""
    office, _ = scene
    out = json.loads(_DocToolExecutor().execute(
        "bash", {"command": command}, _perms(office)))
    assert ".." in out["error"]


@pytest.mark.parametrize("command", [
    "ls -la",
    "for i in {1..3}; do echo $i; done",
    "python3 -c 'print(1.5)'",
    "/usr/bin/env python3 -c 'print(2)'",
    "grep -r Betrag .",
])
def test_ordinary_work_inside_the_folder_is_untouched(scene, command):
    office, _ = scene
    out = _DocToolExecutor().execute(
        "bash", {"command": command}, _perms(office))
    assert not out.lstrip().startswith('{"error"'), command


def test_a_cwd_outside_the_folder_is_refused(scene):
    office, outside = scene
    out = json.loads(_DocToolExecutor().execute(
        "bash", {"command": "ls", "cwd": str(outside)}, _perms(office)))
    assert "error" in out


# ---------------------------------------------------------------------------
# The pure helpers
# ---------------------------------------------------------------------------

def test_outside_path_scanner_exempts_system_directories(tmp_path):
    ws = tmp_path / "OFFICE"
    ws.mkdir()
    assert _bash_paths_outside("/usr/bin/env python3 x.py", ws) == []
    assert _bash_paths_outside("ls /bin", ws) == []
    assert _bash_paths_outside(f"ls {ws}", ws) == []
    assert _bash_paths_outside(f"ls {ws}/unter", ws) == []


def test_outside_path_scanner_finds_a_real_path_outside(tmp_path):
    ws = tmp_path / "OFFICE"
    ws.mkdir()
    other = tmp_path / "privat"
    other.mkdir()
    (other / "data.csv").write_text("x", encoding="utf-8")
    assert _bash_paths_outside(f"cat {other}/data.csv", ws) == [
        f"{other}/data.csv"]
    # …including inside a quoted string, which shell tokenisation would
    # hand over as one token that does not begin with a slash.
    assert _bash_paths_outside(
        f"python3 -c \"print(open('{other}/data.csv').read())\"", ws) == [
        f"{other}/data.csv"]


def test_slash_shaped_arguments_are_not_mistaken_for_paths(tmp_path):
    """A candidate counts only when it, or its parent, exists — otherwise
    `sed 's/a/b/'` would read as a path and block ordinary work."""
    ws = tmp_path / "OFFICE"
    ws.mkdir()
    assert _bash_paths_outside("sed 's/alt/neu/g' liste.csv", ws) == []
    assert _bash_paths_outside("echo a/b/c", ws) == []
    assert _bash_paths_outside("curl https://example.com/x", ws) == []


@pytest.mark.parametrize("command,climbs", [
    ("cd ..", True),
    ("cat ../x", True),
    ("cp '../x' .", True),
    ("for i in {1..10}; do :; done", False),
    ("echo 1.5", False),
    ("ls -la", False),
    ("python -c 'x = 1..2'", False),
])
def test_parent_segment_detection(command, climbs):
    assert _bash_climbs_out(command) is climbs


# ---------------------------------------------------------------------------
# Isolation
# ---------------------------------------------------------------------------

def test_filesystem_isolation_is_forced_for_a_locked_session(scene, monkeypatch):
    """Where bwrap works, take the real containment rather than trusting
    a parser to catch every shape a shell command can take."""
    import delfin.agent.api_client as A

    office, _ = scene
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    argv = A._bash_isolation_argv("ls", office, _perms(office, mode="default"))
    assert argv[0] == "bwrap"


def test_an_unlocked_interactive_session_keeps_plain_bash(scene, monkeypatch):
    import delfin.agent.api_client as A

    office, _ = scene
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    perms = KitToolPermissions(workspace=str(office))
    perms.mode = "default"
    argv = A._bash_isolation_argv("ls", office, perms)
    assert argv[0] != "bwrap"


# ---------------------------------------------------------------------------
# Delegation must not be the way out
# ---------------------------------------------------------------------------

def test_a_subagent_inherits_the_lock(scene):
    from delfin.agent.subagents import _derive_perms

    office, _ = scene
    child = _derive_perms(_perms(office), "acceptEdits")
    assert child.scope_locked is True
    assert child.all_readable_roots() == (office,)


def test_a_subagent_cannot_be_moved_out_of_the_folder(scene, tmp_path):
    """Worktree isolation relocates a child's workspace. Under a lock that
    would contain the child to one folder — just not the one the session
    was confined to."""
    from delfin.agent.subagents import _derive_perms

    office, outside = scene
    child = _derive_perms(_perms(office), "acceptEdits", workspace=outside)
    assert child.workspace == office


def test_an_unlocked_subagent_can_still_be_isolated(scene, outside=None):
    from delfin.agent.subagents import _derive_perms

    office, other = scene
    parent = KitToolPermissions(workspace=str(office))
    parent.agent_role = "solo_agent"
    child = _derive_perms(parent, "acceptEdits", workspace=other)
    assert child.workspace == other
