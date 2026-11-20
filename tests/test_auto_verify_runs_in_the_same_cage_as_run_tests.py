"""Auto-verify runs the workspace's test code in the cage, or not at all.

A turn that edited a ``.py`` file ends with a verification run
(``agent.auto_verify``, default ``smart``). That run executes pytest in
the workspace, and pytest executes the workspace's own ``conftest.py``
before any test -- so the workspace decides what code runs, and nobody
is asked, because auto-verify is not a tool call.

The explicit ``run_tests`` tool wraps its pytest in the filesystem cage
(``_test_run_argv`` -> ``_bash_isolation_argv``). The verification run did
not: it was ``subprocess.run(command, shell=True)``.

Measured against that code, by calling ``_run_auto_verify`` on a workspace
built exactly as ``_workspace`` below builds one: the ``conftest.py``
wrote its file OUTSIDE the workspace and the write landed. Two ways in,
neither of them exotic -- an instruction smuggled into fetched text asks
for a ``conftest.py`` (a file write is not even asked about under
``acceptEdits``), or a cloned repository simply brings one along. The
tests here cannot be run against that code as they stand, because they
hand over the ``perms`` it has no parameter for; what they hold is the
behaviour that must not come back.

Input: a workspace holding a ``conftest.py`` that writes outside it, and
an edited ``mod.py``. Output: whether the outside write landed, and what
the verification status says.

Semantics pinned here: the outside write never lands -- either the cage
stops it, or, on a host that can cage nothing, the command is not run and
the status says so. The second case is the one that keeps unattended runs
working on such a host instead of refusing to verify at all; the syntax
check still runs there.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions

CONFTEST = """\
from pathlib import Path

Path({marker!r}).write_text("escaped")
"""


def _workspace(tmp_path: Path) -> tuple[Path, Path]:
    """A workspace whose ``conftest.py`` writes to a path outside it."""
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    marker = outside / "escaped.txt"
    (ws / "conftest.py").write_text(CONFTEST.format(marker=str(marker)))
    (ws / "test_the_workspace.py").write_text("def test_ok():\n    assert True\n")
    (ws / "mod.py").write_text("VALUE = 1\n")
    return ws, marker


def _perms(ws: Path, mode: str) -> KitToolPermissions:
    return KitToolPermissions(workspace=ws, mode=mode, confirm_callback=None)


@pytest.mark.parametrize("mode", ["bypassPermissions", "acceptEdits", "default"])
def test_the_workspaces_own_code_does_not_write_outside_it(tmp_path, mode):
    """The whole point, in every permission mode.

    ``acceptEdits`` and ``default`` matter as much as the unattended mode:
    the person approves command TEXT, and this command's text is one the
    agent never shows them.
    """
    ws, marker = _workspace(tmp_path)
    status: dict = {}
    A._run_auto_verify([str(ws / "mod.py")], "smart", "", str(ws),
                       status=status, perms=_perms(ws, mode))
    assert not marker.exists(), (
        "a conftest.py in the workspace wrote outside the workspace; "
        f"verification status: {status}")


def test_a_command_the_person_configured_is_caged_too(tmp_path):
    """``auto_verify = "command"`` runs the person's own command -- in the
    workspace, where the ``conftest.py`` is. Their choice of command is not
    a statement about the workspace's contents."""
    ws, marker = _workspace(tmp_path)
    status: dict = {}
    A._run_auto_verify([str(ws / "mod.py")], "command",
                       "python -m pytest -q", str(ws),
                       status=status, perms=_perms(ws, "bypassPermissions"))
    assert not marker.exists(), status


def test_without_permissions_nothing_is_run(tmp_path):
    """No permission bundle means no cage to run in. Verification then
    reports a skip with its reason -- it must not fall back to a plain
    shell, and it must not read as verified."""
    ws, marker = _workspace(tmp_path)
    status: dict = {}
    out = A._run_auto_verify([str(ws / "mod.py")], "smart", "", str(ws),
                             status=status, perms=None)
    assert not marker.exists()
    assert out == ""                      # a skip is not a problem to fix
    assert status.get("skipped") is True
    assert "isolation" in str(status.get("reason", "")).lower()


def test_verification_still_reports_a_real_failure(tmp_path):
    """The cage must not cost the function it protects: a test that fails
    is still reported as a problem to fix.

    On a host that can cage nothing the run is skipped instead, and the
    status says which -- that is the only other outcome allowed here.
    """
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "mod.py").write_text("VALUE = 1\n")
    (ws / "test_mod.py").write_text(
        "import mod\n\n\ndef test_value():\n    assert mod.VALUE == 2\n")
    status: dict = {}
    out = A._run_auto_verify([str(ws / "mod.py")], "smart", "", str(ws),
                             status=status, perms=_perms(ws, "bypassPermissions"))
    if status.get("skipped"):
        assert "isolation" in str(status.get("reason", "")).lower(), status
    else:
        assert "failed" in out, (out, status)


def test_the_syntax_check_still_runs_without_a_cage(tmp_path):
    """A host that can cage nothing still gets the check that needs no
    cage: ``py_compile`` does not import the file, so no workspace code
    runs. Verification degrades to it instead of going silent."""
    ws = tmp_path / "ws"
    ws.mkdir()
    broken = ws / "mod.py"
    broken.write_text("def f(:\n    pass\n")
    status: dict = {}
    out = A._run_auto_verify([str(broken)], "smart", "", str(ws),
                             status=status, perms=None)
    assert out, "a syntax error must be reported with or without a cage"
