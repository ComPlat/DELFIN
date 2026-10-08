"""An auto-allowed command could write outside every workspace root.

The auto-allow list contains commands whose job is to write, and their
containment was expressed inside the patterns themselves, as a blacklist:
``mkdir\\b`` accepted any path at all, ``touch\\s+(?!/)`` accepted every path
that merely did not begin with a slash -- ``~/.local/bin/pytest`` among
them -- and ``cp``/``mv`` refused /etc, /usr, /bin, /lib and /var and
nothing else. None of them looked at a shell redirect, so
``echo ... > ~/.local/bin/pytest`` was auto-allowed on the strength of
``echo``.

That is the route a field report describes as the agent putting a test
runner in the home directory when the environment has none. It mattered
because the filesystem isolation layer is what would otherwise have
stopped the write, and that layer is not available on every host -- which
is the one condition every protection here has to hold without.

Containment is now asked once, centrally, against the roots. A blacklist
is the wrong shape for the question: there is no end of places that are
not the workspace.

The system temp directory is not an exception. Granting one was the first
version of this change and this suite refused it: a workspace under /tmp
made ``chmod 600 ../../outside/file`` land in temp, so a metadata write
that had always been refused became free. The isolation argument for the
exception only holds where that layer runs, and on a shared machine the
real /tmp carries other sessions' lock files.

Universal by construction: the "outside" path is built from
``Path.home()``, never a literal, and nothing is executed -- only the
gate's decision is read.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return KitToolPermissions(workspace=ws, mode="default")


def _home_outside() -> str:
    """A path under this user's home that no workspace root can cover."""
    return str(Path.home() / ".local" / "bin")


# ---------------------------------------------------------------------------
# Routes that must go to the confirm gate
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd_tpl", [
    "python3 -m venv {out}/v",
    "python3 -m venv ~/delfin-venv",
    "virtualenv {out}/v",
    "uv venv {out}/v",
    "python3 -m ensurepip --user",
    "echo 'exec python3 -m pytest \"$@\"' > {out}/pytest",
    "echo x >> {out}/pytest",
    "cat /etc/hostname > {out}/pytest",
    "cp /etc/hostname {out}/pytest",
    "mv x {out}/pytest",
    "install -m 755 x {out}/pytest",
    "mkdir -p {out}",
    "touch {out}/pytest",
    "touch ~/.local/bin/pytest",
    "tee {out}/pytest",
    "ln -s /bin/true {out}/pytest",
    "chmod 777 {out}",
])
def test_a_write_outside_the_roots_is_not_auto_allowed(perms, cmd_tpl):
    cmd = cmd_tpl.format(out=_home_outside())
    assert not perms.matches_bash_auto_allow(cmd), cmd


def test_a_relative_escape_is_outside_too(perms):
    """`..` resolves before the roots are consulted."""
    assert not perms.matches_bash_auto_allow("touch ../escaped")
    assert perms.matches_bash_auto_allow("touch inside.txt")


# ---------------------------------------------------------------------------
# What must stay free — the half that keeps this from being a blanket refusal
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "python3 -m venv .venv",
    "mkdir -p build/out",
    "touch notes.txt",
    "echo hi > out.txt",
    "echo hi >> out.txt",
    "cp a.txt b.txt",
    "mv a.txt b.txt",
    "python3 -m pytest -q",
    "ls -la",
    "git status",
])
def test_work_inside_the_workspace_stays_auto_allowed(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd), cmd


def test_the_system_temp_dir_is_not_an_exception(perms):
    """It is outside the roots, so it is asked about like anywhere else.

    This is the assertion that keeps the earlier, weaker version from
    coming back: with a temp exception, any workspace made under /tmp lost
    its containment for every ``..`` that climbed out of it.
    """
    import tempfile
    out = Path(tempfile.gettempdir()) / "delfin-scratch-x"
    assert not perms.matches_bash_auto_allow(f"mkdir -p {out}")


# ---------------------------------------------------------------------------
# Reading a redirect is not the same as seeing a `>`
# ---------------------------------------------------------------------------

def test_a_quoted_angle_bracket_is_not_a_redirect(perms):
    """`grep '>' .` writes nothing; treating it as a write would charge a
    dialog for a read."""
    assert perms.matches_bash_auto_allow("grep -rn '>' .")
    assert perms._redirect_targets("grep -rn '>' .") == []


def test_a_descriptor_redirect_names_no_file(perms):
    """`2>&1` renames a stream. It must not be read as a path called `&1`."""
    assert perms._redirect_targets("python3 -m pytest -q 2>&1") == []
    assert perms.matches_bash_auto_allow("python3 -m pytest -q 2>&1")


def test_a_dev_null_sink_is_not_the_filesystem(perms):
    assert perms._redirect_targets("cat a > /dev/null") == []
    assert perms.matches_bash_auto_allow("cat a.txt > /dev/null")


def test_a_redirect_is_found_with_and_without_a_space(perms):
    assert perms._redirect_targets("echo x >out.txt") == ["out.txt"]
    assert perms._redirect_targets("echo x > out.txt") == ["out.txt"]
    assert perms._redirect_targets("echo x 2> err.log") == ["err.log"]


def test_an_unparsable_command_goes_to_the_gate(perms):
    """An unbalanced quote cannot be read, so it is not auto-allowed."""
    assert perms._changes_outside_the_workspace("echo 'unterminated")


# ---------------------------------------------------------------------------
# A granted directory is a root, so writing there is work, not an escape
# ---------------------------------------------------------------------------

def test_a_granted_extra_dir_is_inside(tmp_path):
    ws, extra = tmp_path / "ws", tmp_path / "granted"
    ws.mkdir()
    extra.mkdir()
    p = KitToolPermissions(workspace=ws, mode="default",
                           extra_workspace_dirs=(extra,))
    assert p.matches_bash_auto_allow(f"mkdir -p {extra}/sub")
    assert not p.matches_bash_auto_allow(f"mkdir -p {_home_outside()}/sub")
