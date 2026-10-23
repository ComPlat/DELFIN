"""A file git reads as configuration is git configuration on every route.

``git config <program-key> <value>`` is asked in every mode and refused
with nobody to ask (``_git_program_config_write``). The same key reaches
git when a configuration FILE is written directly. The home directory is
never a workspace, so the files a session can reach are those the global
configuration places inside one: the file named by ``GIT_CONFIG_GLOBAL``,
``$XDG_CONFIG_HOME/git/config`` and ``attributes`` (attributes select
filter and diff drivers), a ``~/.gitconfig`` that is a symlink into a
dotfiles checkout, ``core.attributesFile``, and every file an
``include.path`` or ``includeIf.<cond>.path`` names, followed through
includes -- in the global file and in the repository's own config.

Driven through the real executor (``_DocToolExecutor.execute``) in every
mode, attended (a callback that answers no) and headless (no callback),
under a private HOME. Controls: an ordinary file, ``git config user.name``
and reading a configuration file are not asked about.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent import api_client as A
from delfin.agent import job_monitor as jm

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the repository")

MODES = ["default", "acceptEdits", "bypassPermissions"]


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    h = tmp_path / "home"
    work = h / "work"
    (work / "inc").mkdir(parents=True)
    (h / ".gitconfig").write_text(
        "[user]\n\tname = t\n\temail = t@t\n"
        "[include]\n\tpath = ~/work/inc/extra.cfg\n"
        "[includeIf \"gitdir:~/\"]\n\tpath = work/cond.cfg\n"
        "[core]\n\tattributesFile = ~/work/my.attributes ; a comment\n")
    (work / "inc" / "extra.cfg").write_text(
        "[include]\n\tpath = \"nested.cfg\"\n")
    (work / "inc" / "nested.cfg").write_text("[color]\n\tui = auto\n")
    for var in ("XDG_CONFIG_HOME", "GIT_CONFIG_GLOBAL", "GIT_CONFIG_SYSTEM",
                "GIT_CONFIG"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    return h


@pytest.fixture
def work(home):
    return home / "work"


@pytest.fixture
def repo(home, work, tmp_path):
    r = work / "repo"
    r.mkdir()
    env = child_env(tmp_path)
    env["HOME"] = str(home)
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    subprocess.run(["git", "init", "-q", str(r)], check=True,
                   capture_output=True, env=env)
    with open(r / ".git" / "config", "a") as fh:
        fh.write("[include]\n\tpath = ../shared.cfg\n")
    (r / "plain.txt").write_text("a\n")
    A._GIT_ROOTS_CACHE.clear()
    return r


def _perms(ws, mode, answer=None):
    perms = A.KitToolPermissions(workspace=ws, mode=mode)
    perms.asked = []
    if answer is not None:
        def _cb(name, args, preview=""):
            perms.asked.append(str(preview))
            if answer == "all but git":
                return "GIT CONFIGURATION" not in str(preview)
            return answer
        perms.confirm_callback = _cb
    return perms


def _run(perms, tool, args):
    return A._DocToolExecutor().execute(tool, args, perms)


def _read_baseline(perms, ws, read):
    # Writes to an existing file need a read baseline first; without it the
    # refusal would be the baseline's, not the one under test.
    for rel in read:
        if (ws / rel).exists():
            _run(perms, "read_file", {"path": rel})


def _assert_git_config_gate(ws, tool, args, mode, read=()):
    # Attended: every other dialog (the command itself in default mode) is
    # approved; the git-configuration one is shown and refused, and that
    # alone must keep the write out.
    attended = _perms(ws, mode, answer="all but git")
    _read_baseline(attended, ws, read)
    out = _run(attended, tool, args)
    assert any("GIT CONFIGURATION" in p for p in attended.asked), (
        tool, args, mode, attended.asked, out)
    assert "exit_code" not in out, (tool, args, mode, out)
    headless = _perms(ws, mode)
    _read_baseline(headless, ws, read)
    out = _run(headless, tool, args)
    low = out.lower()
    assert "exit_code" not in out and "error" in low, (tool, args, mode, out)
    if tool != "bash" or mode == "bypassPermissions":
        # A headless shell write in default/acceptEdits is refused earlier,
        # as any unlisted command is; elsewhere this gate is the refusal.
        assert "git" in low and "configuration" in low and (
            "refused" in low or "blocked" in low), (tool, args, mode, out)


# Relative to the workspace ``work`` (``repo/`` -- to the repository).
NAMED_FILES = [
    "inc/extra.cfg",                 # include.path in ~/.gitconfig
    "inc/nested.cfg",                # an include inside an include
    "cond.cfg",                      # includeIf, relative to its file
    "my.attributes",                 # core.attributesFile
    "repo/shared.cfg",               # include.path in the repository
]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("rel", NAMED_FILES)
def test_write_file_is_asked_in_every_mode(work, repo, rel, mode):
    ws = repo if rel.startswith("repo/") else work
    rel = rel[len("repo/"):] if rel.startswith("repo/") else rel
    _assert_git_config_gate(
        ws, "write_file",
        {"path": rel, "content": "[core]\n\tpager = ./p.sh\n"}, mode)
    target = ws / rel
    assert not target.exists() or "p.sh" not in target.read_text()


NESTED = ("inc/nested.cfg",)


@pytest.mark.parametrize("mode", MODES)
def test_edit_and_patch_routes_are_asked(work, mode):
    _assert_git_config_gate(
        work, "edit_file",
        {"path": "inc/nested.cfg", "old_string": "ui = auto",
         "new_string": "ui = auto\n[core]\n\tpager = ./p.sh"}, mode, read=NESTED)
    _assert_git_config_gate(
        work, "multi_edit",
        {"path": "inc/nested.cfg",
         "edits": [{"old_string": "ui = auto", "new_string": "ui = never"}]},
        mode, read=NESTED)
    patch = ("--- a/inc/nested.cfg\n+++ b/inc/nested.cfg\n"
             "@@ -1,2 +1,3 @@\n [color]\n \tui = auto\n+[core]\n")
    _assert_git_config_gate(work, "apply_patch", {"diff": patch}, mode,
                            read=NESTED)
    assert (work / "inc" / "nested.cfg").read_text() == \
        "[color]\n\tui = auto\n"


@pytest.mark.parametrize("mode", MODES)
def test_the_global_files_follow_the_environment(
        home, work, monkeypatch, mode):
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(work / "global.cfg"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(work / "xdg"))
    for rel in ("global.cfg", "xdg/git/config", "xdg/git/attributes"):
        _assert_git_config_gate(
            work, "write_file", {"path": rel, "content": "x\n"}, mode)


@pytest.mark.parametrize("mode", MODES)
def test_a_global_file_linked_into_the_workspace_is_asked(home, work, mode):
    (work / "dotfiles").mkdir()
    (work / "dotfiles" / "gitconfig").write_text(
        (home / ".gitconfig").read_text())
    (home / ".gitconfig").unlink()
    os.symlink(work / "dotfiles" / "gitconfig", home / ".gitconfig")
    _assert_git_config_gate(
        work, "write_file", {"path": "dotfiles/gitconfig", "content": "x\n"},
        mode)
    # Its includes are still followed through the link.
    _assert_git_config_gate(
        work, "write_file", {"path": "cond.cfg", "content": "x\n"}, mode)


SHELL_WRITES = [
    "printf '[core]\\n\\tpager = ./p.sh\\n' >> inc/extra.cfg",
    "echo '* filter=x' | tee -a my.attributes",
    "cp plain.cfg cond.cfg",
    "mv plain.cfg inc/extra.cfg",
    "sed -i 's/ui = auto/ui = never/' inc/nested.cfg",
    "python3 -c \"open('cond.cfg','a').write('x')\"",
]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", SHELL_WRITES)
def test_the_shell_routes_are_asked(work, cmd, mode):
    (work / "plain.cfg").write_text("[color]\n\tui = auto\n")
    _assert_git_config_gate(work, "bash", {"command": cmd}, mode,
                            read=("inc/extra.cfg", "inc/nested.cfg"))
    assert (work / "inc" / "extra.cfg").read_text() == \
        "[include]\n\tpath = \"nested.cfg\"\n"
    assert (work / "inc" / "nested.cfg").read_text() == \
        "[color]\n\tui = auto\n"
    assert not (work / "cond.cfg").exists()
    assert not (work / "my.attributes").exists()


GIT_CONFIG_FILE_COMMANDS = [
    "git config --edit",
    "git config -e",
    "git config --global --edit",
    "git config -f other.cfg --edit",
    "git config edit",
    "git config -f other.cfg core.pager ./p.sh",
    "git config --file=other.cfg alias.x '!id'",
    "git config submodule.lib.update '!sh x'",
]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", GIT_CONFIG_FILE_COMMANDS)
def test_git_config_edit_and_file_writes_are_asked(repo, cmd, mode):
    _assert_git_config_gate(repo, "bash", {"command": cmd}, mode)


# --- controls ---------------------------------------------------------------

@pytest.mark.parametrize("mode", ["acceptEdits", "bypassPermissions"])
def test_an_ordinary_file_is_not_asked_about(work, repo, mode):
    for ws, rel in ((work, "notes.txt"), (repo, "new.txt"),
                    (repo, "sub/config"), (work, "inc/other.cfg")):
        perms = _perms(ws, mode, answer=True)
        out = _run(perms, "write_file", {"path": rel, "content": "x\n"})
        assert perms.asked == [], (rel, perms.asked)
        assert (ws / rel).read_text() == "x\n", (rel, out)


def test_a_harmless_git_config_write_is_not_asked_and_the_mount_decides(repo):
    """The gate does not ask about a harmless key. Whether the write lands
    is the sandbox's: bubblewrap and Seatbelt keep .git/config read-only
    (the config names programs git runs outside the sandbox later), so
    there the command fails with git's own message and the file is
    unchanged; on a host without such a layer the write lands."""
    perms = _perms(repo, "bypassPermissions")
    cfg = repo / ".git" / "config"
    before = cfg.read_text()
    out = json.loads(_run(perms, "bash", {"command": "git config user.name x"}))
    assert perms.asked == []
    if out.get("exit_code") == 0:
        assert "name = x" in cfg.read_text()
    else:
        assert "could not write config file" in out.get("stderr", ""), out
        assert cfg.read_text() == before


def test_reading_a_configuration_file_is_not_asked(work):
    perms = _perms(work, "default", answer=False)
    out = _run(perms, "read_file", {"path": "inc/extra.cfg"})
    assert perms.asked == [] and "nested.cfg" in out, out
    out = _run(perms, "bash", {"command": "cat inc/extra.cfg"})
    assert perms.asked == [] and "nested.cfg" in out, out


def test_the_repository_config_stays_denied(repo):
    for rel in (".git/config", ".git/info/attributes"):
        perms = _perms(repo, "bypassPermissions", answer=True)
        out = _run(perms, "write_file", {"path": rel, "content": "x\n"})
        assert "error" in json.loads(out), (rel, out)
    assert "[include]" in (repo / ".git" / "config").read_text()


def test_the_module_under_test_is_this_checkout():
    # A worktree isolates files, not imports.
    here = Path(__file__).resolve().parents[1]
    assert Path(A.__file__).resolve().is_relative_to(here)
