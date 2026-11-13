"""What undo_changes restores from cannot be changed by anything the agent
runs, overlay or not.

The workspace overlay (``overlay_workspace``) keeps the originals a
command replaced in ``<root>/.delfin-trash/.undo/`` and a running
command's layer in ``<root>/.delfin-trash/.layers/``. A background
command, a command in an asking mode, or a file tool could otherwise
delete or rewrite them: ``python -c "shutil.rmtree('.delfin-trash')"``
in ``bash_background`` would take every undo of the session with it.

Contract: every bubblewrap cage binds both directories read-only where
they exist (``api_client._undo_store_binds``); Seatbelt denies writes
there; the file tools' write gate refuses them. The rest of
``.delfin-trash`` (the day folders the trash refusal tells the agent to
move a file into) stays writable. Landlock cannot carve a read-only
subtree out of a writable root without making the root refuse new files,
so a Landlock-only host leaves the store to the file-tool gate and the
command-text checks; that is said in docs/SECURITY_MODEL.md.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    w = tmp_path / "ws"
    store = w / ".delfin-trash" / ".undo" / "2026-11-13" / "cmd-1"
    store.mkdir(parents=True)
    (store / "orig.dat").write_bytes(b"original")
    (w / ".delfin-trash" / ".layers" / "9" / "u").mkdir(parents=True)
    (w / ".delfin-trash" / "2026-11-13").mkdir()
    (w / "big.dat").write_bytes(b"b" * 10)
    return Path(os.path.realpath(w))


def _perms(ws, mode="bypassPermissions"):
    p = KitToolPermissions(mode=mode, workspace=str(ws))
    p.task_session_id = f"us-{os.getpid()}"
    return p


@pytest.mark.parametrize("mode", ["bypassPermissions", "acceptEdits"])
def test_a_cage_cannot_change_the_store(ws, mode, monkeypatch):
    # The argv is run here, outside the executor that starts the session's
    # egress proxy; the network is not what this test is about.
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "open")
    argv = A._bash_isolation_argv(
        "python3 -c \"import shutil; [shutil.rmtree('.delfin-trash/' + d, "
        "ignore_errors=True) for d in ('.undo', '.layers')]\"; "
        "echo planted > .delfin-trash/.layers/9/u/x; "
        "mv big.dat .delfin-trash/2026-11-13/ && echo moved",
        ws, _perms(ws, mode))
    if argv[0] != "bwrap":
        # No bubblewrap here: the binds cannot apply; the file tools'
        # gate (below) is what holds.
        return
    joined = " ".join(argv)
    for sub in (".undo", ".layers"):
        assert f"--ro-bind-try {ws}/.delfin-trash/{sub} " in joined
    proc = subprocess.run(argv, cwd=ws, capture_output=True, text=True)
    store = ws / ".delfin-trash" / ".undo" / "2026-11-13" / "cmd-1"
    assert (store / "orig.dat").read_bytes() == b"original", proc.stderr
    assert not (ws / ".delfin-trash" / ".layers" / "9" / "u" / "x").exists()
    # The day folders stay writable: the trash refusal's remedy works.
    assert "moved" in proc.stdout, proc.stderr
    assert (ws / ".delfin-trash" / "2026-11-13" / "big.dat").exists()


@pytest.mark.parametrize("rel", [".delfin-trash/.undo/2026-11-13/cmd-1/new.dat",
                                 ".delfin-trash/.layers/9/u/x"])
@pytest.mark.parametrize("mode", ["bypassPermissions", "acceptEdits"])
def test_the_file_tools_refuse_the_store(ws, rel, mode):
    perms = _perms(ws, mode)
    out = str(_doc_executor.execute(
        "write_file", {"path": rel, "content": "tampered"}, perms))
    assert "error" in out and "undo" in out, out
    target = ws / rel
    assert not target.exists() or target.read_bytes() == b"original"


def test_a_file_tool_may_still_write_beside_the_store(ws):
    perms = _perms(ws)
    out = str(_doc_executor.execute(
        "write_file", {"path": "notes.txt", "content": "ok"}, perms))
    assert (ws / "notes.txt").read_text() == "ok", out


def test_seatbelt_denies_writes_to_the_store(ws, monkeypatch):
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "none")
    argv = A._seatbelt_argv(["/bin/true"], _perms(ws))
    prof = next(a for a in argv if "(version 1)" in a)
    for sub in (".undo", ".layers"):
        assert f'(subpath "{ws}/.delfin-trash/{sub}")' in prof, prof


def test_the_change_journal_does_not_count_the_store_as_the_commands(ws):
    """A background command's snapshot walk (shell_changes) saw another
    command's layer status file appear and journalled it as created by
    the background job; DELFIN's own store is not the command's change."""
    from delfin.agent import shell_changes as SC
    snap = SC.before([ws])
    (ws / ".delfin-trash" / ".layers" / "9" / "status").write_text("ok")
    (ws / ".delfin-trash" / ".undo" / "2026-11-13" / "cmd-1" / "n").write_text("n")
    (ws / "made.txt").write_text("m")
    names = sorted(Path(c.path).name for c in SC.after(snap))
    assert names == ["made.txt"], names
