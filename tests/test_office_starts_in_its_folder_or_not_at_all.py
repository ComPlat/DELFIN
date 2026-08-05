"""An office session that cannot find its folder must not start anyway.

Office is defined by the folder it works in: the permission layer locks
the role to the workspace, so the workspace IS the boundary the user was
shown. The dashboard resolved that folder like this --

    _office_p = _abs_dir(ctx.office_dir)
    if not _office_p:
        try: ctx.office_dir.mkdir(parents=True, exist_ok=True)
        except Exception: _office_p = None
    if _office_p:
        repo_dir = _office_p           # ... and everything else

-- with no else. A read-only home, a `paths.office_dir` pointing into a
device that is not mounted, a permissions problem: any of them skipped
the whole block, and the session then ran as office_agent with the
LAUNCH directory as its workspace. Locked to it, too, since the lock
follows the role. The user is told "Office works in your documents
folder" while the agent is confined to the DELFIN checkout.

A wrong folder is worse than no session, because nothing about it looks
wrong from the outside. So: try what is configured, then the documented
default, and if neither can exist, refuse and say so.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.dashboard.tab_agent import resolve_office_workspace


def test_the_configured_folder_is_used(tmp_path):
    want = tmp_path / "Dokumente"
    want.mkdir()
    assert resolve_office_workspace(want) == want.resolve()


def test_a_missing_folder_is_created(tmp_path):
    """Asking for office work is asking for the folder to exist."""
    want = tmp_path / "neu" / "Buero"
    assert resolve_office_workspace(want) == want.resolve()
    assert want.is_dir()


def test_no_configuration_falls_back_to_the_documented_default(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    got = resolve_office_workspace(None)
    assert got == (tmp_path / "office").resolve()


def test_an_unusable_folder_does_not_fall_through_to_the_launch_dir(tmp_path, monkeypatch):
    """The failure this exists for: no folder means no office session."""
    blocked = tmp_path / "kein-zugriff"
    blocked.write_text("this is a file, not a directory", encoding="utf-8")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: blocked))
    assert resolve_office_workspace(blocked / "sub") is None


def test_a_file_where_the_folder_should_be_falls_back(tmp_path, monkeypatch):
    """A misconfigured path is not fatal on its own -- the documented
    default is tried next, and only failing BOTH refuses the session."""
    occupied = tmp_path / "Buero"
    occupied.write_text("not a directory", encoding="utf-8")
    home = tmp_path / "home"
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    assert resolve_office_workspace(occupied) == (home / "office").resolve()


def test_the_result_is_absolute(tmp_path, monkeypatch):
    """A relative workspace would be resolved against the process cwd,
    which under the notebook server is not where the user thinks."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "rel").mkdir()
    got = resolve_office_workspace(Path("rel"))
    assert got is not None and got.is_absolute()
