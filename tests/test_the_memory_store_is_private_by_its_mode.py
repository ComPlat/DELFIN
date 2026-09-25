"""Directory modes under ``~/.delfin``: the listing, not the contents.

Input: a directory the memory store wants. Output: that directory, with
every level from ``~/.delfin`` down at mode 0o700.

Semantics. Memory files are written 0600 already, so this is about who can
LIST them. The names are descriptive, plan names are descriptive, and the
store slug under the default ``path`` memory key is the full checkout path.
A directory another account cannot traverse cannot be listed whatever the
modes inside it, so the level that matters most is the topmost one.

Measured, and the reason a one-line fix does not work:

    mkdir(parents=True, mode=0o700)   parents get the umask, not 0o700
    mkdir(mode=0o700, exist_ok=True)  an EXISTING directory keeps its mode

The second is decisive: every installation that has run before already has
these directories, so a fix that only passes ``mode=`` would change nothing
anywhere and still look like a fix.

How the exposure arises: ``~/.delfin`` has no creator that sets its mode.
It reaches 0o700 only as a side effect of credentials.py chmod-ing its own
parent, so an installation that never wrote credentials keeps whatever the
umask gave -- 0o755 under the common 022.
"""

from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import memory_store as M


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """A home of our own, and the permissive umask that causes the defect."""
    old = os.umask(0o022)
    monkeypatch.setattr(M.Path, "home", classmethod(lambda cls: tmp_path))
    yield tmp_path
    os.umask(old)


def _mode(p):
    return stat.S_IMODE(p.stat().st_mode)


# -- the two traps, pinned so a later simplification cannot reintroduce them

def test_plain_mkdir_leaves_the_parents_open(home):
    """The control: what the store did before, and what it produced."""
    d = home / ".delfin" / "projects" / "slug" / "memory"
    d.mkdir(parents=True, mode=0o700)
    assert _mode(d) == 0o700
    assert _mode(d.parent) == 0o755, "parents take the umask, not the mode"
    assert _mode(home / ".delfin") == 0o755


def test_mkdir_does_not_retighten_a_directory_that_exists(home):
    """Why the helper chmods instead of only passing a mode."""
    d = home / ".delfin"
    d.mkdir()
    assert _mode(d) == 0o755
    d.mkdir(mode=0o700, exist_ok=True)
    assert _mode(d) == 0o755, "exist_ok makes the mode a no-op"


# -- what the helper does --------------------------------------------------

def test_a_fresh_store_is_private_at_every_level(home):
    d = M._private_dir(home / ".delfin" / "projects" / "slug" / "memory")
    assert d.is_dir()
    for level in (home / ".delfin", home / ".delfin" / "projects",
                  d.parent, d):
        assert _mode(level) == 0o700, level


def test_an_installation_that_already_ran_is_repaired(home):
    """The case that matters: these directories already exist everywhere."""
    d = home / ".delfin" / "projects" / "slug" / "memory"
    d.mkdir(parents=True)
    assert _mode(home / ".delfin") == 0o755
    M._private_dir(d)
    for level in (home / ".delfin", home / ".delfin" / "projects",
                  d.parent, d):
        assert _mode(level) == 0o700, level


def test_a_home_that_does_not_exist_yet_is_created_too(home):
    """The first run on a machine. Found by a test, not by a probe: a
    manual check used a home tempfile had already created, so the root's
    own parent always existed and the hole stayed invisible."""
    fresh = home / "not-yet"
    d = fresh / ".delfin" / "projects" / "slug" / "memory"
    M.Path.home = classmethod(lambda cls: fresh)
    try:
        M._private_dir(d)
    finally:
        del M.Path.home
    assert d.is_dir(), "a store under a home that did not exist was not made"
    assert _mode(fresh / ".delfin") == 0o700


def test_the_root_itself_can_be_the_target(home):
    M._private_dir(home / ".delfin")
    assert _mode(home / ".delfin") == 0o700


def test_the_home_directory_is_not_touched(home):
    """It is not this tool's directory, so it is not this tool's decision."""
    before = _mode(home)
    M._private_dir(home / ".delfin" / "projects" / "slug" / "memory")
    assert _mode(home) == before


def test_a_path_outside_is_created_but_not_tightened(home):
    """A caller naming somewhere else gets a directory, not a policy."""
    outside = home / "elsewhere" / "dir"
    M._private_dir(outside)
    assert outside.is_dir()
    assert _mode(outside) == 0o755


def test_it_never_raises_when_the_mode_cannot_be_set(home, monkeypatch):
    def _boom(*a, **k):
        raise OSError("read-only")
    monkeypatch.setattr(M.os, "chmod", _boom)
    d = M._private_dir(home / ".delfin" / "projects" / "slug" / "memory")
    assert d.is_dir(), "a store that cannot be tightened is still a store"


def test_it_never_raises_when_the_directory_cannot_be_made(home, monkeypatch):
    def _boom(*a, **k):
        raise OSError("no space")
    monkeypatch.setattr(M.Path, "mkdir", _boom)
    M._private_dir(home / ".delfin" / "projects" / "slug" / "memory")


# -- the store actually uses it --------------------------------------------

def test_the_memory_directory_comes_back_private(home, tmp_path):
    d = M._delfin_memory_dir(tmp_path)
    d.mkdir(parents=True, exist_ok=True)
    M._private_dir(d)
    assert _mode(d) == 0o700
    assert _mode(home / ".delfin") == 0o700


def test_no_site_in_the_store_creates_a_directory_the_plain_way():
    """A new call site that forgets the helper reopens the listing."""
    import inspect
    import re

    # Everything except the helper's own body -- its fallback for a path
    # outside ~/.delfin is the one legitimate plain creation in the module.
    src = inspect.getsource(M).replace(inspect.getsource(M._private_dir), "")
    plain = re.findall(r"\.mkdir\(parents=True,\s*exist_ok=True\)", src)
    assert not plain, (
        f"{len(plain)} directory creations bypass _private_dir")
