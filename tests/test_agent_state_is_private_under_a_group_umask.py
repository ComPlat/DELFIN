"""Under umask 002 the agent's state tree is still owner-only.

Bug hunt 11.11.2026, #7: on a fresh home under the cluster's umask 002,
``~/.delfin`` and ``projects/<slug>/memory/episodes`` came out 0775. The
account's group is shared, so a group member could drop a ``*.md`` into
``episodes/`` and ``recall_episodes`` would put it into the next prompt.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest


@pytest.fixture()
def group_umask():
    old = os.umask(0o002)
    try:
        yield
    finally:
        os.umask(old)


@pytest.fixture()
def fresh_home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    return home


def _mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def test_the_umask_really_is_group_writable(group_umask, tmp_path):
    # The control: without the fix a plain mkdir lands at 0775 here.
    d = tmp_path / "plain"
    d.mkdir(parents=True)
    assert _mode(d) == 0o775


def test_saving_an_episode_makes_every_level_owner_only(
        group_umask, fresh_home, tmp_path):
    from delfin.agent import episodes
    repo = tmp_path / "repo"
    repo.mkdir()
    path = episodes.save_episode(
        "abcdef12", repo_root=repo, goal="g", outcome="o",
        decisions=[], open_items=[])
    root = fresh_home / ".delfin"
    assert path.is_file()
    assert root in path.parents
    for level in [root, *[p for p in path.parents if root in p.parents]]:
        assert _mode(level) == 0o700, (level, oct(_mode(level)))
    assert _mode(path) & 0o077 == 0


def test_an_existing_group_writable_episodes_dir_is_repaired(
        group_umask, fresh_home, tmp_path):
    from delfin.agent import episodes
    repo = tmp_path / "repo"
    repo.mkdir()
    d = episodes._episodes_dir(repo)
    d.mkdir(parents=True)
    os.chmod(fresh_home / ".delfin", 0o775)
    os.chmod(d, 0o775)
    episodes.save_episode("abcdef12", repo_root=repo, goal="g",
                          outcome="o", decisions=[], open_items=[])
    assert _mode(d) == 0o700
    assert _mode(fresh_home / ".delfin") == 0o700


def test_the_first_state_dir_makes_the_root_owner_only(
        group_umask, fresh_home):
    from delfin.agent import state_paths
    target = fresh_home / ".delfin" / "agent_sessions" / "deep"
    state_paths.ensure_dir(target)
    root = fresh_home / ".delfin"
    assert _mode(root) == 0o700
    assert _mode(root / "agent_sessions") == 0o700
    assert _mode(target) == 0o700


def test_a_dir_outside_the_state_root_is_left_to_the_umask(
        group_umask, fresh_home, tmp_path):
    from delfin.agent import state_paths
    outside = tmp_path / "elsewhere" / "inner"
    state_paths.ensure_dir(outside)
    assert _mode(outside) == 0o700          # the level asked for
    assert _mode(outside.parent) == 0o775   # not ours to tighten
