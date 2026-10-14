"""Regression test: dashboard git-lock paths must survive linked worktrees.

In a linked git worktree, ``<worktree>/.git`` is a plain text file containing
``gitdir: <path>``, NOT a directory. The dashboard lock helpers used to build
``Path(repo_dir) / '.git' / 'delfin-dashboard.lock'`` and
``Path(repo_dir) / '.git' / 'index.lock'``, which treats the pointer file as a
directory and fails with ``[Errno 20] Not a directory`` the moment the lock is
opened. The fix resolves the real git dir via :func:`_git_dir_path`.
"""

import os

from delfin.dashboard import _git_dir_path, _git_index_lock_path


def _fake_worktree(tmp_path, real_git_dir):
    """Return a worktree dir whose ``.git`` is a ``gitdir:`` pointer file."""
    worktree = tmp_path / 'wt-some-branch'
    worktree.mkdir()
    (worktree / '.git').write_text(
        f'gitdir: {real_git_dir}\n', encoding='utf-8'
    )
    return worktree


def test_git_dir_path_resolves_linked_worktree_pointer(tmp_path):
    """A .git that is a gitdir: pointer file must resolve to the real git dir."""
    real_git = tmp_path / 'real-git-dir'
    real_git.mkdir()
    worktree = _fake_worktree(tmp_path, real_git)

    assert _git_dir_path(worktree) == real_git
    # the lock lands in the REAL git dir, never "inside" the pointer file
    assert _git_dir_path(worktree) / 'delfin-dashboard.lock' == (
        real_git / 'delfin-dashboard.lock'
    )


def test_git_dir_path_returns_dot_git_for_normal_checkout(tmp_path):
    """A normal checkout (real .git directory) keeps Path(repo)/.git."""
    repo = tmp_path / 'normal'
    repo.mkdir()
    (repo / '.git').mkdir()

    assert _git_dir_path(repo) == repo / '.git'


def test_worktree_lock_path_is_not_built_inside_the_pointer_file(tmp_path):
    """The dashboard lock must not join through a .git pointer file (ENOTDIR)."""
    real_git = tmp_path / 'real-git-dir'
    real_git.mkdir()
    worktree = _fake_worktree(tmp_path, real_git)

    lock_path = _git_dir_path(worktree) / 'delfin-dashboard.lock'

    # Opening this must succeed (the bug produced ENOTDIR here).
    with open(lock_path, 'a+', encoding='utf-8') as handle:
        handle.write(f'{os.getpid()} test-dashboard-lock\n')
        handle.flush()
    assert lock_path.exists()
    assert lock_path.parent == real_git


def test_index_lock_path_resolves_worktree_pointer(tmp_path):
    """index.lock uses the real git dir in a linked worktree too."""
    real_git = tmp_path / 'real-git-dir'
    real_git.mkdir()
    worktree = _fake_worktree(tmp_path, real_git)

    assert _git_index_lock_path(worktree) == real_git / 'index.lock'
