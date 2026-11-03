"""Git writes to a checkout, serialized across processes: UI-free.

The dashboard's header held these helpers (``delfin/dashboard/
__init__.py``) for its Pull and Rollback buttons; the agent's branch
operations (``delfin/agent/branch_ops.py``) make the same kinds of write
from the Pi's surfaces. One checkout, one lock: both call the functions
here, so a dashboard Pull and a service pull cannot interleave, and a
stale ``index.lock`` is cleared the same careful way wherever the write
comes from.

The functions moved here unchanged (names kept, so the dashboard's own
imports and the tests that pin them keep working); the dashboard imports
them back from this module.
"""

from __future__ import annotations

import fcntl
import os
import subprocess
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path


def run_git_capture(repo_dir, *args):
    """Run read-only git in the given repo without taking optional index
    locks, through safe_git (no program the repository names runs);
    stderr is folded into stdout."""
    from delfin.agent import safe_git
    return safe_git.run(args, repo_dir, timeout=60, merge_stderr=True)


@contextmanager
def dashboard_git_write_lock(repo_dir, action_name, timeout=15.0):
    """Serialize repo-changing dashboard Git actions across kernels/processes."""
    lock_path = git_dir(repo_dir) / 'delfin-dashboard.lock'
    with open(lock_path, 'a+', encoding='utf-8') as handle:
        deadline = time.monotonic() + timeout
        while True:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                handle.seek(0)
                handle.truncate()
                handle.write(f'{os.getpid()} {action_name}\n')
                handle.flush()
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise RuntimeError(
                        'Another DELFIN dashboard session is already changing the Git repo. '
                        'Please wait a moment and retry.'
                    )
                time.sleep(0.2)
        try:
            yield
        finally:
            try:
                handle.seek(0)
                handle.truncate()
                handle.flush()
            except Exception:
                pass
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except Exception:
                pass


def git_dir(repo_dir):
    """This checkout's own git directory.

    ``<repo>/.git`` in a main checkout; in a ``git worktree`` checkout
    ``.git`` is a file and the directory (index, index.lock) is
    ``<common>/.git/worktrees/<name>``, which only git can name."""
    plain = Path(repo_dir) / '.git'
    if plain.is_dir():
        return plain
    result = run_git_capture(repo_dir, 'rev-parse', '--absolute-git-dir')
    out = (result.stdout or '').strip()
    if result.returncode == 0 and out and Path(out).is_dir():
        return Path(out)
    return plain


def git_index_lock_path(repo_dir):
    """Return the Git index lock path for the repository."""
    return git_dir(repo_dir) / 'index.lock'


def git_index_lock_message(repo_dir):
    """Describe a present Git index lock, or return an empty string."""
    lock_path = git_index_lock_path(repo_dir)
    if not lock_path.exists():
        return ''
    try:
        stamp = datetime.fromtimestamp(lock_path.stat().st_mtime).strftime('%Y-%m-%d %H:%M:%S')
        return f'Git index lock present: {lock_path} (mtime {stamp})'
    except Exception:
        return f'Git index lock present: {lock_path}'


def git_index_lock_in_use(repo_dir):
    """Return whether another process currently holds the Git index lock."""
    lock_path = git_index_lock_path(repo_dir)
    if not lock_path.exists():
        return False
    try:
        result = subprocess.run(
            ['lsof', str(lock_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        )
    except Exception:
        return False
    return result.returncode == 0


def ensure_git_index_lock_clear(repo_dir, action_name, wait_seconds=10.0):
    """Wait for an active Git index lock or remove it when it is stale."""
    lock_path = git_index_lock_path(repo_dir)
    removed_stale_lock = False
    deadline = time.monotonic() + wait_seconds
    while lock_path.exists():
        if git_index_lock_in_use(repo_dir):
            if time.monotonic() >= deadline:
                raise RuntimeError(
                    f'{action_name} blocked: another git process is still using {lock_path}. '
                    'Please finish that operation or close the other DELFIN session.'
                )
            time.sleep(0.2)
            continue
        try:
            lock_path.unlink()
            removed_stale_lock = True
        except FileNotFoundError:
            break
        except Exception as exc:
            raise RuntimeError(
                f'{action_name} blocked: stale Git lock could not be removed ({exc}).'
            ) from exc
    return removed_stale_lock


def git_output_mentions_index_lock(output):
    """Return whether Git output points to an index lock conflict."""
    message = (output or '').lower()
    return 'index.lock' in message or 'could not write index' in message


def run_git_write(repo_dir, *args, retry_on_index_lock=True):
    """Run a repo-changing Git command with stale-index-lock recovery."""
    from delfin.agent import safe_git
    ensure_git_index_lock_clear(repo_dir, 'Git write')
    # safe_git: no hook, filter or other program the checkout's own config
    # names runs; the user's global credential helper still answers a pull.
    result = safe_git.run(args, repo_dir, timeout=600, merge_stderr=True)
    if result.returncode == 0 or not retry_on_index_lock:
        return result
    if not git_output_mentions_index_lock(result.stdout):
        return result
    ensure_git_index_lock_clear(repo_dir, 'Git write retry', wait_seconds=1.0)
    return run_git_write(repo_dir, *args, retry_on_index_lock=False)
