"""Controls for SIZE and CLEANUP of the shared temp dir (Phase 3, LH).

Red on the previous commit: ``shared_tmp`` has no quota accounting, no
session-end cleanup and no ageing -- a session could fill the runtime
tmpfs without limit, and sessions that crash never clear their
directory.

Contract pinned here:
- ``usage_bytes(dir)``: recursive size of ONE session dir only
  (never walks anything outside it -- no shared-tree traversal);
- ``enforce_quota(dir, max_bytes)``: reports overage, deletes the
  OLDEST files first (mtime), stops as soon as the quota is met, never
  touches anything outside the dir (symlinks are not followed);
- ``cleanup_session(dir)``: removes the whole session dir at session
  end (the one place deletion is allowed);
- ``sweep_stale(base, max_age_s, now=None)``: ageing -- removes
  SESSION dirs under <base>/tmp whose mtime is older than max_age,
  never anything else, never the 'tmp' collection dir itself.
"""

import os
import time
from pathlib import Path

import pytest

from delfin.agent import shared_tmp


def make_session(tmp_path, name="sess", files=()):
    base = tmp_path / "runtime"
    d = shared_tmp.session_dir(name, base=base)
    shared_tmp.ensure_session_dir(d)
    for fname, size in files:
        (d / fname).write_bytes(b"x" * size)
    return base, d


# ---------------------------------------------------------------------------
# usage_bytes
# ---------------------------------------------------------------------------

def test_usage_counts_only_inside_the_session_dir(tmp_path):
    base, d = make_session(tmp_path, files=[("a.bin", 100)])
    other = shared_tmp.session_dir("other", base=base)
    shared_tmp.ensure_session_dir(other)
    (other / "big.bin").write_bytes(b"y" * 1000)
    assert shared_tmp.usage_bytes(d) == 100


def test_usage_is_recursive(tmp_path):
    base, d = make_session(tmp_path)
    sub = d / "deep" / "er"
    sub.mkdir(parents=True)
    (sub / "c.bin").write_bytes(b"z" * 50)
    (d / "a.bin").write_bytes(b"x" * 10)
    assert shared_tmp.usage_bytes(d) == 60


def test_usage_does_not_follow_symlinks_out(tmp_path):
    base, d = make_session(tmp_path)
    secret = tmp_path / "secret.bin"
    secret.write_bytes(b"s" * 5000)
    (d / "link").symlink_to(secret)
    assert shared_tmp.usage_bytes(d) == 0  # symlink counts 0, target
    # is outside the session dir and must not be measured or read


# ---------------------------------------------------------------------------
# enforce_quota
# ---------------------------------------------------------------------------

def test_quota_under_limit_is_a_no_op(tmp_path):
    base, d = make_session(tmp_path, files=[("a.bin", 100)])
    report = shared_tmp.enforce_quota(d, max_bytes=1000)
    assert (d / "a.bin").exists()
    assert report.deleted == []


def test_quota_over_limit_deletes_oldest_first(tmp_path):
    base, d = make_session(tmp_path)
    old = d / "old.bin"
    old.write_bytes(b"x" * 600)
    os.utime(old, (time.time() - 1000, time.time() - 1000))
    new = d / "new.bin"
    new.write_bytes(b"x" * 600)
    report = shared_tmp.enforce_quota(d, max_bytes=700)
    assert not old.exists()      # older one gone
    assert new.exists()          # newer one kept
    assert shared_tmp.usage_bytes(d) <= 700
    assert report.deleted == [old]


def test_quota_never_leaves_the_session_dir(tmp_path):
    base, d = make_session(tmp_path, files=[("a.bin", 10)])
    outside = tmp_path / "outside.bin"
    outside.write_bytes(b"q" * 10)
    (d / "esc").symlink_to(outside)
    report = shared_tmp.enforce_quota(d, max_bytes=0)
    assert outside.exists()      # untouched through the symlink
    assert report.deleted        # something inside WAS deleted


# ---------------------------------------------------------------------------
# cleanup_session
# ---------------------------------------------------------------------------

def test_cleanup_removes_only_the_session_dir(tmp_path):
    base, d = make_session(tmp_path, files=[("a.bin", 10)])
    neighbor = shared_tmp.session_dir("neighbor", base=base)
    shared_tmp.ensure_session_dir(neighbor)
    shared_tmp.cleanup_session(d)
    assert not d.exists()
    assert neighbor.exists()
    assert (base / "tmp").exists()  # collection dir survives


def test_cleanup_of_missing_dir_is_silent(tmp_path):
    shared_tmp.cleanup_session(tmp_path / "runtime" / "tmp" / "gone")


# ---------------------------------------------------------------------------
# sweep_stale (ageing)
# ---------------------------------------------------------------------------

def test_sweep_removes_only_stale_session_dirs(tmp_path):
    base, d_stale = make_session(tmp_path, name="stale")
    base, d_fresh = make_session(tmp_path, name="fresh")
    old = time.time() - 10 * 86400
    os.utime(d_stale, (old, old))
    removed = shared_tmp.sweep_stale(base, max_age_s=86400)
    assert not d_stale.exists()
    assert d_fresh.exists()
    assert removed == [d_stale]


def test_sweep_with_now_argument_is_deterministic(tmp_path):
    base, d = make_session(tmp_path, name="s")
    now = time.time() + 86400 * 3   # three days later
    removed = shared_tmp.sweep_stale(base, max_age_s=86400, now=now)
    assert removed == [d]
    assert not d.exists()
