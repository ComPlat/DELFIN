"""Controls for SECURE CREATION of the shared temp dir (Phase 2, LH run).

Red on the previous commit: ``shared_tmp`` cannot create the session
directory at all -- the operator's second objection to the Welle-5
design was that a pre-created directory or symlink must be DETECTED and
NOT used (owner check, 0700, no symlink anywhere in the path,
O_NOFOLLOW-style lstat/fstat). All scenarios run inside pytest's
``tmp_path`` -- no foreign users involved.

The contract pinned here for ``ensure_session_dir(session_dir_path)``:
- missing: created fresh, every new component 0700, owned by us;
- already ours and 0700: reused (idempotent -- cage and file tools
  both call this);
- pre-placed by "someone else" (wrong owner cannot be simulated
  portably; wrong PERMISSIONS stand in for it, plus chown when root):
  refused, not used;
- a symlink or a file at the target: refused;
- a symlink anywhere in the ancestor path: refused.
"""

import os
import stat
from pathlib import Path

import pytest

from delfin.agent import shared_tmp


def fresh(tmp_path):
    return tmp_path / "runtime" / "tmp" / "sess"


def test_missing_dir_is_created_0700(tmp_path):
    d = fresh(tmp_path)
    shared_tmp.ensure_session_dir(d)
    assert d.is_dir()
    mode = stat.S_IMODE(d.stat().st_mode)
    assert mode == 0o700
    # parent components created by us are 0700 too
    assert stat.S_IMODE(d.parent.stat().st_mode) == 0o700
    assert d.stat().st_uid == os.getuid()


def test_second_call_is_idempotent(tmp_path):
    d = fresh(tmp_path)
    shared_tmp.ensure_session_dir(d)
    shared_tmp.ensure_session_dir(d)  # must not raise, not reset
    assert d.is_dir()
    assert stat.S_IMODE(d.stat().st_mode) == 0o700


def test_preplaced_dir_with_loose_perms_is_refused(tmp_path):
    d = fresh(tmp_path)
    d.mkdir(parents=True)
    d.chmod(0o755)  # "someone else" left it world-visible
    with pytest.raises(shared_tmp.PrePlacedError):
        shared_tmp.ensure_session_dir(d)
    assert d.is_dir()  # untouched, NOT adopted


def test_preplaced_symlink_at_target_is_refused(tmp_path):
    target = fresh(tmp_path)
    target.parent.mkdir(parents=True)
    real = tmp_path / "elsewhere"
    real.mkdir()
    target.symlink_to(real)
    with pytest.raises(shared_tmp.PrePlacedError):
        shared_tmp.ensure_session_dir(target)


def test_regular_file_at_target_is_refused(tmp_path):
    d = fresh(tmp_path)
    d.parent.mkdir(parents=True)
    d.write_text("not a dir", encoding="utf-8")
    with pytest.raises(shared_tmp.PrePlacedError):
        shared_tmp.ensure_session_dir(d)


def test_symlink_in_an_ancestor_is_refused(tmp_path):
    real_parent = tmp_path / "real-tmp"
    real_parent.mkdir()
    link_parent = tmp_path / "runtime"
    link_parent.symlink_to(real_parent)  # ancestor itself is a symlink
    d = link_parent / "tmp" / "sess"
    with pytest.raises(shared_tmp.PrePlacedError):
        shared_tmp.ensure_session_dir(d)


# ---------------------------------------------------------------------------
# Integration through the public call path: session_dir() is how DELFIN
# names the directory; ensure_session_dir is how it comes into being.
# ---------------------------------------------------------------------------

def test_session_dir_can_be_ensured_via_public_path(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_TMP_ROOT", str(tmp_path / "override"))
    sd = shared_tmp.session_dir("pub-path-test")  # base=resolve_location()
    shared_tmp.ensure_session_dir(sd)
    probe = sd / "x.txt"
    probe.write_text("ok", encoding="utf-8")
    host = shared_tmp.map_to_host("/tmp/x.txt", session_dir=sd)
    assert host == probe
    assert shared_tmp.cage_mount(session_dir=sd) == sd
