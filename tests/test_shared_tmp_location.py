"""Controls for the shared-tmp LOCATION decision (Phase 1, LH run).

Red on the previous commit: ``shared_tmp`` has no ``location_candidates``
/ ``resolve_location`` -- the directory was fixed under the HOME state
root, which the operator deferred (quota on the cluster HOME, and a
shared-state-root dir is not cleaned by boot).

The decision these tests pin (reasoning in .gate/ORT.md):
1. ``DELFIN_TMP_ROOT`` -- explicit override, wins unconditionally.
2. ``$XDG_RUNTIME_DIR`` -- per-user tmpfs, wiped at logout, quota-free.
3. ``/run/user/<uid>`` -- the same place when XDG is not exported.
4. state root (``~/.delfin``) -- last resort, always present.
Never the host's ``/tmp``, never ``$TMPDIR`` (here a shared cluster FS).
"""

import os
from pathlib import Path

import pytest

from delfin.agent import shared_tmp


def test_xdg_runtime_dir_is_the_preferred_candidate(monkeypatch, tmp_path):
    xdg = tmp_path / "xdg"
    xdg.mkdir()
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(xdg))
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    cands = shared_tmp.location_candidates()
    assert cands[0] == xdg


@pytest.mark.skipif(
    not Path(f"/run/user/{os.getuid()}").is_dir(),
    reason="/run/user/<uid> does not exist on this machine",
)
def test_without_xdg_the_run_user_dir_is_used(monkeypatch):
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    uid = os.getuid()
    cands = shared_tmp.location_candidates()
    assert cands[0] == Path(f"/run/user/{uid}")


def test_explicit_override_wins_over_everything(monkeypatch, tmp_path):
    xdg = tmp_path / "xdg"
    xdg.mkdir()
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(xdg))
    monkeypatch.setenv("DELFIN_TMP_ROOT", str(tmp_path / "override"))
    assert shared_tmp.location_candidates()[0] == tmp_path / "override"


def test_a_candidate_that_is_not_a_directory_is_skipped(monkeypatch, tmp_path):
    plain = tmp_path / "xdg"  # a FILE, not a directory
    plain.write_text("not a dir", encoding="utf-8")
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(plain))
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    cands = shared_tmp.location_candidates()
    assert plain not in cands


def test_a_symlinked_candidate_is_skipped(monkeypatch, tmp_path):
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "xdg"
    link.symlink_to(real)
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(link))
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    cands = shared_tmp.location_candidates()
    assert link not in cands


def test_home_state_root_is_the_last_candidate_never_the_first(monkeypatch,
                                                               tmp_path):
    xdg = tmp_path / "xdg"
    xdg.mkdir()
    state = tmp_path / "state"
    state.mkdir()
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(xdg))
    monkeypatch.setenv("DELFIN_STATE", str(state))
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    # uid=1: /run/user/1 does not exist, so the chain really is
    # [xdg, state] and the state root is LAST, not first.
    cands = shared_tmp.location_candidates(uid=1)
    assert cands[0] == xdg
    assert cands[-1] == state
    assert shared_tmp.state_root() not in cands[:-1]


def test_the_host_tmp_and_tmpdir_are_never_candidates(monkeypatch, tmp_path):
    other = tmp_path / "tmpdir"  # a perfectly good dir, but only TMPDIR
    other.mkdir()
    xdg = tmp_path / "xdg"
    xdg.mkdir()
    monkeypatch.setenv("TMPDIR", str(other))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(xdg))
    monkeypatch.delenv("DELFIN_TMP_ROOT", raising=False)
    cands = shared_tmp.location_candidates()
    assert Path("/tmp") not in cands
    assert other not in cands  # TMPDIR must not leak into the chain
