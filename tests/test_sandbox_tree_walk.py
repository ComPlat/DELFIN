"""Recursive walks rooted at a whole file system must not reach the shell.

The rules are derived from HOME and the mount table at runtime, so these tests
fake both -- nothing here may depend on the machine they run on.
"""

from pathlib import Path

import pytest

from delfin.agent import sandbox


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    """Pretend HOME is /cluster/site/group/someone, on an unrelated machine."""
    home = tmp_path / "cluster" / "site" / "group" / "someone"
    (home / "jobs" / "run01").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(sandbox.os.path, "ismount", lambda p: False)
    return home


@pytest.mark.parametrize(
    "cmd",
    [
        "du -sb {home}",
        "du -sh {home}",
        "du -s {home}/",
        "find {home} -name '*.out'",
    ],
)
def test_walk_of_home_is_denied(fake_home, cmd):
    result = sandbox.is_allowed(cmd.format(home=fake_home))
    assert not result.allowed
    assert "metadata servers" in result.reason


def test_ncdu_is_denied_even_though_it_never_reaches_the_walk_guard(fake_home):
    """ncdu is not allow-listed at all; it must stay unreachable either way."""
    assert not sandbox.is_allowed(f"ncdu {fake_home}").allowed


def test_walk_of_home_ancestors_is_denied(fake_home):
    for ancestor in (fake_home.parent, fake_home.parent.parent):
        result = sandbox.is_allowed(f"du -sb {ancestor}")
        assert not result.allowed, f"{ancestor} should be denied"


def test_walk_of_root_is_denied(fake_home):
    assert not sandbox.is_allowed("du -sb /").allowed


def test_walk_of_a_job_subdirectory_is_allowed(fake_home):
    """The agent's real use case must keep working."""
    assert sandbox.is_allowed(f"du -sh {fake_home}/jobs/run01").allowed
    assert sandbox.is_allowed(f"find {fake_home}/jobs/run01 -name '*.out'").allowed


def test_tilde_is_expanded_before_the_check(fake_home, monkeypatch):
    monkeypatch.setenv("HOME", str(fake_home))
    assert not sandbox.is_allowed("du -sb ~").allowed
    assert sandbox.is_allowed("du -sh ~/jobs/run01").allowed


def test_mount_points_are_denied_without_hardcoding_their_names(
    monkeypatch, tmp_path, fake_home
):
    """A project or scratch file system is caught by the mount-table signal."""
    scratch = tmp_path / "scratch"
    (scratch / "user" / "job").mkdir(parents=True)
    monkeypatch.setattr(
        sandbox.os.path, "ismount", lambda p: str(p) == str(scratch)
    )

    assert not sandbox.is_allowed(f"du -sb {scratch}").allowed
    # One level down is not a mount point and stays usable.
    assert sandbox.is_allowed(f"du -sb {scratch}/user/job").allowed


def test_bare_du_checks_the_working_directory(fake_home, monkeypatch):
    """`du -sh` with no argument walks CWD; that must be judged too."""
    monkeypatch.chdir(fake_home)
    assert not sandbox.is_allowed("du -sh").allowed

    monkeypatch.chdir(fake_home / "jobs" / "run01")
    assert sandbox.is_allowed("du -sh").allowed


def test_denial_survives_inside_a_pipeline(fake_home):
    """Hiding the walk behind a pipe must not get past the check."""
    assert not sandbox.is_allowed(f"du -sb {fake_home} | tail -1").allowed
    assert not sandbox.is_allowed(f"echo hi && du -sb {fake_home}").allowed


def test_unrelated_allowed_commands_are_unaffected(fake_home):
    for cmd in ("ls -la", "grep -rn foo delfin/", f"cat {fake_home}/.bashrc"):
        assert sandbox.is_allowed(cmd).allowed, cmd


def test_lfs_quota_is_not_blocked(fake_home):
    """The sanctioned replacement must remain reachable if a site allows it."""
    result = sandbox.is_allowed("lfs quota -q -u someone /home/someone")
    # 'lfs' may or may not be on the allow-list; what matters is that the
    # tree-walk guard is not what stops it.
    assert "metadata servers" not in result.reason
