"""Where the agent service runs is the user's own, and holds one service.

tmux's rules for its socket directory: a symlink, another owner or a
non-directory is refused outright, loose permissions on the user's own
directory are tightened. One service per host holds the lock; the run file
says where it runs and is removed only by the service it names.
"""

from __future__ import annotations

import os
import socket
import stat

import pytest

from delfin.agent.service import rundir


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "rt"))
    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    (tmp_path / "rt").mkdir()
    return tmp_path / "rt"


def test_the_runtime_dir_is_created_owner_only(runtime):
    path = rundir.ensure_runtime_dir()
    assert path == runtime / f"delfin-{os.getuid()}"
    assert stat.S_IMODE(path.stat().st_mode) == 0o700


def test_loose_permissions_on_the_users_own_dir_are_tightened(runtime):
    path = runtime / f"delfin-{os.getuid()}"
    path.mkdir(mode=0o777)
    os.chmod(path, 0o777)
    rundir.ensure_runtime_dir()
    assert stat.S_IMODE(path.stat().st_mode) == 0o700


def test_a_symlink_in_its_place_is_refused(runtime, tmp_path):
    target = tmp_path / "elsewhere"
    target.mkdir()
    (runtime / f"delfin-{os.getuid()}").symlink_to(target)
    with pytest.raises(rundir.UnsafeRuntimeDir):
        rundir.ensure_runtime_dir()


def test_a_file_in_its_place_is_refused(runtime):
    (runtime / f"delfin-{os.getuid()}").write_text("not a directory")
    with pytest.raises(rundir.UnsafeRuntimeDir):
        rundir.ensure_runtime_dir()


def test_one_service_per_host_holds_the_lock(runtime):
    first, second = rundir.ServiceLock(), rundir.ServiceLock()
    first.acquire()
    try:
        with pytest.raises(rundir.ServiceBusy):
            second.acquire()
    finally:
        first.release()
    second.acquire()                    # free again once the first lets go
    second.release()


def test_the_run_file_is_the_users_and_only_its_service_removes_it(runtime):
    path = rundir.write_run_file({"pid": 4242, "socket": "/x/agent.sock"})
    assert stat.S_IMODE(path.stat().st_mode) & 0o077 == 0
    assert rundir.read_run_file()["pid"] == 4242
    rundir.remove_run_file(1)           # another pid: left alone
    assert rundir.read_run_file() is not None
    rundir.remove_run_file(4242)
    assert rundir.read_run_file() is None


def test_a_run_file_names_this_shells_service_only_for_its_socket(runtime, monkeypatch):
    # Another DELFIN_RUNTIME_DIR (the same HOME) writes the same run file;
    # its service is not this shell's to wait for, name or stop.
    monkeypatch.setattr(rundir, "_looks_like_the_service", lambda pid: True)
    rundir.write_run_file({"pid": os.getpid(), "socket": "/elsewhere/agent.sock"})
    assert rundir.service_pid() == 0
    rundir.write_run_file({"pid": os.getpid(), "socket": str(rundir.socket_path())})
    assert rundir.service_pid() == os.getpid()


def test_services_on_other_login_nodes_are_found_by_their_run_files(runtime):
    here = rundir.write_run_file({"pid": 1, "host": rundir.hostname()})
    other = here.parent / "login2.json"
    other.write_text('{"pid": 77, "host": "login2"}')
    # Rosters sit beside the run files and are never taken for one: not
    # today's (<host>.roster), nor one an older service left (.roster.json).
    rundir.roster_file().write_text('{"pid": 1, "sessions": []}')
    (here.parent / f"{rundir.hostname()}.roster.json").write_text('{"pid": 1, "sessions": []}')
    assert [h["host"] for h in rundir.other_hosts()] == ["login2"]


def test_the_peer_of_a_socket_is_known_by_its_uid():
    a, b = socket.socketpair(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        assert rundir.peer_uid(a) == os.getuid()
    finally:
        a.close()
        b.close()
