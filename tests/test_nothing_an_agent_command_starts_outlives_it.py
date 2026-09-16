"""Nothing an agent's shell command starts outlives it.

The lifeline ledger ends what DELFIN itself started detached. A command
can start its own: `setsid`, `nohup`, `disown`, a double fork, a tmux or
screen server, `systemd-run --user`, a command over SSH on another machine.
None of those were on any ledger, and in the unattended mode with
filesystem isolation off (a real user's setting) nothing contained them.

Every shell command now runs in a process cage whatever the filesystem
mode: its own PID namespace, ended with the command and with the process
that started it, with the doors to long-lived servers outside masked. The
one way out is a plain push to GitHub, run outside in a hardened form.
"""

from __future__ import annotations

import os
import random
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as A

requires_the_cage = pytest.mark.skipif(
    not A._process_cage_functional(),
    reason="bwrap cannot build the process cage here")


def _perms(tmp_path, mode="bypassPermissions"):
    return A.KitToolPermissions(workspace=str(tmp_path), mode=mode)


def _marks(n: int) -> list[str]:
    base = random.randint(40000, 89999)
    return [str(base + i) for i in range(n)]


def _living(marks: list[str]) -> list[str]:
    """``sleep <mark>`` processes of this user that are still there."""
    found = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            argv = Path(f"/proc/{entry}/cmdline").read_bytes().split(b"\0")
            state = Path(f"/proc/{entry}/stat").read_text().rsplit(")", 1)[1].split()[0]
        except OSError:
            continue
        if state != "Z" and len(argv) >= 2 and argv[0].endswith(b"sleep") \
                and argv[1].decode(errors="replace") in marks:
            found.append(argv[1].decode())
    return found


def _gone(marks: list[str], seconds: float = 10.0) -> bool:
    deadline = time.time() + seconds
    while time.time() < deadline:
        if not _living(marks):
            return True
        time.sleep(0.1)
    return False


@requires_the_cage
def test_every_way_of_leaving_something_behind_ends_with_the_command(tmp_path):
    a, b, c, d = marks = _marks(4)
    cmd = (f"setsid sleep {a} & disown; "
           f"nohup sleep {b} >/dev/null 2>&1 & "
           f"(sleep {c} &); "
           f"python3 -c 'import os,time\n"
           f"if os.fork()==0:\n os.setsid()\n if os.fork()==0: time.sleep({d})\n os._exit(0)'; "
           "sleep 0.5; echo started")
    argv = A._bash_isolation_argv(cmd, tmp_path, _perms(tmp_path), mode="off")
    try:
        done = subprocess.run(argv, cwd=tmp_path, capture_output=True, text=True,
                              timeout=60)
        assert "started" in done.stdout, done.stderr
        assert _gone(marks), f"outlived the command: {_living(marks)}"
    finally:
        subprocess.run(["pkill", "-u", str(os.getuid()), "-f",
                        "sleep (" + "|".join(marks) + ")$"])


@requires_the_cage
def test_the_cage_ends_with_the_process_that_started_it(tmp_path):
    a, b = marks = _marks(2)
    starter = (
        "import subprocess, sys, time\n"
        "subprocess.Popen(sys.argv[1:])\n"
        "print('ready', flush=True)\n"
        "time.sleep(60)\n")
    argv = A._in_process_cage(
        ["/bin/bash", "-c", f"setsid sleep {a} & sleep {b}"], tmp_path)
    kernel = subprocess.Popen([sys.executable, "-c", starter, *argv],
                              stdout=subprocess.PIPE, text=True)
    try:
        assert kernel.stdout.readline().strip() == "ready"
        deadline = time.time() + 10
        while len(_living(marks)) < 2 and time.time() < deadline:
            time.sleep(0.1)
        assert len(_living(marks)) == 2
        kernel.send_signal(signal.SIGKILL)
        kernel.wait()
        assert _gone(marks), f"outlived the kernel: {_living(marks)}"
    finally:
        kernel.kill()
        subprocess.run(["pkill", "-u", str(os.getuid()), "-f",
                        "sleep (" + "|".join(marks) + ")$"])


@requires_the_cage
def test_a_background_job_outlives_the_turn_that_started_it(tmp_path):
    (mark,) = marks = _marks(1)
    started = []
    turn = threading.Thread(target=lambda: started.append(subprocess.Popen(
        A._in_process_cage(["sleep", mark], tmp_path))))
    turn.start()
    turn.join()
    try:
        time.sleep(1.5)
        assert started[0].poll() is None, "ended with the thread of its turn"
    finally:
        started[0].kill()
        started[0].wait()
        assert _gone(marks)


@requires_the_cage
def test_a_background_job_is_caged_and_its_kill_takes_everything(tmp_path):
    from delfin.agent import bash_jobs
    a, b = marks = _marks(2)
    cmd = f"setsid sleep {a} & sleep {b}"
    reg = bash_jobs.get_registry()
    job = reg.start(command=cmd, cwd=str(tmp_path), workspace=str(tmp_path),
                    argv=A._bash_isolation_argv(cmd, tmp_path, _perms(tmp_path),
                                                mode="off"))
    try:
        assert job.proc.args[0] == "bwrap"
        deadline = time.time() + 10
        while len(_living(marks)) < 2 and time.time() < deadline:
            time.sleep(0.1)
        reg.kill(job.job_id)
        assert _gone(marks), f"outlived bash_kill: {_living(marks)}"
    finally:
        reg.kill(job.job_id)


@requires_the_cage
def test_the_doors_to_servers_outside_are_closed(tmp_path, monkeypatch):
    agent_dir = tmp_path / "ssh-probe"
    agent_dir.mkdir()
    (agent_dir / "agent.1").write_text("")
    monkeypatch.setenv("SSH_AUTH_SOCK", str(agent_dir / "agent.1"))
    monkeypatch.setenv("TMUX", "/tmp/tmux-0/default,1,0")
    uid = os.getuid()
    cmd = (f"echo sock=${{SSH_AUTH_SOCK:-unset}} tmux=${{TMUX:-unset}}; "
           f"echo agentdir=$(ls -A {agent_dir} | wc -l); "
           f"echo runuser=$(ls -A /run/user/{uid} 2>/dev/null | wc -l); "
           "echo ssh=$(ls -A ~/.ssh 2>/dev/null | wc -l); "
           "tty >/dev/null 2>&1 && echo tty=yes || echo tty=no; "
           "command -v systemd-run >/dev/null && "
           "{ systemd-run --user --quiet true >/dev/null 2>&1 "
           "&& echo systemd=open || echo systemd=closed; } || echo systemd=closed")
    argv = A._bash_isolation_argv(cmd, tmp_path, _perms(tmp_path), mode="off")
    out = subprocess.run(argv, cwd=tmp_path, capture_output=True, text=True,
                         timeout=60).stdout
    for expected in ("sock=unset", "tmux=unset", "agentdir=0", "runuser=0", "ssh=0",
                     "tty=no", "systemd=closed"):
        assert expected in out, out


def test_filesystem_isolation_off_does_not_turn_the_cage_off(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_process_cage_functional", lambda: True)
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    for mode in ("off", "auto"):
        argv = A._bash_isolation_argv("echo hi", tmp_path,
                                      _perms(tmp_path, "default"), mode=mode)
        assert argv[0] == "bwrap" and "--unshare-pid" in argv
        assert argv[-3:] == ["/bin/bash", "-c", "echo hi"]


def test_only_the_environment_of_the_terminal_turns_the_cage_off(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_process_cage_functional", lambda: True)
    monkeypatch.setenv(A._PROCESS_CAGE_ENV, "off")
    assert A._bash_isolation_argv("echo hi", tmp_path, _perms(tmp_path),
                                  mode="off") == ["/bin/bash", "-c", "echo hi"]


def test_filesystem_isolation_and_the_cage_combine(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_process_cage_functional", lambda: True)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    # The host is supplied, not measured: the forced mode also asks
    # shutil.which, and the CI runner has no bubblewrap installed.
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    argv = A._bash_isolation_argv("echo hi", tmp_path, _perms(tmp_path),
                                  mode="bwrap")
    assert argv[:3] == ["bwrap", "--ro-bind", "/"]
    for option in ("--unshare-pid", "--new-session", "--die-with-parent"):
        assert argv.count(option) == 1, option


def test_where_the_cage_cannot_be_built_it_is_said(tmp_path, monkeypatch):
    events = []
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    monkeypatch.setattr(A, "_PROCESS_CAGE_GAP_ANNOUNCED", False)
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: events.append(a))
    assert A._bash_isolation_argv("echo hi", tmp_path, _perms(tmp_path),
                                  mode="off") == ["/bin/bash", "-c", "echo hi"]
    assert events and "process cage is NOT active" in events[0][2]


def test_the_background_tool_runs_its_command_in_the_cage():
    import ast
    import inspect
    import textwrap
    src = ast.unparse(ast.parse(textwrap.dedent(inspect.getsource(
        A._DocToolExecutor._execute_bash_background))))
    assert "argv=_bash_isolation_argv(cmd, run_cwd, perms)" in src


# -- the push to GitHub ---------------------------------------------------------

@pytest.fixture
def repo(tmp_path):
    def git(*args):
        subprocess.run(["git", *args], cwd=tmp_path, check=True,
                       capture_output=True)
    git("init", "-q", "-b", "main")
    git("remote", "add", "origin", "git@github.com:ComPlat/DELFIN.git")
    git("remote", "add", "cluster", "ka_user@uc3n991:repo.git")
    return tmp_path, git


@pytest.mark.parametrize("cmd", [
    "git push origin main",
    "git push -u origin agent/topic",
    "git push origin main 2>&1 | tail -3",
    "git push origin HEAD:refs/heads/agent/topic",
])
def test_a_plain_push_to_github_is_the_one_way_out(repo, cmd):
    path, _git = repo
    plan = A._github_push_plan(cmd, path)
    assert plan is not None and plan["url"] == "git@github.com:ComPlat/DELFIN.git"
    argv = plan["argv"]
    assert argv[0] == "git" and "push" in argv
    for setting in ("core.hooksPath=/dev/null", "protocol.allow=never",
                    "protocol.ssh.allow=always", "push.recurseSubmodules=no"):
        assert setting in argv
    ssh = next(a for a in argv if a.startswith("core.sshCommand="))
    for option in ("ControlPath=none", "ProxyCommand=none", "ProxyJump=none",
                   "PermitLocalCommand=no", "BatchMode=yes"):
        assert option in ssh


@pytest.mark.parametrize("cmd", [
    "git push origin main; touch x",
    "git push origin main && echo done",
    "cd sub && git push origin main",
    "git push origin +main",
    "git push origin :main",
    "git push --force origin main",
    "git push --receive-pack=/tmp/x origin main",
    "git push --exec=/tmp/x origin main",
    "git push cluster main",
    "git push git@github.com:ComPlat/DELFIN.git main",
    "git push origin main | sh",
    "git push origin main 2>&1 | tail -3 | sh",
])
def test_anything_else_stays_in_the_cage(repo, cmd):
    path, _git = repo
    assert A._github_push_plan(cmd, path) is None


def test_a_remote_rewritten_away_from_github_stays_in_the_cage(repo):
    path, git = repo
    git("config", "url.ka_user@uc3n991:.pushInsteadOf", "git@github.com:ComPlat/")
    assert A._github_push_plan("git push origin main", path) is None


def test_a_remote_with_a_helper_stays_in_the_cage(repo):
    path, git = repo
    git("config", "remote.origin.vcs", "evil")
    assert A._github_push_plan("git push origin main", path) is None


def test_the_push_is_reported_the_way_the_shell_would(monkeypatch, tmp_path):
    seen = {}

    def fake_run(argv, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(argv, 1, "a\nb\nc\nd\n", "rejected\n")

    monkeypatch.setattr(A.subprocess, "run", fake_run)
    plan = {"argv": ["git", "push"], "merge_stderr": True, "filter": "tail",
            "lines": 2}
    done = A._run_github_push(plan, tmp_path, {"GIT_SSH_COMMAND": "evil",
                                               "PATH": "/usr/bin"}, 30)
    assert done.stdout == "d\nrejected\n" and done.stderr == ""
    assert done.returncode == 0, "a pipeline's status is its last command's"
    assert "GIT_SSH_COMMAND" not in seen["env"]
    assert seen["start_new_session"] and seen["stdin"] is subprocess.DEVNULL
