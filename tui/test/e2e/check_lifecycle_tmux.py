"""The window across the service's life, live in tmux, against a private
agent service:

- kill -9 of the service under an open window: the window says it is
  reconnecting, starts the next service, and goes on with its session
  (woken from where the dead service left it): a prompt is answered;
- `delfin-agent service stop` from another shell: the window says the
  service stopped and that the sessions are saved, and leaves cleanly;
  the next service offers the session again, asleep;
- a hangup of a window (its SSH connection dropped): the window ends, its
  session runs on in the service.

Only processes this check started are signalled: the service by the pid
in its private run file, the window by its tmux pane.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_lifecycle_tmux.py
"""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[2]
REPO = TUI.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TUI / "spike"))

from check_tmux import Tmux  # noqa: E402

from delfin.agent.service.client import ServiceClient  # noqa: E402

PYTHON = os.environ.get("DELFIN_PYTHON") or sys.executable


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-life-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtl", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO), "DELFIN_TUI_SCREEN": "inline",
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    full_env = {**os.environ, **env}
    run_dir = tmp / "scratch" / ".delfin" / "agent_service" / "run"
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(100, 40)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def service_pid() -> int:
        for path in run_dir.glob("*.json"):
            try:
                pid = int(json.loads(path.read_text()).get("pid") or 0)
            except (OSError, ValueError):
                continue
            try:
                if b"delfin.agent.service.supervisor" in Path(f"/proc/{pid}/cmdline").read_bytes():
                    return pid
            except OSError:
                pass
        return 0

    def expect(pattern: str, what: str, timeout: float = 15, where: str = "history") -> bool:
        found = tmux.wait_for(pattern, timeout) if where == "history" else _screen_wait(tmux, pattern, timeout)
        if found:
            return True
        problems.append(f"{what}: no {pattern!r}\n{tmux.screen()}")
        return False

    def exited(pattern_count: int) -> str | None:
        deadline = time.monotonic() + 40
        while time.monotonic() < deadline:
            codes = re.findall(r"TUI_EXIT=(\d+)", tmux.history())
            if len(codes) >= pattern_count:
                return codes[pattern_count - 1]
            time.sleep(0.2)
        return None

    try:
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        tmux.type("before the crash")
        tmux.key("Enter")
        expect(r"\(re: before the crash\)", "the first turn", 30)
        # kill -9 of the service under the open window.
        old = service_pid()
        if not old:
            problems.append("no private service found")
            raise SystemExit
        os.kill(old, signal.SIGKILL)
        expect(r"reconnecting to the agent service", "the window's reconnecting note", 15, where="screen")
        deadline = time.monotonic() + 60
        while time.monotonic() < deadline and (not service_pid() or service_pid() == old):
            time.sleep(0.2)
        new = service_pid()
        if not new or new == old:
            problems.append("the window did not start the next service")
            raise SystemExit
        _screen_gone(tmux, r"reconnecting to the agent service", 60)
        time.sleep(1.0)
        tmux.type("after the crash")
        tmux.key("Enter")
        expect(r"\(re: after the crash\)", "a prompt after the crash", 60)
        with ServiceClient(sock) as c:
            sessions = c.request("sessions.list")["sessions"]
            said = [e["text"] for e in c.request("session.replay", {"sessionId": sessions[0]["sessionId"]})["events"]
                    if e.get("type") == "user.message"] if sessions else []
        if said != ["before the crash", "after the crash"]:
            problems.append(f"the session after the crash holds: {said}")
        sid = sessions[0]["sessionId"] if sessions else ""
        # `service stop` from another shell: the window leaves, saying why.
        stop = subprocess.run([PYTHON, "-m", "delfin.agent.cli", "service", "stop"], cwd=str(REPO),
                              env=full_env, capture_output=True, text=True, timeout=120)
        if stop.returncode != 0:
            problems.append(f"service stop failed: {stop.stdout} {stop.stderr}")
        expect(r"the agent service stopped \(stop\)", "the window's note about the stop", 30)
        code = exited(1)
        if code != "0":
            problems.append(f"the window did not leave cleanly after the stop: {code}")
        if service_pid():
            problems.append("a service runs again after `service stop` (the window must not start one)")
        # The next service offers the session again, asleep.
        start = subprocess.run([PYTHON, "-m", "delfin.agent.service.launcher"], cwd=str(REPO),
                               env=full_env, capture_output=True, text=True, timeout=120)
        if start.returncode != 0:
            problems.append(f"the launcher could not start a service: {start.stderr}")
        with ServiceClient(sock) as c:
            row = next((s for s in c.hello["sessions"] if s["sessionId"] == sid), None)
        if not row or row["state"] != "sleeping":
            problems.append(f"the next service should offer the session asleep: {row}")
        # A window that hangs up (its SSH connection dropped): it ends, its
        # session runs on.
        tmux.key("C-l")
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the second window", 60):
            raise SystemExit
        tmux.type("from the second window")
        tmux.key("Enter")
        expect(r"\(re: from the second window\)", "the second window's turn", 30)
        bun = _window_of(tmux, runtime)
        if not bun:
            problems.append("the second window's process was not found")
            raise SystemExit
        os.kill(bun, signal.SIGHUP)
        if exited(2) is None:
            problems.append("the window did not end on a hangup")
        with ServiceClient(sock) as c:
            sessions = c.request("sessions.list")["sessions"]
            alive = [s for s in sessions if s.get("lastPrompt") == "from the second window"]
            if not alive or alive[0]["state"] == "sleeping":
                problems.append(f"the hung-up window's session should run on: {sessions}")
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        pid = service_pid()
        if pid:
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux lifecycle check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the window outlives a killed service, leaves on a stop, and a hangup leaves its session running")
    return 0


def _pane_pid(tmux: Tmux) -> int:
    """The shell of this check's own tmux pane (its private server)."""
    try:
        out = tmux.run("list-panes", "-t", tmux.name, "-F", "#{pane_pid}")
        return int(out.split()[0])
    except (subprocess.CalledProcessError, IndexError, ValueError):
        return 0


def _window_of(tmux: Tmux, runtime: Path) -> int:
    """The window this check started: a direct child of its own pane's
    shell, running main.tsx or the built window (the launcher prefers a
    fresh build), with this check's private runtime directory. Nothing
    else is ever signalled (0 when it is not found)."""
    shell = _pane_pid(tmux)
    if shell <= 1:
        return 0
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
            if int(stat.rsplit(")", 1)[1].split()[1]) != shell:
                continue
            cmdline = (entry / "cmdline").read_bytes()
            if b"main.tsx" not in cmdline and b"/dist/delfin-tui-" not in cmdline:
                continue
            environ = (entry / "environ").read_bytes().split(b"\0")
            if f"DELFIN_RUNTIME_DIR={runtime}".encode() in environ:
                return int(entry.name)
        except (OSError, IndexError, ValueError):
            continue
    return 0


def _screen_wait(tmux: Tmux, pattern: str, timeout: float):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        m = re.search(pattern, tmux.screen())
        if m:
            return m
        time.sleep(0.1)
    return None


def _screen_gone(tmux: Tmux, pattern: str, timeout: float) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not re.search(pattern, tmux.screen()):
            return True
        time.sleep(0.1)
    return False


if __name__ == "__main__":
    sys.exit(main())
