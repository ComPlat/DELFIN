"""Switching sessions, live: `delfin-agent tui` in tmux against a private
agent service in which a second session runs (started from outside, as
another window would). ctrl+p → resume… → that session: the window shows
a divider, the session's last turns, and from then on its live events; a
prompt goes to it; the first session keeps running in the service.

Also: a bare /model opens the model picker, and esc closes it.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_switch_tmux.py
"""

from __future__ import annotations

import os
import shutil
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
    tmp = Path(tempfile.mkdtemp(prefix="tui-switch-"))
    runtime = Path(tempfile.mkdtemp(prefix="dts", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO), "DELFIN_TUI_SCREEN": "inline",
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(100, 40)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def expect(pattern: str, what: str, timeout: float = 15, where: str = "history") -> bool:
        found = tmux.wait_for(pattern, timeout) if where == "history" else _screen_wait(tmux, pattern, timeout)
        if found:
            return True
        problems.append(f"{what}: no {pattern!r}\n{tmux.screen()}")
        return False

    try:
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        tmux.type("hello from A")
        tmux.key("Enter")
        expect(r"\(re: hello from A\)", "A's turn", 30)
        # A second session, from outside: as another window would start it.
        with ServiceClient(sock) as other:
            b = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                                                 "name": "second-session"}, timeout=60)["sessionId"]
            mark = len(other.events)
            other.request("session.prompt", {"sessionId": b, "text": "hello from B"})
            other.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == b,
                           after=mark, timeout=30)
        # A bare /model opens the picker; esc closes it.
        tmux.type("/model")
        time.sleep(0.4)
        tmux.key("Escape")                      # the completion list
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"fake-model.*\(now\)", "the model picker", 15, where="screen")
        tmux.key("Escape")
        time.sleep(0.5)
        # ctrl+p → resume… → the second session.
        tmux.key("C-p")
        expect(r"resume…", "the palette", 10, where="screen")
        tmux.type("resume")
        tmux.key("Enter")
        expect(r"second-session", "the resume picker", 15, where="screen")
        tmux.type("second")
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"now: second-session", "the divider", 15)
        expect(r"❯ hello from B", "B's replay", 15)
        tmux.type("and again")
        tmux.key("Enter")
        expect(r"\(re: and again\)", "a prompt to B", 30)
        time.sleep(0.5)
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        with ServiceClient(sock) as c:
            sessions = {s["sessionId"]: s for s in c.request("sessions.list")["sessions"]}
            if len(sessions) != 2:
                problems.append(f"both sessions should run on: {list(sessions)}")
            events = c.request("session.replay", {"sessionId": b})["events"]
            said = [e["text"] for e in events if e.get("type") == "user.message"]
            if said != ["hello from B", "and again"]:
                problems.append(f"B got: {said}")
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux switch check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the resume picker switches the window to another session, live")
    return 0


def _screen_wait(tmux: Tmux, pattern: str, timeout: float):
    import re
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        m = re.search(pattern, tmux.screen())
        if m:
            return m
        time.sleep(0.1)
    return None


if __name__ == "__main__":
    sys.exit(main())
