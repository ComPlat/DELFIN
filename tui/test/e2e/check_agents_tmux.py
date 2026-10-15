"""The agent view, live: `delfin-agent tui` in tmux against a private agent
service in which two more sessions run (started from outside, as other
windows would): one done, one waiting on a question.

The hints say a session waits; ← opens the agent view over the terminal
(its alternate screen: esc gives the screen back exactly as it was),
grouped; space peeks at a session and a reply goes to it without opening
it; enter opens the waiting one (divider, its question) and y answers it;
← again, typing and enter start a new session in this window's folder,
which takes the prompt. Leaving the window leaves all four running, and
no row of the scrollback is lost or written twice by the view opening and
closing over it.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_agents_tmux.py
"""

from __future__ import annotations

import os
import re
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
    tmp = Path(tempfile.mkdtemp(prefix="tui-agents-"))
    runtime = Path(tempfile.mkdtemp(prefix="dta", dir="/tmp"))
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
        found = (tmux.wait_for(pattern, timeout) if where == "history"
                 else _screen_wait(tmux, pattern, timeout, tmux.history if where == "scrollback" else tmux.screen))
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
        with ServiceClient(sock) as other:
            b = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                                                 "name": "second-session"}, timeout=60)["sessionId"]
            mark = len(other.events)
            other.request("session.prompt", {"sessionId": b, "text": "hello from B"})
            other.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == b,
                           after=mark, timeout=30)
            c = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                                                 "name": "third-asks", "permissionMode": "default"},
                              timeout=60)["sessionId"]
            mark = len(other.events)
            other.request("session.prompt", {"sessionId": c, "text": "please edit the parser"})
            other.wait_for(lambda e: e.get("type") == "permission.asked" and e.get("session") == c,
                           after=mark, timeout=30)
        # At the right end of the hints line: first the session that asks
        # names itself for a few seconds, then the count takes its place.
        expect(r"◉ third-asks asks|agents\s+◉ 1\b", "the note at the hints' right", 15, where="screen")
        expect(r"agents\s+◉ 1\b", "the count at the hints' right", 15, where="screen")
        # Opened and closed again, the view leaves the screen as it was.
        time.sleep(0.5)
        before = tmux.screen()
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view", 15, where="screen")
        tmux.key("Escape")
        time.sleep(1.0)
        if tmux.screen() != before:
            problems.append(f"closing the view changed the screen:\n--- before\n{before}\n--- after\n{tmux.screen()}")
        # ← opens the view over the terminal.
        tmux.key("Left")
        expect(r"◆ AGENTS\s+3 sessions · 1 needs you", "the agent view", 15, where="screen")
        screen = tmux.screen()
        for title in ("NEEDS YOU", "third-asks", "DONE", "second-session"):
            if title not in screen:
                problems.append(f"the agent view lacks {title!r}\n{screen}")
        if screen.index("NEEDS YOU") > screen.index("DONE"):
            problems.append("NEEDS YOU should come before DONE")
        # Peek at B (one up from this window's row) and reply to it.
        tmux.key("Up")
        time.sleep(0.3)
        tmux.key("Space")
        expect(r"│ ❯ hello from B", "the peek", 15, where="screen")
        tmux.type("a reply to B")
        tmux.key("Enter")
        expect(r"sent to second-session|queued for second-session", "the reply", 15, where="screen")
        tmux.key("Escape")                     # the peek
        time.sleep(0.4)
        # Open the waiting one: its divider and its question, answered with y.
        tmux.key("Up")
        time.sleep(0.3)
        tmux.key("Enter")
        # Written into the scrollback, ahead of C's banner and replay: on a
        # 40-row pane the colour logo's 15 rows and the open diff can
        # push it above the visible screen.
        expect(r"now: third-asks", "the divider", 15, where="scrollback")
        expect(r"please edit the parser", "C's replay", 15, where="scrollback")
        expect(r"\byes\b", "C's question", 15, where="screen")
        time.sleep(0.6)                        # the pause before a key answers
        tmux.type("y")
        expect(r"approved: edit_file", "C's answer", 30)
        expect(r"The parser now reads the energy", "C's turn after the answer", 30)
        # ← again: a new session from what is typed.
        time.sleep(0.5)
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view again", 15, where="screen")
        tmux.type("a brand new job")
        tmux.key("Enter")
        expect(r"a new session", "the new session's divider", 30)
        expect(r"\(re: a brand new job\)", "the new session's turn", 30)
        time.sleep(0.5)
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        history = tmux.history()
        for line in ("❯ hello from A", "(re: hello from A)"):
            n = history.count(line)
            if n != 1:
                problems.append(f"{line!r} is in the scrollback {n} times, not once")
        with ServiceClient(sock) as cl:
            sessions = {s["sessionId"]: s for s in cl.request("sessions.list")["sessions"]}
            if len(sessions) != 4:
                problems.append(f"four sessions should run on: {[s.get('title') for s in sessions.values()]}")
            said = [e["text"] for e in cl.request("session.replay", {"sessionId": b})["events"]
                    if e.get("type") == "user.message"]
            if said != ["hello from B", "a reply to B"]:
                problems.append(f"B got: {said}")
            new = [s for s in sessions.values() if s["sessionId"] not in (b, c)
                   and s.get("lastPrompt") == "a brand new job"]
            if not new or new[0].get("workspace") != str(ws):
                problems.append(f"the new session should run in {ws}: {new}")
            cl.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux agent view check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the agent view lists, peeks, replies, opens and starts sessions, live")
    return 0


def _screen_wait(tmux: Tmux, pattern: str, timeout: float, read=None):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        m = re.search(pattern, (read or tmux.screen)())
        if m:
            return m
        time.sleep(0.1)
    return None


if __name__ == "__main__":
    sys.exit(main())
