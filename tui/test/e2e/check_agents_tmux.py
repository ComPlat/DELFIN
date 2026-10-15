"""The agent view, live: `delfin-agent tui` in tmux against a private agent
service in which two more sessions run (started from outside, as other
windows would): one done, one waiting on a question.

The waiting session says so once in the message row, and the hints line
counts the other two ("2 more sessions (1 waits for you) ← agents");
`delfin-agent service status` lists all three, with the window on this
one. ← opens the agent view over the terminal (its alternate screen: esc
gives the screen back exactly as it was), grouped, its keys named in its
head; space peeks at a session and a reply goes to it without opening it;
enter opens the waiting one (divider, its question) and y answers it; ←
again, typing and enter start a new session in this window's folder, which
takes the prompt. `delfin-agent service end <prefix> --yes` ends one of the
others and the count in the hints line follows; `--after-turn` on a session
in a turn of some seconds waits for that turn, which finishes (not stopped). Leaving the window leaves
the other three running, and no row of the scrollback is lost or written
twice by the view opening and closing over it.

DELFIN_CHECK_CAPTURES=<dir>: the screens and command outputs checked, as
text files there.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_agents_tmux.py
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[2]
REPO = TUI.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TUI / "spike"))

from check_tmux import Tmux, scrubbed_env  # noqa: E402

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
    captures = Path(os.environ["DELFIN_CHECK_CAPTURES"]) if os.environ.get("DELFIN_CHECK_CAPTURES") else None

    def capture(name: str, text: str) -> None:
        if captures is not None:
            captures.mkdir(parents=True, exist_ok=True)
            (captures / f"{name}.txt").write_text(text)

    def service(*args: str) -> subprocess.CompletedProcess:
        """`delfin-agent service ...` as the user types it, against this
        check's service (no credentials in its environment)."""
        return subprocess.run([PYTHON, "-m", "delfin.agent.cli", "service", *args],
                              env=scrubbed_env(**env), capture_output=True, text=True, timeout=120,
                              cwd=str(REPO))

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
            # What the window finds in its first 2.5 s is no news (ARRIVAL_QUIET_MS).
            time.sleep(3.0)
            c = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                                                 "name": "third-asks", "permissionMode": "default"},
                              timeout=60)["sessionId"]
            mark = len(other.events)
            other.request("session.prompt", {"sessionId": c, "text": "please edit the parser"})
            other.wait_for(lambda e: e.get("type") == "permission.asked" and e.get("session") == c,
                           after=mark, timeout=30)
        # The session that asks says so once in the message row; the hints
        # line counts the other two, and at its right what waits.
        if expect(r"◉ session 'third-asks' waits for you — ← to open", "the note of a new wait", 15,
                  where="screen"):
            capture("1-notice-and-count", tmux.screen())
        expect(r"· 2 more sessions \(1 waits for you\) ← agents", "the count of the others", 15, where="screen")
        expect(r"agents\s+◉ 1\b", "the count at the hints' right", 15, where="screen")
        # `service status` lists the three, the window on this one.
        status = service("status")
        capture("2-service-status", status.stdout + status.stderr)
        rows = {line.split()[0]: line for line in status.stdout.splitlines()[2:] if line.startswith("  ")
                and not line.strip().startswith(("`", "sleeping:"))}
        if status.returncode != 0 or not re.search(r"3 sessions, 1 window", status.stdout):
            problems.append(f"service status: {status.returncode}\n{status.stdout}{status.stderr}")
        if "second-session" not in rows.get(b[:8], "") or "third-asks" not in rows.get(c[:8], ""):
            problems.append(f"service status should list B and C by id:\n{status.stdout}")
        if "waiting for you" not in rows.get(c[:8], ""):
            problems.append(f"C should be waiting for you:\n{status.stdout}")
        listed = json.loads(service("status", "--json").stdout or "{}")
        windows = {r["id"]: r["windows"] for r in listed.get("sessions", [])}
        if sorted(windows.values()) != [0, 0, 1]:
            problems.append(f"one window, on this window's session: {windows}")
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
        if expect(r"◆ AGENTS.*→/enter open · ctrl\+x ×2 stop · esc back", "the view's keys in its head", 15,
                  where="screen"):
            capture("3-agent-view", tmux.screen())
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
        # The other three, counted; `service end` ends B and the count follows.
        expect(r"· 3 more sessions ← agents", "the count of the others", 15, where="screen")
        with ServiceClient(sock) as cl:
            said = [e["text"] for e in cl.request("session.replay", {"sessionId": b})["events"]
                    if e.get("type") == "user.message"]
        if said != ["hello from B", "a reply to B"]:
            problems.append(f"B got: {said}")
        ended = service("end", b[:6], "--yes")
        capture("4-service-end", ended.stdout + ended.stderr)
        # A scripted session started without fakeSave is not saved, and is not said to be.
        if ended.returncode != 0 or "second-session" not in ended.stdout or "not saved" not in ended.stdout:
            problems.append(f"service end: {ended.returncode}\n{ended.stdout}{ended.stderr}")
        if expect(r"· 2 more sessions ← agents", "the count after an end", 15, where="screen"):
            capture("5-count-after-end", tmux.screen())
        # --after-turn: a session in a turn of some seconds ends after it, not in it.
        with ServiceClient(sock) as other:
            d = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.12,
                                                 "name": "fourth-slow"}, timeout=60)["sessionId"]
            mark = len(other.events)
            other.request("session.prompt", {"sessionId": d, "text": "a slow one"})
            other.wait_for(lambda e: e.get("type") == "turn.started" and e.get("session") == d,
                           after=mark, timeout=30)
            expect(r"· 3 more sessions ← agents", "the count with the slow one", 15, where="screen")
            after = service("end", d[:8], "--after-turn", "--yes")
            capture("6-service-end-after-turn", after.stdout + after.stderr)
            if after.returncode != 0 or "waiting for" not in after.stdout or "after its turn" not in after.stdout:
                problems.append(f"service end --after-turn: {after.returncode}\n{after.stdout}{after.stderr}")
            turn_end = next((e for e in other.events[mark:] if e.get("type") == "turn.ended"
                             and e.get("session") == d), None)
            if not turn_end or turn_end.get("stopped"):
                problems.append(f"--after-turn should let the turn finish: {turn_end}")
        expect(r"· 2 more sessions ← agents", "the count after --after-turn", 15, where="screen")
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
            if len(sessions) != 3 or b in sessions:
                problems.append(f"three sessions should run on, B ended: "
                                f"{[s.get('title') for s in sessions.values()]}")
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
