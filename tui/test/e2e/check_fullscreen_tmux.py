"""The window full screen (the default), live in tmux:

- it takes the whole terminal: the shell's lines are gone while it runs,
  the composer is at the bottom, the conversation above it;
- a tool call is its headline and what came back; ctrl+o shows its
  details in place and hides them again;
- PgUp scrolls back to the session's banner, ctrl+end to its end;
- the agent view takes the same screen and gives it back;
- after a resize every row fits the new width and the composer is still
  at the bottom;
- leaving gives the shell its lines back as they were, says the session
  keeps running, and it does, and prints the command that comes back to it;
- `--inline` prints that command too, and running it opens that session;
- `--inline` draws under the shell's lines instead, as before.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_fullscreen_tmux.py
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
MARK = "before-the-window-was-opened"


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-full-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtf", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    (ws / "README.md").write_text("# tui_demo\n")
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    for name in ("DELFIN_TUI_SCREEN", "STY"):
        env.pop(name, None)

    def command(extra: str = "") -> str:
        # Neither a screen setting nor GNU screen's $STY (which picks inline)
        # may come from the shell the check runs in.
        return (f"cd {REPO} && env -u DELFIN_TUI_SCREEN -u STY {' '.join(f'{k}={v}' for k, v in env.items())} "
                f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws} -n 'scan of benzene' {extra}; echo TUI_EXIT=$?")

    problems: list[str] = []
    tmux = Tmux(100, 30)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def expect(pattern: str, what: str, timeout: float = 15) -> bool:
        if tmux.wait_for(pattern, timeout):
            return True
        problems.append(f"{what}: no {pattern!r}\n{tmux.screen()}")
        return False

    def rows() -> list[str]:
        return tmux.screen().rstrip("\n").split("\n")

    try:
        tmux.type(f"clear; echo {MARK}")
        tmux.key("Enter")
        time.sleep(0.5)
        tmux.type(command())
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window", 60):
            raise SystemExit
        screen = tmux.screen()
        if MARK in screen:
            problems.append("the shell's lines are still on the screen: the window did not take it")
        composer = [i for i, line in enumerate(rows()) if line.lstrip().startswith("❯ ask DELFIN")]
        if not composer or composer[-1] < len(rows()) - 4:
            problems.append(f"the composer is not at the bottom: rows {composer} of {len(rows())}")
        for n in range(4):
            tmux.type(f"turn number {n}")
            tmux.key("Enter")
            expect(rf"\(re: turn number {n}\)", f"turn {n}", 30)
            time.sleep(0.4)
        if "(re: turn number 3)" not in tmux.screen():
            problems.append("the latest answer is not in view")
        # A tool call is its headline and what came back; ctrl+o shows what
        # it printed, in place, and again hides it.
        if "# tui_demo" in tmux.screen():
            problems.append("a tool call's details were shown before ctrl+o")
        tmux.key("C-o")
        if not tmux.wait_for(r"# tui_demo", 10):
            problems.append(f"ctrl+o did not show the tool calls' details:\n{tmux.screen()}")
        tmux.key("C-o")
        time.sleep(0.8)
        if "# tui_demo" in tmux.screen():
            problems.append("ctrl+o again did not hide them")
        # Back to the start, and to the end again.
        tmux.key("C-Home")
        time.sleep(0.6)
        if "DELFIN  fake/fake-model" not in tmux.screen():
            problems.append(f"ctrl+home did not show the banner:\n{tmux.screen()}")
        tmux.key("C-End")
        time.sleep(0.6)
        tmux.key("PageUp")
        time.sleep(0.6)
        if "(re: turn number 3)" in tmux.screen():
            problems.append("PgUp did not scroll back")
        tmux.key("C-End")
        time.sleep(0.6)
        if "(re: turn number 3)" not in tmux.screen():
            problems.append("ctrl+end did not come back to the end")
        # The agent view in the same screen.
        tmux.key("Left")
        if expect(r"◆ AGENTS", "the agent view", 10):
            if "scan of benzene" not in tmux.screen():
                problems.append("the agent view does not list the session")
            tmux.key("Escape")
            time.sleep(1.0)
            if "◆ AGENTS" in tmux.screen() or "(re: turn number 3)" not in tmux.screen():
                problems.append(f"esc did not give the conversation back:\n{tmux.screen()}")
        # A resize: every row fits, the composer stays at the bottom.
        tmux.run("resize-window", "-t", tmux.name, "-x", "70", "-y", "24")
        time.sleep(1.5)
        too_wide = [line for line in rows() if len(line.rstrip()) > 70]
        if too_wide:
            problems.append(f"after the resize, rows wider than 70: {too_wide[:3]}")
        composer = [i for i, line in enumerate(rows()) if line.lstrip().startswith("❯ ")]
        if not composer or composer[-1] < len(rows()) - 4:
            problems.append(f"after the resize the composer is not at the bottom: {composer} of {len(rows())}")
        # Leaving (at the first size, so the shell's lines fit): they come
        # back, and a line says the session runs on.
        tmux.run("resize-window", "-t", tmux.name, "-x", "100", "-y", "30")
        time.sleep(1.0)
        tmux.key("C-d")
        if expect(r"TUI_EXIT=(\d+)", "leaving", 30):
            after = tmux.screen()
            if MARK not in after:
                problems.append(f"the shell's lines did not come back:\n{after}")
            if "keeps running in the agent service" not in after:
                problems.append("no line said the session keeps running")
            if "TUI_EXIT=0" not in after:
                problems.append("the window did not leave cleanly")
            if not re.search(r"resume with: delfin-agent tui --resume [0-9A-Za-z_.:-]+", after):
                problems.append(f"no line said how to come back:\n{after}")
        with ServiceClient(sock) as c:
            if not any(s.get("title") == "scan of benzene" for s in c.request("sessions.list")["sessions"]):
                problems.append("the session did not keep running")
        # --inline: under the shell's lines, as before.
        tmux.type(f"clear; echo {MARK}")
        tmux.key("Enter")
        time.sleep(0.5)
        tmux.type(command("--inline"))
        tmux.key("Enter")
        if expect(r"ask DELFIN", "the inline window", 60):
            # Under the shell's lines: not the alternate screen, and the line
            # printed before is still there (scrolled up by the long command
            # line and the banner on a short pane, so in the history).
            alternate = tmux.run("display-message", "-p", "-t", tmux.name, "#{alternate_on}").strip()
            if alternate == "1" or MARK not in tmux.history():
                problems.append(f"--inline: the shell's lines are gone (it took the screen):\n{tmux.screen()}")
            tmux.key("C-d")
            if expect(r"TUI_EXIT=0", "leaving the inline window", 30):
                # The printed command, typed as printed, opens that session.
                back = re.findall(r"resume with: (delfin-agent tui --resume (\S+))", tmux.history())
                if not back:
                    problems.append(f"--inline: no line said how to come back:\n{tmux.screen()}")
                else:
                    sid = back[-1][1]
                    tmux.type("clear; " + command(f"--resume {sid}").replace(f"tui --fake --cwd {ws} -n 'scan of benzene'",
                                                                   f"tui --fake --cwd {ws}"))
                    tmux.key("Enter")
                    until = time.time() + 60
                    while time.time() < until and not ("ask DELFIN" in tmux.screen()
                                                       and " scan of benzene " in tmux.screen()):
                        time.sleep(0.3)
                    if " scan of benzene " not in tmux.screen():
                        problems.append(f"the printed command did not open that session:\n{tmux.screen()}")
                    tmux.key("C-d")
                    time.sleep(2)
        with ServiceClient(sock) as c:
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        time.sleep(1)
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux full-screen check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the window takes the screen, scrolls inside it, gives the shell back, and --inline does not")
    return 0


if __name__ == "__main__":
    sys.exit(main())
