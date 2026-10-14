"""A resize storm, live in tmux: the window of a session with a dozen turns
is resized thirty times in three seconds (the way a tiling window manager
or a dragged corner does it), ends at another width, and then

- writes the session again once, at the final width: every turn of the
  replay is in the scrollback exactly once;
- shows its footer once, at the bottom;
- and goes on: a prompt after the storm is answered, once.

Also with the agent view open during the storm: closing it gives a screen
drawn at the final width.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_resize_tmux.py
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
TURNS = 12
WIDTHS = [80, 130, 90, 140, 100, 70, 120, 85, 135, 95] * 3


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-resize-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtr", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(110, 40)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def expect(pattern: str, what: str, timeout: float = 15) -> bool:
        if tmux.wait_for(pattern, timeout):
            return True
        problems.append(f"{what}: no {pattern!r}\n{tmux.screen()}")
        return False

    def storm(final: int) -> None:
        for width in WIDTHS:
            tmux.run("resize-window", "-t", tmux.name, "-x", str(width), "-y", "40")
            time.sleep(0.1)
        tmux.run("resize-window", "-t", tmux.name, "-x", str(final), "-y", "40")

    def once_each(history: str, when: str) -> None:
        for n in range(1, TURNS + 1):
            count = len(re.findall(rf"❯ turn {n}\b", history))
            if count != 1:
                problems.append(f"{when}: '❯ turn {n}' is in the scrollback {count} times")
        footers = history.count("ask DELFIN")
        if footers != 1:
            problems.append(f"{when}: the footer is in the scrollback {footers} times")

    try:
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        for n in range(1, TURNS + 1):
            tmux.type(f"turn {n}")
            tmux.key("Enter")
            if not expect(rf"\(re: turn {n}\)", f"turn {n}", 30):
                raise SystemExit
        time.sleep(0.5)
        storm(100)
        time.sleep(3.0)                         # the replay after it settles
        once_each(tmux.history(), "after the storm")
        width = max(len(line) for line in tmux.screen().splitlines())
        if width > 100:
            problems.append(f"after the storm a row is {width} wide at 100 columns")
        tmux.type("after the storm")
        tmux.key("Enter")
        expect(r"\(re: after the storm\)", "a prompt after the storm", 30)
        time.sleep(0.5)
        if tmux.history().count("❯ after the storm") != 1:
            problems.append("the prompt after the storm is not in the scrollback once")
        # A storm while the agent view is open.
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view", 15)
        storm(90)
        time.sleep(1.0)
        tmux.key("Escape")
        time.sleep(3.0)
        screen = tmux.screen()
        if "◆ AGENTS" in screen:
            problems.append(f"the agent view is still drawn after it closed:\n{screen}")
        width = max(len(line) for line in screen.splitlines())
        if width > 90:
            problems.append(f"after the view closed a row is {width} wide at 90 columns")
        once_each(tmux.history(), "after the storm with the view open")
        tmux.key("C-d")
        code = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not code or code.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {code and code.group(1)}")
        with ServiceClient(sock) as c:
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux resize storm check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: a resize storm ends in one replay at the final width, also with the agent view open")
    return 0


if __name__ == "__main__":
    sys.exit(main())
