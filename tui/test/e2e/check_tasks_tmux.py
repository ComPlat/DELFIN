"""The task list during a run, live: `delfin-agent tui --fake` in tmux
against a private agent service, in both screens.

The fake engine's "tasks" turn writes three tasks into the workspace task
store (as the task tools do) and works through them over a few seconds.
While the turn runs, the list must be on screen above the prompt with the
task in hand marked (◼, under its active form), and move on as the store
does: each of the three is seen in hand, the earlier ones done (☒).

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_tasks_tmux.py
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
IN_HAND = ("reading the input files", "fixing the energy parser", "running the tests")


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-tasks-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtt", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    problems: list[str] = []
    tmux = Tmux(100, 34)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def run(screen: str) -> None:
        extra = "--inline" if screen == "inline" else ""
        tmux.type(f"clear; cd {REPO} && env -u DELFIN_TUI_SCREEN -u STY "
                  f"{' '.join(f'{k}={v}' for k, v in env.items())} "
                  f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws} {extra}; echo TUI_EXIT=$?")
        tmux.key("Enter")
        if not tmux.wait_for(r"ask DELFIN", 60):
            problems.append(f"{screen}: the window did not come up\n{tmux.screen()}")
            return
        tmux.type("work through the tasks")
        tmux.key("Enter")
        seen: list[str] = []
        done_seen = False
        captures: list[str] = []
        until = time.monotonic() + 30
        while time.monotonic() < until:
            rows = tmux.screen().split("\n")
            prompt = max((i for i, r in enumerate(rows) if "❯ " in r), default=len(rows))
            above = rows[max(0, prompt - 6):prompt]
            for label in IN_HAND:
                if any(f"◼ {label}" in r for r in above) and label not in seen:
                    seen.append(label)
                    captures.append("\n".join(above))
            done_seen = done_seen or any("☒ read the input files" in r for r in above)
            if "All three tasks are done." in "\n".join(rows):
                break
            time.sleep(0.1)
        if seen != list(IN_HAND):
            problems.append(f"{screen}: the tasks in hand seen above the prompt were {seen}\n{tmux.screen()}")
        if not done_seen:
            problems.append(f"{screen}: a finished task was never shown done (☒)")
        if captures:
            print(f"--- {screen}, while the turn ran:\n{captures[-1]}")
        tmux.key("C-d")
        if not tmux.wait_for(r"TUI_EXIT=0", 30):
            problems.append(f"{screen}: the window did not leave cleanly")

    try:
        run("full")
        run("inline")
        with ServiceClient(sock) as c:
            c.request("service.shutdown")
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        time.sleep(1)
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux task-list check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the task list is drawn above the prompt while the turn works through it, in both screens")
    return 0


if __name__ == "__main__":
    sys.exit(main())
