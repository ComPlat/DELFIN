"""Clicks, live: `delfin-agent tui --fake` full screen in tmux against a
private agent service. The mouse reports are sent as the terminal sends
them (SGR, ESC [ < b ; x ; y M/m) into the pane.

- a click on "[y] yes" of a shell approval answers it: the command runs;
- a click on a question's choice answers it with that choice;
- a click that arrives as the panel appears (within the pause) does not.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_mouse_tmux.py
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
    tmp = Path(tempfile.mkdtemp(prefix="tui-mouse-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtm", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO), "TERM": "xterm-256color",
           "LANG": "C.UTF-8", "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    problems: list[str] = []
    tmux = Tmux(100, 34)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"

    def find(text: str, timeout: float = 15) -> tuple[int, int] | None:
        until = time.monotonic() + timeout
        while time.monotonic() < until:
            rows = tmux.screen().split("\n")
            for y, row in enumerate(rows):
                if text in row:
                    return row.index(text), y
            time.sleep(0.1)
        return None

    def click(x: int, y: int) -> None:
        # 1-based cells; press then release on the same cell.
        tmux.type(f"\x1b[<0;{x + 1};{y + 1}M")
        time.sleep(0.05)
        tmux.type(f"\x1b[<0;{x + 1};{y + 1}m")

    try:
        # A clean environment: nothing of the shell's (credentials) reaches the session.
        tmux.type(f"clear; cd {REPO} && env -i {' '.join(f'{k}={v}' for k, v in env.items())} "
                  f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
        tmux.key("Enter")
        if not find("ask DELFIN", 60):
            problems.append(f"the window did not come up\n{tmux.screen()}")
            raise SystemExit
        tmux.type("run the bash tests")
        tmux.key("Enter")
        at = find("[y] yes")
        if not at:
            problems.append(f"no approval panel\n{tmux.screen()}")
            raise SystemExit
        click(at[0] + 1, at[1])                 # as it appears: within the pause
        time.sleep(0.6)
        if "All three tests pass." in tmux.screen():
            problems.append("a click as the panel appeared answered it")
        at = find("[y] yes", 2)
        if at:
            click(at[0] + 1, at[1])
        if not find("All three tests pass.", 15):
            problems.append(f"a click on [y] yes did not answer it\n{tmux.screen()}")
        tmux.type("ask about the solvent")
        tmux.key("Enter")
        at = find("density-based")
        if not at:
            problems.append(f"no question panel\n{tmux.screen()}")
            raise SystemExit
        time.sleep(0.6)
        click(at[0], at[1])
        if not find("Using SMD.", 15):
            problems.append(f"a click on a choice did not answer with it\n{tmux.screen()}")
        tmux.key("C-d")
        if not find("TUI_EXIT=0", 30):
            problems.append("the window did not leave cleanly")
        with ServiceClient(sock) as c:
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        time.sleep(1)
        shutil.rmtree(runtime, ignore_errors=True)
        shutil.rmtree(tmp / "home", ignore_errors=True)
    print(f"tmux mouse check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: a click answers an approval and a question, and not as the panel appears")
    return 0


if __name__ == "__main__":
    sys.exit(main())
