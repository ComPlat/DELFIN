"""The composer, live: `delfin-agent tui` in tmux against a private agent
service, typed into like a person would.

  /st             the command list opens, from the session's own commands
  /dep            a workspace command (.delfin/commands/deploy.md) is offered
  @del + tab      paths complete from the workspace, folder by folder
  alt+enter       a prompt of two lines goes out as one message
  a large paste   stands in as [Pasted text #1 +30 lines], arrives whole
  up              the history brings the last prompt back

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_composer_tmux.py
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
    tmp = Path(tempfile.mkdtemp(prefix="tui-composer-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtc", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws / "delfin", ws / "docs",
              ws / ".delfin" / "commands"):
        d.mkdir(parents=True)
    (ws / "delfin" / "xtb.py").write_text("")
    (ws / "docs" / "guide.md").write_text("")
    (ws / ".delfin" / "commands" / "deploy.md").write_text(
        "---\ndescription: Ship the build\n---\nDeploy $ARGUMENTS\n")
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO), "DELFIN_TUI_SCREEN": "inline",
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(100, 35)

    def expect(pattern: str, what: str, timeout: float = 15) -> bool:
        if tmux.wait_for(pattern, timeout):
            return True
        problems.append(f"{what}: no {pattern!r} on screen\n{tmux.screen()}")
        return False

    def clear() -> None:
        tmux.key("C-u")
        time.sleep(0.3)

    try:
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        tmux.type("/st")
        expect(r"/status\s+Model, mode, tokens, cost", "the command list")
        clear()
        tmux.type("/dep")
        expect(r"/deploy\s+Ship the build", "a workspace command")
        clear()
        tmux.type("look at @del")
        expect(r"@delfin/", "the path list")
        tmux.key("Tab")
        expect(r"@delfin/xtb\.py", "the list inside the folder")
        tmux.key("Tab")
        expect(r"❯ look at @delfin/xtb\.py", "the completed path")
        clear()
        # Two lines, one message.
        tmux.type("first line")
        tmux.key("M-Enter")
        tmux.type("second line")
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"\(re: first line", "the two-line turn", 30)
        # A paste the size of a log.
        log = "\n".join(f"log line {i:02d}" for i in range(30))
        tmux.run("set-buffer", "-b", "log", log)
        tmux.run("paste-buffer", "-p", "-b", "log", "-t", tmux.name)
        expect(r"\[Pasted text #1 \+30 lines\]", "the paste placeholder")
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"❯ \[Pasted text #1 \+30 lines\]", "the placeholder in the scrollback", 30)
        expect(r"\(re: log line 00", "the pasted turn", 30)
        time.sleep(1.0)
        # The history: the paste comes back as its placeholder.
        tmux.key("Up")
        expect(r"❯ \[Pasted text #1 \+30 lines\]\s*$|❯ \[Pasted text #1 \+30 lines\]",
               "the history")
        tmux.key("C-c")                     # clears what is typed
        time.sleep(0.3)
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"
        with ServiceClient(sock) as c:
            sessions = c.request("sessions.list")["sessions"]
            sid = sessions[0]["sessionId"]
            events = c.request("session.replay", {"sessionId": sid})["events"]
            said = [e["text"] for e in events if e.get("type") == "user.message"]
            if "first line\nsecond line" not in said:
                problems.append(f"the two-line prompt did not arrive as one: {said}")
            if log not in said:
                problems.append("the paste did not arrive whole")
            c.request("service.shutdown")
        history = tmp / "scratch" / ".delfin" / "tui_history.jsonl"
        if not history.exists() or "Pasted text #1" not in history.read_text():
            problems.append("the history file does not hold the prompts")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux composer check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: completion, new lines, pastes and the history work live")
    return 0


if __name__ == "__main__":
    sys.exit(main())
