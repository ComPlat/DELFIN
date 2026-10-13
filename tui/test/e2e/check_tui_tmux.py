"""The window, live: `delfin-agent tui` in tmux against a private agent
service, typed into like a person would, must show what the snapshot test
renders from the recorded conversation (test/fixtures/).

Then it leaves (ctrl+d), and the session must still be running in the
service: closing the window never ends a session.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_tui_tmux.py
"""

from __future__ import annotations

import json
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
WORKSPACE = "/home/user/project"


def normalise(text: str) -> str:
    """The snapshot test's normalise(), line for line."""
    lines = [re.sub(r"\s+$", "", line) for line in text.split("\n")]
    lines = [re.sub(r"[✓✗] \d+(\.\d+)?s", lambda m: m.group(0)[0] + " N.Ns", line) for line in lines]
    lines = [re.sub(r"↑[\d.]+[kM]? ↓[\d.]+[kM]?", "↑N ↓N", line) for line in lines]
    return re.sub(r"\n{3,}", "\n\n", "\n".join(lines)).strip()


def main() -> int:
    fixture = json.loads((TUI / "test" / "fixtures" / "fake-conversation.json").read_text())
    expected = (TUI / "test" / "fixtures" / "fake-conversation.expected.txt").read_text().strip()
    tmp = Path(tempfile.mkdtemp(prefix="tui-e2e-"))
    runtime = Path(tempfile.mkdtemp(prefix="dte", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --permission-mode default "
               f"--cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(100, 35)
    try:
        tmux.type(command)
        tmux.key("Enter")
        if not tmux.wait_for(r"ask DELFIN", 60):
            problems.append("the window did not come up")
            raise SystemExit
        turns = 0
        for step in fixture["script"]:
            tmux.type(step["prompt"])
            tmux.key("Enter")
            reply = step.get("reply") or {}
            # A person reads the question before answering it; a key in the
            # same instant it opens is typeahead and not taken (approval.tsx).
            if "choice" in reply:
                if not tmux.wait_for(r"\[y\] ", 30):
                    problems.append(f"no approval for {step['prompt']!r}")
                time.sleep(0.6)
                tmux.key(reply["choice"])
            elif "answers" in reply:
                found = tmux.wait_for(r"(\d)\. " + re.escape(reply["answers"][0]), 30)
                if not found:
                    problems.append(f"no question for {step['prompt']!r}")
                else:
                    time.sleep(0.6)
                    tmux.key(found.group(1))
            turns += 1
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline:
                if len(re.findall(r"✓ \d+\.\ds · ↑", tmux.history())) >= turns:
                    break
                time.sleep(0.2)
            else:
                problems.append(f"turn {turns} did not end")
        time.sleep(0.5)
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        history = tmux.history().replace(str(ws), WORKSPACE)
        (tmp / "history.txt").write_text(history)
        start = history.find("◆ DELFIN")
        ends = [m.end() for m in re.finditer(r"[✓✗] \d+\.\ds · ↑\S+ ↓\S+", history)]
        shown = normalise(history[start:ends[-1]]) if start >= 0 and ends else ""
        if shown != expected:
            import difflib
            diff = "\n".join(difflib.unified_diff(expected.split("\n"), shown.split("\n"),
                                                  "snapshot", "live", lineterm="", n=1))
            problems.append(f"the live window differs from the snapshot:\n{diff}")
        sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"
        with ServiceClient(sock) as c:
            running = c.request("sessions.list")["sessions"]
            if len(running) != 1 or running[0]["state"] != "idle":
                problems.append(f"after the window closed: {running}")
            c.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux live window vs snapshot; history in {tmp / 'history.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: the live window shows the snapshot, and the session runs on")
    return 0


if __name__ == "__main__":
    sys.exit(main())
