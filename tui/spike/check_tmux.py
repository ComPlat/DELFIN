"""Run the spike inside tmux and check the scrollback it leaves.

tmux is a terminal emulator with its own history, so what it keeps after the
run is what a user could scroll back to. Each run uses a private tmux server
(``tmux -L``, no config file), so the user's own tmux sessions and settings
are never touched.

Modes:

  default        the spike's --auto script; the history must hold every
                 row the spike committed, exactly once and in order
                 (scrollback_check.py), and the spike must exit 0 with no
                 errors and nothing on stderr;
  --prefill N    print N shell lines first, so the spike starts N rows down
                 the screen; those lines must survive above its output;
  --resize-to C  resize the pane to C columns mid-run: the spike must replay
                 exactly once, at the new width;
  --ssh-key K    run the spike through ``ssh -tt localhost`` (a real hop:
                 its own pty, window-size forwarding, a reduced environment);
  --keys         interactive spike; send keys and check what each one did
                 (← opens the agent view, Esc closes it, shift+tab moves
                 the approval posture, ctrl+c quits).

Usage (from tui/):  python3 -I spike/check_tmux.py [--cols 110] [--rows 35] ...
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrollback_check import check, summary  # noqa: E402

TUI = Path(__file__).resolve().parents[1]
BUN = shutil.which("bun") or str(Path.home() / ".bun" / "bin" / "bun")


class Tmux:
    """A private tmux server with one session in it."""

    def __init__(self, cols: int, rows: int):
        self.socket = f"spikecheck{os.getpid()}"
        self.name = "spike"
        # One command: history-limit applies to panes created after it is
        # set, so it must be in place before the session's pane exists.
        self.run("-f", "/dev/null", "start-server", ";",
                 "set-option", "-g", "history-limit", "200000", ";",
                 "set-option", "-g", "allow-rename", "off", ";",
                 "new-session", "-d", "-s", self.name, "-x", str(cols),
                 "-y", str(rows), "bash --noprofile --norc")

    def run(self, *args: str) -> str:
        return subprocess.run(["tmux", "-L", self.socket, *args], check=True,
                              capture_output=True, text=True).stdout

    def type(self, text: str) -> None:
        self.run("send-keys", "-t", self.name, "-l", text)

    def key(self, *keys: str) -> None:
        self.run("send-keys", "-t", self.name, *keys)

    def screen(self) -> str:
        return self.run("capture-pane", "-p", "-t", self.name)

    def history(self) -> str:
        return self.run("capture-pane", "-p", "-S", "-", "-E", "-", "-t", self.name)

    def wait_for(self, pattern: str, timeout: float) -> re.Match | None:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            found = re.search(pattern, self.screen())
            if found:
                return found
            time.sleep(0.2)
        return None

    def close(self) -> None:
        subprocess.run(["tmux", "-L", self.socket, "kill-server"],
                       capture_output=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=110)
    ap.add_argument("--rows", type=int, default=35)
    ap.add_argument("--prefill", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--resize-to", type=int, default=0)
    ap.add_argument("--resize-after", type=float, default=6.0)
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--keep", action="store_true",
                    help="leave the tmux server running for inspection")
    ap.add_argument("--ssh-key", default="")
    ap.add_argument("--keys", action="store_true")
    args = ap.parse_args()

    work = Path(tempfile.mkdtemp(prefix="spike-tmux-"))
    report_path = work / "report.json"
    stderr_path = work / "stderr.txt"
    mode = "" if args.keys else f" --auto --seed {args.seed}"
    inner = (f"cd {shlex.quote(str(TUI))} && SPIKE_REPORT={report_path} "
             f"{BUN} spike/inline.tsx{mode} 2>{stderr_path}")
    if args.ssh_key:
        inner = (f"ssh -tt -i {args.ssh_key} -o BatchMode=yes -o LogLevel=ERROR "
                 f"localhost {shlex.quote(inner)}")
    command = f"{inner}; echo SPIKE_EXIT=$?"
    if args.prefill:
        command = f"seq -f 'shell line %g' {args.prefill}; {command}"

    tmux = Tmux(args.cols, args.rows)
    problems: list[str] = []
    resized = False
    try:
        tmux.type(command)
        tmux.key("Enter")
        if args.keys:
            problems += drive_keys(tmux)
        started = time.monotonic()
        exit_code = None
        while time.monotonic() - started < args.timeout:
            time.sleep(0.3)
            if args.resize_to and not resized and \
                    time.monotonic() - started > args.resize_after:
                tmux.run("resize-window", "-t", tmux.name, "-x", str(args.resize_to))
                resized = True
            found = re.search(r"SPIKE_EXIT=(\d+)", tmux.screen())
            if found:
                exit_code = int(found.group(1))
                break
        history = tmux.history()
    finally:
        if not args.keep:
            tmux.close()
    (work / "history.txt").write_text(history)

    if exit_code is None:
        problems.append(f"spike did not finish within {args.timeout:.0f}s")
    elif exit_code != 0:
        problems.append(f"spike exited with {exit_code}")
    report = json.loads(report_path.read_text()) if report_path.exists() else {}
    if not report:
        problems.append("no report written")
    stderr = stderr_path.read_text() if stderr_path.exists() else ""
    if stderr.strip():
        problems.append(f"stderr: {stderr.strip()[:300]}")
    if report and not args.keys:
        problems += check(history, report,
                          prefill=0 if resized else args.prefill,
                          replays=1 if resized else 0,
                          width=args.resize_to if resized else args.cols)
    if report and args.keys:
        want = ["agents:open", "agents:close", "posture:default", "posture:acceptEdits"]
        if report.get("keyLog") != want:
            problems.append(f"keys did: {report.get('keyLog')} (expected {want})")

    where = f"tmux {tmux_version()} {args.cols}x{args.rows}"
    if args.prefill:
        where += f", from row {args.prefill + 1}"
    if args.ssh_key:
        where += ", over ssh"
    if resized:
        where += f", resized to {args.resize_to} cols"
    print(f"{where}: {summary(history, report) if report else 'no report'}; "
          f"work dir {work}")
    if report.get("replay1000Ms") is not None:
        print(f"  1000-line replay {report['replay1000Ms']} ms, "
              f"30 realistic turns {report.get('realReplay30Ms')} ms")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS")
    return 0


def drive_keys(tmux: Tmux) -> list[str]:
    """Keys a user presses, each after the screen shows its effect."""
    if not tmux.wait_for(r"ask DELFIN", 30):
        return ["the interactive spike did not come up"]
    steps = [(["Left"], r"DELFIN agents"),          # ← opens the agent view
             (["Escape"], r"ask DELFIN"),           # Esc closes it
             (["BTab"], r"default · esc"),          # shift+tab: plan → default
             (["BTab"], r"acceptEdits · esc")]      # → acceptEdits
    for keys, effect in steps:
        tmux.key(*keys)
        if not tmux.wait_for(effect, 10):
            return [f"after {keys} the screen never showed {effect!r}"]
        time.sleep(0.3)
    tmux.key("C-c")
    return []


def tmux_version() -> str:
    out = subprocess.run(["tmux", "-V"], capture_output=True, text=True).stdout
    return out.strip().removeprefix("tmux ")


if __name__ == "__main__":
    sys.exit(main())
