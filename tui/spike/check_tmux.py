"""Run the spike's --auto script inside tmux and check the scrollback it leaves.

tmux is a terminal emulator with its own history, so what it keeps after the
run is what a user could scroll back to. The checks are the phase-0 risks:

  * every marker the spike reported appears, in the order it was written;
  * all replay lines are present exactly once and in sequence (a footer that
    grows must not delete or duplicate rows: OpenTUI issue #1630);
  * the spike reported no errors.

With --resize the pane is resized mid-run, which makes the spike clear the
screen and scrollback and reprint the session; the same checks then apply
to the reprinted history.

Usage (from tui/):  python3 -I spike/check_tmux.py [--cols 110] [--rows 35]
                                                    [--resize-to 80 --resize-after 6]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[1]
BUN = shutil.which("bun") or str(Path.home() / ".bun" / "bin" / "bun")


def tmux(*args: str) -> str:
    return subprocess.run(["tmux", *args], check=True, capture_output=True,
                          text=True).stdout


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=110)
    ap.add_argument("--rows", type=int, default=35)
    ap.add_argument("--resize-to", type=int, default=0)
    ap.add_argument("--resize-after", type=float, default=6.0)
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--keep", action="store_true",
                    help="leave the tmux session open for inspection")
    ap.add_argument("--ssh-key", default="",
                    help="run the spike through `ssh -tt localhost` with this key")
    args = ap.parse_args()

    name = f"spikecheck{os.getpid()}"
    work = Path(tempfile.mkdtemp(prefix="spikecheck-"))
    report_path = work / "report.json"
    tmux("new-session", "-d", "-s", name, "-x", str(args.cols), "-y",
         str(args.rows), "bash --noprofile --norc")
    tmux("set-option", "-t", name, "history-limit", "50000")
    inner = (f"cd {TUI} && SPIKE_REPORT={report_path} {BUN} spike/inline.tsx"
             f" --auto 2>{work / 'stderr.txt'}")
    if args.ssh_key:
        # A real SSH hop: its own pty, window-size forwarding, and only the
        # environment SSH passes (COLORTERM, for one, is not forwarded).
        inner = (f"ssh -tt -i {args.ssh_key} -o BatchMode=yes localhost "
                 f"\"{inner}\"")
    command = f"{inner}; echo SPIKE_EXIT=$?"
    tmux("send-keys", "-t", name, command, "Enter")

    started = time.monotonic()
    resized = False
    exit_code = None
    while time.monotonic() - started < args.timeout:
        time.sleep(0.5)
        if args.resize_to and not resized and \
                time.monotonic() - started > args.resize_after:
            tmux("resize-window", "-t", name, "-x", str(args.resize_to))
            resized = True
        screen = tmux("capture-pane", "-p", "-t", name)
        found = re.search(r"SPIKE_EXIT=(\d+)", screen)
        if found:
            exit_code = int(found.group(1))
            break
    history = tmux("capture-pane", "-p", "-J", "-S", "-", "-E", "-", "-t", name)
    (work / "history.txt").write_text(history)
    if not args.keep:
        tmux("kill-session", "-t", name)

    problems: list[str] = []
    if exit_code is None:
        problems.append(f"spike did not finish within {args.timeout:.0f}s")
    elif exit_code != 0:
        problems.append(f"spike exited with {exit_code}")
    report = json.loads(report_path.read_text()) if report_path.exists() else {}
    if not report:
        problems.append("no report written")
    for err in report.get("errors", []):
        problems.append(f"spike error: {err}")
    stderr = (work / "stderr.txt").read_text() if (work / "stderr.txt").exists() else ""
    if stderr.strip():
        problems.append(f"stderr: {stderr.strip()[:300]}")

    # Markers in order (the session is linear, so first occurrences after the
    # previous marker must exist).
    pos = 0
    for marker in report.get("markers", []):
        at = history.find(marker, pos)
        if at < 0:
            problems.append(f"marker missing or out of order: {marker!r}")
        else:
            pos = at

    numbers = [int(n) for n in re.findall(r"replay line (\d+) —", history)]
    if numbers != list(range(1, len(numbers) + 1)) or len(numbers) != 1000:
        dupes = len(numbers) - len(set(numbers))
        missing = sorted(set(range(1, 1001)) - set(numbers))
        problems.append(f"replay lines: {len(numbers)} seen, {dupes} duplicated,"
                        f" {len(missing)} missing (first missing: {missing[:5]})")

    print(f"history lines: {history.count(chr(10))}  work dir: {work}")
    print(f"replay of 1000 lines: {report.get('replay1000Ms', float('nan')):.0f} ms")
    if resized:
        print(f"resized to {args.resize_to} columns mid-run")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: markers in order, 1000 replay lines exactly once, no errors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
