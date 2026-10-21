"""Run the spike inside GNU screen and check screen's own scrollback.

GNU screen is common on HPC login nodes and differs from the other
terminals in the ways that matter here: it has no ED3 (clear scrollback),
no truecolor, and its own idea of scroll regions. It runs in a private tmux
pane (check_tmux.py's server) so it has a real terminal, with its socket
directory under the work directory and no screenrc; its scrollback is read
back with ``hardcopy -h`` and checked with scrollback_check.py.

Usage (from tui/):  python3 -I spike/check_screen.py [--cols 110] [--rows 35]
                                                     [--prefill N] [--seed N]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_tmux import BUN, TUI, Tmux  # noqa: E402
from scrollback_check import check, summary  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=110)
    ap.add_argument("--rows", type=int, default=35)
    ap.add_argument("--prefill", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--resize-to", type=int, default=0)
    ap.add_argument("--resize-after", type=float, default=6.0)
    args = ap.parse_args()

    work = Path(tempfile.mkdtemp(prefix="spike-screen-"))
    screendir = work / "screens"
    screendir.mkdir(mode=0o700)
    report_path, stderr_path, dump = work / "report.json", work / "stderr.txt", work / "hardcopy.txt"
    spike = (f"cd {shlex.quote(str(TUI))} && SPIKE_REPORT={report_path} {BUN} "
             f"spike/inline.tsx --auto --seed {args.seed} 2>{stderr_path}; echo SPIKE_EXIT=$?")
    if args.prefill:
        spike = f"seq -f 'shell line %g' {args.prefill}; {spike}"
    name = f"spike{os.getpid()}"
    env = {**os.environ, "SCREENDIR": str(screendir)}
    inner = (f"SCREENDIR={shlex.quote(str(screendir))} screen -c /dev/null -h 200000 "
             f"-S {name} bash -c {shlex.quote(spike + '; sleep 600')}")

    tmux = Tmux(args.cols, args.rows)
    problems: list[str] = []
    exit_code = None
    try:
        tmux.type(inner)
        tmux.key("Enter")
        started = time.monotonic()
        resized = False
        while time.monotonic() - started < args.timeout:
            time.sleep(0.5)
            if args.resize_to and not resized and \
                    time.monotonic() - started > args.resize_after:
                tmux.run("resize-window", "-t", tmux.name, "-x", str(args.resize_to))
                resized = True
            found = re.search(r"SPIKE_EXIT=(\d+)", tmux.screen())
            if found:
                exit_code = int(found.group(1))
                break
        subprocess.run(["screen", "-S", name, "-X", "hardcopy", "-h", str(dump)],
                       env=env, capture_output=True, timeout=30)
        for _ in range(20):
            if dump.exists() and dump.stat().st_size:
                break
            time.sleep(0.25)
        subprocess.run(["screen", "-S", name, "-X", "quit"], env=env,
                       capture_output=True, timeout=30)
    finally:
        tmux.close()

    history = dump.read_text(errors="replace") if dump.exists() else ""
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
    if not history:
        problems.append("screen wrote no hardcopy")
    elif report:
        # GNU screen has no ED3: a replay cannot erase what it reprints.
        problems += check(history, report, prefill=0 if args.resize_to else args.prefill,
                          replays=1 if args.resize_to else 0)
    version = subprocess.run(["screen", "--version"], capture_output=True,
                             text=True).stdout.split()[2]
    where = f"GNU screen {version} {args.cols}x{args.rows}, from row {args.prefill + 1}"
    if args.resize_to:
        where += f", resized to {args.resize_to} cols"
    print(f"{where}: {summary(history, report) if report else 'no report'}; work dir {work}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
