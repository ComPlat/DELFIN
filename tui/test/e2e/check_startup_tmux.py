"""Launch to footer, live in tmux: with the agent service running, the
window's footer is on the screen within 300 ms of pressing enter on
`delfin-agent tui` (the median of five launches) -- with the window as one
executable, which this check builds first (scripts/build-binary.ts). From
source, Bun compiles the TypeScript at every start: that time is reported,
not held to the budget.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_startup_tmux.py
"""

from __future__ import annotations

import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[2]
REPO = TUI.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TUI / "spike"))

from check_tmux import BUN, Tmux  # noqa: E402

from delfin.agent.service.client import ServiceClient  # noqa: E402

PYTHON = os.environ.get("DELFIN_PYTHON") or sys.executable
BUDGET_S = 0.300
RUNS = 5


def launches(tmux: Tmux, env: dict, ws: Path) -> list[float]:
    command = (f"cd {REPO} && env -u DELFIN_SCRATCH_STATE "
               f"{' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    times = []
    for _ in range(RUNS):
        tmux.key("C-l")
        time.sleep(0.3)
        tmux.type(command)
        t0 = time.monotonic()
        tmux.key("Enter")
        while time.monotonic() - t0 < 30:
            screen = tmux.screen()
            if "ask DELFIN" in screen or "starting the session" in screen:
                break
            time.sleep(0.005)
        times.append(time.monotonic() - t0)
        time.sleep(0.8)
        tmux.key("C-d")
        tmux.wait_for(r"TUI_EXIT=\d+", 30)
        time.sleep(0.3)
    return times


def main() -> int:
    problems: list[str] = []
    built = subprocess.run([BUN, "scripts/build-binary.ts"], cwd=str(TUI), capture_output=True, text=True)
    if built.returncode != 0:
        print(f"FAIL: the binary did not build: {built.stdout} {built.stderr}")
        return 1
    tmp = Path(tempfile.mkdtemp(prefix="tui-start-"))
    runtime = Path(tempfile.mkdtemp(prefix="dts", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", ws):
        d.mkdir(parents=True)
    # A user's environment: a home of its own, no scratch state.
    env = {"HOME": str(tmp / "home"), "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    base = {k: v for k, v in os.environ.items() if k != "DELFIN_SCRATCH_STATE"}
    subprocess.run([PYTHON, "-m", "delfin.agent.service.launcher"], env={**base, **env},
                   check=True, capture_output=True, timeout=120)
    tmux = Tmux(100, 30)
    try:
        binary = launches(tmux, env, ws)
        source = launches(tmux, {**env, "DELFIN_TUI_BINARY": "/nonexistent"}, ws)
        with ServiceClient(runtime / f"delfin-{os.getuid()}" / "agent.sock") as c:
            c.request("service.shutdown")
    finally:
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    b, s = statistics.median(binary), statistics.median(source)
    print("binary: " + " ".join(f"{x * 1000:.0f}" for x in binary) + f" ms, median {b * 1000:.0f} ms")
    print("source: " + " ".join(f"{x * 1000:.0f}" for x in source) + f" ms, median {s * 1000:.0f} ms (not held to the budget)")
    if b > BUDGET_S:
        problems.append(f"launch to footer with the binary: median {b * 1000:.0f} ms, over {BUDGET_S * 1000:.0f} ms")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"PASS: the footer is on the screen {b * 1000:.0f} ms after enter (median, binary)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
