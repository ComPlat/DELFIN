"""Run the spike in a real VTE terminal widget (the library GNOME Terminal is
built on) and check the scrollback through VTE's own API.

Runs on a private Xvfb display this script starts on a free display number
(``Xvfb -displayfd``), so nothing appears on the user's desktop, no
keystrokes are sent anywhere and the clipboard is never touched. The
scrollback is checked with scrollback_check.py: every row the spike
committed, exactly once and in order, and with --prefill the shell lines
printed before it still above it.

Usage (from tui/, with the system python that has python3-gi + VTE):
    /usr/bin/python3 -I spike/check_vte.py <out-dir> [--cols 110] [--rows 35]
                                                     [--prefill N] [--seed N]
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrollback_check import check, summary  # noqa: E402

TUI = Path(__file__).resolve().parents[1]
BUN = shutil.which("bun") or str(Path.home() / ".bun" / "bin" / "bun")


def private_display() -> tuple[subprocess.Popen, str]:
    read_fd, write_fd = os.pipe()
    xvfb = subprocess.Popen(["Xvfb", "-displayfd", str(write_fd), "-screen", "0",
                             "1920x1200x24", "-nolisten", "tcp"],
                            pass_fds=(write_fd,), stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    os.close(write_fd)
    with os.fdopen(read_fd) as fh:
        number = fh.readline().strip()
    return xvfb, f":{number}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--cols", type=int, default=110)
    ap.add_argument("--rows", type=int, default=35)
    ap.add_argument("--prefill", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    for stale in ("report.json", "stderr.txt", "scrollback.txt"):
        (out / stale).unlink(missing_ok=True)

    xvfb, display = private_display()
    os.environ["DISPLAY"] = display
    os.environ["NO_AT_BRIDGE"] = "1"
    os.environ.pop("WAYLAND_DISPLAY", None)
    try:
        return run(out, args, display)
    finally:
        xvfb.terminate()


def run(out: Path, args, display: str) -> int:
    import gi

    gi.require_version("Gtk", "3.0")
    gi.require_version("Vte", "2.91")
    from gi.repository import GLib, Gtk, Vte

    term = Vte.Terminal()
    term.set_scrollback_lines(-1)
    term.set_size(args.cols, args.rows)
    window = Gtk.Window()
    window.add(term)
    window.show_all()

    state = {"status": None}

    def exited(_term, status):
        state["status"] = status
        GLib.timeout_add(1500, Gtk.main_quit)

    term.connect("child-exited", exited)
    command = (f"cd {shlex.quote(str(TUI))} && SPIKE_REPORT={shlex.quote(str(out))}/report.json "
               f"{shlex.quote(BUN)} spike/inline.tsx --auto --seed {args.seed} "
               f"2>{shlex.quote(str(out))}/stderr.txt; echo SPIKE_EXIT=$?")
    if args.prefill:
        command = f"seq -f 'shell line %g' {args.prefill}; {command}"
    env = [f"{k}={v}" for k, v in os.environ.items()]
    term.spawn_async(Vte.PtyFlags.DEFAULT, str(TUI), ["/bin/bash", "-c", command], env,
                     GLib.SpawnFlags.DEFAULT, None, None, -1, None, None)
    GLib.timeout_add_seconds(180, Gtk.main_quit)  # backstop
    Gtk.main()

    adj = term.get_vadjustment()
    first, last = int(adj.get_lower()), int(adj.get_upper())
    text, _ = term.get_text_range_format(Vte.Format.TEXT, first, 0, last, args.cols)
    history = text or ""
    (out / "scrollback.txt").write_text(history)
    subprocess.run(["import", "-display", display, "-window", "root", "-crop",
                    f"{term.get_allocated_width()}x{term.get_allocated_height()}+0+0",
                    str(out / "vte-final.png")], check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    version = f"{Vte.get_major_version()}.{Vte.get_minor_version()}.{Vte.get_micro_version()}"
    problems: list[str] = []
    report = json.loads((out / "report.json").read_text()) if (out / "report.json").exists() else {}
    if state["status"] is None:
        problems.append("the spike did not finish")
    if "SPIKE_EXIT=0" not in history:
        problems.append("no SPIKE_EXIT=0 in the scrollback")
    stderr = (out / "stderr.txt").read_text() if (out / "stderr.txt").exists() else ""
    if stderr.strip():
        problems.append(f"stderr: {stderr.strip()[:300]}")
    if not report:
        problems.append("no report written")
    else:
        # No width check: VTE joins soft-wrapped lines in its text.
        problems += check(history, report, prefill=args.prefill, replays=0)
    where = f"VTE {version} {args.cols}x{args.rows}, from row {args.prefill + 1}, seed {args.seed}"
    print(f"{where}: {summary(history, report) if report else 'no report'}; out {out}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
