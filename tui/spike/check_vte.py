"""Run the spike's --auto script inside a real VTE terminal widget (the
library GNOME Terminal is built on) and read the scrollback back through
VTE's own API.

Runs on a private Xvfb display started here, so nothing appears on the
user's desktop, no keystrokes are sent anywhere and the clipboard is never
touched. The captured text is checked with check_text.py's rules.

Usage (from tui/, with the system python that has python3-gi + VTE):
    /usr/bin/python3 -I spike/check_vte.py <out-dir> [cols] [rows]
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[1]
BUN = shutil.which("bun") or str(Path.home() / ".bun" / "bin" / "bun")


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    cols = int(sys.argv[2]) if len(sys.argv) > 2 else 110
    rows = int(sys.argv[3]) if len(sys.argv) > 3 else 35
    out.mkdir(parents=True, exist_ok=True)

    display = ":97"
    xvfb = subprocess.Popen(["Xvfb", display, "-screen", "0", "1920x1200x24", "-nolisten", "tcp"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(1.0)
    os.environ["DISPLAY"] = display
    os.environ["NO_AT_BRIDGE"] = "1"
    os.environ.pop("WAYLAND_DISPLAY", None)
    try:
        return run(out, cols, rows, display)
    finally:
        xvfb.terminate()


def run(out: Path, cols: int, rows: int, display: str) -> int:
    import gi

    gi.require_version("Gtk", "3.0")
    gi.require_version("Vte", "2.91")
    from gi.repository import GLib, Gtk, Vte

    term = Vte.Terminal()
    term.set_scrollback_lines(-1)
    term.set_size(cols, rows)
    window = Gtk.Window()
    window.add(term)
    window.show_all()

    state = {"status": None}

    def exited(_term, status):
        state["status"] = status
        GLib.timeout_add(1500, Gtk.main_quit)

    term.connect("child-exited", exited)
    command = (f"cd '{TUI}' && SPIKE_REPORT='{out}/report.json' '{BUN}' spike/inline.tsx --auto"
               f" 2>'{out}/stderr.txt'; echo SPIKE_EXIT=$?")
    env = [f"{k}={v}" for k, v in os.environ.items()]
    term.spawn_async(Vte.PtyFlags.DEFAULT, str(TUI), ["/bin/bash", "-c", command], env,
                     GLib.SpawnFlags.DEFAULT, None, None, -1, None, None)
    GLib.timeout_add_seconds(150, Gtk.main_quit)  # backstop
    Gtk.main()

    adj = term.get_vadjustment()
    first, last = int(adj.get_lower()), int(adj.get_upper())
    text, _ = term.get_text_range_format(Vte.Format.TEXT, first, 0, last, cols)
    (out / "scrollback.txt").write_text(text or "")
    subprocess.run(["import", "-display", display, "-window", "root", "-crop",
                    f"{term.get_allocated_width()}x{term.get_allocated_height()}+0+0",
                    str(out / "vte-final.png")], check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"VTE {Vte.get_major_version()}.{Vte.get_minor_version()}.{Vte.get_micro_version()}"
          f" {cols}x{rows}, child status {state['status']}, rows {first}..{last}")
    stderr = (out / "stderr.txt").read_text() if (out / "stderr.txt").exists() else ""
    if stderr.strip():
        print("stderr:", stderr.strip()[:300])
    if state["status"] is None or not (out / "report.json").exists():
        print("FAIL: the spike did not finish")
        return 1
    return subprocess.run([sys.executable, "-I", str(TUI / "spike" / "check_text.py"),
                           str(out / "scrollback.txt"), str(out / "report.json")]).returncode


if __name__ == "__main__":
    sys.exit(main())
