#!/usr/bin/env bash
# Run the spike's --auto script in a real GNOME Terminal (VTE) window on the
# desktop: screenshots while it runs, then select-all + copy of the whole
# scrollback, checked with the same rules as the tmux and xterm.js checks.
#
# Needs the desktop session's DISPLAY, XAUTHORITY and DBUS_SESSION_BUS_ADDRESS
# in the environment, plus xdotool, xclip and ImageMagick's `import`.
# Usage (from tui/): spike/check_gnome.sh <output-dir>
set -u
cd "$(dirname "$0")/.."
out=${1:?output directory}
mkdir -p "$out"
bun=$HOME/.bun/bin/bun

before=$(xdotool search --class gnome-terminal 2>/dev/null | sort)
gnome-terminal --geometry=110x35 -- bash -c \
  "cd '$PWD' && SPIKE_REPORT='$out/report.json' '$bun' spike/inline.tsx --auto 2>'$out/stderr.txt'; echo SPIKE_EXIT=\$?; sleep 600"

wid=""
for _ in $(seq 1 50); do
  wid=$(comm -13 <(echo "$before") <(xdotool search --class gnome-terminal 2>/dev/null | sort) | tail -1)
  [ -n "$wid" ] && break
  sleep 0.2
done
[ -n "$wid" ] || { echo "no new gnome-terminal window"; exit 1; }
echo "window $wid"

start=$(date +%s.%N)
for t in 3 5.5 8 11 14 18; do
  while [ "$(echo "$(date +%s.%N) - $start < $t" | bc)" = 1 ]; do sleep 0.1; done
  import -window "$wid" "$out/shot-$t.png" 2>/dev/null
done

for _ in $(seq 1 120); do [ -s "$out/report.json" ] && break; sleep 0.5; done
sleep 1.5
import -window "$wid" "$out/shot-final.png" 2>/dev/null

# Select-all sets the X11 PRIMARY selection. Read that, never CLIPBOARD:
# the clipboard is shared with the host (VirtualBox), so copying would
# overwrite the user's clipboard and reading it could pick up their data.
# A marker in PRIMARY first tells a real selection from a stale one.
printf 'no-selection-yet' | xclip -selection primary -i
xdotool windowactivate --sync "$wid"
sleep 0.4
xdotool key ctrl+shift+a
sleep 0.8
xclip -o -selection primary > "$out/scrollback.txt" 2>/dev/null
if [ "$(head -c 16 "$out/scrollback.txt")" = "no-selection-yet" ]; then
  echo "select-all did not reach the window"; : > "$out/scrollback.txt"
fi
xdotool windowclose "$wid" 2>/dev/null || xdotool windowkill "$wid" 2>/dev/null

[ -s "$out/stderr.txt" ] && { echo "stderr:"; head -20 "$out/stderr.txt"; }
python3 -I spike/check_text.py "$out/scrollback.txt" "$out/report.json"
