#!/usr/bin/env bash
# Run the spike's --auto script in a real xterm on a private Xvfb display,
# then have xterm dump its screen and scrollback (print-everything, bound to
# F12) and check the dump with check_text.py. The key press goes to the
# private display only; nothing touches the user's desktop or clipboard.
# Usage (from tui/): spike/check_xterm_real.sh <out-dir> [cols] [rows]
set -u
cd "$(dirname "$0")/.."
out=$(realpath -m "${1:?output dir}"); cols=${2:-110}; rows=${3:-35}
mkdir -p "$out"; rm -f "$out"/dump.txt "$out"/report.json
display=:96
Xvfb $display -screen 0 1920x1200x24 -nolisten tcp >/dev/null 2>&1 &
xvfb=$!
trap 'kill $xvfb 2>/dev/null' EXIT
sleep 1
export DISPLAY=$display

xterm -geometry ${cols}x${rows} -sl 100000 \
  -xrm "XTerm*printerCommand: cat > $out/dump.txt" \
  -xrm 'XTerm*printAttributes: 0' \
  -xrm 'XTerm*VT100.translations: #override <Key>F12: print-everything()' \
  -e bash -c "cd '$PWD' && SPIKE_REPORT='$out/report.json' $HOME/.bun/bin/bun spike/inline.tsx --auto 2>'$out/stderr.txt'; echo SPIKE_EXIT=\$?; sleep 300" &
xpid=$!

for _ in $(seq 1 300); do [ -s "$out/report.json" ] && break; sleep 0.5; done
sleep 2
import -window root -crop 1000x800+0+0 "$out/xterm-final.png" 2>/dev/null
# No window manager on the private display: give xterm input focus
# explicitly, otherwise the key goes nowhere.
wid=$(xdotool search --class xterm | head -1)
xdotool windowfocus --sync "$wid" 2>/dev/null
xdotool mousemove --window "$wid" 20 20
xdotool key F12
for _ in $(seq 1 20); do [ -s "$out/dump.txt" ] && break; sleep 0.5; done
kill $xpid 2>/dev/null

[ -s "$out/stderr.txt" ] && { echo "stderr:"; head -20 "$out/stderr.txt"; }
echo "xterm $(xterm -version) ${cols}x${rows}"
python3 -I spike/check_text.py "$out/dump.txt" "$out/report.json"
