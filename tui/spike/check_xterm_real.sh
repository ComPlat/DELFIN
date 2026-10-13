#!/usr/bin/env bash
# Run the spike's --auto script in a real XTerm on a private Xvfb display
# (a free display number, picked by xvfb-run -a), then have XTerm dump its
# screen and scrollback (print-everything, bound to F12) and check the dump
# with scrollback_check.py. The key press goes to the private display only,
# to the XTerm this script started (found by its pid); nothing touches the
# user's desktop or clipboard.
# Usage (from tui/): spike/check_xterm_real.sh <out-dir> [cols] [rows] [prefill]
set -u
cd "$(dirname "$0")/.."
out=$(realpath -m "${1:?output dir}"); cols=${2:-110}; rows=${3:-35}; prefill=${4:-0}
mkdir -p "$out"; rm -f "$out"/dump.txt "$out"/report.json "$out"/stderr.txt

if [ -z "${SPIKE_XTERM_INNER:-}" ]; then
  exec env SPIKE_XTERM_INNER=1 xvfb-run -a -s "-screen 0 1920x1200x24 -nolisten tcp" \
    "$0" "$out" "$cols" "$rows" "$prefill"
fi

pre=""
[ "$prefill" -gt 0 ] && pre="seq -f 'shell line %g' $prefill; "
xterm -geometry ${cols}x${rows} -sl 200000 \
  -xrm "XTerm*printerCommand: cat > $out/dump.txt" \
  -xrm 'XTerm*printAttributes: 0' \
  -xrm 'XTerm*VT100.translations: #override <Key>F12: print-everything()' \
  -e bash -c "${pre}cd '$PWD' && SPIKE_REPORT='$out/report.json' $HOME/.bun/bin/bun spike/inline.tsx --auto 2>'$out/stderr.txt'; echo SPIKE_EXIT=\$?; sleep 300" &
xpid=$!

for _ in $(seq 1 360); do [ -s "$out/report.json" ] && break; sleep 0.5; done
# No window manager on the private display: give the XTerm this script
# started input focus explicitly, otherwise the key goes nowhere.
wid=$(timeout 10 xdotool search --sync --pid "$xpid" 2>/dev/null | head -1)
# This display holds nothing but that XTerm, so its class finds the same one.
[ -n "$wid" ] || wid=$(xdotool search --class xterm 2>/dev/null | head -1)
xdotool windowfocus --sync "$wid" 2>/dev/null
xdotool mousemove --window "$wid" 20 20
# The report is written before XTerm has drawn everything (a loaded machine
# draws half a megabyte of output slowly): dump until the last line is in.
for _ in $(seq 1 60); do
  rm -f "$out/dump.txt"
  xdotool key --window "$wid" F12
  for _ in $(seq 1 20); do [ -s "$out/dump.txt" ] && break; sleep 0.25; done
  grep -q "SPIKE_EXIT=" "$out/dump.txt" 2>/dev/null && break
  sleep 1
done
import -window root -crop 1000x800+0+0 "$out/xterm-final.png" 2>/dev/null
kill $xpid 2>/dev/null

fail=0
[ -s "$out/stderr.txt" ] && { echo "stderr:"; head -20 "$out/stderr.txt"; fail=1; }
grep -q "SPIKE_EXIT=0" "$out/dump.txt" 2>/dev/null || { echo "no SPIKE_EXIT=0 in the dump"; fail=1; }
echo "XTerm $(xterm -version) ${cols}x${rows}, from row $((prefill + 1))"
python3 -I spike/scrollback_check.py "$out/dump.txt" "$out/report.json" --prefill "$prefill" --replays 0 || fail=1
exit $fail
