#!/usr/bin/env bash
# Phase-0 spike: every automated terminal check in one go.
#   1. tmux, no resize          scrollback integrity after footer growth
#   2. tmux, resize to 80 cols  clear-and-replay at the new width
#   3. xterm.js (headless)      the emulator family of OpenTUI issue #1630,
#                               fed the raw bytes of a real pseudo-terminal run
set -u
cd "$(dirname "$0")/.."
export PATH="$HOME/.bun/bin:$PATH"
fail=0

echo "== 1. tmux 110x35, no resize"
python3 -I spike/check_tmux.py --cols 110 --rows 35 || fail=1

echo "== 2. tmux 110x35, resized to 80 columns mid-run"
python3 -I spike/check_tmux.py --cols 110 --rows 35 --resize-to 80 --resize-after 4 || fail=1

echo "== 3. xterm.js 110x35 (recorded in a pty)"
work=$(mktemp -d)
script -q -e -f -c "stty cols 110 rows 35; SPIKE_REPORT=$work/report.json bun spike/inline.tsx --auto 2>$work/stderr.txt" \
  "$work/rec.raw" >/dev/null 2>&1 || fail=1
if [ -s "$work/stderr.txt" ]; then echo "stderr:"; head -20 "$work/stderr.txt"; fail=1; fi
bun spike/check_xterm.ts "$work/rec.raw" "$work/report.json" 110 35 || fail=1
clears=$(python3 -I -c 'import re,sys; print(len(re.findall(rb"\x1b\[3J", open(sys.argv[1],"rb").read())))' "$work/rec.raw")
echo "clear-scrollback sequences without a resize: $clears"
[ "$clears" = "0" ] || fail=1

if [ "$fail" = 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; fi
exit "$fail"
