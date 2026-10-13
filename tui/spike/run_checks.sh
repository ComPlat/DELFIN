#!/usr/bin/env bash
# The inline-mode regression suite: every terminal the TUI must work in,
# from the start rows and link speeds that broke it before. Each check
# compares the terminal's real scrollback with the rows the spike committed
# (scrollback_check.py): every row exactly once, in order, nothing above
# the spike overwritten.
#
#   spike/run_checks.sh          everything available on this machine
#   spike/run_checks.sh quick    unit tests, typecheck, xterm.js and tmux only
#
# Terminals that are not installed are skipped and say so. Window-based
# checks (VTE, XTerm) run on private Xvfb displays; nothing touches the
# user's desktop. Results and logs go to $SPIKE_CHECKS_OUT (default: a new
# temporary directory).
set -u
cd "$(dirname "$0")/.."
export PATH="$HOME/.bun/bin:$PATH"
out=${SPIKE_CHECKS_OUT:-$(mktemp -d -t spike-checks-XXXX)}
mkdir -p "$out"
quick=${1:-}
passed=0; failed=0; skipped=0; known=0
failures=()

run() {   # run <label> <command...>
  local label=$1; shift
  if "$@" >"$out/$(echo "$label" | tr ' /' '__').log" 2>&1; then
    passed=$((passed + 1)); echo "  PASS  $label"
  else
    failed=$((failed + 1)); failures+=("$label"); echo "  FAIL  $label"
    tail -5 "$out/$(echo "$label" | tr ' /' '__').log" | sed 's/^/          /'
  fi
}
known_limit() {   # known_limit <label> <why> <command...>: expected to fail
  local label=$1 why=$2; shift 2
  if "$@" >"$out/$(echo "$label" | tr ' /' '__').log" 2>&1; then
    echo "  PASS  $label (no longer a known limitation: $why)"; passed=$((passed + 1))
  else
    known=$((known + 1)); echo "  KNOWN $label: $why"
  fi
}
skip() { skipped=$((skipped + 1)); echo "  SKIP  $1"; }

# The Python that has delfin (the protocol types and the live window need it).
export DELFIN_PYTHON=${DELFIN_PYTHON:-$( [ -x "$HOME/software/delfin/.venv/bin/python" ] && echo "$HOME/software/delfin/.venv/bin/python" || command -v python3)}

echo "== unit tests and types"
run "bun test" bun test
run "typecheck" bun x tsc --noEmit
run "protocol types match the Python definition" bun scripts/gen-protocol.ts --check
run "native build is the patched one" scripts/build-opentui-native.sh --check

echo "== the window, live"
if command -v tmux >/dev/null; then
  run "delfin-agent tui in tmux shows its snapshot" "$DELFIN_PYTHON" test/e2e/check_tui_tmux.py
else
  skip "the live window (no tmux)"
fi

echo "== xterm.js (headless, live; replies delayed to play a slow link)"
for latency in 0 300 2000; do
  for prefill in 0 34; do
    run "xterm.js latency ${latency}ms from row $((prefill + 1))" \
      bun spike/check_xterm.ts --cols 110 --rows 35 --latency "$latency" --prefill "$prefill"
  done
done

echo "== tmux $(tmux -V 2>/dev/null | cut -d' ' -f2)"
if command -v tmux >/dev/null; then
  run "tmux from row 1" python3 -I spike/check_tmux.py
  run "tmux from the bottom row" python3 -I spike/check_tmux.py --prefill 40
  run "tmux resized 110 to 80" python3 -I spike/check_tmux.py --resize-to 80 --resize-after 6
  run "tmux keys" python3 -I spike/check_tmux.py --keys
  key=$HOME/.ssh/delfin_tui_localhost_test
  if [ -r "$key" ] && ssh -i "$key" -o BatchMode=yes -o ConnectTimeout=5 localhost true 2>/dev/null; then
    run "tmux over ssh from the bottom row" python3 -I spike/check_tmux.py --ssh-key "$key" --prefill 40
  else
    skip "tmux over ssh (no key for ssh to localhost)"
  fi
else
  skip "tmux (not installed)"
fi
[ "$quick" = quick ] && { echo "(quick: window terminals and GNU screen left out)"; }

if [ "$quick" != quick ]; then
  echo "== GNU screen"
  if command -v screen >/dev/null; then
    run "screen from row 1" python3 -I spike/check_screen.py
    run "screen from the bottom row" python3 -I spike/check_screen.py --prefill 34
    known_limit "screen resized 110 to 80" \
      "no ED3: a replay cannot erase screen's scrollback (phase 4 prints a divider instead)" \
      python3 -I spike/check_screen.py --resize-to 80 --resize-after 6
  else
    skip "GNU screen (not installed)"
  fi

  echo "== VTE (GNOME Terminal's widget), private display"
  if command -v Xvfb >/dev/null && /usr/bin/python3 -c "import gi; gi.require_version('Vte', '2.91')" 2>/dev/null; then
    run "VTE from row 1" /usr/bin/python3 -I spike/check_vte.py "$out/vte-1"
    run "VTE from the bottom row" /usr/bin/python3 -I spike/check_vte.py "$out/vte-35" --prefill 34
    run "VTE 80x24 from the bottom row" /usr/bin/python3 -I spike/check_vte.py "$out/vte-80" --cols 80 --rows 24 --prefill 23
  else
    skip "VTE (needs Xvfb and python3-gi with VTE 2.91)"
  fi

  echo "== XTerm, private display"
  if command -v xterm >/dev/null && command -v xvfb-run >/dev/null && command -v xdotool >/dev/null; then
    run "XTerm from row 1" spike/check_xterm_real.sh "$out/xterm-1" 110 35 0
    run "XTerm from the bottom row" spike/check_xterm_real.sh "$out/xterm-35" 110 35 34
  else
    skip "XTerm (needs xterm, xvfb-run and xdotool)"
  fi
fi

echo
echo "$passed passed, $failed failed, $known known limitations, $skipped skipped; logs in $out"
[ "$failed" = 0 ] || { printf '  failed: %s\n' "${failures[@]}"; exit 1; }
