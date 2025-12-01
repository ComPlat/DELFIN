#!/bin/bash
# Reliably stop DELFIN and all ORCA processes in current directory

set -e

CURRENT_DIR=$(pwd)
echo "Looking for DELFIN process in: $CURRENT_DIR"

# Find DELFIN process running in this directory
DELFIN_PID=""
for pid in $(pgrep -f "python.*delfin" 2>/dev/null || true); do
    if [ -d "/proc/$pid" ]; then
        cwd=$(readlink -f "/proc/$pid/cwd" 2>/dev/null || echo "")
        if [ "$cwd" == "$CURRENT_DIR" ]; then
            DELFIN_PID=$pid
            break
        fi
    fi
done

if [ -z "$DELFIN_PID" ]; then
    echo "No DELFIN process found in $CURRENT_DIR"
    exit 1
fi

echo "Found DELFIN process: PID $DELFIN_PID"
echo "Killing all ORCA child processes..."

# Kill all ORCA processes that are children of this DELFIN
pkill -9 -P "$DELFIN_PID" 2>/dev/null || true

# Also kill ORCA processes in this directory
for orca_pid in $(pgrep -f "orca.*\.inp" 2>/dev/null || true); do
    if [ -d "/proc/$orca_pid" ]; then
        orca_cwd=$(readlink -f "/proc/$orca_pid/cwd" 2>/dev/null || echo "")
        if [[ "$orca_cwd" == "$CURRENT_DIR"* ]]; then
            echo "Killing ORCA PID $orca_pid"
            kill -9 "$orca_pid" 2>/dev/null || true
        fi
    fi
done

echo "Killing DELFIN process $DELFIN_PID..."
kill -9 "$DELFIN_PID" 2>/dev/null || true

sleep 0.5

if ps -p "$DELFIN_PID" >/dev/null 2>&1; then
    echo "✗ Process still running!"
    exit 1
else
    echo "✓ DELFIN stopped successfully"
fi
