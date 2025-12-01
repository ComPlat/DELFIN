#!/bin/bash
# Test if signals work in tmux

echo "Testing signal handling in tmux..."
echo "Current session: $TMUX"

if [ -n "$TMUX" ]; then
    echo "✓ Running in tmux"
else
    echo "✗ Not in tmux"
fi

# Run a simple Python script that should catch SIGINT
python3.13 - << 'EOF'
import signal
import sys
import time

def handler(sig, frame):
    print("\n✓✓✓ SIGNAL RECEIVED!", file=sys.stderr, flush=True)
    print("\n✓✓✓ SIGNAL RECEIVED!", file=sys.stdout, flush=True)
    with open("/tmp/test_signal_received.txt", "w") as f:
        f.write("Signal received!\n")
    sys.exit(0)

signal.signal(signal.SIGINT, handler)
signal.signal(signal.SIGTERM, handler)

print("Signal handler installed. Press Ctrl+C...")
print("Waiting 30 seconds...")

try:
    time.sleep(30)
except KeyboardInterrupt:
    print("KeyboardInterrupt (this means handler didn't work)")
    sys.exit(1)

print("Timeout - no signal received")
EOF

# Check if file was created
if [ -f "/tmp/test_signal_received.txt" ]; then
    echo "✓ Signal was received!"
    cat /tmp/test_signal_received.txt
    rm /tmp/test_signal_received.txt
else
    echo "✗ Signal was NOT received"
fi
