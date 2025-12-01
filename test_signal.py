#!/usr/bin/env python3
"""Test if SIGINT signals are received."""
import signal
import sys
import time

def handle_sigint(signum, frame):
    print(f"\n🛑 SIGINT received! (signal {signum})", file=sys.stderr, flush=True)
    sys.exit(130)

def handle_sigterm(signum, frame):
    print(f"\n🛑 SIGTERM received! (signal {signum})", file=sys.stderr, flush=True)
    sys.exit(143)

# Install handlers
signal.signal(signal.SIGINT, handle_sigint)
signal.signal(signal.SIGTERM, handle_sigterm)

print("✓ Signal handlers installed")
print("✓ Press Ctrl+C to test...")
print("✓ Waiting (will timeout after 30 seconds)...")

# Wait
try:
    time.sleep(30)
    print("Timeout - no signal received")
except KeyboardInterrupt:
    print("KeyboardInterrupt caught (this shouldn't happen if handler works)")
