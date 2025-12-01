#!/usr/bin/env python3
"""Minimal test for Ctrl+C in tmux"""
import signal
import sys
import time

caught = False

def handler(sig, frame):
    global caught
    caught = True
    print("\n✓✓✓ GOT IT! Signal caught!", flush=True)
    sys.exit(0)

signal.signal(signal.SIGINT, handler)
print("Ready. Press Ctrl+C (you have 10 seconds)...")
sys.stdout.flush()

time.sleep(10)

if not caught:
    print("✗ No signal received after 10 seconds")
    sys.exit(1)
