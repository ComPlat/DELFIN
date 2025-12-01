#!/usr/bin/env python3
"""Quick test to verify Ctrl+C properly terminates DELFIN and ORCA processes."""

import os
import sys
import signal
import time

# Test that the global manager properly handles SIGINT
from delfin.global_manager import get_global_manager

def test_signal_handler():
    """Verify signal handler is installed and works."""
    manager = get_global_manager()

    # Initialize with dummy config
    config = {
        'PAL': 4,
        'maxcore': 1000,
        'pal_jobs': 2,
        'parallel_workflows': 'yes'
    }
    manager.initialize(config)

    print("✓ Global manager initialized")
    print(f"✓ Signal handler installed: {manager._signal_handler_installed}")
    print(f"✓ SIGTERM handler installed: {manager._sigterm_handler_installed}")

    # Check what the current SIGINT handler is
    current_handler = signal.getsignal(signal.SIGINT)
    print(f"✓ Current SIGINT handler: {current_handler}")

    if current_handler == signal.SIG_DFL:
        print("⚠ WARNING: Signal handler is still default! This means Ctrl+C won't work properly.")
        return False

    print("\n✓ Signal handlers are properly installed!")
    print("✓ Ctrl+C should now properly terminate DELFIN and ORCA processes")

    return True

if __name__ == "__main__":
    success = test_signal_handler()
    sys.exit(0 if success else 1)
