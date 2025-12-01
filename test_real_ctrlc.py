#!/usr/bin/env python3
"""Test Ctrl+C behavior with real ORCA subprocess."""

import os
import sys
import signal
import time
import subprocess

from delfin.global_manager import get_global_manager

def test_ctrlc_with_subprocess():
    """Test that Ctrl+C properly kills a subprocess."""
    print("=" * 70)
    print("CTRL+C TEST WITH SUBPROCESS")
    print("=" * 70)

    # Initialize global manager
    manager = get_global_manager()
    config = {
        'PAL': 12,
        'maxcore': 6000,
        'pal_jobs': 4,
        'parallel_workflows': 'yes'
    }
    manager.initialize(config)

    print(f"\n✓ Global manager initialized")
    print(f"✓ SIGINT handler: {manager._signal_handler_installed}")
    print(f"✓ SIGTERM handler: {manager._sigterm_handler_installed}")

    # Get current SIGINT handler
    current_sigint = signal.getsignal(signal.SIGINT)
    print(f"✓ Current SIGINT handler: {current_sigint}")
    print(f"✓ Expected handler: {manager._handle_sigint}")

    if current_sigint != manager._handle_sigint:
        print("\n⚠ WARNING: SIGINT handler was not installed correctly!")
        print(f"  Expected: {manager._handle_sigint}")
        print(f"  Got: {current_sigint}")
        return False

    # Start a dummy long-running subprocess (simulating ORCA)
    print("\n▶ Starting dummy subprocess (sleep 300)...")
    try:
        # Use process group like ORCA subprocess
        proc = subprocess.Popen(
            ["sleep", "300"],
            start_new_session=True,  # Create new process group
        )

        # Register the subprocess with the manager
        token = manager.register_subprocess(
            proc,
            label="test_sleep",
            cwd=os.getcwd()
        )
        print(f"✓ Subprocess started (PID={proc.pid}, token={token})")
        print(f"✓ Tracked processes: {len(manager._tracked_processes)}")

        # Give it a moment to start
        time.sleep(0.5)

        # Check if subprocess is running
        if proc.poll() is not None:
            print("⚠ Subprocess exited early!")
            return False

        print(f"\n✓ Subprocess is running (PID={proc.pid})")
        print("\n" + "=" * 70)
        print("NOW PRESS Ctrl+C TO TEST SIGNAL HANDLING")
        print("The subprocess should be killed within ~1 second")
        print("=" * 70)

        # Wait for Ctrl+C or subprocess to finish
        try:
            proc.wait()
        except KeyboardInterrupt:
            print("\n🛑 Ctrl+C received! Cleanup should happen automatically...")
            raise

        print("\n⚠ Subprocess finished without Ctrl+C")
        return False

    except KeyboardInterrupt:
        # This is expected!
        print("\n✓ KeyboardInterrupt caught as expected")

        # Check if subprocess was killed
        time.sleep(1)  # Give cleanup time to complete
        returncode = proc.poll()

        if returncode is not None:
            print(f"✓ Subprocess was killed (exit code: {returncode})")
            print(f"✓ Tracked processes remaining: {len(manager._tracked_processes)}")
            return True
        else:
            print(f"⚠ Subprocess is still running! (PID={proc.pid})")
            print("This means signal handling did NOT work correctly!")
            # Clean up manually
            try:
                proc.kill()
                proc.wait(timeout=2)
            except:
                pass
            return False
    finally:
        # Cleanup
        manager.shutdown()

if __name__ == "__main__":
    try:
        success = test_ctrlc_with_subprocess()
        if success:
            print("\n" + "=" * 70)
            print("✓ TEST PASSED: Ctrl+C correctly kills subprocesses")
            print("=" * 70)
            sys.exit(0)
        else:
            print("\n" + "=" * 70)
            print("✗ TEST FAILED: Ctrl+C did not work correctly")
            print("=" * 70)
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ TEST ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
