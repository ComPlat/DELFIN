#!/usr/bin/env python3
"""Check what signal handler is installed for SIGINT"""
import signal
import sys

# Import delfin to trigger handler installation
from delfin.global_manager import get_global_manager

# Initialize the manager (this should install handlers)
mgr = get_global_manager()

# Check what handler is installed
handler = signal.getsignal(signal.SIGINT)
print(f"SIGINT handler: {handler}")
print(f"Handler type: {type(handler)}")

if hasattr(handler, '__name__'):
    print(f"Handler name: {handler.__name__}")
if hasattr(handler, '__module__'):
    print(f"Handler module: {handler.__module__}")

# Check if it's the default
if handler == signal.SIG_DFL:
    print("❌ Using DEFAULT handler (not our custom handler!)")
elif handler == signal.SIG_IGN:
    print("❌ Signal is IGNORED")
elif callable(handler):
    print("✓ Custom handler is installed")
else:
    print(f"⚠️  Unknown handler type: {handler}")
