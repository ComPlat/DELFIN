"""Acceptance self-test: always passes. Proves the runner reaches it."""
import sys
from pathlib import Path

ws = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
print(f"workspace seen: {ws}")
sys.exit(0)
