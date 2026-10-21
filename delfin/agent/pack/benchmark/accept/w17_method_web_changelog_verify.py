#!/usr/bin/env python3
"""Acceptance entry: web changelog task (delegates to the shared verifier)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from w17_method_web_verify import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main(sys.argv, variant="changelog"))
