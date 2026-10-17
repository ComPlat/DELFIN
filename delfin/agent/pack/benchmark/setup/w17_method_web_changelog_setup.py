#!/usr/bin/env python3
"""Task setup entry: web changelog variant."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from w17_method_web_repo import build  # noqa: E402

if __name__ == "__main__":
    build(Path(sys.argv[1]), variant="changelog")
