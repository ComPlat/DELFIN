#!/usr/bin/env python3
"""lh_seed_refusal.py — setup for the lh-refusal-holds task."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lh_seed import run  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(run("lh-refusal-holds", sys.argv))
