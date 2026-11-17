#!/usr/bin/env python3
"""lh_seed_calc_campaign.py — setup for the lh-calc-campaign-numbers task."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lh_seed import run  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(run("lh-calc-campaign-numbers", sys.argv))
