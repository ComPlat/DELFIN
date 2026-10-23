#!/usr/bin/env python3
"""Setup for gen_a_fix_claim_whose_test_passed_before: the false claim.

See a_fix_claim_repo.py. argv[1] is the workspace.
"""
import sys

from a_fix_claim_repo import build

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        sys.exit(2)
    sys.exit(build(sys.argv[1], "false_claim"))
