#!/usr/bin/env python3
"""Acceptance for gen_a_fix_claim_whose_test_passed_before: expects
``refuted``. See fix_claim_verdict.py. argv[1] is the workspace."""
import sys

from fix_claim_verdict import check

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        sys.exit(2)
    sys.exit(check(sys.argv[1], "refuted"))
