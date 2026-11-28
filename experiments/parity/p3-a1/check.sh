#!/bin/sh
# One command for the check this package repeats: the personal-data scanner
# (the guard that caught the home path in after.json) plus the two red/green
# review tests. Run from anywhere; the tree is the one this file lives in.
set -e
cd "$(dirname "$0")/../../.." || exit 1
PY="${PY:-python}"
PYTHONPATH=. "$PY" -m pytest \
  tests/test_no_personal_data_is_committed.py \
  tests/test_a_review_finding_carries_its_location_and_its_failure.py \
  tests/test_a_review_reports_only_findings_it_could_verify.py \
  -q --timeout=120 -p no:randomly
